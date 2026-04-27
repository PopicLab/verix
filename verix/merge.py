from collections import defaultdict
from dataclasses import dataclass
import networkx as nx
from networkx.algorithms.components import connected_components
from collections import Counter

from verix.sv import MatchCallset, SV


@dataclass
class SVConsensus:
    merge_id: str
    representative_sv: SV
    support_vec: str
    support_count: int
    support_count_vec: str
    stringified_recs: dict[str, str]


def stringify_merge_rec(sv):
    # Process the records to be written in the merged record
    stringify_rec = []
    for rec in sv.records:
        chrom2 = rec.chrom
        if 'CHROM2' in rec.info:
            chrom2 = rec.info['CHROM2']
        stringify_rec.append(
            f"{rec.id}|{rec.chrom}:{rec.pos}|{chrom2}:{rec.stop}|{sv.parent_type}|{sv.parent_id}")
    return stringify_rec


def join_stringified_records(record_by_source):
    return {source: ','.join(record_list) for source, record_list in record_by_source.items()}


def compute_support(support, ordered_sources):
    # Compute support vector
    supp_vec = "".join(['1' if src in support else '0' for src in ordered_sources])
    supp_count = len(support)

    source_tally = Counter(support)
    supp_counts_vec = "-".join([str(source_tally[src]) for src in ordered_sources])
    return {'supp_vec': supp_vec, 'supp_count_vec': supp_counts_vec, 'supp_count': supp_count}


class ConsensusCallset(MatchCallset):
    def __init__(self, svs, thresh, ordered_sources, unmatched_bnds_thresh=None, match_ratio_thresh=None,
                 enforce_svtype=False, enforce_genotype=False):
        super().__init__(svs, thresh=thresh, unmatched_bnds_thresh=unmatched_bnds_thresh, match_ratio_thresh=match_ratio_thresh,
                         enforce_svtype=enforce_svtype, enforce_genotype=enforce_genotype)
        self.ordered_sources = ordered_sources
        self.merge_list = []

    def build_merge_graph(self):
        # Build a graph of overlapping SV fulfilling the merge constraints
        edge_list = []
        for sv in self.svs:
            self.find_candidates(sv, self.pid2sv, self.all_bnds, self.all_bnds_pos)

            for candidate_sv_id in sv.candidates:
                candidate_sv = self.pid2sv[candidate_sv_id]
                # prevent double counting
                if candidate_sv_id <= sv.parent_id:
                    continue

                matched_bnds = self.align_candidates(sv, candidate_sv)

                total_num_bnds_call = len(sv.breakends)
                total_num_bnds_candidate = len(candidate_sv.breakends)
                diff_bnd = total_num_bnds_call + total_num_bnds_candidate - 2 * len(matched_bnds)

                if self.filter_candidate(diff_bnd, len(matched_bnds), diff_bnd + len(matched_bnds)): continue
                edge_list.append((sv, candidate_sv))
        merge_graph = nx.Graph()
        merge_graph.add_nodes_from(self.svs)
        merge_graph.add_edges_from(edge_list)
        return merge_graph

    def compute_merge(self, merge_graph):
        merge_components = connected_components(merge_graph)

        # To ensure determinism of the representative_sv
        sorted_components = [
            sorted(list(comp), key=lambda sv: sv.parent_id)
            for comp in merge_components
        ]

        # To ensure determinism of the components order
        sorted_components.sort(key=lambda comp: comp[0].parent_id)

        # Define the representative SV for each connected component
        for component_id, component in enumerate(sorted_components):
            merge_id = f"merge_{component_id}"
            sources = []
            record_by_source = defaultdict(list)

            max_sv = None
            for sv in component:
                sources.append(sv.source)
                record_by_source[sv.source].extend(stringify_merge_rec(sv))

                if not max_sv or len(max_sv.breakends) < len(sv.breakends):
                    max_sv = sv

            stringified_recs = join_stringified_records(record_by_source)
            support = compute_support(sources, self.ordered_sources)
            self.merge_list.append(SVConsensus(
                merge_id=merge_id,
                representative_sv=max_sv,
                stringified_recs=stringified_recs,
                support_vec=support['supp_vec'],
                support_count=support['supp_count'],
                support_count_vec=support['supp_count_vec'],
            ))

    def merge(self):
        merge_graph = self.build_merge_graph()
        self.compute_merge(merge_graph)
