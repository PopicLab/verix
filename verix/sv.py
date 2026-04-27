import bisect
import pandas as pd
from collections import defaultdict
from dataclasses import dataclass
from typing import Tuple
from networkx import Graph
from networkx.algorithms.matching import min_weight_matching

from verix.vcf_writer import stringify_matched_bnds, stringify_bnd


class Breakend:
    # Memory optimization
    __slots__ = ('chrom', 'pos', 'parent_id', 'bnd_matches', 'match_bp')

    def __init__(self, chrom, pos, parent_id):
        # initialized during parsing of caller output
        self.chrom = chrom
        self.pos = int(pos)
        self.parent_id = parent_id  # <- store PID in the breakend obj so we can look up the SV it belongs to
        self.bnd_matches = defaultdict(list)
        self.match_bp = None

    def __eq__(self, other):
        if not isinstance(other, Breakend):
            return False
        return (self.chrom == other.chrom and
                self.pos == other.pos and
                self.parent_id == other.parent_id)

    def __hash__(self):
        return hash((self.chrom, self.pos, self.parent_id))

    def __repr__(self):
        return f"Breakend({self.chrom}, {self.pos}, {self.parent_id})"


class SV:
    def __init__(self, parent_id, parent_type, breakends, genotype={}, bp_merge_threshold=1, source=None, records=()):
        # initialized for each output record of a comparison caller output callset
        self.parent_id = parent_id
        self.parent_type = parent_type
        self.bp_merge_threshold = bp_merge_threshold
        self.breakends = breakends
        self.start = self.breakends[0]
        self.end = self.breakends[-1]
        self.genotype = genotype

        # default detection label to be overwritten if complete match
        self.detection_mode = 'miss'
        self.match_type = []  # <- if the SV finds a match, store the type
        self.match_id = []  # <- if the SV finds a match, store the ID
        self.match_score = []
        self.optimal_candidate = None
        self.match_list = []
        self.spurious = 0
        self.optimal_score = ()
        self.candidates = []

        self.bnd2idx = {}

        # Keep track of the VCF records
        self.records = records

        # Keep track of the source VCF for merge
        self.source = source

    def update_breakends(self):
        self.breakends = sorted(self.breakends, key=lambda b: (b.chrom, b.pos))
        self.deduplicate_breakends()
        self.end = self.breakends[-1]
        self.start = self.breakends[0]

    def deduplicate_breakends(self):
        # remove near-duplicate breakends reported within the same SV
        # -> needed esp. for Mako, which often reports multiple fragments with
        # -> near-duplicate breakends for complex SVs
        groups = []
        prev_bp = None
        gp = []
        # collect groups of breakend separated by a distance < threshold
        for bp in self.breakends:  # <- traverse sorted self.breakends in order of position
            if not gp or (bp.chrom == prev_bp.chrom and bp.pos - prev_bp.pos < self.bp_merge_threshold):
                gp.append(bp)
            else:
                groups.append(gp)
                gp = [bp]
            prev_bp = bp
        groups.append(gp)
        # for each group, retain the middle-most breakend
        filtered_bps = [g[len(g) // 2] for g in groups]
        self.breakends = filtered_bps


@dataclass
class SVMatch:
    call_sv: SV
    gt_sv: SV
    mode: str
    contiguity: str
    diff_bnd: int
    total_dist: int
    matched_bnds: Tuple
    fragmented: bool = False


class Callset:
    def __init__(self, svs):
        self.svs = svs
        self.pid2sv = {sv.parent_id: sv for sv in self.svs}

        self.all_bnds = defaultdict(list)
        self.all_bnds_pos = defaultdict(list)  # <- dict of the form {chromosome: [breakends]}
        self.chrom_set = set()
        self.sample_names = set()

        self.get_all_bnds()

    def callset2dataframe(self, is_gt=False):
        # Create dataframe for downstream analysis
        rows = []

        for sv in self.svs:
            span = sv.end.pos - sv.start.pos if sv.breakends[0].chrom == sv.breakends[-1].chrom else -1

            row = {
                'parent_id': str(sv.parent_id),
                'parent_type': str(sv.parent_type),
                'num_breakends': len(sv.breakends),
                'detection_mode': sv.detection_mode,
                'sv_span': span,
                'spurious': sv.spurious,
                'num_candidates': len(sv.candidates),
                'breakends': "|".join([stringify_bnd(bnd) for bnd in sv.breakends]),
                'matched_breakends': "|".join(
                    [stringify_bnd(bnd) for bnd in sv.breakends if bnd.bnd_matches]) or None,
                'match_type': ";".join(map(str, sv.match_type)) if sv.match_type else None,
                'match_id': ";".join(map(str, sv.match_id)) if sv.match_id else None,
                'match_num_unmatched_bnds': ";".join(str(score[0]) for score in sv.match_score) if sv.match_score else None,
                'match_total_distance': ";".join(str(score[1]) for score in sv.match_score) if sv.match_score else None,
            }

            opt = sv.optimal_candidate
            if opt:
                other_sv = opt.call_sv if is_gt else opt.gt_sv
                row.update({
                    'optimal_candidate_id': str(other_sv.parent_id),
                    'optimal_contiguity': opt.contiguity,
                    'fragmented': opt.fragmented,
                    'optimal_type': other_sv.parent_type,
                    'optimal_breakends': stringify_matched_bnds(opt.matched_bnds),
                    'optimal_num_unmatched': opt.diff_bnd,
                    'optimal_total_distance': opt.total_dist,
                    'optimal_full': len(opt.matched_bnds) == len(opt.gt_sv.breakends),

                    # Store raw counts temporarily for vectorized ratio calculation
                    '_opt_match_count': len(opt.matched_bnds),
                    '_gt_bnd_count': len(opt.gt_sv.breakends),
                    '_call_bnd_count': len(opt.call_sv.breakends)
                })
            else:
                row.update({
                    'optimal_candidate_id': None,
                    'optimal_contiguity': None,
                    'fragmented': False,
                    'optimal_type': None,
                    'optimal_mode': None,
                    'optimal_breakends': None,
                    'optimal_num_unmatched': None,
                    'optimal_total_distance': None,
                    '_opt_match_count': 0, '_gt_bnd_count': 1, '_call_bnd_count': 1  # Avoid div by zero
                })

            rows.append(row)

        df = pd.DataFrame(rows)
        if not df.empty:
            df['optimal_miss'] = 1 - df['_opt_match_count'] / df['_gt_bnd_count']
            df['optimal_extra'] = 1 - df['_opt_match_count'] / df['_call_bnd_count']
            df.drop(columns=['_opt_match_count', '_gt_bnd_count', '_call_bnd_count'], inplace=True)
        return df

    def get_all_bnds(self):
        for sv in self.svs:
            sv.breakends.sort(key=lambda bnd: (bnd.chrom, bnd.pos))
            for bnd_idx, bnd in enumerate(sv.breakends):
                self.all_bnds[bnd.chrom].append(bnd)
                sv.bnd2idx[bnd] = bnd_idx
            self.sample_names = self.sample_names.union(sv.genotype.keys())
        self.chrom_set = set(self.all_bnds.keys())

        for chrom, chrom_bnds in self.all_bnds.items():
            chrom_bnds.sort(key=lambda bnd: bnd.pos)
            self.all_bnds_pos[chrom] = [bnd.pos for bnd in chrom_bnds]


class MatchCallset(Callset):
    def __init__(self, svs, thresh=500, unmatched_bnds_thresh=None, match_ratio_thresh=None, enforce_svtype=False,
                 enforce_genotype=False):
        super().__init__(svs)

        self.thresh = thresh
        self.unmatched_bnds_thresh = unmatched_bnds_thresh
        self.match_ratio_thresh = match_ratio_thresh
        self.enforce_svtype = enforce_svtype
        self.enforce_genotype = enforce_genotype

    def find_candidates(self, sv, pid2sv, all_bnds, all_bnds_pos):
        candidates = set()
        for bnd in sv.breakends:
            # for each breakend, find the overlapping calls breakends
            list_bnds_chrom = all_bnds_pos[bnd.chrom]
            lwr = bisect.bisect_left(list_bnds_chrom, bnd.pos - self.thresh)
            upr = bisect.bisect_right(list_bnds_chrom, bnd.pos + self.thresh)
            bnd_matches = all_bnds[bnd.chrom][lwr:upr]
            for match_bnd in bnd_matches:
                candidate_sv = pid2sv[match_bnd.parent_id]

                if self.enforce_svtype and candidate_sv.parent_type != sv.parent_type: continue
                if self.enforce_genotype and candidate_sv.genotype != sv.genotype: continue

                bnd.bnd_matches[match_bnd.parent_id].append(match_bnd)
                candidates.add(match_bnd.parent_id)
        sv.candidates = candidates

    @staticmethod
    def align_candidates(sv, candidate_sv):
        # Use minimum weight matching to find the best alignments
        edges = []
        for bnd in sv.breakends:
            for match in bnd.bnd_matches[candidate_sv.parent_id]:
                # Remove ambiguities if the breakends are the same
                edges.append((('sv', bnd), ('candidate', match), abs(bnd.pos - match.pos) + 1))

        graph = Graph()
        graph.add_weighted_edges_from(edges)
        matching = min_weight_matching(graph)

        # Reformat the matching output and remove duplicates
        clean_matching = []
        seen = set()
        for node_1, node_2 in matching:
            n1_type, n1_bnd = node_1
            _, n2_bnd = node_2
            if n1_type == 'sv':
                if n1_bnd in seen: continue
                clean_matching.append((n1_bnd, n2_bnd))
                seen.add(n1_bnd)
            else:
                if n2_bnd in seen: continue
                clean_matching.append((n2_bnd, n1_bnd))
                seen.add(n2_bnd)
        return clean_matching

    def filter_candidate(self, unmatched_bnds, num_matched, total_num_bnds):
        if self.unmatched_bnds_thresh is not None and unmatched_bnds > self.unmatched_bnds_thresh:
            return True

        if self.match_ratio_thresh is not None and num_matched / total_num_bnds < self.match_ratio_thresh:
            return True
        return False