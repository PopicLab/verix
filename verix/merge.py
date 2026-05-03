from collections import defaultdict
import json
import logging
import networkx as nx
from networkx.algorithms.components import connected_components

from verix.io import write_merge_vcf
from verix.sv import BreakpointAlignment, Callset


class SVCluster:
    def __init__(self, cluster_id, svs_by_sample, samples):
        self.cluster_id = cluster_id
        self.svs_by_sample = svs_by_sample
        self.samples = samples
        self.representative = self.get_representative()
        self.support_vec = self.compute_support()
        self.size = sum(len(v) for v in svs_by_sample.values())

    def get_representative(self):
        return max((sv for s in self.svs_by_sample for sv in self.svs_by_sample[s]), key=lambda x: len(x.bkps))

    def compute_support(self):
        return [len(self.svs_by_sample[src]) if src in self.svs_by_sample else 0 for src in self.samples]

    def serialize_records(self, sample):
        if sample not in self.svs_by_sample: return "."
        return "|".join(f'{sv.id},{sv.type},{sv.breakpoints2str()}' for sv in self.svs_by_sample[sample])

    def to_vcf_info_dict(self):
        info = {'SUPPORT': len(self.svs_by_sample),
                'SUPPORT_COUNT': ",".join(str(c) for c in self.support_vec),
                'SUPPORT_BINARY': "".join(str(int(c > 0)) for c in self.support_vec)}
        return info

    vcf_info_fields = {
        'SUPPORT': ('1', 'Integer', 'Number of distinct input samples supporting this SV'),
        'SUPPORT_COUNT': ('1', 'String', 'Count vector indicating how many variants were merged from each sample'),
        'SUPPORT_BINARY': ('1', 'String', 'Binary vector indicating sample support'),
    }

class MergeEngine:
    def __init__(self, svs, aligner):
        self.callsets = svs
        self.aligner = aligner
        self.sv_clusters = []
        self.sample_names = [c.sample_name for c in self.callsets]

    def find_sv_clusters(self):
        graph = self.build_consensus_graph()
        merge_components = connected_components(graph)
        sorted_components = [sorted(list(comp)) for comp in merge_components]
        #sorted_components.sort()
        for component_id, component in enumerate(sorted_components):
            sample2svs = defaultdict(list)
            for callset_id, sv_id in component:
                sv = self.callsets[callset_id].id2sv[sv_id]
                sample2svs[self.sample_names[callset_id]].append(sv)
            self.sv_clusters.append(SVCluster(component_id, sample2svs, self.sample_names))

    def write_vcf(self, filepath):
        write_merge_vcf(self.callsets, self.sv_clusters, self.sample_names, SVCluster.vcf_info_fields, filepath)
        logging.info(f"Wrote consensus VCF to {filepath}")

    def write_stats(self, filepath):
        stats = {
            "n_total_variants": sum(len(c) for c in self.callsets),
            "n_variants_in_sample": {c.sample_name: len(c) for c in self.callsets},
            "n_clusters": len(self.sv_clusters),
            "support_vec_types": list(set(",".join(str(v) for v in c.support_vec) for c in self.sv_clusters)),
        }
        if self.sv_clusters:
           stats.update({
               "max_cluster_size": max(c.size for c in self.sv_clusters),
               "min_cluster_size": min(c.size for c in self.sv_clusters)}),
        logging.info("Results:\n" + json.dumps(stats, indent=4))
        json.dump(stats, open(filepath, "w"), indent=4)
        logging.info(f"Wrote benchmark report to {filepath}")

    def build_consensus_graph(self):
        graph = nx.Graph()
        for i, callset1 in enumerate(self.callsets):
            graph.add_nodes_from((i, k) for k in callset1.id2sv)
            for j, callset2 in enumerate(self.callsets[i+1:], start=i+1):
                for sv in callset1.svs:
                    for tid, bnd_matches in self.aligner.find_candidates(sv, callset2).items():
                        alignment = BreakpointAlignment(sv, callset2.id2sv[tid], self.aligner.align(bnd_matches))
                        if alignment.num_unmatched() != 0: continue
                        graph.add_edge((i, sv.id), (j, tid))
        return graph
