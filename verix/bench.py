from enum import Enum
import logging
from collections import defaultdict, Counter
import json
import numpy as np
import pandas as pd
from tqdm import tqdm
from typing import Any

from verix.plot import BenchPlotter
from verix.sv import *
from verix.io import write_bench_vcf


class MatchType(str, Enum):
    COMPLETE = "complete"
    PARTIAL = "partial"
    AGGREGATE = "aggregate"
    SPURIOUS = "spurious"

class MatchAnnotation:
    def __init__(self, sv):
        self.sv = sv
        self.alignments = []
        self.optimal = None
        self.match_class = None
        self.spurious = 0
        self.fragmented = False

    def add_alignment(self, aln):
        self.alignments.append(aln)
        if self.optimal is None or aln.score(side="query") < self.optimal.score(side="query"):
            self.optimal = aln

    def classify(self):
        if self.optimal is None:
            self.match_class = MatchType.SPURIOUS
            self.spurious = len(self.sv.bkps)
            return
        if self.optimal.num_unmatched() == 0:
            self.match_class = MatchType.COMPLETE
            return

        bnd_in_opt = set(self.optimal.assignment)
        bnd_in_other = {
            b for aln in self.alignments
            if aln.sv_target.id != self.optimal.sv_target.id
            for b in aln.assignment
        }
        self.match_class = MatchType.AGGREGATE if (bnd_in_other - bnd_in_opt) else MatchType.PARTIAL
        self.spurious = sum(1 for b in self.sv.bkps if b not in (bnd_in_other | bnd_in_opt))

    def to_vcf_info_dict(self):
        info = {'BEST_MATCH_CLASS': self.match_class.value}
        if self.spurious: info['SPURIOUS'] = self.spurious
        if self.fragmented: info['FRAGMENTED'] = self.fragmented
        if self.optimal:
            info['BEST_MATCH_ID'] = str(self.optimal.sv_target.id)
            info['BEST_MATCH_COV'] = self.optimal.coverage(side="target")
            info['BEST_MATCH_TYPE'] = self.optimal.sv_target.type
            info['BEST_N_MATCHED'] = self.optimal.num_matched
            info['BEST_BND_DIST'] = self.optimal.distance
        if self.alignments:
            info['MATCHES'] = "|".join(f'{aln.sv_target.id},{aln.sv_target.type},{str(aln)}' for aln in self.alignments)
        return info

    def to_stats_info_dict(self):
        info = self.to_vcf_info_dict()
        info.update(QID=self.sv.id, QTYPE=self.sv.type, QNBKPS=len(self.sv.bkps),
                    TNBKPS=len(self.optimal.sv_target.bkps) if self.optimal else None,
                    NTARGETS=len({aln.sv_target.id for aln in self.alignments}))
        return info

    vcf_info_fields = {
        'BEST_MATCH_CLASS': ('1', 'String', 'Match classification: complete | partial | aggregate | spurious'),
        'BEST_MATCH_ID': ('1', 'String', 'ID of the best matching target event'),
        'BEST_MATCH_COV': ('1', 'String', 'full iff all target event breakpoints are covered, else partial'),
        'BEST_MATCH_TYPE': ('1', 'String', 'SV type of the best matching target event'),
        'BEST_N_MATCHED': ('1', 'Integer', 'Number of breakpoints matched in the best alignment'),
        'BEST_BND_DIST': ('1', 'Integer', 'Total breakpoints distance in the best alignment'),
        'MATCHES': ('.', 'String', 'All candidate alignments'),
        'SPURIOUS': ('1', 'Integer', 'Number of query breakpoints that match no target event'),
        'FRAGMENTED': ('0', 'Flag', 'Best target match was covered by other query calls'),
    }
    stats_info_fields = list(vcf_info_fields.keys()) + ["QID", "QTYPE", "QNBKPS", "TNBKPS", "NTARGETS"]

class BenchmarkEngine:
    # Compares a query SV callset against a target truthset
    def __init__(self, query, target, aligner):
        self.query = query
        self.target = target
        self.matches = {sv.id: MatchAnnotation(sv) for sv in query.svs}
        self.aligner = aligner
        self.target_matches = defaultdict(set)

    def get_match(self, svid):
        return self.matches[svid]

    def find_matches(self):
        logging.info("Finding matches...")
        for sv in tqdm(self.query.svs):
            self.match_to_target(sv)
        self.global_annotation()

    def match_to_target(self, sv):
        match = self.matches[sv.id]
        # 1. find candidate matches in target truthset
        candidates = self.aligner.find_candidates(sv, self.target)
        # 2. compute candidate alignments
        for target_id, bp_matches in candidates.items():
            aln = BreakpointAlignment(sv, self.target.id2sv[target_id], self.aligner.align(bp_matches))
            match.add_alignment(aln)
            self.target_matches[target_id].update(aln.assignment.values())
        # 3. classify the match based on optimal alignments + other matches
        match.classify()

    def global_annotation(self):
        for match in self.matches.values():
            if match.match_class in [MatchType.SPURIOUS, MatchType.COMPLETE]: continue
            optimal = set(match.optimal.assignment.values())
            extra_matches = self.target_matches[match.optimal.sv_target.id] - optimal
            if extra_matches:
                match.fragmented = True

    def compute_stats(self):
        n_query = len(self.query.svs)
        n_target = len(self.target.svs)
        stats: dict[str, Any] = {
            "n_query": n_query,
            "n_target": n_target,
        }
        matches = self.matches.values()
        # ---- per-class stats
        class2query = {cls: [m for m in matches if m.match_class == cls] for cls in MatchType}
        class2targets = {cls: {m.optimal.sv_target for m in ms}
                         for cls, ms in class2query.items() if cls != MatchType.SPURIOUS}
        query_type_totals = Counter(s.type for s in self.query.svs)
        target_type_totals = Counter(s.type for s in self.target.svs)

        # ---- TP, FN, FP for complete matches
        if class2query[MatchType.COMPLETE]:
            tp_query = len(class2query[MatchType.COMPLETE])
            tp_target = len(class2targets[MatchType.COMPLETE])
            fp = n_query - tp_query
            fn = n_target - tp_target
            precision = tp_query / n_query
            recall = tp_target / n_target
            stats.update(tp_query=tp_query, tp_target=tp_target, fp=fp, fn=fn,
                         precision=precision, recall=recall,
                         f1=2 * precision * recall / (precision + recall))

        stats_by_class = {}
        for cls in class2query:
            if cls in [MatchType.SPURIOUS]: continue
            if not class2query[cls]: continue
            cls_matches = class2query[cls]
            cls_targets = class2targets[cls]
            qtype_counts = Counter(m.sv.type for m in cls_matches)
            ttype_counts = Counter(s.type for s in cls_targets)
            stats_by_class[cls] = {
                "num_matches": len(cls_matches),
                "num_unique_targets": len(cls_targets),
                "mean_breakpoint_distance": float(np.mean([m.optimal.distance / m.optimal.num_matched for m in cls_matches])),
                "mean_breakpoint_hit_rate": float(np.mean([m.optimal.num_matched / len(m.optimal.sv_target.bkps) for m in cls_matches])),
                "mean_spurious_breakpoint_rate": float(np.mean([m.spurious / len(m.sv.bkps) for m in cls_matches])),
                "mean_targets_per_record": float(np.mean([len({aln.sv_target.id for aln in m.alignments}) for m in cls_matches])),
                "query_type_counts": dict(qtype_counts.most_common()),
                "query_type_proportions": {t: c / query_type_totals[t] for t, c in qtype_counts.most_common()},
                "target_type_counts": dict(ttype_counts.most_common()),
                "target_type_proportions": {t: c / target_type_totals[t] for t, c in ttype_counts.most_common()},
            }
        stats["class_proportions"] = {cls.value: len(class2query[cls]) / n_query for cls in MatchType}
        stats["by_class"] = {cls.value: stats_by_class[cls]
                             for cls in (MatchType.COMPLETE, MatchType.PARTIAL, MatchType.AGGREGATE)
                             if cls in stats_by_class}
        return stats

    def write_stats(self, filepath):
        stats = self.compute_stats()
        logging.info("Results:\n" + json.dumps(stats, indent=4))
        json.dump(stats, open(filepath, "w"), indent=4)
        logging.info(f"Wrote benchmark report to {filepath}")

    def write_vcf(self, filepath):
        write_bench_vcf(self.query, self.matches.values(), MatchAnnotation.vcf_info_fields, filepath)
        logging.info(f"Wrote benchmark VCF to {filepath}")

    def generate_plots(self, filepath):
        match_df = pd.DataFrame([m.to_stats_info_dict() for m in self.matches.values()],
                          columns=MatchAnnotation.stats_info_fields)
        target_hits = {m.optimal.sv_target.id: m.optimal.coverage()
                       for m in self.matches.values() if m.match_class != MatchType.SPURIOUS}
        target_df = pd.DataFrame([
            {"ID": str(s.id), "TYPE": s.type, "COV": target_hits.get(s.id, "miss"), "NBP": len(s.bkps),
             "UNION_COV": "miss" if s.id not in self.target_matches else \
                 "full" if len(self.target_matches[s.id]) == len(s.bkps) else "partial"}
            for s in self.target.svs])
        plotter = BenchPlotter(filepath, match_df, target_df, self.compute_stats())
        plotter.make_plots()
        logging.info(f"Generated plots in: {filepath}")



