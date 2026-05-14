import pandas as pd
from collections import defaultdict
from networkx import Graph
from networkx.algorithms.matching import min_weight_matching
import bisect

class Breakpoint:
    __slots__ = ('chrom', 'pos', 'id')
    def __init__(self, chrom, pos, svid):
        self.chrom = chrom
        self.pos = pos
        self.id = svid

    def __eq__(self, other):
        return self.chrom == other.chrom and self.pos == other.pos and self.id == other.id

    def __hash__(self):
        return hash((self.chrom, self.pos, self.id))

    def __str__(self):
        return f"{self.chrom}:{self.pos}"

class SV:
    def __init__(self, svid, svtype, bkps, genotype, sample_name, records):
        self.id = svid
        self.genotype = genotype
        self.sample_name = sample_name  # VCF file
        self.records = records  # linked VCF records
        # consolidate breakpoints and types
        self.bkps = bkps
        self.type = svtype
        self.start = self.bkps[0]
        self.end = self.bkps[-1]

    @classmethod
    def from_records(cls, svid, records, vcf_bp, types, genotype, sample_name, bp_merge_threshold):
        return cls(
            svid=svid,
            svtype=cls.consolidate_type(types),
            bkps=cls.consolidate_breakpoints(vcf_bp, bp_merge_threshold),
            genotype=genotype,
            sample_name=sample_name,
            records=records,
        )

    @staticmethod
    def consolidate_type(types):
        return '+'.join(sorted(types))

    @staticmethod
    def consolidate_breakpoints(bps, bp_merge_threshold):
        groups = []
        for bp in sorted(bps, key=lambda b: (b.chrom, b.pos)):
            if (groups and bp.chrom == groups[-1][-1].chrom and
                    bp.pos - groups[-1][-1].pos <= bp_merge_threshold): groups[-1].append(bp)
            else: groups.append([bp])
        return [g[len(g) // 2] for g in groups]

    def get_min_max_size(self):
        min_size, max_size = 0, 0
        for b1, b2 in zip(self.bkps, self.bkps[1:]):
            if b1.chrom != b2.chrom: continue
            interval_len = abs(b1.pos - b2.pos)
            if not min_size or interval_len < min_size:
                min_size = interval_len
            if not max_size or interval_len >= max_size:
                max_size = interval_len
        return min_size, max_size

    def breakpoints2str(self):
        return ','.join([f'{bnd.chrom}:{bnd.pos}' for bnd in self.bkps])

    def __str__(self):
        return f"SV({self.sample_name}, {self.id}, {self.type}, {self.breakpoints2str()})"

class Callset:
    def __init__(self, svs, name=None, has_gt=False):
        self.svs = svs
        self.id2sv = {sv.id: sv for sv in self.svs}
        self.all_bnds = defaultdict(list)  # {chrom: [Breakend, ...] sorted by pos}
        self.all_bnds_pos = defaultdict(list)  # {chrom: [pos, ...]} parallel to all_bnds
        self.chrom_set = set()
        self.has_gt = has_gt
        self.sample_name = name
        self._index_breakpoints()

    def __len__(self):
        return len(self.svs)

    def _index_breakpoints(self):
        for sv in self.svs:
            for b in sv.bkps: self.all_bnds[b.chrom].append(b)
        self.chrom_set = set(self.all_bnds.keys())
        for chrom, chrom_bkps in self.all_bnds.items():
            chrom_bkps.sort(key=lambda b: b.pos)
            self.all_bnds_pos[chrom] = [b.pos for b in chrom_bkps]

    def lookup_breakpoints(self, chrom, start, end):
        chrom_positions = self.all_bnds_pos[chrom]
        lwr = bisect.bisect_left(chrom_positions, start)
        upr = bisect.bisect_right(chrom_positions, end)
        return self.all_bnds[chrom][lwr:upr]

class BreakpointAlignment:
    def __init__(self, sv_query, sv_target, assignment):
        self.sv_query = sv_query
        self.sv_target = sv_target
        self.assignment = assignment

    @property
    def num_matched(self):
        return len(self.assignment)

    @property
    def distance(self):
        return sum(abs(b1.pos - b2.pos) for b1, b2 in self.assignment.items())

    def num_unmatched(self, side=None):
        if side == "query": return len(self.sv_query.bkps) - self.num_matched
        if side == "target": return len(self.sv_target.bkps) - self.num_matched
        return len(self.sv_query.bkps) + len(self.sv_target.bkps) - 2 * self.num_matched

    def coverage(self, side=None):
        return "full" if self.num_unmatched(side) == 0 else "partial"

    def contiguous(self, side=None):
        all_bp = self.sv_target.bkps if side == "target" else self.sv_query.bkps
        matched_bp = self.assignment.values() if side == "target" else self.assignment.keys()
        if len(matched_bp) < 2: return False
        by_chrom = defaultdict(list)
        matched_by_chrom = defaultdict(list)
        for bp in all_bp: by_chrom[bp.chrom].append(bp)
        for bp in matched_bp: matched_by_chrom[bp.chrom].append(bp)
        multi = {c: m for c, m in matched_by_chrom.items() if len(m) >= 2}
        if not multi: return False
        for chrom, m in multi.items():
            idx = sorted(by_chrom[chrom].index(bp) for bp in m)
            if idx != list(range(idx[0], idx[0] + len(idx))): return False
        return True

    def score(self, side=None):
        return self.num_unmatched(side), self.distance

    def __str__(self, side=None):
        return ",".join(str(x) for x in (self.num_matched, self.distance,
                                         self.coverage(side="target"),
                                         self.contiguous(side="target"),
                                         ",".join(f"{p}-{t}" for p, t in self.assignment.items())))


class BreakpointAligner:
    def __init__(self, thresh=500, enforce_type=False, enforce_genotype=False):
        self.thresh = thresh
        self.enforce_type = enforce_type
        self.enforce_genotype = enforce_genotype

    def find_candidates(self, sv_query, target_callset):
        candidates = defaultdict(list)
        for b in sv_query.bkps:
            for target_b in target_callset.lookup_breakpoints(b.chrom, b.pos - self.thresh, b.pos + self.thresh):
                gt_sv = target_callset.id2sv[target_b.id]
                if self.enforce_type and gt_sv.type != sv_query.type: continue
                if self.enforce_genotype and gt_sv.genotype != sv_query.genotype: continue
                candidates[target_b.id].append((b, target_b))
        return candidates

    @staticmethod
    def align(bnd_matches):
        # min-weight matching over the given candidate breakend pairs
        edges = [(('a', c), ('b', g), abs(c.pos - g.pos) + 1)  # +1 keeps weights > 0
                 for c, g in bnd_matches]
        graph = Graph()
        graph.add_weighted_edges_from(edges)
        assignment = {}
        for n1, n2 in min_weight_matching(graph):
            call_node, gt_node = (n1, n2) if n1[0] == 'a' else (n2, n1)
            assignment[call_node[1]] = gt_node[1]
        return assignment