import pytest
from verix.sv import Breakpoint, SV, Callset
from verix.bench import *
from verix.merge import *

def sv_factory(svid, svtype, chrom, positions):
    return SV(svid, svtype, [Breakpoint(chrom, p, svid) for p in positions], None, None, None)

def test_consolidate_types():
    assert SV.consolidate_type({"INV", "DUP"}) == "DUP+INV"
    assert SV.consolidate_type({"DEL", "DEL"}) == "DEL"

def test_consolidate_breakpoints():
    bps = [
        Breakpoint("chr1", 100, "x"),
        Breakpoint("chr1", 101, "x"),  # within thr=2 of 100
        Breakpoint("chr1", 102, "x"),  # within thr=2 of 100
        Breakpoint("chr1", 200, "x"),
        Breakpoint("chr2", 200, "x"),  # different chrom -> separate group
    ]
    merged = SV.consolidate_breakpoints(bps, bp_merge_threshold=2)
    # 100/101 collapse, chr1:200 alone, chr2:200 -> 3 breakpoints
    assert len(merged) == 3
    assert (merged[0].chrom, merged[0].pos) == ("chr1", 101)

def test_min_max_size():
    sv = sv_factory("s", "BND", "chr1", [100, 110, 200])
    sv.bkps.append(Breakpoint("chr2", 50, "s"))
    sv.end = sv.bkps[0]
    sv.end = sv.bkps[-1]
    mn, mx = sv.get_min_max_size()
    assert (mn, mx) == (10, 90)

def test_callset_lookup_breakpoints():
    sv = sv_factory("s", "DEL", "chr1", [100, 200, 300])
    cs = Callset([sv])
    assert {b.pos for b in cs.lookup_breakpoints("chr1", 150, 250)} == {200}
    assert {b.pos for b in cs.lookup_breakpoints("chr1", 100, 300)} == {100, 200, 300}
    assert cs.lookup_breakpoints("chrX", 0, 10 ** 9) == []

def test_alignment_distance_and_coverage():
    q = sv_factory("q", "DEL", "chr1", [100, 200])
    t = sv_factory("t", "DEL", "chr1", [105, 195])
    assignment = {q.bkps[0]: t.bkps[0], q.bkps[1]: t.bkps[1]}
    aln = BreakpointAlignment(q, t, assignment)
    assert aln.num_matched == 2
    assert aln.distance == 10
    assert aln.coverage("target") == "full"
    assert aln.coverage("query") == "full"

def test_classification(tmp_path):
    t1 = sv_factory("t1", "BND", "chr1", [100, 200, 300])
    t2 = sv_factory("t2", "BND", "chr1", [400, 500])
    q1 = sv_factory("q1", "BND", "chr1", [100, 200, 300])  # complete t1
    q2 = sv_factory("q2", "BND", "chr1", [100, 200, 300, 600])  # partial t1, due to 1 spurious
    q3 = sv_factory("q3", "BND", "chr1", [100, 200]) # partial t1
    q4 = sv_factory("q4", "BND","chr1", [100, 200, 600]) # partial t1 with 1 spurious
    q5 = sv_factory("q5", "BND", "chr1", [100, 200, 300, 400]) # aggregate (full, partial)
    q6 = sv_factory("q6", "BND", "chr1", [100, 200, 300, 400, 500]) # aggregate (full, full)
    q7 = sv_factory("q7", "BND", "chr1", [100, 200, 400]) # aggregate (partial, partial)
    q8 = sv_factory("q8", "BND", "chr1", [600, 700]) # spurious
    q9 = sv_factory("q9", "BND", "chr1", [100]) # partial, fragmented with q10
    q10 = sv_factory("q10", "BND", "chr1", [200, 300])  # partial, fragmented with q10
    engine = BenchmarkEngine(Callset([q1, q2, q3, q4, q5, q6, q7, q8, q9, q10]), Callset([t1, t2]),
                             BreakpointAligner(50, False, False))
    engine.find_matches()
    assert engine.get_match("q1").sv.id == "q1"
    assert len(engine.get_match("q1").sv.bkps) == 3
    assert engine.get_match("q1").match_class == MatchType.COMPLETE
    assert engine.get_match("q1").optimal.sv_target.id == "t1"
    assert engine.get_match("q1").optimal.distance == 0
    assert engine.get_match("q1").optimal.num_matched == 3
    assert engine.get_match("q1").spurious == 0
    assert engine.get_match("q1").fragmented is False
    assert len(engine.get_match("q1").alignments) == 1
    assert engine.get_match("q2").match_class == MatchType.PARTIAL
    assert engine.get_match("q2").optimal.sv_target.id == "t1"
    assert engine.get_match("q2").optimal.coverage("target") == "full"
    assert engine.get_match("q2").optimal.coverage("query") == "partial"
    assert engine.get_match("q2").spurious == 1
    assert engine.get_match("q2").fragmented is False
    assert len(engine.get_match("q2").alignments) == 1
    assert engine.get_match("q3").match_class == MatchType.PARTIAL
    assert engine.get_match("q3").optimal.sv_target.id == "t1"
    assert engine.get_match("q3").spurious == 0
    assert engine.get_match("q3").fragmented is True
    assert engine.get_match("q4").match_class == MatchType.PARTIAL
    assert engine.get_match("q4").optimal.sv_target.id == "t1"
    assert engine.get_match("q4").spurious == 1
    assert engine.get_match("q4").fragmented is True
    assert engine.get_match("q5").match_class == MatchType.AGGREGATE
    assert engine.get_match("q5").optimal.sv_target.id == "t1"
    assert engine.get_match("q5").optimal.coverage("target") == "full"
    assert engine.get_match("q5").spurious == 0
    assert engine.get_match("q5").fragmented is False
    assert len(engine.get_match("q5").alignments) == 2
    assert engine.get_match("q6").match_class == MatchType.AGGREGATE
    assert engine.get_match("q6").optimal.sv_target.id == "t1" # more breakpoints
    assert engine.get_match("q6").optimal.coverage("target") == "full"
    assert engine.get_match("q7").match_class == MatchType.AGGREGATE
    assert engine.get_match("q7").optimal.sv_target.id == "t1"
    assert engine.get_match("q7").optimal.coverage("target") == "partial"
    assert engine.get_match("q8").match_class == MatchType.SPURIOUS
    assert engine.get_match("q8").optimal is None
    assert engine.get_match("q9").match_class == MatchType.PARTIAL
    assert engine.get_match("q9").optimal.sv_target.id == "t1"
    assert engine.get_match("q9").fragmented is True
    assert engine.get_match("q10").match_class == MatchType.PARTIAL
    assert engine.get_match("q10").optimal.sv_target.id == "t1"
    assert engine.get_match("q10").fragmented is True

    stats_path = tmp_path / "report.json"
    engine.write_stats(stats_path)
    stats = json.loads(stats_path.read_text())
    assert stats["n_query"] == 10
    assert stats["n_target"] == 2
    assert stats["tp"] == 1
    assert stats["fp"] == 9
    assert stats["fn"] == 1
    assert stats["class_proportions"]["complete"] == pytest.approx(1 / 10)
    assert stats["class_proportions"]["partial"] == pytest.approx(5 / 10)
    assert stats["class_proportions"]["spurious"] == pytest.approx(1 / 10)
    assert stats["class_proportions"]["aggregate"] == pytest.approx(3 / 10)


def test_merge_single_cluster(tmp_path):
    a1 = sv_factory("a1", "BND", "chr1", [100, 200])
    a1_1 = sv_factory("a1_1", "BND", "chr1", [101, 205])
    a2 = sv_factory("a2", "BND", "chr1", [1000, 1100])
    a3 = sv_factory("a3", "BND", "chr1", [5000, 5100])  # singleton
    b1 = sv_factory("b1", "BND", "chr1", [105, 198])
    b2 = sv_factory("b2", "BND", "chr1", [1010, 1095])
    c1 = sv_factory("c1", "BND", "chr1", [102, 203])
    c2 = sv_factory("c2", "BND", "chr2", [100, 200])  # singleton (diff chrom)

    engine = MergeEngine([Callset([a1, a1_1, a2, a3], name="A"),
                          Callset([b1, b2], name="B"),
                          Callset([c1, c2], name="C")],
                         BreakpointAligner(thresh=50))
    engine.find_sv_clusters()

    # 4 clusters: {a1,a1_1,b1,c1}, {a2,b2}, {a3}, {c2}
    assert len(engine.sv_clusters) == 4
    a, b, *c = sorted(engine.sv_clusters, key=lambda c: -c.size)
    assert a.size == 4 and a.support_vec == [2, 1, 1]
    assert a.to_vcf_info_dict() == {
        "SUPPORT": 3, "SUPPORT_BINARY": "111", "SUPPORT_COUNT": "2,1,1",
    }
    assert b.serialize_records("B") == "b2,BND,chr1:1010,chr1:1095"
    assert b.size == 2 and b.support_vec == [1, 1, 0]
    assert b.to_vcf_info_dict()["SUPPORT_BINARY"] == "110"
    assert b.serialize_records("C") == "."
    assert {tuple(e.support_vec) for e in c} == {(1, 0, 0), (0, 0, 1)}

    stats_path = tmp_path / "report.json"
    engine.write_stats(stats_path)
    stats = json.loads(stats_path.read_text())
    assert stats["n_total_variants"] == 8
    assert stats["n_variants_in_sample"] == {"A": 4, "B": 2, "C": 2}
    assert stats["n_clusters"] == 4
    assert (stats["min_cluster_size"], stats["max_cluster_size"]) == (1, 4)
    assert set(stats["support_vec_types"]) == {"2,1,1", "1,1,0", "1,0,0", "0,0,1"}

