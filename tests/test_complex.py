import pytest
from verix.sv import Breakend, SV, Callset
from verix.bench import BenchCallset


# Helper to create a list of Breakends
def bnd_factory(chrom, positions, parent_id):
    return [Breakend(chrom, p, parent_id) for p in positions]


def test_defragmentation():
    """
    Scenario 1: One GT SV is split into three separate Call SVs.
    Tests: fragmentation flags.
    """
    gt_bnds = bnd_factory("chr1", [100, 200, 300, 400, 500, 600], "GT_1")
    gt_sv = SV("GT_1", "BND", gt_bnds)

    call1 = SV("C1", "BND", bnd_factory("chr1", [100, 200], "C1"))
    call2 = SV("C2", "BND", bnd_factory("chr1", [300, 400], "C2"))
    call3 = SV("C3", "BND", bnd_factory("chr1", [500, 600], "C3"))

    comp_callset = Callset([call1, call2, call3])
    gt_callset = BenchCallset([gt_sv], comp_callset, thresh=10)

    gt_callset.collect_candidate_detections()

    # GT_1 should be a 'fragmented' match because it's covered by multiple calls
    assert gt_sv.detection_mode == 'partial'
    assert gt_sv.optimal_candidate.fragmented

    # Each call should see GT_1 as fragmented relative to the total match
    for call in [call1, call2, call3]:
        assert call.optimal_candidate.gt_sv.parent_id == "GT_1"
        assert call.optimal_candidate.fragmented is True


def test_aggregation():
    """
    Scenario 1.5: One GT SV is split into three separate Call SVs.
    Tests: Aggregate detection mode.
    """
    call_bnds = bnd_factory("chr1", [100, 200, 300, 400, 500, 600], "C1")
    call_sv = SV("C1", "BND", call_bnds)

    gt1 = SV("GT1", "BND", bnd_factory("chr1", [100, 200], "GT1"))
    gt2 = SV("GT2", "BND", bnd_factory("chr1", [300, 400], "GT2"))
    gt3 = SV("GT3", "BND", bnd_factory("chr1", [500, 600], "GT3"))

    comp_callset = Callset([call_sv])
    gt_callset = BenchCallset([gt1, gt2, gt3], comp_callset, thresh=10)

    gt_callset.collect_candidate_detections()

    # GT_1 should be a 'fragmented' match because it's covered by multiple calls
    assert call_sv.detection_mode == 'aggregate'
    assert not call_sv.optimal_candidate.fragmented

    # Each call should see call_1 as fragmented relative to the total match
    for gt in gt_callset.svs:
        assert gt.optimal_candidate.call_sv.parent_id == "C1"
        assert gt.optimal_candidate.mode == 'aggregate'
        assert not gt.optimal_candidate.fragmented


def test_spurious_breakends_and_jitter():
    """
    Scenario 2: Caller adds an extra breakend and has slight coordinate shift.
    Tests: Spurious count and jitter calculation.
    """
    gt_sv = SV("GT_INV", "INV", bnd_factory("chr1", [1000, 2000], "GT_INV"))
    # Call has jitter (+5) and an extra BND at 1500
    call_sv = SV("CALL_INV", "INV", bnd_factory("chr1", [1005, 1500, 2005], "CALL_INV"))

    comp_callset = Callset([call_sv])
    gt_callset = BenchCallset([gt_sv], comp_callset, thresh=20)

    gt_callset.collect_candidate_detections()

    # Total distance should be |1000-1005| + |2000-2005| = 10
    assert call_sv.optimal_candidate.total_dist == 10
    # The breakend at 1500 doesn't match anything in GT
    assert call_sv.spurious == 1
    assert call_sv.detection_mode == 'partial'


def test_discontiguous_skipping():
    """
    Scenario 3: Caller captures the ends but misses the middle of a complex SV.
    Tests: The 'is_contiguous' logic when breakends are skipped.
    """
    # GT has 4 breakends
    gt_bnds = bnd_factory("chr1", [100, 200, 300, 400], "GT_CX")
    gt_sv = SV("GT_CX", "BND", gt_bnds)

    # Call captures 100 and 400, skipping the middle two
    call_sv = SV("CALL_CX", "BND", bnd_factory("chr1", [100, 400], "CALL_CX"))

    comp_callset = Callset([call_sv])
    gt_callset = BenchCallset([gt_sv], comp_callset, thresh=5)

    gt_callset.collect_candidate_detections()

    # Since call_idx 0 (pos 100) matches and call_idx 1 (pos 400) matches,
    # but in GT they are indices 0 and 3, next_match_idx != prev_match_idx + 1
    # and call_idx != prev_gt + 1 logic should trigger is_contiguous = False.
    assert call_sv.optimal_candidate.contiguity == "Discontiguous"


def test_nested_variants_best_match():
    """
    Scenario 4: Two GT variants close together, one Call matches one better.
    Tests: Branch and bound/scoring picking the better candidate.
    """
    gt_del = SV("GT_DEL", "DEL", bnd_factory("chr2", [500, 600], "GT_DEL"))
    gt_ins = SV("GT_INS", "INS", bnd_factory("chr2", [510, 520], "GT_INS"))

    # Call matches the Deletion perfectly, but is near the Insertion
    call_sv = SV("CALL_DEL", "DEL", bnd_factory("chr2", [500, 600], "CALL_DEL"))

    comp_callset = Callset([call_sv])
    gt_callset = BenchCallset([gt_del, gt_ins], comp_callset, thresh=50)

    gt_callset.collect_candidate_detections()

    # The score (matched_bnds, -dist) should favor GT_DEL (2 matches, 0 dist)
    # over GT_INS (0 matches or high dist)
    assert call_sv.optimal_candidate.gt_sv.parent_id == "GT_DEL"
    assert call_sv.detection_mode == "complete"


def test_deduplicate_near_breakends():
    """
    Tests the SV class internal logic for merging near-duplicate breakends.
    """
    # Two breakends at 100 and 101 with threshold 5 should merge
    bnds = bnd_factory("chr1", [100, 101, 500], "SV_1")
    sv = SV("SV_1", "BND", bnds, bp_merge_threshold=5)
    sv.update_breakends()

    # Should result in 2 breakends (one from the [100, 101] group and 500)
    assert len(sv.breakends) == 2
    assert sv.breakends[0].pos == 101  # middle-most of [100, 101] is index 1


def test_thresholds_both_pass():
    """
    Scenario 5: Call matches GT partially, but comfortably inside both thresholds.
    GT: 5 breakends. Call: 3 breakends.
    Matched: 3. Unmatched: 2 (missing from GT).
    Ratio: 3 / 5 = 0.6
    """
    gt_sv = SV("GT_LARGE", "BND", bnd_factory("chr1", [100, 200, 300, 400, 500], "GT_LARGE"))
    call_sv = SV("C_PART", "BND", bnd_factory("chr1", [200, 300, 400], "C_PART"))

    comp_callset = Callset([call_sv])

    # We tolerate up to 2 unmatched breakends and require at least 50% match ratio
    gt_callset = BenchCallset([gt_sv], comp_callset, thresh=10, unmatched_bnds_thresh=2, match_ratio_thresh=0.5)
    gt_callset.collect_candidate_detections()

    assert call_sv.optimal_candidate is not None
    assert call_sv.detection_mode == 'partial'


def test_thresholds_unmatched_fails():
    """
    Scenario 6: Match ratio is acceptable, but absolute unmatched count is too high.
    GT: 5 breakends. Call: 3 breakends.
    Matched: 3. Unmatched: 2. Ratio: 0.6.
    Requirement: max 1 unmatched breakend.
    """
    gt_sv = SV("GT_LARGE", "BND", bnd_factory("chr1", [100, 200, 300, 400, 500], "GT_LARGE"))
    call_sv = SV("C_PART", "BND", bnd_factory("chr1", [200, 300, 400], "C_PART"))

    comp_callset = Callset([call_sv])

    # Pruning should trigger due to unmatched_bnds_thresh=1
    gt_callset = BenchCallset([gt_sv], comp_callset, thresh=10, unmatched_bnds_thresh=1, match_ratio_thresh=0.5)
    gt_callset.collect_candidate_detections()

    assert call_sv.optimal_candidate is None
    assert call_sv.detection_mode == 'miss'


def test_thresholds_ratio_fails():
    """
    Scenario 7: Absolute unmatched count is fine, but match ratio is too low.
    GT: 5 breakends. Call: 3 breakends.
    Matched: 3. Unmatched: 2. Ratio: 0.6.
    Requirement: minimum 75% match ratio.
    """
    gt_sv = SV("GT_LARGE", "BND", bnd_factory("chr1", [100, 200, 300, 400, 500], "GT_LARGE"))
    call_sv = SV("C_PART", "BND", bnd_factory("chr1", [200, 300, 400], "C_PART"))

    comp_callset = Callset([call_sv])

    # Pruning should trigger due to match_ratio_thresh=0.75
    gt_callset = BenchCallset([gt_sv], comp_callset, thresh=10, unmatched_bnds_thresh=5,
                              match_ratio_thresh=0.75)
    gt_callset.collect_candidate_detections()

    assert call_sv.optimal_candidate is None
    assert call_sv.detection_mode == 'miss'


def test_thresholds_spurious_penalty():
    """
    Scenario 8: Call contains too many spurious breakends, violating the thresholds.
    GT: 3 breakends. Call: 5 breakends.
    Matched: 3. Unmatched: 2 (spurious in Call).
    Ratio: 3 / 5 = 0.6
    Requirement: max 1 unmatched breakend.
    """
    gt_sv = SV("GT_SMALL", "BND", bnd_factory("chr1", [100, 200, 300], "GT_SMALL"))
    call_sv = SV("C_LARGE", "BND", bnd_factory("chr1", [100, 200, 300, 400, 500], "C_LARGE"))

    comp_callset = Callset([call_sv])

    # The 2 spurious breakends should cause diff_bnds to exceed the threshold of 1
    gt_callset = BenchCallset([gt_sv], comp_callset, thresh=10, unmatched_bnds_thresh=1)
    gt_callset.collect_candidate_detections()

    assert call_sv.optimal_candidate is None
    assert call_sv.detection_mode == 'miss'


def test_thresholds_perfect_match_required():
    """
    Scenario 9: Strict requirement for perfect matching (ratio = 1.0, unmatched = 0).
    """
    gt_sv = SV("GT_EXACT", "BND", bnd_factory("chr1", [100, 200, 300], "GT_EXACT"))
    call_perfect = SV("C_PERF", "BND", bnd_factory("chr1", [100, 200, 300], "C_PERF"))
    call_flawed = SV("C_FLAW", "BND", bnd_factory("chr1", [100, 200], "C_FLAW"))

    comp_callset = Callset([call_perfect, call_flawed])

    # Require 100% match
    gt_callset = BenchCallset([gt_sv], comp_callset, thresh=10, unmatched_bnds_thresh=0, match_ratio_thresh=1.0)
    gt_callset.collect_candidate_detections()

    assert call_perfect.optimal_candidate is not None
    assert call_perfect.detection_mode == 'complete'

    assert call_flawed.optimal_candidate is None
    assert call_flawed.detection_mode == 'miss'
