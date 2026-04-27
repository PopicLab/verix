from collections import Counter
from copy import deepcopy

from verix.sv import SV, Breakend
from verix.merge import ConsensusCallset


def create_sv(sv_id, sv_type, positions, source, chrom="chr1"):
    """
    Helper function to quickly generate SVs for testing.
    """
    bnds = [Breakend(chrom, pos, sv_id) for pos in positions]
    return SV(parent_id=sv_id, parent_type=sv_type, breakends=bnds, source=source)


def test_jittered_clustering_and_median():
    """
    Scenario 1: 3 callers find the same DEL, but with slight coordinate jitter.
    Tests: Do they merge into 1 cluster? Is the median correctly computed?
    """
    sv1 = create_sv("SV1", "DEL", [100, 200], "Severus")
    sv2 = create_sv("SV2", "DEL", [105, 195], "SVision")
    sv3 = create_sv("SV3", "DEL", [95, 205], "Sniffles")

    # Threshold 15 covers the max distance (105 to 95 = 10)
    merge = ConsensusCallset([sv1, sv2, sv3], thresh=15, ordered_sources=['Sniffles', 'Severus', 'SVision'])
    merge.merge()

    # All 3 should merge into exactly 1 merge SV
    assert len(merge.merge_list) == 1
    merge_sv = merge.merge_list[0]
    representative_sv = merge_sv.representative_sv
    # The medians of the groups:
    # Group 1: [95, 100, 105] -> Median is 100
    # Group 2: [195, 200, 205] -> Median is 200
    assert len(representative_sv.breakends) == 2
    assert representative_sv.breakends[0].pos == 100
    assert representative_sv.breakends[1].pos == 200

    # Check support tracking
    assert merge_sv.support_count == 3
    assert merge_sv.support_count_vec == '1-1-1'
    assert representative_sv.parent_type == "DEL"


def test_disjoint_svs_do_not_cluster():
    """
    Scenario 2: 2 callers find SVs, but they are too far apart to be the same event.
    Tests: Do they stay separated?
    """
    sv1 = create_sv("SV1", "DEL", [100, 200], "Severus")
    sv2 = create_sv("SV2", "DEL", [5000, 6000], "SVision")

    merge = ConsensusCallset([sv1, sv2], thresh=50, ordered_sources=['Severus', 'SVision'])
    merge.merge()

    # Should result in 2 separate merge events
    assert len(merge.merge_list) == 2


def test_complex_sv_clustering():
    """
    Scenario 3: Complex SVs with 3 breakends.
    Tests: Can the clustering and median logic handle N-breakend graphs?
    """
    sv1 = create_sv("SV_C1", "BND", [1000, 2000, 3000], "Severus")
    sv2 = create_sv("SV_C2", "BND", [1010, 1990, 3010], "SVision")

    # Threshold 20
    merge = ConsensusCallset([sv1, sv2], thresh=20, unmatched_bnds_thresh=0, ordered_sources=['Severus', 'SVision'])
    merge.merge()

    assert len(merge.merge_list) == 1
    merge_sv = merge.merge_list[0]
    representative_sv = merge_sv.representative_sv

    assert len(representative_sv.breakends) == 3
    assert representative_sv.breakends[0].pos == 1000
    assert representative_sv.breakends[1].pos == 2000
    assert representative_sv.breakends[2].pos == 3000

    # Check support tracking
    assert merge_sv.support_count == 2
    assert merge_sv.support_count_vec == '1-1'
    assert representative_sv.parent_type == "BND"


def test_multiple_svs_from_same_caller():
    """
    Scenario 5: One caller reports a fragmented SV as two records, another reports it as one.
    Tests: Counter functionality and partial matching threshold.
    """
    # Severus finds it as one big SV
    sv_severus = create_sv("SV_M", "DEL", [100, 200, 300], "Severus")

    # SVision finds it broken up
    sv_svision_1 = create_sv("SV_D1", "DEL", [100, 200], "SVision")
    sv_svision_2 = create_sv("SV_D2", "DEL", [200, 300], "SVision")

    # Allow 1 unmatched breakend (diff_bnd <= 1)
    merge = ConsensusCallset([sv_severus, sv_svision_1, sv_svision_2],
                                 thresh=10, unmatched_bnds_thresh=2, ordered_sources=['Severus', 'SVision'])
    merge.merge()

    # They should all merge into 1 cluster
    assert len(merge.merge_list) == 1
    merge_sv = merge.merge_list[0]
    representative_sv = merge_sv.representative_sv

    assert len(representative_sv.breakends) == 3
    assert representative_sv.breakends[0].pos == 100
    assert representative_sv.breakends[1].pos == 200
    assert representative_sv.breakends[2].pos == 300

    # Check support tracking
    assert merge_sv.support_count == 3
    assert merge_sv.support_count_vec == '1-2'
    assert representative_sv.parent_type == "DEL"
    representative_sv.parent_id = 'SV_M'


def test_large_sv_overlapping_multiple_smaller_svs():
    """
    Scenario 6: A massive complex SV overlaps with two smaller fragmentary SVs.
    Tests: unmatched_bnds_thresh and match_ratio_thresh correctly gate clustering.
    """
    # SV_LARGE: 7 breakends
    sv_large = create_sv("SV_LARGE", "BND", [100, 200, 300, 400, 500, 600, 700], "Manta")

    # SV_SMALL_1: 2 breakends.
    # Matched = 2.
    # diff_bnd = (7 + 2) - (2 * 2) = 5 unmatched breakends.
    # Jaccard ratio = matched / (diff_bnd + matched) = 2 / (5 + 2) = ~0.285
    sv_small_1 = create_sv("SV_S1", "BND", [200, 300], "Severus")

    # SV_SMALL_2: 4 breakends.
    # Matched = 4.
    # diff_bnd = (7 + 4) - (2 * 4) = 3 unmatched breakends.
    # Jaccard ratio = 4 / (3 + 4) = ~0.571
    sv_small_2 = create_sv("SV_S2", "BND", [400, 500, 600, 700], "SVision")

    # ---------------------------------------------------------
    # Configuration 1: Strict on unmatched_bnds_thresh (max 3)
    # SV_SMALL_1 (diff=5) -> REJECTED
    # SV_SMALL_2 (diff=3) -> ACCEPTED
    # ---------------------------------------------------------
    merge_1 = ConsensusCallset([deepcopy(sv_large), deepcopy(sv_small_1), deepcopy(sv_small_2)],
                                   thresh=10, unmatched_bnds_thresh=3, ordered_sources=['Manta', 'Severus', 'SVision'])
    merge_1.merge()
    # Expect 2 final SVs: {SV_LARGE + SV_S2} (size 2) and {SV_S1} (size 1)
    assert len(merge_1.merge_list) == 2
    cluster_sizes_1 = [c.support_count for c in merge_1.merge_list]
    assert sorted(cluster_sizes_1) == [1, 2]

    # ---------------------------------------------------------
    # Configuration 2: Strict on match_ratio_thresh (min 50% match)
    # SV_SMALL_1 (ratio 0.28) -> REJECTED
    # SV_SMALL_2 (ratio 0.57) -> ACCEPTED
    # ---------------------------------------------------------
    merge_2 = ConsensusCallset([deepcopy(sv_large), deepcopy(sv_small_1), deepcopy(sv_small_2)],
                                   thresh=10, match_ratio_thresh=0.50, ordered_sources=['Manta', 'Severus', 'SVision'])
    merge_2.merge()

    # Expect 2 final SVs: {SV_LARGE + SV_S2} (size 2) and {SV_S1} (size 1)
    assert len(merge_2.merge_list) == 2
    cluster_sizes_2 = [c.support_count for c in merge_2.merge_list]
    assert sorted(cluster_sizes_2) == [1, 2]

    # ---------------------------------------------------------
    # Configuration 3: Loose thresholds (accept both fragments)
    # Max diff 5, min ratio 0.20
    # ---------------------------------------------------------
    merge_3 = ConsensusCallset([deepcopy(sv_large), deepcopy(sv_small_1), deepcopy(sv_small_2)],
                                   thresh=10,
                                   unmatched_bnds_thresh=5, match_ratio_thresh=0.20, ordered_sources=['Manta', 'Severus', 'SVision'])
    merge_3.merge()

    # We expect 1 big cluster containing all 3 SVs
    assert len(merge_3.merge_list) == 1

    merge_3_sv = merge_3.merge_list[0]
    assert merge_3_sv.support_count == 3
    # Ensure all 7 breakends survived the median joining process
    assert len(merge_3_sv.representative_sv.breakends) == 7
