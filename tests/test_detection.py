import pytest
from copy import deepcopy

from verix.bench import BenchCallset
from verix.sv import Callset, Breakend, SV
from test_parsers import check_performance_self


@pytest.fixture
def long_sv():
    breakends = [Breakend(chrom='chr1', pos=pos * 500, parent_id=0) for pos in range(7)]
    gt_sv = SV(parent_id=0, parent_type='ABCDEF', breakends=breakends, bp_merge_threshold=1)

    return gt_sv


def detection(call_svs, gt_svs, match_threshold=500, pick='single'):
    callset = Callset(call_svs)
    gt_callset = BenchCallset(gt_svs, callset, match_threshold)

    gt_callset.collect_candidate_detections()
    gt_callset.select_best_matches(pick=pick)
    return gt_callset


def change_id(call_sv, new_id):
    call_sv.parent_id = new_id
    for bp in call_sv.breakends:
        bp.parent_id = new_id


def check_optimal_sv(call, detection_mode, distance, diff_bnd, pids, contig, fragmented, is_gt=False):
    optimal = call.optimal_candidate
    other_sv = optimal.call_sv if is_gt else optimal.gt_sv
    print(optimal)
    assert ((optimal.mode == detection_mode) and (optimal.contiguity == contig) and
            (optimal.diff_bnd == diff_bnd and (optimal.total_dist == distance) and
             (other_sv.parent_id in pids) and (optimal.fragmented == fragmented))), 'Wrong optimal sv'


def move_bps_limits(call_sv):
    add = True
    for bp in call_sv.breakends:
        if add:
            bp.pos += 249
        else:
            bp.pos -= 249
        add = not add


def test_complete(long_sv):
    gt_sv = long_sv
    call_sv = deepcopy(gt_sv)
    num_bp = len(gt_sv.breakends)

    move_bps_limits(call_sv)

    gt_callset = detection([call_sv], [gt_sv])
    check_performance_self(gt_callset.compute_metrics(gt_callset.pred_callset.callset2dataframe(), gt_callset.callset2dataframe()))

    assert gt_sv.detection_mode == 'complete', 'Ground truth not matched'
    assert call_sv.detection_mode == 'complete', 'Callset not matched'
    check_optimal_sv(gt_sv, 'complete', 249 * num_bp, 0, [call_sv.parent_id], 'Contiguous', False, is_gt=True)
    assert gt_sv.match_score[0][1] == call_sv.match_score[0][1] == 249 * num_bp, 'Wrong distance computation'


def check_incomplete(gt_sv, call_sv, detection_mode, num_matches, contiguity, fragmented, dist=None, match_threshold=500):
    gt_callset = detection([call_sv], [gt_sv], match_threshold=match_threshold)
    performance = gt_callset.compute_metrics(gt_callset.pred_callset.callset2dataframe(), gt_callset.callset2dataframe())
    assert performance['f1'] == performance['recall'] == performance['precision'] == 0.0, 'Wrong metrics'
    assert performance['Valid Calls'] == performance['Truth Valid SVs'] == 1, 'Missing Calls'

    diff_bp = len(gt_sv.breakends) + len(call_sv.breakends) - 2 * num_matches

    assert gt_sv.detection_mode != 'complete', 'Ground truth not matched'
    assert call_sv.detection_mode != 'complete', 'Callset not matched'
    check_optimal_sv(gt_sv, detection_mode, dist if dist is not None else 249 * num_matches, diff_bp, [call_sv.parent_id], contiguity, fragmented, is_gt=True)

    assert len(gt_callset.incomplete_detections[gt_sv.parent_id]) == 1, 'Wrong incomplete detections'
    assert not gt_callset.complete_detections[gt_sv.parent_id], 'Wrong complete detections'


def test_contiguous_subset(long_sv):
    gt_sv = long_sv
    call_sv = deepcopy(gt_sv)

    call_sv.breakends = call_sv.breakends[2:][:-2]
    move_bps_limits(call_sv)
    check_incomplete(gt_sv, call_sv, 'partial', len(call_sv.breakends), 'Contiguous', False)


def test_discontiguous_subset(long_sv):
    gt_sv = long_sv
    call_sv = deepcopy(gt_sv)
    call_sv.breakends = [bp for idx_bp, bp in enumerate(call_sv.breakends) if idx_bp % 2 == 0]
    move_bps_limits(call_sv)
    check_incomplete(gt_sv, call_sv, 'partial', len(call_sv.breakends), 'Discontiguous', False)
    assert call_sv.spurious == 0
    assert gt_sv.spurious == 3


def test_spurious_superset(long_sv):
    gt_sv = long_sv
    call_sv = deepcopy(gt_sv)
    move_bps_limits(call_sv)
    call_sv.breakends.append(Breakend(chrom='chr1', pos=25000, parent_id=0))
    check_incomplete(gt_sv, call_sv, 'partial', len(call_sv.breakends) - 1, 'Contiguous', False)
    assert call_sv.spurious == 1
    assert gt_sv.spurious == 0


def test_discontiguous_superset(long_sv):
    gt_sv = long_sv
    call_sv = deepcopy(gt_sv)
    call_sv.breakends.append(Breakend(chrom='chr1', pos=250, parent_id=0))
    check_incomplete(gt_sv, call_sv, 'partial', len(call_sv.breakends) - 1, 'Discontiguous', False, dist=0)
    assert call_sv.spurious == 1
    assert gt_sv.spurious == 0


def test_partial(long_sv):
    gt_sv = long_sv
    call_sv = deepcopy(gt_sv)
    call_sv.breakends = call_sv.breakends[2:]
    move_bps_limits(call_sv)
    call_sv.breakends.append(Breakend(chrom='chr1', pos=25000, parent_id=0))
    check_incomplete(gt_sv, call_sv, 'partial', len(call_sv.breakends) - 1, 'Contiguous', False)
    assert call_sv.spurious == 1
    assert gt_sv.spurious == 2


def test_discontiguous_partial(long_sv):
    gt_sv = long_sv
    call_sv = deepcopy(gt_sv)
    call_sv.breakends = [bp for idx_bp, bp in enumerate(call_sv.breakends) if idx_bp % 2 == 1]
    move_bps_limits(call_sv)
    call_sv.breakends.append(Breakend(chrom='chr1', pos=25000, parent_id=0))
    check_incomplete(gt_sv, call_sv, 'partial', len(call_sv.breakends) - 1, 'Discontiguous', False)
    assert call_sv.spurious == 1
    assert gt_sv.spurious == 4


def test_complete_incomplete(long_sv):
    gt_sv = long_sv
    call_sv = deepcopy(gt_sv)
    num_bp = len(gt_sv.breakends)
    move_bps_limits(call_sv)

    incomplete_call_sv = deepcopy(gt_sv)
    change_id(incomplete_call_sv, 1)
    incomplete_call_sv.breakends.append(Breakend(chrom='chr1', pos=25000, parent_id=1))

    gt_callset = detection([incomplete_call_sv, call_sv], [gt_sv])
    performance = gt_callset.compute_metrics(gt_callset.pred_callset.callset2dataframe(), gt_callset.callset2dataframe())

    assert performance['Valid Calls'] == 2 and performance['Truth Valid SVs'] == 1, 'Missing Calls'
    assert performance['TP-pred'] == performance['TP-truth'] == 1, 'Missing TP'
    assert performance['FP'] == 1, 'Missing FP'
    assert performance['f1'] == 1 / 1.5 and performance['recall'] == 1.0 and performance['precision'] == 0.5, 'Wrong metrics'

    assert gt_sv.detection_mode == 'complete', 'Ground truth not matched'
    assert call_sv.detection_mode == 'complete', 'Callset not matched'
    assert incomplete_call_sv.detection_mode != 'complete', 'Incomplete detection missing'
    assert gt_sv.match_score[0][1] == call_sv.match_score[0][1] == 249 * num_bp, 'Wrong distance computation'

    diff_bp = len(gt_sv.breakends) + len(incomplete_call_sv.breakends) - 2 * num_bp
    check_optimal_sv(gt_sv, 'complete', 249 * num_bp, 0, [call_sv.parent_id], 'Contiguous', False, is_gt=True)
    check_optimal_sv(incomplete_call_sv, 'partial', 0, diff_bp, [gt_sv.parent_id], 'Contiguous', False)
    assert incomplete_call_sv.spurious == 1
    assert gt_sv.spurious == 0
    assert call_sv.spurious == 0


def test_multiple_incomplete(long_sv):
    gt_sv = long_sv
    num_bp = len(gt_sv.breakends)

    superset_call_sv = deepcopy(gt_sv)
    move_bps_limits(superset_call_sv)
    change_id(superset_call_sv, 1)
    superset_call_sv.breakends.append(Breakend(chrom='chr1', pos=25000, parent_id=1))

    subset_call_sv = deepcopy(gt_sv)
    change_id(subset_call_sv, 2)
    subset_call_sv.breakends = subset_call_sv.breakends[2:]

    partial_call_sv = deepcopy(gt_sv)
    change_id(partial_call_sv, 3)
    partial_call_sv.breakends = partial_call_sv.breakends[2:]
    partial_call_sv.breakends.append(Breakend(chrom='chr1', pos=25000, parent_id=3))

    no_detection = SV(breakends=[Breakend(chrom='chr1', pos=25000, parent_id=4),
                                   Breakend(chrom='chr1', pos=30000, parent_id=4)],
                      parent_id=4, parent_type='Test')

    discontiguous_aggregate = deepcopy(gt_sv)
    change_id(discontiguous_aggregate, 5)
    discontiguous_aggregate.breakends.append(Breakend(chrom='chr1', pos=250, parent_id=5))

    incomplete_calls = [superset_call_sv, subset_call_sv, partial_call_sv, discontiguous_aggregate]
    gt_callset = detection(incomplete_calls + [no_detection], [gt_sv], match_threshold=249)
    performance = gt_callset.compute_metrics(gt_callset.pred_callset.callset2dataframe(), gt_callset.callset2dataframe())

    assert performance['Valid Calls'] == 5 and performance['Truth Valid SVs'] == 1, 'Missing Calls'
    assert performance['TP-pred'] == performance['TP-truth'] == 0, 'Missing TP'
    assert performance['FP'] == 5, 'Missing FP'
    assert performance['f1'] == performance['recall'] == performance['precision'] == 0., 'Wrong metrics'

    assert gt_sv.detection_mode != 'complete', 'Ground truth matched or no detection'
    assert all(call.detection_mode != 'complete' for call in incomplete_calls), 'Incomplete detection missing'

    assert not gt_sv.match_score, 'Wrong distance computation'
    diff_bp = len(gt_sv.breakends) + len(discontiguous_aggregate.breakends) - 2 * num_bp
    check_optimal_sv(gt_sv, 'partial', 0, diff_bp, [discontiguous_aggregate.parent_id], 'Discontiguous', False, is_gt=True)
    check_optimal_sv(superset_call_sv, 'partial', 249*num_bp, diff_bp, [gt_sv.parent_id], 'Contiguous', False)
    diff_bp = len(gt_sv.breakends) + len(discontiguous_aggregate.breakends) - 2 * num_bp
    check_optimal_sv(discontiguous_aggregate, 'partial', 0, diff_bp, [gt_sv.parent_id], 'Discontiguous', False)
    diff_bp = len(gt_sv.breakends) - len(subset_call_sv.breakends)
    check_optimal_sv(subset_call_sv, 'partial', 0, diff_bp, [gt_sv.parent_id], 'Contiguous', True)
    diff_bp = len(gt_sv.breakends) - len(partial_call_sv.breakends) + 2
    check_optimal_sv(partial_call_sv, 'partial', 0, diff_bp, [gt_sv.parent_id], 'Contiguous', True)

    assert no_detection.detection_mode == 'miss', 'Wrong detection'
    assert no_detection.optimal_candidate is None, 'Wrong optimal candidate'


def test_close_gt_complete(long_sv):
    gt_sv = long_sv
    num_bp = len(gt_sv.breakends)

    gt2 = SV(breakends=[Breakend(pos=gt_sv.breakends[-1].pos, parent_id=1, chrom='chr1'),
                        Breakend(pos=gt_sv.breakends[-1].pos+50, parent_id=1, chrom='chr1')],
                        parent_id=1, parent_type='A')

    call_complete = deepcopy(gt_sv)
    move_bps_limits(call_complete)
    call_complete_2 = deepcopy(gt2)
    call_complete_2.breakends[0].pos -= 490
    call_complete_2.breakends[1].pos -= 50

    gt_callset = detection([call_complete, call_complete_2], [gt_sv, gt2])
    print('gt callset, imcomplete', gt_callset.incomplete_detections)
    performance = gt_callset.compute_metrics(gt_callset.pred_callset.callset2dataframe(), gt_callset.callset2dataframe())

    assert performance['Valid Calls'] == 2 and performance['Truth Valid SVs'] == 2, 'Missing Calls'
    assert performance['TP-pred'] == performance['TP-truth'] == 2, 'Missing TP'
    assert performance['FP'] == 0, 'Missing FP'
    assert performance['f1'] == performance['recall'] == performance['precision'] == 1., 'Wrong metrics'

    check_optimal_sv(gt_sv, 'complete', 249*num_bp, 0, [call_complete.parent_id], 'Contiguous', False, is_gt=True)
    check_optimal_sv(gt2, 'complete', 540, 0, [call_complete_2.parent_id], 'Contiguous', False, is_gt=True)
    print('gt callset, complete', gt_callset.complete_detections)
    assert len(gt_callset.complete_detections[gt_sv.parent_id]) == len(gt_callset.complete_detections[gt2.parent_id]) == 1

def test_close_gt_incomplete(long_sv):
    gt_sv = long_sv
    num_bp = len(gt_sv.breakends)

    gt2 = SV(breakends=[Breakend(pos=gt_sv.breakends[-1].pos, parent_id=1, chrom='chr1'),
                          Breakend(pos=gt_sv.breakends[-1].pos+50, parent_id=1, chrom='chr1')],
                          parent_id=1, parent_type='A')
    print('gt', gt_sv.breakends)
    print('gt2', gt2.breakends)
    call_complete = deepcopy(gt_sv)
    move_bps_limits(call_complete)
    call_complete_2 = deepcopy(gt2)
    call_complete_2.breakends[0].pos -= 510
    call_complete_2.breakends[1].pos -= 50

    gt_callset = detection([call_complete, call_complete_2], [gt_sv, gt2])
    performance = gt_callset.compute_metrics(gt_callset.pred_callset.callset2dataframe(), gt_callset.callset2dataframe())
    print('DETECTION')
    assert performance['Valid Calls'] == 2 and performance['Truth Valid SVs'] == 2, 'Missing Calls'
    assert performance['TP-pred'] == performance['TP-truth'] == 1, 'Missing TP'
    assert performance['FP'] == 1, 'Missing FP'
    assert performance['f1'] == 0.5 / 1 and performance['recall'] == .5 and performance['precision'] == .5, 'Wrong metrics'

    check_optimal_sv(gt_sv, 'complete', 249*num_bp, 0, [call_complete.parent_id], 'Contiguous', False, is_gt=True)
    diff_bp = len(gt2.breakends) + len(call_complete_2.breakends) - 2 * 1
    check_optimal_sv(gt2, 'aggregate', 0, diff_bp, [call_complete_2.parent_id], 'Contiguous', True, is_gt=True)

    assert len(gt_callset.incomplete_detections[gt_sv.parent_id]) == 1, 'Wrong number of incomplete detections for gt'
    assert len(gt_callset.incomplete_detections[gt2.parent_id]) == 2, 'Missing incomplete detection for gt2'


def test_close_gt_aggregate(long_sv):
    gt_sv = long_sv
    num_bp = len(gt_sv.breakends)
    print('gt', gt_sv.breakends)
    gt2_bps = [Breakend(pos=gt_sv.breakends[-1].pos, parent_id=1, chrom='chr1'),
                          Breakend(pos=gt_sv.breakends[-1].pos+50, parent_id=1, chrom='chr1')]
    gt2 = SV(breakends=gt2_bps, parent_id=1, parent_type='A')
    gt2.update_breakends()
    print('gt2',gt2.breakends)
    call_aggregate = deepcopy(gt_sv)
    call_aggregate.breakends += deepcopy(gt2_bps)
    change_id(call_aggregate, 0)
    call_aggregate.update_breakends()
    print('agg', call_aggregate.breakends)
    gt_callset = detection([call_aggregate], [gt_sv, gt2])
    performance = gt_callset.compute_metrics(gt_callset.pred_callset.callset2dataframe(), gt_callset.callset2dataframe())

    assert performance['Valid Calls'] == 1 and performance['Truth Valid SVs'] == 2, 'Missing Calls'
    assert performance['TP-pred'] == performance['TP-truth'] == 0, 'Missing TP'
    assert performance['FP'] == 1, 'Missing FP'
    assert performance['f1'] == performance['recall'] == performance['precision'] == 0., 'Wrong metrics'
    print(gt_callset.incomplete_detections)
    diff_bp = len(gt_sv.breakends) + len(call_aggregate.breakends) - 2 * num_bp
    check_optimal_sv(gt_sv, 'aggregate', 0, diff_bp, [call_aggregate.parent_id], 'Contiguous', False, is_gt=True)
    diff_bp = len(gt2.breakends) + len(call_aggregate.breakends) - 2 * len(gt2_bps)
    check_optimal_sv(gt2, 'aggregate', 0, diff_bp, [call_aggregate.parent_id], 'Contiguous', False, is_gt=True)

    assert len(gt_callset.incomplete_detections[gt_sv.parent_id]) == 1, 'Wrong number of incomplete detections for gt'
    assert len(gt_callset.incomplete_detections[gt2.parent_id]) == 1, 'Wrong number of  incomplete detection for gt2'


def test_multiple_matches():
    gt_bps = [Breakend(pos=0, parent_id='A', chrom='chr1'), Breakend(pos=100, parent_id='A', chrom='chr1')]
    gt_sv = SV(breakends=gt_bps, parent_id='A', parent_type='A')

    gt_bps_2 = [Breakend(pos=150, parent_id='B', chrom='chr1'), Breakend(pos=300, parent_id='B', chrom='chr1')]
    gt_sv_2 = SV(breakends=gt_bps_2, parent_id='B', parent_type='A')

    call_bps = [Breakend(pos=50, parent_id='test', chrom='chr1'), Breakend(pos=200, parent_id='test', chrom='chr1')]
    call_sv = SV(breakends=call_bps, parent_id='test', parent_type='A')

    call_sv_copy = deepcopy(call_sv)

    gt_callset = detection([call_sv, call_sv_copy], [gt_sv, gt_sv_2], pick='multi')

    check_optimal_sv(gt_sv, 'complete', 150, 0, ['test'], 'Contiguous', False, is_gt=True)
    check_optimal_sv(gt_sv_2, 'complete', 200, 0, ['test'], 'Contiguous', False, is_gt=True)

    performance = gt_callset.compute_metrics(gt_callset.pred_callset.callset2dataframe(), gt_callset.callset2dataframe())
    assert performance['Valid Calls'] == 2 and performance['Truth Valid SVs'] == 2, 'Missing Calls'
    assert performance['TP-pred'] == 2 and performance['TP-truth'] == 1, 'Missing TP'
    assert performance['FP'] == 0, 'Missing FP'
    assert performance['f1'] == 1 / 1.5 and performance['recall'] == .5 and performance['precision'] == 1., 'Wrong metrics'

    assert all([not gt_callset.incomplete_detections[pid] for pid in
                ['A', 'B']]), 'Wrong number of incomplete detections for gt'
    assert all([len(gt_callset.complete_detections[pid]) == 2 for pid in
                ['A', 'B']]), 'Wrong number of complete detections for gt'


def test_multiple_incomplete_interchr(long_sv):
    gt_sv = long_sv
    num_bp = len(gt_sv.breakends)

    gt_sv_2_bps = [
        Breakend(chrom='chr2', pos=30000, parent_id=1),
        Breakend(chrom='chr1', pos=0, parent_id=1), Breakend(chrom='chr1', pos=500, parent_id=1),
        Breakend(chrom='chr1', pos=1000, parent_id=1), Breakend(chrom='chr1', pos=1500, parent_id=1),
        Breakend(chrom='chr2', pos=25000, parent_id=1),
    ]
    gt_sv_2 = SV(breakends=gt_sv_2_bps, parent_id=1, parent_type='AB_C')

    complete_match = deepcopy(gt_sv)
    change_id(complete_match, -1)

    complete_match_2 = deepcopy(gt_sv_2)
    change_id(complete_match_2, 0)

    aggregate_call_sv = deepcopy(gt_sv)
    move_bps_limits(aggregate_call_sv)
    change_id(aggregate_call_sv, 1)
    aggregate_call_sv.breakends.append(Breakend(chrom='chr2', pos=25000, parent_id=1))

    subset_call_sv = deepcopy(gt_sv)
    change_id(subset_call_sv, 2)
    subset_call_sv.breakends = subset_call_sv.breakends[2:]

    partial_call_sv = deepcopy(gt_sv)
    change_id(partial_call_sv, 3)
    partial_call_sv.breakends = partial_call_sv.breakends[2:]
    partial_call_sv.breakends.append(Breakend(chrom='chr1', pos=25000, parent_id=3))

    no_detection = SV(breakends=[Breakend(chrom='chr1', pos=25000, parent_id=4),
                                 Breakend(chrom='chr1', pos=30000, parent_id=4)],
                      parent_id=4, parent_type='Test')

    discontiguous_aggregate = deepcopy(gt_sv)
    change_id(discontiguous_aggregate, 5)
    move_bps_limits(discontiguous_aggregate)
    discontiguous_aggregate.breakends.append(Breakend(chrom='chr1', pos=0, parent_id=5))
    discontiguous_aggregate.breakends.append(Breakend(chrom='chr2', pos=25000, parent_id=5))

    contiguous_aggregate_interchr = SV(breakends=[Breakend(chrom='chr2', pos=25000, parent_id='A')] + deepcopy(gt_sv.breakends),
                                       parent_id='A', parent_type='A_B')
    change_id(contiguous_aggregate_interchr, 'A')

    incomplete_calls = [aggregate_call_sv, subset_call_sv, partial_call_sv, discontiguous_aggregate, contiguous_aggregate_interchr,]
    complete_calls = [complete_match, complete_match_2]
    gt_callset = detection(incomplete_calls + complete_calls + [no_detection], [gt_sv, gt_sv_2], match_threshold=249)
    performance = gt_callset.compute_metrics(gt_callset.pred_callset.callset2dataframe(), gt_callset.callset2dataframe())

    assert performance['Valid Calls'] == 8 and performance['Truth Valid SVs'] == 2, 'Missing Calls'
    assert performance['TP-pred'] == performance['TP-truth'] == 2, 'Missing TP'
    assert performance['FP'] == 6, 'Missing FP'
    assert performance['f1'] == 4 / 10 and performance['recall'] == 1 and performance['precision'] == 2/8, 'Wrong metrics'

    assert gt_sv.detection_mode == 'complete', 'Ground truth matched or no detection'
    assert gt_sv_2.detection_mode == 'complete', 'Ground truth matched or no detection'
    assert all(call.detection_mode != 'complete' for call in incomplete_calls), 'Incomplete detection missing'

    check_optimal_sv(gt_sv, 'complete', 0, 0, [complete_match.parent_id], 'Contiguous', False, is_gt=True)
    check_optimal_sv(gt_sv_2, 'complete', 0, 0, [complete_match_2.parent_id], 'Contiguous', False, is_gt=True)
    diff_bp = len(aggregate_call_sv.breakends) - num_bp
    check_optimal_sv(aggregate_call_sv, 'aggregate', 249*num_bp, diff_bp, [gt_sv.parent_id], 'Contiguous', False)
    diff_bp = len(discontiguous_aggregate.breakends) - num_bp
    check_optimal_sv(discontiguous_aggregate, 'aggregate', 6*249, diff_bp, [gt_sv.parent_id], 'Discontiguous', False)
    diff_bp = len(contiguous_aggregate_interchr.breakends) - num_bp
    check_optimal_sv(contiguous_aggregate_interchr, 'aggregate', 0, diff_bp, [gt_sv.parent_id], 'Contiguous', False)
    diff_bp = num_bp + len(subset_call_sv.breakends) - 2 * len(subset_call_sv.breakends)
    check_optimal_sv(subset_call_sv, 'partial', 0, diff_bp, [gt_sv.parent_id], 'Contiguous', True)
    diff_bp = num_bp + len(partial_call_sv.breakends) - 2 * (len(partial_call_sv.breakends)-1)
    check_optimal_sv(partial_call_sv, 'partial', 0, diff_bp, [gt_sv.parent_id], 'Contiguous', True)

    assert gt_sv.match_score[0] == (0, 0), 'Wrong distance computation'
    assert no_detection.detection_mode == 'miss', 'Wrong detection'
    assert no_detection.optimal_candidate is None, 'Wrong optimal candidate'

    assert gt_sv.match_id[0] == complete_match.parent_id, 'Wrong match interchromosomal'
    assert gt_sv_2.match_id[0] == complete_match_2.parent_id, 'Wrong match interchromosomal'
