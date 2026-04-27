from verix.parser import parse
from verix.main import default_csv_info
from verix.sv import Callset
from verix.bench import BenchCallset
import os

FILE_PATH = os.path.dirname(os.path.abspath(__file__))


def run_benchmark(compare, base, format='standard', baseformat='standard', typeignore=True, pick='single',
                  match_threshold=500, bp_merge_treshold=2, chr_list=None, sizemin=0, sizemax=None):

    pred_svs = parse(compare, vcf_format=format, csv_info=default_csv_info(None, format),
                         merge_threshold=bp_merge_treshold,
                         chr_list=chr_list,
                         sizemin=sizemin, sizemax=sizemax)
    gt_svs = parse(base,vcf_format=baseformat, csv_info=default_csv_info(None, baseformat),
                   merge_threshold=bp_merge_treshold, chr_list=chr_list,
                   sizemin=sizemin, sizemax=sizemax)

    pred_callset = Callset(pred_svs)
    gt_callset = BenchCallset(gt_svs, pred_callset, match_threshold)

    gt_callset.collect_candidate_detections()
    gt_callset.select_best_matches(pick=pick)

    return gt_callset.compute_metrics(pred_callset.callset2dataframe(), gt_callset.callset2dataframe())


def check_performance_self(performance):
    assert performance['Valid Calls'] > 0
    assert performance['Valid Calls'] == performance['Truth Valid SVs'], 'Missing Calls'
    print(performance)
    assert performance['TP-pred'] == performance['TP-truth'] == performance['Valid Calls'], 'Missing TP'
    assert performance['FP'] == performance['FN'] == 0, 'Unmatched calls'
    assert performance['f1'] == performance['recall'] == performance['precision'] == 1.0, 'Wrong metrics'


def test_insilico2insilico():
    insilicosv_vcf_path = FILE_PATH + '/test_vcf/insilicosv.vcf'

    performance = run_benchmark(insilicosv_vcf_path, insilicosv_vcf_path, format='multi_rec', baseformat='multi_rec')
    check_performance_self(performance)


def test_standard2standard():
    insilicosv_vcf_path = FILE_PATH + '/test_vcf/insilicosv.vcf'

    performance = run_benchmark(insilicosv_vcf_path, insilicosv_vcf_path)
    check_performance_self(performance)


def test_svim2svim():
    vcf_path = FILE_PATH + '/test_vcf/svim.vcf.gz'

    performance = run_benchmark(vcf_path, vcf_path, format='bnd', baseformat='bnd', match_threshold=0)
    check_performance_self(performance)


def test_svision2svision():
    vcf_path = FILE_PATH + '/test_vcf/svision.vcf'

    performance = run_benchmark(vcf_path, vcf_path, format='single_rec', baseformat='single_rec', match_threshold=0)
    check_performance_self(performance)


def test_severus2severus():
    vcf_path = FILE_PATH + '/test_vcf/severus.vcf'

    performance = run_benchmark(vcf_path, vcf_path, format='multi_rec', baseformat='multi_rec', match_threshold=0)
    check_performance_self(performance)


def test_sniffles2sniffles():
    vcf_path = FILE_PATH + '/test_vcf/sniffles.vcf'

    performance = run_benchmark(vcf_path, vcf_path, format='bnd', baseformat='bnd', match_threshold=0)
    check_performance_self(performance)
