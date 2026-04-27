import pytest
import pandas as pd
import json

from verix.sv import Callset
from verix.bench import BenchCallset
from verix.parser import parse

@pytest.fixture
def run_benchmarking(tmp_path):
    """
    Setup function to parse VCFs and run the BenchCallset logic.
    """
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    base_path = "tests/test_vcf/base.vcf"
    comp_path = "tests/test_vcf/comp.vcf"

    base_svs = parse(base_path, csv_info='EVENT', vcf_format='multi_rec')
    comp_svs = parse(comp_path, csv_info='EVENT', vcf_format='multi_rec')

    comp_callset = Callset(comp_svs)
    gt_callset = BenchCallset(base_svs, comp_callset, thresh=50)

    gt_callset.collect_candidate_detections()
    gt_callset.select_best_matches()
    gt_callset.create_outputs(str(output_dir) + "/")
    return output_dir


def test_detection_modes(run_benchmarking):
    """Verify that each SV was assigned the correct detection mode in truth.csv"""
    df_base = pd.read_csv(run_benchmarking / "truth.csv")

    # GT_COMP should be complete
    assert df_base.loc[df_base['parent_id'] == 'GT_COMP', 'detection_mode'].iloc[0] == 'complete'

    # GT_PART should be partial (Call only has 2 of the 3 BNDs)
    assert df_base.loc[df_base['parent_id'] == 'GT_PART', 'detection_mode'].iloc[0] == 'partial'

    # GT_AGGR should be aggregate (Matches multiple call SVs)
    assert df_base.loc[df_base['parent_id'] == 'GT_AGGR', 'detection_mode'].iloc[0] == 'aggregate'

    # GT_MISS should be miss
    assert df_base.loc[df_base['parent_id'] == 'GT_MISS', 'detection_mode'].iloc[0] == 'miss'


def test_summary_json_counts(run_benchmarking):
    """Verify the aggregate stats in summary.json"""
    with open(run_benchmarking / "summary.json") as f:
        summary = json.load(f)

    # We expect 1 complete, 1 partial, 1 aggregate, 1 miss in 'base'
    assert summary['complete-truth'] == 1
    assert summary['partial-truth'] == 1
    assert summary['aggregate-truth'] == 1
    assert summary['miss-truth'] == 1

    # Check Precision/Recall logic
    # TP (for precision) = 1 (C_COMP)
    # TP (for recall) = 1 (GT_COMP)
    # Note: Your logic currently counts only 'complete' for TPs in the select_best_matches
    assert summary['TP-truth'] == 1


def test_optimal_mapping(run_benchmarking):
    """Verify the optimal.csv links IDs correctly"""
    df_optimal = pd.read_csv(run_benchmarking / "optimal.csv")

    # GT_PART should have C_PART as its optimal match
    part_match = df_optimal[df_optimal['GT ID'] == 'GT_PART']
    assert part_match['Call ID'].iloc[0] == 'C_PART'
    assert part_match['detection mode'].iloc[0] == 'partial'

    # C_SPUR should have NA as GT ID
    spur_match = df_optimal[df_optimal['Call ID'] == 'C_SPUR']
    assert spur_match['detection mode'].iloc[0] == 'miss'