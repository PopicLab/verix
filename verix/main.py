import argparse
import logging
import sys
from pathlib import Path

from verix import __version__, __toolname__
from verix.parser import parse
from verix.sv import Callset
from verix.bench import BenchCallset
from verix.plot import Plotter
from verix.merge import ConsensusCallset
from verix.vcf_writer import write_merge_vcf


FORMAT_CHOICES = ['default', 'single_rec', 'multi_rec']


def parse_args():
    parser = argparse.ArgumentParser(description=f"{__toolname__} {__version__}")
    subparsers = parser.add_subparsers(dest='command', required=True)

    shared = argparse.ArgumentParser(add_help=False)
    shared.add_argument('-o', '--output', required=True, help='Output directory.')
    shared.add_argument('-st', '--svtype', default='SVTYPE', help='INFO field used as SV type.')
    shared.add_argument('-mt', '--match_threshold', default=500, type=int, help='Max distance between breakends (500).')
    shared.add_argument('-mrt', '--match_ratio_thresh', default=None, type=float,
                        help='Min ratio of matched breakends.')
    shared.add_argument('-bpt', '--bp_merge_threshold', default=2, type=int, help='Distance to merge breakends (2).')
    shared.add_argument('-s', '--sizemin', default=0.0, type=float, help='Filter SVs smaller than this (0).')
    shared.add_argument('-S', '--sizemax', default=None, type=float, help='Filter SVs larger than this.')
    shared.add_argument('-chr', '--chr_list', default=None, nargs='+', help='Filter SV records not in this list.')
    shared.add_argument('-q', '--qual', default=0.0, type=float, help='Filter SV records of insufficient quality (0).')
    shared.add_argument('--enforce_svtype', action='store_true', help='Whether the SV type should be ignored when finding candidates (False).')
    shared.add_argument('--enforce_genotype', action='store_true',
                        help='Whether the genotype should be ignored when finding candidates (False).')

    # Benchmarking parameters
    bench = subparsers.add_parser('bench', parents=[shared], help='Benchmark VCF files against a truth set.')
    bench.set_defaults(func=benchmark)

    bench.add_argument('-p', '--pred', required=True, help='VCF file to compare.')
    bench.add_argument('-t', '--truth', required=True, help='Truth set VCF file.')
    bench.add_argument('-ft', '--format_truth', required=False, default='default', choices=FORMAT_CHOICES, help='Truth format (default).')
    bench.add_argument('-fp', '--format_pred', required=False, default='default', choices=FORMAT_CHOICES, help='Compare format (default).')
    bench.add_argument('-csvt', '--csv_info_truth', default=None, help='INFO field for base SVID/BKPS.')
    bench.add_argument('-csvp', '--csv_info_pred', default=None, help='INFO field for compare SVID/BKPS.')
    #bench.add_argument('--pick', default='single', choices=['single', 'multi'], help='Match allowance (single).')
    bench.add_argument('-ubt', '--unmatched_bnds_thresh', default=None, type=int, help='Max unmatched breakends allowed (None).')
    bench.add_argument('--plot', action='store_true',
                       help='Whether benchmarking figures should be generated (False).')

    # Merge parameters
    merge = subparsers.add_parser('merge', parents=[shared], help='Merge VCF files into a merge.')
    merge.set_defaults(func=run_merge)

    merge.add_argument('-i', '--inputs', nargs='+', required=True, help='List of VCF files to merge.')
    merge.add_argument('-f', '--formats', nargs='+', required=True, choices=FORMAT_CHOICES,
                           help='Formats for inputs.')
    merge.add_argument('-csv', '--csv_info_list', nargs='+', default=None, help='INFO field list for SVID/BKPS.')
    merge.add_argument('-src', '--sources', nargs='+', help='Ordered list of names for merging.')
    merge.add_argument('-ubt', '--unmatched_bnds_thresh', default=0, type=int, help='Max unmatched breakends allowed (0).')

    return parser.parse_args()


def default_csv_info(csv_info, vcf_format):
    if csv_info is None:
        csv_info = 'SVID' if vcf_format == 'multi_rec' else 'BKPS'
    return csv_info


def benchmark(pred, truth, output, format_truth='bnd', format_pred='bnd', csv_info_truth=None, csv_info_pred=None,
              svtype='SVTYPE', pick='single', qual=None, match_threshold=500, bp_merge_threshold=2, chr_list=None,
              sizemin=0, sizemax=None, unmatched_bnds_thresh=None, match_ratio_thresh=None, enforce_svtype=True, enforce_genotype=True,
              plot=True):
    output_folder = output + '/output/'

    csv_info_truth = default_csv_info(csv_info_truth, format_truth)
    csv_info_pred = default_csv_info(csv_info_pred, format_pred)

    pred_svs = parse(pred, vcf_format=format_pred, csv_info=csv_info_pred, svtype_name=svtype,
                         merge_threshold=bp_merge_threshold, chr_list=chr_list, sizemin=sizemin, sizemax=sizemax, qual=qual)
    assert pred_svs, 'No SV found in pred'
    logging.info(f"Loaded a total of {len(pred_svs)} call SVs.")

    pred_callset = Callset(pred_svs)

    gt_svs = parse(truth, vcf_format=format_truth, csv_info=csv_info_truth, svtype_name=svtype, merge_threshold=bp_merge_threshold, chr_list=chr_list,
                   sizemin=sizemin, sizemax=sizemax, qual=qual)
    assert gt_svs, 'No SV found in truth'
    logging.info(f"Loaded a total of {len(gt_svs)} truth SVs.")

    gt_callset = BenchCallset(gt_svs, pred_callset=pred_callset, thresh=match_threshold, unmatched_bnds_thresh=unmatched_bnds_thresh,
                              match_ratio_thresh=match_ratio_thresh, enforce_svtype=enforce_svtype, enforce_genotype=enforce_genotype)

    logging.info("Matching SVs")
    gt_callset.collect_candidate_detections()
    gt_callset.select_best_matches(pick=pick)

    # Write VCF output
    logging.info("Generating outputs")
    results_dict = gt_callset.create_outputs(output_folder)

    if plot:
        Plotter(results_dict, output_folder).plot()


def run_merge(inputs, output, formats, csv_info_list=None, svtype='SVTYPE', qual=0, sources=None,
              match_threshold=500, bp_merge_threshold=1, chr_list=None, sizemin=0, sizemax=None,
              unmatched_bnds_thresh=None, match_ratio_thresh=None, enforce_svtype=True, enforce_genotype=True,):
    output_folder = Path(output) / 'output'
    output_folder.mkdir(parents=True, exist_ok=True)

    all_svs = []
    ordered_sources = [] if sources is None else sources

    logging.info(f"Parsing {len(inputs)} VCF files for merge building...")

    for idx_vcf, vcf_path in enumerate(inputs):
        # Use the filename (without extension) as the source identifier

        if sources is None:
            source_name = Path(vcf_path).stem + '_' + str(idx_vcf)
            ordered_sources.append(source_name)
        else:
            source_name = sources[idx_vcf]

        assert idx_vcf < len(formats), f'Please provide the format for each file provided {len(formats)} for {len(inputs)} files.'
        vcf_format = formats[idx_vcf]
        csv_info = None
        if csv_info_list is not None:
            assert idx_vcf < len(
                csv_info_list), ('When provided csv_info_list must define the csv info name for each file incluing the ones in format bnd. '
                                 f'Provided {len(csv_info_list)} for {len(inputs)} files.')
            csv_info = csv_info_list[idx_vcf]
        csv_info = default_csv_info(csv_info, vcf_format)
        svs = parse(vcf_path, vcf_format=vcf_format, csv_info=csv_info, svtype_name=svtype, merge_threshold=bp_merge_threshold,
                    chr_list=chr_list, sizemin=sizemin, sizemax=sizemax, qual=qual, source=source_name)
        if not svs:
            logging.warning(f'No SV found in {vcf_path}')
        all_svs.extend(svs)
    assert all_svs, 'No SV found in all VCF files'
    logging.info(f"Loaded a total of {len(all_svs)} SVs. Building merge...")

    merge_callset = ConsensusCallset(
        svs=all_svs,
        thresh=match_threshold,
        unmatched_bnds_thresh=unmatched_bnds_thresh,
        match_ratio_thresh=match_ratio_thresh,
        ordered_sources=ordered_sources,
        enforce_svtype=enforce_svtype,
        enforce_genotype=enforce_genotype
    )

    merge_callset.merge()

    logging.info(f"Generated {len(merge_callset.merge_list)} merge SV clusters.")

    write_merge_vcf(merge_callset, ordered_sources, output_folder)
    logging.info(f"Consensus outputs written to {output_folder}/merge.vcf")


def run_main():
    args = parse_args()
    out_path = Path(args.output)
    log_dir = out_path / "logs"
    out_dir = out_path / "output"

    # setup the experiment directory structure
    Path(log_dir).mkdir(parents=True, exist_ok=True)
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    # logging
    logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s',
                        handlers=[logging.FileHandler(log_dir / 'main.log', mode='w'),
                                  logging.StreamHandler(sys.stdout)])
    
    # Prevents plotly to fill the terminal when generating sankey plots
    logging.getLogger("kaleido").setLevel(logging.WARNING)
    logging.getLogger("choreographer").setLevel(logging.WARNING)

    logging.info(f"{__toolname__} version {__version__}")
    logging.info(args)

    args_dict = vars(args)
    command = args_dict.pop('command')
    func = args_dict.pop('func')
    func(**args_dict)


if __name__ == '__main__':
    run_main()
