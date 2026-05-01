import argparse
import logging
from pathlib import Path
import sys

from verix import __version__, __toolname__
from verix.io import parse_vcf, VCFFormat, write_csv_vcf
from verix.sv import BreakpointAligner
from verix.bench import BenchmarkEngine
from verix.merge import MergeEngine


def parse_args():
    parser = argparse.ArgumentParser(description=f"{__toolname__} {__version__}",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(dest='command', required=True)
    shared = argparse.ArgumentParser(add_help=False)
    shared.add_argument('-o', '--output_dir', metavar='', required=True, help='Output directory')
    shared.add_argument('-d', '--match_thr', metavar='', default=500, type=int, help='Max distance between matching breakpoints')
    shared.add_argument('-s', '--sizemin', metavar='', default=0, type=int, help='Minimum SV interval size')
    shared.add_argument('-S', '--sizemax', metavar='', default=None, type=int, help='Maximum SV interval size')
    shared.add_argument('-b', '--merge_thr', metavar='', default=1, type=int,
                        help='Collapse breakends in a CSV within this distance into a single breakpoint')
    shared.add_argument('--enforce_type', action='store_true', help='Require SV types to match')
    shared.add_argument('--enforce_genotype', action='store_true', help='Require SV genotypes to match')

    # Benchmarking parameters
    bench = subparsers.add_parser('bench', parents=[shared], help='Compare two VCF files (query and target/truthset)',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    bench.set_defaults(func=benchmark)
    bench.add_argument('-p', '--query', metavar='', required=True, help='VCF file with query CSVs')
    bench.add_argument('-t', '--target', metavar='', required=True, help='VCF file with target CSVs')
    bench.add_argument('-fq', '--format_q', required=False, default='default', choices=[e.value for e in VCFFormat], help='VCF format for query')
    bench.add_argument('-ft', '--format_t', required=False, default='default', choices=[e.value for e in VCFFormat], help='VCF format for target')
    bench.add_argument('-lq', '--csv_link_q', metavar='', default=None, help='INFO field for CSV linking in query VCF')
    bench.add_argument('-lt', '--csv_link_t', metavar='', default=None, help='INFO field for CSV linking in target VCF')
    bench.add_argument('--plot', action='store_true', help='Generate benchmarking figures')

    # Consensus parameters
    merge = subparsers.add_parser('consensus', parents=[shared], help='Merge a single or multiple VCF files into a single consensus VCF',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    merge.set_defaults(func=consensus)
    merge.add_argument('-i', '--inputs', nargs='+', metavar='VCF', required=True, help='List of VCF files to merge')
    merge.add_argument('-n', '--names', nargs='+', metavar='NAME', default=[], help='Ordered list of names for each VCF')
    merge.add_argument('-f', '--formats', nargs='+', required=True, choices=[e.value for e in VCFFormat], help='Format type for each VCF')
    merge.add_argument('-csv', '--csv_links', metavar='LINK', nargs='+', default=None, help='INFO field for CSV linking in each VCF')
    merge.add_argument('-ubt', '--unmatched_thr', metavar='', default=0, type=int, help='Max number of unmatched breakpoints allowed for a merge')
    return parser.parse_args()


def benchmark(args):
    if args.csv_link_q is None and args.format_q != VCFFormat.DEFAULT:
        raise ValueError(f"--csv_link_q has to be specified if --format_q is {args.format_q}")
    if args.csv_link_t is None and args.format_t != VCFFormat.DEFAULT:
        raise ValueError(f"--csv_link_t has to be specified if --format_t is {args.format_t}")
    query_svs = parse_vcf(args.query, args.format_q, args.csv_link_q, args.sizemin, args.sizemax, args.merge_thr)
    target_svs = parse_vcf(args.target, args.format_t, args.csv_link_t, args.sizemin, args.sizemax, args.merge_thr)
    logging.info(f"Loaded {len(query_svs)} query SVs")
    logging.info(f"Loaded {len(target_svs)} target/truthset SVs")
    engine = BenchmarkEngine(query_svs, target_svs,
                             BreakpointAligner(args.match_thr, args.enforce_type, args.enforce_genotype))
    engine.find_matches()
    engine.write_stats(Path(args.output_dir) / "report.json")
    engine.write_vcf(Path(args.output_dir) / "matches.vcf")
    write_csv_vcf(query_svs, Path(args.output_dir) / "query.vcf")
    write_csv_vcf(target_svs, Path(args.output_dir) / "target.vcf")
    #if args.plot:
    #    engine.make_plots(args.output_dir / "report.pdf")) #Plotter(args.out_dir).plot()


def consensus(args):
    if len(args.formats) != len(args.inputs):
        raise ValueError(f"--formats has {len(args.formats)} entries for {len(args.inputs)} inputs")
    if args.csv_links is not None and len(args.csv_links) != len(args.inputs):
        raise ValueError(f"--csv_links has {len(args.csv_links)} entries for {len(args.inputs)} inputs")
    if args.names is not None and (len(args.names) != len(set(args.names)) or len(args.names) != len(args.inputs)):
        raise ValueError(f"--names cannot have duplicate entries and must match the length of VCF inputs")
    logging.info(f"Parsing {len(args.inputs)} VCF files for consensus generation...")
    names = args.names or range(len(args.inputs))

    svs = []
    for i, vcf_path in enumerate(args.inputs):
        svs.append(parse_vcf(vcf_path, vcf_format=args.formats[i], csv_link_name=args.csv_links[i] if args.csv_links else None,
                        merge_threshold=args.merge_thr, sizemin=args.sizemin, sizemax=args.sizemax, name=names[i]))
    engine = MergeEngine(svs, BreakpointAligner(args.match_thr, args.enforce_type, args.enforce_genotype))
    engine.find_sv_clusters()
    engine.write_stats(Path(args.output_dir) / "report.json")
    engine.write_vcf(Path(args.output_dir) / "merged.vcf")

def main():
    args = parse_args()
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s',
                        handlers=[logging.FileHandler(args.output_dir + '/main.log', mode='w'),
                                  logging.StreamHandler(sys.stdout)])
    logging.info(f"{__toolname__} version {__version__}")
    logging.info("\n ***Params***\n" + "\n".join(f" {k}: {v}" for k, v in vars(args).items() if not callable(v)))
    args.func(args)


if __name__ == '__main__':
    main()
