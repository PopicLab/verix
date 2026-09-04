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
    shared.add_argument('-b', '--merge_thr', metavar='', default=2, type=int,
                        help='Collapse breakends in a CSV within this distance into a single breakpoint')
    shared.add_argument('--passonly', action='store_true', help='Only keep records whose FILTER is PASS or unset')
    shared.add_argument('--enforce_type', action='store_true', help='Require SV types to match')
    shared.add_argument('--enforce_genotype', action='store_true', help='Require SV genotypes to match')
    shared.add_argument('-f', '--formats', nargs='+', default=[], choices=[e.value for e in VCFFormat],
                        help='Format type for each VCF (expected order for bench: query, target)')
    shared.add_argument('-l', '--csv_links', metavar='LINK', nargs='+', default=[],
                       help='INFO field for CSV linking in each VCF (expected order for bench: query, target)')
    shared.add_argument('-svt', '--types', nargs='+', default=[], help='INFO field for SV type extraction (default SVTYPE)')

    # Benchmarking parameters
    bench = subparsers.add_parser('bench', parents=[shared], help='Compare two VCF files (query and target/truthset)',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    bench.set_defaults(func=benchmark)
    bench.add_argument('-q', '--query', metavar='', required=True, help='VCF file with query CSVs')
    bench.add_argument('-t', '--target', metavar='', required=True, help='VCF file with target CSVs')
    bench.add_argument('--plot', action='store_true', help='Generate benchmarking figures')

    # Consensus parameters
    merge = subparsers.add_parser('consensus', parents=[shared], help='Merge a single or multiple VCF files into a single consensus VCF',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    merge.set_defaults(func=consensus)
    merge.add_argument('-i', '--inputs', nargs='+', metavar='VCF', required=True, help='List of VCF files to merge')
    merge.add_argument('-n', '--names', nargs='+', metavar='NAME', default=[], help='Ordered list of names for each VCF')
    return parser.parse_args()


def benchmark(args):
    for param_name, param in [("formats", args.formats), ("types", args.types)]:
        if param and len(param) != 2:
            raise ValueError(f"--{param_name} must have exactly 2 entries: query, target")
    if args.formats and len(args.csv_links) != 2:
        raise ValueError("--csv_links must have 2 entries when --formats is set")
    if not args.formats and args.csv_links:
        raise ValueError("--csv_links requires --formats to be set")
    fq, ft = args.formats if args.formats else (VCFFormat.DEFAULT, VCFFormat.DEFAULT)
    lq, lt = args.csv_links if args.csv_links else (None, None)
    tq, tt = args.types if args.types else ("SVTYPE", "SVTYPE")
    query_svs = parse_vcf(args.query, fq, lq, tq, args.sizemin, args.sizemax, args.merge_thr, passonly=args.passonly)
    target_svs = parse_vcf(args.target, ft, lt, tt, args.sizemin, args.sizemax, args.merge_thr, passonly=args.passonly)
    logging.info(f"Loaded {len(query_svs)} query SVs")
    logging.info(f"Loaded {len(target_svs)} target/truthset SVs")
    engine = BenchmarkEngine(query_svs,
                             target_svs,
                             BreakpointAligner(args.match_thr, args.enforce_type, args.enforce_genotype))
    engine.find_matches()
    engine.write_stats(Path(args.output_dir) / "report.json")
    engine.write_vcf(Path(args.output_dir) / "matches.vcf")
    write_csv_vcf(query_svs, Path(args.output_dir) / "query.vcf")
    write_csv_vcf(target_svs, Path(args.output_dir) / "target.vcf")
    if args.plot:
        (Path(args.output_dir) / "plots").mkdir(parents=True, exist_ok=True)
        engine.generate_plots(Path(args.output_dir) / "plots")


def consensus(args):
    n_inputs = len(args.inputs)
    for param_name, param in [("formats", args.formats), ("types", args.types),
                              ("csv_links", args.csv_links), ("names", args.names)]:
        if param and len(param) != n_inputs:
            raise ValueError(f"--{param_name} must have exactly {n_inputs} entries")
        if param_name == "names" and len(param) != len(set(param)):
            raise ValueError(f"--names cannot have duplicate entries")
    if not args.formats and args.csv_links:
        raise ValueError("--csv_links requires --formats to be set")
    logging.info(f"Parsing {len(args.inputs)} VCF files for consensus generation...")
    names = args.names or range(len(args.inputs))
    svs = []
    for i, vcf_path in enumerate(args.inputs):
        svs.append(parse_vcf(vcf_path,
                             vcf_format=args.formats[i] if args.formats else VCFFormat.DEFAULT,
                             csv_link_name=args.csv_links[i] if args.csv_links else None,
                             type_name=args.types[i] if args.types else "SVTYPE",
                             merge_threshold=args.merge_thr, sizemin=args.sizemin, sizemax=args.sizemax,
                             name=names[i], passonly=args.passonly))
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
