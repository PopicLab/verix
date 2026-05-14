from pysam import VariantFile, VariantHeader

from verix.sv import SV, Breakpoint, Callset
from verix import __version__, __toolname__
import re
from enum import Enum
import logging

# ----- Input -----

class VCFFormat(str, Enum):
    MULTI = "multi"
    SINGLE = "single"
    DEFAULT = "default"

def find_mate(rec):
    match = re.search(r'[\[\]](.+):(\d+)[\[\]]', rec.alts[0])
    if not match: return None
    return match.group(1), int(match.group(2))

def process_bkps(bkps_raw, chrom):
    breakends = []
    for segment in (bkps_raw if isinstance(bkps_raw, (list, tuple)) else [bkps_raw]):
        for bp in segment.split('-')[1:]:
            chrom, pos = bp.split(':') if ':' in bp else (chrom, bp)
            breakends.append((chrom, int(pos)))
    return breakends

def extract_breakends(rec, svid, vcf_format, bkp_link_field):
    breakends = set()
    if vcf_format == VCFFormat.SINGLE and bkp_link_field in rec.info:
        breakends.update(process_bkps(rec.info[bkp_link_field], rec.chrom))
    else:
        breakends.update({(rec.chrom, rec.pos), (rec.chrom, rec.stop)})
    if 'TARGET' in rec.info:
        breakends.add((rec.info.get('TARGET_CHROM', rec.chrom), rec.info['TARGET']))
    return [Breakpoint(chrom=c, pos=p, svid=svid) for c, p in breakends]

def parse_vcf(input_file, vcf_format, csv_link_name, type_name, sizemin, sizemax, merge_threshold, name=None):
    sv_groups = {}
    made2id = {}
    # --- parse and group all CSV records
    vcf = VariantFile(input_file)
    sample_name = name or next(iter(vcf.header.samples), None)
    has_gt = 'GT' in vcf.header.formats
    for rec_idx, rec in enumerate(vcf):
        svid = rec_idx # by default, use the index in the VCF as the SV ID
        sv_type = rec.info.get(type_name, 'NA')
        genotype = next(iter(rec.samples.values()), {}).get('GT', (None, None))
        if vcf_format == VCFFormat.MULTI and csv_link_name in rec.info: # use the provided link ID
            svid = rec.info[csv_link_name]
        # check for BND links
        mate_info = find_mate(rec)
        if mate_info:
            if mate_info in made2id: svid = made2id[mate_info]
            else: made2id[(rec.chrom, rec.pos)] = svid
        bkps = extract_breakends(rec, svid, vcf_format, csv_link_name)
        if svid in sv_groups:
            if sv_groups[svid]['genotype'] != genotype: raise ValueError(f"Genotype conflict grouping record {rec}")
            sv_groups[svid]['types'].add(sv_type)
            sv_groups[svid]['vcf_bp'].extend(bkps)
            sv_groups[svid]['records'].append(rec)
        else:
            sv_groups[svid] = {'genotype': genotype, 'types': {sv_type}, 'vcf_bp': list(bkps), 'records': [rec]}
    if not sv_groups: logging.warning(f'No SVs found in {input_file}')

    svs = []
    for svid, g in sv_groups.items():
        sv = SV.from_records(svid=svid, sample_name=sample_name, bp_merge_threshold=merge_threshold, **g)
        min_size, max_size = sv.get_min_max_size()
        if sizemin <= min_size and max_size <= (sizemax or float('inf')): svs.append(sv)
    if not svs and sv_groups: logging.warning(f'No SVs left after size filtering in {input_file}')
    return Callset(svs, name=sample_name, has_gt=has_gt)

# ----- Output -----

def create_base_header(chrom_set, source_name):
    header = VariantHeader()
    header.add_line('##fileformat=VCFv4.2')
    header.add_line(f'##source={source_name}')
    for chrom in sorted(chrom_set):
        header.add_line(f'##contig=<ID={chrom}>')
    header.add_line('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the SV">')
    header.add_line('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">')
    header.add_line('##INFO=<ID=CHROM2,Number=1,Type=String,Description="Chr of the last breakpoint, if different">')
    header.add_line('##INFO=<ID=BKPS,Number=.,Type=String,Description="List of all SV breakpoints">')
    return header

def populate_rec(rec, sv):
    rec.chrom = sv.start.chrom
    rec.pos = sv.start.pos
    rec.stop = max(sv.start.pos + 1, sv.end.pos) # if the SV is a single breakend
    rec.id = str(sv.id)
    rec.alleles = ("N", "<SV>")
    if sv.end.chrom != rec.chrom: rec.info['CHROM2'] = sv.end.chrom
    rec.info['SVTYPE'] = sv.type
    rec.info['BKPS'] = sv.breakpoints2str()

def write_csv_vcf(callset, vcf_path):
    header = create_base_header(callset.chrom_set, f'{__toolname__}_{__version__}')
    header.add_line('##INFO=<ID=ID,Number=.,Type=String,Description="Unique CSV identifier">')
    if callset.has_gt:
        header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    if callset.sample_name:
        header.add_sample(callset.sample_name)
    with VariantFile(vcf_path, 'w', header=header) as csv_vcf:
        for sv in callset.svs:
            rec = csv_vcf.new_record()
            rec.info['ID'] = str(sv.id)
            populate_rec(rec, sv)
            if callset.has_gt: rec.samples[callset.sample_name]['GT'] = sv.genotype
            csv_vcf.write(rec)

def write_bench_vcf(callset, matches, info_fields, vcf_path):
    header = create_base_header(callset.chrom_set, f'{__toolname__}_bench_{__version__}')
    header.add_line('##INFO=<ID=ID,Number=.,Type=String,Description="Unique CSV identifier">')
    if callset.has_gt:
        header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    for name, (number, typ, desc) in info_fields.items():
        header.add_line(f'##INFO=<ID={name},Number={number},Type={typ},Description="{desc}">')
    if callset.sample_name:
        header.add_sample(callset.sample_name)
    with VariantFile(vcf_path, 'w', header=header) as bench_vcf:
        for sv, match in zip(callset.svs, matches):
            rec = bench_vcf.new_record()
            rec.info['ID'] = str(sv.id)
            populate_rec(rec, sv)
            if callset.has_gt: rec.samples[callset.sample_name]['GT'] = sv.genotype
            for key, val in match.to_vcf_info_dict().items():
                rec.info[key] = val
            bench_vcf.write(rec)

def write_merge_vcf(callsets, sv_clusters, samples, info_fields, vcf_path):
    contigs = set()
    for c in callsets:
        contigs.update(c.chrom_set)
    header = create_base_header(contigs, f'{__toolname__}_merge-{__version__}')
    has_gt = any(c.has_gt for c in callsets)
    if has_gt:
        header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    for name, (number, typ, desc) in info_fields.items():
        header.add_line(f'##INFO=<ID={name},Number={number},Type={typ},Description="{desc}">')
    header.add_line(
        '##FORMAT=<ID=REC,Number=1,Type=String,Description="Merged records from the corresponding sample">')
    for src in samples:
        header.add_sample(src)
    with VariantFile(vcf_path, 'w', header=header) as vcf_out:
        for cluster in sv_clusters:
            rec = vcf_out.new_record()
            populate_rec(rec, cluster.representative)
            for key, val in cluster.to_vcf_info_dict().items():
                rec.info[key] = val
            for sample in samples:
                if has_gt:
                    rec.samples[sample]['GT'] = cluster.representative.genotype
                rec.samples[sample]['REC'] = cluster.serialize_records(sample)
            vcf_out.write(rec)