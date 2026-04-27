from pysam import VariantFile, VariantHeader

from verix import __version__, __toolname__


def create_base_header(chrom_set, source_name):
    """Generates the common VCF header used by both benchmarking and merge."""
    header = VariantHeader()
    header.add_line('##fileformat=VCFv4.2')
    header.add_line(f'##source={source_name}')
    contigs = chrom_set

    for chrom in sorted(contigs):
        header.add_line(f'##contig=<ID={chrom}>')

    # Truth INFO fields shared by both formats
    header.add_line('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the SV">')
    header.add_line('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">')
    header.add_line('##INFO=<ID=CHROM2,Number=1,Type=String,Description="Chromosome of the end breakend">')
    header.add_line('##INFO=<ID=BKPS,Number=.,Type=String,Description="Internal Breakpoints of a complex SV">')
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')

    return header


def populate_bkps_rec(rec, sv):
    # Generate a single record with BKPS notation to represent an SV
    rec.contig = sv.breakends[0].chrom
    rec.start = sv.start.pos
    rec.stop = max(sv.start.pos + 1, sv.end.pos) # If the SV is a single breakend, ensure the position satisfies VCF format
    rec.id = str(sv.parent_id)
    rec.alleles = ('N', f'<{sv.parent_type}>')
    rec.chrom = sv.start.chrom

    rec.info['CHROM2'] = sv.breakends[-1].chrom
    rec.info['SVTYPE'] = sv.parent_type

    if len(sv.breakends) > 2:
        bkps = f"{sv.parent_id}-" + '-'.join([f"{bnd.chrom}:{bnd.pos}" for bnd in sv.breakends[1:-1]])
        rec.info['BKPS'] = bkps


def add_rec_info(rec, sv, is_gt=False):
    formatted_matches = []
    for match in sv.match_list:
        bnd_str = stringify_matched_bnds(match.matched_bnds)
        other_sv = match.call_sv if is_gt else match.gt_sv
        formatted_matches.append(f"{other_sv.parent_id}|{match.mode}|{match.diff_bnd}|{match.total_dist}|"
                                 f"{int(match.fragmented)}|{match.contiguity}|{bnd_str}")
    if formatted_matches:
        rec.info['MATCH_LIST'] = tuple(formatted_matches)

    if sv.optimal_candidate:
        other_sv = sv.optimal_candidate.call_sv if is_gt else sv.optimal_candidate.gt_sv
        rec.info['OPTIMAL_MATCH_MODE'] = sv.optimal_candidate.mode
        rec.info['OPTIMAL_MATCH_TYPE'] = other_sv.parent_type
        rec.info['OPTIMAL_MATCH_ID'] = str(other_sv.parent_id)
        rec.info['OPTIMAL_MATCH_SCORE'] = (f"{sv.optimal_candidate.diff_bnd}-{sv.optimal_candidate.total_dist}-"
                                           f"{int(sv.optimal_candidate.fragmented)}-{sv.optimal_candidate.contiguity}")

    for sample_name, gt in sv.genotype.items():
        rec.samples[sample_name]['GT'] = gt

    best_matches = []
    for bnd in sv.breakends:
        best_match = None
        best_dist = None
        for bnd_list in bnd.bnd_matches.values():
            for match_bnd in bnd_list:
                dist = abs(match_bnd.pos - bnd.pos)
                if best_dist is None or best_dist > dist:
                    best_match = (match_bnd.parent_id, match_bnd.pos)
                    best_dist = dist
        if not best_match:
            best_match = '.'
        else:
            best_match = str(best_match[0]) + ':' + str(best_match[1])
        best_matches.append(best_match)
    if best_matches:
        rec.info['BND_MATCH'] = tuple(best_matches)


def stringify_matched_bnds(matched_bnds):
    return "|".join(f"{stringify_bnd(b1)}_{stringify_bnd(b2)}" for b1, b2 in matched_bnds) if matched_bnds else "None"

def stringify_bnd(bnd):
    return f"{bnd.chrom}:{bnd.pos}"

def get_rec_parent_id(rec, pid_name):
    if pid_name not in rec.info:
        return rec.id
    return rec.info[pid_name]


def add_header_info(header, infos_fields):
    for info_field in infos_fields:
        header.info.add(*info_field.values())


def write_bench_vcf(callset, fp_path, tp_path, is_gt):
    header = create_base_header(callset.chrom_set, f'{__toolname__}_bench_{__version__}')

    VCF_BENCH_INFOS = [
        ('DETECTION_MODE', '1', 'String', 'Mode of detection in {complete, partial, aggregate, miss}'),
        ('MATCH_LIST', '.', 'String', 'List of matches'),
        ('OPTIMAL_MATCH_MODE', '1', 'String', "Mode of detection of the SV's best match"),
        ('OPTIMAL_MATCH_TYPE', '1', 'String', "Type of the SV's best match"),
        ('OPTIMAL_MATCH_ID', '1', 'String', "ID of the SV's best match"),
        ('OPTIMAL_MATCH_SCORE', '1', 'String', 'Score of the best match'),
        ('BND_MATCH', '.', 'String', 'SV ID: distance of the optimal matched bnd for each bnd'),
    ]
    for uid, num, typ, desc in VCF_BENCH_INFOS:
        header.add_line(f'##INFO=<ID={uid},Number={num},Type={typ},Description="{desc}">')

    for sample in callset.sample_names:
        header.add_sample(sample)

    with VariantFile(fp_path, 'w', header=header) as fp_vcf, \
         VariantFile(tp_path, 'w', header=header) as tp_vcf:
        for sv in callset.svs:
            rec = fp_vcf.new_record()
            populate_bkps_rec(rec, sv)
            add_rec_info(rec, sv, is_gt)

            if sv.match_id:
                tp_vcf.write(rec)
            else:
                fp_vcf.write(rec)


def write_merge_vcf(merge_callset, ordered_sources, output_folder):
    vcf_path = f"{output_folder}/merge.vcf"

    header = create_base_header(merge_callset.chrom_set, f'{__toolname__}_merge-{__version__}')

    sources_str = ", ".join(ordered_sources)
    header.add_line(
        '##INFO=<ID=SUPP,Number=1,Type=Integer,Description="Number of distinct callers supporting the variant">')
    header.add_line(
        f'##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="Binary vector of supporting samples. Order: {sources_str}">')
    header.add_line(
        f'##INFO=<ID=SUPP_COUNTS,Number=1,Type=String,Description="Counts contributed by each caller. Order: {sources_str}">')
    header.add_line(
        '##FORMAT=<ID=SRC,Number=1,Type=String,Description="Original records from this source joined by semi-colons">')

    for src in ordered_sources:
        header.add_sample(src)

    with VariantFile(vcf_path, 'w', header=header) as vcf_out:
        for merge_sv in merge_callset.merge_list:
            # Write a VCF record for each merge SV
            rec = vcf_out.new_record()
            populate_bkps_rec(rec, merge_sv.representative_sv)

            # Add INFO tags
            rec.info['SUPP'] = merge_sv.support_count
            rec.info['SUPP_VEC'] = merge_sv.support_vec
            rec.info['SUPP_COUNTS'] = merge_sv.support_count_vec

            for src in ordered_sources:
                rec.samples[src]['SRC'] = merge_sv.stringified_recs[src] if src in merge_sv.stringified_recs else '.'

            vcf_out.write(rec)
