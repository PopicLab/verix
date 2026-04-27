from pysam import VariantFile

from verix.sv import SV, Breakend


def find_mate(rec, bnd_match2idx, index):
    # Pair a BND with its mate if already found
    parent_id = index
    pos_mate = rec.alts[0]
    if ':' not in rec.alts[0]:
        pos_mate = rec.ref

    if ':' not in pos_mate or '[' not in pos_mate:
        # No mate BND specified with a standard notation
        return None

    pair_pos = int(pos_mate.split(':')[1].split('[')[0].split(']')[0]) - 1
    pair_chrom = pos_mate.split(':')[0]
    bracket = ']'
    if bracket not in pair_chrom:
        bracket = '['
    pair_chrom = pair_chrom.split(bracket)[1]
    pair_coord = (pair_pos, pair_chrom, rec.start, rec.chrom)

    current_mate = (rec.start, rec.chrom, pair_pos, pair_chrom)
    if pair_coord in bnd_match2idx:
        parent_id = bnd_match2idx[pair_coord]
    bnd_match2idx[current_mate] = parent_id
    return parent_id


def process_bkps(bkps_raw, chrom):
    bnd_positions = []
    if isinstance(bkps_raw, str):
        bkps_raw = [bkps_raw]

    for segment in bkps_raw:
        # Segment format is 'TAG-bp-list'
        for bp in segment.split('-')[1:]:
            # bp can have the format pos or chrom:pos
            pos_bp = bp
            chrom_bp = chrom
            if ':' in bp:
                chrom_bp, pos_bp = bp.split(':')
            bnd_positions.append((chrom_bp, pos_bp))
    return bnd_positions


def parse(input_file, vcf_format, csv_info, svtype_name='SVTYPE', merge_threshold=1,
          chr_list=None, sizemin=0, sizemax=None, qual=None, source=''):
    svs = {}
    bnd_match2idx = {}
    # Keep track of the source file and ensure the SV_id is unique for merge
    source_cpt = source + '_' if source else ''
    for rec_idx, rec in enumerate(VariantFile(input_file)):
        chrom = rec.chrom

        if chr_list and chrom not in chr_list: continue
        if isinstance(rec.qual, int) and rec.qual < qual: continue

        start = rec.start
        stop = rec.stop
        svtype = rec.info[svtype_name] if svtype_name in rec.info else \
            (rec.info['SVTYPE'] if 'SVTYPE' in rec.info else 'unspecified')

        genotype = {
            sample_name: sample_obj.get('GT', (None, None))
            for sample_name, sample_obj in rec.samples.items()
        }

        chrom2 = chrom
        if 'CHROM2' in rec.info:
            chrom2 = rec.info['CHROM2']
        bnd_positions = [(chrom, start)]

        svid = f"{source_cpt}{rec_idx}"

        # CSVs are represented as multiple records linked through the parent ID csv_info field
        if vcf_format == 'multi_rec' and csv_info in rec.info:
            svid = source_cpt + str(rec.info[csv_info])

        if start != stop:
            # If the record wasn't a BND
            bnd_positions.append((chrom2, stop))
            if chrom == chrom2:
                sv_length = stop - start
                if sv_length and (sv_length < sizemin or (sizemax and sv_length > sizemax)): continue
        else:
            # If the format is default or multi_rec, the breakend notation is allowed to represent SVs
            mate_id = find_mate(rec, bnd_match2idx, svid)

            if vcf_format == 'multi_rec' and csv_info in rec.info:
                assert mate_id == svid, ('If the breakend notation is used jointly to a '
                                         f'multi record notation through cluster ID, mate '
                                         f'breakends must present the same parent ID {csv_info}')
            if mate_id is not None:
                svid = mate_id

        # Support dispersion with the TARGET INFO field
        if 'TARGET' in rec.info:
            target = rec.info['TARGET']
            target_chrom = rec.info.get('TARGET_CHROM', chrom)
            bnd_positions.append((target_chrom, target))

        # All the internal breakends are reported in a same line using the csv_info field
        if vcf_format == 'single_rec' and csv_info in rec.info:
            bnd_positions += process_bkps(rec.info[csv_info], chrom)

        breakends = [Breakend(chrom=chrom, pos=pos, parent_id=svid) for chrom, pos in bnd_positions]

        # Recover other potential breakends of a CSV
        if svid not in svs:
            svs[svid] = SV(parent_id=svid, parent_type=svtype,
                           breakends=breakends, bp_merge_threshold=merge_threshold,
                           source=source, records=[rec], genotype=genotype)
        else:
            if vcf_format == 'multi_rec':
                assert svs[svid].genotype == genotype, (
                    f"Records with the same parent ID must have the same genotype. "
                    f"Conflict on SV {svid} (Record ID {rec_idx}): "
                    f"Expected {svs[svid].genotype}, but got {genotype}."
                )

            parent_bp = svs[svid].breakends
            svs[svid].breakends = parent_bp + breakends
            svs[svid].records.append(rec)

    # Sort and deduplicates the breakends
    svs = list(svs.values())
    for sv in svs:
        sv.update_breakends()
    return svs
