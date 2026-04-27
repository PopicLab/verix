# Consensus
The consensus command merges CSV callsets from multiple VCF files (e.g., different callers or samples) 
into a single consensus VCF.

`Verix` identifies clusters of SVs that share structurally similar breakend topologies, and handles fragmented calls.
The representative SV selected for each cluster is the one with the most breakends.

## Usage
```bash
verix consensus -i <vcf1> <vcf2> ... -f <format1> <format2> ... -o <out_dir> [OPTIONS]
```

## Required Arguments

* **`-i, --inputs`** *(List of strings)*: Space-separated paths to the input VCF files to be merged.
* **`-f, --formats`** *(List of choices: `bnd`, `single_rec`, `multi_rec`)*: The VCF representation format for each input file. The number of formats must match the number of inputs.
* **`-o, --output`** *(String)*: Path to the output directory. Results will be written to `<output>/output/`.

## Clustering & Matching Thresholds

These parameters govern how strict the algorithm is when deciding if two SVs belong to the same consensus event.

* **`-mt, --match_threshold`** *(Integer, default: 500)*: 
    The maximum base-pair distance between two breakends for them to be considered a match.
* **`-ubt, --unmatched_bnds_thresh`** *(Integer, default: 0)*: 
    The maximum absolute number of unmatched breakends tolerated between two SVs to merge them.
* **`-mrt, --match_ratio_thresh`** *(Float, default: None)*: 
    The minimum ratio of matched breakends required to merge two SVs: $\frac{\text{matched_breakends}}{\text{total_unique_breakends}}$.
* **`-bpt, --bp_merge_threshold`** *(Integer, default: 1)*: 
    Distance threshold to merge near-duplicate breakends within the same input SV before cross-matching.
* **`-src, --sources`** *(List of strings, default: None)*: 
    Custom names for each input file to be used in the output VCF metadata. If not provided, the filename and the position of the file in the list are used.
* **`--ignore_svtype`** *(Bool, default: True)*: If True, SVs can be clustered even if their SVTYPE labels differ.
* **`--ignore_genotype`** *(Bool, default: True)*: If True, matching ignores genotype information.

## Input Parsing & Filtering Options

* **`-csv, --csv_info_list`** *(List of strings)*: The `INFO` field tags used for SVID/BKPS for each input file. 
If omitted, defaults are chosen based on the file format (`SVID` for `multi_rec`, `BKPS` for `single_rec`).
* **`-st, --svtype`** *(String, default: SVTYPE)*: The `INFO` field tag indicating the variant type.
* **`-s, --sizemin`** *(Float, default: 0)*: Filter out input SVs with a span smaller than this threshold.
* **`-S, --sizemax`** *(Float, default: None)*: Filter out input SVs with a span larger than this threshold.
* **`-chr, --chr_list`** *(List of strings)*: Restrict analysis to specific chromosomes (e.g., `-chr chr1 chr2`).
* **`-q, --qual`** *(Float, default: 0)*: Minimum `QUAL` score to retain an input record.

---

## Output Files

The tool generates two primary outputs in the specified `<output>/output/` directory:

### `consensus_merged.vcf`
A representative VCF for each consensus cluster. 
The VCF format follows:

**Key INFO fields:**
* `SVTYPE` (or `--svtype` if provided): The lexicographically sorted, hyphen-joined types (e.g., `BND-DEL`).
* **`SUPP`**: The total number of input files that support this consensus SV.
* **`SUPP_VEC`**: A binary string (e.g., `101`) indicating which input sources contributed to the cluster.
* **`SUPP_COUNTS`**: A hyphen-separated string of the exact number of SV records each caller contributed to the cluster (e.g., `2-0-1`).

**Samples Field from original records:**
* `SRC`: One column per source, details the original records from this source (ID|chrom:pos-chrom2:end|SVTYPE|PARENTID) joined by commas.

```angular2html
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  severus svision sniffles
chr1    87751608        severus_1013    N       <INV>   .       .       END=87753971;CHROM2=chr1;SVTYPE=INV;SUPP=3;SUPP_VEC=101;SUPP_COUNTS=1-0-2       SRC     severus_INV688|chr1:87751608|chr1:87753971|INV|severus_1013     .       Sniffles2.INV.8CAS0|chr1:87751609|chr1:87753971|INV|sniffles_659,Sniffles2.DUP.60BS0|chr1:87751615|chr1:87753974|DUP|sniffles_660
```
---

## Examples

### Basic Usage
Merge three callers using a standard 500bp distance between breakends. The SVs must have all their breakends matching 1-to-1 to be merged.
```bash
verix consensus \
    -i sniffles.vcf svision.vcf severus.vcf \
    -f bnd single_rec multi_rec \
    -o ./consensus_results \
    -mt 500
```

### High-Fidelity Consensus with Custom Names
Require breakends to be within 100bp, allow a maximum of 2 unmatched breakends between callers, and demand at least 75% 
of the breakends to be matched to accept the merge.
```bash
verix consensus \
    -i caller1.vcf caller2.vcf caller3.vcf \
    -f format1 format2 format3 \
    -src ONT_Caller Illumina_Caller \
    -o ./strict_consensus \
    -mt 100 \
    -ubt 2 \
    -mrt 0.75 
```
