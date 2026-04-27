# Bench

The `bench` command is the evaluation module of Verix. 
It benchmarks a set of predicted CSV calls against a ground truth (GT) VCF. 

It resolves topological matches between complex events, calculates standard 
performance metrics (Precision, Recall, F1), and provides diagnostic insights into fragmentation, contiguity, and spurious breakends
(refer to [Bench Output](bench_outputs.md) for a detailled presentation of the command output).

---

## Usage
```bash
verix bench -p <pred.vcf> -fp <format_pred> -t <truth.vcf> -ft <format_truth> -o <out_dir> [OPTIONS]
```

## Required Arguments

* **`-p, --pred`** *(String)*: Path to the caller's VCF file to be evaluated.
* **`-t, --truth`** *(String)*: Path to the ground truth VCF file.
* **`-o, --output`** *(String)*: Path to the output directory. Results will be written to `<output>/output/`.
* **`-ft, --format_truth`** *(Choice: `bnd`, `single_rec`, `multi_rec`)*: The VCF representation format used in the truth set (see [VCF Format Support](vcf_format.md) for a detailed description of each format).
* **`-fp, --format_pred`** *(Choice: `bnd`, `single_rec`, `multi_rec`)*: The VCF representation format used in the caller's file.

---

## Matching & Assignment Parameters

These parameters govern how predicted breakends are mapped to ground truth events.

* **`-mt, --match_threshold`** *(Integer, default: 500)*: 
    Maximum genomic distance (in base pairs) between two breakends for them to be considered a valid match. 
* **`-ubt, --unmatched_bnds_thresh`** *(Integer, default: None)*: 
    The maximum absolute number of unmatched breakends tolerated between two SVs to consider them a candidate match.
* **`-mrt, --match_ratio_thresh`** *(Float, default: None)*: 
    The minimum Jaccard similarity ratio required to merge two SVs: $\frac{\text{matched_breakends}}{\text{total_unique_breakends_in_both_svs}}$.
* **`-bpt, --bp_merge_threshold`** *(Integer, default: 2)*: 
    Breakends within the same SV that are closer than this distance are merged to account for redundancy.
* **`--enforce_svtype`** *(Bool, default: True)*: If True, the algorithm will match breakends even if their SVTYPE (e.g., DEL vs INV) differs.
* **`--enforce_genotype`** *(Bool, default: True)*: If True, genotype discrepancies are ignored during matching.

---

## Parsing & Filtering Options

* **`-csvt, --csv_info_truth`** *(String)*: The `INFO` field tag used to group breakends in the truth VCF. 
    Defaults to `SVID` for `multi_rec` and `BKPS` for `single_rec`.
* **`-csvp, --csv_info_pred`** *(String)*: The `INFO` field tag used to group breakends in the predicted calls VCF.
* **`-st, --svtype`** *(String, default: SVTYPE)*: The `INFO` field tag indicating the structural variant type.
* **`-s, --sizemin`** *(Float, default: 0)*: Filters out SV records with a length strictly smaller than this threshold.
* **`-S, --sizemax`** *(Float, default: None)*: Filters out SV records with a length strictly larger than this threshold.
* **`-chr, --chr_list`** *(List of strings)*: Retains only SV records located on the specified chromosomes.
* **`-q, --qual`** *(Float, default: 0)*: Filters out SV records with a `QUAL` score below this threshold.

> **Note:** Segments are evaluated using their absolute length (the distance between the start and end breakends) disregarding the reported INFO length.

> **Note:** All VCF record filters (`sizemin`, `chr_list`, `qual`) are applied at the record level, not the global CSV level. 
>For example, in a transchromosomal CSV spanning `chr1` and `chr2`, if `chr2` is excluded by `chr_list`, only the breakend on `chr2` is dropped.

---

## Output Files

The tool generates a comprehensive suite of outputs in the `<output>/output/` directory for downstream analysis.
More details and examples are provided at [Bench Outputs](bench_outputs.md).

### Summary Metrics
`summary.json`: The primary metrics report. It contains global Precision, Recall, and F1 scores. 
It also includes average breakend accuracy, total fragmented calls, and spurious breakend counts.

### Relational Mappings
* `closest.csv`: A table identifying the single optimal candidate match for every SV (both from the GT and the Callset 
perspectives), including the detection mode (Complete, Partial, Aggregate, Miss) and the distance scores.
* `detections.csv`: A comprehensive log of all valid candidate overlaps found, not just the optimal ones.
* `truth.csv` & `pred.csv`: DataFrames containing structural characteristics (span, number of breakends, spurious counts) for every SV evaluated.

### Classified VCFs
Verix physically separates the variants into distinct VCF files based on their validation status, appending the `DETECTION_MODE` and match information directly to the `INFO` field of each record:
* `tp-truth.vcf`: Ground truth SVs successfully reconstructed and matched (`Complete` match).
* `fn.vcf`: Ground truth SVs missed or only partially reconstructed (`False Negatives`).
* `tp-pred.vcf`: Caller SVs that perfectly matched a ground truth event (`True Positives`).
* `fp.vcf`: Caller SVs that failed to perfectly match a ground truth event, including aggregates, partials, and total misses (`False Positives`).

---

## Examples

### Standard Benchmarking
Evaluate a caller (in bnd format) against a truth set (multi_rec format) using a standard 500bp distance window for breakend matching.
```bash
verix bench \
    -p caller_output.vcf \
    -fp bnd \
    -t ground_truth.vcf \
    -ft multi_rec \
    -o ./bench_results 
```

### Strict Matching with Quality Filtering
Require breakends to be within 50bp, enforce a minimum quality score of 20, and restrict the evaluation to Chromosome 1 and Chromosome 2.
```bash
verix bench \
    -p caller_output.vcf \
    -fp bnd \
    -t ground_truth.vcf \
    -ft multi_rec \
    -o ./strict_bench_results \
    -mt 50 \
    -q 20 \
    -chr chr1 chr2
```
