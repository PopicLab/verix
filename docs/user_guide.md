## User Guide

### Table of Contents
[CLI](#cli)  
[VCF Input Formats](#inputs)  
[Outputs](#outputs)  
&emsp;&emsp;[Benchmarking](#bench)  
&emsp;&emsp;&emsp;&emsp;[Annotated VCF](#benchvcf)  
&emsp;&emsp;&emsp;&emsp;[Stats](#benchstats)  
&emsp;&emsp;&emsp;&emsp;[Plots](#benchstats)  
&emsp;&emsp;[Consensus](#merge)  
&emsp;&emsp;&emsp;&emsp;[Annotated VCF](#mergevcf)  
&emsp;&emsp;&emsp;&emsp;[Stats](#mergestats)
 
<a name="cli"></a>
### CLI

The full list of `verix` parameters is provided below.

`verix bench`: compares a query VCF against a target (truth) VCF
```bash
usage: verix bench [options] -q QUERY_VCF -t TARGET_VCF -o OUTPUT_DIR
params:
-h, --help            show this help message and exit
-q, --query           VCF file with query SVs
-t, --target          VCF file with target SVs
--plot                Generate PDF report with summary plots (default: False)
-o, --output_dir      Output directory
-d, --match_thr       Max distance between matching breakpoints (default: 500)
-s, --sizemin         Minimum SV interval size (default: 0)
-S, --sizemax         Maximum SV interval size (default: None)
-b, --merge_thr       Collapse breakends in a CSV within this distance into a single breakpoint (default: 1)
--enforce_type        Require SV types to match (default: False)
--enforce_genotype    Require SV genotypes to match (default: False)
-f, --formats {multi,single,default} [{multi,single,default} ...]
                      Format type for each VCF (for bench: query, target) (default: default[default...])
-l, --csv_links LINK [LINK ...]
                      INFO field for CSV linking in each VCF (for bench: query, target) (default: [])
-svt, --types TYPES [TYPES ...]
                      INFO field for SV type extraction (default: SVTYPE[SVTYPE...])
```

`verix consensus`: merges multiple VCF files into a single consensus callset
```bash
usage: verix consensus [options] -i VCF [VCF ...] -o OUTPUT_DIR
params:
-h, --help            show this help message and exit
-i, --inputs VCF [VCF ...]  List of VCF files to merge
-n, --names NAME [NAME ...]  Ordered list of names for each VCF (default: [])
-o, --output_dir      Output directory
-d, --match_thr       Max distance between matching breakpoints (default: 500)
-s, --sizemin         Minimum SV interval size (default: 0)
-S, --sizemax         Maximum SV interval size (default: None)
-b, --merge_thr       Collapse breakends in a CSV within this distance into a single breakpoint (default: 1)
--enforce_type        Require SV types to match (default: False)
--enforce_genotype    Require SV genotypes to match (default: False)
-f, --formats {multi,single,default} [{multi,single,default} ...]
                      Format type for each VCF (for bench: query, target) (default: default[default...])
-l, --csv_links LINK [LINK ...]
                      INFO field for CSV linking in each VCF (for bench: query, target) (default: [])
-svt, --types TYPES [TYPES ...]
                      INFO field for SV type extraction (default: SVTYPE[SVTYPE...])
```

#### Core parameters


- `--match_thr`: a pair of breakpoints is considered a candidate match if they are at most this many base pairs apart

- `--merge_thr`: breakends from the same CSV record within this distance are collapsed into a single breakpoint at the median position

- `--sizemin` / `--sizemax`: size filters applied to intervals between consecutive breakpoints on the same chromosome within a CSV; 
events with at least one interval smaller than `--sizemin` or larger than `--sizemax` are removed


<a name="inputs"></a>
### VCF Inputs

`verix` assembles CSVs by linking related VCF records. The `--formats` parameter can be used to specify 
how to group input records in each VCF: 

- **`default`**: collects breakpoints from `POS`, `END`, the `TARGET` INFO field, and BND mates in the ALT field

- **`single`**: uses `default` fields and a custom INFO field containing internal breakpoints specified using `--csv_links` 
(see details on the field format below)

- **`multi`**:uses `default` fields and links records that share the same ID specified as a custom INFO field using 
`--csv_links`

**Note**: formats can be mixed across inputs but not within the same VCF file.

**Genotype parsing**: if each input VCF file contains multiple samples, only the genotype of the first sample is retained.

**SV type parsing**: the `SVTYPE` field (or the field specified via `--types`) is used to extract the type from each 
linked record; if a CSV spans records with multiple distinct types, the values are joined with `+` in sorted order 
to form a consolidated type string.

###### Expected INFO field structure for the `single` record VCF format 

The INFO field provided using `--csv_links` (storing the internal breakpoints of a CSV record) is expected in the
following format:

```<prefix>-<bp_1>[-<bp_2>...]```

where each `<bp_i>` is either `<chr>:<pos>` or just `<pos>` (defaults to the record's `CHROM`). 
Note: the leading `<prefix>` is ignored by the `verix` parser (anything before the first `-` is discarded )

### Outputs

All `verix` outputs are written to the `--output_dir` folder, which is created if it doesn't exist. 
Each run also produces a `main.log` file recording input parameters and execution details.
The outputs of each command are described below.

<a name="bench"></a>
#### Benchmarking

The `bench` command generates the following files:

- `matches.vcf`: query SVs annotated with match information (see below)
- `report.json`: summary statistics
- `report.pdf`: diagnostic plots (only with `--plot`)
- `query.vcf` and `target.vcf`: SVs assembled from the query and target input VCFs, respectively 
(in a `single` record VCF format, with all breakpoints listed in the `BKPS` field)

<a name="benchvcf"></a>
##### Annotated VCF (`matches.vcf`)

Contains the assembled query SVs (one line per SV). Each record includes:

- `SVTYPE`: consolidated SV type
- `BKPS`: comma-separated list of all consolidated breakpoints (```chr:pos,...```)
- `CHROM2`: chromosome of the last breakpoint, set only when it differs from the record's `CHROM`

Match-related INFO fields:

- `BEST_MATCH_CLASS`: overall match classification, one of:
  - `complete`: every query breakpoint matched every target breakpoint
  - `partial`: some query breakpoints matched and either the rest matched no target or target has extra points
  - `aggregate`: distinct subsets of query breakpoints matched multiple distinct target events
  - `spurious`: no query breakpoints matched any target
- `SPURIOUS`: number of query breakpoints that matched no target (omitted when 0)
- `FRAGMENTED`: flag set when the optimal target was also matched by other queries (only for target breakpoints that 
were not matched by this query)
- `MATCHES`: pipe-separated list of all candidate alignments considered (see format below)

For records not classified as `spurious`, the following fields describe the optimal alignment:
- `BEST_MATCH_ID`, `BEST_MATCH_TYPE`: target event ID and SV type
- `BEST_N_MATCHED`, `BEST_BND_DIST`: number of matched breakpoints and total breakpoint distance
- `BEST_MATCH_COV`: `full` if every target breakpoint was matched, `partial` otherwise

###### Candidate alignment format
Each alignment entry in `MATCHES` is a comma-separated tuple:
```
<target_ID>,<target_SVTYPE>,<n_matched>,<bp_dist>,<target_coverage>,<bp_pair_1>,<bp_pair_2>,...
```
where: 
- `<target_ID>`, `<target_SVTYPE>`: ID and SV type of the matched target event
- `<n_matched>`: number of matched breakpoints
- `<bp_dist>`: total breakpoint distance
- `<coverage>`: `full` if every target breakpoint was matched, `partial` otherwise
- `<bp_pair_i>`: matched breakpoint pair `<query_chr:pos>-<target_chr:pos>` (e.g. `chr1:817452-chr1:817452`)

###### Example record

A query CSV on chr2 with four breakpoints, classified as `aggregate` because two distinct subsets of 
its breakpoints align to two different target events, the best alignment has the smaller breakpoint distance:

```
chr2  1450200  q_42  END=1492800;SVTYPE=INVDUP;BKPS=chr2:1450200,chr2:1471050,chr2:1488300,chr2:1492800;
BEST_MATCH_CLASS=aggregate;BEST_MATCH_ID=truth_88;BEST_MATCH_TYPE=INV;BEST_MATCH_COV=full;BEST_N_MATCHED=2;BEST_BND_DIST=15;
MATCHES=truth_88,INV,2,15,full,chr2:1450200-chr2:1450195,chr2:1471050-chr2:1471060|truth_91,DUP,2,42,full,chr2:1488300-chr2:1488280,chr2:1492800-chr2:1492840
```

<a name="benchstats"></a>
##### Stats (`report.json`)

A report with various summary statistics, including: 

- `n_query`: number of query SVs after parsing and size filtering
- `n_target`: number of target SVs after parsing and size filtering  
- `tp-query`, `tp-target`, `fp`, `fn`, `precision`, `recall`, `f1`; note: only `complete` matches count as TP 
- `class_proportions`: shows the fraction of query events in each match category
- `by_class` maps each non-spurious match category to metrics computed over the optimal alignment of each query in this category:
  - `num_matches`: number of query records 
  - `mean_breakpoint_distance`: mean per-breakpoint distance 
  - `mean_breakpoint_hit_rate`: mean fraction of target breakpoints matched 
  - `mean_spurious_breakpoint_rate`: mean fraction of query breakpoints that matched no target 
  - `mean_targets_per_record`: mean number of distinct target events appearing in the candidate alignments of each query 
  - `query_type_counts`: count of query records in this category broken down by predicted SV type 
  - `query_type_proportions`: per SV type, the fraction of all query records of that type that fell into this category 
  - `target_type_counts`: count of unique target events in this category broken down by truthset SV type 
  - `target_type_proportions`: per SV type, the fraction of all target events of that type that fell into this category 
<a name="benchplots"></a>
##### Plots

With `--plot`, `verix bench` writes a multi-page PDF with various diagnostic plots summarizing how query events
aligned to the target, which includes breakpoint-level accuracy and coverage, the prevalence of spurious and fragmented
calls, the correspondence between query and target SV types -- stratified by SV type and match category.

<a name="merge"></a>
#### Consensus

The `consensus` command generates the following files: 

- `merged.vcf`: the consensus callset (see below)
- `report.json`: summary statistics

<a name="mergevcf"></a>
##### Annotated VCF (`merged.vcf`)

Contains one record per SV cluster (group of matching SVs); the SV with the most breakpoints 
is used as the cluster representative (its `CHROM`, `POS`, `END`, `SVTYPE`, `BKPS`, and optionally `CHROM2` 
are used to write the record).

Consensus-related INFO fields:

- `SUPPORT`: number of distinct input callsets contributing SVs to this cluster
- `SUPPORT_COUNT`: comma-separated count vector, one entry per input VCF, giving the number of SVs contributed to this cluster by this input
- `SUPPORT_BINARY`: concatenated 0/1 string, where 1 indicates that the input VCF has SV in this cluster, 0 otherwise

A per-sample `REC` FORMAT field lists the input records merged into each consensus call. 
For each input sample, the value is either `.` (no contribution) or a `|`-separated list of records, each formatted
as `<sv_id>,<sv_type>,<bp1>,<bp2>,...`.

<a name="mergestats"></a>
##### Stats (`report.json`)

A report with various input and consensus-based stats, including: 
- `n_total_variants`: total number of variants summed across all input callsets (after parsing and size filtering)
- `n_variants_in_sample`: maps each input sample name to its variant count
- `n_clusters`: number of consensus clusters
- `support_vec_types`: list of distinct support patterns observed across clusters
- `max_cluster_size`, `min_cluster_size`: maximum and minimum total number of input SVs contained in any cluster


