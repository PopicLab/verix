## Outputs

### Performance metrics
In `summary.json`, `Verix` reports standard evaluation metrics alongside advanced diagnostic measurements:
**Global Summary:**
* `Valid Calls`: Total number of predicted SVs analyzed.
* `Truth Valid SVs`: Total number of ground truth SVs analyzed.
* `TP-pred`, `TP-truth`, `FP`, `FN`.
* `precision`, `recall`, and `f1`.

**Breakend-Level Fidelity:**
* `Proportion truth breakend matched (all matches)`: The ratio of breakends that found a valid partner over the total number of ground truth breakends.
* `Proportion truth breakend (complete matches)`: The ratio of breakends belonging to SVs with a `complete` detection mode over the total number of ground truth breakends.
* `Mean breakend accuracy (optimal matches)`: The average genomic distance (bp) between partner breakends across all optimal matches.
* `Mean breakend accuracy (complete matches)`: The average genomic distance (bp) between partner breakends across all `complete` matches.

**Error Diagnostics:**
* `Mean GT breakend miss rate (optimal matches)`: Average number of unmatched GT breakends in the optimal match.
* `Mean extra breakend rate (optimal matches)`: Average number of call breakends that are unmatched in the optimal match.
* `Average number spurious breakend rate (all matches)`: Average number of spurious predicted breakends (not matched with any ground-truth breakend).
* `Average number fragmented calls (optimal matches)`: Proportion of optimal matches where a single GT SV is shattered into multiple call SVs.
* `Average number discontiguous calls (optimal matches)`: Proportion of predicted calls with an optimal match skipping internal breakends.
* `Average number full matches (optimal matches)`: Proportion of predicted calls with a full optimal match (all the ground truth breakends are matched but the call might have extra breakends).
* Count of the different modes of the optimal matches for the predicted calls and ground-truth events.

> **Note:** Detailed logic for these metrics can be found in the [Benchmarking Methodology](methodology.md).

Example
```angular2html
----------------------------------------------------------------------------------------------------
            Verix Benchmarking Summary            
----------------------------------------------------------------------------------------------------
Valid Calls                                                            :        880
Truth Valid SVs                                                        :       1000
TP-pred                                                                :        400
TP-truth                                                               :        400
FP                                                                     :        480
FN                                                                     :        600
precision                                                              : 0.45454545
recall                                                                 : 0.40000000
f1                                                                     : 0.42553191
Proportion truth breakend matched (all matches)                        : 0.91892857
Proportion truth breakend (complete matches)                           : 0.29000000
Mean breakend accuracy in bp (optimal matches)                         : 1.94077449
Mean breakend accuracy in bp (complete matches)                        : 1.79821429
Mean GT breakend miss rate (optimal matches)                           : 0.13566667
Mean extra breakend rate (optimal matches)                             : 0.07088156
Average number spurious breakend rate (all matches)                    : 0.00438920
Average number of GT SV single match (all matches)                     : 0.66100000
Average number of GT SV single match (complete matches)                : 0.23400000
Average number fragmented calls (optimal matches)                      : 0.15500000
Average number discontiguous calls (optimal matches)                   : 0.01590909
Average number full matches (optimal matches)                          : 0.79414733
complete-pred                                                          :        560
partial-pred                                                           :        237
aggregate-pred                                                         :         81
miss-pred                                                              :          2
complete-truth                                                         :        400
partial-truth                                                          :        195
aggregate-truth                                                        :        396
miss-truth                                                             :          9
----------------------------------------------------------------------------------------------------
```

### Detection JSONs
Verix generates four nested JSON files with a breakdown of the types of SV matched:
* `count_detections_gt.json` / `count_detections_call.json`: Statistics for all possible overlaps found.
* `count_optimal_gt.json` / `count_optimal_call.json`: Statistics restricted to the optimal match for each SV.

These files use a nested structure: **[Primary Type] $\rightarrow$ [Detection Mode] $\rightarrow$ [Secondary Type]**.
Example of a count_optimal_call.json file. Here the primary type will be the predicted type and the secondary type the truth type the call is matched to.
```angular2html
{
    "DEL": {
        "complete": {
            "DEL": 198
        },
        "miss": {
            "NA": 2
        }
    },
    "DUP_INV": {
        "aggregate": {
            "BND": 2
        },
        "complete": {
            "DUP": 100,
            "INV": 93
        },
        "miss": {
            "NA": 2
        },
        "partial": {
            "BND": 2,
            "INS": 1
        }
    }
}
```

### Annotated VCFs
`Verix` generates annotated VCF files (`fn.vcf`, `fp.vcf`, `tp-pred.vcf`, and `tp-truth.vcf` ) containing the records 
classified into each respective performance category.
These output records contain custom `INFO` fields that describe their specific detection status.

The appended `INFO` fields include:
* `MATCH_LIST`: A list of all matches that overlapped with this record. The format is `match_SV_ID|detection_mode|#unmatched_bnds|Distance_matched_bnds|is_fragmented|Is_contiguous|matched_bnds_positions`
* `OPTIMAL_MATCH_MODE`: The assigned mode of the highest-scoring match (`complete`, `partial`, `aggregate`, or `miss`).
* `OPTIMAL_MATCH_TYPE`: The `SVTYPE` of the partner in the optimal match.
* `OPTIMAL_MATCH_ID`: The ID of the specific partner SV used for the optimal match.
* `OPTIMAL_MATCH_SCORE`: Formatted as `(#BND unmatched - total distance - fragmented flag - Contiguity)`.
* `BND_MATCH`: A map of internal breakend matching, formatted as `Local_BND_ID:Partner_BND_ID:Distance|...`.

```angular2html
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  example.30x
chr1    1464461 6       N       <INV+DUP>       .       .       END=1483489;CHROM2=chr1;SVTYPE=INV+DUP;BKPS=6-chr1:1466741-chr1:1468755;MATCH_LIST=sv878|partial|4|15|1|Discontiguous|chr1:1468755_chr1:1468766|chr1:1483489_chr1:1483493;OPTIMAL_MATCH_MODE=partial;OPTIMAL_MATCH_TYPE=delINVdel;OPTIMAL_MATCH_ID=sv878;OPTIMAL_MATCH_SCORE=4-15-1-Discontiguous;BND_MATCH=.,.,sv878:1468766,sv878:1483493     GT      ./.
chr1    5523565 59      N       <INV>   .       .       END=5532963;CHROM2=chr1;SVTYPE=INV;MATCH_LIST=sv273|complete|0|3|0|Contiguous|chr1:5532963_chr1:5532962|chr1:5523564_chr1:5523566;OPTIMAL_MATCH_MODE=complete;OPTIMAL_MATCH_TYPE=DUP_INV;OPTIMAL_MATCH_ID=sv273;OPTIMAL_MATCH_SCORE=0-3-0-Contiguous;BND_MATCH=sv273:5523566,sv273:5532962      GT      0/1
```

### Tabular Data (CSVs)
For custom analysis, Verix exports all internal DataFrames.

#### `truth.csv` & `pred.csv`
These files contain one row per SV event, summarizing its structural properties and matching:
* `parent_id`
* `parent_type`
* `num_breakends`: Total number of breakends in the event.
* `sv_span`: Total genomic span (bp). `-1` for inter-chromosomal SVs.
* `breakends`: Coordinates of all breakends in the event.
* `matched_breakends`: Coordinates of breakends that successfully matched.
* `num_candidates`: Number of overlapping SVs.
* `spurious`: Number of unmatched breakends.
* `match_type` / `match_id`: Attributes of the matched SV(s).
* `match_num_unmatched_bnds` / `match_total_distance`: Score components of the matched SV(s).
* `optimal_candidate_id`: The ID of the best-scoring matched event.
* `optimal_contiguity`: Contiguity of the best match (`Contiguous`, `Discontiguous`, or `None`).
* `fragmented`: Boolean; true if, in the optimal match, the unmatched ground-truth breakends overlap another call.
* `detection_mode`: Detection mode of the optimal match (`complete`, `partial`, `aggregate`, `miss`).
* `optimal_type`: Parent type of the optimally matched SV.
* `optimal_breakends`: List of matched breakends in the optimal match.
* `optimal_num_unmatched`: Number of breakends that are not matched between the two optimally matched SVs.
* `optimal_total_distance`: Sum of the genomic distances between matched breakends.
* `optimal_full`: Boolean; true if the ground-truth SV's breakends are all matched in the optimal match.

```csv
parent_id,parent_type,num_breakends,detection_mode,sv_span,spurious,num_candidates,breakends,matched_breakends,match_type,match_id,match_num_unmatched_bnds,match_total_distance,optimal_candidate_id,optimal_contiguity,fragmented,optimal_type,optimal_breakends,optimal_num_unmatched,optimal_total_distance,optimal_full,optimal_mode,optimal_miss,optimal_extra
sv561,dDUP,3,partial,96691,1,2,chr1:74129|chr1:168161|chr1:170820,chr1:74129|chr1:168161,,,,,1,Contiguous,TRUE,INS,chr1:74129_chr1:74129,2,0,FALSE,,0.666666667,0
sv799,INVdel,3,aggregate,7472,0,1,chr1:2659885|chr1:2662360|chr1:2667357,chr1:2659885|chr1:2662360|chr1:2667357,,,,,severus_36,Discontiguous,FALSE,BND,chr1:2662360_chr1:2662360|chr1:2659885_chr1:2659885|chr1:2667357_chr1:2667357,12,0,TRUE,,0,0.8
```

#### `detections.csv` & `optimal.csv`
These files describe the matches between SVs. `detections.csv` lists all valid overlaps, 
while `optimal.csv` lists only the optimal matches.
A ground-truth or call SV can participate to several matches and be the optimal match to multiple SVs.
* `GT ID` / `Call ID`: The IDs of the two matched SVs.
* `GT type` / `Call type`: The SV types of the pair.
* `is_optimal_GT` / `is_optimal_call`: Whether the match is the optimal match for the GT / call SV.
* `#BND unmatched`: Number of breakends that are not matched for the two SVs.
* `Total dist BND`: Sum of genomic distances between matched breakends.
* `detection mode`: The detection mode of this specific match among (`complete`, `partial`, `aggregate`, `miss`).
* `fragmented`: Whether this specific match is considered fragmented.
* `Matched BNDS`: A detailed string mapping GT breakend positions to Call breakend positions.

```csv
GT ID,GT type,is_optimal_GT,Call ID,Call type,is_optimal_call,#BND unmatched,Total dist BND,detection mode,fragmented,Matched BNDS
sv121,DEL,TRUE,2,DEL,TRUE,0,1,complete,FALSE,chr1:441653_chr1:441654|chr1:446718_chr1:446718
sv878,delINVdel,TRUE,severus_35,DEL,FALSE,12,5,aggregate,FALSE,chr1:1487562_chr1:1487562|chr1:1474629_chr1:1474629|chr1:1483493_chr1:1483493|chr1:1468771_chr1:1468766
```
