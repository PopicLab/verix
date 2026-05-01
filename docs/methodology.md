## Benchmarking Methodology

This section describes the logic used to evaluate the accuracy of complex structural variant (CSV) calls against a ground truth (GT) set.

### SV Alignment and Matching Logic

The benchmarking process begins by aggregating breakends, potentially across multiple VCF records, to recover the complexity of each CSV.
Please refer to the [Inputs Documentation](https://www.google.com/search?q=inputs.md) for more details about the VCF format.

#### Breakend Matching

Two breakends are considered a match if their genomic coordinates fall within a user-defined `matching_threshold` (default: 500bp).
A ground-truth call is considered a match for a call SV if they share at least one matching breakend.

* **SV Type Matching:** By default, SV type is disregarded when identifying candidates. This can be enforced by using the flag `--enforce_svtype`.
* **Genotype Matching:** Genotype is also disregarded by default. Enforce matching by using the flag `--enforce_genotype`.
* **Candidate Constraints:** While a single matching breakend qualifies for matching, stricter filters can be applied 
using `--unmatched_bnds_thresh` (maximum number of unmatched breakends) and `--match_ratio_thresh` (minimum ratio of matched breakends).

#### Match Scoring
The quality of a candidate is ranked using in order:
1.  **Unmatched Breakends:** The primary score is the total number of unmatched breakends between the two SVs.
2.  **Total Distance:** The secondary score is the sum of genomic distances across all aligned breakend pairs.

In particular, this scoring is used to identify the optimal match for each call and ground-truth SV.

#### Classification of Detection Types

The matches are classified according to the cardinal of the overlap and whether the call overlaps other ground-truth SVs.

| Detection Type | Definition                                                                                                                                                                   |
|----------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Complete**   | Perfect correspondence. The number of predicted breakends equals the number of GT breakends, and all are matched                                                             |
| **Aggregate**  | The predicted call maps to multiple distinct ground-truth events.                                                                                                            |
| **Partial**    | The call represents an incomplete match (e.g., covering only a subset of the GT breakends, or containing extra unmatched breakends) that does not bridge multiple GT events. |
| **Miss**       | No candidate breakend matches were found within the threshold.                                                                                                               |

#### Error diagnostic 

To further stratify reconstruction error types among incomplete calls, `Verix` uses four diagnostic properties:

  * **Fragmented:** Highlights ground-truth CSV split into multiple pieces.
  * **Spurious:** Calculates the number of predicted breakends within a call that fail to map to any ground truth event.
  * **Full:** Characterizes whether all the breakends of the ground-truth events were matched.
  * **Contiguity:** Whether the breakends are matched out of order (i.e., skipping an internal breakend in either the GT or the call).


### Definitions: True / False Positive and False Negative

Only `Complete` calls can be considered TP.
To defined TP, FP and FN, a greedy 1-to-1 assignment is performed.
Each GT SV is matched to its optimal available `Complete` call which is then removed from the pool of available predictions.

Hence, the definitions:
  * **True Positive (TP):** A `Complete` call that is successfully matched to a ground-truth SV.
  * **False Positive (FP):** Any other call.
  * **False Negative (FN):** Any ground-truth SV that fails to match a `Complete` predicted call.

> **Note:** While incomplete matches (`Partial`, `Aggregate`) provide valuable diagnostic information, they are strictly penalized as FPs/FNs to maintain a rigorous standard for complex variant reconstruction.
