# Verix: a toolkit for benchmarking and harmonization of complex structural variants

## Table of Contents

* [Overview](#overview)
* [Key Functionality](#func)
* [Installation](#install)
* [Quick Start](#start)
* [User Guide](#guide)


<a name="overview"></a>
## Overview
`verix` is a lightweight toolkit for benchmarking and integrating complex structural variant (CSV) callsets. 
Each CSV is represented as a set of breakpoints, and comparing two events is formulated as min-weight bipartite 
matching between these breakpoint sets that (1) maximizes the number of matched breakpoint pairs and (2) minimizes the total
distance between them. Two breakpoints are eligible to match only if they fall within a configurable distance threshold, 
and matching can be additionally restricted to event pairs whose SV type and genotype agree. In benchmarking mode, 
each query event is classified according to how its breakpoints line up with the target truthset, distinguishing 
between complete reconstruction, partial capture, over-aggregation (i.e., call that collapsed multiple target events 
into one), and fully spurious calls. To handle diverse callsets, `verix` supports multiple input VCF record 
linking conventions to group related records into a single complex event. `verix` outputs annotated VCFs with 
per-event match details, alongside comprehensive summary statistics and diagnostic plots.


<a name="func"></a>
### Key Functionality

* `bench`: Compares a query VCF against a truth VCF and assigns each query event a match class 
(complete, partial, aggregate, or spurious). Outputs precision, recall, and F1 based on complete matches, along with 
per-class breakpoint accuracy and hit-rate metrics, an annotated VCF with detailed match information 
(e.g. breakpoint alignment) for every event, and summary plots.
* `consensus`: Collapses matching events from one or more VCFs (e.g., from different callers or samples) into a single 
representative call. Generates an integrated VCF with per-sample support annotations and summary statistics on 
call concordance and support patterns.

<a name="install"></a>
## Installation

* Clone the repository: `git clone git@github.com:PopicLab/verix`
* Navigate into the top-level folder: `cd verix`
* Install the framework: `pip install .`
* Set the PYTHONPATH: `export PYTHONPATH=$PYTHONPATH:/path/to/verix/`

<a name="start"></a>
## Quick Start
1. Compare two VCFs: 

`verix bench --query predictions.vcf --target truth.vcf --output_dir results/`

Key outputs: `results/matches.vcf` with per-CSV match annotations and `results/report.json` with summary statistics. 

2. Integrate multiple VCFs into a single consensus callset:

`verix consensus --inputs caller_a.vcf caller_b.vcf caller_c.vcf --names A B C --output_dir results/`

Key outputs: `results/merged.vcf` with matching calls collapsed into one record and `results/report.json` 
with summary statistics.

<a name="guide"></a>
## User Guide

For detailed information about how to use `verix`, please refer to the [User Guide](docs/user_guide.md).