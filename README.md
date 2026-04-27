# Verix: a toolkit for benchmarking and harmonization of complex structural variants

**Verix** is a principled toolkit designed to benchmark, merge, and evaluate **Complex Structural Variants (CSV)** calls.

## Table of Contents

* [Overview](#Overview)
* [Key Functionality](#Key-Functionality)
* [Installation](#Installation)
* [Quick Start](#Quick-Start)
* [Documentation](#Documentation)
* [Citation](#Citation)
* [Authors](#Authors)

-----

## Overview
`Verix` CSV evaluation and merging tool.
It uses a set-of-breakends matching paradigm and a graph-based algorithm designed to handle the complexity of CSVs. 

In `bench` mode, `Verix` evaluates the accuracy of predicted SVs against a ground-truth dataset.
It provides detailed performance statistics, annotated VCF/CSV files, and diagnostic plots.

In `merge` mode, `Verix` consolidates redundant SV calls within or across multiple VCF files into a single consensus set.

## Key Functionality

  * **`bench`**: Compares a predicted call VCF `pred` against a ground-truth VCF `truth`. 
It characterizes SV matches in several categories for a thorough diagnostic.
  * **`merge`**: Combines complex SV callsets from a single or multiple VCF files (e.g., across different callers or samples) 
into a representative SV relying on the connected components in an SV overlap graph.


-----

## Installation

1.  **Clone the repository:**

<!-- end list -->

```bash
git clone git@github.com:PopicLab/verix
cd verix
```

2.  **Install the package:**

<!-- end list -->

```bash
pip install .
```

3.  **Set your PYTHONPATH:**

<!-- end list -->

```bash
export PYTHONPATH=$PYTHONPATH:/path/to/verix/
```

-----

## Quick Start
* **Format your VCFs:** Ensure your VCFs follow the specifications in the [Inputs Documentation](https://www.google.com/search?q=Documentation/vcf_format.md).
* **Run the benchmark:**



```bash
verix bench \
  --pred caller_output.vcf \
  --format_pred single_rec \
  --truth ground_truth.vcf \
  --format_truth multi_rec \
  --output ./results_path/ \
  --plot
```

* **Run the merge:**
```bash
verix merge \
  --inputs vcf1.vcf vcf2.vcf vcf3.vcf\
  --formats format1 format2 format3 \
  --output ./results_path/ 
```

<!-- end list -->
-----

## Documentation

For detailed information about how to use `Verix`, please refer to the following documentation:

* [Methodology](Documentation/methodology.md)
* [VCF Format](Documentation/vcf_format.md)
* [Benchmarking](Documentation/bench.md)
* [Benchmarking Outputs](Documentation/bench_outputs)
* [Merge](Documentation/merge)
* [Plotting](Documentation/plot.md)

-----

## Authors

  * **Enzo Battistella** 
  * **Chris Rohlicek** 
  * **Tony Cui** 
  * **Bert Huang** 
  * **Yueyao Gao** 
  * **Anant Maheshwari** 
  * **Victoria Popic** 
