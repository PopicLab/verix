This demo provides several examples of how to run Verix on a small simulated dataset 
and shows the expected results for each. The dataset (`demo/sim.vcf`) is a simulated SV 
truthset produced by insilicoSV (using the `demo/case_study.yml` config), containing 1,000 events drawn from five SV classes 
(DEL, DUP_INV, dDUP, INVdel, delINVdel), with 200 instances per class.

### Test 1: Benchmarking with multi-record linking (perfect-match self-test)
In this test we simply compare the truthset against itself using the multi-record 
parsing mode, where records sharing a common SVID are 
aggregated into a single CSV (matching the insilicoSV output format). 

Every record should be matched to itself, such that all 1,000 CSVs are recovered 
as complete matches with zero breakpoint distance, yielding perfect precision and recall. 

```
verix bench \
    -q demo/sim.vcf \
    -t demo/sim.vcf \
    -o demo/test_bench_multi \
    -f multi multi \
    -l SVID SVID \
    --plot
```

The contents of `demo/test_bench_multi/report.json` should match:
```json{
    "n_query": 1000,
    "n_target": 1000,
    "tp_query": 1000,
    "tp_target": 1000,
    "fp": 0,
    "fn": 0,
    "precision": 1.0,
    "recall": 1.0,
    "f1": 1.0,
    "class_proportions": {
        "complete": 1.0,
        "partial": 0.0,
        "aggregate": 0.0,
        "spurious": 0.0
    }
    "by_class": {
        "complete": {
            "num_matches": 1000,
            "num_contiguous_matches": 1000,
            "num_unique_targets": 1000,
            "mean_breakpoint_distance": 0.0,
            "mean_breakpoint_hit_rate": 1.0,
            "mean_spurious_breakpoint_rate": 0.0,
            "mean_targets_per_record": 1.0,
            "num_matches_by_query_target_type": {
                "dDUP/dDUP": 200,
                "DEL/DEL": 200,
                "INVdel/INVdel": 200,
                "delINVdel/delINVdel": 200,
                "DUP_INV/DUP_INV": 200
            },
            "num_contiguous_matches_by_query_target_type": {
                "dDUP/dDUP": 200,
                "DEL/DEL": 200,
                "INVdel/INVdel": 200,
                "delINVdel/delINVdel": 200,
                "DUP_INV/DUP_INV": 200
            },
            "query_type_counts": {
                "dDUP": 200,
                "DEL": 200,
                "INVdel": 200,
                "delINVdel": 200,
                "DUP_INV": 200
            },
            "query_type_proportions": {
                "dDUP": 1.0,
                "DEL": 1.0,
                "INVdel": 1.0,
                "delINVdel": 1.0,
                "DUP_INV": 1.0
            },
            "target_type_counts": {
                "INVdel": 200,
                "dDUP": 200,
                "DUP_INV": 200,
                "DEL": 200,
                "delINVdel": 200
            },
            "target_type_proportions": {
                "INVdel": 1.0,
                "dDUP": 1.0,
                "DUP_INV": 1.0,
                "DEL": 1.0,
                "delINVdel": 1.0
            }
        }
    }
}
```

### Test 2: Benchmarking with mixed parsing modes (aggregation demo)
The same truthset can be parsed with different conventions on each side to illustrate 
how CSV assembly affects evaluation. Here the query is parsed in multi-record mode 
(CSVs aggregated by SVID), while the target is parsed using default mode 
(pairs of records linked via the standard POS/END/BND notation and TARGET for single dispersion loci). 

In this scenario, DELs, dDUPs, and DUP_INVs are assembled identically under both modes, 
but INVdel and delINVdel are split into multiple smaller target events. 
The asymmetric parsing produces 1,000 query CSVs and 1,800 target events. 
Verix classifies the query CSVs whose breakpoints span multiple target events as aggregate matches.
The resulting 400 aggregate matches correspond to INVdel and delINVdel query events which match multiple corresponding 
target events. 

```
    verix bench \
    -q demo/sim.vcf \
    -t demo/sim.vcf \
    -o demo/test_bench_mixed \
    -f multi default \
    -l SVID default \
    --plot
 ```

The contents of `demo/test_bench_mixed/report.json` should match:
```json{
    "n_query": 1000,
    "n_target": 1800,
    "tp_query": 600,
    "tp_target": 600,
    "fp": 400,
    "fn": 1200,
    "precision": 0.600,
    "recall": 0.333,
    "f1": 0.429,
    "class_proportions": {
        "complete": 0.600,
        "partial": 0.000,
        "aggregate": 0.400,
        "spurious": 0.000
    }
    "by_class": {
        "complete": {
            "num_matches": 600,
            "num_contiguous_matches": 600,
            "num_unique_targets": 600,
            "mean_breakpoint_distance": 0.0,
            "mean_breakpoint_hit_rate": 1.0,
            "mean_spurious_breakpoint_rate": 0.0,
            "mean_targets_per_record": 1.333333,
            "num_matches_by_query_target_type": {
                "dDUP/dDUP": 200,
                "DEL/DEL": 200,
                "DUP_INV/DUP_INV": 200
            },
            "num_contiguous_matches_by_query_target_type": {
                "dDUP/dDUP": 200,
                "DEL/DEL": 200,
                "DUP_INV/DUP_INV": 200
            },
            "query_type_counts": {
                "dDUP": 200,
                "DEL": 200,
                "DUP_INV": 200
            },
            "query_type_proportions": {
                "dDUP": 1.0,
                "DEL": 1.0,
                "DUP_INV": 1.0
            },
            "target_type_counts": {
                "dDUP": 200,
                "DUP_INV": 200,
                "DEL": 200
            },
            "target_type_proportions": {
                "dDUP": 1.0,
                "DUP_INV": 0.5,
                "DEL": 1.0
            }
        },
        "aggregate": {
            "num_matches": 400,
            "num_contiguous_matches": 400,
            "num_unique_targets": 400,
            "mean_breakpoint_distance": 0.0,
            "mean_breakpoint_hit_rate": 1.0,
            "mean_spurious_breakpoint_rate": 0.0,
            "mean_targets_per_record": 2.5,
            "num_matches_by_query_target_type": {
                "INVdel/INVdel": 200,
                "delINVdel/delINVdel": 200
            },
            "num_contiguous_matches_by_query_target_type": {
                "INVdel/INVdel": 200,
                "delINVdel/delINVdel": 200
            },
            "query_type_counts": {
                "INVdel": 200,
                "delINVdel": 200
            },
            "query_type_proportions": {
                "INVdel": 1.0,
                "delINVdel": 1.0
            },
            "target_type_counts": {
                "INVdel": 200,
                "delINVdel": 200
            },
            "target_type_proportions": {
                "INVdel": 0.5,
                "delINVdel": 0.3333333333333333
            }
        }
    }

}
```

### Test 3: Harmonization (self-test)

A similar sanity check for the harmonization command merges two copies of `demo/sim.vcf`. 
Each input record should be paired with its counterpart from the second copy, 
producing 1,000 consensus clusters of size 2.

```
    verix consensus \
    -i demo/sim.vcf demo/sim.vcf \
    -f multi multi \
    -l SVID SVID \
    -n sample1 sample2 \
    -o demo/test_merge
```

The contents of `demo/test_merge/report.json` should match:
```json{
    "n_total_variants": 2000,
    "n_variants_in_sample": {
        "sample1": 1000,
        "sample2": 1000
    },
    "n_clusters": 1000,
    "support_vec_types": [
        "1,1"
    ],
    "n_singleton_clusters": 0,
    "n_multi_caller_clusters": 1000,
    "n_full_support_clusters": 1000,
    "support_vec_counts": {
        "1,1": 1000
    },
    "cluster_size_counts": {
        "2": 1000
    },
    "max_cluster_size": 2,
    "min_cluster_size": 2
}
```

