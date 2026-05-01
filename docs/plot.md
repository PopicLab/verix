# Plot

---

## Plot Types Reference

* **`sankey`**: Generates a Plotly Sankey diagram showing the flow of SVs from Ground Truth pools (Simple vs. Complex) into their respective detection modes (Complete, Aggregate, Partial, Miss).

<div style="text-align: center;">
<img src="figures/sankey.svg" width="500" alt="Sankey">
</div>

* **`stacked_barplot`**: Displays log-scaled counts per called SV type, stacked and color-coded by detection mode (Complete, Partial, Aggregate, Miss).

<div style="text-align: center;">
<img src="figures/stacked_barplot.svg" width="500" alt="barplot">
</div>

* **`breakends_dist`**: Generates boxplots showing the distribution of matched called and GT breakends per SV.

<div style="text-align: center;">
<img src="figures/breakends_distribution.svg" width="500" alt="breakends">
</div>

* **`fragmentation`**: Shows the total number of fragmented predictions, broken down by the GT SV type they attempted to cover.

<div style="text-align: center;">
<img src="figures/fragmented_calls.svg" width="500" alt="frag">
</div>

