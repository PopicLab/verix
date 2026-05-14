from plotnine import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import warnings
warnings.filterwarnings("ignore", message="divide by zero encountered in log10", category=RuntimeWarning)
plt.rcParams["svg.fonttype"] = "path"

class BenchPlotter:
    MATCH_COLORS = {
        "complete":  '#3E8F6E',
        "partial":   '#B8956A',
        "aggregate": '#B5A8D1',
        "spurious":  '#6B1F2E',
    }
    TARGET_COLORS = {
        'full': '#5A8F5A',
        'partial': '#D4A84A',
        'miss': '#A04545',
    }

    def __init__(self, out_dir, match_df, target_df, stats):
        self.out_dir = out_dir
        self.match_df = match_df
        self.target_df = target_df
        self.stats = stats
        self.n_query_types = self.match_df["QTYPE"].nunique()
        self.n_target_types = self.target_df["TYPE"].nunique()
        self.n_match_types =  self.match_df["BEST_MATCH_TYPE"].nunique()

    @staticmethod
    def base_theme(figure_size):
        return (
            theme_minimal() +
            theme(axis_line=element_line(size=2, color="black"),
            text=element_text(size=12, color="black"),
            axis_title=element_text(size=12, color="black"),
            plot_title=element_text(size=12, color="black"),
            strip_text=element_text(size=10, color="black"),
            legend_title=element_blank(),
            panel_grid_major_x=element_blank(),
            panel_grid_minor_y=element_blank(),
            panel_grid_major_y=element_line(size=1, color="grey", linetype="dotted"),
            figure_size=figure_size)
        )

    @staticmethod
    def scale_dim(n_categories, facets=1, extra=0):
        return max(6, 0.8 * n_categories * facets + 3 + extra)

    # how many target records were fully missed or fully/partially covered by an optimal query alignment
    def plot_target_capture_by_type(self, fraction=False):
        return (
                ggplot(self.target_df, aes(x="TYPE", fill="COV"))
                + self.base_theme(figure_size=(self.scale_dim(self.n_target_types), 5))
                + geom_bar(position=("fill" if fraction else "stack"))
                + scale_fill_manual(values=self.TARGET_COLORS)
                + labs(title="Target breakpoint capture across optimal query alignments",
                       x="Target SV type", y=("Fraction" if fraction else "Count"))
                + theme(axis_text_x=element_text(rotation=30, ha="right")))

    # how many target records were fully missed or fully/partially covered by the union of all query matches
    def plot_target_capture_by_union(self, fraction=False):
        return (
                ggplot(self.target_df, aes(x="TYPE", fill="UNION_COV"))
                + self.base_theme(figure_size=(self.scale_dim(self.n_target_types), 5))
                + geom_bar(position=("fill" if fraction else "stack"))
                + scale_fill_manual(values=self.TARGET_COLORS)
                + labs(title="Target breakpoint capture across all candidate query alignments",
                       x="Target SV type", y="Number of targets")
                + theme(axis_text_x=element_text(rotation=30, ha="right")))

    # how many target breakpoints were matched in optimal alignments, stratified by match category and target type
    def plot_target_breakpoint_hit_rate(self):
        df = self.match_df.dropna(subset=["BEST_N_MATCHED", "BEST_MATCH_TYPE"]).copy()
        df = df[df["BEST_MATCH_CLASS"] != "spurious"]
        return (
                ggplot(df, aes(x="BEST_N_MATCHED"))
                + self.base_theme(figure_size=(self.scale_dim(self.n_match_types), self.scale_dim(self.n_target_types)))
                + geom_bar(aes(fill="..x..", group=1), color="white", show_legend=False)
                + scale_fill_gradient(low="#F7FCB9", high="#1B5E3F")
                + scale_y_log10(minor_breaks=[])
                + scale_x_continuous(breaks=range(0, int(df["BEST_N_MATCHED"].max()) + 1), minor_breaks=[])
                + facet_grid("BEST_MATCH_TYPE~BEST_MATCH_CLASS")
                + labs(title="Target breakpoints matched by target SV type",
                       x="Number of target breakpoints matched", y="Count")
                + theme(axis_text_x=element_text(rotation=0, ha="right")))

    # how many query events are in each match category, stratified by query type
    def plot_match_category_by_type(self):
        return (
                ggplot(self.match_df, aes(x="QTYPE", fill="BEST_MATCH_CLASS"))
                + self.base_theme(figure_size=(self.scale_dim(self.n_query_types), 5))
                + geom_bar(position=position_dodge2(preserve="single"))
                + scale_y_log10(minor_breaks=[])
                + scale_fill_manual(values=self.MATCH_COLORS)
                + labs(title="Match category split by query SV type", x="Query SV type", y="Count")
                + theme(axis_text_x=element_text(rotation=30, ha="right")))

    # how many different targets were matched by the query (on at least 1 breakpoint)
    def plot_targets_per_query(self):
        df = self.match_df.copy()
        df = df[df["BEST_MATCH_CLASS"] != "spurious"]
        return (
                ggplot(df, aes(x="NTARGETS"))
                + self.base_theme(figure_size=(self.scale_dim(self.n_match_types, extra=3),
                                               self.scale_dim(self.n_query_types)))
                + geom_bar(aes(fill="..x.."), color="white", show_legend=False)
                + scale_fill_gradientn(colors=["#1B7837", "#A50F15"])
                + scale_y_log10(minor_breaks=[])
                + scale_x_continuous(breaks=range(0, int(df["NTARGETS"].max()) + 1), minor_breaks=[])
                + facet_grid("QTYPE~BEST_MATCH_CLASS")
                + theme(axis_text_x=element_text(rotation=90, ha="right", size=8))
                + labs(title="Distinct targets matched per query record",
                       x="Number of distinct targets in candidate alignments",
                       y="Number of query records"))

    # total breakpoint distance to optimal match
    def plot_optimal_distance_histogram(self):
        df = self.match_df.dropna(subset=["BEST_BND_DIST"]).copy()
        df = df[df["BEST_MATCH_CLASS"] != "spurious"]
        return (
                ggplot(df, aes(x="BEST_BND_DIST"))
                + self.base_theme(figure_size=(self.scale_dim(self.n_match_types), self.scale_dim(self.n_query_types)))
                + geom_histogram(aes(fill="..x.."), bins=10, color="white", show_legend=False)
                + scale_fill_gradient(low="#FEE5D9", high="#A50F15")
                + facet_grid("QTYPE~BEST_MATCH_CLASS", scales="free")
                + scale_y_log10(minor_breaks=[])
                + labs(title="Optimal alignment breakpoint distance",
                       x="Sum of distances across all matched breakpoints", y="Count")
                + theme(axis_text_x=element_text(rotation=30, ha="right")))

    # number of query records with at least 1 spurious breakpoint, stratified by query type
    def plot_records_with_spurious(self):
        df = self.match_df[self.match_df["SPURIOUS"] >= 1]
        return (
                ggplot(df, aes(x="QTYPE", fill="BEST_MATCH_CLASS"))
                + self.base_theme(figure_size=(self.scale_dim(self.n_query_types), 5))
                + geom_bar(position=position_dodge2(preserve="single"))
                + scale_y_log10(minor_breaks=[])
                + scale_fill_manual(values=self.MATCH_COLORS)
                + labs(title="Records with at least one spurious breakpoint",
                       x="Query SV type", y="Number of records")
                + theme(axis_text_x=element_text(rotation=30, ha="right")))

    # number of query records that are fragmented, stratified by query type
    def plot_fragmented_records(self):
        df = self.match_df[self.match_df["FRAGMENTED"] >= 1]
        return (
                ggplot(df, aes(x="QTYPE", fill="BEST_MATCH_CLASS"))
                + self.base_theme(figure_size=(self.scale_dim(self.n_query_types), 5))
                + geom_bar(position=position_dodge2(preserve="single"))
                + scale_y_log10(minor_breaks=[])
                + scale_fill_manual(values=self.MATCH_COLORS)
                + labs(title="Fragmented records",
                       x="Query SV type", y="Number of records")
                + theme(axis_text_x=element_text(rotation=30, ha="right")))

    def plot_optimal_type_correspondence(self):
        df = self.match_df.dropna(subset=["BEST_MATCH_TYPE"])
        counts = (df.groupby(["BEST_MATCH_CLASS", "QTYPE", "BEST_MATCH_TYPE"], observed=True)
                  .agg(n=("BEST_IS_CONTIGUOUS", "size"),
                       n_contig=("BEST_IS_CONTIGUOUS", "sum"))
                  .reset_index())
        query_types = sorted(self.match_df["QTYPE"].dropna().unique())
        target_types = sorted(self.match_df["BEST_MATCH_TYPE"].dropna().unique())
        full_index = pd.MultiIndex.from_product(
            [counts["BEST_MATCH_CLASS"].unique(), query_types, target_types],
            names=["BEST_MATCH_CLASS", "QTYPE", "BEST_MATCH_TYPE"])
        counts = (counts.set_index(["BEST_MATCH_CLASS", "QTYPE", "BEST_MATCH_TYPE"])
                  .reindex(full_index, fill_value=0).reset_index())
        counts["n_contig"] = counts["n_contig"].astype(int)
        counts["label"] = counts.apply(
            lambda r: f"{r['n']}\n({int(r['n_contig'])})" if r["n_contig"] > 0 else f"{r['n']}",
            axis=1)
        n_facets = counts["BEST_MATCH_CLASS"].nunique()
        return (
                ggplot(counts, aes(x="BEST_MATCH_TYPE", y="QTYPE", fill="n"))
                + theme_minimal()
                + theme(figure_size=(max(6.0, 0.6 * len(target_types) + 2),
                                     max(6.0, 0.6 * len(query_types) * n_facets + 2)),
                        axis_text_x=element_text(rotation=45, ha="right", size=14),
                        axis_text_y=element_text(rotation=0, size=14),
                        strip_text=element_text(size=24),
                        panel_border=element_rect(color="black", size=1, fill=None),
                        panel_grid_major=element_blank(),
                        panel_grid_minor=element_blank(),
                        axis_ticks=element_blank())
                + facet_wrap("~BEST_MATCH_CLASS", ncol=1, scales="free_x")
                + geom_tile(color="white", size=0.4)
                + geom_text(aes(label="label"), data=counts[counts["n"] > 0], size=14)
                + scale_fill_gradient(low="#f7f7f7", high="#1f77b4", name="Count")
                + labs(title="Predicted versus true SV type correspondence for optimal matches",
                       x="True SV type", y="Predicted SV type"))

    def plot_summary(self):
        s = self.stats
        classes = [c for c in ("complete", "partial", "aggregate") if c in s.get("by_class", {})]
        n_classes = len(classes)
        fig = plt.figure(figsize=(5 * max(n_classes, 1) + 2, 12))
        gs = fig.add_gridspec(4, n_classes, height_ratios=[0.5, 0.4, 1, 1], hspace=0.5, wspace=0.4)
        count_fields = ["tp_query", "tp_target", "fp", "fn"]
        rate_fields = ["precision", "recall", "f1"]
        fmt = lambda t: f"{t}: {s[t]:.3f}" if isinstance(s[t], float) else f"{t}: {s[t]}"
        counts_line = " , ".join(fmt(t) for t in count_fields if t in s)
        rates_line = " , ".join(fmt(t) for t in rate_fields if t in s)
        headline_lines = [line for line in (counts_line, rates_line) if line]
        ax_head = fig.add_subplot(gs[0, :])
        ax_head.axis("off")
        ax_head.set_title(f"Benchmark summary: n_query={s['n_query']}, n_target={s['n_target']}",
                          fontsize=24, pad=12)
        ax_head.text(0.5, 0.3, "\n\n".join(headline_lines), ha="center", va="center", fontsize=20,
                     color="#34495E", transform=ax_head.transAxes)

        for col, cls in enumerate(classes):
            cs = s["by_class"][cls]
            color = self.MATCH_COLORS[cls]
            ax_label = fig.add_subplot(gs[1, col])
            ax_label.axis("off")
            ax_label.text(0.5, 0.9, cls.upper(), ha="center", va="top", fontsize=20,
                          color=color, transform=ax_label.transAxes)
            ax_label.text(0.5, 0.3,
                          f"num_matches: {cs['num_matches']}\nnum_unique_targets: {cs['num_unique_targets']}",
                          color="#34495E", ha="center", va="center", fontsize=20, transform=ax_label.transAxes)
            for row, key_count, key_frac in [(2, "query_type_counts", "query_type_proportions"),
                                             (3, "target_type_counts", "target_type_proportions")]:
                ax = fig.add_subplot(gs[row, col])
                items = list(cs[key_count].items())[::-1]
                labels = [k if len(k) <= 12 else k[:11] + "…" for k, _ in items]
                values = [v for _, v in items]
                bars = ax.barh(labels, values, color=color)
                ax.set_title(key_count, fontsize=20)
                ax.tick_params(labelsize=20)
                if values:
                    if max(values) > 1000:
                        ax.set_xscale("log")
                        ax.set_xlim(0.08, max(values) * 3)
                    else:
                        ax.set_xlim(0, max(values) * 1.2)
                for bar, v, frac in zip(bars, values, list(cs[key_frac].values())[::-1]):
                    bw = bar.get_width()
                    if bw > max(values) * 0.2:
                        ax.text(bw - max(values) * 0.01, bar.get_y() + bar.get_height() / 2,
                                f"{100 * frac:.1f}%", va="center", ha="right", fontsize=16, color="black")
                    else:
                        ax.text(bw + max(values) * 0.01, bar.get_y() + bar.get_height() / 2,
                                f"{100 * frac:.1f}%", va="center", ha="left", fontsize=16, color="black")
        return fig

    def make_plots(self):
        plots = {
            "optimal_match_category_by_type": self.plot_match_category_by_type(),
            "records_with_spurious": self.plot_records_with_spurious(),
            "records_fragmented": self.plot_fragmented_records(),
            "optimal_distance": self.plot_optimal_distance_histogram(),
            "targets_per_query": self.plot_targets_per_query(),
            "target_capture_by_type": self.plot_target_capture_by_type(),
            "target_capture_by_union": self.plot_target_capture_by_union(),
            "target_breakpoint_hit_rate": self.plot_target_breakpoint_hit_rate(),
            "optimal_type_correspondence": self.plot_optimal_type_correspondence(),
        }
        with PdfPages(self.out_dir / f"report.pdf") as pdf:
            summary = self.plot_summary()
            pdf.savefig(summary, bbox_inches="tight")
            plt.close(summary)
            for name, p in plots.items():
                fig = p.draw()
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)
