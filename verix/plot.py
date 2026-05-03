from plotnine import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import warnings
warnings.filterwarnings("ignore", message="divide by zero encountered in log10", category=RuntimeWarning)


class BenchPlotter:
    MATCH_COLORS = {
        "complete":  '#1B5E3F',  # dark green
        "partial":   '#B8956A',  # warm tan
        "aggregate": '#B5A8D1',  # soft lavender
        "spurious":  '#6B1F2E',  # dark wine
    }
    TARGET_COLORS = {
        'full': '#5A8F5A',  # muted sage green
        'partial': '#D4A84A',  # soft amber
        'miss': '#A04545',  # muted brick red
    }

    def __init__(self, out_dir, match_df, target_df):
        self.out_dir = out_dir
        self.match_df = match_df
        self.target_df = target_df
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
                       x="Number of records", y="Number of target breakpoints matched")
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
                + self.base_theme(figure_size=(self.scale_dim(self.n_match_types), self.scale_dim(self.n_query_types)))
                + geom_bar(aes(fill="..x.."), color="white", show_legend=False)
                + scale_fill_gradientn(colors=["#1B7837", "#A50F15"])
                + scale_y_log10(minor_breaks=[])
                + scale_x_continuous(breaks=range(0, int(df["NTARGETS"].max()) + 1), minor_breaks=[])
                + facet_grid("QTYPE~BEST_MATCH_CLASS")
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
                  .size().reset_index(name="n"))
        all_types = sorted(set(self.match_df["QTYPE"].dropna().unique())
                           | set(self.match_df["BEST_MATCH_TYPE"].dropna().unique()))
        full_index = pd.MultiIndex.from_product(
            [counts["BEST_MATCH_CLASS"].unique(), all_types, all_types],
            names=["BEST_MATCH_CLASS", "QTYPE", "BEST_MATCH_TYPE"])
        counts = (counts.set_index(["BEST_MATCH_CLASS", "QTYPE", "BEST_MATCH_TYPE"])
                  .reindex(full_index, fill_value=0).reset_index())
        n_facets = counts["BEST_MATCH_CLASS"].nunique()
        return (
                ggplot(counts, aes(x="BEST_MATCH_TYPE", y="QTYPE", fill="n"))
                + theme_minimal()
                + theme(figure_size=(max(6.0, 0.4 * len(all_types) + 2.0),
                                     max(4, 0.4 * len(all_types) * n_facets + 2)),
                        axis_text_x=element_text(rotation=45, ha="right"),
                        strip_text=element_text(size=14),
                        panel_border=element_rect(color="black", size=1, fill=None),
                        panel_grid_major=element_blank(),
                        panel_grid_minor=element_blank(),
                        axis_ticks=element_blank())
                + facet_wrap("~BEST_MATCH_CLASS", ncol=1, scales="free_x")
                + geom_tile(color="white", size=0.4)
                + geom_text(aes(label="n"), data=counts[counts["n"] > 0], size=8)
                + scale_fill_gradient(low="#f7f7f7", high="#1f77b4", name="Count")
                + labs(title="Query versus target SV type correspondence for optimal matches",
                       x="Target SV type", y="Query SV type"))

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
            for name, p in plots.items():
                fig = p.draw()
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)
