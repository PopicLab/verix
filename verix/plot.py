import textwrap
from plotnine import *
import pandas as pd
import plotly.graph_objects as go
import matplotlib.colors as mcolors
import logging
import os
import numpy as np
from mizani.transforms import pseudo_log_trans

NORMAL_FONT = 11
AXIS_FONT = 13
TITLE_FONT = 17


def process_counts(counts, is_gt=True):
    records = []
    for t1, modes in counts.items():
        for mode, pred_types in modes.items():
            for t2, val in pred_types.items():
                gt_type = t1 if is_gt else t2
                p_type = t2 if is_gt else t1
                records.append([gt_type, p_type, val, mode])
    return records


class Plotter:
    COLORS = {
        'complete': '#00A138',  # Green
        'partial': '#FFD700',  # Gold
        'aggregate': '#FF7800',  # Orange
        'miss': '#BF001E',  # Red
        'grey': '#787878',  # Neutral for GT/Background
        'FP': '#BF001E',  # Alias for Red
        'FN': '#BF001E'  # Alias for Red
    }

    def __init__(self, results_dict, output):
        self.plot_function_dict = {'sankey': self.plot_sankey,
                                   'breakends_dist': self.plot_breakends_distribution,
                                   'stacked_barplot': self.plot_detection_barplot,
                                   'fragmentation': self.plot_fragmentation_barplot
                                   }

        self.output = output

        self.truth_df = results_dict['truth_df']
        self.pred_df = results_dict['pred_df']
        self.optimal_df = results_dict['optimal_df']

        self.optimal_counts_gt = pd.DataFrame(
            process_counts(results_dict['count_optimal_gt'], is_gt=True),
            columns=['gt_type', 'predicted_sv_type', 'counts', 'detection_mode'])
        self.optimal_counts_calls = pd.DataFrame(
            process_counts(results_dict['count_optimal_calls'], is_gt=True),
            columns=['predicted_sv_type', 'gt_type', 'counts', 'detection_mode'])

        self.plot_folder = self.output + '/plots/'
        os.makedirs(self.plot_folder, exist_ok=True)

    def plot(self):
        for plottype, plot_function in self.plot_function_dict.items():
            logging.info(f'Plotting {plottype} scatter')
            plot_function()

    def plot_sankey(self):
        def hex_to_rgba(hex_code, alpha):
            rgb = mcolors.to_rgb(hex_code)
            return f'rgba({int(rgb[0] * 255)}, {int(rgb[1] * 255)}, {int(rgb[2] * 255)}, {alpha})'

        df = self.optimal_counts_gt.copy()
        modes = ['complete', 'aggregate', 'partial', 'miss']
        masks = {m: df['detection_mode'].str.contains(m, na=False, regex=False) for m in modes}

        compare_col = 'predicted_sv_type'
        base_col = 'gt_type'

        fp_mask = (df[compare_col] == 'miss')
        simple_list = ['INV', 'DUP', 'DEL', 'INS', 'DUP:TANDEM']
        simple_mask = df[base_col].isin(simple_list) & ~fp_mask
        complex_mask = ~df[base_col].isin(simple_list) & ~fp_mask

        gt_label = 'Ground Truth'

        def get_w(type_mask, mode_key):
            return df[type_mask & masks[mode_key]]['counts'].sum()

        weights = {}
        for prefix, m_mask in [('SMPL', simple_mask), ('CPLX', complex_mask)]:
            weights[f'gt->{prefix}'] = df[m_mask]['counts'].sum()
            for m in modes:
                weights[f'{prefix}->{m}'] = get_w(m_mask, m)

        cat_counts = {
            gt_label: df[~fp_mask]['counts'].sum(),
            'Simple': weights['gt->SMPL'],
            'Complex': weights['gt->CPLX'],
            'Complete': weights['SMPL->complete'] + weights['CPLX->complete'],
            'Aggregate': weights['SMPL->aggregate'] + weights['CPLX->aggregate'],
            'Partial': weights['SMPL->partial'] + weights['CPLX->partial'],
            'Miss': weights['SMPL->miss'] + weights['CPLX->miss']
        }

        nodes = [gt_label, "Simple", "Complex", "Complete", "Aggregate", "Partial", "Miss"]
        active_nodes = [n for n in nodes if cat_counts.get(n, 0) > 0]
        ind_map = {name: i for i, name in enumerate(active_nodes)}

        cat_colors = {
            gt_label: hex_to_rgba(self.COLORS['grey'], 0.8),
            'Simple': hex_to_rgba(self.COLORS['grey'], 0.8),
            'Complex': hex_to_rgba(self.COLORS['grey'], 0.8),
            'Complete': hex_to_rgba(self.COLORS['complete'], 0.8),
            'Aggregate': hex_to_rgba(self.COLORS['aggregate'], 0.8),
            'Partial': hex_to_rgba(self.COLORS['partial'], 0.8),
            'Miss': hex_to_rgba(self.COLORS['miss'], 0.8)
        }

        raw_conns = [
            (gt_label, 'Simple', weights['gt->SMPL'], hex_to_rgba(self.COLORS['grey'], 0.3)),
            (gt_label, 'Complex', weights['gt->CPLX'], hex_to_rgba(self.COLORS['grey'], 0.3)),
        ]

        mode_to_node = {'complete': 'Complete', 'aggregate': 'Aggregate', 'partial': 'Partial', 'miss': 'Miss'}
        for prefix in ['SMPL', 'CPLX']:
            label = 'Simple' if prefix == 'SMPL' else 'Complex'
            for m_key, target_node in mode_to_node.items():
                w_key = f'{prefix}->{m_key}'
                raw_conns.append((label, target_node, weights[w_key], hex_to_rgba(self.COLORS[m_key], 0.4)))

        sources, targets, values, colors = [], [], [], []
        for s, t, v, c in raw_conns:
            if v > 0 and s in ind_map and t in ind_map:
                sources.append(ind_map[s])
                targets.append(ind_map[t])
                values.append(v)
                colors.append(c)

        labels = [f"{n}<br>({int(cat_counts[n])})" for n in active_nodes]
        if 'Complete' in ind_map and cat_counts[gt_label] > 0:
            metric_val = cat_counts['Complete'] / cat_counts[gt_label]
            metric_name = "Recall"
            labels[ind_map['Complete']] += f"<br>{metric_name}: {metric_val:.1%}"

        fig = go.Figure(data=[go.Sankey(
            node=dict(pad=20, thickness=20, line=dict(color="black", width=0.5), label=labels, color=[cat_colors[n] for n in active_nodes]),
            link=dict(source=sources, target=targets, value=values, color=colors)
        )])

        fig.update_layout(
            title=dict(text=f"Sankey Diagram", x=0.5, font=dict(size=TITLE_FONT, color="black")),
            font=dict(size=NORMAL_FONT),
            paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)'
        )
        self.save_plot(fig, 'sankey', is_plotly=True)

    def plot_breakends_distribution(self):
        """
        Plots the distribution of Total Predicted Breakends, Spurious Breakends,
        and Missing (Unmatched) Ground Truth Breakends per SV.
        """
        figsize = (8, 6)

        data_map = [
            (self.truth_df, 'num_breakends', 'GT Breakends (All)'),
            (self.truth_df, 'spurious', 'Missed GT Breakends'),
            (self.truth_df.query("detection_mode != 'complete'"), 'num_breakends', 'GT Breakends (Incomplete)'),
            (self.pred_df, 'num_breakends', 'Call Breakends (All)'),
            (self.pred_df, 'spurious', 'Spurious Call Breakends'),
            (self.pred_df.query("detection_mode != 'complete'"), 'num_breakends', 'Call Breakends (Incomplete)')
        ]

        plot_list = []
        for df, col, label in data_map:
            if not df.empty:
                temp_df = pd.DataFrame({'value': df[col], 'Metric': label})
                plot_list.append(temp_df)

        if not plot_list:
            return None

        plot_df = pd.concat(plot_list).dropna(subset=['value'])
        metric_order = [m[2] for m in data_map]
        plot_df['Metric'] = pd.Categorical(plot_df['Metric'], categories=metric_order, ordered=True)

        p = (ggplot(plot_df, aes(x='Metric', y='value', fill='Metric'))
             + geom_boxplot(outlier_size=0.5, outlier_alpha=0.5, alpha=0.8)
             + scale_y_continuous(trans=pseudo_log_trans(base=10))
             + labs(title=f'Breakend Distributions per SV', x='', y='Number Breakends (pseudo-log scale)')
             + theme_minimal()
             + theme(figure_size=figsize, legend_position='none',
                     axis_text_x=element_text(rotation=25, hjust=1, size=NORMAL_FONT, color="black"),
                     axis_text_y=element_text(size=NORMAL_FONT, color="black"),
                     plot_title=element_text(size=TITLE_FONT, weight='bold', ha='center'),
                     panel_border=element_rect(color="black", size=1))
             )

        self.save_plot(p, 'breakends_distribution')
        return p

    def plot_detection_barplot(self):
        """
        Plot a stacked barplot with the type of detection for each ground truth type
        """
        figsize = (6, 8)
        df_plot = self.optimal_counts_calls.copy()

        df_plot = df_plot[df_plot['predicted_sv_type'] != 'NA']
        if df_plot.empty: return

        df_plot['predicted_sv_type'] = df_plot['predicted_sv_type'].apply(lambda x: textwrap.fill(str(x), width=35))
        df_plot['detection_mode'] = pd.Categorical(df_plot['detection_mode'],
                                                   categories=['complete', 'aggregate', 'partial', 'miss'],
                                                   ordered=True)

        df_plot = df_plot.groupby(['predicted_sv_type', 'detection_mode'], observed=True)['counts'].sum().reset_index()
        df_totals = df_plot.groupby(['predicted_sv_type'], observed=True)['counts'].sum().reset_index()

        df_plot = df_plot.merge(df_totals.rename(columns={'counts': 'total_counts'}), on=['predicted_sv_type'],
                                how='left')
        #Manually compute log scale to allow correct propotions for the different stacks
        df_plot['log_total'] = np.log1p(df_plot['total_counts'])
        df_plot['proportion'] = df_plot['counts'] / df_plot['total_counts']
        # Ensures the proportions of the stacked bar are correct despite logscale
        df_plot['pseudo_count'] = df_plot['proportion'] * df_plot['log_total']

        df_totals['log_total'] = np.log1p(df_totals['counts'])
        df_plot['inner_label'] = df_plot.apply(lambda row: f"{row['proportion'] * 100:.0f}%" if (
                    row['counts'] > 0 and row['proportion'] >= 0.05 and row['proportion'] != 1.) else "", axis=1)

        log_breaks_raw = [0, 10, 100, 1000, 10000, 100000]
        log_breaks_mapped = [np.log1p(v) for v in log_breaks_raw]

        p = (ggplot(df_plot, aes(x='predicted_sv_type', y='pseudo_count', fill='detection_mode'))
             + geom_col(width=0.6, position='stack')
             + scale_fill_manual(values=self.COLORS, name="Detection Mode")
             + scale_y_continuous(breaks=log_breaks_mapped, labels=['0', '10', '100', '1,000', '10,000', '100,000'],
                                  expand=(0, 0, 0.15, 0.02))
             + geom_text(aes(label='inner_label'), position=position_stack(vjust=0.5),
                         size=NORMAL_FONT - 5, color="black")
             + geom_text(data=df_totals, mapping=aes(label='counts', x='predicted_sv_type', y='log_total'),
                         inherit_aes=False, angle=45, va='bottom', ha='center', format_string="{:,.0f}",
                         size=NORMAL_FONT - 1)
             + labs(title=f'Counts per Detection Mode', x='Predicted SV type', y='SV count')
             + theme_minimal()
             + theme(figure_size=figsize, legend_position='top', axis_line=element_line(size=1, color="black"),
                     panel_grid_major_x=element_blank(),
                     axis_text_x=element_text(size=NORMAL_FONT, rotation=70, hjust=1, color="black"),
                     axis_text_y=element_text(size=NORMAL_FONT, color="black"),
                     axis_title=element_text(size=AXIS_FONT, weight="bold"),
                     plot_title=element_text(size=TITLE_FONT, weight='bold', ha='center'))
             )

        self.save_plot(p, 'stacked_barplot')
        return p

    def plot_fragmentation_barplot(self):
        """
        Plots a bar chart showing the total count of fragmented predictions
        broken down by the Ground Truth SV Type they attempted to cover.
        """
        figsize = (6, 5)
        df = self.optimal_df.copy()
        df = df[(df['fragmented'] == True) & (df['is_optimal_GT'] == True)]
        if df.empty: return

        counts = df.groupby(['GT type']).size().reset_index(name='fragment_count')

        p = (ggplot(counts, aes(x='GT type', y='fragment_count', fill='GT type'))
             + geom_col(width=0.6, color="black", size=0.2)
             + geom_text(aes(label='fragment_count'), va='bottom', size=NORMAL_FONT - 1)
             + scale_y_continuous(expand=(0, 0, 0.15, 0))
             + labs(title='Fragmented Calls per GT SV Type', x='Ground Truth SV Type',
                    y='Fragmented Calls')
             + theme_minimal()
             + theme(figure_size=figsize, legend_position='none',
                     axis_text_x=element_text(rotation=35, hjust=1, size=NORMAL_FONT, color="black"),
                     axis_text_y=element_text(size=NORMAL_FONT, color="black"),
                     plot_title=element_text(size=TITLE_FONT, weight='bold', ha='center'),
                     panel_grid_major_x=element_blank())
             )

        self.save_plot(p, 'fragmented_calls')
        return p

    def save_plot(self, plot_obj, name, is_plotly=False):
        filename = f"{name}.svg"
        full_path = os.path.join(self.plot_folder, filename)

        if is_plotly:
            plot_obj.write_image(full_path, format='svg')
        else:
            plot_obj.save(full_path, format='svg', verbose=False, dpi=500, limitsize=False)

        logging.info(f"Saved: {full_path}")
