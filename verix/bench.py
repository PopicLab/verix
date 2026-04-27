import json
import logging
from collections import defaultdict
from pathlib import Path

import pandas as pd

from verix.sv import MatchCallset, SVMatch
from verix.vcf_writer import write_bench_vcf, stringify_matched_bnds


class BenchCallset(MatchCallset):
    def __init__(self, svs, pred_callset=None, thresh=500,
                 unmatched_bnds_thresh=None, match_ratio_thresh=None, enforce_svtype=False, enforce_genotype=False):
        super().__init__(svs, thresh=thresh, unmatched_bnds_thresh=unmatched_bnds_thresh,
                         match_ratio_thresh=match_ratio_thresh, enforce_svtype=enforce_svtype, enforce_genotype=enforce_genotype)
        # dicts of the form {gt_pid: SV_Match}
        self.complete_detections = defaultdict(list)
        self.incomplete_detections = defaultdict(list)

        # Counts of each detection mode
        self.matched_calls = set()

        # Match the bnds of all call events to be used to find the match candidates and in the output information
        self.pred_callset = pred_callset
        self.pred_callset.get_all_bnds()
        for sv in self.pred_callset.svs:
            self.find_candidates(sv, self.pid2sv, self.all_bnds, self.all_bnds_pos)

        for sv in self.svs:
            self.find_candidates(sv, self.pred_callset.pid2sv, self.pred_callset.all_bnds, self.pred_callset.all_bnds_pos)

    def callset2dataframe(self, is_gt=False):
        return super().callset2dataframe(is_gt=True)

    @staticmethod
    def is_aggregate(sv, matched_bnds, candidate_idx):
        # Determine if a match is an aggregate
        for bnd_idx, bnd in enumerate(sv.breakends):
            if bnd in matched_bnds: continue
            if any(bnd_id != candidate_idx for bnd_id in bnd.bnd_matches):
                return True
        return False

    @staticmethod
    def is_contiguous(matched_bnds, bnd2idx, bnd2idx_candidate):
        # Determine if a match is contiguous
        list_bnd_idx = defaultdict(list)
        list_candidate_idx = defaultdict(list)

        for bnd, candidate_bnd in matched_bnds:
            # The BNDs have to match in order if on the same chrom
            bnd_idx = bnd2idx[bnd]
            candidate_bnd_idx = bnd2idx_candidate[candidate_bnd]
            list_bnd_idx[bnd.chrom].append(bnd_idx)
            list_candidate_idx[candidate_bnd.chrom].append(candidate_bnd_idx)

        for list_idx_chr in list(list_bnd_idx.values()) + list(list_candidate_idx.values()):
            # If one idx is missing between min and max, it is not consecutive
            if max(list_idx_chr) - min(list_idx_chr) != len(list_idx_chr) - 1:
                return False

        return True

    def find_detections(self, call_sv):
        total_matched_bnds = set()

        for gt_pid in call_sv.candidates:
            gt_sv = self.pid2sv[gt_pid]

            matched_bnds = self.align_candidates(call_sv, gt_sv)

            # Compute number of unmatched bnds
            total_num_bnds_call = len(call_sv.breakends)
            total_num_bnds_gt = len(gt_sv.breakends)
            diff_bnd = total_num_bnds_call + total_num_bnds_gt - 2 * len(matched_bnds)

            # Apply filters
            if self.filter_candidate(diff_bnd, len(matched_bnds), diff_bnd + len(matched_bnds)): continue

            total_distance = sum(abs(match[0].pos - match[1].pos) for match in matched_bnds)

            # To determine later if the SV is spurious
            call_matched_bnds = set(call_bnd for call_bnd, _ in matched_bnds)
            total_matched_bnds = total_matched_bnds.union(call_matched_bnds)

            # Identify the type of match
            if len(matched_bnds) == total_num_bnds_call == total_num_bnds_gt:
                detection_mode = 'complete'
            else:
                detection_mode = 'partial'

                if self.is_aggregate(call_sv, call_matched_bnds, gt_pid):
                    detection_mode = 'aggregate'

            contig = 'Contiguous' if self.is_contiguous(matched_bnds, call_sv.bnd2idx, gt_sv.bnd2idx) else 'Discontiguous'
            match = SVMatch(call_sv, gt_sv, detection_mode, contig, diff_bnd, total_distance, matched_bnds)

            gt_sv.match_list.append(match)
            call_sv.match_list.append(match)
            if detection_mode == 'complete':
                self.complete_detections[gt_pid].append(match)
            else:
                self.incomplete_detections[gt_pid].append(match)

            # update the optimal candidate for both the call and GT SVs
            for obj, other_sv in [(call_sv, gt_sv), (gt_sv, call_sv)]:
                score = (-diff_bnd, -total_distance)
                if (obj.optimal_candidate is None or (detection_mode == 'complete' and obj.detection_mode != 'complete')
                        or (score > obj.optimal_score)):
                    obj.optimal_candidate = match
                    obj.optimal_score = score
                    obj.detection_mode = detection_mode

        call_sv.spurious = len(call_sv.breakends) - len(total_matched_bnds)

    def collect_candidate_detections(self):
        # for a given comparison callset, we identify which SVs match with the various groundtruth SVs
        for sv in self.pred_callset.svs:
            self.find_detections(sv)

        # Find the number of spurious breakends for the ground truth and if the matches are fragmented
        for sv in self.svs:
            matched_bnds = set(bnds[1] for match in sv.match_list for bnds in match.matched_bnds)
            sv.spurious = len(sv.breakends) - len(matched_bnds)
            for match in sv.match_list:
                if len(match.matched_bnds) < len(matched_bnds):
                    match.fragmented = True

    def select_best_matches(self, pick='single'):
        # for each groundtruth SV we identify and mark the best-match call SV
        for candidate_list in self.complete_detections.values():
            # sort complete matches by increasing breakend distance
            candidate_list.sort(key=lambda m: m.total_dist)

        for gt_sv in self.svs:
            for match in self.complete_detections[gt_sv.parent_id]:
                call_sv = match.call_sv
                if call_sv.match_id: continue
                call_sv.match_type.append(gt_sv.parent_type)
                call_sv.match_id.append(gt_sv.parent_id)
                call_sv.match_score.append((0, match.total_dist))

                gt_sv.match_type.append(call_sv.parent_type)
                gt_sv.match_id.append(call_sv.parent_id)
                gt_sv.match_score.append((0, match.total_dist))

                if pick == 'single':
                    break

    @staticmethod
    def compute_metrics(pred_df, truth_df):
        pred_count = len(pred_df)
        truth_count = len(truth_df)

        tp_mask = pred_df['match_type'].notna()
        gt_tp_mask = truth_df['match_type'].notna()

        gt_complete_mask = truth_df['detection_mode'] == 'complete'

        # ============================== type-aggregate precision/recall/f1 ===============================
        tp = int(tp_mask.sum())
        fp = int((~tp_mask).sum())
        fn = int((~gt_tp_mask).sum())
        tp_truth = int(gt_tp_mask.sum())

        precision = float(tp / pred_count) if pred_count != 0 else 0
        recall = float(tp_truth / truth_count) if truth_count != 0 else 0
        f1 = (2 * precision * recall) / (precision + recall) if precision + recall != 0 else 0

        total_n_bnds = truth_df['num_breakends'].sum()
        num_matched_bnds = total_n_bnds - truth_df['spurious'].sum()

        spurious_rate = (pred_df['spurious'] / pred_df['num_breakends']).mean()

        num_complete_matched_bnds = truth_df[gt_complete_mask]['num_breakends'].sum()
        p_bnds = num_matched_bnds / total_n_bnds if total_n_bnds else 0
        complete_p_bnds = num_complete_matched_bnds / total_n_bnds if total_n_bnds else 0

        complete_comp = pred_df[pred_df['detection_mode'] == 'complete']
        avg_jitter = complete_comp["optimal_total_distance"].mean()

        avg_single_match_gt = (truth_df['num_candidates'] == 1).mean() if truth_count else 0.0
        avg_single_match_gt_complete = ((truth_df['num_candidates'] == 1) & gt_tp_mask).mean() if truth_count else 0.0

        truth_modes = truth_df['detection_mode'].value_counts().to_dict()
        pred_modes = pred_df['detection_mode'].value_counts().to_dict()

        return {
            'Valid Calls': pred_count,
            'Truth Valid SVs': truth_count,
            'TP-pred': tp,
            'TP-truth': tp_truth,
            'FP': fp,
            'FN': fn,
            'precision': precision,
            'recall': recall,
            'f1': f1,

            'Proportion truth breakend matched (all matches)': float(p_bnds),
            'Proportion truth breakend (complete matches)': float(complete_p_bnds),
            'Mean breakend accuracy in bp (optimal matches)': pred_df["optimal_total_distance"].mean(),
            'Mean breakend accuracy in bp (complete matches)': avg_jitter,
            'Mean GT breakend miss rate (optimal matches)': truth_df["optimal_miss"].mean(),
            'Mean extra breakend rate (optimal matches)': pred_df["optimal_extra"].mean(),
            'Average number spurious breakend rate (all matches)': spurious_rate,
            'Average number of GT SV single match (all matches)': avg_single_match_gt,
            'Average number of GT SV single match (complete matches)': avg_single_match_gt_complete,
            'Average number fragmented calls (optimal matches)': truth_df['fragmented'].mean(),
            'Average number discontiguous calls (optimal matches)': (pred_df['optimal_contiguity'] == 'Discontiguous').mean(),
            'Average number full matches (optimal matches)': truth_df['optimal_full'].mean(),

            'complete-pred': pred_modes.get('complete', 0),
            'partial-pred': pred_modes.get('partial', 0),
            'aggregate-pred': pred_modes.get('aggregate', 0),
            'miss-pred': pred_modes.get('miss', 0),
            'complete-truth': truth_modes.get('complete', 0),
            'partial-truth': truth_modes.get('partial', 0),
            'aggregate-truth': truth_modes.get('aggregate', 0),
            'miss-truth': truth_modes.get('miss', 0),
        }

    def create_outputs(self, output_folder):
        out = Path(output_folder)

        # Get Performance Metrics
        pred_df = self.pred_callset.callset2dataframe()
        gt_df = self.callset2dataframe()

        pred_df.to_csv(out / 'pred.csv', index=False)
        gt_df.to_csv(out / 'truth.csv', index=False)

        performance = self.compute_metrics(pred_df, gt_df)
        readable_perf = format_performance_log(performance)
        logging.info(f"Performance Metrics:\n{readable_perf}")
        with open(out / 'summary.json', 'w') as json_file:
            json.dump(performance, json_file, indent=4)

        write_bench_vcf(self.pred_callset, out / 'fp.vcf', out / 'tp-pred.vcf', is_gt=False)
        write_bench_vcf(self, out / 'fn.vcf', out / 'tp-truth.vcf', is_gt=True)

        all_matches = [m for sublist in self.complete_detections.values() for m in sublist] + \
                      [m for sublist in self.incomplete_detections.values() for m in sublist]

        df_columns = ['GT ID', 'GT type', 'is_optimal_GT', 'Call ID', 'Call type', 'is_optimal_call', '#BND unmatched',
                      'Total dist BND', 'detection mode', 'fragmented', 'Matched BNDS']

        detections_data = [
            (m.gt_sv.parent_id, m.gt_sv.parent_type, m == m.gt_sv.optimal_candidate, m.call_sv.parent_id, m.call_sv.parent_type,
             m == m.call_sv.optimal_candidate, m.diff_bnd, m.total_dist, m.mode, m.fragmented, stringify_matched_bnds(m.matched_bnds))
            for m in all_matches
        ]
        df_detections = pd.DataFrame(detections_data, columns=df_columns)
        df_detections.to_csv(out / 'detections.csv', index=False)

        optimal_matches = [m for m in all_matches if m.gt_sv.optimal_candidate == m or m.call_sv.optimal_candidate == m]
        optimal_data = [
            (m.gt_sv.parent_id, m.gt_sv.parent_type, m == m.gt_sv.optimal_candidate, m.call_sv.parent_id, m.call_sv.parent_type,
             m == m.call_sv.optimal_candidate, m.diff_bnd, m.total_dist, m.mode, m.fragmented, stringify_matched_bnds(m.matched_bnds))
            for m in optimal_matches
        ]
        df_optimal = pd.DataFrame(optimal_data, columns=df_columns)

        gt_misses = gt_df[gt_df['detection_mode'] == 'miss']
        call_misses = pred_df[pred_df['detection_mode'] == 'miss']

        miss_rows = []
        for _, row in gt_misses.iterrows():
            miss_rows.append([row['parent_id'], row['parent_type'], True, 'NA', 'NA', False, 'NA', 'NA', 'miss', 'NA', 'NA'])
        for _, row in call_misses.iterrows():
            miss_rows.append(['NA', 'NA', False, row['parent_id'], row['parent_type'], True, 'NA', 'NA', 'miss', 'NA', 'NA'])

        if miss_rows:
            df_optimal = pd.concat([df_optimal, pd.DataFrame(miss_rows, columns=df_columns)], ignore_index=True)

        df_optimal.to_csv(out / 'optimal.csv', index=False)

        def save_json_counts(df, filter_column, primary_col, secondary_col, filename):
            """Uses Pandas to group data into a nested dictionary structure."""
            if df.empty:
                counts_dict = {}
            else:
                filtered_df = df[df[filter_column] == True] if filter_column is not None else df
                # Aggregate counts: [Primary Type][Detection Mode][Secondary Type]
                counts = filtered_df.groupby([primary_col, 'detection mode', secondary_col]).size()

                # Convert MultiIndex Series to nested dictionary
                counts_dict = {}
                for (p_type, mode, s_type), val in counts.items():
                    counts_dict.setdefault(p_type, {}).setdefault(mode, {})[s_type] = int(val)

            with open(out / f"{filename}.json", 'w') as f:
                json.dump(counts_dict, f, indent=4)
            return counts_dict

        count_detections_gt = save_json_counts(df_detections, None, 'GT type', 'Call type', 'count_detections_gt')
        count_detections_call = save_json_counts(df_detections, None, 'Call type', 'GT type', 'count_detections_call')
        count_optimal_gt = save_json_counts(df_optimal, 'is_optimal_GT', 'GT type', 'Call type', 'count_optimal_gt')
        count_optimal_call = save_json_counts(df_optimal, 'is_optimal_call', 'Call type', 'GT type', 'count_optimal_call')

        return {
            'pred_df': pred_df,
            'truth_df': gt_df,
            'optimal_df': df_optimal,
            'count_optimal_gt': count_optimal_gt,
            'count_optimal_calls': count_optimal_call,
        }


def format_performance_log(perf_dict):
    lines = ["\n" + "-" * 100, f"{'Verix Benchmarking Summary':^50}", "-" * 100]

    for key, value in perf_dict.items():
        val_str = f"{value:>10}"
        if isinstance(value, float):
            val_str = f"{value:>10.8f}"
        lines.append(f"{key:<70} : {val_str}")
    lines.append("-" * 100)
    return "\n".join(lines)
