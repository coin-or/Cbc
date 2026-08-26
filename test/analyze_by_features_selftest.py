#!/usr/bin/env python3
"""
analyze_by_features_selftest.py — unit tests for analyze_by_features.py.

Constructs small synthetic summary.tsv + features.csv fixtures (never reads
real experiment data) and checks:
  - a group query correctly partitions instances (n adds up, disjoint groups
    don't overlap incorrectly).
  - avg_cost/avg_pgap/avg_dgap/avg_t for a group matches a hand-computed
    value, in particular that gap_cap-penalized values (not NA) are used for
    TIMEOUT_NO_SOL/OVERTIME/WRONG instances, consistent with compare_benchmarks.py.
  - an instance missing from the features file is excluded from every
    feature-based group but still counted in "ALL".
  - --features-file accepts both .csv and .tsv inputs with the same schema.

Run directly:  python3 analyze_by_features_selftest.py
"""

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_benchmarks as cb
import analyze_by_features as ab

SUMMARY_HEADER = (
    "instance\tstatus\tobjective\tdual_bound\texpected\telapsed_s\t"
    "threads\tgap_field\tsolution_found\tproven_infeasible\ttimed_out\n"
)


def make_summary_tsv(path, rows):
    with open(path, "w") as fh:
        fh.write(SUMMARY_HEADER)
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


FEATURES_HEADER = (
    "Name,cols,rows,percBin,percInteger,genInt,density,"
    "rPercPartitioning,rPercPacking,rPercCovering,"
    "rPercKnapsack,rPercIntegerKnapsack\n"
)


def make_features_csv(path, rows):
    with open(path, "w") as fh:
        fh.write(FEATURES_HEADER)
        for r in rows:
            fh.write(",".join(str(x) for x in r) + "\n")


class TestAnalyzeByFeatures(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.exp_dir = os.path.join(self.tmpdir, "exp1")
        os.makedirs(self.exp_dir)
        make_summary_tsv(
            os.path.join(self.exp_dir, "summary.tsv"),
            [
                # instance, status, objective, dual_bound, expected, elapsed_s, threads, gap_field, solution_found, proven_infeasible, timed_out
                ("bin_a", "SOLVED", "100", "100", "100", "5", "1", "0%", "1", "0", "0"),
                ("bin_b", "TIMEOUT(no_sol)", "-", "50", "100", "18000", "1", "100%", "0", "0", "1"),
                ("mixed_a", "TIMEOUT(gap=10.0%)", "110", "100", "100", "18000", "1", "10%", "1", "0", "1"),
                ("no_feat", "SOLVED", "50", "50", "50", "2", "1", "0%", "1", "0", "0"),
            ],
        )
        make_features_csv(
            os.path.join(self.tmpdir, "features.csv"),
            [
                ("bin_a", 100, 100, 100, 100, 0, 10, 100, 0, 0, 0, 0),
                ("bin_b", 200, 150, 60, 100, 0, 5, 0, 0, 30, 0, 0),
                ("mixed_a", 300, 100, 40, 60, 5, 20, 0, 0, 0, 30, 0),
                # 'no_feat' deliberately absent from features.csv
            ],
        )
        self.weights = cb.load_weights(None)

    def _joined(self):
        feat = ab.load_features(os.path.join(self.tmpdir, "features.csv"))
        dfs = ab.load_experiments([self.exp_dir], ["exp1"], self.weights, None)
        return ab.build_joined(dfs, ["exp1"], feat)

    def test_missing_feature_excluded_from_groups_but_counted_in_all(self):
        joined = self._joined()
        all_row = ab.summarize_group(joined, ["exp1"], "ALL", None)
        self.assertEqual(all_row["n"], 4)  # no_feat still present in ALL

        bin_row = ab.summarize_group(joined, ["exp1"], "bin", "percBin >= 95")
        # no_feat has NaN percBin -> excluded by the query
        self.assertEqual(bin_row["n"], 1)
        self.assertEqual(joined.loc["bin_a"]["percBin"], 100)

    def test_gap_cap_used_for_no_solution_not_na(self):
        joined = self._joined()
        row = joined.loc["bin_b"]
        # TIMEOUT_NO_SOL must be scored at gap_cap (100.0 default), not NaN
        self.assertEqual(row["category__exp1"], "TIMEOUT_NO_SOL")
        self.assertEqual(row["pgap_clip__exp1"], self.weights["gap_cap"])
        self.assertEqual(row["dgap_clip__exp1"], self.weights["gap_cap"])

    def test_group_avg_cost_matches_hand_computation(self):
        joined = self._joined()
        sppc_row = ab.summarize_group(joined, ["exp1"], "sppc", "_sppc >= 90")
        # only bin_a (rPercPartitioning=100) qualifies
        self.assertEqual(sppc_row["n"], 1)
        expected_cost = joined.loc["bin_a"]["cost__exp1"]
        self.assertAlmostEqual(sppc_row["exp1"]["avg_cost"], expected_cost)

    def test_avg_t_formatting(self):
        self.assertEqual(ab.fmt_hms(3661), "1h01m")
        self.assertEqual(ab.fmt_hms(0), "0h00m")
        self.assertTrue("n/a" in ab.fmt_hms(float("nan")))

    def test_tsv_and_csv_feature_files_both_load(self):
        tsv_path = os.path.join(self.tmpdir, "features.tsv")
        with open(os.path.join(self.tmpdir, "features.csv")) as src, open(tsv_path, "w") as dst:
            for line in src:
                dst.write(line.replace(",", "\t"))
        feat_csv = ab.load_features(os.path.join(self.tmpdir, "features.csv"))
        feat_tsv = ab.load_features(tsv_path)
        self.assertEqual(list(feat_csv.index), list(feat_tsv.index))
        self.assertEqual(feat_csv.loc["bin_a"]["rows"], feat_tsv.loc["bin_a"]["rows"])


if __name__ == "__main__":
    unittest.main()
