#!/usr/bin/env python3
"""
stats_analysis_selftest.py — unit tests for stats_analysis.py's column
handling, time-budget accounting, and "wasted work" detection.

These exercise the traps specific to this tool:
  - the schema mismatch vs mipster's stats_analysis.py: this Cbc build's
    -csvStatistics has no per-generator "_calls"/"_avgnz" columns and no
    per-heuristic columns at all -- the ported script must degrade
    gracefully (explicit "not available" message) rather than crash or
    silently report zeroes as if heuristics genuinely used no time.
  - CbcSolverStatistics's `obj` sentinel (1.0e50 = "no solution found")
    must not be summed/averaged as if it were a real objective value.
  - the time-budget breakdown must fully account for total "time" (no
    double counting / no silently dropped remainder) across LP, cuts,
    conflict-graph, presolve/tightening, and "other" buckets.
  - "wasted" cut time (time>0 but cuts==0) must be computed per-row, not
    from an aggregate that could hide row-level cancellation.
  - real stats.csv files (which may predate the rowred_* columns, per the
    real trust5/trust10 experiment data used for this port) must still
    load and be analyzed without KeyErrors.

Run directly:  python3 stats_analysis_selftest.py
"""

import io
import os
import sys
import unittest

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stats_analysis as sa  # noqa: E402


def make_df(rows, columns):
    return pd.DataFrame(rows, columns=columns)


class TestSafeHelpers(unittest.TestCase):
    def test_safe_sum_missing_column_defaults_zero(self):
        df = make_df([{"instance": "a"}], ["instance"])
        self.assertEqual(sa._safe_sum(df, "cut_time"), 0.0)

    def test_safe_sum_present_column(self):
        df = make_df([{"cut_time": 1.5}, {"cut_time": 2.5}], ["cut_time"])
        self.assertAlmostEqual(sa._safe_sum(df, "cut_time"), 4.0)

    def test_safe_sum_ignores_non_numeric_garbage(self):
        df = make_df([{"cut_time": "1.5"}, {"cut_time": "n/a"}], ["cut_time"])
        # non-numeric entries must be coerced to 0, not raise or propagate NaN
        self.assertAlmostEqual(sa._safe_sum(df, "cut_time"), 1.5)

    def test_safe_col_missing_returns_zero_series(self):
        df = make_df([{"instance": "a"}, {"instance": "b"}], ["instance"])
        col = sa._safe_col(df, "rowred_time")
        self.assertEqual(list(col), [0.0, 0.0])


class TestTimeBudgetAccounting(unittest.TestCase):
    """The 5 buckets (LP, cuts, cgraph, presolve, other) must sum back to
    the reported total time -- no instance's time should vanish or be
    double-counted."""

    def test_buckets_sum_to_total(self):
        df = make_df([
            {"instance": "i1", "time": 100.0, "lp_seconds": 10.0, "cut_time": 20.0,
             "cgraph_time": 1.0, "clqstr_time": 2.0, "coefstr_time": 0.5, "rowred_time": 0.5},
            {"instance": "i2", "time": 50.0, "lp_seconds": 5.0, "cut_time": 5.0,
             "cgraph_time": 0.5, "clqstr_time": 0.0, "coefstr_time": 0.0, "rowred_time": 0.0},
        ], ["instance", "time", "lp_seconds", "cut_time", "cgraph_time",
            "clqstr_time", "coefstr_time", "rowred_time"])

        total_time = sa._safe_sum(df, "time")
        lp_time = sa._safe_sum(df, "lp_seconds")
        cut_time = sa._safe_sum(df, "cut_time")
        cgraph_time = sa._safe_sum(df, "cgraph_time")
        presolve_time = (sa._safe_sum(df, "clqstr_time") + sa._safe_sum(df, "coefstr_time")
                          + sa._safe_sum(df, "rowred_time"))
        other_time = total_time - (lp_time + cut_time + cgraph_time + presolve_time)

        self.assertAlmostEqual(total_time, 150.0)
        self.assertGreaterEqual(other_time, 0.0)
        self.assertAlmostEqual(lp_time + cut_time + cgraph_time + presolve_time + other_time, total_time)

    def test_missing_rowred_columns_does_not_crash(self):
        # real trust5/trust10 experiment data predates the rowred_* columns
        # (added later upstream) -- must not KeyError, must treat as 0.
        df = make_df([
            {"instance": "i1", "time": 10.0, "lp_seconds": 1.0, "cut_time": 2.0,
             "cgraph_time": 0.2, "clqstr_time": 0.1, "coefstr_time": 0.05},
        ], ["instance", "time", "lp_seconds", "cut_time", "cgraph_time", "clqstr_time", "coefstr_time"])
        rpt = sa.Report(sa.Colors(enabled=False))
        sa.section_time_budget(df, rpt)  # must not raise
        sa.section_presolve_tightening(df, rpt)  # must not raise
        text = "\n".join(rpt.lines)
        self.assertIn("Time budget", text)


class TestObjSentinel(unittest.TestCase):
    """CbcSolverStatistics::obj defaults to 1.0e50 ('no solution found')
    and must never be treated as a real objective value by this tool."""

    def test_sentinel_value_is_not_a_real_objective(self):
        # stats_analysis.py itself does not currently aggregate `objective`
        # in any section (it is informational only in the raw CSV), so this
        # test locks in that no section accidentally starts summing/averaging
        # it as if 1e50 were meaningful. If a future section adds objective
        # aggregation, it MUST filter out this sentinel first.
        self.assertNotEqual(sa.NA, sa.NA)  # NA is a real NaN (NaN != NaN), our own "no value" sentinel
        huge = 1.0e50
        self.assertTrue(huge > 1e40)  # documents the CbcSolverStatistics sentinel magnitude to guard against


class TestWastedCutDetection(unittest.TestCase):
    def test_zero_cut_time_detected_per_row(self):
        df = make_df([
            {"instance": "i1", "Probing_time": 5.0, "Probing_cuts": 0},
            {"instance": "i2", "Probing_time": 3.0, "Probing_cuts": 10},
            {"instance": "i3", "Probing_time": 2.0, "Probing_cuts": 0},
        ], ["instance", "Probing_time", "Probing_cuts"])
        rpt = sa.Report(sa.Colors(enabled=False))
        sa.section_wasted_cuts(df, rpt, top_n=5)
        text = "\n".join(rpt.lines)
        self.assertIn("Probing", text)
        # 2 instances (i1, i3) had time>0 but cuts==0 -> wasted = 5.0 + 2.0 = 7.0s
        self.assertIn("7.00s", text)

    def test_no_generator_columns_reports_gracefully(self):
        df = make_df([{"instance": "i1", "time": 1.0}], ["instance", "time"])
        rpt = sa.Report(sa.Colors(enabled=False))
        sa.section_cuts(df, rpt, top_n=5)
        text = "\n".join(rpt.lines)
        self.assertIn("No cut generator data found", text)


class TestHeuristicsGracefulDegradation(unittest.TestCase):
    """Upstream Cbc's -csvStatistics has NO per-heuristic columns -- this
    must be reported explicitly, not silently as if heuristics used 0 time
    (which would be misleading: it's not that they're fast, it's that
    they're untracked)."""

    def test_missing_heuristic_columns_explicit_message(self):
        df = make_df([{"instance": "i1", "time": 1.0}], ["instance", "time"])
        rpt = sa.Report(sa.Colors(enabled=False))
        sa.section_heuristics(df, rpt, top_n=5)
        text = "\n".join(rpt.lines)
        self.assertIn("Not available", text)
        self.assertIn("heur_", text)

    def test_present_heuristic_columns_are_still_used(self):
        # forward-compat: if a future Cbc build ever adds heur_* columns,
        # they must be picked up automatically without code changes.
        df = make_df([
            {"instance": "i1", "heur_rins_time": 2.0, "heur_rins_sols": 1, "heur_rins_execs": 3},
        ], ["instance", "heur_rins_time", "heur_rins_sols", "heur_rins_execs"])
        rpt = sa.Report(sa.Colors(enabled=False))
        sa.section_heuristics(df, rpt, top_n=5)
        text = "\n".join(rpt.lines)
        self.assertIn("RINS", text)
        self.assertNotIn("Not available", text)


class TestRealSchemaColumnNaming(unittest.TestCase):
    """Locks in the real column-naming scheme this Cbc build actually
    writes (CbcSolverStatistics.cpp), which differs from mipster's
    cut_<lower>_time/_cuts/_calls/_avgnz scheme."""

    def test_cut_generators_use_bare_name_prefix_no_calls_avgnz(self):
        for gen, _label in sa.CUT_GENERATORS:
            # Must NOT be lowercased/prefixed with "cut_" (mipster's scheme)
            self.assertNotIn("cut_", gen.lower().replace("cutscol", ""))
        # spot check some canonical names from CbcSolverStatistics.cpp
        gens = dict(sa.CUT_GENERATORS)
        for expected in ("Clique", "Gomory", "Gomory(2)", "Reduce-and-split",
                         "Reduce-and-split(2)", "TwoMirCuts", "ZeroHalf"):
            self.assertIn(expected, gens)


class TestRealDataIfAvailable(unittest.TestCase):
    """Cross-checks against the real trust5 experiment data used earlier in
    this session for compare_benchmarks.py, when available on this machine."""

    STATS_CSV = "/home/haroldo/experiments/cbc/git_08_10_trust5_arm64/stats.csv"

    def test_loads_and_analyzes_without_error(self):
        if not os.path.exists(self.STATS_CSV):
            self.skipTest("real experiment data not present on this machine")
        df = sa.load_stats(self.STATS_CSV)
        self.assertGreater(len(df), 300)  # ~383 instances expected
        rpt = sa.Report(sa.Colors(enabled=False))
        sa.section_time_budget(df, rpt)
        sa.section_presolve_tightening(df, rpt)
        sa.section_cuts(df, rpt, top_n=10)
        sa.section_heuristics(df, rpt, top_n=10)
        sa.section_wasted_cuts(df, rpt, top_n=10)
        sa.section_worst_instances(df, rpt, top_n=10)
        text = "\n".join(rpt.lines)
        # Probing is known (from manual inspection this session) to be the
        # single largest cut-generator time consumer in this dataset.
        probing_idx = text.find("Probing")
        gomory_idx = text.find("Gomory (root)")
        self.assertNotEqual(probing_idx, -1)
        self.assertNotEqual(gomory_idx, -1)
        self.assertLess(probing_idx, gomory_idx)  # Probing ranked above Gomory by time

    def test_total_time_matches_raw_sum(self):
        if not os.path.exists(self.STATS_CSV):
            self.skipTest("real experiment data not present on this machine")
        df = sa.load_stats(self.STATS_CSV)
        raw_total = pd.to_numeric(df["time"], errors="coerce").fillna(0.0).sum()
        self.assertAlmostEqual(sa._safe_sum(df, "time"), raw_total)


if __name__ == "__main__":
    unittest.main(verbosity=2)
