#!/usr/bin/env python3
"""
compare_benchmarks_selftest.py — unit tests for compare_benchmarks.py's status
classification, gap computation, and cross-check logic.

These exercise the specific traps this tool exists to avoid:
  - a "TIMEOUT(no_sol)" row must never be treated as having a real solution,
    even though Cbc's log/.sol machinery can print a fractional-relaxation
    "objective value" in that situation (summary.tsv is expected to leave
    objective="-" for these, which we assert on with --verify).
  - a proven-infeasible run must not be scored via primal/dual gap vs BKS.
  - OVERTIME must be treated as *worse* than TIMEOUT_NO_SOL even when the
    Cbc log shows a transient incumbent existed, as long as nothing was
    persisted (solution_found=0) -- see the chromaticindex1024-7 case
    investigated in the real experiment data this tool targets.
  - a validator-caught wrong answer (WRONG_OBJ / INFEASIBLE_WRONG / FAIL)
    must cost far more than a plain missing solution, regardless of weights.
  - results.tsv (cbc-workspace ./test) rows must classify the same way as
    an equivalent summary.tsv row once remapped.

Run directly:  python3 compare_benchmarks_selftest.py
"""

import io
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_benchmarks as cb


SUMMARY_HEADER = (
    "instance\tstatus\tobjective\tdual_bound\texpected\telapsed_s\t"
    "threads\tgap_field\tsolution_found\tproven_infeasible\ttimed_out\n"
)


def make_summary_tsv(rows):
    fh = tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False)
    fh.write(SUMMARY_HEADER)
    for r in rows:
        fh.write("\t".join(str(x) for x in r) + "\n")
    fh.close()
    return fh.name


class TestClassification(unittest.TestCase):
    def test_optimal(self):
        cat = cb.classify_status("SOLVED", 1, 0, 0)
        self.assertEqual(cat, "OPTIMAL")

    def test_infeasible_confirmed(self):
        cat = cb.classify_status("SOLVED(inf)", 0, 1, 0)
        self.assertEqual(cat, "INFEASIBLE_CONFIRMED")

    def test_timeout_with_real_solution(self):
        cat = cb.classify_status("TIMEOUT(gap=15.9%)", 1, 0, 1)
        self.assertEqual(cat, "TIMEOUT_WITH_SOL")

    def test_timeout_no_solution_is_never_a_solution(self):
        # The exact trap: Cbc's log/.sol can carry a fractional relaxation
        # "objective", but solution_found=0 must still classify as no-sol.
        cat = cb.classify_status("TIMEOUT(no_sol)", 0, 0, 1)
        self.assertEqual(cat, "TIMEOUT_NO_SOL")

    def test_overtime_with_flag_zero_even_if_log_shows_transient_incumbent(self):
        # Mirrors chromaticindex1024-7: cbc's log shows "BestSol: 4" moments
        # before being hard-killed, but nothing was ever persisted to a .sol
        # file, so solution_found is (correctly) 0 -- must classify OVERTIME,
        # not TIMEOUT_WITH_SOL, regardless of what the raw log looked like.
        cat = cb.classify_status("OVERTIME", 0, 0, 1)
        self.assertEqual(cat, "OVERTIME")

    def test_memlimit_no_sol(self):
        cat = cb.classify_status("MEMLIMIT(no_sol)", 0, 0, 0)
        self.assertEqual(cat, "TIMEOUT_NO_SOL")

    def test_wrong_obj_detected(self):
        cat = cb.classify_status("WRONG_OBJ", 1, 0, 0)
        self.assertEqual(cat, "WRONG")

    def test_infeasible_wrong_detected(self):
        cat = cb.classify_status("INFEASIBLE_WRONG", 0, 1, 0)
        self.assertEqual(cat, "WRONG")

    def test_invalid_sol_detected(self):
        cat = cb.classify_status("INVALID_SOL", 1, 0, 0)
        self.assertEqual(cat, "WRONG")

    def test_crash_and_asan(self):
        self.assertEqual(cb.classify_status("CRASH(exit=134)", 0, 0, 0), "ERROR")
        self.assertEqual(cb.classify_status("ASAN_ERROR", 0, 0, 0), "ERROR")


class TestGapsAndCost(unittest.TestCase):
    def setUp(self):
        self.weights = dict(cb.DEFAULT_WEIGHTS)

    def _row(self, **kw):
        base = dict(
            instance="x", status="SOLVED", objective=cb.NA, dual_bound=cb.NA,
            bks=cb.NA, elapsed_s=1.0, gap_field=cb.NA,
            solution_found=cb.NA, proven_infeasible=cb.NA, timed_out=cb.NA,
        )
        base.update(kw)
        import pandas as pd
        return pd.DataFrame([base])

    def test_optimal_cost_zero(self):
        df = self._row(status="SOLVED", objective=10, dual_bound=10, bks=10,
                        solution_found=1, proven_infeasible=0, timed_out=0)
        e = cb.enrich(df, self.weights)
        self.assertEqual(e["cost"].iloc[0], 0.0)
        self.assertTrue(e["concluded"].iloc[0])
        self.assertAlmostEqual(e["primal_gap_bks"].iloc[0], 0.0)
        self.assertAlmostEqual(e["dual_gap_bks"].iloc[0], 0.0)

    def test_infeasible_confirmed_no_gap(self):
        df = self._row(status="SOLVED(inf)", bks=cb.NA,
                        solution_found=0, proven_infeasible=1, timed_out=0)
        e = cb.enrich(df, self.weights)
        self.assertEqual(e["cost"].iloc[0], 0.0)
        self.assertTrue(e["concluded"].iloc[0])
        self.assertTrue(e["primal_gap_bks"].isna().iloc[0])
        self.assertTrue(e["dual_gap_bks"].isna().iloc[0])

    def test_timeout_no_sol_cost_matches_weight(self):
        df = self._row(status="TIMEOUT(no_sol)", bks=100,
                        solution_found=0, proven_infeasible=0, timed_out=1)
        e = cb.enrich(df, self.weights)
        self.assertEqual(e["cost"].iloc[0], self.weights["cost_no_solution"])
        self.assertFalse(e["concluded"].iloc[0])
        self.assertFalse(bool(e["has_solution"].iloc[0]))

    def test_overtime_cost_worse_than_no_sol(self):
        df = self._row(status="OVERTIME", bks=100,
                        solution_found=0, proven_infeasible=0, timed_out=1)
        e = cb.enrich(df, self.weights)
        self.assertEqual(e["cost"].iloc[0], self.weights["cost_overtime"])
        self.assertGreater(self.weights["cost_overtime"], self.weights["cost_no_solution"])

    def test_wrong_cost_dominates_by_default(self):
        df = self._row(status="WRONG_OBJ", objective=5, bks=100,
                        solution_found=1, proven_infeasible=0, timed_out=0)
        e = cb.enrich(df, self.weights)
        self.assertEqual(e["cost"].iloc[0], self.weights["cost_wrong"])
        self.assertGreater(self.weights["cost_wrong"], self.weights["cost_overtime"])

    def test_timeout_with_sol_gap_used_and_capped(self):
        df = self._row(status="TIMEOUT(gap=250%)", objective=5, dual_bound=1, bks=100,
                        gap_field=250.0, solution_found=1, proven_infeasible=0, timed_out=1)
        e = cb.enrich(df, self.weights)
        self.assertEqual(e["cost"].iloc[0], self.weights["gap_cap"])

    def test_primal_gap_sign_worse_than_bks_is_positive(self):
        # minimization: objective worse (higher) than bks -> positive primal gap
        df = self._row(status="TIMEOUT(gap=10%)", objective=110, dual_bound=90, bks=100,
                        gap_field=10.0, solution_found=1, proven_infeasible=0, timed_out=1)
        e = cb.enrich(df, self.weights)
        self.assertAlmostEqual(e["primal_gap_bks"].iloc[0], 10.0)
        self.assertAlmostEqual(e["dual_gap_bks"].iloc[0], 10.0)


class TestSummaryLoaderAndVerify(unittest.TestCase):
    def test_load_and_classify_realistic_rows(self):
        path = make_summary_tsv([
            ["inst_opt", "SOLVED", 1698, 1698, 1698, 0.02, 1, "0%", 1, 0, 0],
            ["inst_inf", "SOLVED(inf)", "-", "-", "-", 4, 1, "-", 0, 1, 0],
            ["inst_to_sol", "TIMEOUT(gap=15.9%)", 3650.13, 3149.07, 3311.18, 20019, 1, "15.9%", 1, 0, 1],
            ["inst_to_nosol", "TIMEOUT(no_sol)", "-", 202, 302, 20001, 1, "200%", 0, 0, 1],
            ["inst_overtime", "OVERTIME", "-", "-", 4, 20600, 1, "300%", 0, 0, 1],
            ["inst_wrong", "WRONG_OBJ", 100, 100, 90, 60, 1, "-", 1, 0, 0],
        ])
        try:
            raw = cb.load_summary_tsv(path)
            e = cb.enrich(raw, cb.DEFAULT_WEIGHTS)
            cats = dict(zip(e["instance"], e["category"]))
            self.assertEqual(cats["inst_opt"], "OPTIMAL")
            self.assertEqual(cats["inst_inf"], "INFEASIBLE_CONFIRMED")
            self.assertEqual(cats["inst_to_sol"], "TIMEOUT_WITH_SOL")
            self.assertEqual(cats["inst_to_nosol"], "TIMEOUT_NO_SOL")
            self.assertEqual(cats["inst_overtime"], "OVERTIME")
            self.assertEqual(cats["inst_wrong"], "WRONG")

            issues = cb.verify_consistency(e, "test")
            self.assertEqual(issues, [], f"unexpected inconsistencies: {issues}")
        finally:
            os.unlink(path)

    def test_verify_catches_inconsistent_flags(self):
        # Broken input: solution_found=1 but objective is empty ("-").
        path = make_summary_tsv([
            ["inst_bad", "TIMEOUT(gap=5%)", "-", 90, 100, 10, 1, "5%", 1, 0, 1],
        ])
        try:
            raw = cb.load_summary_tsv(path)
            e = cb.enrich(raw, cb.DEFAULT_WEIGHTS)
            issues = cb.verify_consistency(e, "test")
            self.assertTrue(any("solution_found=1 but objective is empty" in i for i in issues))
        finally:
            os.unlink(path)


class TestResultsTsvLoader(unittest.TestCase):
    def test_pass_optimal_row(self):
        fh = tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False)
        fh.write("instance\tstatus\telapsed_s\tnodes\tgap_pct\tis_optimal\tobj\tbound\n")
        fh.write("i1\tPASS\t0.1\t0\t\t1\t55\t55\n")
        fh.write("i2\tPASS\t10.0\t100\t12.3\t0\t60\t53\n")
        fh.write("i3\tOVERTIME\t180.0\t0\t\t0\t\t\n")
        fh.close()
        try:
            raw = cb.load_results_tsv(fh.name, bks_lookup={})
            e = cb.enrich(raw, cb.DEFAULT_WEIGHTS)
            cats = dict(zip(e["instance"], e["category"]))
            self.assertEqual(cats["i1"], "OPTIMAL")
            self.assertEqual(cats["i2"], "TIMEOUT_WITH_SOL")
            self.assertEqual(cats["i3"], "OVERTIME")
        finally:
            os.unlink(fh.name)


class TestRealExperimentData(unittest.TestCase):
    """
    Cross-checks against /home/haroldo/experiments/cbc real benchmark data
    when present on this machine (skipped otherwise), tying the tool's
    n_solved count back to the independently-computed count already printed
    in the harness's own report.txt ("Results: N/M solved/proved...").
    """
    DIR = "/home/haroldo/experiments/cbc/git_08_10_trust5_arm64"

    def test_n_solved_matches_report_txt(self):
        if not os.path.exists(os.path.join(self.DIR, "summary.tsv")):
            self.skipTest("real experiment data not present on this machine")
        raw = cb.load_summary_tsv(os.path.join(self.DIR, "summary.tsv"))
        e = cb.enrich(raw, cb.DEFAULT_WEIGHTS)
        n_solved = int(e["concluded"].sum())
        report_path = os.path.join(self.DIR, "report.txt")
        with open(report_path) as f:
            text = f.read()
        import re
        m = re.search(r"Results:\s+(\d+)/(\d+)\s+solved/proved", text)
        self.assertIsNotNone(m)
        self.assertEqual(n_solved, int(m.group(1)))


if __name__ == "__main__":
    unittest.main(verbosity=2)
