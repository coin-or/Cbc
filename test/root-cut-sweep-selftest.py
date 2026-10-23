#!/usr/bin/env python3
"""Unit tests plus optional fixture-backed checks of root-cut replay settings."""

import importlib.util
import re
from pathlib import Path
import subprocess
import sys
import unittest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("root_cut_sweep", HERE / "root-cut-sweep.py")
sweep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sweep)
FIXTURES = Path(sys.argv.pop(1)).resolve() if len(sys.argv) > 1 else None


class SweepTests(unittest.TestCase):
    def test_selection_is_deterministic_and_stratified(self):
        population = {
            f"{tier}{i}": {"tier": tier}
            for tier in ["small", "medium", "large"] for i in range(20)
        }
        first = sweep.select_instances(population, 15)
        self.assertEqual(first, sweep.select_instances(population, 15))
        self.assertEqual(len(set(first)), 15)
        for tier in ["small", "medium", "large"]:
            self.assertEqual(sum(name.startswith(tier) for name in first), 5)

    def test_missing_result_is_not_a_bound(self):
        self.assertEqual(sweep.fields("unrelated output", "[root-result]"), {})
        self.assertFalse(sweep.finite(float("nan")))
        self.assertFalse(sweep.finite(1e50))


@unittest.skipUnless(FIXTURES, "pass a root-fixture directory for replay checks")
class ReplayTests(unittest.TestCase):
    def run_replay(self, *options, instance="pk1"):
        run = subprocess.run(
            [str(HERE / "mip-root-replay"), instance,
             "--fixture-dir=" + str(FIXTURES), "--no-heur",
             "--root-only", "--adaptive-skip", "--sec=10", *options],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        self.assertEqual(run.returncode, 0, run.stderr)
        return run.stdout

    def test_default_small_tier_and_completed_root(self):
        text = self.run_replay()
        settings = sweep.fields(text, "[root-settings]")
        root = sweep.fields(text, "[root-result]")
        self.assertEqual(settings["passes"], "-100")
        self.assertEqual(settings["dropPatience"], "1")
        self.assertEqual(root["completed"], "1")
        self.assertEqual(root["timeLimit"], "0")
        self.assertIn("nodes=0 ", text)

    def test_explicit_drop_and_medium_tier(self):
        text = self.run_replay(
            "--param=sizeSmallRows=0", "--param=sizeSmallCols=0",
            "--param=sizeLargeCols=100000", "--param=cutPassMedium=7",
            "--param=minDrop=0.123", "--param=passCuts=auto",
            "--param=passTreeCuts=2",
        )
        settings = sweep.fields(text, "[root-settings]")
        self.assertEqual(settings["passes"], "7")
        self.assertAlmostEqual(float(settings["minimumDrop"]), 0.123)
        self.assertEqual(settings["treePasses"], "2")

    def test_large_tier_and_auto_drop(self):
        text = self.run_replay(
            "--param=sizeSmallRows=0", "--param=sizeSmallCols=0",
            "--param=sizeLargeCols=0", "--param=cutPassLarge=9",
            "--param=minDrop=auto",
        )
        settings = sweep.fields(text, "[root-settings]")
        self.assertEqual(settings["passes"], "9")
        self.assertAlmostEqual(float(settings["minimumDrop"]), 1e-5)

    def test_dedicated_override_and_scaling(self):
        text = self.run_replay(
            "--param=passCuts=8", "--pass-cuts=-3",
            "--param=minDrop=0.02", "--min-drop-scale=0.5",
        )
        settings = sweep.fields(text, "[root-settings]")
        self.assertEqual(settings["passes"], "-3")
        self.assertAlmostEqual(float(settings["minimumDrop"]), 0.01)

    def test_time_limit_is_not_a_completed_root(self):
        root = sweep.fields(self.run_replay("--sec=0.000001"), "[root-result]")
        self.assertEqual(root["completed"], "0")
        self.assertEqual(root["timeLimit"], "1")

    def test_patience_one_preserves_default(self):
        options = ("--pass-cuts=50", "--param=minDrop=0.05")
        default = sweep.fields(self.run_replay(*options), "[root-result]")
        explicit = sweep.fields(self.run_replay(
            *options, "--param=rootCutDropPatience=1"), "[root-result]")
        for key in ("completed", "bound", "passes", "optimal", "infeasible"):
            self.assertEqual(default[key], explicit[key])

    def test_patience_is_ignored_for_negative_budget(self):
        default = sweep.fields(self.run_replay("--pass-cuts=-15"), "[root-result]")
        text = self.run_replay("--pass-cuts=-15", "--param=rootCutDropPatience=5",
                               "--log=2")
        patience = sweep.fields(text, "[root-result]")
        self.assertNotIn("[root-drop]", text)
        for key in ("completed", "bound", "passes"):
            self.assertEqual(default[key], patience[key])

    def test_patience_counts_consecutive_misses_and_resets(self):
        text = self.run_replay(
            "--pass-cuts=100", "--param=minDrop=0.5",
            "--param=rootCutDropPatience=5", "--log=2", instance="binkar10_1")
        records = [sweep.fields(line, "[root-drop]") for line in text.splitlines()
                   if line.startswith("[root-drop]")]
        self.assertTrue(records, text)
        misses = 0
        reset = False
        for record in records:
            previous = misses
            misses = misses + 1 if float(record["drop"]) < float(record["minimumDrop"]) else 0
            reset |= previous > 0 and misses == 0
            self.assertEqual(int(record["misses"]), misses)
            self.assertEqual(record["patience"], "5")
        self.assertTrue(any(int(record["misses"]) > 1 for record in records))
        self.assertTrue(reset, "fixture must exercise a successful round after a stall")

    def test_patience_rejects_zero(self):
        run = subprocess.run(
            [str(HERE / "mip-root-replay"), "pk1",
             "--fixture-dir=" + str(FIXTURES),
             "--param=rootCutDropPatience=0"],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertNotEqual(run.returncode, 0)

    def test_zerohalf_cadence_and_initial_call(self):
        for period in (1, 2, 4, 8):
            text = self.run_replay(
                "--pass-cuts=-20", "--param=cutSwitchOff=zeroHalf:0",
                "--param=zeroHalfRootPeriod=" + str(period), "--log=2")
            passes = [int(value) for value in re.findall(
                r"\[root-cadence\] generator=ZeroHalf pass=(\d+)", text)]
            self.assertTrue(passes, text)
            self.assertEqual(passes[0], 1)
            self.assertTrue(all((value - 1) % period == 0 for value in passes))
            generator = sweep.fields(text, "[root-generator] name=ZeroHalf")
            self.assertEqual(int(generator["calls"]), len(passes))
            self.assertEqual(int(generator["rootPeriod"]), period)
            root = sweep.fields(text, "[root-result]")
            self.assertEqual(passes, list(range(1, int(root["passes"]) + 1, period)))

    def test_zerohalf_zero_period_preserves_default(self):
        baseline = sweep.fields(self.run_replay(), "[root-result]")
        zero = sweep.fields(self.run_replay("--param=zeroHalfRootPeriod=0"), "[root-result]")
        for key in ("completed", "bound", "passes"):
            self.assertEqual(baseline[key], zero[key])

    def test_combined_round_selector(self):
        for scope in (1, 2):
            text = self.run_replay(
                "--param=rootCutFilterScope=" + str(scope),
                "--param=rootCutFilterScore=1", "--param=cutFilterAlways=on",
                "--param=rootCutDropPatience=2", "--log=2", instance="binkar10_1")
            selections = [sweep.fields(line, "[root-selection]")
                          for line in text.splitlines()
                          if line.startswith("[root-selection]")]
            self.assertTrue(selections, text)
            self.assertTrue(any(int(r["kept"]) < int(r["candidates"]) for r in selections))
            self.assertEqual(sweep.fields(text, "[root-result]")["completed"], "1")

    def test_scheduled_reuse_respects_period_and_cap(self):
        for period in (1, 2, 4, 8):
            text = self.run_replay(
                "--param=rootCutReusePeriod=" + str(period),
                "--param=rootCutReuseMaxCuts=2", "--param=rootCutDropPatience=2",
                "--log=2", instance="p200x1188c")
            records = [sweep.fields(line, "[root-reuse]") for line in text.splitlines()
                       if line.startswith("[root-reuse]")]
            self.assertTrue(records, text)
            self.assertTrue(any(int(r["added"]) > 0 for r in records))
            for record in records:
                self.assertEqual(int(record["pass"]) % period, 0)
                self.assertLessEqual(int(record["added"]), 2)
                self.assertLessEqual(int(record["added"]), int(record["violated"]))
            self.assertEqual(sweep.fields(text, "[root-result]")["completed"], "1")

    def test_reuse_zero_preserves_default(self):
        default = sweep.fields(self.run_replay(), "[root-result]")
        zero = sweep.fields(self.run_replay("--param=rootCutReusePeriod=0"), "[root-result]")
        for key in ("completed", "bound", "passes"):
            self.assertEqual(default[key], zero[key])


if __name__ == "__main__":
    unittest.main()
