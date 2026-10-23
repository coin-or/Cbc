#!/usr/bin/env python3
"""Aggregate work caps must survive artificial-penalty epoch restarts."""
import argparse
import csv
import importlib.util
from pathlib import Path
import shutil
import subprocess

HERE = Path(__file__).resolve().parent
WORK = HERE / ".spp-grasp-total-work-selftest"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bench", default=str(HERE / "spp-grasp-bench"))
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location(
        "features", HERE / "spp-grasp-adaptive-features-selftest.py")
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    if WORK.exists():
        raise RuntimeError("refusing to overwrite " + str(WORK))
    WORK.mkdir()
    try:
        model = WORK / "odd-cycle.mps"
        helper.write_mps(model, ["E"] * 3, [
            (1, [(0, 1), (1, 1)]),
            (1, [(1, 1), (2, 1)]),
            (1, [(0, 1), (2, 1)]),
        ])
        def run(label, total):
            progress = WORK / (label + ".progress.tsv")
            trace = WORK / (label + ".trace.tsv")
            command = [
                args.bench, str(model), "--nolp", "--artificial=2",
                "--artgrow=2", "--artrounds=20", "--artwork=100",
                "--maxwork=0", "--noimp=0", "--stall=0", "--lns=100000",
                "--iters=1", "--time=1e9", "--ls=0", "--lls=0", "--lp=0",
                "--seed=17", f"--totalwork={total}", f"--progress={progress}",
                f"--trace={trace}",
            ]
            result = subprocess.run(command, capture_output=True, text=True,
                                    timeout=15, check=True)
            assert not result.stderr, result.stderr
            fields = result.stdout.strip().split("\t")
            events = list(csv.DictReader(progress.open(), delimiter="\t"))
            work = float(trace.read_text().split("\t")[1].split(":")[1])
            return fields, events, work
        capped, events, work = run("capped", 500)
        uncapped, longer, long_work = run("uncapped", 0)
        repeated, again, repeat_work = run("repeated", 500)
        # The core checks at move boundaries. An indivisible construction/move
        # may overshoot, but must not acquire a fresh aggregate cap each epoch.
        assert 500 <= work < 1000, (work, capped, events)
        assert long_work > 2 * work, (work, long_work)
        assert len({event["penalty_round"] for event in events}) > 1, events
        signature = lambda es: [
            tuple((key, value) for key, value in sorted(event.items())
                  if key not in ("time", "seconds")) for event in es]
        assert signature(events) == signature(again), (events, again)
        assert work == repeat_work, (work, repeat_work)
        for invalid in ("-1", "nan", "inf", "garbage", ""):
            result = subprocess.run(
                [args.bench, str(model), "--nolp", f"--totalwork={invalid}"],
                capture_output=True, text=True, timeout=15)
            assert result.returncode != 0, (invalid, result.stdout)
        print("PASS: aggregate cap survives penalty epochs, includes overshoot, "
              "preserves seeded trajectory and rejects invalid budgets")
    finally:
        shutil.rmtree(WORK)


if __name__ == "__main__":
    main()
