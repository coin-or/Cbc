#!/usr/bin/env python3
"""Repair phases depend on genuine feasibility, without changing default RNG."""
import argparse
import csv
import importlib.util
from pathlib import Path
import shutil
import subprocess

HERE = Path(__file__).resolve().parent
WORK = HERE / ".spp-grasp-phase-search-selftest"


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
        native = WORK / "native.mps"
        helper.write_mps(native, ["E"] * 2, [
            (2, [(0, 1)]), (2, [(1, 1)]), (1, [(0, 1), (1, 1)]),
        ])
        impossible = WORK / "half-native-impossible.mps"
        triangle = [
            (1, [(0, 1), (1, 1)]), (1, [(1, 1), (2, 1)]),
            (1, [(0, 1), (2, 1)]),
        ]
        private = [(0, [(row, 1)]) for row in range(3, 6)]
        helper.write_mps(impossible, ["E"] * 6, triangle + private)
        switch = WORK / "half-native-feasible.mps"
        helper.write_mps(switch, ["E"] * 6,
                         triangle + private + [(4, [(0, 1)])])
        serial = 0

        def run(model, mode, seed=17, initial=2):
            nonlocal serial
            serial += 1
            progress = WORK / f"{serial}.progress.tsv"
            trace = WORK / f"{serial}.trace.tsv"
            command = [
                args.bench, str(model), "--nolp", "--artificial=2",
                "--maxwork=0", "--noimp=0", "--stall=0", "--lns=200",
                "--iters=1", "--time=1e9", f"--ls={initial}",
                f"--lls={mode}", "--lp=0", f"--seed={seed}",
                "--row=1", "--col=0", "--calpha=0", "--totalwork=100000",
                f"--progress={progress}", f"--trace={trace}",
            ]
            result = subprocess.run(command, capture_output=True, text=True,
                                    timeout=15, check=True)
            assert not result.stderr, result.stderr
            events = list(csv.DictReader(progress.open(), delimiter="\t"))
            work = float(trace.read_text().split("\t")[1].split(":")[1])
            signature = [
                tuple((key, value) for key, value in sorted(event.items())
                      if key not in ("time", "seconds")) for event in events]
            return signature, events, work

        def first_real(events):
            return next((i for i, event in enumerate(events)
                         if int(event["artificial_used"]) == 0), None)

        for seed in (1, 17, 38):
            inherited = run(native, -1, seed)
            fixed_two = run(native, 2, seed)
            assert inherited[0] == fixed_two[0]
            assert inherited[2] == fixed_two[2]
            for mode, fixed in ((3, 2), (4, 2), (5, 1)):
                phased = run(native, mode, seed)
                reference = run(native, fixed, seed)
                assert phased[0] == reference[0], (mode, phased, reference)
                assert phased[2] == reference[2], (mode, phased, reference)
            for mode, fixed in ((3, 1), (4, 0), (5, 2)):
                phased = run(impossible, mode, seed)
                reference = run(impossible, fixed, seed)
                assert first_real(phased[1]) is None, phased
                assert phased[0] == reference[0], (mode, phased, reference)
                assert phased[2] == reference[2], (mode, phased, reference)

        exercised = set()
        for seed in range(1, 65):
            for mode, fixed in ((3, 1), (4, 0), (5, 2)):
                if mode in exercised:
                    continue
                phased = run(switch, mode, seed, initial=0)
                index = first_real(phased[1])
                if index is None or int(phased[1][index]["iteration"]) == 0:
                    continue
                reference = run(switch, fixed, seed, initial=0)
                assert phased[0][:index + 1] == reference[0][:index + 1]
                assert int(phased[1][-1]["artificial_used"]) == 0, phased
                exercised.add(mode)
            if len(exercised) == 3:
                break
        assert exercised == {3, 4, 5}, exercised
        print("PASS: native seeded work/RNG unchanged, artificial vectors never "
              "trigger phase switches, and all three genuine-feasibility "
              "switches retain the exact pre-switch trajectory")
    finally:
        shutil.rmtree(WORK)


if __name__ == "__main__":
    main()
