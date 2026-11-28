#!/usr/bin/env python3
"""Check opt-in stage-local native coverage and explicit policy precedence."""
import argparse
import importlib.util
from pathlib import Path
import re
import shutil
import subprocess

HERE = Path(__file__).resolve().parent
WORK = HERE / ".spp-grasp-native-coverage-selftest"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cbc", default=str(HERE.parent / "src/cbc"))
    parser.add_argument("--reference-cbc")
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location(
        "feature_fixtures", HERE / "spp-grasp-adaptive-features-selftest.py")
    fixtures = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fixtures)
    if WORK.exists():
        raise RuntimeError(f"Refusing to overwrite {WORK}")
    WORK.mkdir()
    try:
        def model(name, rows, cols, **kwargs):
            file = WORK / (name + ".mps")
            fixtures.write_mps(file, rows, cols, **kwargs)
            return file

        full = model("full", ["E", "E"], [
            (-2, [(0, 1)]), (1, [(1, 1)]), (0, [(0, 1), (1, 1)])])
        partial = model("partial", ["E", "E"], [
            (2, [(0, 1)]), (1, [(0, 1), (1, 1)]),
            (3, [(0, 1), (1, 1)])])
        artificial = model("artificial", ["E", "E"], [
            (1, [(0, 1), (1, 1)]), (2, [(0, 1), (1, 1)])])
        packed = model("packed", ["E", "E", "L"], [
            (1, [(0, 1), (2, 1)]), (1, [(1, 1), (2, 1)]),
            (3, [(0, 1), (1, 1)])])
        forced = model("forced", ["E", "E"], [
            (0, [(0, 1)]), (1, [(0, 1), (1, 1)]),
            (2, [(1, 1)])], fixed=(0,))
        complemented = model("complemented", ["E"], [
            (-2, [(0, 1)]), (3, [(0, -1)])], rhs_values=[0])

        def run(file, threshold, controls=()):
            coverage_controls = [] if args.cbc == args.reference_cbc else [
                "-sppGraspMinNativeCoverage", str(threshold)]
            command = [args.cbc, str(file), "-log", "2", "-threads", "1",
                       "-randomSeed", "71", "-randomCbcSeed", "71",
                       "-sec", "5", "-preprocess", "off",
                       "-preRootLPStrenghtening", "off",
                       "-sppGrasp", "on", "-sppGraspArtificial", "2",
                       *coverage_controls,
                       "-sppGraspWork", "0.001", "-sppGraspNoImprove", "1",
                       *controls, "-solve", "-quit"]
            result = subprocess.run(command, capture_output=True, text=True,
                                    check=True, timeout=15)
            text = result.stdout + result.stderr
            assert "No match for" not in text, text
            return text

        for file, threshold, fragment in (
            (full, 1, "native coverage 2/2 (1), forced 0, minimum 1: run"),
            (partial, .5, "native coverage 1/2 (0.5), forced 0, minimum 0.5: run"),
            (partial, .5001, "native coverage 1/2 (0.5), forced 0, minimum 0.5001: skip"),
            (artificial, .5, "native coverage 0/2 (0), forced 0, minimum 0.5: skip"),
            (packed, .5, "native coverage 0/2 (0), forced 0, minimum 0.5: skip"),
            (complemented, 1, "native coverage 1/1 (1), forced 0, minimum 1: run"),
        ):
            text = run(file, threshold)
            assert "SppGrasp original-entry " + fragment in text, text
        text = run(forced, 1)
        assert "minimum 1: skip" in text and "forced 0," not in text, text
        text = run(partial, .5, ("-sppGraspArtificial", "0"))
        assert "incomplete coverage and artificial slacks disabled: skip" in text, text
        text = run(artificial, .5, ("-feas", "on"))
        assert "feasibility pump switched off" not in text, text
        text = run(artificial, .5)
        assert "feasibility pump switched off" not in text and "FPump invocation:" in text, text
        default = run(partial, 0)
        assert "original-entry native coverage" not in default, default
        off = run(full, 1, ("-sppGrasp", "off"))
        assert "SppGrasp original-entry" not in off and "SppGrasp: first" not in off, off
        explicit_pump = run(full, 1, ("-feas", "on"))
        assert "feasibility pump switched off" not in explicit_pump, explicit_pump
        if args.reference_cbc:
            def trajectory(text):
                return re.findall(r"SppGrasp: (?:first .*?|no feasible solution .*?)after \d+ LNS iterations"
                                  r"|Objective value:\s+\S+|Enumerated nodes:\s+\d+"
                                  r"|LP-solves=\d+ passes=\d+", text)
            for file in (full, partial, artificial):
                current = run(file, 0)
                executable = args.cbc
                args.cbc = args.reference_cbc
                reference = run(file, 0)
                args.cbc = executable
                assert trajectory(current) == trajectory(reference), (current, reference)
        print("PASS: native coverage 0/partial/1, no artificial or packing-coupled credit, forced/complement normalization, legacy default and explicit overrides")
    finally:
        shutil.rmtree(WORK)


if __name__ == "__main__":
    main()
