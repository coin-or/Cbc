#!/usr/bin/env python3
"""Check strict native completion, saved solutions and pre-LP progress."""
import argparse
import csv
import importlib.util
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile

import cbc_pdi


HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cbc", type=Path, default=HERE.parent / "src/cbc")
    parser.add_argument("--reference-cbc", type=Path)
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location(
        "native_fixtures", HERE / "spp-grasp-adaptive-features-selftest.py")
    fixtures = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fixtures)
    with tempfile.TemporaryDirectory(prefix="spp-native-prelp-") as directory:
        work = Path(directory)
        cases = [
            ("full", ["E", "E"],
             [(2, [(0, 1)]), (3, [(1, 1)]), (1, [(0, 1), (1, 1)])],
             {}, 5, ()),
            ("negative", ["E", "E"],
             [(-2, [(0, 1)]), (1, [(1, 1)]), (0, [(0, 1), (1, 1)])],
             {}, -1, ()),
            ("competing-private", ["E", "E"],
             [(7, [(0, 1)]), (2, [(0, 1)]), (3, [(1, 1)]),
              (1, [(0, 1), (1, 1)])], {}, 5, ()),
            ("complemented", ["E"],
             [(-2, [(0, 1)]), (3, [(0, -1)])],
             {"rhs_values": [0]}, 0, ()),
            ("offset", ["E", "E"],
             [(2, [(0, 1)]), (3, [(1, 1)]), (1, [(0, 1), (1, 1)])],
             {}, -2, ()),
            ("partial", ["E", "E"],
             [(2, [(0, 1)]), (1, [(0, 1), (1, 1)])],
             {}, None, ("-sppGraspArtificial", "5")),
            ("packed", ["E", "E", "L"],
             [(1, [(0, 1), (2, 1)]), (1, [(1, 1), (2, 1)]),
              (3, [(0, 1), (1, 1)])], {}, None, ()),
            ("benign-not-private", ["E", "L"],
             [(1, [(0, 1), (1, -1)])],
             {"rhs_values": [1, 0]}, None, ()),
            ("choose-strict-not-cheaper-coupled", ["E", "E", "L"],
             [(2, [(0, 1)]), (3, [(1, 1)]),
              (1, [(0, 1), (1, 1), (2, 1)]), (0, [(0, 1), (2, -1)])],
             {"rhs_values": [1, 1, 0]}, 5, ()),
            ("forced", ["E", "E"],
             [(0, [(0, 1)]), (1, [(0, 1), (1, 1)]), (2, [(1, 1)])],
             {"fixed": (0,)}, None, ()),
            ("spp-off", ["E", "E"],
             [(2, [(0, 1)]), (3, [(1, 1)]), (1, [(0, 1), (1, 1)])],
             {}, None, ("-sppGrasp", "off")),
            ("cutoff", ["E", "E"],
             [(2, [(0, 1)]), (3, [(1, 1)]), (1, [(0, 1), (1, 1)])],
             {}, None, ("-cutoff", "2")),
            ("maximize", ["E", "E"],
             [(2, [(0, 1)]), (3, [(1, 1)]), (1, [(0, 1), (1, 1)])],
             {}, None, ("-maximize",)),
            ("staged", ["E", "E"],
             [(2, [(0, 1)]), (3, [(1, 1)]), (1, [(0, 1), (1, 1)])],
             {}, None, ("-sppGraspRootSequence", "spp-fpump-staged")),
            ("ordered", ["E", "E"],
             [(2, [(0, 1)]), (3, [(1, 1)]), (1, [(0, 1), (1, 1)])],
             {}, None, ("-sppGraspRootSequence", "fpump-spp")),
            ("mipstart", ["E", "E"],
             [(2, [(0, 1)]), (3, [(1, 1)]), (1, [(0, 1), (1, 1)])],
             {}, None, ("-mipstart", str(work / "start.mst"))),
            ("expired", ["E", "E"],
             [(2, [(0, 1)]), (3, [(1, 1)]), (1, [(0, 1), (1, 1)])],
             {}, None, ("-sec", "0")),
        ]
        (work / "start.mst").write_text("0 X0 0\n1 X1 0\n2 X2 1\n")

        def solve(model, name, controls=(), binary=None):
            progress = work / (name + ".csv")
            solution = work / (name + ".sol")
            command = [
                str(binary or args.cbc), str(model), "-log", "2", "-threads", "1",
                "-randomSeed", "71", "-randomCbcSeed", "71", "-sec", "5",
                "-presolve", "off", "-preprocess", "off",
                "-preRootLPStrenghtening", "off", "-sppGrasp", "on",
                "-sppGraspFpump", "on", "-sppGraspWork", ".001",
                "-sppGraspNoImprove", "1", *controls,
                "-progressCsv", str(progress), "-solve", "-solu", str(solution),
                "-quit"]
            result = subprocess.run(command, capture_output=True, text=True,
                                    check=True, timeout=20)
            text = result.stdout + result.stderr
            assert not any(error in text for error in (
                "No match for", "Illegal keyword", "String of", "unrecognized")), text
            events = list(csv.DictReader(progress.open())) if progress.exists() else []
            return text, events, solution

        for name, rows, columns, kwargs, objective, extra in cases:
            model = work / (name + ".mps")
            fixtures.write_mps(model, rows, columns, **kwargs)
            offset = 7 if name == "offset" else 0
            if offset:
                model.write_text(model.read_text().replace(
                    "RHS\n", "RHS\n    RHS1      OBJ       7\n"))
            text, events, solution = solve(
                model, name, ("-sppGraspNativePreLp", "on", *extra))
            accepted = "SppGrasp native pre-LP: real feasible certificate" in text
            assert accepted == (objective is not None), (name, text)
            native = [event for event in events
                      if event["source"] == "SppGraspNativePreLp"]
            if not accepted:
                assert "SppGrasp native pre-LP: skipped:" in text, (name, text)
                assert not native, (name, events)
                if name == "mipstart":
                    assert "MIPStart provided solution with cost 1" in text, text
                continue
            assert len(native) == 1, (name, events)
            assert float(native[0]["time"]) > 0, native
            assert abs(float(native[0]["primal"]) - objective) < 1e-7, native
            assert float(native[0]["dual"]) <= -1e49, native
            assert events[-1]["source"] == "final", (name, events, text)
            times = [float(event["time"]) for event in events]
            assert times == sorted(times), events
            lines = solution.read_text().splitlines()
            assert "Optimal" in lines[0], (name, lines)
            values = {}
            for line in lines[1:]:
                fields = line.split()
                if len(fields) >= 3 and fields[1].startswith("X"):
                    values[int(fields[1][1:])] = float(fields[2])
            activity = [0.0] * len(rows)
            cost = -offset
            for j, (coefficient, entries) in enumerate(columns):
                value = values.get(j, 0.0)
                assert value in (0, 1), (name, j, value)
                if j in kwargs.get("fixed", ()):
                    assert value == 1, (name, j, value)
                cost += coefficient * value
                for i, coefficient in entries:
                    activity[i] += coefficient * value
            for i, kind in enumerate(rows):
                rhs = kwargs.get("rhs_values", [1] * len(rows))[i]
                assert (abs(activity[i] - rhs) < 1e-7 if kind == "E"
                        else activity[i] <= rhs + 1e-7), (name, activity)
            saved_cost = float(re.search(r"objective value\s+(\S+)", lines[0])[1])
            assert abs(saved_cost - cost) < 1e-7, (name, lines[0], cost)

        model = work / "full.mps"
        default, events, _ = solve(model, "default")
        automatic, auto_events, _ = solve(model, "auto", ("-sppGraspNativePreLp", "auto"))
        off, off_events, _ = solve(model, "off", ("-sppGraspNativePreLp", "off"))
        assert "real feasible certificate" in default, default
        assert "real feasible certificate" in automatic, automatic
        assert "native pre-LP:" not in off, off
        assert len([event for event in events
                    if event["source"] == "SppGraspNativePreLp"]) == 1, events
        assert len([event for event in auto_events
                    if event["source"] == "SppGraspNativePreLp"]) == 1, auto_events
        assert not any(event["source"] == "SppGraspNativePreLp" for event in off_events)
        assert events[-1]["source"] == "final", events
        for name, controls in (
                ("heur-off", ("-heur", "off")),
                ("spp-off", ("-sppGrasp", "off"))):
            disabled, disabled_events, _ = solve(model, name, controls)
            assert "real feasible certificate" not in disabled, disabled
            assert not any(event["source"] == "SppGraspNativePreLp"
                           for event in disabled_events), disabled_events

        def trajectory(text):
            return re.findall(
                r"Objective value:\s+\S+|Enumerated nodes:\s+\d+"
                r"|SppGrasp:.*?after \d+ LNS iterations", text)

        assert trajectory(default) == trajectory(automatic), (default, automatic)
        if args.reference_cbc:
            reference, _, _ = solve(model, "reference", binary=args.reference_cbc)
            assert trajectory(off) == trajectory(reference), (off, reference)
        flags = next(line.split("=", 1)[1] for line in
                     (HERE / "Makefile").read_text().splitlines()
                     if line.startswith("CBCLIB_CFLAGS ="))
        guards = work / "guards"
        subprocess.run([
            str(HERE.parent / "libtool"), "--quiet", "--tag=CXX", "--mode=link",
            *shlex.split(os.environ.get("CXX", "c++")), "-std=c++11",
            "-DHAVE_CONFIG_H", "-I" + str(HERE.parent / "src"), *shlex.split(flags),
            str(HERE / "spp-grasp-native-prelp-guards.cpp"),
            str(HERE.parent / "src/libCbc.la"), "-o", str(guards),
        ], check=True)
        for mode in ("callback", "event", "cuts", "sos", "better", "worse", "hotstart"):
            trace = work / (mode + "-guard.csv")
            result = subprocess.run([str(guards), str(model), mode, str(trace)],
                                    capture_output=True, text=True, timeout=20)
            assert result.returncode == 0, (mode, result.stdout, result.stderr)
            assert "real feasible certificate" not in result.stdout, (mode, result.stdout)
            if trace.exists():
                assert "SppGraspNativePreLp" not in trace.read_text(), mode
        timeout_model = work / "preprocess-timeout.mps"
        timeout_columns = [
            (2, [(0, 1)]), (2, [(1, 1)]), (2, [(2, 1)]),
            (1, [(0, 1), (1, 1)]), (1, [(1, 1), (2, 1)]),
            (1, [(0, 1), (2, 1)]),
        ]
        fixtures.write_mps(timeout_model, ["E"] * 3, timeout_columns)
        for mode in ("preprocess-empty", "preprocess-original", "preprocess-existing"):
            trace = work / (mode + ".csv")
            result = subprocess.run(
                [str(guards), str(timeout_model), mode, str(trace)],
                capture_output=True, text=True, timeout=15, check=True)
            assert "Pre-processing stopped on time limit" in result.stdout, result.stdout
            log = work / (mode + ".log")
            log.write_text(result.stdout)
            events = list(csv.DictReader(trace.open()))
            assert events[0]["source"] == "start" and events[-1]["source"] == "final", events
            assert sum(e["source"] == "final" for e in events) == 1, events
            times = [float(e["time"]) for e in events]
            assert times == sorted(times) and times[-1] > 0, events
            lines = Path(str(trace) + ".sol").read_text().splitlines()
            assert lines[0].startswith("Stopped on time"), lines[0]
            parsed = cbc_pdi.build_trace(cbc_pdi.run_events(log, trace, "csv"), False)
            metrics = cbc_pdi.integrate(parsed, 10.0)
            assert parsed.status == "stopped", parsed.status
            if mode == "preprocess-empty":
                assert "no integer solution" in lines[0], lines[0]
                assert float(events[-1]["primal"]) >= 1e49, events
                assert metrics["pdi"] == 2.0 and metrics["primal"] is None, metrics
                continue
            assert "no integer solution" not in lines[0], lines[0]
            values = {}
            for line in lines[1:]:
                fields = line.split()
                if len(fields) >= 3:
                    values[int(fields[0])] = float(fields[1])
            activity = [0.0] * 3
            objective = 0.0
            for j, (cost, entries) in enumerate(timeout_columns):
                value = values.get(j, 0.0)
                assert value in (0, 1), (mode, j, value)
                objective += cost * value
                for i, coefficient in entries:
                    activity[i] += coefficient * value
            assert activity == [1, 1, 1], (mode, activity, lines)
            saved = float(re.search(r"objective value\s+(\S+)", lines[0])[1])
            assert abs(saved - objective) < 1e-7, (mode, saved, objective)
            assert abs(float(events[-1]["primal"]) - objective) < 1e-7, events
            assert metrics["primal"] == objective, metrics
            if mode == "preprocess-original":
                original = [e for e in events if e["source"] == "SppGraspOriginal"]
                assert len(original) == 1 and 0 < float(original[0]["time"]) < times[-1], events
        print("PASS: strict actual-row native completion, cheapest private choice, "
              "complements, offsets, saved-vector feasibility, progress ordering, "
              "cutoff, MIP starts, opt-out, callbacks, custom objects/cuts, "
              "existing incumbents, unsupported modes, default auto and off parity; "
              "preprocessing timeouts preserve progress and saved integer vectors")


if __name__ == "__main__":
    main()
