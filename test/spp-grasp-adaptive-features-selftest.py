#!/usr/bin/env python3
"""Check structural feature summaries on tiny SPP models."""
import csv
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys


HERE = Path(__file__).resolve().parent
BENCH = Path(sys.argv[1] if len(sys.argv) > 1 else HERE / "spp-grasp-bench").resolve()
WORK = HERE / ".spp-grasp-adaptive-feature-selftest"


def write_mps(path, rows, columns, rhs_values=None, fixed=()):
    rhs_values = rhs_values or [1] * len(rows)
    lines = ["NAME          FEATURES", "ROWS", " N  OBJ"]
    lines += [f" {kind}  R{i}" for i, kind in enumerate(rows)]
    lines.append("COLUMNS")
    for j, (cost, entries) in enumerate(columns):
        name = f"X{j}"
        lines.append(f"    {name:<8}  {'OBJ':<8}  {cost:g}")
        lines += [f"    {name:<8}  {'R' + str(i):<8}  {value:g}"
                  for i, value in entries]
    lines.append("RHS")
    lines += [f"    {'RHS1':<8}  {'R' + str(i):<8}  {rhs_values[i]:g}"
              for i in range(len(rows))]
    lines.append("BOUNDS")
    for j in range(len(columns)):
        if j in fixed:
            lines.append(f" LI {'BND1':<8}  {'X' + str(j):<8}  1")
            lines.append(f" UI {'BND1':<8}  {'X' + str(j):<8}  1")
        else:
            lines.append(f" BV {'BND1':<8}  X{j}")
    lines.append("ENDATA")
    path.write_text("\n".join(lines) + "\n")


def features(model, name, extra=()):
    path = WORK / f"{name}.tsv"
    command = [str(BENCH), str(model), "--info", f"--features={path}", *extra]
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    assert not result.stderr, result.stderr
    with path.open() as stream:
        row = next(csv.DictReader(stream, delimiter="\t"))
    for key in ("load_seconds", "feature_seconds"):
        value = float(row[key])
        assert math.isfinite(value) and value >= 0.0, row
    return row


def main():
    if WORK.exists():
        shutil.rmtree(WORK)
    WORK.mkdir()
    try:
        simple = WORK / "simple.mps"
        write_mps(simple, ["E", "E"], [
            (-5, [(0, 1)]),
            (1, [(0, 1), (1, 1)]),
            (3, [(1, 1)]),
            (0, [(0, 1), (1, 1)]),
        ])
        row = features(simple, "simple")
        assert row["native_singleton_rows"] == "2", row
        assert row["native_singleton_missing"] == "0", row
        assert float(row["row_options_min"]) == 3.0, row
        assert abs(float(row["row_options_mean"]) - 3.0) < 1e-12, row
        assert abs(float(row["singleton_premium_q25"]) + 0.5) < 1e-9, row
        assert abs(float(row["singleton_premium_median"])) < 1e-9, row
        assert abs(float(row["singleton_premium_q75"]) - 0.5) < 1e-9, row
        assert abs(float(row["column_degree_mean"]) - 1.5) < 1e-12, row
        assert float(row["lp_available"]) == 0.0, row
        native = WORK / "native-trivial.mps"
        write_mps(native, ["E", "E"], [(0, [(0, 1)]), (0, [(1, 1)])])
        native_progress = WORK / "native-trivial.progress.tsv"
        result = subprocess.run([
            str(BENCH), str(native), "--nolp", "--iters=1", "--lns=0",
            f"--progress={native_progress}",
        ], capture_output=True, text=True, check=True)
        assert not result.stderr, result.stderr
        with native_progress.open() as stream:
            events = list(csv.DictReader(stream, delimiter="\t"))
        assert len(events) == 1, events
        event = events[0]
        assert int(event["kind"]) == 1 and int(event["artificial_used"]) == 0, event
        assert float(event["value"]) == 0 and float(event["work"]) == 6, event

        artificial = WORK / "artificial.mps"
        write_mps(artificial, ["E", "E"], [
            (1, [(0, 1), (1, 1)]),
            (2, [(0, 1), (1, 1)]),
        ])
        row = features(artificial, "artificial", ("--artificial=2",))
        assert row["native_singleton_rows"] == "0", row
        assert row["native_singleton_missing"] == "2", row
        assert float(row["row_options_min"]) == 2.0, row

        # Each equality row has one coefficient, but its only variable also
        # appears in the packing row. These are row singletons, not independent
        # singleton columns and therefore cannot certify full row coverage.
        packed_singletons = WORK / "packed-singletons.mps"
        write_mps(packed_singletons, ["E", "E", "L"], [
            (1, [(0, 1), (2, 1)]),
            (1, [(1, 1), (2, 1)]),
        ], rhs_values=[1, 1, 1])
        row = features(packed_singletons, "packed-singletons", ("--artificial=2",))
        assert row["native_singleton_rows"] == "0", row
        assert row["native_singleton_missing"] == "2", row
        assert float(row["row_options_min"]) == 1.0, row

        complemented = WORK / "complemented.mps"
        write_mps(complemented, ["E"], [
            (2, [(0, 1)]),
            (3, [(0, -1)]),
        ], rhs_values=[0])
        row = features(complemented, "complemented")
        assert row["complemented_columns"] == "1", row
        assert float(row["complemented_fraction"]) == 0.5, row

        restricted = WORK / "restricted.mps"
        write_mps(restricted, ["E", "L"], [
            (1, [(0, 1), (1, 1)]),
            (2, [(0, 1)]),
        ], rhs_values=[1, 0])
        row = features(restricted, "restricted")
        assert int(row["restricting_rows"]) == 1, row
        assert int(row["restricted_columns"]) == 1, row
        assert float(row["restricted_column_fraction"]) == 0.5, row

        forced = WORK / "forced.mps"
        write_mps(forced, ["E", "E"], [
            (0, [(0, 1)]),
            (1, [(0, 1), (1, 1)]),
            (2, [(1, 1)]),
        ], fixed=(0,))
        row = features(forced, "forced")
        assert int(row["forced_columns"]) == 1, row
        assert abs(float(row["forced_fraction"]) - 1.0 / 3.0) < 1e-9, row

        # A normal run solves the LP before extracting features; info-only does not.
        lp_features = WORK / "lp.tsv"
        result = subprocess.run([
            str(BENCH), str(simple), "--iters=1", "--lns=0",
            f"--features={lp_features}",
        ], capture_output=True, text=True, check=True)
        assert not result.stderr, result.stderr
        with lp_features.open() as stream:
            row = next(csv.DictReader(stream, delimiter="\t"))
        assert row["lp_available"] == "1", row
        assert float(row["lp_fractionality_mean"]) >= 0.0, row

        # Conflicting real choices make the tiny model infeasible; artificial
        # rows permit deterministic penalty escalations, which are logged as
        # epochs and never mistaken for objective improvements.
        conflict = WORK / "artificial-conflict.mps"
        write_mps(conflict, ["E", "E", "L"], [
            (1, [(0, 1), (2, 1)]),
            (1, [(1, 1), (2, 1)]),
        ], rhs_values=[1, 1, 1])
        progress_signatures = []
        for label in ("first", "repeat"):
            path = WORK / f"{label}.progress.tsv"
            result = subprocess.run([
                str(BENCH), str(conflict), "--nolp", "--iters=1",
                "--lns=20", "--maxwork=0", "--noimp=1", "--time=1e9",
                "--artificial=2", "--artgrow=2", "--artrounds=2",
                "--artwork=1", f"--progress={path}",
            ], capture_output=True, text=True, check=True)
            assert not result.stderr, result.stderr
            with path.open() as stream:
                events = list(csv.DictReader(stream, delimiter="\t"))
            epochs = [int(event["penalty_round"]) for event in events
                      if int(event["kind"]) == 3]
            assert epochs and epochs == list(range(1, len(epochs) + 1)), events
            assert all(int(event["artificial_used"]) > 0 for event in events
                       if int(event["kind"]) in (1, 3)), events
            progress_signatures.append([
                tuple(event[key] for key in
                      ("iteration", "kind", "penalty_round", "artificial_used",
                       "value", "work"))
                for event in events
            ])
        assert progress_signatures[0] == progress_signatures[1], progress_signatures
        print("PASS: native trivial feasibility is logged without requiring improvement; coverage excludes artificials/packing-coupled columns; optional LP features")
    finally:
        shutil.rmtree(WORK, ignore_errors=True)


if __name__ == "__main__":
    try:
        main()
    except (AssertionError, subprocess.CalledProcessError) as error:
        print(error, file=sys.stderr)
        sys.exit(1)
