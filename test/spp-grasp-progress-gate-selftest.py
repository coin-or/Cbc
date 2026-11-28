#!/usr/bin/env python3
"""Deterministic checks for the opt-in real-incumbent progress gate."""
import csv
from pathlib import Path
import shutil
import subprocess
import sys

bench = Path(sys.argv[1] if len(sys.argv) > 1 else "./spp-grasp-bench").resolve()
root = Path(__file__).resolve().parent
work = root / ".spp-grasp-progress-gate-selftest"
shutil.rmtree(work, ignore_errors=True)
work.mkdir(exist_ok=True)


def write_model(name, body):
    path = work / f"{name}.mps"
    path.write_text(body)
    return path


def run(label, model, extra):
    progress = work / f"{label}.progress.tsv"
    trace = work / f"{label}.trace.tsv"
    command = [
        str(bench), str(model), "--seed=271", "--iters=1", "--lns=100",
        "--noimp=0", "--maxwork=1000000", "--drows=1", "--stall=0",
        f"--progress={progress}", f"--trace={trace}", *extra,
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode != 0 or result.stderr:
        raise AssertionError(
            f"{label}: exit {result.returncode}\n{result.stdout}\n{result.stderr}"
        )
    fields = result.stdout.strip().split("\t")
    if len(fields) < 11:
        raise AssertionError(f"{label}: malformed result {result.stdout!r}")
    with progress.open(newline="") as stream:
        events = list(csv.DictReader(stream, delimiter="\t"))
    return fields, events, trace.read_text()


def deterministic_trace(trace):
    fields = trace.strip().split("\t")
    result = [":".join(fields[1].split(":")[:2])]
    result.extend(
        ":".join(event.split(":")[:2] + event.split(":")[3:])
        for event in fields[2:]
    )
    return result


def deterministic_events(events):
    return [
        tuple((key, event[key]) for key in event if key != "time")
        for event in events
    ]


feasible = write_model(
    "feasible",
    """NAME          GATE
ROWS
 N  OBJ
 E  R1
 E  R2
COLUMNS
    X1        OBJ       5       R1        1
    X2        OBJ       5       R2        1
    X3        OBJ       7       R1        1
    X3        R2        1
RHS
    RHS1      R1        1       R2        1
BOUNDS
 BV BND1      X1
 BV BND1      X2
 BV BND1      X3
ENDATA
""",
)
full, full_events, _ = run("full", feasible, [])
gated1, events1, trace1 = run("gated1", feasible, ["--realstallwork=1"])
gated2, events2, trace2 = run("gated2", feasible, ["--realstallwork=1"])
if any(int(event["kind"]) == 4 for event in full_events):
    raise AssertionError("disabled progress gate unexpectedly fired")
if not any(int(event["kind"]) == 4 for event in events1):
    raise AssertionError("progress gate did not stop a feasible stalled run")
if ([v for i, v in enumerate(gated1) if i not in (5, 6)]
    != [v for i, v in enumerate(gated2) if i not in (5, 6)]
    or deterministic_trace(trace1) != deterministic_trace(trace2)
    or deterministic_events(events1) != deterministic_events(events2)):
    raise AssertionError("same-seed gate trajectory is not deterministic")
if float(gated1[3]) > float(full[3]) + 1.0e-8:
    raise AssertionError("gate discarded a better solution already found")
if int(gated1[10]) != 0:
    raise AssertionError("gate returned a solution using artificial columns")

infeasible = write_model(
    "artificial_only",
    """NAME          NOREAL
ROWS
 N  OBJ
 E  R1
 E  R2
 E  R3
COLUMNS
    X1        OBJ       1       R1        1
    X1        R2        1
    X2        OBJ       1       R2        1
    X2        R3        1
    X3        OBJ       1       R1        1
    X3        R3        1
RHS
    RHS1      R1        1       R2        1
    RHS1      R3        1
BOUNDS
 BV BND1      X1
 BV BND1      X2
 BV BND1      X3
ENDATA
""",
)
artificial, art_events, _ = run(
    "artificial", infeasible,
    ["--artificial=1.1", "--artwork=1000", "--artgrow=10", "--artrounds=2",
     "--realstallwork=1", "--maxwork=1000000"],
)
if not any(int(event["kind"]) == 3 for event in art_events):
    raise AssertionError("artificial-only run did not exercise penalty resets")
if any(int(event["kind"]) == 4 for event in art_events):
    raise AssertionError("gate stopped a run without a real incumbent")
if int(artificial[10]) <= 0:
    raise AssertionError("infeasible real model did not retain its artificial result")

shutil.rmtree(work)
print("PASS: disabled/default behavior, deterministic stop, best preservation, "
      "and no-real-incumbent penalty escalation")
