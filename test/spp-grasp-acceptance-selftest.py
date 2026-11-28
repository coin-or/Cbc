#!/usr/bin/env python3
"""Compile and run direct checks for SppGrasp's interchangeable LNS policies."""
import csv
import os
from pathlib import Path
import shlex
import subprocess

test_dir = Path(__file__).resolve().parent
source_dir = test_dir.parent / "src"
binary = test_dir / ".spp-grasp-acceptance-selftest"
compiler = shlex.split(os.environ.get("CXX", "c++"))
bench = Path(os.environ.get("SPP_GRASP_BENCH", test_dir / "spp-grasp-bench")).resolve()
work_dir = test_dir / ".spp-grasp-acceptance-selftest-work"


def write_cycle_model(path):
    rows = 8
    edges = [(i, (i + 1) % rows) for i in range(rows)]
    costs = [0, 4, -1, 3, 0, 2, 1, -2]
    columns = list(zip(costs, edges)) + [(6 + 0.5 * i, (i,)) for i in range(rows)]
    lines = ["NAME          ACCEPT", "ROWS", " N  OBJ"]
    lines += [f" E  R{i}" for i in range(rows)]
    lines += ["COLUMNS"]
    for j, (cost, covered) in enumerate(columns):
        name = f"X{j}"
        lines.append(f"    {name:<8}  {'OBJ':<8}  {cost:g}")
        lines += [f"    {name:<8}  {'R' + str(i):<8}  1" for i in covered]
    lines += ["RHS"]
    lines += [f"    {'RHS1':<8}  {'R' + str(i):<8}  1" for i in range(rows)]
    lines += ["BOUNDS"]
    lines += [f" BV {'BND1':<8}  X{j}" for j in range(len(columns))]
    lines += ["ENDATA"]
    path.write_text("\n".join(lines) + "\n")


def run_sa_integration(model, repetition):
    trace = work_dir / f"{repetition}.tr"
    btstats = work_dir / f"{repetition}.bt"
    acceptstats = work_dir / f"{repetition}.acc"
    command = [
        str(bench), str(model), "--nolp", "--iters=1", "--lns=200",
        "--time=1e9", "--row=6", "--col=0", "--lp=0", "--ls=0",
        "--lrow=6", "--lls=0", "--drows=1", "--dcols=1", "--dfrac=.25",
        "--rrows=1", "--stall=0", "--noimp=0", "--maxwork=0",
        "--amode=3", "--satemp=.5", "--sacool=1e9", "--sareheat=1",
        "--saschedule=1", "--sacalibrate=40", "--satarget=.3",
        "--saquantile=.5",
        "--seed=62", "--btprob=.25", f"--trace={trace}",
        f"--btstats={btstats}", f"--acceptstats={acceptstats}",
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    assert not result.stderr, result.stderr
    fields = result.stdout.strip().split("\t")
    assert abs(float(fields[3])) < 1.0e-7, fields
    with acceptstats.open() as stream:
        acceptance = next(csv.DictReader(stream, delimiter="\t"))
    trials = int(acceptance["trials"])
    accepted = int(acceptance["accepted"])
    rejected = int(acceptance["rejected"])
    assert trials == 200 and accepted + rejected == trials, acceptance
    assert int(acceptance["accepted_worse"]) > 0, acceptance
    assert int(acceptance["reheats"]) > 0 and int(acceptance["invalid"]) == 0, acceptance
    assert int(acceptance["calibration_trials"]) == 40, acceptance
    assert int(acceptance["calibration_positive"]) > 0, acceptance
    assert float(acceptance["initial_temperature"]) > 0.0, acceptance
    with btstats.open() as stream:
        repair = next(csv.DictReader(stream, delimiter="\t"))
    calls = int(repair["calls"])
    exits = sum(int(repair[key]) for key in
                ("work_exits", "node_exits", "exhausted", "invalid_partial"))
    assert calls > 0 and calls == exits, repair
    path = trace.read_text().strip().split("\t")
    signature = (path[0], path[1].split(":")[:2],
                 [(row.split(":")[0], row.split(":")[1], row.split(":")[3])
                  for row in path[2:]])
    return signature, acceptance, repair


try:
    subprocess.run(
        compiler + [
            "-std=c++11", "-Wall", "-Wextra", "-pedantic",
            "-I", str(source_dir),
            str(test_dir / "spp-grasp-acceptance-selftest.cpp"),
            "-o", str(binary),
        ],
        check=True,
    )
    subprocess.run([str(binary)], check=True)
    assert bench.is_file(), f"build the benchmark first: {bench}"
    work_dir.mkdir(exist_ok=True)
    model = work_dir / "acceptance-cycle.mps"
    write_cycle_model(model)
    first = run_sa_integration(model, "first")
    second = run_sa_integration(model, "second")
    assert first == second, "same-seed SA trajectory or diagnostic counts changed"
    print("PASS: seeded SA integration accepts uphill trials, reheats, preserves best and checks --btstats")
finally:
    binary.unlink(missing_ok=True)
    if work_dir.exists():
        for path in work_dir.iterdir():
            if path.is_file():
                path.unlink()
        work_dir.rmdir()
