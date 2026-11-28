#!/usr/bin/env python3
"""Small exhaustive checks for spp-grasp-bench's bounded repair."""
import itertools
import csv
from pathlib import Path
import subprocess
import sys
import tempfile

bench = Path(sys.argv[1] if len(sys.argv) > 1 else "./spp-grasp-bench").resolve()


def check(directory, name, rows, columns, fixed=(), extra=()):
    model = directory / f"{name}.mps"
    lines = ["NAME          REPAIR", "ROWS", " N  OBJ"]
    lines += [f" {kind}  R{i}" for i, kind in enumerate(rows)]
    lines += ["COLUMNS"]
    for j, (cost, indices) in enumerate(columns):
        lines.append(f"    {'X' + str(j):<8}  {'OBJ':<8}  {cost}")
        lines += [f"    {'X' + str(j):<8}  {'R' + str(i):<8}  1" for i in indices]
    lines += ["RHS"]
    lines += [f"    {'RHS1':<8}  {'R' + str(i):<8}  1" for i in range(len(rows))]
    lines += ["BOUNDS"]
    for j in range(len(columns)):
        if j in fixed:
            lines.append(f" LI {'BND1':<8}  {'X' + str(j):<8}  1")
            lines.append(f" UI {'BND1':<8}  {'X' + str(j):<8}  1")
        else:
            lines.append(f" BV {'BND1':<8}  X{j}")
    lines += ["ENDATA"]
    model.write_text("\n".join(lines) + "\n")

    def feasible(x):
        if any(x[j] != 1 for j in fixed):
            return False
        activity = [sum(x[j] for j, (_, indices) in enumerate(columns) if i in indices)
                    for i in range(len(rows))]
        return all(a == 1 if kind == "E" else a <= 1
                   for a, kind in zip(activity, rows))

    optimum = min(sum(x[j] * c for j, (c, _) in enumerate(columns))
                  for x in itertools.product((0, 1), repeat=len(columns)) if feasible(x))
    objectives = {}
    greedy_trace = None
    for probability, nodes, work in ((0, 2000, "1e7"), (0.25, 2000, "1e7"),
                                    (1, 2000, "1e7"), (1, 1, "1e7"), (1, 2000, "1")):
        solution = directory / f"{name}-{probability}-{nodes}-{work}.sol"
        diagnostics = solution.with_suffix(".tsv")
        trace = solution.with_suffix(".tr")
        command = [
            str(bench), str(model), "--nolp", "--iters=1", "--lns=30",
            "--noimp=0", "--maxwork=1e8", "--time=100", "--artificial=2",
            "--ls=0", "--row=4", "--col=0", "--lp=0", "--lrow=4",
            "--calpha=0", "--drows=20", "--dcols=20", "--stall=0",
            f"--btprob={probability}", "--btwidth=3", "--btdepth=8",
            f"--btnodes={nodes}", f"--btwork={work}", f"--solfile={solution}",
            f"--btstats={diagnostics}",
            f"--trace={trace}",
            *extra,
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode == 0, result.stdout + result.stderr
        with diagnostics.open() as stream:
            stats = next(csv.DictReader(stream, delimiter="\t"))
        calls = int(stats["calls"])
        exits = sum(int(stats[key]) for key in
                    ("work_exits", "node_exits", "exhausted", "invalid_partial"))
        assert calls == exits, stats
        assert int(stats["no_leaf"]) <= calls
        assert int(stats["real_leaves"]) <= int(stats["leaves"])
        assert int(stats["feasibility_gains"]) <= int(stats["improved"]) <= calls
        if probability == 0:
            assert calls == 0
        else:
            assert calls + int(stats["skipped_feasible"]) > 0
        if work == "1":
            assert int(stats["work_exits"]) == calls, stats
        fields = trace.read_text().split()
        path = (fields[1].split(":")[:2],
                [(q[0], q[1], q[3]) for q in
                 (event.split(":") for event in fields[2:])])
        if probability == 0:
            greedy_trace = path
        if name == "greedy-trap" and "--btonlyart=1" in extra:
            assert calls == 0 and float(stats["extra_work"]) == 0
            assert path == greedy_trace, "skips must preserve RNG and work path"
        if result.stderr:
            raise AssertionError(result.stderr)
        values = [0] * len(columns)
        for line in solution.read_text().splitlines()[1:]:
            fields = line.split()
            values[int(fields[0])] = int(fields[2])
        assert feasible(values), (name, probability, values)
        objective = sum(values[j] * c for j, (c, _) in enumerate(columns))
        if nodes == 2000 and work == "1e7":
            objectives[probability] = objective
        assert objective >= optimum - 1e-8
        if probability == 1 and nodes == 2000 and work == "1e7" and "--btonlyart=1" not in extra:
            assert abs(objective - optimum) < 1e-8, (name, objective, optimum)
        if nodes == 1 or work == "1":
            assert objective == objectives[0], "budget exhaustion must keep greedy"
    if name == "greedy-trap" and "--btonlyart=1" not in extra:
        assert objectives[0] > objectives[1], objectives


with tempfile.TemporaryDirectory(prefix="spp-backtracking-") as tmp:
    directory = Path(tmp)
    check(directory, "greedy-trap", ["E"] * 3,
          [(1, [0, 1]), (2, [1, 2]), (1, [0]), (10, [2])])
    check(directory, "artificial", ["E"] * 4,
          [(1, [0, 1]), (3, [1, 2]), (1, [2, 3]), (3, [3, 0])])
    check(directory, "fixed-negative-packing", ["E", "E", "L"],
          [(4, [0]), (2, [1]), (-3, [2]), (-5, [1, 2]), (8, [1])],
          fixed=(0,))
    check(directory, "greedy-trap", ["E"] * 3,
          [(1, [0, 1]), (2, [1, 2]), (1, [0]), (10, [2])],
          extra=("--btonlyart=1", "--btrealfirst=1"))
    check(directory, "artificial", ["E"] * 4,
          [(1, [0, 1]), (3, [1, 2]), (1, [2, 3]), (3, [3, 0])],
          extra=("--btrealfirst=1",))
print("PASS: greedy/mixed/backtracking feasibility, exhaustive objectives, "
      "artificials, fixed columns and negative packing costs")
