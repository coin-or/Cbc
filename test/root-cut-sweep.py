#!/usr/bin/env python3
"""Screen root cut budgets on identical warm-started fixtures via GNU parallel.

Timing under parallel dispatch is exploratory, not a serial speedup claim.
Incomplete roots, cold starts and contradictory bounds are reported separately.
"""

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import random
import re
import resource
import statistics
import subprocess
import sys

SCRIPT = Path(__file__).resolve()
CONFIGS = {
    "baseline": [],
    "passes50": ["--pass-cuts=50"],
    "passes100": ["--pass-cuts=100"],
    "passes200": ["--pass-cuts=200"],
    "forced50": ["--pass-cuts=-50"],
    "forced100": ["--pass-cuts=-100"],
    "forced200": ["--pass-cuts=-200"],
    "drop_half": ["--min-drop-scale=0.5"],
    "drop_fifth": ["--min-drop-scale=0.2"],
    "skip_off": ["--no-adaptive-skip"],
    "shutdown_off": ["--param=cutSwitchOff=0"],
    "patience2": ["--param=rootCutDropPatience=2"],
    "patience3": ["--param=rootCutDropPatience=3"],
    "patience5": ["--param=rootCutDropPatience=5"],
    "half_patience3": ["--min-drop-scale=0.5", "--param=rootCutDropPatience=3"],
    "cap_double": ["--param=cutPassMedium=200", "--param=cutPassLarge=100"],
    "cap_large100": ["--param=cutPassLarge=100"],
    "skip_delay": ["--param=cutSkipMinTries=20"],
    "skip_short": ["--param=cutSkipInitialPeriod=1", "--param=cutSkipMaxPeriod=4"],
    "shutdown_zerohalf": ["--param=cutSwitchOff=zeroHalf:0"],
    "shutdown_lift": ["--param=cutSwitchOff=liftAndProject:0"],
    "zerohalf_k1": ["--param=cutSwitchOff=zeroHalf:0", "--param=zeroHalfRootPeriod=1"],
    "zerohalf_k2": ["--param=cutSwitchOff=zeroHalf:0", "--param=zeroHalfRootPeriod=2"],
    "zerohalf_k4": ["--param=cutSwitchOff=zeroHalf:0", "--param=zeroHalfRootPeriod=4"],
    "zerohalf_k8": ["--param=cutSwitchOff=zeroHalf:0", "--param=zeroHalfRootPeriod=8"],
}
PRIORITIES = [
    "neos-3381206-awhea", "p200x1188c", "n5-3", "swath3", "30n20b8",
    "binkar10_1", "rococoC10-001000", "ns1644855", "s250r10", "netdiversion",
    "roi2alpha3n4", "mas74", "gen-ip002", "gen-ip054", "pk1",
    "neos-3988577-wolgan", "mcsched", "pg5_34", "comp07-2idx",
]


def finite(value):
    return math.isfinite(value) and abs(value) < 1e40


def metadata(path):
    return dict(line.split(maxsplit=1) for line in path.read_text().splitlines()
                if len(line.split(maxsplit=1)) == 2)


def fixture_population(directory):
    records = {}
    for path in sorted(directory.glob("*.root.meta")):
        values = metadata(path)
        rows, cols = int(values["rows"]), int(values["cols"])
        name = path.name.removesuffix(".root.meta")
        if values.get("lpOptimal") != "1":
            continue
        records[name] = {
            "instance": name, "rows": rows, "cols": cols,
            "elements": int(values["elements"]),
            "integer_columns": int(values["integerColumns"]),
            "tier": "small" if rows < 500 or cols < 500
            else ("medium" if cols < 5000 else "large"),
        }
    return records


def select_instances(population, sample):
    if sample <= 0 or sample >= len(population):
        return sorted(population)
    selected = [name for name in PRIORITIES if name in population][:sample]
    rng = random.Random(20261005)
    remaining = {}
    for tier in ["small", "medium", "large"]:
        names = [name for name, record in population.items()
                 if record["tier"] == tier and name not in selected]
        rng.shuffle(names)
        remaining[tier] = names
    while len(selected) < sample:
        for tier in ["small", "medium", "large"]:
            if remaining[tier] and len(selected) < sample:
                selected.append(remaining[tier].pop())
    return sorted(selected)


def fields(text, prefix):
    line = next((line for line in text.splitlines() if line.startswith(prefix)), None)
    return dict(re.findall(r"(\w+)=([^\s]+)", line)) if line else {}


def worker(path):
    job = json.loads(path.read_text())
    output = Path(job["output"])
    command = [job["replay"], job["instance"],
               "--fixture-dir=" + job["fixture_dir"],
               "--data-dir=" + job["data_dir"],
               "--root-only", "--no-heur", "--adaptive-skip", "--log=0",
               "--sec=" + str(job["seconds"])] + job["options"]
    result = {
        "instance": job["instance"], "config": job["config"],
        "status": "ERROR", "completed": False, "bound": None,
        "seconds": None, "passes": None, "warm_iterations": None,
        "initial_bound": None, "reference": job["reference"],
        "repeat": job.get("repeat", 0),
    }
    with output.with_suffix(".log").open("w") as log:
        try:
            process = subprocess.run(
                command, stdout=log, stderr=subprocess.STDOUT,
                timeout=job["seconds"] + 30,
                preexec_fn=lambda: resource.setrlimit(
                    resource.RLIMIT_AS, (job["memory_bytes"], job["memory_bytes"])),
            )
            result["exit_code"] = process.returncode
        except subprocess.TimeoutExpired:
            result["status"] = "HARD_TIMEOUT"
            result["exit_code"] = None
    text = output.with_suffix(".log").read_text(errors="replace")
    load = fields(text, "[replay]")
    root = fields(text, "[root-result]")
    if root and result.get("exit_code") == 0:
        bound = float(root["bound"])
        result.update(
            completed=root["completed"] == "1",
            seconds=float(root["seconds"]), passes=int(root["passes"]),
            warm_iterations=int(load["warmStartIters"]),
            initial_bound=float(load["lpBound"]),
            optimal=root["optimal"] == "1",
            infeasible=root["infeasible"] == "1",
        )
        if not result["completed"]:
            result["status"] = "INFEASIBLE" if result["infeasible"] else "INCOMPLETE"
        elif result["warm_iterations"] != 0:
            result["status"] = "COLD_START"
        elif not finite(bound):
            result["status"] = "INVALID_BOUND"
        elif job["reference"] is not None and (
                bound > job["reference"] + max(1e-4, abs(job["reference"]) * 1e-4)):
            result["status"] = "INVALID_BOUND"
        else:
            result["status"] = "OK"
            result["bound"] = bound
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    if result["status"] not in ("OK", "INCOMPLETE", "INFEASIBLE"):
        print(f"{job['config']} {job['instance']}: {result['status']}; "
              f"see {output.with_suffix('.log')}", file=sys.stderr)
    return 0


def summarize(out):
    results = [json.loads(path.read_text()) for path in sorted((out / "results").glob("*.json"))]
    if not results:
        raise ValueError("No result files")
    with (out / "results.tsv").open("w") as handle:
        columns = sorted(set().union(*(result.keys() for result in results)))
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)
    grouped = {}
    for row in results:
        grouped.setdefault((row["config"], row["instance"]), []).append(row)
    indexed = {}
    for key, samples in grouped.items():
        if all(row["status"] == "OK" for row in samples):
            representative = dict(samples[0])
            representative["seconds"] = statistics.median(row["seconds"] for row in samples)
            representative["bound"] = min(row["bound"] for row in samples)
            indexed[key] = representative
    names = sorted({row["instance"] for row in results})
    print("Config          OK Incomplete Other Paired Better Worse Tie  SumTime(s) PairedTimeRatio")
    for config in sorted({row["config"] for row in results}):
        rows = [row for row in results if row["config"] == config]
        counts = {status: sum(row["status"] == status for row in rows)
                  for status in ["OK", "INCOMPLETE"]}
        better = worse = ties = 0
        base_time = candidate_time = 0.0
        for name in names:
            base, candidate = indexed.get(("baseline", name)), indexed.get((config, name))
            if not base or not candidate or base["status"] != "OK" or candidate["status"] != "OK":
                continue
            scale = max(abs(base["initial_bound"]), abs(base["bound"]),
                        abs(base["reference"] or 0.0), 1.0)
            delta = candidate["bound"] - base["bound"]
            tolerance = max(1e-6, scale * 1e-6)
            better += delta > tolerance
            worse += delta < -tolerance
            ties += abs(delta) <= tolerance
            base_time += base["seconds"]
            candidate_time += candidate["seconds"]
        total_time = sum(row["seconds"] or 0 for row in rows)
        ratio = candidate_time / base_time if base_time else float("nan")
        print(f"{config:15} {counts['OK']:2} {counts['INCOMPLETE']:10} "
              f"{len(rows)-sum(counts.values()):5} {better+worse+ties:6} "
              f"{better:6} {worse:5} {ties:3} {total_time:11.2f} {ratio:15.3f}")
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture-dir", type=Path)
    parser.add_argument("--data-dir", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--replay", type=Path, default=SCRIPT.parent / "mip-root-replay")
    parser.add_argument("--sample", type=int, default=70, help="0 means all fixtures")
    parser.add_argument("--instances", type=Path)
    parser.add_argument("--sec", type=float, default=30)
    parser.add_argument("--jobs", type=int, default=os.cpu_count())
    parser.add_argument("--memory-gb", type=float, default=2)
    parser.add_argument("--configs", default=",".join(CONFIGS))
    parser.add_argument("--base-param", action="append", default=[],
                        help="NAME=VALUE override shared by every configuration")
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--worker", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--summarize", type=Path)
    args = parser.parse_args()
    if args.worker:
        return worker(args.worker)
    if args.summarize:
        summarize(args.summarize.resolve())
        return 0
    if not args.fixture_dir or not args.data_dir or not args.out:
        parser.error("--fixture-dir, --data-dir and --out are required")
    if args.jobs < 1 or args.sec <= 0 or args.memory_gb <= 0 or args.repeats < 1:
        parser.error("jobs, sec, memory-gb and repeats must be positive")
    config_names = args.configs.split(",")
    if "baseline" not in config_names or any(name not in CONFIGS for name in config_names):
        parser.error("--configs must include baseline and only known configurations")
    if len(set(config_names)) != len(config_names):
        parser.error("--configs must not contain duplicates")
    configs = {name: ["--param=" + value for value in args.base_param] + CONFIGS[name]
               for name in config_names}
    if not args.replay.is_file() or not os.access(args.replay, os.X_OK):
        parser.error("Build mip-root-replay before running this sweep")
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    (out / "jobs").mkdir()
    (out / "results").mkdir()
    population = fixture_population(args.fixture_dir)
    names = (args.instances.read_text().splitlines() if args.instances
             else select_instances(population, args.sample))
    if not names or len(set(names)) != len(names) or any(name not in population for name in names):
        parser.error("Instance list must contain unique, available fixture names")
    with (args.data_dir / "bks.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        references = {}
        for row in reader:
            raw = row.get("bks", row.get("objective", "")) or ""
            try:
                value = float(raw)
            except ValueError:
                continue
            if finite(value):
                references[row["instance"].strip()] = value
    manifest = {
        "replay": str(args.replay.resolve()),
        "binary_sha256": hashlib.sha256(args.replay.read_bytes()).hexdigest(),
        "fixture_dir": str(args.fixture_dir.resolve()),
        "data_dir": str(args.data_dir.resolve()),
        "configs": configs, "repeats": args.repeats,
        "instances": [population[name] for name in names],
        "seconds": args.sec, "jobs": args.jobs, "memory_gb": args.memory_gb,
        "timing": "isolated dispatch" if args.jobs == 1 else "parallel screening only",
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (out / "instances.txt").write_text("\n".join(names) + "\n")
    jobs = []
    for repeat in range(args.repeats):
        for config, options in configs.items():
            for name in names:
                path = out / "jobs" / f"{config}--{name}--r{repeat}.json"
                path.write_text(json.dumps({
                    "instance": name, "config": config, "options": options,
                    "reference": references.get(name), "repeat": repeat,
                    "fixture_dir": str(args.fixture_dir.resolve()),
                    "data_dir": str(args.data_dir.resolve()),
                    "replay": str(args.replay.resolve()),
                    "output": str(out / "results" / path.name),
                    "seconds": args.sec, "memory_bytes": int(args.memory_gb * 1024**3),
                }) + "\n")
                jobs.append(str(path))
    random.Random(20261005).shuffle(jobs)
    print(f"Dispatching {len(jobs)} replays: {len(names)} fixtures, "
          f"{len(configs)} configs, {args.repeats} repeats, {args.jobs} jobs; {out}", flush=True)
    subprocess.run(
        ["parallel", "--jobs", str(args.jobs), "--halt", "soon,fail=1",
         "--joblog", str(out / "parallel.tsv"),
         sys.executable, str(SCRIPT), "--worker", "{}"],
        input="\n".join(jobs) + "\n", text=True, check=True,
    )
    results = summarize(out)
    if len(results) != len(jobs):
        raise RuntimeError(f"Expected {len(jobs)} results, found {len(results)}")
    return int(any(row["status"] in ("ERROR", "HARD_TIMEOUT") for row in results))


if __name__ == "__main__":
    sys.exit(main())
