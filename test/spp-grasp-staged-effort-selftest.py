#!/usr/bin/env python3
"""Focused checks for the opt-in SppGrasp staged-work policy."""

import argparse
import pathlib
import re
import subprocess
import sys


HERE = pathlib.Path(__file__).resolve().parent
MODEL = HERE / "spp-grasp-stage-policy.mps"
AIR04 = HERE / "mip-sanity-data" / "mips" / "air04.mps.gz"


def solve(cbc, model, extra, timeout=20):
    command = [
        str(cbc),
        str(model),
        "-log",
        "2",
        "-preprocess",
        "off",
        "-presolve",
        "off",
        "-cuts",
        "off",
        "-sppGraspDestroyRows",
        "2",
        "-sppGraspStall",
        "0",
        "-sppGraspDestroyGrow",
        "0.1",
        "-sppGraspRelatedRows",
        "0.75",
        "-sppGraspRealStallWork",
        "0",
        *extra,
        "-sec",
        "5",
        "-maxNodes",
        "100",
        "-solve",
        "-quit",
    ]
    result = subprocess.run(
        command,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=timeout,
    )
    if result.returncode:
        raise AssertionError(
            f"CBC exited {result.returncode}\ncommand: {command}\n{result.stdout}"
        )
    if "No match for" in result.stdout or "Unknown keyword" in result.stdout:
        raise AssertionError(f"A staged-effort option was rejected:\n{result.stdout}")
    return result.stdout


def policy(cbc, sequence, work="0.01", stage_parameter=True, fpump=False):
    options = [
        "-feas",
        "on" if fpump else "off",
        "-sppGraspRootSequence",
        sequence,
        "-sppGraspWork",
        work,
    ]
    if stage_parameter:
        options.extend(["-sppGraspStageFirstPercent", "25"])
    return solve(
        cbc,
        MODEL,
        options,
    )


def spp_result(output):
    match = re.search(
        r"SppGrasp: first ([^,]+), best ([^ ]+) .* after (\d+) LNS iterations",
        output,
    )
    return match.groups() if match else None


def full_air04(cbc):
    if not AIR04.is_file():
        raise AssertionError(f"Missing optional sanity fixture: {AIR04}")
    command = [
        str(cbc),
        str(AIR04),
        "-randomSeed",
        "271",
        "-randomCbcSeed",
        "271",
        "-threads",
        "1",
        "-sec",
        "30",
        "-log",
        "2",
        "-preprocess",
        "sos",
        "-sppGraspDestroyRows",
        "5",
        "-sppGraspStall",
        "200",
        "-sppGraspDestroyGrow",
        "0",
        "-sppGraspRelatedRows",
        "0",
        "-sppGraspArtificial",
        "1.1",
        "-sppGraspWork",
        "500",
        "-sppGraspNoImprove",
        "1000",
        "-sppGraspRealStallWork",
        "0",
        "-sppGrasp",
        "on",
        "-feas",
        "on",
        "-sppGraspFpump",
        "on",
        "-sppGraspRootSequence",
        "spp-fpump-staged",
        "-sppGraspStageFirstPercent",
        "25",
        "-maxNodes",
        "1",
        "-solve",
        "-quit",
    ]
    result = subprocess.run(
        command,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=45,
    )
    if result.returncode:
        raise AssertionError(
            f"Air04 integration run exited {result.returncode}\n{result.stdout}"
        )
    first = re.search(r"stage 1 used ([0-9.]+) of ([0-9.]+) work", result.stdout)
    second = re.search(r"stage 2 used ([0-9.]+) of ([0-9.]+) total work", result.stdout)
    if not first or not second:
        raise AssertionError(f"Both staged calls must be exercised:\n{result.stdout}")
    first_used, total = map(float, first.groups())
    second_used, second_total = map(float, second.groups())
    if total != second_total:
        raise AssertionError("The two stages reported inconsistent total budgets")
    if first_used + second_used > total + 5.0e6:
        raise AssertionError(
            f"Staged counted work exceeded its total by more than one bounded "
            f"iteration: {first_used:g} + {second_used:g} > {total:g}"
        )
    return first_used, second_used, total, result.stdout


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cbc", default=str(HERE.parent / "src" / "cbc")
    )
    parser.add_argument("--baseline", help="optional frozen legacy CBC binary")
    parser.add_argument(
        "--air04",
        action="store_true",
        help="also run the 30-second root-budget integration check",
    )
    args = parser.parse_args()
    cbc = pathlib.Path(args.cbc).resolve()

    staged = policy(cbc, "spp-fpump-staged")
    if "SppGrasp stage 1 used" not in staged:
        raise AssertionError("the initial staged SppGrasp pass did not run")
    if "staged SppGrasp with a post-cuts improvement pass" not in staged:
        raise AssertionError("the staged root policy was not configured")
    if "Root heuristics (after cuts)" not in staged or "calling SppGrasp" not in staged:
        raise AssertionError("the post-cuts SppGrasp hook did not execute")
    if "SppGrasp stage 2 used" not in staged:
        raise AssertionError("the staged second-pass accounting was not reported")
    if "Result - Optimal solution found" not in staged:
        raise AssertionError("the staged policy changed the fixture's valid optimum")

    for sequence in ("spp-fpump-upfront", "spp-fpump-short"):
        output = policy(cbc, sequence)
        if "SppGrasp initial effort pass used" not in output:
            raise AssertionError(f"{sequence} did not enforce its initial work cap")
        if "Root sequence: calling SppGrasp" in output:
            raise AssertionError(f"{sequence} unexpectedly registered a root SppGrasp pass")
        with_fpump = policy(cbc, sequence, fpump=True)
        if "Heuristic feasibility pump:" not in with_fpump:
            raise AssertionError(f"{sequence} unexpectedly suppressed the feasibility pump")
        if "Root sequence: calling SppGrasp" in with_fpump:
            raise AssertionError(f"{sequence} registered a root SppGrasp pass with FPump on")

    legacy = policy(cbc, "legacy")
    if "SppGrasp stage " in legacy:
        raise AssertionError("legacy mode unexpectedly entered staged accounting")
    if spp_result(legacy) is None:
        raise AssertionError("legacy SppGrasp did not produce its expected fixture result")
    if args.baseline:
        baseline = policy(
            pathlib.Path(args.baseline).resolve(), "legacy", stage_parameter=False
        )
        if spp_result(legacy) != spp_result(baseline):
            raise AssertionError(
                "default-off legacy SppGrasp trajectory differs from the frozen "
                f"baseline: {spp_result(legacy)} != {spp_result(baseline)}"
            )

    for sequence in (
        "spp-fpump-staged",
        "spp-fpump-upfront",
        "spp-fpump-short",
    ):
        unlimited = policy(cbc, sequence, "0")
        if "requires a positive sppGraspWork limit" not in unlimited:
            raise AssertionError(f"{sequence} accepted an unbounded work setting")
        if "SppGrasp stage 1 used" in unlimited or "initial effort pass used" in unlimited:
            raise AssertionError(f"{sequence} incorrectly split an unlimited run")

    print("SppGrasp staged-effort selftest: PASS")
    print("  staged mode calls the original and post-cuts hooks; legacy remains default-off")
    print("  an unbounded sppGraspWork setting falls back without splitting the run")
    if args.baseline:
        print("  legacy construction/objective/iteration trajectory matches frozen CBC")
    if args.air04:
        first, second, total, output = full_air04(cbc)
        print(
            f"  air04 actual deterministic work: stage 1 {first:.0f}, "
            f"stage 2 {second:.0f}, total {total:.0f}"
        )
        if "SppGrasp staged incumbent map rejected" in output:
            print("  incompatible root incumbent was rejected rather than used as an SPP start")


if __name__ == "__main__":
    try:
        main()
    except (AssertionError, subprocess.TimeoutExpired) as error:
        print(error, file=sys.stderr)
        sys.exit(1)
