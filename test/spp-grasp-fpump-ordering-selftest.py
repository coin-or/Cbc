#!/usr/bin/env python3
"""Verify ordered FPump runs can do inner work with an incumbent present."""

import argparse
import pathlib
import re
import subprocess
import sys


HERE = pathlib.Path(__file__).resolve().parent
MODEL = HERE / "mip-sanity-data" / "mips" / "air04.mps.gz"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cbc", default=str(HERE.parent / "src" / "cbc"))
    parser.add_argument("--model", default=str(MODEL))
    args = parser.parse_args()
    model = pathlib.Path(args.model)
    if not model.is_file():
        raise AssertionError(f"FPump ordering fixture is unavailable: {model}")
    command = [
        str(pathlib.Path(args.cbc).resolve()),
        str(model),
        "-randomSeed",
        "271",
        "-randomCbcSeed",
        "271",
        "-threads",
        "1",
        "-sec",
        "20",
        "-log",
        "2",
        "-preprocess",
        "sos",
        "-sppGraspWork",
        "100",
        "-sppGraspRealStallWork",
        "50",
        "-sppGraspNoImprove",
        "1000",
        "-sppGrasp",
        "on",
        "-feas",
        "on",
        "-sppGraspFpump",
        "auto",
        "-sppGraspRootSequence",
        "spp-fpump",
        "-sppGraspRootKeepOriginal",
        "on",
        "-maxNodes",
        "10",
        "-solve",
        "-quit",
    ]
    result = subprocess.run(
        command,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=30,
    )
    if result.returncode:
        raise AssertionError(f"CBC exited {result.returncode}\n{result.stdout}")
    if "Ordered sequence FPump root placements: LC" not in result.stdout:
        raise AssertionError(f"Ordered sequence did not enable LC:\n{result.stdout}")
    incumbent_pump = False
    awaiting_incumbent_pump = False
    for line in result.stdout.splitlines():
        if "Root sequence: calling feasibility pump (incumbent available)" in line:
            awaiting_incumbent_pump = True
        elif "FPump execution:" in line and awaiting_incumbent_pump:
            match = re.search(r"LP-solves=(\d+) passes=(\d+)", line)
            if match and int(match.group(1)) > 0 and int(match.group(2)) > 0:
                incumbent_pump = True
            awaiting_incumbent_pump = False
    if not incumbent_pump:
        raise AssertionError(
            "FPump was not shown completing LP solves and passes with an "
            f"incumbent present:\n{result.stdout}"
        )
    print("SppGrasp/FPump ordering selftest: PASS")
    print("  post-incumbent FPump reached LP solves and completed pump passes")


if __name__ == "__main__":
    try:
        main()
    except (AssertionError, subprocess.TimeoutExpired) as error:
        print(error, file=sys.stderr)
        sys.exit(1)
