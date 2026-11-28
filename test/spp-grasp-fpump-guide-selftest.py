#!/usr/bin/env python3
"""Compile and run the focused FPump-to-SppGrasp guide-selection checks."""
import os
from pathlib import Path
import shlex
import subprocess

test_dir = Path(__file__).resolve().parent
source_dir = test_dir.parent / "src"
binary = test_dir / ".spp-grasp-fpump-guide-selftest"
compiler = shlex.split(os.environ.get("CXX", "c++"))

try:
    subprocess.run(
        compiler + [
            "-std=c++11", "-Wall", "-Wextra", "-pedantic",
            "-I", str(source_dir),
            str(test_dir / "spp-grasp-fpump-guide-selftest.cpp"),
            "-o", str(binary),
        ],
        check=True,
    )
    subprocess.run([str(binary)], check=True)
    print("PASS: pump-guide selection validates quality, bounds, finiteness, and mapping")
finally:
    binary.unlink(missing_ok=True)
