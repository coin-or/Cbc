#!/usr/bin/env python3
"""The validator must retain objective constants in recomputation and BKS."""
import argparse
from pathlib import Path
import shutil
import subprocess

HERE = Path(__file__).resolve().parent
WORK = HERE / ".cbc-validate-sol-offset-selftest"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--validator", default=str(HERE / "cbc_validate_sol"))
    args = parser.parse_args()
    if WORK.exists():
        raise RuntimeError("refusing to overwrite " + str(WORK))
    WORK.mkdir()
    try:
        for cost in (2, -2):
            for offset in (0, 5, -5, 1000000):
                model = WORK / "constant.mps"
                model.write_text(
                    "NAME          CONSTANT\nROWS\n N  OBJ\n E  R0\n"
                    f"COLUMNS\n    X0        OBJ       {cost}\n"
                    "    X0        R0        1\n"
                    f"RHS\n    RHS1      OBJ       {offset}\n"
                    "    RHS1      R0        1\n"
                    "BOUNDS\n BV BND1      X0\nENDATA\n")
                solution = WORK / "constant.sol"
                objective = cost - offset
                solution.write_text(
                    f"Optimal - objective value {objective}\n0 X0 1 {cost}\n")
                command = [args.validator, "--expected-status", "optimal",
                           "--expected-obj", str(objective),
                           str(model), str(solution)]
                good = subprocess.run(command, capture_output=True, text=True,
                                      timeout=15)
                assert good.returncode == 0, (cost, offset, good.stdout, good.stderr)
                if offset:
                    bad = subprocess.run(
                        [*command[:-3], str(cost), str(model), str(solution)],
                        capture_output=True, text=True, timeout=15)
                    assert bad.returncode == 1, bad.stdout
                    solution.write_text(
                        f"Optimal - objective value {cost}\n0 X0 1 {cost}\n")
                    bad = subprocess.run(command, capture_output=True, text=True,
                                         timeout=15)
                    assert bad.returncode == 1, bad.stdout
        print("PASS: positive/negative/large objective constants and negative "
              "costs preserved; unshifted claims and references rejected")
    finally:
        shutil.rmtree(WORK)


if __name__ == "__main__":
    main()
