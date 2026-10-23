#!/usr/bin/env python3
"""Exercise solver-facing SppGrasp controls and conditional FPump policy."""

import argparse
import pathlib
import re
import subprocess
import sys


HERE = pathlib.Path(__file__).resolve().parent
MODEL = HERE / "spp-grasp-policy.mps"
DISABLED = "SppGrasp: set partitioning with slacks - feasibility pump switched off"


def solve(cbc, name, options, model=MODEL):
    command = [
        cbc,
        str(model),
        "-log",
        "2",
        "-preprocess",
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
        "-sppGraspArtificial",
        "2",
        *options,
        "-sec",
        "5",
        "-solve",
        "-quit",
    ]
    result = subprocess.run(
        command,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=20,
    )
    if result.returncode:
        raise AssertionError(
            f"{name}: Cbc exited {result.returncode}\n{result.stdout}"
        )
    if "No match for" in result.stdout or "Unknown keyword" in result.stdout:
        raise AssertionError(f"{name}: a new parameter was rejected\n{result.stdout}")
    return result.stdout


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cbc", default=str(pathlib.Path.home() / "prog/cbc/bin/cbc")
    )
    args = parser.parse_args()
    cbc = str(pathlib.Path(args.cbc).resolve())

    auto = solve(cbc, "auto", ["-sppGraspFpump", "auto"])
    keep = solve(cbc, "on", ["-sppGraspFpump", "on"])
    after = solve(cbc, "after", ["-sppGraspFpump", "after"])
    explicit = solve(
        cbc, "explicit-feas", ["-sppGraspFpump", "after", "-feas", "on"]
    )
    no_improvement = solve(
        cbc, "no-spp-improvement",
        ["-sppGraspFpump", "after", "-cutoff", "0"],
    )
    gated = solve(
        cbc, "real-stall-gate",
        ["-sppGraspRealStallWork", "0.000001"],
    )
    root_spp_then_pump = solve(
        cbc, "root-spp-then-pump",
        ["-sppGraspRootSequence", "spp-fpump", "-feas", "on"],
    )
    root_spp_then_pump_fallback = solve(
        cbc, "root-spp-then-pump-if-no-incumbent",
        ["-sppGraspRootSequence", "spp-fpump-if-no-feasible", "-feas", "on"],
    )
    root_keep_original = solve(
        cbc, "root-sequence-keeping-original-spp",
        [
            "-sppGraspRootSequence", "spp-fpump-if-no-feasible",
            "-sppGraspRootKeepOriginal", "on", "-feas", "on",
        ],
    )
    root_mipstart_fallback = solve(
        cbc, "root-spp-then-pump-with-mipstart",
        [
            "-mipstart", str(HERE / "spp-grasp-policy.mst"),
            "-sppGraspRootSequence", "spp-fpump-if-no-feasible",
            "-feas", "on",
        ],
    )
    root_spp_then_pump_no_integer = solve(
        cbc, "root-spp-then-pump-without-integer-incumbent",
        ["-sppGraspRootSequence", "spp-fpump-if-no-feasible", "-feas", "on"],
        HERE / "spp-grasp-no-integer-policy.mps",
    )
    root_pump_then_spp = solve(
        cbc, "root-pump-then-spp",
        ["-sppGraspRootSequence", "fpump-spp", "-feas", "on"],
    )
    root_pump_guide_spp = solve(
        cbc, "root-pump-guide-then-spp",
        ["-sppGraspRootSequence", "fpump-spp-guide", "-feas", "on"],
    )
    root_explicit_pump_off = solve(
        cbc, "root-sequence-explicit-pump-off",
        ["-sppGraspRootSequence", "spp-fpump", "-feas", "off"],
    )

    if "SppGrasp: first 1, best 1" not in after:
        raise AssertionError("fixture did not produce an original-model SppGrasp solution")
    if DISABLED not in auto or DISABLED not in after:
        raise AssertionError("legacy auto/after-success policies did not disable FPump")
    if DISABLED in keep or DISABLED in explicit or DISABLED in no_improvement:
        raise AssertionError("FPump was disabled under a retain/override/failure policy")
    if "SppGrasp: progress gate stopped" not in gated:
        raise AssertionError("solver-facing real-stall gate was not activated")
    expected_sequences = (
        (root_spp_then_pump, "SppGrasp then FPump"),
        (
            root_spp_then_pump_fallback,
            "SppGrasp then FPump only without an incumbent",
        ),
        (
            root_pump_then_spp,
            "FPump then SppGrasp seeded from the incumbent",
        ),
        (
            root_pump_guide_spp,
            "FPump then SppGrasp with an optional less-fractional pump guide",
        ),
    )
    for output, description in expected_sequences:
        if not re.search(
            r"SppGrasp root sequence:\s+" + re.escape(description), output
        ):
            raise AssertionError(f"root sequence was not configured: {description}")
        if "SppGrasp not applicable:" in output:
            raise AssertionError(f"SppGrasp failed applicability: {description}")
    if "No match for" in root_explicit_pump_off:
        raise AssertionError("explicit -feas off failed to parse with root sequence")
    if "SppGrasp: first 1, best 1" not in root_explicit_pump_off:
        raise AssertionError("explicit -feas off suppressed the legacy SppGrasp call")
    if "calling feasibility pump (incumbent available)" not in root_spp_then_pump:
        raise AssertionError("unconditional sequence did not call FPump after an incumbent")
    if "Root sequence: completed feasibility pump (status " not in root_spp_then_pump:
        raise AssertionError("unconditional sequence did not report FPump execution time")
    if "Root sequence: completed SppGrasp (status " not in root_spp_then_pump:
        raise AssertionError("root sequence did not report SppGrasp execution time")
    if "skipped feasibility pump (incumbent already present)" not in root_spp_then_pump_fallback:
        raise AssertionError("conditional sequence did not skip FPump after a genuine incumbent")
    original_spp = root_keep_original.find("SppGrasp: first 1, best 1")
    ordered_root = root_keep_original.find("SppGrasp root sequence:")
    if original_spp < 0 or ordered_root < 0 or original_spp > ordered_root:
        raise AssertionError(
            "root keep-original did not run original-model SppGrasp before the ordered root calls"
        )
    if "Root sequence: calling feasibility pump" in root_keep_original:
        raise AssertionError("root keep-original unexpectedly invoked FPump after finding a solution")
    if "MIPStart provided solution with cost 2" not in root_mipstart_fallback:
        raise AssertionError("MIP start did not establish an incumbent before root heuristics")
    if "calling SppGrasp (incumbent available)" not in root_mipstart_fallback:
        raise AssertionError("SppGrasp was not called to improve the MIP-start incumbent")
    if "skipped feasibility pump (incumbent already present)" not in root_mipstart_fallback:
        raise AssertionError("conditional sequence did not honor the MIP-start incumbent")
    if "calling feasibility pump (no incumbent)" not in root_spp_then_pump_no_integer:
        raise AssertionError("conditional sequence did not run FPump after SppGrasp failed")
    print("SppGrasp solver parameter and FPump policy selftest: PASS")
    print("  auto preserves legacy native-SPP behavior")
    print("  on retains FPump; after disables only after a real SppGrasp solution")
    print("  explicit -feas on and no-improvement cases retain FPump")
    print("  root orderings call/skip FPump from genuine incumbent availability")
    print("  explicit -feas off retains the original-model SppGrasp fallback")


if __name__ == "__main__":
    try:
        main()
    except (AssertionError, subprocess.TimeoutExpired) as error:
        print(error, file=sys.stderr)
        sys.exit(1)
