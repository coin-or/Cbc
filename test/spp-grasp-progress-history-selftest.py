#!/usr/bin/env python3
"""Check pre-search incumbents in the common-origin progress CSV."""

import argparse
import csv
import pathlib
import re
import subprocess
import sys


HERE = pathlib.Path(__file__).resolve().parent
MODEL = HERE / "spp-grasp-progress-max-offset.mps"
STARTS = (
    HERE / "spp-grasp-progress.mst",
    HERE / "spp-grasp-progress.before.mst",
)
TRACE = HERE / "spp-grasp-progress-history-selftest.csv"
EMPTY_TRACE = HERE / "spp-grasp-progress-history-selftest-empty.csv"
EMPTY_LOG = HERE / "spp-grasp-progress-history-selftest-empty.log"
PARTIAL_TRACE = HERE / "spp-grasp-progress-history-selftest-partial.csv"
PARTIAL_LOG = HERE / "spp-grasp-progress-history-selftest-partial.log"
INFEASIBLE_TRACE = HERE / "spp-grasp-progress-history-selftest-infeasible.csv"
INFEASIBLE_LOG = HERE / "spp-grasp-progress-history-selftest-infeasible.log"
sys.path.insert(0, str(HERE))
import cbc_pdi  # noqa: E402


def run(cbc, mipstart=None):
    TRACE.unlink(missing_ok=True)
    command = [
        str(cbc),
        str(MODEL),
        "-maximize",
        "-timeMode",
        "elapsed",
        "-threads",
        "1",
        "-preprocess",
        "off",
        "-presolve",
        "off",
        "-cuts",
        "off",
        "-progressCsv",
        str(TRACE),
        "-sec",
        "5",
        "-solve",
        "-quit",
    ]
    if mipstart is not None:
        command[3:3] = ["-mipstart", str(mipstart)]
        command[3:3] = ["-sppGrasp", "off"]
    else:
        command[3:3] = ["-sppGrasp", "on"]
    result = subprocess.run(
        command,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=15,
    )
    if result.returncode:
        raise AssertionError(f"CBC exited {result.returncode}\n{result.stdout}")
    if "Objective value:                8" not in result.stdout:
        raise AssertionError(f"Maximization/objective offset changed:\n{result.stdout}")
    with TRACE.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    final_rows = [row for row in rows if row["source"] == "final"]
    if mipstart is None and len(final_rows) != 1:
        raise AssertionError(f"Expected one final progress event, got {final_rows!r}")
    if len(final_rows) > 1:
        raise AssertionError(f"Expected at most one final progress event, got {final_rows!r}")
    if mipstart is None:
        original_events = [
            row for row in rows if row["source"] == "SppGraspOriginal"
        ]
        root_events = [row for row in rows if row["source"] == "SppGrasp"]
        if not original_events or not root_events:
            raise AssertionError(
                f"Expected original-model and root SppGrasp events: {rows!r}"
            )
        event = original_events[0]
        if float(event["primal"]) != 2.0 or float(event["dual"]) > -1.0e49:
            raise AssertionError(f"Original-model SppGrasp event is incorrect: {event!r}")
    else:
        if "MIPStart ..." not in result.stdout or "Solution with cost 8" not in result.stdout:
            raise AssertionError(f"MIP start was not accepted:\n{result.stdout}")
        if "found by Reduced search" in result.stdout:
            raise AssertionError(f"Duplicate MIP-start incumbent message:\n{result.stdout}")
        if not re.search(r"\n\n[^\n]*MIPStart \.\.\.", result.stdout):
            raise AssertionError(f"Missing blank line before MIPStart:\n{result.stdout}")
        events = [row for row in rows if row["source"] == "mipstart"]
        if len(events) != 1:
            raise AssertionError(f"Expected one MIP-start event, got {events!r}")
        event = events[0]
    expected = 2.0 if mipstart is None else 8.0
    if float(event["primal"]) != expected or float(event["dual"]) > -1.0e49:
        raise AssertionError(f"Incorrect true objective or unknown-bound sentinel: {event!r}")
    if float(event["time"]) <= 0.0:
        raise AssertionError(
            f"Pre-search incumbent was assigned a zero/fabricated origin: {event!r}"
        )
    times = [float(row["time"]) for row in rows]
    if times != sorted(times):
        raise AssertionError(f"Progress events are not in global-time order: {times!r}")
    if rows.index(event) >= next(
        (i for i, row in enumerate(rows) if row["source"] in ("rootLP", "tree", "node")),
        len(rows),
    ):
        raise AssertionError("The pre-search MIP start was not recorded before root progress")
    if cbc_pdi.gap_value(float(event["primal"]), None) != 1.0:
        raise AssertionError("An unknown dual bound must count as gap 1, not as proof")
    parsed = {
        "prim": [],
        "bnb": [],
        "root_lp": None,
        "cuts": [],
        "final": {"kind": "Optimal", "time": None, "obj": None, "bound": None},
        "infeasible": False,
    }
    cbc_pdi.apply_progress_csv(parsed, TRACE)
    if final_rows:
        final = final_rows[0]
        if parsed["final"]["time"] != float(final["time"]):
            raise AssertionError("PDI final time did not use the common-origin CSV clock")
        if parsed["final"]["obj"] != float(final["primal"]):
            raise AssertionError("PDI final primal did not use the final progress event")
        if parsed["final"]["bound"] != float(final["dual"]):
            raise AssertionError("PDI final dual did not use the final progress event")
    total = re.search(r"Total time \(Wallclock seconds\):\s*([0-9.eE+-]+)", result.stdout)
    if total and times[-1] > float(total.group(1)) + 0.01:
        raise AssertionError("Progress clock exceeds the whole-solve wall-clock duration")


def check_empty_progress_no_incumbent_fallback():
    EMPTY_TRACE.write_text("")
    EMPTY_LOG.write_text(
        '# command: ["cbc", "-threads", "1"]\n'
        "✔ LP Stopped — Frac: 795/174997 (0.5%) Obj: 89.4439 "
        "Iters: 92974 Time: 120s\n"
        "Result - Stopped on time limit\n"
        "No feasible solution found\n"
        "Total time (Wallclock seconds): 120.1\n"
    )
    try:
        events = cbc_pdi.run_events(EMPTY_LOG, EMPTY_TRACE, "csv")
        if events["event_source"] != "timestampedLog (no progress CSV)":
            raise AssertionError(f"Empty progress source is mislabeled: {events!r}")
        trace = cbc_pdi.build_trace(events, False)
        metrics = cbc_pdi.integrate(trace, 120.0)
        if metrics["pdi"] != 2.0 or metrics["primal"] is not None:
            raise AssertionError(
                f"No-incumbent full horizon must have PDI 2 without a primal: {metrics!r}"
            )
    finally:
        EMPTY_TRACE.unlink(missing_ok=True)
        EMPTY_LOG.unlink(missing_ok=True)


def check_start_only_csv_falls_back_to_log():
    PARTIAL_TRACE.write_text(
        "time,primal,dual,source,nodes\n"
        "0.000,1e+50,-1e+50,start,0\n"
    )
    PARTIAL_LOG.write_text(
        '# command: ["cbc", "-threads", "1"]\n'
        "Result - Optimal solution found\n"
        "Objective value: 7\n"
        "Total time (Wallclock seconds): 0.1\n"
    )
    try:
        events = cbc_pdi.run_events(PARTIAL_LOG, PARTIAL_TRACE, "csv")
        if events["event_source"] != (
            "timestampedLog (progress CSV missing final event)"
        ):
            raise AssertionError(f"Start-only progress source is mislabeled: {events!r}")
        trace = cbc_pdi.build_trace(events, False)
        metrics = cbc_pdi.integrate(trace, 120.0)
        if trace.status != "optimal" or metrics["primal"] != 7.0:
            raise AssertionError(f"Early proof was lost in log fallback: {metrics!r}")
        if abs(metrics["pdi"] - (2.0 * 0.1 / 120.0)) > 1e-12:
            raise AssertionError(
                f"Early presolve proof did not close the gap at log end time: {metrics!r}"
            )
        PARTIAL_LOG.write_text(
            '# command: ["cbc", "-threads", "1"]\n'
            "Result - Optimal solution found\n"
            "Objective value: 7\n"
        )
        events = cbc_pdi.run_events(PARTIAL_LOG, PARTIAL_TRACE, "csv")
        trace = cbc_pdi.build_trace(events, False)
        metrics = cbc_pdi.integrate(trace, 120.0)
        if trace.primal or trace.dual or metrics["pdi"] is not None \
                or metrics["primal"] != 7.0:
            raise AssertionError(
                f"An untimed proof must not fabricate a time-zero event or PDI: {metrics!r}"
            )
        PARTIAL_LOG.write_text(
            '# command: ["cbc", "-threads", "1"]\n'
            "Result - Optimal solution found\n"
            "Objective value: 7\n"
            "Total time (Wallclock seconds): 0.000\n"
        )
        events = cbc_pdi.run_events(PARTIAL_LOG, PARTIAL_TRACE, "csv")
        trace = cbc_pdi.build_trace(events, False)
        metrics = cbc_pdi.integrate(trace, 120.0)
        if trace.primal or trace.dual or metrics["pdi"] is not None \
                or metrics["primal"] != 7.0:
            raise AssertionError(
                f"A rounded-zero proof time must not be treated as an exact time-zero event: {metrics!r}"
            )
    finally:
        PARTIAL_TRACE.unlink(missing_ok=True)
        PARTIAL_LOG.unlink(missing_ok=True)


def check_infeasible_proof_closes_pdi_gap():
    INFEASIBLE_TRACE.write_text("")
    try:
        cases = (
            ("Problem is infeasible - 0.284917 seconds\n"
             "Total time (Wallclock seconds): 0.306276\n", 0.306276),
            ("Bound tightening: INFEASIBLE in round 2 — row 4423, 0.014 s.\n"
             "Bound tightening: infeasibility proved — skipping solve.\n"
             "Total time (Wallclock seconds): 1.17179\n", 1.17179),
        )
        for output, end_time in cases:
            INFEASIBLE_LOG.write_text(
                '# command: ["cbc", "-threads", "1"]\n' + output
            )
            events = cbc_pdi.run_events(INFEASIBLE_LOG, INFEASIBLE_TRACE, "csv")
            trace = cbc_pdi.build_trace(events, False)
            metrics = cbc_pdi.integrate(trace, 120.0)
            if events["event_source"] != "timestampedLog (no progress CSV)":
                raise AssertionError(f"Early infeasibility log source is mislabeled: {events!r}")
            if trace.status != "infeasible":
                raise AssertionError(f"Infeasibility proof was not recognized: {trace!r}")
            if abs(metrics["pdi"] - (2.0 * end_time / 120.0)) > 1e-12 \
                    or metrics["final_gap"] != 0.0:
                raise AssertionError(
                    f"Infeasible proof must count gap 2 until proof and 0 thereafter: {metrics!r}"
                )
            if metrics["primal"] is not None:
                raise AssertionError(f"Proof of infeasibility fabricated a primal: {metrics!r}")
        INFEASIBLE_LOG.write_text(
            '# command: ["cbc", "-threads", "1"]\n'
            "Problem is infeasible - 0.1 seconds\n"
            "Result - Stopped on time limit\n"
            "No feasible solution found\n"
            "Total time (Wallclock seconds): 120\n"
        )
        events = cbc_pdi.run_events(INFEASIBLE_LOG, INFEASIBLE_TRACE, "csv")
        trace = cbc_pdi.build_trace(events, False)
        metrics = cbc_pdi.integrate(trace, 120.0)
        if trace.status != "stopped" or metrics["pdi"] != 2.0 \
                or metrics["final_gap"] != 2.0:
            raise AssertionError(
                f"An intermediate infeasible relaxation must not become an infeasibility proof: {metrics!r}"
            )
    finally:
        INFEASIBLE_TRACE.unlink(missing_ok=True)
        INFEASIBLE_LOG.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cbc", default=str(HERE.parent / "src" / "cbc"))
    args = parser.parse_args()
    cbc = pathlib.Path(args.cbc).resolve()
    try:
        for mipstart in STARTS:
            run(cbc, mipstart)
        run(cbc)
        check_empty_progress_no_incumbent_fallback()
        check_start_only_csv_falls_back_to_log()
        check_infeasible_proof_closes_pdi_gap()
    finally:
        TRACE.unlink(missing_ok=True)
        EMPTY_TRACE.unlink(missing_ok=True)
        EMPTY_LOG.unlink(missing_ok=True)
        PARTIAL_TRACE.unlink(missing_ok=True)
        PARTIAL_LOG.unlink(missing_ok=True)
        INFEASIBLE_TRACE.unlink(missing_ok=True)
        INFEASIBLE_LOG.unlink(missing_ok=True)
    print("SppGrasp progress-history selftest: PASS")
    print("  standard/.before MIP starts and original-model SppGrasp precede root progress")
    print("  unknown dual remains gap 1; maximization and the MPS objective offset are preserved")
    print("  empty/start-only CSVs use an explicitly labeled log fallback")
    print("  an infeasibility proof closes PDI at proof time without inventing a primal")
    print("  missing/rounded-zero proof times do not fabricate time-zero events or PDI values")


if __name__ == "__main__":
    try:
        main()
    except (AssertionError, subprocess.TimeoutExpired) as error:
        print(error, file=sys.stderr)
        sys.exit(1)
