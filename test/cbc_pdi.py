#!/usr/bin/env python3
"""Primal-dual gap integral (PDI) for Cbc runs: how early and how far a run
closes the gap.

Sub-commands
  run      run cbc on instances, one log per instance (parallel)
  analyze  parse logs -> <dir>/pdi.tsv (one row per instance)
  compare  rank N>=2 experiment directories on the PDI metrics
  plot     mean gap-vs-time curves of experiments (needs matplotlib)
  explain  one figure per run(s): bounds, gap(t) and the shaded PDI area
  credit   which method (cuts, heuristics, tree) produced the improvements

Metric (all minimisation internally; max problems are negated)
  gap(t) = 2                      no incumbent yet            (worst case)
         = 1                      incumbent but no dual bound, or p*d <= 0
         = |p-d| / max(|p|,|d|)   otherwise, in [0,1]; 0 when p == d
  PDI    = (1/T) * integral_0^T gap(t) dt      in [0,2]; MINIMISE it
  T is the time horizon (default: the run's time limit; --sec overrides it).
  After a proven-optimal finish the gap stays 0; after an early stop (node
  limit, hard kill) the last gap is held until T, and a run longer than T is
  truncated at T. To compare runs with different time limits, analyze all of
  them with the same --sec. Besides the combined
  PDI the tool reports, when a best-known objective (bks.tsv) is available,
  the primal integral and dual integral against it (each in [0,1] with 1 =
  no incumbent / no bound), and t_first, t_gap10 and t_gap1 (time to first
  incumbent, to 10% and to 1% gap).

Parsed log events (Cbc "next" console format; exact times, 6 significant
digits): root LP objective, cut-pass objectives, every starred incumbent
(feasibility pump, heuristics, B&B), B&B progress rows (incumbent + bound),
and the final summary line. Incumbents printed without a timestamp in the
root-heuristics table are placed at the end of the preceding timed row, which
is approximate (a few ms to a few hundred ms).

For single-threaded runs, `auto` uses the complete progress CSV, including
pre-search incumbents and its common-origin final timestamp. For runs with
multiple CBC threads, `auto` uses the timestamped log because worker clones
may omit tree events from the CSV; the output `events` column labels that
different, approximate source. Do not pool those sources as if their clocks
were identical.

Examples
  cbc_pdi.py run --out /tmp/base --sec 60 'bpc_*' 'cvrp_*'
  cbc_pdi.py run --out /tmp/new  --sec 60 --cbc ~/prog/cbc-new/bin/cbc 'bpc_*' 'cvrp_*'
  cbc_pdi.py analyze /tmp/base /tmp/new
  cbc_pdi.py compare base=/tmp/base new=/tmp/new
  cbc_pdi.py plot --out gap.png base=/tmp/base new=/tmp/new
"""
import argparse
import concurrent.futures as cf
import csv
import fnmatch
import json
import math
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "mip-sanity-data")
INF = float("inf")
NOGAP = 2.0

# ---------------------------------------------------------------- parsing


def num(tok):
    """Float from a log token; None for '-', '—', '' or unparsable."""
    try:
        v = float(tok)
    except (TypeError, ValueError):
        return None
    if v >= 1e49 or v <= -1e49:
        return None
    return v


class Trace:
    """Step functions of primal and dual bound (minimisation sign)."""

    def __init__(self):
        self.primal = []   # (t, value), improving
        self.dual = []     # (t, value), improving
        self.end_time = None
        self.status = "unknown"   # optimal | stopped | infeasible | unknown
        self.final_primal = None
        self.final_dual = None
        self.maximize = None

    def add_primal(self, t, v):
        if v is None:
            return
        if not self.primal or v < self.primal[-1][1] - 1e-12 * max(1, abs(v)):
            self.primal.append((t, v))

    def add_dual(self, t, v):
        if v is None:
            return
        if not self.dual or v > self.dual[-1][1] + 1e-12 * max(1, abs(v)):
            self.dual.append((t, v))


def parse_log(path):
    """Parse a Cbc log into raw events (values in the log's own sign)."""
    prim, dual = [], []          # (t, v)
    section = ""
    root_lp = None               # (t, obj)
    cut_passes = []              # (t, obj)
    bnb_rows = []                # (t, best, bound)
    cursor = 0.0                 # running time estimate for untimed rows
    final = {}
    res = {}
    infeasible = False
    with open(path, errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n")
            star = "\u2605" in line or re.match(r"^\s*\*\s", line) is not None
            line = line.replace("\u2605", " ").replace("\u2502", " ").replace("|", " ")
            if re.match(r"^\s*\*\s", raw):
                line = re.sub(r"^\s*\*", " ", line)
            s = line.strip()
            if s.startswith("\u25b6") or re.match(r"^>\s", s):
                section = s.lstrip("\u25b6> ").strip()
                continue
            if s.startswith("\u2714") or s.startswith("\u2718"):
                m = re.search(r"LP Optimal.*Obj:\s*(\S+).*Time:\s*([\d.eE+-]+)s", s)
                if m and root_lp is None:
                    root_lp = (float(m.group(2)), num(m.group(1)))
                    cursor = root_lp[0]
                m = re.search(r"best\s+\S+\s+in\s+([\d.eE+-]+)s", s)
                if m:
                    cursor += float(m.group(1))
                m = re.search(r"Preprocessing complete.*Time:\s*([\d.eE+-]+)s", s)
                if m:
                    cursor = max(cursor, float(m.group(1)))
                m = re.match(r"\S+\s+(Optimal|Stopped|Infeasible)\b(.*)", s)
                if m and "Time:" in s and ("Obj:" in s or "Bound:" in s or "BestSol" in s
                                           or "No integer" in s or "Infeasible" in s):
                    final["kind"] = m.group(1)
                    tm = re.search(r"Time:\s*([\d.eE+-]+)s", s)
                    final["time"] = float(tm.group(1)) if tm else None
                    for key, pat in (("obj", r"(?:Obj|BestSol):\s*(\S+)"), ("bound", r"Bound:\s*(\S+)")):
                        mm = re.search(pat, s)
                        final[key] = num(mm.group(1)) if mm else None
                    if m.group(1) == "Infeasible":
                        infeasible = True
                continue
            if s.startswith("Result - "):
                res["kind"] = ("Optimal" if "Optimal solution" in s else
                               "Infeasible" if "infeasible" in s.lower() else "Stopped")
                if res["kind"] == "Infeasible":
                    infeasible = True
            elif s.lower().startswith("problem is infeasible"):
                res["kind"] = "Infeasible"
                infeasible = True
                m = re.search(r"-\s*([\d.eE+-]+)\s*seconds", s)
                if m:
                    res["time"] = num(m.group(1))
            elif "bound tightening: infeasibility proved" in s.lower():
                res["kind"] = "Infeasible"
                infeasible = True
            m = re.match(r"Objective value:\s*(\S+)", s)
            if m:
                res["obj"] = num(m.group(1))
            m = re.match(r"Lower bound:\s*(\S+)", s)
            if m:
                res["bound"] = num(m.group(1))
            m = re.match(r"Total time \(Wallclock seconds\):\s*(\S+)", s)
            if m:
                res["time"] = num(m.group(1))
            toks = s.split()
            heur_row = (len(toks) >= 3 and toks[0][:1].isalpha() and num(toks[-1]) is not None
                        and "solution" in s and section.startswith(("Root node", "Feasibility"))
                        and not toks[0].startswith(("Round", "Pass", "Nodes", "Heuristic", "Effort", "Phase")))
            if not toks or not re.match(r"^[\d.]", toks[0]) and not toks[0][0].isalpha():
                continue

            if section.startswith("Cut generation"):
                # Pass Rows Tight Frac Suminf Objective Time
                if len(toks) == 7 and toks[0].isdigit():
                    t, o = num(toks[-1]), num(toks[-2])
                    if t is not None and o is not None:
                        cut_passes.append((t, o))
                        cursor = max(cursor, t)
            elif section.startswith("Feasibility pump"):
                # Round Pass Frac Suminf Status BestSol Time
                if len(toks) >= 7 and toks[0].isdigit():
                    t = num(toks[-1])
                    if t is not None:
                        cursor = max(cursor, t)
                        v = num(toks[-2])
                        if star and v is not None:
                            prim.append((t, v))
            elif heur_row and not section.startswith(("Branch and bound", "Cut generation")):
                # Heuristic table row ("Name  status  BestSol  Time"); the table
                # continues under the Feasibility Jump header in the log.
                d = num(toks[-1])
                cursor += d
                v = num(toks[-2])
                if v is not None:
                    prim.append((cursor, v))
            elif star and not section.startswith(("Branch and bound", "Root node heuristics")):
                # Starred incumbent in any other table (Feasibility Jump,
                # proximity search, ...): no timestamp, so place it at the
                # running time cursor and take the last numeric column.
                nums = [num(x) for x in toks if num(x) is not None]
                if nums:
                    prim.append((cursor, nums[-1]))
            elif section.startswith("Root node heuristics"):
                # Name... Status BestSol Time   (Time = duration of the heuristic)
                d = num(toks[-1])
                if d is not None and len(toks) >= 3:
                    cursor += d
                    v = num(toks[-2])
                    if v is not None and "solution" not in toks[-2]:
                        prim.append((cursor, v))
            elif section.startswith("Branch and bound"):
                # Nodes OnTree Depth BestSol [Method] BestBound Gap% Time
                if len(toks) >= 7 and toks[0].isdigit() and toks[1].isdigit() and toks[2].isdigit():
                    t = num(toks[-1])
                    if t is None:
                        continue
                    best, bound = num(toks[3]), num(toks[-3])
                    bnb_rows.append((t, best, bound))
                    if best is not None:
                        prim.append((t, best))
                    cursor = max(cursor, t)
    if not final.get("kind") and res.get("kind"):
        # no B&B summary (solved before the tree): use the closing Result block
        final = dict(res)
        if final["kind"] == "Optimal" and final.get("obj") is not None:
            final["bound"] = final["obj"]
    return dict(prim=prim, root_lp=root_lp, cuts=cut_passes, bnb=bnb_rows,
                final=final, infeasible=infeasible)


def apply_progress_csv(ev, path):
    """Replace log events by progress CSV times, including the global final clock."""
    prim, bnb = [], []
    final_event = None
    with open(path) as f:
        for r in csv.DictReader(f):
            t, p, d = float(r["time"]), float(r["primal"]), float(r["dual"])
            pv = p if p < 1e49 else None
            dv = d if d > -1e49 else None
            if r["source"] == "final":
                final_event = (t, pv, dv)
            if pv is not None:
                prim.append((t, pv))
            if dv is not None:
                bnb.append((t, pv, dv))
    if final_event is None:
        ev["event_source"] = "timestampedLog (progress CSV missing final event)"
        return ev
    ev["prim"], ev["bnb"], ev["root_lp"], ev["cuts"] = prim, bnb, None, []
    ev["event_source"] = "progressCsv"
    t, primal, dual = final_event
    ev["final"]["time"] = t
    ev["final"]["obj"] = primal
    ev["final"]["bound"] = dual
    return ev


def build_trace(ev, sense_max):
    """Combine raw events into monotone primal/dual step functions (min sign)."""
    sg = -1.0 if sense_max else 1.0
    tr = Trace()
    tr.maximize = sense_max
    final = ev["final"]
    for t, v in sorted(ev["prim"]):
        tr.add_primal(t, sg * v)
    # dual bound: root LP, cut passes (offset-calibrated), B&B rows
    bb_dual = [(t, sg * b) for t, _, b in ev["bnb"] if b is not None]
    cuts = [(t, sg * o) for t, o in ev["cuts"]]
    if cuts:
        # Cut-pass objectives live in the preprocessed space: shift them so
        # that the last pass matches the first B&B bound (or root LP if no tree).
        anchor = bb_dual[0][1] if bb_dual else (sg * ev["root_lp"][1] if ev["root_lp"] and ev["root_lp"][1] is not None else None)
        off = (anchor - cuts[-1][1]) if anchor is not None else 0.0
        cuts = [(t, o + off) for t, o in cuts]
    pts = []
    if ev["root_lp"] and ev["root_lp"][1] is not None:
        pts.append((ev["root_lp"][0], sg * ev["root_lp"][1]))
    pts += cuts + bb_dual
    for t, v in sorted(pts):
        tr.add_dual(t, v)
    tr.end_time = final.get("time")
    kind = final.get("kind")
    if kind == "Infeasible" or (ev["infeasible"] and kind is None):
        tr.status = "infeasible"
    elif kind == "Optimal":
        tr.status = "optimal"
    elif kind == "Stopped":
        tr.status = "stopped"
    tr.final_primal = None if final.get("obj") is None else sg * final["obj"]
    tr.final_dual = None if final.get("bound") is None else sg * final["bound"]
    if tr.status == "stopped" and tr.final_dual is not None and tr.end_time:
        tr.add_dual(tr.end_time, tr.final_dual)
    if tr.status == "optimal" and tr.final_primal is not None \
            and tr.end_time is not None and tr.end_time > 0.0:
        tr.add_primal(tr.end_time, tr.final_primal)
        # a proof closes the gap exactly at the end
        tr.dual.append((tr.end_time, tr.final_primal))
        tr.dual.sort()
    elif tr.final_dual is not None and tr.dual and tr.dual[-1][1] > tr.final_dual:
        # never claim a bound stronger than the one Cbc ends with
        tr.dual = [(t, min(v, tr.final_dual)) for t, v in tr.dual]
    return tr


# ---------------------------------------------------------------- metric


def gap_value(p, d):
    if p is None:
        return NOGAP
    if d is None or d == -INF:
        return 1.0
    if p == d:
        return 0.0
    if p * d <= 0:
        return 1.0
    return min(1.0, abs(p - d) / max(abs(p), abs(d)))


def one_sided(v, ref):
    """Gap of a value to a reference (best-known) objective; 1 if unknown."""
    if v is None:
        return 1.0
    if v == ref:
        return 0.0
    if v * ref <= 0:
        return 1.0
    return min(1.0, abs(v - ref) / max(abs(v), abs(ref)))


def step_at(points, t):
    """Value of a step function (list of (time, value) sorted) at time t."""
    cur = None
    for tt, v in points:
        if tt <= t:
            cur = v
        else:
            break
    return cur


def integrate(trace, T, bks=None):
    """Return dict of metrics over horizon T."""
    solved = trace.status in ("optimal", "infeasible")
    if solved and (trace.end_time is None or trace.end_time <= 0.0):
        out = dict(pdi=None, t_first=None, t_gap10=None, t_gap1=None,
                   final_gap=0.0, primal=trace.final_primal, dual=trace.final_dual)
        if bks is not None:
            out["primal_int"] = None
            out["dual_int"] = None
        return out
    times = {0.0, T} | {t for t, _ in trace.primal if t <= T}
    times |= {t for t, _ in trace.dual if t <= T}
    if solved and trace.end_time is not None and 0.0 <= trace.end_time <= T:
        times.add(trace.end_time)
    times = sorted(times)
    end = trace.end_time if trace.end_time is not None else (times[-2] if len(times) > 1 else 0.0)
    pdi = pint = dint = 0.0
    t_first = t10 = t1 = None
    for a, b in zip(times[:-1], times[1:]):
        p, d = step_at(trace.primal, a), step_at(trace.dual, a)
        g = gap_value(p, d)
        if solved and a >= end:
            g = 0.0
        pdi += g * (b - a)
        if bks is not None:
            pint += one_sided(p, bks) * (b - a)
            dint += one_sided(d, bks) * (b - a)
        if t_first is None and p is not None:
            t_first = a
        if t10 is None and g <= 0.10:
            t10 = a
        if t1 is None and g <= 0.01:
            t1 = a
    last_p, last_d = step_at(trace.primal, T), step_at(trace.dual, T)
    final_gap = 0.0 if solved else gap_value(last_p, last_d)
    out = dict(pdi=pdi / T, t_first=t_first, t_gap10=t10, t_gap1=t1, final_gap=final_gap,
               primal=last_p, dual=last_d)
    if bks is not None:
        out["primal_int"] = pint / T
        out["dual_int"] = dint / T
    if solved and trace.status == "optimal" and t1 is None:
        out["t_gap1"] = end
    return out


def gap_curve(trace, T, grid):
    """Gap at each grid time (for plotting)."""
    end = trace.end_time
    solved = trace.status in ("optimal", "infeasible")
    out = []
    for t in grid:
        if solved and end is not None and t >= end:
            out.append(0.0)
        else:
            out.append(gap_value(step_at(trace.primal, t), step_at(trace.dual, t)))
    return out


# ---------------------------------------------------------------- data


def read_tsv(path):
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return {r["instance"]: r for r in csv.DictReader(f, delimiter="\t")}


def load_meta(data_dir):
    return read_tsv(os.path.join(data_dir, "bks.tsv")), read_tsv(os.path.join(data_dir, "limits.tsv"))


def bks_value(bks, name):
    r = bks.get(name)
    if not r or r.get("status") == "infeasible":
        return None, None
    try:
        return float(r["objective"]), r.get("sense", "min")
    except (KeyError, ValueError):
        return None, None


# ---------------------------------------------------------------- run


def find_instances(inst_dir, patterns):
    names = sorted(f[:-len(".mps.gz")] for f in os.listdir(inst_dir) if f.endswith(".mps.gz"))
    if not patterns:
        return names
    return [n for n in names if any(fnmatch.fnmatch(n, p) for p in patterns)]


def cmd_run(a):
    bks, limits = load_meta(a.data_dir)
    inst_dir = a.inst_dir or os.path.join(a.data_dir, "mips")
    names = find_instances(inst_dir, a.patterns)
    os.makedirs(a.out, exist_ok=True)
    extra = a.args.split() if a.args else []

    def one(name):
        sec = a.sec or float(limits.get(name, {}).get("time_limit_sec") or 120)
        log = os.path.join(a.out, name + ".log")
        cmd = [a.cbc, os.path.join(inst_dir, name + ".mps.gz"), "-threads", str(a.threads),
               "-sec", str(sec),
               "-progressCsv", os.path.join(a.out, name + ".progress.csv")] + ([] if a.nodes is None else ["-maxNodes", str(a.nodes)]) + extra + ["-solve"]
        with open(log, "w") as f:
            f.write("# cmd: %s\n# limit: %g\n" % (" ".join(cmd), sec))
            f.flush()
            try:
                subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, timeout=sec * 1.5 + 30)
            except subprocess.TimeoutExpired:
                f.write("\n# HARD KILL\n")
        return name

    with cf.ThreadPoolExecutor(a.jobs) as ex:
        for i, n in enumerate(ex.map(one, names), 1):
            print("[%d/%d] %s" % (i, len(names), n), file=sys.stderr)
    print("logs in", a.out)


# ---------------------------------------------------------------- analyze

COLS = ["instance", "status", "events", "T", "pdi", "primal_int", "dual_int", "t_first", "t_gap10",
        "t_gap1", "final_gap", "primal", "dual", "bks"]


def log_limit(path):
    with open(path, errors="replace") as f:
        for _ in range(3):
            m = re.match(r"# limit:\s*(\S+)", f.readline())
            if m:
                return float(m.group(1))
    return None


def log_threads(path):
    """Value of -threads in a cbc_pdi command header, or 0 if absent."""
    with open(path, errors="replace") as f:
        for line in f:
            if line.startswith("# command:"):
                try:
                    command = json.loads(line.split(":", 1)[1].strip())
                except (ValueError, TypeError):
                    command = []
                if "-threads" in command:
                    try:
                        return int(command[command.index("-threads") + 1])
                    except (ValueError, IndexError):
                        pass
            if line.startswith("# cmd:"):
                match = re.search(r"\s-threads\s+(\d+)", line)
                if match:
                    return int(match.group(1))
            if line.startswith("CBC ") or line.startswith("CBC version"):
                break
    return 0


def run_events(path, pcsv, mode):
    """Choose one event source; parallel log fallback has a distinct label."""
    ev = parse_log(path)
    have = os.path.exists(pcsv) and os.path.getsize(pcsv) > 0
    threads = log_threads(path)
    if mode == "csv" or (mode == "auto" and have and threads <= 1):
        if have:
            ev = apply_progress_csv(ev, pcsv)
        else:
            ev["event_source"] = "timestampedLog (no progress CSV)"
    elif mode == "auto":
        ev["event_source"] = (
            "timestampedLog (parallel CSV incomplete)"
            if have else "timestampedLog (no progress CSV)"
        )
    else:
        ev["event_source"] = "timestampedLog"
    return ev


def analyze_dir(d, a):
    bks, limits = load_meta(a.data_dir)
    rows = []
    for fn in sorted(os.listdir(d)):
        if not fn.endswith(".log"):
            continue
        name = fn[:-4]
        path = os.path.join(d, fn)
        ev = run_events(path, os.path.join(d, name + ".progress.csv"), a.events)
        ref, sense = bks_value(bks, name)
        maximize = (sense == "max") if sense else False
        tr = build_trace(ev, maximize)
        T = a.sec or log_limit(path) or float((limits.get(name) or {}).get("time_limit_sec") or 0) or (tr.end_time or 1.0)
        T = max(T, 1e-9)
        ref_min = None if ref is None else (-ref if maximize else ref)
        if ref_min is not None and tr.primal and tr.primal[-1][1] < ref_min - 1e-4 * max(1.0, abs(ref_min)):
            print("WARNING %s: incumbent %g beats best-known %g (parse error or stale bks)" % (
                name, tr.primal[-1][1], ref_min), file=sys.stderr)
        m = integrate(tr, T, ref_min)
        sg = -1.0 if maximize else 1.0
        row = dict(instance=name, status=tr.status,
                   events=ev.get("event_source", "unknown"),
                   T=T, bks="" if ref is None else ref, **m)
        row["primal"] = "" if m["primal"] is None else sg * m["primal"]
        row["dual"] = "" if m["dual"] is None else sg * m["dual"]
        rows.append(row)
    out = os.path.join(d, "pdi.tsv")
    with open(out, "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(COLS)
        for r in rows:
            w.writerow([fmt(r.get(c)) for c in COLS])
    return rows, out


def fmt(v):
    if v is None or v == "":
        return ""
    if isinstance(v, float):
        return "%.6g" % v
    return str(v)


def cmd_analyze(a):
    for d in a.dirs:
        rows, out = analyze_dir(d, a)
        rows = [r for r in rows if r["pdi"] is not None]
        n = len(rows)
        if not n:
            print(d, ": no logs")
            continue
        print("%s: %d instances, mean PDI %.4f, no-solution %d, optimal %d -> %s" % (
            d, n, sum(r["pdi"] for r in rows) / n, sum(1 for r in rows if r["primal"] == ""),
            sum(1 for r in rows if r["status"] == "optimal"), out))


# ---------------------------------------------------------------- compare


def load_pdi(d):
    p = os.path.join(d, "pdi.tsv")
    if not os.path.exists(p):
        sys.exit("%s missing: run 'analyze %s' first" % (p, d))
    return read_tsv(p)


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def parse_exps(items):
    out = []
    for it in items:
        if "=" in it:
            n, p = it.split("=", 1)
        else:
            n, p = os.path.basename(os.path.normpath(it)), it
        out.append((n, p))
    return out


def cmd_compare(a):
    exps = [(n, load_pdi(p)) for n, p in parse_exps(a.exps)]
    if len(exps) < 2:
        sys.exit("need at least two experiments")
    common = sorted(i for i in set.intersection(*[set(r) for _, r in exps])
                    if all(fnum(r[i]["pdi"]) is not None for _, r in exps))
    if not common:
        sys.exit("no common instances")
    print("%d common instances\n" % len(common))
    mism = [i for i in common if len({fnum(r[i]["T"]) for _, r in exps}) > 1]
    if mism:
        print("WARNING: %d instances have different time horizons T across experiments (e.g. %s); "
              "PDI is only comparable on one horizon: re-run 'analyze --sec <T>' on every directory.\n"
              % (len(mism), mism[0]), file=sys.stderr)
    hdr = "%-14s %8s %8s %8s %8s %8s %8s %8s %6s" % (
        "experiment", "PDI", "primInt", "dualInt", "t_first", "t_gap1", "finalGap", "noSol", "opt")
    print(hdr)
    print("-" * len(hdr))
    stats = {}
    for n, r in exps:
        def mean(col):
            v = [fnum(r[i][col]) for i in common if fnum(r[i][col]) is not None]
            return sum(v) / len(v) if v else float("nan")
        # shifted geometric mean of t_first (shift 1s); unsolved-to-incumbent -> T
        tf = [fnum(r[i]["t_first"]) if fnum(r[i]["t_first"]) is not None else fnum(r[i]["T"]) for i in common]
        tg = [fnum(r[i]["t_gap1"]) if fnum(r[i]["t_gap1"]) is not None else fnum(r[i]["T"]) for i in common]
        sgm = lambda xs: math.exp(sum(math.log(x + 1) for x in xs) / len(xs)) - 1
        nosol = sum(1 for i in common if r[i]["primal"] == "")
        opt = sum(1 for i in common if r[i]["status"] == "optimal")
        stats[n] = mean("pdi")
        print("%-14s %8.4f %8.4f %8.4f %8.2f %8.2f %8.4f %8d %6d" % (
            n, mean("pdi"), mean("primal_int"), mean("dual_int"), sgm(tf), sgm(tg),
            mean("final_gap"), nosol, opt))
    print("\n(t_first / t_gap1: shifted geometric mean, shift 1s; a missing event counts as T)")
    rank = sorted(stats, key=stats.get)
    print("Ranking by PDI (lower is better): " + " < ".join("%s (%.4f)" % (n, stats[n]) for n in rank))

    base_n, base = exps[0]
    for n, r in exps[1:]:
        diffs = [(fnum(r[i]["pdi"]) - fnum(base[i]["pdi"]), i) for i in common]
        better = [x for x in diffs if x[0] < -a.tol]
        worse = [x for x in diffs if x[0] > a.tol]
        print("\n%s vs %s: better %d, worse %d, equal %d (tol %g)" % (
            n, base_n, len(better), len(worse), len(diffs) - len(better) - len(worse), a.tol))
        for title, lst, rev in (("Biggest PDI improvements", better, False), ("Biggest PDI regressions", worse, True)):
            lst = sorted(lst, reverse=rev)[: a.top]
            if lst:
                print(" " + title + ":")
                for dlt, i in lst:
                    print("   %-52s %7.4f -> %7.4f  (%+.4f)" % (
                        i, fnum(base[i]["pdi"]), fnum(r[i]["pdi"]), dlt))


# ---------------------------------------------------------------- plot


def cmd_credit(a):
    """Which method produced the improvements (from *.progress.csv files)."""
    prim = {}   # source -> [count, sum relative improvement]
    dual = {}
    first = {}
    n = 0
    for d in a.dirs:
        for fn in sorted(os.listdir(d)):
            if not fn.endswith(".progress.csv"):
                continue
            n += 1
            lp, ld = None, None
            with open(os.path.join(d, fn)) as f:
                for r in csv.DictReader(f):
                    p_, d_, src = float(r["primal"]), float(r["dual"]), r["source"]
                    if src in ("start", "final"):
                        continue
                    if p_ < 1e49 and (lp is None or p_ < lp - 1e-9 * max(1, abs(p_))):
                        if lp is None:
                            first[src] = first.get(src, 0) + 1
                        else:
                            e = prim.setdefault(src, [0, 0.0])
                            e[0] += 1
                            e[1] += (lp - p_) / max(1.0, abs(lp))
                        lp = p_
                    if d_ > -1e49 and (ld is None or d_ > ld + 1e-9 * max(1, abs(d_))):
                        if ld is not None:
                            e = dual.setdefault(src, [0, 0.0])
                            e[0] += 1
                            e[1] += (d_ - ld) / max(1.0, abs(ld))
                        ld = d_
    print("%d runs" % n)
    for title, tab in (("first incumbent", {k: [v, 0.0] for k, v in first.items()}),
                       ("primal improvements", prim), ("dual improvements", dual)):
        print("\n%-22s %8s %12s" % (title, "count", "sum rel.chg"))
        for k, (c, v) in sorted(tab.items(), key=lambda kv: -kv[1][0]):
            print("  %-20s %8d %12.4f" % (k, c, v))


def cmd_plot(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    bks, limits = load_meta(a.data_dir)
    exps = parse_exps(a.exps)
    grid_n = 200
    fig, ax = plt.subplots(figsize=(8, 5))
    for n, d in exps:
        curves = []
        for fn in sorted(os.listdir(d)):
            if not fn.endswith(".log"):
                continue
            name = fn[:-4]
            path = os.path.join(d, fn)
            ref, sense = bks_value(bks, name)
            tr = build_trace(parse_log(path), sense == "max")
            T = a.sec or log_limit(path) or float((limits.get(name) or {}).get("time_limit_sec") or 0) or (tr.end_time or 1.0)
            # normalised time axis: 0..1 of the limit; x ticks relabelled below
            grid = [T * k / (grid_n - 1) for k in range(grid_n)]
            curves.append(gap_curve(tr, T, grid))
        if not curves:
            continue
        mean = [sum(c[k] for c in curves) / len(curves) for k in range(grid_n)]
        ax.plot([k / (grid_n - 1) for k in range(grid_n)], mean,
                label="%s (PDI %.3f)" % (n, sum(mean) / grid_n))
    ax.set_xlabel("fraction of time limit")
    ax.set_ylabel("mean primal-dual gap (2 = no incumbent)")
    ax.set_ylim(0, 2.05)
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(a.out, dpi=130)
    print("wrote", a.out)


def gap_steps(trace, T):
    """Breakpoints and the gap held from each one (for plotting/inspection)."""
    times = sorted({0.0, T} | {t for t, _ in trace.primal if t <= T} | {t for t, _ in trace.dual if t <= T})
    end = trace.end_time
    solved = trace.status in ("optimal", "infeasible")
    gaps = []
    for a in times[:-1]:
        g = gap_value(step_at(trace.primal, a), step_at(trace.dual, a))
        gaps.append(0.0 if solved and end is not None and a >= end else g)
    return times, gaps


def cmd_explain(a):
    """One figure that shows what the PDI is: bounds on top, gap(t) below, and
    the area under gap(t) (= PDI * T) shaded. Several runs are overlaid."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    bks, limits = load_meta(a.data_dir)
    runs = []
    for n, path in parse_exps(a.runs):
        name = os.path.basename(path)[:-4]
        ref, sense = bks_value(bks, name.split(".")[0] if not os.path.exists(path) else name)
        maximize = (sense == "max") if sense else False
        tr = build_trace(run_events(path, path[:-4] + ".progress.csv", a.events), maximize)
        T = a.sec or log_limit(path) or (tr.end_time or 1.0)
        runs.append((n, tr, T, maximize, ref))
    T = max(r[2] for r in runs) if not a.sec else a.sec
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True, gridspec_kw=dict(height_ratios=[1.1, 1]))
    for k, (n, tr, _, maximize, ref) in enumerate(runs):
        c = colors[k % len(colors)]
        sg = -1.0 if maximize else 1.0
        for pts, ls, lab in ((tr.primal, "-", "primal (incumbent)"), (tr.dual, "--", "dual (bound)")):
            if not pts:
                continue
            xs = [t for t, _ in pts if t <= T]
            ys = [sg * v for t, v in pts if t <= T]
            if not xs:
                continue
            xs, ys = [0.0] + xs + [T], [ys[0]] + ys + [ys[-1]]
            ax1.step(xs, ys, where="post", color=c, ls=ls, lw=1.6, label="%s: %s" % (n, lab))
        if tr.end_time is not None and tr.end_time < T:
            ax1.axvline(tr.end_time, color=c, ls=":", lw=1)
            ax2.axvline(tr.end_time, color=c, ls=":", lw=1)
            ax2.annotate("%s ends (%s)\ngap is held" % (n, tr.status), (tr.end_time, 1.75), color=c,
                         fontsize=8, ha="right")
        times, gaps = gap_steps(tr, T)
        m = integrate(tr, T)
        ax2.fill_between(times, gaps + [gaps[-1]], step="post", color=c, alpha=0.25)
        ax2.step(times, gaps + [gaps[-1]], where="post", color=c, lw=1.6,
                 label="%s: PDI = %.4f (area / T)" % (n, m["pdi"]))
    ref = next((r[4] for r in runs if r[4] is not None), None)
    if ref is not None:
        ax1.axhline(ref, color="gray", ls=":", lw=1, label="best known")
    ax1.set_ylabel("objective value")
    ax1.legend(fontsize=8, ncol=2)
    ax1.grid(alpha=0.3)
    ax2.set_ylabel("gap(t)   (2 = no incumbent)")
    ax2.set_xlabel("time (s)   - horizon T = %g s" % T)
    ax2.set_ylim(0, 2.05)
    ax2.grid(alpha=0.3)
    ax2.legend(loc="center right", fontsize=9)
    if a.logx:
        ax2.set_xscale("symlog", linthresh=1.0)
    ax2.set_xlim(0, T)
    fig.suptitle("Primal-dual integral: PDI = shaded area under gap(t) divided by T  (smaller is better)", fontsize=10)
    fig.tight_layout()
    fig.savefig(a.out, dpi=130)
    print("wrote", a.out)


# ---------------------------------------------------------------- main


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", default=DATA, help="mip-sanity-data directory (bks.tsv, limits.tsv, mips/)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="run cbc on instances")
    r.add_argument("--out", required=True)
    r.add_argument("--cbc", default=os.path.expanduser("~/prog/cbc/bin/cbc"))
    r.add_argument("--sec", type=float, help="time limit (default: limits.tsv, else 120)")
    r.add_argument("--nodes", type=int, help="optional -maxNodes")
    r.add_argument("--threads", type=int, default=1)
    r.add_argument("--jobs", type=int, default=os.cpu_count())
    r.add_argument("--args", help="extra cbc arguments, one quoted string")
    r.add_argument("--inst-dir", help="directory of <name>.mps.gz (default: <data-dir>/mips)")
    r.add_argument("patterns", nargs="*", help="instance name globs")
    r.set_defaults(fn=cmd_run)

    z = sub.add_parser("analyze", help="logs -> pdi.tsv")
    z.add_argument("dirs", nargs="+")
    z.add_argument("--sec", type=float, help="override the time horizon T (use the same value for every "
                   "experiment you compare; shorter runs hold their last gap, longer ones are truncated)")
    z.add_argument("--events", choices=("auto", "csv", "log"), default="auto",
                   help="event source: auto = *.progress.csv unless the run used -threads N>0 (then the log)")
    z.set_defaults(fn=cmd_analyze)

    c = sub.add_parser("compare", help="rank experiments (pdi.tsv files)")
    c.add_argument("exps", nargs="+", help="[name=]dir ...; the first is the baseline")
    c.add_argument("--tol", type=float, default=0.005)
    c.add_argument("--top", type=int, default=10)
    c.set_defaults(fn=cmd_compare)

    cr = sub.add_parser("credit", help="improvements per source (needs *.progress.csv)")
    cr.add_argument("dirs", nargs="+")
    cr.set_defaults(fn=cmd_credit)

    p = sub.add_parser("plot", help="mean gap curves")
    p.add_argument("exps", nargs="+")
    p.add_argument("--out", default="gap-curves.png")
    p.add_argument("--sec", type=float)
    p.set_defaults(fn=cmd_plot)

    ex = sub.add_parser("explain", help="figure of bounds, gap(t) and the PDI area for one or more runs")
    ex.add_argument("runs", nargs="+", help="[name=]path/to/<instance>.log ...")
    ex.add_argument("--out", default="pdi-explained.png")
    ex.add_argument("--sec", type=float, help="common time horizon T (default: longest limit)")
    ex.add_argument("--events", choices=("auto", "csv", "log"), default="auto")
    ex.add_argument("--logx", action="store_true", help="symlog time axis")
    ex.set_defaults(fn=cmd_explain)

    a = ap.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
