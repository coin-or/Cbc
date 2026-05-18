#!/usr/bin/env python3
"""Generate a Markdown parameter reference from 'mipster -dumpParameters' JSON.

Usage:
    python3 doc/generate_markdown_reference.py [path/to/mipster] > doc/mipster-parameters.md
"""

import subprocess, sys, os, json, re
from datetime import datetime


def load_params(cbc):
    out = subprocess.run([cbc, "-dumpParameters"],
                         capture_output=True, text=True, timeout=10).stdout
    start = out.index("["); end = out.rindex("]") + 1
    return json.loads(out[start:end])


def get_version(cbc):
    out = subprocess.run([cbc, "--version"],
                         capture_output=True, text=True, timeout=10).stdout
    m = re.search(r"(?:MIPster|CBC)\s+(\S+)", out, re.IGNORECASE)
    return m.group(1) if m else "devel"


def group_by_topic(params):
    topics = {}
    for p in params:
        t = p.get("topic", "")
        if not t or p.get("displayPriority", 0) == 0:
            continue
        topics.setdefault(t, []).append(p)
    return topics


TYPE_LABELS = {
    "integer": "Integer", "double": "Double", "keyword": "Keyword",
    "string": "String", "action": "Action", "file": "File",
    "directory": "Directory",
}

ORDER = [
    "Stopping", "Cuts", "Heuristics", "Preprocessing", "Branching",
    "Tolerances", "Conflict Graph", "Strategy", "Solving",
    "Simplex", "Barrier", "Scaling", "Output", "I/O", "Parallelism",
]


def emit(version, topics, out):
    date = datetime.now().strftime("%B %Y")
    out.write(f"# MIPster Parameter Reference\n\n")
    out.write(f"*MIPster {version} — {date}*\n\n")
    out.write("Parameters are specified on the command line **before** `-solve`:\n")
    out.write("```\nmipster model.mps -sec 300 -cuts ifmove -solve\n```\n\n")
    out.write("Both single-dash (`-sec`) and double-dash (`--sec`) styles are accepted.\n\n")

    # Table of contents
    out.write("## Contents\n\n")
    for topic in ORDER:
        if topic in topics:
            anchor = topic.lower().replace(" ", "-")
            out.write(f"- [{topic}](#{anchor}) ({len(topics[topic])} parameters)\n")
    for topic in topics:
        if topic not in ORDER:
            anchor = topic.lower().replace(" ", "-")
            out.write(f"- [{topic}](#{anchor}) ({len(topics[topic])} parameters)\n")
    out.write("\n---\n\n")

    for topic in ORDER:
        if topic in topics:
            emit_topic(topic, topics[topic], out)
    for topic in topics:
        if topic not in ORDER:
            emit_topic(topic, topics[topic], out)


def emit_topic(topic, params, out):
    out.write(f"## {topic}\n\n")

    if topic == "Heuristics":
        emit_heuristics_topic(params, out)
        return

    for p in params:
        emit_param(p, out, level="###")


# Heuristic parameters classified by category.
# Constructive: do not need an existing feasible solution.
# Improvement: need at least one feasible solution.
# Improvement (2+ solutions): need at least two feasible solutions.
CONSTRUCTIVE_HEURISTICS = {
    "feasibilityPump", "DivingCoefficient", "DivingFractional",
    "DivingGuided", "DivingLineSearch", "DivingPseudocost",
    "DivingVectorLength", "DivingSome", "roundingHeuristic",
    "greedyHeuristic", "naiveHeuristics", "pivotAndFix",
    "randomizedRounding", "Rens",
}

IMPROVEMENT_HEURISTICS = {
    "Rins", "VndVariableNeighborhoodSearch", "Dins",
    "proximitySearch", "dwHeuristic", "localTreeSearch",
}

IMPROVEMENT_2_HEURISTICS = {
    "combineSolutions", "combine2Solutions",
}


def emit_heuristics_topic(params, out):
    """Emit heuristic parameters grouped by category."""
    constructive = []
    improvement = []
    improvement2 = []
    general = []

    for p in params:
        name = p["name"]
        if name in CONSTRUCTIVE_HEURISTICS:
            constructive.append(p)
        elif name in IMPROVEMENT_HEURISTICS:
            improvement.append(p)
        elif name in IMPROVEMENT_2_HEURISTICS:
            improvement2.append(p)
        else:
            general.append(p)

    if constructive:
        out.write("### Constructive Heuristics\n\n")
        out.write("These heuristics do **not** require an existing feasible solution. "
                  "They attempt to construct a feasible solution from scratch.\n\n")
        for p in constructive:
            emit_param(p, out, level="####")

    if improvement:
        out.write("### Improvement Heuristics\n\n")
        out.write("These heuristics require **at least one** existing feasible solution. "
                  "They attempt to improve upon the incumbent.\n\n")
        for p in improvement:
            emit_param(p, out, level="####")

    if improvement2:
        out.write("### Improvement Heuristics (2+ solutions)\n\n")
        out.write("These heuristics require **at least two** existing feasible solutions. "
                  "They combine or crossover multiple solutions.\n\n")
        for p in improvement2:
            emit_param(p, out, level="####")

    if general:
        out.write("### General Heuristic Settings\n\n")
        for p in general:
            emit_param(p, out, level="####")


def humanize_bound(v):
    """Replace sentinel values with human-readable labels."""
    s = str(v)
    if s == "2147483647":
        return "INT_MAX"
    if s == "-2147483647":
        return "-INT_MAX"
    if s in ("1e+60", "1e+20"):
        return "∞"
    if s in ("-1e+60", "-1e+20"):
        return "-∞"
    return s


def emit_param(p, out, level="###"):
    name = p["name"]
    ptype = p.get("type", "")
    short_help = p.get("shortHelp", "")
    long_help = p.get("longHelp", "")
    keywords = p.get("keywords", [])

    out.write(f"{level} `-{name}`\n\n")
    out.write(f"{short_help}\n\n")

    if long_help and long_help != short_help:
        # Replace literal \n with actual newlines
        text = long_help.replace("\\n", "\n")
        out.write(f"{text}\n\n")

    if keywords:
        dv = p.get("defaultValue", "")
        kw_str = ", ".join(f"`{k}`" for k in keywords)
        out.write(f"**Values:** {kw_str}")
        if dv:
            out.write(f" (default: `{dv}`)")
        out.write("\n\n")

    if ptype == "integer":
        lo, hi = p.get("lowerInt", ""), p.get("upperInt", "")
        dv = p.get("defaultValue", "")
        out.write(f"**Range:** {humanize_bound(lo)} to {humanize_bound(hi)}")
        if dv != "":
            out.write(f" (default: {dv})")
        out.write("\n\n")
    elif ptype == "double":
        lo, hi = p.get("lowerDbl", ""), p.get("upperDbl", "")
        dv = p.get("defaultValue", "")
        out.write(f"**Range:** {humanize_bound(lo)} to {humanize_bound(hi)}")
        if dv != "":
            out.write(f" (default: {dv})")
        out.write("\n\n")


def main():
    cbc = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/prog/cbc/bin/mipster")
    params = load_params(cbc)
    version = get_version(cbc)
    topics = group_by_topic(params)
    total = sum(len(v) for v in topics.values())
    print(f"Parsed {len(topics)} topics, {total} parameters", file=sys.stderr)
    emit(version, topics, sys.stdout)


if __name__ == "__main__":
    main()
