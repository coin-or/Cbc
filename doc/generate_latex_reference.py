#!/usr/bin/env python3
"""Generate a LaTeX parameter reference from 'mipster -dumpParameters' JSON.

Usage:
    python3 doc/generate_latex_reference.py [path/to/mipster] > doc/mipster-parameters.tex
"""

import subprocess, sys, os, json, re
from datetime import datetime


def load_params(cbc):
    """Run cbc -dumpParameters and parse the JSON array."""
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


def tex_escape(s):
    s = s.replace('\\', r'\textbackslash{}')
    for old, new in [
        ('&', r'\&'), ('%', r'\%'), ('$', r'\$'), ('#', r'\#'),
        ('_', r'\_'), ('{', r'\{'), ('}', r'\}'),
        ('~', r'\textasciitilde{}'), ('^', r'\textasciicircum{}'),
        ('\u2014', '---'), ('\u2013', '--'), ('\u2018', '`'), ('\u2019', "'"),
        ('\u201c', '``'), ('\u201d', "''"), ('\u2208', r'$\in$'),
        ('\u03ba', r'$\kappa$'), ('\u2264', r'$\le$'), ('\u2265', r'$\ge$'),
    ]:
        s = s.replace(old, new)
    s = s.encode('ascii', 'replace').decode('ascii')
    return s


TYPE_LABELS = {
    "integer": "Integer", "double": "Double", "keyword": "Keyword",
    "string": "String", "action": "Action", "file": "File",
    "directory": "Directory",
}


def emit(version, topics, out):
    date = datetime.now().strftime("%B %Y")
    order = [
        "Stopping", "Cuts", "Heuristics", "Preprocessing", "Branching",
        "Tolerances", "Conflict Graph", "Strategy", "Solving",
        "Simplex", "Barrier", "Scaling", "Output", "I/O", "Parallelism",
    ]

    out.write(r"""\documentclass[11pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage[margin=2.5cm]{geometry}
\usepackage{xcolor}
\usepackage{hyperref}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{enumitem}
\usepackage{parskip}

\hypersetup{colorlinks=true, linkcolor=blue!60!black, urlcolor=blue!60!black}

\title{MIPster Parameter Reference}
\author{MIPster Development Team}
""")
    out.write(f"\\date{{{date} --- MIPster {tex_escape(version)}}}\n")
    out.write(r"""
\begin{document}
\maketitle
\tableofcontents
\newpage

\section*{Introduction}

MIPster is an open-source mixed-integer programming (MIP) solver.
This document provides a complete reference for all command-line
parameters, organized by functional topic.

Parameters are specified on the command line \emph{before} \texttt{-solve}:
\begin{verbatim}
  mipster model.mps -sec 300 -cuts ifmove -solve
\end{verbatim}

Both single-dash (\texttt{-sec}) and double-dash (\texttt{--sec}) styles
are accepted.  In interactive mode, omit the leading dash.

""")

    for topic in order:
        if topic not in topics:
            continue
        emit_topic(topic, topics[topic], out)
    for topic in topics:
        if topic in order:
            continue
        emit_topic(topic, topics[topic], out)

    out.write("\\end{document}\n")


def emit_topic(topic, params, out):
    out.write(f"\\section{{{tex_escape(topic)}}}\n\n")

    if topic == "Heuristics":
        emit_heuristics_topic(params, out)
        return

    for p in params:
        emit_param(p, out, sub="subsection")


# Heuristic parameters classified by category.
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
        out.write("\\subsection{Constructive Heuristics}\n\n")
        out.write("These heuristics do \\textbf{not} require an existing feasible solution. "
                  "They attempt to construct a feasible solution from scratch.\n\n")
        for p in constructive:
            emit_param(p, out, sub="subsubsection")

    if improvement:
        out.write("\\subsection{Improvement Heuristics}\n\n")
        out.write("These heuristics require \\textbf{at least one} existing feasible solution. "
                  "They attempt to improve upon the incumbent.\n\n")
        for p in improvement:
            emit_param(p, out, sub="subsubsection")

    if improvement2:
        out.write("\\subsection{Improvement Heuristics (2+ solutions)}\n\n")
        out.write("These heuristics require \\textbf{at least two} existing feasible solutions. "
                  "They combine or crossover multiple solutions.\n\n")
        for p in improvement2:
            emit_param(p, out, sub="subsubsection")

    if general:
        out.write("\\subsection{General Heuristic Settings}\n\n")
        for p in general:
            emit_param(p, out, sub="subsubsection")


def humanize_bound(v):
    """Replace sentinel values with human-readable labels."""
    s = str(v)
    if s == "2147483647":
        return "INT\\_MAX"
    if s == "-2147483647":
        return "$-$INT\\_MAX"
    if s in ("1e+60", "1e+20"):
        return "$\\infty$"
    if s in ("-1e+60", "-1e+20"):
        return "$-\\infty$"
    return s


def emit_param(p, out, sub="subsection"):
    name = p["name"]
    ptype = p.get("type", "")
    short_help = p.get("shortHelp", "")
    long_help = p.get("longHelp", "")
    keywords = p.get("keywords", [])

    out.write(f"\\{sub}*{{\\texttt{{-{tex_escape(name)}}}}}\n")
    out.write(f"\\addcontentsline{{toc}}{{{sub}}}"
              f"{{\\texttt{{-{tex_escape(name)}}}}}\n\n")

    # Short description
    out.write(f"{tex_escape(short_help)}\n\n")

    # Long description
    if long_help and long_help != short_help:
        out.write(f"{tex_escape(long_help)}\n\n")

    # Keywords
    if keywords:
        out.write("\\textbf{Values:} ")
        out.write(", ".join(f"\\texttt{{{tex_escape(k)}}}" for k in keywords))
        dv = p.get("defaultValue", "")
        if dv:
            out.write(f" (default: \\texttt{{{tex_escape(str(dv))}}})")
        out.write("\n\n")

    # Numeric range
    if ptype == "integer":
        lo = p.get("lowerInt", "")
        hi = p.get("upperInt", "")
        dv = p.get("defaultValue", "")
        out.write(f"\\textbf{{Range:}} {humanize_bound(lo)} to {humanize_bound(hi)}")
        if dv != "":
            out.write(f" (default: {dv})")
        out.write("\n\n")
    elif ptype == "double":
        lo = p.get("lowerDbl", "")
        hi = p.get("upperDbl", "")
        dv = p.get("defaultValue", "")
        out.write(f"\\textbf{{Range:}} {humanize_bound(lo)} to {humanize_bound(hi)}")
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
