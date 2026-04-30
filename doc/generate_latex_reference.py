#!/usr/bin/env python3
"""Generate a LaTeX parameter reference from 'cbc -dumpParameters' JSON.

Usage:
    python3 doc/generate_latex_reference.py [path/to/cbc] > doc/cbc-parameters.tex
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
    m = re.search(r"CBC\s+(\S+)", out, re.IGNORECASE)
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

\title{CBC Parameter Reference}
\author{COIN-OR Cbc Development Team}
""")
    out.write(f"\\date{{{date} --- CBC {tex_escape(version)}}}\n")
    out.write(r"""
\begin{document}
\maketitle
\tableofcontents
\newpage

\section*{Introduction}

CBC (COIN-OR Branch and Cut) is an open-source mixed-integer programming
solver.  This document provides a complete reference for all command-line
parameters, organized by functional topic.

Parameters are specified on the command line \emph{before} \texttt{-solve}:
\begin{verbatim}
  cbc model.mps -sec 300 -cuts ifmove -solve
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

    for p in params:
        name = p["name"]
        ptype = p.get("type", "")
        short_help = p.get("shortHelp", "")
        long_help = p.get("longHelp", "")
        source = p.get("source", "")
        keywords = p.get("keywords", [])

        out.write(f"\\subsection*{{\\texttt{{-{tex_escape(name)}}}}}\n")
        out.write(f"\\addcontentsline{{toc}}{{subsection}}"
                  f"{{\\texttt{{-{tex_escape(name)}}}}}\n\n")

        # Type and source badge
        tlabel = TYPE_LABELS.get(ptype, ptype)
        badge = f"\\textit{{{tlabel}}}"
        if source:
            badge += f" ({source})"
        out.write(f"{badge}\n\\medskip\n\n")

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
            out.write(f"\\textbf{{Range:}} {lo} to {hi}")
            if dv != "":
                out.write(f" (default: {dv})")
            out.write("\n\n")
        elif ptype == "double":
            lo = p.get("lowerDbl", "")
            hi = p.get("upperDbl", "")
            dv = p.get("defaultValue", "")
            out.write(f"\\textbf{{Range:}} {lo} to {hi}")
            if dv != "":
                out.write(f" (default: {dv})")
            out.write("\n\n")


def main():
    cbc = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/prog/cbc/bin/cbc")
    params = load_params(cbc)
    version = get_version(cbc)
    topics = group_by_topic(params)
    total = sum(len(v) for v in topics.values())
    print(f"Parsed {len(topics)} topics, {total} parameters", file=sys.stderr)
    emit(version, topics, sys.stdout)


if __name__ == "__main__":
    main()
