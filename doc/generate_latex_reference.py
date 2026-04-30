#!/usr/bin/env python3
"""Generate a LaTeX parameter reference from the cbc binary.

Usage:
    python3 doc/generate_latex_reference.py [path/to/cbc] > doc/cbc-parameters.tex

Runs 'cbc --help' to get the topic-grouped parameter list, then queries
each parameter with 'param??' to get its long description, type info,
and valid keywords/ranges.
"""

import subprocess
import sys
import os
import re
from datetime import datetime


def run_cbc(cbc, args):
    """Run cbc with given args, return stdout."""
    result = subprocess.run(
        [cbc] + args,
        capture_output=True, text=True, timeout=10
    )
    return result.stdout


def get_version(cbc):
    for line in run_cbc(cbc, ["--version"]).splitlines():
        m = re.search(r"CBC\s+(\S+)", line, re.IGNORECASE)
        if m:
            return m.group(1)
    return "devel"


def parse_topics(help_text):
    """Parse --help output into {topic: [(name, short_desc)]}."""
    topics = {}
    current = None
    in_params = False
    for line in help_text.splitlines():
        if line.strip() == "Parameters by topic:":
            in_params = True
            continue
        if not in_params:
            continue
        m = re.match(r"^  (\S.*):$", line)
        if m:
            current = m.group(1)
            topics[current] = []
            continue
        m = re.match(r"^    -(\S+)\s{2,}(.+)$", line)
        if m and current:
            topics[current].append((m.group(1), m.group(2).strip()))
    return topics


def query_long_help(cbc, param_name):
    """Query 'cbc -param_name??' and parse the long help text."""
    out = run_cbc(cbc, [f"-{param_name}??"])
    lines = out.splitlines()

    # Skip header lines (CBC version, args, default strategy, match line)
    help_lines = []
    started = False
    keywords = None
    param_range = None

    for line in lines:
        # Skip cbc header
        if line.startswith("CBC ") or line.startswith("  args:") or line.startswith("  default strategy"):
            continue
        if line.startswith("Total time"):
            continue
        # "Match for ..." or "Multiple matches" header
        if line.startswith("Match for") or line.startswith("Multiple matches"):
            started = True
            continue
        if not started:
            continue
        # "Possible options for X are:" line
        m = re.match(r"Possible options for \S+ are:", line)
        if m:
            continue
        # Keyword list (indented, space-separated)
        if keywords is None and re.match(r"^  \S", line) and not help_lines:
            # This might be keywords, but could also be help text
            pass
        # Range line: "Range of values is -1e+08 to 1e+08"
        m = re.match(r".*[Rr]ange.*?(-?[\d.e+]+)\s+to\s+(-?[\d.e+]+)", line)
        if m:
            param_range = (m.group(1), m.group(2))
            continue
        # Current value line
        if re.match(r"^Current (value|option)", line):
            continue
        # Keyword options line (e.g., "  off  on  root  ifmove  forceon")
        if re.match(r"^  \w+(\s{2,}\w+)+$", line):
            keywords = line.split()
            continue
        # Help text
        if line.strip():
            help_lines.append(line)

    long_help = " ".join(l.strip() for l in help_lines).strip()
    return long_help, keywords, param_range


def tex_escape(s):
    """Escape LaTeX special characters."""
    # Replace backslash first
    s = s.replace('\\', r'\textbackslash{}')
    for old, new in [
        ('&', r'\&'), ('%', r'\%'), ('$', r'\$'), ('#', r'\#'),
        ('_', r'\_'), ('{', r'\{'), ('}', r'\}'),
        ('~', r'\textasciitilde{}'), ('^', r'\textasciicircum{}'),
    ]:
        s = s.replace(old, new)
    # Replace common UTF-8 characters that may appear in help text
    for old, new in [
        ('\u2014', '---'), ('\u2013', '--'), ('\u2018', '`'), ('\u2019', "'"),
        ('\u201c', '``'), ('\u201d', "''"), ('\u2208', r'$\in$'),
        ('\u03ba', r'$\kappa$'), ('\u2264', r'$\le$'), ('\u2265', r'$\ge$'),
    ]:
        s = s.replace(old, new)
    # Strip any remaining non-ASCII that LaTeX can't handle
    s = s.encode('ascii', 'replace').decode('ascii')
    return s


def emit_latex(version, topics, param_details, out):
    """Write the LaTeX document."""
    date = datetime.now().strftime("%B %Y")

    order = [
        "Stopping", "Cuts", "Heuristics", "Preprocessing", "Branching",
        "Tolerances", "Conflict Graph", "Strategy", "Solving",
        "Simplex", "Barrier", "Scaling",
        "Output", "I/O", "Parallelism",
    ]

    out.write(r"""\documentclass[11pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage[margin=2.5cm]{geometry}
\usepackage{hyperref}
\usepackage{xcolor}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{enumitem}
\usepackage{parskip}

\hypersetup{colorlinks=true, linkcolor=blue!60!black, urlcolor=blue!60!black}

\title{CBC Parameter Reference}
\author{COIN-OR Cbc Development Team}
""")
    out.write(f"\\date{{{date} — CBC {tex_escape(version)}}}\n")
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
        emit_topic(topic, topics[topic], param_details, out)

    for topic in topics:
        if topic in order:
            continue
        emit_topic(topic, topics[topic], param_details, out)

    out.write(r"""
\end{document}
""")


def emit_topic(topic, params, param_details, out):
    out.write(f"\\section{{{tex_escape(topic)}}}\n\n")

    for name, short_desc in params:
        long_help, keywords, param_range = param_details.get(name, ("", None, None))

        out.write(f"\\subsection*{{\\texttt{{-{tex_escape(name)}}}}}\n")
        out.write("\\addcontentsline{toc}{subsection}")
        out.write(f"{{\\texttt{{-{tex_escape(name)}}}}}\n\n")

        # Short description
        out.write(f"{tex_escape(short_desc)}\n\n")

        # Long description (if different from short)
        if long_help and long_help != short_desc:
            out.write(f"{tex_escape(long_help)}\n\n")

        # Keywords
        if keywords:
            out.write("\\textbf{Values:} ")
            out.write(", ".join(f"\\texttt{{{tex_escape(k)}}}" for k in keywords))
            out.write("\n\n")

        # Range
        if param_range:
            lo, hi = param_range
            out.write(f"\\textbf{{Range:}} {tex_escape(lo)} to {tex_escape(hi)}\n\n")


def main():
    cbc = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/prog/cbc/bin/cbc")

    print("Parsing --help output...", file=sys.stderr)
    help_text = run_cbc(cbc, ["--help"])
    version = get_version(cbc)
    topics = parse_topics(help_text)

    if not topics:
        print(f"Error: no topics found in '{cbc} --help' output", file=sys.stderr)
        sys.exit(1)

    total = sum(len(v) for v in topics.values())
    print(f"Found {len(topics)} topics, {total} parameters", file=sys.stderr)

    # Query long help for each parameter
    seen = set()
    param_details = {}
    count = 0
    for topic, params in topics.items():
        for name, short_desc in params:
            if name in seen:
                continue
            seen.add(name)
            count += 1
            if count % 20 == 0:
                print(f"  querying parameter {count}/{total}...", file=sys.stderr)
            long_help, keywords, param_range = query_long_help(cbc, name)
            param_details[name] = (long_help, keywords, param_range)

    print(f"Queried {len(param_details)} unique parameters", file=sys.stderr)

    emit_latex(version, topics, param_details, sys.stdout)


if __name__ == "__main__":
    main()
