#!/bin/bash
# Generate all MIPster documentation from the installed binary.
#
# Usage:
#   ./doc/generate_docs.sh [path/to/mipster]
#
# Generates:
#   doc/mipster.1              — man page
#   doc/mipster-parameters.tex — LaTeX reference
#   doc/mipster-parameters.pdf — PDF reference (if pdflatex available)
#   doc/mipster-parameters.md  — Markdown reference

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
MIPSTER="${1:-$HOME/prog/cbc/bin/mipster}"

if [ ! -x "$MIPSTER" ]; then
    echo "Error: mipster binary not found at $MIPSTER" >&2
    echo "Usage: $0 [path/to/mipster]" >&2
    exit 1
fi

echo "Using: $MIPSTER"

echo "Generating man page..."
python3 "$SCRIPT_DIR/generate_manpage.py" "$MIPSTER" > "$SCRIPT_DIR/mipster.1"

echo "Generating LaTeX reference..."
python3 "$SCRIPT_DIR/generate_latex_reference.py" "$MIPSTER" > "$SCRIPT_DIR/mipster-parameters.tex"

echo "Generating Markdown reference..."
python3 "$SCRIPT_DIR/generate_markdown_reference.py" "$MIPSTER" > "$SCRIPT_DIR/mipster-parameters.md"

if command -v pdflatex >/dev/null 2>&1; then
    echo "Compiling PDF..."
    TMPDIR=$(mktemp -d)
    PDFLOG="$TMPDIR/pdflatex.log"
    if pdflatex -interaction=nonstopmode -output-directory="$TMPDIR" "$SCRIPT_DIR/mipster-parameters.tex" >"$PDFLOG" 2>&1 && \
       pdflatex -interaction=nonstopmode -output-directory="$TMPDIR" "$SCRIPT_DIR/mipster-parameters.tex" >>"$PDFLOG" 2>&1; then
        mv "$TMPDIR/mipster-parameters.pdf" "$SCRIPT_DIR/mipster-parameters.pdf"
        echo "Done: $SCRIPT_DIR/mipster-parameters.pdf"
    else
        echo "Warning: pdflatex failed (non-fatal) — skipping PDF. Last 40 lines of LaTeX log:"
        tail -40 "$PDFLOG" || true
    fi
    rm -rf "$TMPDIR"
else
    echo "pdflatex not found — skipping PDF generation."
    echo "Compile manually: pdflatex doc/mipster-parameters.tex"
fi

echo ""
echo "Generated:"
ls -lh "$SCRIPT_DIR/mipster.1" "$SCRIPT_DIR/mipster-parameters.tex" "$SCRIPT_DIR/mipster-parameters.md" "$SCRIPT_DIR/mipster-parameters.pdf" 2>/dev/null
