#!/bin/bash
# Generate all CBC documentation from the installed binary.
#
# Usage:
#   ./doc/generate_docs.sh [path/to/cbc]
#
# Generates:
#   doc/cbc.1              — man page
#   doc/cbc-parameters.tex — LaTeX reference
#   doc/cbc-parameters.pdf — PDF reference (if pdflatex available)

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CBC="${1:-$HOME/prog/cbc/bin/cbc}"

if [ ! -x "$CBC" ]; then
    echo "Error: cbc binary not found at $CBC" >&2
    echo "Usage: $0 [path/to/cbc]" >&2
    exit 1
fi

echo "Using: $CBC"

echo "Generating man page..."
python3 "$SCRIPT_DIR/generate_manpage.py" "$CBC" > "$SCRIPT_DIR/cbc.1"

echo "Generating LaTeX reference..."
python3 "$SCRIPT_DIR/generate_latex_reference.py" "$CBC" > "$SCRIPT_DIR/cbc-parameters.tex"

if command -v pdflatex >/dev/null 2>&1; then
    echo "Compiling PDF..."
    TMPDIR=$(mktemp -d)
    pdflatex -interaction=nonstopmode -output-directory="$TMPDIR" "$SCRIPT_DIR/cbc-parameters.tex" >/dev/null 2>&1
    pdflatex -interaction=nonstopmode -output-directory="$TMPDIR" "$SCRIPT_DIR/cbc-parameters.tex" >/dev/null 2>&1
    mv "$TMPDIR/cbc-parameters.pdf" "$SCRIPT_DIR/cbc-parameters.pdf"
    rm -rf "$TMPDIR"
    echo "Done: $SCRIPT_DIR/cbc-parameters.pdf"
else
    echo "pdflatex not found — skipping PDF generation."
    echo "Compile manually: pdflatex doc/cbc-parameters.tex"
fi

echo ""
echo "Generated:"
ls -lh "$SCRIPT_DIR/cbc.1" "$SCRIPT_DIR/cbc-parameters.tex" "$SCRIPT_DIR/cbc-parameters.pdf" 2>/dev/null
