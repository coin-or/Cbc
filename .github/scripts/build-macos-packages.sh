#!/usr/bin/env bash
# build-macos-packages.sh — Create a macOS .pkg installer from a dist tarball.
#
# Usage: build-macos-packages.sh <tarball.tar.gz> <arch>
# e.g.:  build-macos-packages.sh mipster-macos-arm64.tar.gz arm64
#        build-macos-packages.sh mipster-macos-x86_64.tar.gz x86_64
#
# Output: mipster-<version>-macos-<arch>.pkg (in cwd)
#
# Installs to /usr/local by default (standard Homebrew prefix on Intel Macs and
# the traditional CLI prefix on Apple Silicon).  Uses only Apple-supplied tools
# (pkgbuild) — no third-party dependencies required.

set -euo pipefail

TARBALL="${1:?Usage: $0 <tarball.tar.gz> <arch>}"
ARCH="${2:?Usage: $0 <tarball.tar.gz> <arch>}"

SRC_DIR="$(cd "$(dirname "$0")/../.." && pwd)"

if [ ! -f "$TARBALL" ]; then
  echo "Error: tarball not found: $TARBALL" >&2
  exit 1
fi

TARBALL_BASENAME="$(basename "$TARBALL" .tar.gz)"

echo "==> MIPster macOS .pkg builder"
echo "    Tarball: $TARBALL"
echo "    Arch:    $ARCH"

# ── Extract tarball ───────────────────────────────────────────────────────────
WORK_DIR="/tmp/mipster-pkg-work-${ARCH}"
rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR"
tar -xzf "$TARBALL" -C "$WORK_DIR"
DIST_DIR="$WORK_DIR/$TARBALL_BASENAME"

# ── Detect version ────────────────────────────────────────────────────────────
# Try launcher or any real binary (not symlinks)
MIPSTER_BIN="$(find "$DIST_DIR/bin" -type f -name 'mipster*' | head -1)"
VERSION=$("$MIPSTER_BIN" --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || echo "0.0.0")
echo "    Version: $VERSION"

# ── Stage in FHS layout ───────────────────────────────────────────────────────
STAGING="$WORK_DIR/staging"
mkdir -p \
  "$STAGING/usr/local/bin" \
  "$STAGING/usr/local/lib" \
  "$STAGING/usr/local/include" \
  "$STAGING/usr/local/share/doc/mipster" \
  "$STAGING/usr/local/share/man/man1"

# Binaries (launcher + variants)
cp -r "$DIST_DIR/bin/"* "$STAGING/usr/local/bin/"

# Shared libraries — keep variant subdirs (generic/, haswell/) if present
if ls "$DIST_DIR/lib/"*/ >/dev/null 2>&1; then
  # x86_64: lib/generic/, lib/haswell/
  mkdir -p "$STAGING/usr/local/lib/mipster"
  cp -rP "$DIST_DIR/lib/"*/ "$STAGING/usr/local/lib/mipster/"
else
  # arm64: flat lib/
  cp -rP "$DIST_DIR/lib/"* "$STAGING/usr/local/lib/" 2>/dev/null || true
fi

# Headers
[ -d "$DIST_DIR/include/mipster" ] && \
  cp -r "$DIST_DIR/include/mipster" "$STAGING/usr/local/include/"

# Docs (from tarball if present, else from source tree)
DOC_SRC="${SRC_DIR}/doc"
for f in mipster-parameters.pdf mipster-parameters.md; do
  if [ -f "$DIST_DIR/share/doc/mipster/$f" ]; then
    cp "$DIST_DIR/share/doc/mipster/$f" "$STAGING/usr/local/share/doc/mipster/"
  elif [ -f "${DOC_SRC}/$f" ]; then
    cp "${DOC_SRC}/$f" "$STAGING/usr/local/share/doc/mipster/"
  fi
done

# Man page
if [ -f "$DIST_DIR/share/man/man1/mipster.1.gz" ]; then
  cp "$DIST_DIR/share/man/man1/mipster.1.gz" "$STAGING/usr/local/share/man/man1/"
elif [ -f "${DOC_SRC}/mipster.1" ]; then
  gzip -c "${DOC_SRC}/mipster.1" > "$STAGING/usr/local/share/man/man1/mipster.1.gz"
fi

echo "    Staging tree:"
find "$STAGING" -type f | sed "s|$STAGING||" | sort | sed 's/^/      /'

# ── Build .pkg ────────────────────────────────────────────────────────────────
PKG_NAME="mipster-${VERSION}-macos-${ARCH}.pkg"
echo ""
echo "==> Building ${PKG_NAME}..."

pkgbuild \
  --root "$STAGING" \
  --identifier "com.mipster.solver" \
  --version "$VERSION" \
  --install-location "/" \
  "$PKG_NAME"

echo "==> Done: $PKG_NAME  ($(du -sh "$PKG_NAME" | cut -f1))"
