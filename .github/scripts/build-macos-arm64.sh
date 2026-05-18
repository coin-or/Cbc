#!/usr/bin/env bash
# build-macos-arm64.sh — Build MIPster on macOS ARM64 (Apple Silicon).
#
# Produces: dist/mipster-macos-arm64.tar.gz
#   bin/mipster              — single ARM64 binary (all M1/M2/M3/M4)
#   lib/libmipster.dylib     — shared library
#   include/coin-or/         — public C API headers
#
# Single variant: all Apple Silicon (M1 through M4+) shares the same
# relevant ISA. Generational differences are microarchitectural, not
# exploitable via -march for a MIP solver.
#
# Uses Apple's Accelerate framework for BLAS/LAPACK (no Homebrew deps).
# Requires Xcode Command Line Tools.

set -euo pipefail

DIST_NAME="mipster-macos-arm64"
INSTALL_DIR="$(pwd)/dist/${DIST_NAME}"
BUILD_DIR="/tmp/mipster-build-arm64"
SRC_DIR="$(cd "$(dirname "$0")/../.." && pwd)"

# macOS 11.0 Big Sur: first macOS with Apple Silicon support
export MACOSX_DEPLOYMENT_TARGET=11.0

echo "==> MIPster macOS arm64 build"
echo "    Clang: $(clang --version | head -1)"
echo "    Source: ${SRC_DIR}"

mkdir -p "${INSTALL_DIR}/bin" "${INSTALL_DIR}/include" "${INSTALL_DIR}/lib"

# ── Build ─────────────────────────────────────────────────────────────────────
echo ""
echo "==> Building (arm64, ARMv8.5-A / M1 baseline)  CXXFLAGS='-O3 -march=armv8.5-a -ffp-contract=off'"

rm -rf "${BUILD_DIR}"
mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

# Accelerate provides BLAS + LAPACK, always available on macOS
"${SRC_DIR}/configure" \
  --prefix="${BUILD_DIR}/install" \
  --enable-shared \
  --enable-static \
  --without-amd \
  '--with-lapack-lflags=-framework Accelerate' \
  CXXFLAGS="-O3 -march=armv8.5-a -ffp-contract=off" \
  LDFLAGS="-mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}" \
  2>&1 | tail -3

make -j"$(sysctl -n hw.logicalcpu)" 2>&1 | tail -3
make install 2>&1 | tail -2
echo "    Build: OK"

# ── Test ──────────────────────────────────────────────────────────────────────
cd "${BUILD_DIR}/test"
make -j"$(sysctl -n hw.logicalcpu)" CInterfaceTest 2>&1 | tail -2
./CInterfaceTest
echo "    CInterfaceTest: PASSED"

# ── Binary ────────────────────────────────────────────────────────────────────
strip "${BUILD_DIR}/src/.libs/mipster"
cp "${BUILD_DIR}/src/.libs/mipster" "${INSTALL_DIR}/bin/mipster"
echo "    bin/mipster: $(du -sh "${INSTALL_DIR}/bin/mipster" | cut -f1)"

# ── Shared library ────────────────────────────────────────────────────────────
cp -P "${BUILD_DIR}/install/lib"/libmipster*.dylib "${INSTALL_DIR}/lib/" 2>/dev/null || true

real_dylib=$(find "${INSTALL_DIR}/lib" -name 'libmipster.*.dylib' ! -type l | head -1)
if [ -n "${real_dylib}" ]; then
  install_name_tool -id "@loader_path/$(basename "${real_dylib}")" "${real_dylib}"
  strip -x "${real_dylib}"
fi

echo "    Shared lib deps:"
if [ -n "${real_dylib}" ]; then
  otool -L "${real_dylib}" | grep -v "libmipster\|/usr/lib\|/System\|^${INSTALL_DIR}" \
    || echo "      (none — only system libs)"
fi

# ── Headers ───────────────────────────────────────────────────────────────────
cp -r "${BUILD_DIR}/install/include/coin-or" "${INSTALL_DIR}/include/"
echo "    include/coin-or: copied"

# ── Smoke test ────────────────────────────────────────────────────────────────
echo ""
echo "==> Smoke test (--version):"
"${INSTALL_DIR}/bin/mipster" --version 2>&1 | head -1 || true
echo "    PASSED"

# ── Package ───────────────────────────────────────────────────────────────────
echo ""
echo "==> Packaging..."
cd "${SRC_DIR}/dist"
tar -czf "${DIST_NAME}.tar.gz" "${DIST_NAME}/"
echo "==> Done: dist/${DIST_NAME}.tar.gz  ($(du -sh "${DIST_NAME}.tar.gz" | cut -f1))"
tar -tzf "${DIST_NAME}.tar.gz"
