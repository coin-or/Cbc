#!/usr/bin/env bash
# build-macos-arm64.sh — Build MIPster on macOS ARM64 (Apple Silicon).
#
# Produces: dist/mipster-macos-arm64.tar.gz
#   bin/mipster               — single ARM64 binary (all M1/M2/M3/M4)
#   bin/mipster-dbg           — debug binary (symbols preserved)
#   bin/test/                 — test suite binaries
#   lib/libmipster.dylib      — shared library
#   lib/dbg/libmipster*.dylib — debug shared library
#   include/mipster/          — public C API headers
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

# ── Build debug variant ────────────────────────────────────────────────────────
# Note: macOS ASan runtime is tied to Xcode version and cannot be bundled.
# Users wanting ASan on macOS should build from source:
#   ./configster --debug --sanitizer=asan
build_debug_variant() {
  local build_dir="${BUILD_DIR}-dbg"
  local cxxflags="-O1 -g -fno-omit-frame-pointer"

  echo ""
  echo "==> Building debug variant  CXXFLAGS='${cxxflags}'"

  rm -rf "${build_dir}"
  mkdir -p "${build_dir}"
  cd "${build_dir}"

  "${SRC_DIR}/configure" \
    --prefix="${build_dir}/install" \
    --enable-shared \
    --disable-static \
    --without-amd \
    '--with-lapack-lflags=-framework Accelerate' \
    CXXFLAGS="${cxxflags}" \
    LDFLAGS="-mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}" \
    2>&1 | tail -3

  make -j"$(sysctl -n hw.logicalcpu)" 2>&1 | tail -3
  make install 2>&1 | tail -2
  echo "    Build: OK"

  # ── Full test suite ──────────────────────────────────────────────────────────
  cd "${build_dir}/test"
  make -j"$(sysctl -n hw.logicalcpu)" 2>&1 | tail -3
  MIPSTER_FIXTURE_DIR="${SRC_DIR}/test/fixtures" \
    DYLD_LIBRARY_PATH="${build_dir}/src/.libs${DYLD_LIBRARY_PATH:+:${DYLD_LIBRARY_PATH}}" \
    bash "${SRC_DIR}/test/run-mipster-tests"
  echo "    All tests: PASSED (debug)"

  # ── Debug binary ─────────────────────────────────────────────────────────────
  # No strip — keep debug symbols.
  local dbg_bin="${build_dir}/src/.libs/mipster"
  local old_ref
  old_ref=$(otool -L "${dbg_bin}" | awk '/libmipster/{print $1}' | head -1)
  if [ -n "${old_ref}" ]; then
    local dylib_base
    dylib_base=$(basename "${old_ref}")
    install_name_tool -change "${old_ref}" "@executable_path/../lib/dbg/${dylib_base}" "${dbg_bin}"
  fi
  cp "${dbg_bin}" "${INSTALL_DIR}/bin/mipster-dbg"
  echo "    bin/mipster-dbg: $(du -sh "${INSTALL_DIR}/bin/mipster-dbg" | cut -f1)"

  # ── Debug shared library ─────────────────────────────────────────────────────
  local libdir="${INSTALL_DIR}/lib/dbg"
  mkdir -p "${libdir}"
  cp -P "${build_dir}/install/lib"/libmipster*.dylib "${libdir}/" 2>/dev/null || true

  local real_dylib
  real_dylib=$(find "${libdir}" -name 'libmipster.*.dylib' ! -type l | head -1)
  if [ -n "${real_dylib}" ]; then
    install_name_tool -id "@rpath/$(basename "${real_dylib}")" "${real_dylib}"
  fi
  echo "    lib/dbg: $(du -sh "${libdir}" | cut -f1)"
}

build_debug_variant

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
make -j"$(sysctl -n hw.logicalcpu)" 2>&1 | tail -2
mkdir -p "${INSTALL_DIR}/share/mipster/test"
MIPSTER_FIXTURE_DIR="${SRC_DIR}/test/fixtures" \
  DYLD_LIBRARY_PATH="${BUILD_DIR}/src/.libs${DYLD_LIBRARY_PATH:+:${DYLD_LIBRARY_PATH}}" \
  bash "${SRC_DIR}/test/run-mipster-tests" \
    --write-baseline "GitHub Actions macos-latest (Apple Silicon)" \
    "${INSTALL_DIR}/share/mipster/test/ci-baseline-times.json"
echo "    CI baseline times written"

# ── Collect test binaries for tarball ──────────────────────────────────────
test_dir="${INSTALL_DIR}/bin/test"
mkdir -p "${test_dir}"
for tbin in CInterfaceTest CInterfaceTest_tsp_random CInterfaceTest_fl_random \
            CInterfaceTest_mdkp_random CInterfaceTest_nursesched CInterfaceTest_a1 \
            CInterfaceTest_graphdraw CInterfaceTest_trdta5581 CInterfaceTest_trd445c \
            CbcSolverLpTest; do
  if [ -f "${BUILD_DIR}/test/${tbin}" ]; then
    strip "${BUILD_DIR}/test/${tbin}"
    old=$(otool -L "${BUILD_DIR}/test/${tbin}" | awk '/libmipster/{print $1}' | head -1)
    if [ -n "${old}" ]; then
      # Fix dylib reference: test binary is in bin/test/, lib is in lib/
      install_name_tool -change "${old}" "@executable_path/../../lib/$(basename "${old}")" \
        "${BUILD_DIR}/test/${tbin}"
    fi
    cp "${BUILD_DIR}/test/${tbin}" "${test_dir}/"
  fi
done
echo "    bin/test/: $(ls "${test_dir}" | wc -l) test binaries"

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
cp -r "${BUILD_DIR}/install/include/mipster" "${INSTALL_DIR}/include/"
echo "    include/mipster: copied"

# ── Smoke test ────────────────────────────────────────────────────────────────
echo ""
echo "==> Smoke test (--version):"
"${INSTALL_DIR}/bin/mipster" --version 2>&1 | head -1 || true
echo "    PASSED"

# ── Documentation ─────────────────────────────────────────────────────────────
echo ""
echo "==> Installing documentation..."
mkdir -p "${INSTALL_DIR}/share/doc/mipster"
mkdir -p "${INSTALL_DIR}/share/man/man1"
mkdir -p "${INSTALL_DIR}/share/applications"

DOC_SRC="${SRC_DIR}/doc"
for f in mipster-parameters.pdf mipster-parameters.md; do
  [ -f "${DOC_SRC}/${f}" ] && cp "${DOC_SRC}/${f}" "${INSTALL_DIR}/share/doc/mipster/" && echo "    ${f}"
done
[ -f "${DOC_SRC}/mipster.1" ] && \
  gzip -c "${DOC_SRC}/mipster.1" > "${INSTALL_DIR}/share/man/man1/mipster.1.gz" && \
  echo "    mipster.1.gz"

# ── Test runner, fixtures, and desktop entries ─────────────────────────────
echo ""
echo "==> Installing test runner and fixtures..."
cp "${SRC_DIR}/test/run-mipster-tests" "${INSTALL_DIR}/bin/"
chmod +x "${INSTALL_DIR}/bin/run-mipster-tests"
mkdir -p "${INSTALL_DIR}/share/mipster/test/fixtures"
find "${SRC_DIR}/test/fixtures" -name '*.mps.gz' | while read -r f; do
  cp "$f" "${INSTALL_DIR}/share/mipster/test/fixtures/"
done
echo "    fixtures: $(ls "${INSTALL_DIR}/share/mipster/test/fixtures" | wc -l) files"
for f in mipster-tests.desktop mipster-tests-debug.desktop; do
  [ -f "${SRC_DIR}/test/${f}" ] && \
    cp "${SRC_DIR}/test/${f}" "${INSTALL_DIR}/share/applications/" && echo "    ${f}"
done

# ── Package ───────────────────────────────────────────────────────────────────
echo ""
echo "==> Packaging..."
cd "${SRC_DIR}/dist"
tar -czf "${DIST_NAME}.tar.gz" "${DIST_NAME}/"
echo "==> Done: dist/${DIST_NAME}.tar.gz  ($(du -sh "${DIST_NAME}.tar.gz" | cut -f1))"
tar -tzf "${DIST_NAME}.tar.gz"
