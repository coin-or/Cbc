#!/usr/bin/env bash
# build-macos-x86_64.sh — Build MIPster on macOS x86_64 (Intel).
#
# Produces: dist/mipster-macos-x86_64.tar.gz
#   bin/mipster              — CPU-dispatch launcher (generic vs haswell)
#   bin/mipster-generic      — baseline x86_64 binary
#   bin/mipster-haswell      — AVX2/FMA binary (Haswell 2013+)
#   lib/generic/libCbc.dylib — shared library, baseline
#   lib/haswell/libCbc.dylib — shared library, AVX2
#   include/coin-or/         — public C API headers
#
# Uses Apple's Accelerate framework for BLAS/LAPACK (no Homebrew deps).
# Requires Xcode Command Line Tools.

set -euo pipefail

DIST_NAME="mipster-macos-x86_64"
INSTALL_DIR="$(pwd)/dist/${DIST_NAME}"
BUILD_BASE="/tmp/mipster-build"
SRC_DIR="$(cd "$(dirname "$0")/../.." && pwd)"

# macOS 10.15 Catalina: last macOS to drop 32-bit, reasonable baseline
export MACOSX_DEPLOYMENT_TARGET=10.15

echo "==> MIPster macOS x86_64 build"
echo "    Clang: $(clang --version | head -1)"
echo "    Source: ${SRC_DIR}"

mkdir -p "${INSTALL_DIR}/bin" "${INSTALL_DIR}/include"

# ── Build one variant ─────────────────────────────────────────────────────────
build_variant() {
  local name="$1"
  local cxxflags="$2"
  local build_dir="${BUILD_BASE}-${name}"

  echo ""
  echo "==> Building variant: ${name}  CXXFLAGS='${cxxflags}'"

  rm -rf "${build_dir}"
  mkdir -p "${build_dir}"
  cd "${build_dir}"

  # Accelerate provides BLAS + LAPACK, always available on macOS
  "${SRC_DIR}/configure" \
    --prefix="${build_dir}/install" \
    --enable-shared \
    --enable-static \
    --without-amd \
    '--with-lapack-lflags=-framework Accelerate' \
    CXXFLAGS="${cxxflags}" \
    LDFLAGS="-mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}" \
    2>&1 | tail -3

  make -j"$(sysctl -n hw.logicalcpu)" 2>&1 | tail -3
  make install 2>&1 | tail -2
  echo "    Build: OK"

  # ── Test ────────────────────────────────────────────────────────────────────
  cd "${build_dir}/test"
  make -j"$(sysctl -n hw.logicalcpu)" CInterfaceTest 2>&1 | tail -2
  ./CInterfaceTest
  echo "    CInterfaceTest: PASSED"

  # ── Binary ──────────────────────────────────────────────────────────────────
  strip "${build_dir}/src/.libs/mipster"
  cp "${build_dir}/src/.libs/mipster" "${INSTALL_DIR}/bin/mipster-${name}"
  echo "    bin/mipster-${name}: $(du -sh "${INSTALL_DIR}/bin/mipster-${name}" | cut -f1)"

  # ── Shared library ──────────────────────────────────────────────────────────
  local libdir="${INSTALL_DIR}/lib/${name}"
  mkdir -p "${libdir}"
  # Copy versioned .dylib and symlinks
  cp -P "${build_dir}/install/lib"/libCbc*.dylib "${libdir}/" 2>/dev/null || true

  # Fix install name and RPATH so dylib is self-contained when placed next to the binary
  local real_dylib
  real_dylib=$(find "${libdir}" -name 'libCbc.*.dylib' ! -type l | head -1)
  if [ -n "${real_dylib}" ]; then
    install_name_tool -id "@loader_path/$(basename "${real_dylib}")" "${real_dylib}"
    strip -x "${real_dylib}"
  fi

  # Verify only system deps
  echo "    Shared lib deps:"
  if [ -n "${real_dylib}" ]; then
    otool -L "${real_dylib}" | grep -v "libCbc\|/usr/lib\|/System\|^${libdir}" \
      || echo "      (none — only system libs)"
  fi

  # ── Headers (once) ──────────────────────────────────────────────────────────
  if [ "${name}" = "generic" ]; then
    cp -r "${build_dir}/install/include/coin-or" "${INSTALL_DIR}/include/"
    echo "    include/coin-or: copied"
  fi
}

# -ffp-contract=off: prevent FMA fusion from changing LP floating-point results
build_variant "generic" "-O3 -ffp-contract=off"
# x86-64-v3 = AVX2 + BMI1/BMI2 + FMA (all Intel Macs from 2013+)
build_variant "haswell" "-O3 -march=x86-64-v3 -ffp-contract=off"

# ── CPU-dispatch launcher ─────────────────────────────────────────────────────
echo ""
echo "==> Compiling launcher..."
cat > /tmp/mipster-launcher.c << 'EOF'
/*
 * mipster launcher (macOS x86_64) — detects AVX2 via CPUID, then exec's
 * the appropriate variant binary from the same directory.
 */
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <stdint.h>
#include <mach-o/dyld.h>

static int has_avx2(void)
{
    unsigned int ebx;
    __asm__ volatile (
        "movl $7,  %%eax\n\t"
        "xorl %%ecx, %%ecx\n\t"
        "cpuid"
        : "=b"(ebx) : : "eax", "ecx", "edx"
    );
    return (ebx >> 5) & 1;  /* EBX bit 5 = AVX2 */
}

int main(int argc, char *argv[])
{
    char self[PATH_MAX], target[PATH_MAX];
    uint32_t size = sizeof(self);
    char *slash;

    if (_NSGetExecutablePath(self, &size) != 0) {
        fputs("mipster: executable path too long\n", stderr);
        return 1;
    }

    slash = strrchr(self, '/');
    if (!slash) { fputs("mipster: cannot find binary dir\n", stderr); return 1; }
    *slash = '\0';

    snprintf(target, sizeof(target), "%s/mipster-%s",
             self, has_avx2() ? "haswell" : "generic");

    execv(target, argv);
    perror(target);
    return 1;
}
EOF

clang -O2 -mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET} \
  -o "${INSTALL_DIR}/bin/mipster" /tmp/mipster-launcher.c
echo "    bin/mipster (launcher): $(du -sh "${INSTALL_DIR}/bin/mipster" | cut -f1)"

# ── Smoke test ────────────────────────────────────────────────────────────────
echo ""
echo "==> Smoke test (launcher → --version):"
"${INSTALL_DIR}/bin/mipster" --version 2>&1 | head -1
echo "    PASSED"

# ── Package ───────────────────────────────────────────────────────────────────
echo ""
echo "==> Packaging..."
cd "$(pwd)/dist"
tar -czf "${DIST_NAME}.tar.gz" "${DIST_NAME}/"
echo "==> Done: dist/${DIST_NAME}.tar.gz  ($(du -sh "${DIST_NAME}.tar.gz" | cut -f1))"
tar -tzf "${DIST_NAME}.tar.gz"
