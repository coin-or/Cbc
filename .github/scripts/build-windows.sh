#!/usr/bin/env bash
# build-windows.sh — Build MIPster on Windows using MSYS2/MinGW64.
#
# Run via: C:\msys64\usr\bin\bash.exe -l .github/scripts/build-windows.sh
#
# Produces:
#   dist/mipster-windows-x86_64.zip              — portable archive
#   dist/mipster-<ver>-setup-windows-x86_64.exe  — NSIS installer (updates PATH)
#
# Contents:
#   bin/mipster.exe              — CPU-dispatch launcher (CPUID → generic or avx2)
#   bin/mipster-dbg.exe          — debug binary (symbols preserved)
#   bin/mipster-generic.exe      — baseline x86_64 binary (static, no DLL deps)
#   bin/mipster-avx2.exe         — AVX2/FMA binary, Haswell 2013+ / Zen2 2019+
#   bin/test/                    — test suite binaries (generic build)
#   bin/libmipster-dbg.dll       — debug shared library
#   lib/generic/libmipster.dll.a — import library, baseline
#   lib/avx2/libmipster.dll.a    — import library, AVX2
#   include/mipster/             — public C API headers
#
# No OpenBLAS dependency: built without BLAS/LAPACK.
# Binaries are fully self-contained: libgcc, libstdc++, libwinpthread are
# linked statically so no extra DLLs need to be distributed.
#
# Requirements: MSYS2 with mingw-w64-x86_64-gcc, make, autoconf, automake,
#               libtool, pkg-config, nsis, zip (installed by the CI step).

set -euo pipefail

export PATH="/mingw64/bin:/usr/bin:$PATH"

DIST_NAME="mipster-windows-x86_64"
# Script lives at .github/scripts/; source root is two levels up.
SRC_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
INSTALL_DIR="${SRC_DIR}/dist/${DIST_NAME}"
BUILD_BASE="/tmp/mipster-build-win"

echo "==> MIPster Windows/MinGW64 build"
echo "    GCC:    $(gcc --version | head -1)"
echo "    Source: ${SRC_DIR}"

# bin/ holds both executables AND the DLLs they depend on, so that everything
# works from any directory on a clean Windows machine without PATH changes.
mkdir -p "${INSTALL_DIR}/bin" "${INSTALL_DIR}/include"

# ── Build debug variant ────────────────────────────────────────────────────────
# Note: ASan is not supported with MinGW. Debug binary uses -O1 -g only.
build_debug_variant() {
  local name="dbg"
  local cxxflags="-O1 -g -fno-omit-frame-pointer"
  local build_dir="${BUILD_BASE}-${name}"

  echo ""
  echo "==> Building debug variant: ${name}  CXXFLAGS='${cxxflags}'"

  rm -rf "${build_dir}"
  mkdir -p "${build_dir}"
  cd "${build_dir}"

  "${SRC_DIR}/configure" \
    --prefix="${build_dir}/install" \
    --enable-shared \
    --disable-static \
    --without-amd \
    --without-lapack \
    --disable-bzlib \
    CXXFLAGS="${cxxflags}" \
    LDFLAGS="-static-libgcc -static-libstdc++ -Wl,-Bstatic,-lwinpthread,-Bdynamic -Wl,--export-all-symbols" \
    2>&1 | tail -3

  make -j"$(nproc)" 2>&1 | tail -3
  make install 2>&1 | tail -2
  echo "    Build: OK"

  # ── Test ──────────────────────────────────────────────────────────────────
  cd "${build_dir}/test"
  make -j"$(nproc)" CInterfaceTest.exe 2>&1 | tail -2
  ./CInterfaceTest.exe
  echo "    CInterfaceTest: PASSED (debug)"

  # ── Debug executable ───────────────────────────────────────────────────────
  # No strip — keep debug symbols.
  cp "${build_dir}/src/.libs/mipster.exe" "${INSTALL_DIR}/bin/mipster-dbg.exe"
  echo "    bin/mipster-dbg.exe: $(du -sh "${INSTALL_DIR}/bin/mipster-dbg.exe" | cut -f1)"

  # ── Debug DLL ─────────────────────────────────────────────────────────────
  find "${build_dir}/install/bin" -name 'libmipster*.dll' | while read -r f; do
    # No strip — keep debug symbols.
    local base
    base=$(basename "${f}")
    cp "${f}" "${INSTALL_DIR}/bin/${base%.dll}-${name}.dll"
  done
  echo "    bin/libmipster-${name}.dll: done"
}

# ── Build one variant (exe + DLL) ─────────────────────────────────────────────
build_variant() {
  local name="$1"
  local cxxflags="$2"
  local build_dir="${BUILD_BASE}-${name}"

  echo ""
  echo "==> Building variant: ${name}  CXXFLAGS='${cxxflags}'"

  rm -rf "${build_dir}"
  mkdir -p "${build_dir}"
  cd "${build_dir}"

  # Shared build (DLL): CoinUtils/Clp/Cgl are noinst convenience libs so libtool
  # merges them all into libmipster.dll — one self-contained DLL.
  # MinGW runtime (libgcc, libstdc++, libwinpthread) is embedded statically via
  # LDFLAGS so libmipster.dll and the exe need only standard Windows system DLLs.
  "${SRC_DIR}/configure" \
    --prefix="${build_dir}/install" \
    --enable-shared \
    --disable-static \
    --without-amd \
    --without-lapack \
    --disable-bzlib \
    CXXFLAGS="${cxxflags}" \
    LDFLAGS="-static-libgcc -static-libstdc++ -Wl,-Bstatic,-lwinpthread,-Bdynamic -Wl,--export-all-symbols" \
    2>&1 | tail -3

  make -j"$(nproc)"
  make install 2>&1 | tail -2
  echo "    Build: OK"

  # ── Test ──────────────────────────────────────────────────────────────────
  cd "${build_dir}/test"
  make -j"$(nproc)" CInterfaceTest.exe
  ./CInterfaceTest.exe
  echo "    CInterfaceTest: PASSED"

  # ── Collect test binaries for tarball (generic variant only) ──────────────
  if [ "${name}" = "generic" ]; then
    make -j"$(nproc)"
    local test_dir="${INSTALL_DIR}/bin/test"
    mkdir -p "${test_dir}"
    for tbin in CInterfaceTest CInterfaceTest_tsp_random CInterfaceTest_fl_random \
                CInterfaceTest_mdkp_random CInterfaceTest_nursesched CInterfaceTest_a1 \
                CInterfaceTest_graphdraw CInterfaceTest_trdta5581 CInterfaceTest_trd445c \
                CbcSolverLpTest; do
      if [ -f "${build_dir}/test/${tbin}.exe" ]; then
        strip "${build_dir}/test/${tbin}.exe"
        cp "${build_dir}/test/${tbin}.exe" "${test_dir}/"
      fi
    done
    echo "    bin/test/: $(ls "${test_dir}" | wc -l) test binaries"
  fi

  # ── Executable ────────────────────────────────────────────────────────────
  # On Windows, libtool places the real exe at src/.libs/mipster.exe
  strip "${build_dir}/src/.libs/mipster.exe"
  cp "${build_dir}/src/.libs/mipster.exe" "${INSTALL_DIR}/bin/mipster-${name}.exe"
  echo "    bin/mipster-${name}.exe: $(du -sh "${INSTALL_DIR}/bin/mipster-${name}.exe" | cut -f1)"

  # ── Shared library (.dll) ─────────────────────────────────────────────────
  # libtool installs DLLs into $prefix/bin on Windows. Copy libmipster-N.dll
  # alongside the executables so they are found without any PATH change.
  find "${build_dir}/install/bin" -name 'libmipster*.dll' | while read -r f; do
    strip --strip-unneeded "${f}" 2>/dev/null || true
    # Suffix dll with variant name to avoid collision between generic/avx2 builds
    local base
    base=$(basename "${f}")
    cp "${f}" "${INSTALL_DIR}/bin/${base%.dll}-${name}.dll"
  done
  # Also copy the import lib (.dll.a) to lib/ for developers linking against the DLL
  mkdir -p "${INSTALL_DIR}/lib"
  find "${build_dir}/install/lib" -name 'libmipster*.dll.a' | while read -r f; do
    cp "${f}" "${INSTALL_DIR}/lib/"
  done
  echo "    bin/libmipster-${name}.dll: done"
}

build_debug_variant
build_variant "generic" "-O3 -ffp-contract=off"
build_variant "avx2"    "-O3 -march=x86-64-v3 -ffp-contract=off"

# ── Headers ───────────────────────────────────────────────────────────────────
cp -r "${BUILD_BASE}-generic/install/include/mipster" "${INSTALL_DIR}/include/"
echo "    include/mipster: copied"

# ── Compile CPU-dispatch launcher ─────────────────────────────────────────────
echo ""
echo "==> Compiling launcher..."
cat > /tmp/mipster_launcher.c << 'EOF2'
/*
 * mipster_launcher.c — Windows CPU-dispatch launcher for MIPster.
 *
 * Uses CPUID leaf 7, EBX bit 5 to detect AVX2 support, then exec-replaces
 * itself with mipster-avx2.exe or mipster-generic.exe from the same directory.
 */
#include <windows.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <process.h>

static int has_avx2(void)
{
    unsigned int ebx = 0;
    __asm__ volatile (
        "movl $7, %%eax\n\t"
        "xorl %%ecx, %%ecx\n\t"
        "cpuid"
        : "=b"(ebx) : : "eax", "ecx", "edx"
    );
    return (ebx >> 5) & 1;
}

int main(int argc, char *argv[])
{
    char exe_path[MAX_PATH];
    char dir[MAX_PATH];
    char target[MAX_PATH];

    if (!GetModuleFileNameA(NULL, exe_path, MAX_PATH)) {
        fprintf(stderr, "mipster: GetModuleFileNameA failed (%lu)\n", GetLastError());
        return 1;
    }

    char *sep = strrchr(exe_path, '\\');
    if (!sep) sep = strrchr(exe_path, '/');
    if (sep) {
        size_t n = (size_t)(sep - exe_path);
        memcpy(dir, exe_path, n);
        dir[n] = '\0';
    } else {
        dir[0] = '.'; dir[1] = '\0';
    }

    const char *variant = has_avx2() ? "mipster-avx2.exe" : "mipster-generic.exe";
    snprintf(target, MAX_PATH, "%s\\%s", dir, variant);

    _execv(target, (const char * const *)argv);
    fprintf(stderr, "mipster: could not exec '%s': %s\n", target, strerror(errno));
    return 1;
}
EOF2

gcc -O2 -o "${INSTALL_DIR}/bin/mipster.exe" /tmp/mipster_launcher.c
strip "${INSTALL_DIR}/bin/mipster.exe"
echo "    bin/mipster.exe (launcher): $(du -sh "${INSTALL_DIR}/bin/mipster.exe" | cut -f1)"

# ── Smoke test ────────────────────────────────────────────────────────────────
echo ""
echo "==> Smoke test (launcher → --version):"
"${INSTALL_DIR}/bin/mipster.exe" --version 2>&1 | head -1 || true
echo "    PASSED"

# ── Verify no unexpected DLL dependencies ─────────────────────────────────────
echo ""
echo "==> DLL dependencies (mipster-generic.exe):"
objdump -p "${INSTALL_DIR}/bin/mipster-generic.exe" 2>/dev/null \
  | grep 'DLL Name' \
  | grep -iv 'kernel32\|ntdll\|msvcrt\|user32\|advapi32\|ws2_32\|libmipster\|api-ms' \
  || echo "    (only system DLLs + libmipster — OK)"
echo ""
echo "==> DLL dependencies (libmipster-N-generic.dll):"
find "${INSTALL_DIR}/bin" -name 'libmipster*-generic.dll' | head -1 | xargs objdump -p 2>/dev/null \
  | grep 'DLL Name' \
  | grep -iv 'kernel32\|ntdll\|msvcrt\|user32\|advapi32\|ws2_32\|api-ms' \
  || echo "    (only system DLLs — self-contained OK)"

# ── NSIS installer ────────────────────────────────────────────────────────────
echo ""
echo "==> Building NSIS installer..."

# Determine version (from configure.ac or a fallback).
# AC_INIT([Name],[version],...) — version is the second bracketed argument.
VERSION=$(grep '^AC_INIT' "${SRC_DIR}/configure.ac" 2>/dev/null \
  | sed 's/AC_INIT(\[[^]]*\],\[\([^]]*\)\].*/\1/' | head -1)
VERSION=${VERSION:-devel}
echo "    Version: ${VERSION}"

# Substitute placeholders in the .nsi template.
NSIS_SCRIPT="/tmp/mipster-installer.nsi"
# cygpath -m gives Windows-style paths with forward slashes (e.g. D:/a/foo/dist).
# Forward slashes are safe in sed replacements (no backslash escaping needed)
# and NSIS accepts them in File/SetOutPath commands.
DIST_WIN=$(cygpath -m "${SRC_DIR}/dist" 2>/dev/null || echo "${SRC_DIR}/dist")
DOC_WIN=$(cygpath -m "${SRC_DIR}/doc" 2>/dev/null || echo "${SRC_DIR}/doc")

sed \
  -e "s|@VERSION@|${VERSION}|g" \
  -e "s|@DIST_DIR@|${DIST_WIN}|g" \
  -e "s|@DIST_SUBDIR@|${DIST_NAME}|g" \
  -e "s|@DOC_DIR@|${DOC_WIN}|g" \
  "${SRC_DIR}/.github/scripts/mipster-installer.nsi" > "${NSIS_SCRIPT}"

# makensis is available from the NSIS package in MSYS2
makensis "${NSIS_SCRIPT}"

INSTALLER_FILE="${SRC_DIR}/dist/mipster-${VERSION}-setup-windows-x86_64.exe"
if [ -f "${INSTALLER_FILE}" ]; then
  echo "    Installer: $(du -sh "${INSTALLER_FILE}" | cut -f1)"
else
  echo "    Warning: installer not found at expected path, searching..."
  find "${SRC_DIR}/dist" -name '*.exe' | head -5
fi

# ── Package as zip ────────────────────────────────────────────────────────────
echo ""
echo "==> Packaging..."
cd "${SRC_DIR}/dist"
zip -r "${DIST_NAME}.zip" "${DIST_NAME}/"
echo "==> Done: dist/${DIST_NAME}.zip  ($(du -sh "${DIST_NAME}.zip" | cut -f1))"
echo "    Installer: dist/mipster-${VERSION}-setup-windows-x86_64.exe"
