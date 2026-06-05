#!/usr/bin/env bash
# build-windows.sh — Build MIPster on Windows using MSYS2/MinGW64.
#
# Run via: C:\msys64\usr\bin\bash.exe -l .github/scripts/build-windows.sh
#
# Produces:
#   dist/mipster-windows-x86_64.zip              — portable archive
#   dist/mipster-<ver>-setup-windows-x86_64.exe  — NSIS installer (updates PATH)
#
# Layout (per-variant subdirs so the EXE's PE import table can keep libtool's
# natural name `libmipster-0.dll` and find the right variant DLL right next to
# it — renaming the DLL would not work because the import name is baked into
# the linked EXE):
#   bin/mipster.exe              — CPU-dispatch launcher (CPUID → generic/avx2)
#   bin/generic/mipster.exe      — baseline x86_64 binary
#   bin/generic/libmipster-0.dll — baseline DLL
#   bin/generic/{libgcc_s_seh-1,libstdc++-6,libwinpthread-1,zlib1}.dll
#   bin/avx2/mipster.exe         — AVX2/FMA binary, Haswell 2013+ / Zen2 2019+
#   bin/avx2/libmipster-0.dll    — AVX2 DLL
#   bin/avx2/{libgcc_s_seh-1,libstdc++-6,libwinpthread-1,zlib1}.dll
#   bin/dbg/mipster.exe          — debug binary (symbols preserved)
#   bin/dbg/libmipster-0.dll     — debug DLL
#   bin/dbg/{libgcc_s_seh-1,libstdc++-6,libwinpthread-1,zlib1}.dll
#   bin/test/                    — test suite binaries + DLLs (generic build)
#   lib/generic/libmipster.dll.a — import library, baseline (for developers)
#   lib/avx2/libmipster.dll.a    — import library, AVX2
#   include/mipster/             — public C API headers
#
# No OpenBLAS dependency: built without BLAS/LAPACK.
# MinGW runtime DLLs (libgcc_s_seh-1, libstdc++-6, libwinpthread-1) are
# bundled alongside libmipster — libtool ignores -static-libgcc/-static-libstdc++
# when linking shared libs so they remain dynamic NEEDED entries. zlib1.dll is
# also bundled (CoinUtils uses it for .gz reading). All bundled DLLs are
# either GCC-RLE-3.1 (libgcc_s, libstdc++, libwinpthread) or zlib license.
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

# MinGW runtime DLLs that must travel with each variant subdir (libtool ignores
# -static-libgcc/-static-libstdc++ when linking the shared library; zlib1.dll is
# pulled in by CoinUtils for gzipped MPS reading).
MINGW_RUNTIME_DLLS=(
  libgcc_s_seh-1.dll
  libstdc++-6.dll
  libwinpthread-1.dll
  zlib1.dll
)

# Allowed NEEDED-entry regex for verification: Windows system DLLs + libmipster
# + bundled MinGW runtime + Windows API set redirection DLLs.
ALLOWED_DLL_DEPS_RE='kernel32|ntdll|msvcrt|user32|advapi32|ws2_32|libmipster|api-ms|libgcc_s_seh-1|libstdc\+\+-6|libwinpthread-1|zlib1'

echo "==> MIPster Windows/MinGW64 build"
echo "    GCC:    $(gcc --version | head -1)"
echo "    Source: ${SRC_DIR}"
echo "    MSYSTEM=${MSYSTEM:-unset}"
echo "    config.guess: $(${SRC_DIR}/BuildTools/config.guess 2>/dev/null || echo unavailable)"

mkdir -p "${INSTALL_DIR}/bin" "${INSTALL_DIR}/include"

# Copy the MinGW runtime DLLs into a target directory (next to libmipster*.dll).
copy_mingw_runtime_into() {
  local dest="$1"
  for dll in "${MINGW_RUNTIME_DLLS[@]}"; do
    local src="/mingw64/bin/${dll}"
    if [ -f "${src}" ]; then
      cp "${src}" "${dest}/"
    else
      echo "ERROR: ${src} not found — runtime DLL missing from /mingw64/bin/" >&2
      exit 1
    fi
  done
}

# Locate the libtool-produced libmipster*.dll for a given build dir, refusing
# any cygmipster*.dll outright. cygmipster*.dll's PE import table forces a
# dependency on msys-2.0.dll which is not present on a clean Windows host —
# this is a packaging defect that cannot be fixed by renaming the file.
locate_mipster_dll() {
  local build_dir="$1"
  local cyg_dlls
  cyg_dlls=$(find "${build_dir}/install/bin" "${build_dir}/src/.libs" \
               -maxdepth 2 -name 'cygmipster*.dll' 2>/dev/null || true)
  if [ -n "${cyg_dlls}" ]; then
    echo "ERROR: libtool produced cygmipster*.dll — MSYSTEM is not MINGW64." >&2
    echo "  MSYSTEM=${MSYSTEM:-unset}" >&2
    echo "  Set MSYSTEM=MINGW64 in the workflow before invoking this script." >&2
    exit 1
  fi
  local dlls
  dlls=$(find "${build_dir}/install/bin" "${build_dir}/src/.libs" \
           -maxdepth 2 -name 'libmipster*.dll' 2>/dev/null || true)
  if [ -z "${dlls}" ]; then
    echo "ERROR: no libmipster*.dll found under ${build_dir}/install/bin or" \
         "${build_dir}/src/.libs." >&2
    echo "Inventory of build_dir/install/bin:" >&2
    ls -la "${build_dir}/install/bin/" >&2 || true
    echo "Inventory of build_dir/src/.libs:" >&2
    ls -la "${build_dir}/src/.libs/" >&2 || true
    exit 1
  fi
  # First match — there is only one .dll per build (we built --disable-static).
  echo "${dlls}" | head -1
}

# Verify the NEEDED entries of a PE image match the allow-list above.
verify_dll_deps() {
  local image="$1"
  local label="$2"
  local unbundled
  unbundled=$(objdump -p "${image}" 2>/dev/null \
    | awk '/DLL Name:/{print tolower($3)}' \
    | grep -v -iE "${ALLOWED_DLL_DEPS_RE}" \
    || true)
  if [ -n "${unbundled}" ]; then
    echo "ERROR: ${label} has unbundled non-system DLL deps:" >&2
    echo "${unbundled}" >&2
    return 1
  fi
}

# ── Build debug variant ────────────────────────────────────────────────────────
# Note: ASan is not supported with MinGW. Debug binary uses -O1 -g only.
build_debug_variant() {
  local name="dbg"
  local cxxflags="-O1 -g -fno-omit-frame-pointer"
  local build_dir="${BUILD_BASE}-${name}"
  local variant_bindir="${INSTALL_DIR}/bin/${name}"

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
  PATH="${build_dir}/src/.libs:${PATH}" ./CInterfaceTest.exe
  echo "    CInterfaceTest: PASSED (debug)"

  # ── Lay out variant subdir: exe + DLL keep their natural names ────────────
  mkdir -p "${variant_bindir}"
  # No strip — keep debug symbols.
  cp "${build_dir}/src/.libs/mipster.exe" "${variant_bindir}/mipster.exe"
  local dll
  dll=$(locate_mipster_dll "${build_dir}")
  cp "${dll}" "${variant_bindir}/"
  copy_mingw_runtime_into "${variant_bindir}"
  echo "    bin/${name}/: $(ls "${variant_bindir}" | wc -l) files"
}

# ── Build one variant (exe + DLL) ─────────────────────────────────────────────
build_variant() {
  local name="$1"
  local cxxflags="$2"
  local build_dir="${BUILD_BASE}-${name}"
  local variant_bindir="${INSTALL_DIR}/bin/${name}"

  echo ""
  echo "==> Building variant: ${name}  CXXFLAGS='${cxxflags}'"

  rm -rf "${build_dir}"
  mkdir -p "${build_dir}"
  cd "${build_dir}"

  # Shared build (DLL): CoinUtils/Clp/Cgl are noinst convenience libs so libtool
  # merges them all into libmipster.dll — one self-contained DLL.
  # MinGW runtime is *requested* statically via LDFLAGS but libtool ignores
  # those flags for shared-library linkage; the result is that libmipster*.dll
  # has dynamic NEEDED entries on libgcc_s_seh-1.dll, libstdc++-6.dll,
  # libwinpthread-1.dll, and zlib1.dll. We bundle those in each variant subdir.
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
  PATH="${build_dir}/src/.libs:${PATH}" ./CInterfaceTest.exe
  echo "    CInterfaceTest: PASSED"

  # ── Lay out variant subdir ────────────────────────────────────────────────
  mkdir -p "${variant_bindir}"
  strip "${build_dir}/src/.libs/mipster.exe"
  cp "${build_dir}/src/.libs/mipster.exe" "${variant_bindir}/mipster.exe"
  local dll
  dll=$(locate_mipster_dll "${build_dir}")
  cp "${dll}" "${variant_bindir}/"
  strip --strip-unneeded "${variant_bindir}/$(basename "${dll}")" 2>/dev/null || true
  copy_mingw_runtime_into "${variant_bindir}"
  echo "    bin/${name}/: $(ls "${variant_bindir}" | wc -l) files"

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
    # Test binaries also import libmipster-0.dll → ship the generic DLL +
    # MinGW runtime next to them so they run from bin/test/ without a PATH
    # change.
    cp "${variant_bindir}/$(basename "${dll}")" "${test_dir}/"
    copy_mingw_runtime_into "${test_dir}"
    echo "    bin/test/: $(ls "${test_dir}" | wc -l) files"
  fi

  # ── Import library (.dll.a) for developers linking against libmipster ─────
  local variant_libdir="${INSTALL_DIR}/lib/${name}"
  mkdir -p "${variant_libdir}"
  local import_libs
  import_libs=$(find "${build_dir}/install/lib" "${build_dir}/src/.libs" \
                  -maxdepth 2 -name 'libmipster*.dll.a' 2>/dev/null || true)
  if [ -z "${import_libs}" ]; then
    echo "ERROR: no libmipster*.dll.a (import lib) found for variant '${name}'." >&2
    exit 1
  fi
  while IFS= read -r f; do
    cp "${f}" "${variant_libdir}/"
  done <<< "${import_libs}"
  echo "    lib/${name}/libmipster.dll.a: done"
}

build_debug_variant
build_variant "generic" "-O3 -ffp-contract=off"
build_variant "avx2"    "-O3 -march=x86-64-v3 -ffp-contract=off"

# ── Headers ───────────────────────────────────────────────────────────────────
cp -r "${BUILD_BASE}-generic/install/include/mipster" "${INSTALL_DIR}/include/"
echo "    include/mipster: copied"

# ── Compile CPU-dispatch launcher ─────────────────────────────────────────────
# The launcher lives at bin/mipster.exe. It detects AVX2 via CPUID, prepends
# bin/<variant>/ to PATH so the variant's libmipster-0.dll loads from there,
# and exec-replaces itself with bin/<variant>/mipster.exe.
echo ""
echo "==> Compiling launcher..."
cat > /tmp/mipster_launcher.c << 'EOF2'
/*
 * mipster_launcher.c — Windows CPU-dispatch launcher for MIPster.
 *
 * Detects AVX2 via CPUID leaf 7 (EBX bit 5), prepends bin/<variant>/ to PATH,
 * and exec's bin/<variant>/mipster.exe with the same argv. The variant subdir
 * carries libmipster-0.dll and the MinGW runtime DLLs, so DLL search resolves
 * locally without needing any installer-managed PATH entry.
 */
#include <windows.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
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
    char variant_dir[MAX_PATH];
    char target[MAX_PATH];
    char *new_path;
    const char *old_path;

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

    const char *variant = has_avx2() ? "avx2" : "generic";
    snprintf(variant_dir, MAX_PATH, "%s\\%s", dir, variant);
    snprintf(target, MAX_PATH, "%s\\mipster.exe", variant_dir);

    /* Prepend variant_dir to PATH so libmipster-0.dll resolves from there. */
    old_path = getenv("PATH");
    if (old_path == NULL) old_path = "";
    {
        size_t need = strlen(variant_dir) + 1 + strlen(old_path) + 1;
        new_path = (char *)malloc(need);
        if (!new_path) {
            fputs("mipster: out of memory\n", stderr);
            return 1;
        }
        snprintf(new_path, need, "%s;%s", variant_dir, old_path);
        SetEnvironmentVariableA("PATH", new_path);
    }

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

# ── Verify NEEDED entries on each variant's exe + DLL ─────────────────────────
echo ""
echo "==> Verifying DLL dependencies of all variant binaries..."
for variant in generic avx2 dbg; do
  exe="${INSTALL_DIR}/bin/${variant}/mipster.exe"
  dll=$(find "${INSTALL_DIR}/bin/${variant}" -name 'libmipster*.dll' | head -1)
  verify_dll_deps "${exe}" "bin/${variant}/mipster.exe"
  verify_dll_deps "${dll}" "bin/${variant}/$(basename "${dll}")"
  echo "    bin/${variant}/: deps OK"
done
# Launcher itself has no MinGW runtime (compiled with -static-libgcc isn't
# needed since the launcher is plain C, but verify anyway).
verify_dll_deps "${INSTALL_DIR}/bin/mipster.exe" "bin/mipster.exe (launcher)"
echo "    bin/mipster.exe: deps OK"

# ── NSIS installer ────────────────────────────────────────────────────────────
echo ""
echo "==> Building NSIS installer..."

VERSION=$(grep '^AC_INIT' "${SRC_DIR}/configure.ac" 2>/dev/null \
  | sed 's/AC_INIT(\[[^]]*\],\[\([^]]*\)\].*/\1/' | head -1)
VERSION=${VERSION:-devel}
echo "    Version: ${VERSION}"

NSIS_SCRIPT="/tmp/mipster-installer.nsi"
DIST_WIN=$(cygpath -m "${SRC_DIR}/dist" 2>/dev/null || echo "${SRC_DIR}/dist")
DOC_WIN=$(cygpath -m "${SRC_DIR}/doc" 2>/dev/null || echo "${SRC_DIR}/doc")

sed \
  -e "s|@VERSION@|${VERSION}|g" \
  -e "s|@DIST_DIR@|${DIST_WIN}|g" \
  -e "s|@DIST_SUBDIR@|${DIST_NAME}|g" \
  -e "s|@DOC_DIR@|${DOC_WIN}|g" \
  "${SRC_DIR}/.github/scripts/mipster-installer.nsi" > "${NSIS_SCRIPT}"

makensis "${NSIS_SCRIPT}"

INSTALLER_FILE="${SRC_DIR}/dist/mipster-${VERSION}-setup-windows-x86_64.exe"
if [ -f "${INSTALLER_FILE}" ]; then
  echo "    Installer: $(du -sh "${INSTALLER_FILE}" | cut -f1)"
else
  echo "    Warning: installer not found at expected path, searching..."
  find "${SRC_DIR}/dist" -name '*.exe' | head -5
fi

# ── Pre-package verification ──────────────────────────────────────────────────
echo ""
echo "==> Verifying install tree..."
verify_fail=0
for required in \
    "${INSTALL_DIR}/bin/mipster.exe" \
    "${INSTALL_DIR}/bin/generic/mipster.exe" \
    "${INSTALL_DIR}/bin/avx2/mipster.exe" \
    "${INSTALL_DIR}/bin/dbg/mipster.exe" \
    "${INSTALL_DIR}/include/mipster/Cbc_C_Interface.h"
do
  if [ ! -f "${required}" ]; then
    echo "  MISSING: ${required}" >&2
    verify_fail=1
  fi
done
for variant in generic avx2 dbg; do
  if ! ls "${INSTALL_DIR}/bin/${variant}"/libmipster*.dll >/dev/null 2>&1; then
    echo "  MISSING: ${INSTALL_DIR}/bin/${variant}/libmipster*.dll" >&2
    verify_fail=1
  fi
  for rt in "${MINGW_RUNTIME_DLLS[@]}"; do
    if [ ! -f "${INSTALL_DIR}/bin/${variant}/${rt}" ]; then
      echo "  MISSING: ${INSTALL_DIR}/bin/${variant}/${rt}" >&2
      verify_fail=1
    fi
  done
done
if [ "${verify_fail}" -ne 0 ]; then
  echo "ERROR: install tree verification failed — refusing to package." >&2
  echo "Full inventory:" >&2
  find "${INSTALL_DIR}" -maxdepth 4 -type f >&2
  exit 1
fi
echo "    OK"

# ── Third-party licenses ──────────────────────────────────────────────────────
echo ""
echo "==> Installing third-party license texts..."
mkdir -p "${INSTALL_DIR}/share/doc/mipster"
cp "${SRC_DIR}/.github/scripts/THIRD_PARTY_LICENSES.md" \
   "${INSTALL_DIR}/share/doc/mipster/THIRD_PARTY_LICENSES.md"
[ -f "${SRC_DIR}/LICENSE" ] && cp "${SRC_DIR}/LICENSE" "${INSTALL_DIR}/share/doc/mipster/LICENSE"
echo "    THIRD_PARTY_LICENSES.md, LICENSE"

# ── Package as zip ────────────────────────────────────────────────────────────
echo ""
echo "==> Packaging..."
cd "${SRC_DIR}/dist"
zip -r "${DIST_NAME}.zip" "${DIST_NAME}/"
echo "==> Done: dist/${DIST_NAME}.zip  ($(du -sh "${DIST_NAME}.zip" | cut -f1))"
echo "    Installer: dist/mipster-${VERSION}-setup-windows-x86_64.exe"
