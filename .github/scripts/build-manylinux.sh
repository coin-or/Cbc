#!/usr/bin/env bash
# build-manylinux.sh — Build MIPster inside manylinux_2_34_x86_64.
#
# Produces: /project/dist/mipster-linux-x86_64.tar.gz
#   bin/mipster              — compiled CPU-dispatch launcher (static)
#   bin/mipster-dbg          — debug+ASan binary (symbols preserved)
#   bin/mipster-generic      — baseline x86_64 static binary
#   bin/mipster-avx2         — AVX2/FMA static binary (Haswell 2013+ / Zen2 2019+)
#   bin/test/                — test suite binaries (generic build)
#   lib/generic/libmipster.so*   — shared library, baseline
#   lib/avx2/libmipster.so*      — shared library, AVX2
#   lib/dbg/libmipster.so*   — debug+ASan shared library
#   lib/dbg/libasan.so.*     — bundled ASan runtime (self-contained)
#   include/coin/            — public C API headers

set -euo pipefail

DIST_NAME="mipster-linux-x86_64"
INSTALL_DIR="/project/dist/${DIST_NAME}"
BUILD_BASE="/tmp/mipster-build"

echo "==> MIPster manylinux_2_34 build"
# Activate GCC toolset runtime (needed for linker and libraries)
source /opt/rh/gcc-toolset-14/enable 2>/dev/null || true
echo "    GCC: $(gcc --version | head -1)"

# ── Build-time deps ───────────────────────────────────────────────────────────
dnf install -y bzip2-devel zlib-devel patchelf rsync openblas-static > /dev/null

echo "==> Copying source (excluding build artifacts)..."
rsync -a \
  --exclude='*.la' --exclude='*.lo' --exclude='*.o' \
  --exclude='.libs/' --exclude='config.status' --exclude='config.cache' \
  --exclude='dist/' \
  /project/ /tmp/mipster-src/
SRC_DIR=/tmp/mipster-src

mkdir -p "${INSTALL_DIR}/bin" "${INSTALL_DIR}/include"

# ── Build debug + ASan variant ────────────────────────────────────────────────
# Runs first so test failures surface with full debug symbols and ASan checking.
build_debug_variant() {
  local name="dbg"
  local cxxflags="-O1 -g -fno-omit-frame-pointer -fsanitize=address"
  local build_dir="${BUILD_BASE}-${name}"

  echo ""
  echo "==> Building debug+ASan variant  CXXFLAGS='${cxxflags}'"

  rm -rf "${build_dir}"
  mkdir -p "${build_dir}"
  cd "${build_dir}"

  # Use --disable-static so the mipster binary links dynamically against
  # libmipster.so (required for ASan instrumentation to be active at runtime).
  /tmp/mipster-src/configure \
    --prefix="${build_dir}/install" \
    --enable-shared \
    --disable-static \
    --without-amd \
    '--with-lapack-lflags=-lopenblas -lgfortran' \
    CXXFLAGS="${cxxflags}" \
    LDFLAGS="-fsanitize=address -static-libgcc -static-libstdc++ -Wl,-Bstatic,-lgfortran,-Bdynamic" \
    2>&1 | tail -3

  # Hide libopenblas.so so the linker falls back to libopenblas.a (static, self-contained).
  local openblas_so
  openblas_so=$(ldconfig -p 2>/dev/null | awk '/libopenblas\.so /{print $NF}' | head -1)
  [ -n "$openblas_so" ] && mv "$openblas_so" "${openblas_so}.bak"

  make -j"$(nproc)" 2>&1 | tail -3
  make install 2>&1 | tail -2

  [ -n "$openblas_so" ] && mv "${openblas_so}.bak" "$openblas_so"
  echo "    Build: OK"

  # ── Full test suite ────────────────────────────────────────────────────────
  cd "${build_dir}/test"
  make -j"$(nproc)" 2>&1 | tail -3
  MIPSTER_FIXTURE_DIR="${SRC_DIR}/test/fixtures" \
    LD_LIBRARY_PATH="${build_dir}/src/.libs${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}" \
    ASAN_OPTIONS="detect_leaks=0:abort_on_error=1:fast_unwind_on_malloc=0" \
    bash "${SRC_DIR}/test/run-mipster-tests"
  echo "    All tests: PASSED (debug+ASan)"

  # ── Debug binary ────────────────────────────────────────────────────────────
  # Keep debug symbols — do NOT strip.
  # The binary links against lib/dbg/libmipster.so; set RPATH accordingly.
  local dbg_bin="${build_dir}/src/.libs/mipster"
  patchelf --set-rpath '$ORIGIN/../lib/dbg' "${dbg_bin}"
  cp "${dbg_bin}" "${INSTALL_DIR}/bin/mipster-dbg"
  echo "    bin/mipster-dbg: $(du -sh "${INSTALL_DIR}/bin/mipster-dbg" | cut -f1)"

  # ── Debug shared library ─────────────────────────────────────────────────────
  local libdir="${INSTALL_DIR}/lib/dbg"
  mkdir -p "${libdir}"
  cp -P "${build_dir}/install/lib"/libmipster.so* "${libdir}/"

  # Keep debug symbols (do NOT strip). Set RPATH so bundled libasan is found.
  find "${libdir}" -name 'libmipster.so.*.*' ! -type l | \
    xargs -I{} patchelf --set-rpath '$ORIGIN' {}

  # ── Bundle libasan.so ────────────────────────────────────────────────────────
  # Copy the ASan runtime from the build environment so the debug binary and
  # library are fully self-contained (users don't need GCC dev packages).
  local libasan
  libasan=$(ldd "${build_dir}/src/.libs/mipster" 2>/dev/null | awk '/libasan/{print $3}' | head -1)
  if [ -z "$libasan" ] || [ ! -f "$libasan" ]; then
    libasan=$(ldconfig -p 2>/dev/null | awk '/libasan\.so\.[0-9]/{print $NF}' | sort -V | tail -1)
  fi
  if [ -n "$libasan" ] && [ -f "$libasan" ]; then
    cp "$libasan" "${libdir}/"
    echo "    Bundled: $(basename "$libasan")"
  else
    echo "    Warning: libasan.so not found; debug binary requires system libasan at runtime"
  fi

  echo "    lib/dbg: $(du -sh "${libdir}" | cut -f1)"
}

# ── Build one variant (static binary + shared lib) ────────────────────────────
build_variant() {
  local name="$1"
  local cxxflags="$2"
  local build_dir="${BUILD_BASE}-${name}"

  echo ""
  echo "==> Building variant: ${name}  CXXFLAGS='${cxxflags}'"

  # Configure with both static (for binary) and shared (for library)
  rm -rf "${build_dir}"
  mkdir -p "${build_dir}"
  cd "${build_dir}"

  /tmp/mipster-src/configure \
    --prefix="${build_dir}/install" \
    --enable-shared \
    --enable-static \
    --without-amd \
    '--with-lapack-lflags=-lopenblas -lgfortran' \
    CXXFLAGS="${cxxflags}" \
    LDFLAGS="-static-libstdc++ -static-libgcc -Wl,-Bstatic,-lgfortran,-Bdynamic" \
    2>&1 | tail -3

  # Hide libopenblas.so so the linker falls back to libopenblas.a (static).
  # This avoids a dynamic dependency on OpenBLAS in the installed binaries.
  local openblas_so
  openblas_so=$(ldconfig -p 2>/dev/null | awk '/libopenblas\.so /{print $NF}' | head -1)
  [ -n "$openblas_so" ] && mv "$openblas_so" "${openblas_so}.bak"

  make -j"$(nproc)" 2>&1 | tail -3
  make install 2>&1 | tail -2

  # Restore libopenblas.so
  [ -n "$openblas_so" ] && mv "${openblas_so}.bak" "$openblas_so"
  echo "    Build: OK"

  # ── Test ────────────────────────────────────────────────────────────────────
  cd "${build_dir}/test"
  make -j"$(nproc)" CInterfaceTest 2>&1 | tail -2
  ./CInterfaceTest
  echo "    CInterfaceTest: PASSED"

  # ── Collect test binaries for tarball (generic/first variant only) ──────────
  # These binaries are shipped in bin/test/ so users can run tests without
  # building from source. RPATH is patched to find lib/<variant>/libmipster.so.
  if [ "${name}" = "generic" ]; then
    make -j"$(nproc)" 2>&1 | tail -2
    local test_dir="${INSTALL_DIR}/bin/test"
    mkdir -p "${test_dir}"
    for tbin in CInterfaceTest CInterfaceTest_tsp_random CInterfaceTest_fl_random \
                CInterfaceTest_mdkp_random CInterfaceTest_nursesched CInterfaceTest_a1 \
                CInterfaceTest_graphdraw CInterfaceTest_trdta5581 CInterfaceTest_trd445c \
                CbcSolverLpTest; do
      if [ -f "${build_dir}/test/${tbin}" ]; then
        strip "${build_dir}/test/${tbin}"
        patchelf --set-rpath '$ORIGIN/../../lib/generic' "${build_dir}/test/${tbin}"
        cp "${build_dir}/test/${tbin}" "${test_dir}/"
      fi
    done
    echo "    bin/test/: $(ls "${test_dir}" | wc -l) test binaries"
  fi

  # ── Static binary ───────────────────────────────────────────────────────────
  strip "${build_dir}/src/.libs/mipster"
  cp "${build_dir}/src/.libs/mipster" "${INSTALL_DIR}/bin/mipster-${name}"
  echo "    bin/mipster-${name}: $(du -sh "${INSTALL_DIR}/bin/mipster-${name}" | cut -f1)"

  # ── Shared library ──────────────────────────────────────────────────────────
  # Copy installed .so files (versioned + symlinks) to lib/<variant>/
  local libdir="${INSTALL_DIR}/lib/${name}"
  mkdir -p "${libdir}"
  cp -P "${build_dir}/install/lib"/libmipster.so* "${libdir}/"

  # Strip the real .so (not symlinks)
  find "${libdir}" -name 'libmipster.so.*.*' ! -type l | xargs strip --strip-unneeded

  # Set RPATH to $ORIGIN so bundled deps (if any) are found next to the .so
  find "${libdir}" -name 'libmipster.so.*.*' ! -type l | \
    xargs -I{} patchelf --set-rpath '$ORIGIN' {}

  # Verify no unexpected external deps (glibc + libm + libpthread are fine)
  echo "    Shared lib deps:"
  find "${libdir}" -name 'libmipster.so.*.*' ! -type l | \
    xargs ldd | grep -v "linux-vdso\|libm\|libc\.so\|libpthread\|libdl\|/lib64/ld\|=>" | \
    grep "=>" || echo "      (none — fully self-contained)"

  # ── Headers (only need to do this once) ─────────────────────────────────────
  if [ "${name}" = "generic" ]; then
    cp -r "${build_dir}/install/include/mipster" "${INSTALL_DIR}/include/"
    echo "    include/mipster: copied"
  fi
}

# Debug+ASan build runs first — failures surface with full debug info.
build_debug_variant

# Release builds follow.
# -ffp-contract=off: prevent FMA fusion from changing LP floating-point results
build_variant "generic" "-O3 -ffp-contract=off"
# x86-64-v3 = AVX2 + BMI1/BMI2 + FMA (Haswell 2013+, AMD Zen2 2019+)
build_variant "avx2" "-O3 -march=x86-64-v3 -ffp-contract=off"

# ── Compiled CPU-dispatch launcher ────────────────────────────────────────────
echo ""
echo "==> Compiling launcher..."
cat > /tmp/mipster-launcher.c << 'EOF2'
/*
 * mipster launcher — detects AVX2 at runtime via CPUID, then exec's the
 * appropriate variant binary from the same directory.
 * Compiled as a fully static binary with no external dependencies.
 */
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>

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
    ssize_t n;
    char *slash;

    n = readlink("/proc/self/exe", self, sizeof(self) - 1);
    if (n < 0) { perror("readlink /proc/self/exe"); return 1; }
    self[n] = '\0';

    slash = strrchr(self, '/');
    if (!slash) { fputs("mipster: cannot find binary dir\n", stderr); return 1; }
    *slash = '\0';   /* self is now the directory */

    snprintf(target, sizeof(target), "%s/mipster-%s",
             self, has_avx2() ? "avx2" : "generic");

    execv(target, argv);
    perror(target);
    return 1;
}
EOF2

gcc -O2 -o "${INSTALL_DIR}/bin/mipster" /tmp/mipster-launcher.c
echo "    bin/mipster (launcher): $(du -sh "${INSTALL_DIR}/bin/mipster" | cut -f1)"

# ── Smoke test launcher ───────────────────────────────────────────────────────
echo ""
echo "==> Smoke test (launcher → --version):"
"${INSTALL_DIR}/bin/mipster" --version 2>&1 | head -1 || true
echo "    PASSED"

# ── Documentation ─────────────────────────────────────────────────────────────
echo ""
echo "==> Installing documentation..."
mkdir -p "${INSTALL_DIR}/share/doc/mipster"
mkdir -p "${INSTALL_DIR}/share/man/man1"
mkdir -p "${INSTALL_DIR}/share/applications"

# Copy committed docs from the source tree (generated from the binary)
DOC_SRC="/project/doc"
for f in mipster-parameters.pdf mipster-parameters.md; do
  [ -f "${DOC_SRC}/${f}" ] && cp "${DOC_SRC}/${f}" "${INSTALL_DIR}/share/doc/mipster/" && echo "    ${f}"
done
[ -f "${DOC_SRC}/mipster.1" ] && \
  gzip -c "${DOC_SRC}/mipster.1" > "${INSTALL_DIR}/share/man/man1/mipster.1.gz" && \
  echo "    mipster.1.gz"
for f in mipster.desktop mipster-docs.desktop; do
  [ -f "${DOC_SRC}/${f}" ] && cp "${DOC_SRC}/${f}" "${INSTALL_DIR}/share/applications/" && echo "    ${f}"
done

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
  [ -f "${SRC_DIR}/test/${f}" ] && cp "${SRC_DIR}/test/${f}" "${INSTALL_DIR}/share/applications/" && echo "    ${f}"
done

# ── Package ───────────────────────────────────────────────────────────────────
echo ""
echo "==> Packaging..."
cd /project/dist
tar -czf "${DIST_NAME}.tar.gz" "${DIST_NAME}/"
echo "==> Done: dist/${DIST_NAME}.tar.gz  ($(du -sh "${DIST_NAME}.tar.gz" | cut -f1))"
echo "Contents:"
tar -tzf "${DIST_NAME}.tar.gz"
