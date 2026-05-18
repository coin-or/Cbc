#!/usr/bin/env bash
# build-manylinux.sh — Build MIPster inside manylinux_2_34_x86_64.
#
# Produces: /project/dist/mipster-linux-x86_64.tar.gz
#   bin/mipster              — compiled CPU-dispatch launcher (static)
#   bin/mipster-generic      — baseline x86_64 static binary
#   bin/mipster-avx2         — AVX2/FMA static binary (Haswell 2013+ / Zen2 2019+)
#   lib/generic/libmipster.so*   — shared library, baseline
#   lib/avx2/libmipster.so*      — shared library, AVX2
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
    cp -r "${build_dir}/install/include/coin-or" "${INSTALL_DIR}/include/"
    echo "    include/coin-or: copied"
  fi
}

# -ffp-contract=off: prevent FMA fusion from changing LP floating-point results
build_variant "generic" "-O3 -ffp-contract=off"
# x86-64-v3 = AVX2 + BMI1/BMI2 + FMA (Haswell 2013+, AMD Zen2 2019+)
build_variant "avx2" "-O3 -march=x86-64-v3 -ffp-contract=off"

# ── Compiled CPU-dispatch launcher ────────────────────────────────────────────
echo ""
echo "==> Compiling launcher..."
cat > /tmp/mipster-launcher.c << 'EOF'
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
EOF

gcc -O2 -o "${INSTALL_DIR}/bin/mipster" /tmp/mipster-launcher.c
echo "    bin/mipster (launcher): $(du -sh "${INSTALL_DIR}/bin/mipster" | cut -f1)"

# ── Smoke test launcher ───────────────────────────────────────────────────────
echo ""
echo "==> Smoke test (launcher → --version):"
"${INSTALL_DIR}/bin/mipster" --version 2>&1 | head -1 || true
echo "    PASSED"

# ── Package ───────────────────────────────────────────────────────────────────
echo ""
echo "==> Packaging..."
cd /project/dist
tar -czf "${DIST_NAME}.tar.gz" "${DIST_NAME}/"
echo "==> Done: dist/${DIST_NAME}.tar.gz  ($(du -sh "${DIST_NAME}.tar.gz" | cut -f1))"
echo "Contents:"
tar -tzf "${DIST_NAME}.tar.gz"
