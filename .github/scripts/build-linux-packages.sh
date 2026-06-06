#!/usr/bin/env bash
# build-linux-packages.sh — Build DEB, RPM, and Arch pkg.tar.zst from a
# pre-built MIPster Linux tarball.
#
# Usage:
#   build-linux-packages.sh <tarball> [arch]
#
#   <tarball>  Path to a mipster-linux-<arch>.tar.gz from build-manylinux.sh
#   [arch]     Package architecture (default: x86_64)
#
# Produces in the current working directory:
#   mipster_<ver>_<arch>.deb
#   mipster-<ver>-1.<arch>.rpm
#   mipster-<ver>-1-<arch>.pkg.tar.zst  (only if makepkg is available)
#
# Requirements:
#   DEB/RPM: fpm  (gem install fpm)
#   Arch:    makepkg + base-devel (Arch Linux only)
set -euo pipefail

TARBALL="${1:?Usage: $0 <tarball> [arch]}"
PKG_ARCH="${2:-x86_64}"

echo "==> MIPster Linux package builder"
echo "    tarball : ${TARBALL}"
echo "    arch    : ${PKG_ARCH}"

# ── Extract tarball ───────────────────────────────────────────────────────────
WORK_DIR="$(mktemp -d)"
trap 'rm -rf "${WORK_DIR}"' EXIT

echo ""
echo "==> Extracting tarball..."
tar -xzf "${TARBALL}" -C "${WORK_DIR}"
DIST_DIR="$(find "${WORK_DIR}" -maxdepth 1 -mindepth 1 -type d | head -1)"
echo "    dist dir: ${DIST_DIR}"

# ── Detect version ─────────────────────────────────────────────────────────────
VERSION="$("${DIST_DIR}/bin/mipster-generic" --version 2>&1 \
  | grep -Eo '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || echo '')"
if [ -z "${VERSION}" ]; then
  # Dev build: use latest git tag + short commit hash for traceability
  LATEST_TAG="$(git describe --tags --abbrev=0 2>/dev/null | sed 's/^v//' || echo '0.0.0')"
  GIT_HASH="$(git rev-parse --short HEAD 2>/dev/null || echo 'unknown')"
  VERSION="${LATEST_TAG}+devel.${GIT_HASH}"
fi
echo "    version : ${VERSION}"

# ── Build staging tree (FHS layout) ──────────────────────────────────────────
STAGE="${WORK_DIR}/stage"
mkdir -p \
  "${STAGE}/usr/bin" \
  "${STAGE}/usr/lib/mipster/generic" \
  "${STAGE}/usr/lib/mipster/avx2" \
  "${STAGE}/usr/lib/mipster/neon" \
  "${STAGE}/usr/include/mipster" \
  "${STAGE}/usr/share/doc/mipster" \
  "${STAGE}/usr/share/man/man1" \
  "${STAGE}/usr/share/applications"

# Binaries
install -m 755 "${DIST_DIR}/bin/mipster"         "${STAGE}/usr/bin/"
install -m 755 "${DIST_DIR}/bin/mipster-generic" "${STAGE}/usr/bin/"
[ -f "${DIST_DIR}/bin/mipster-avx2" ]  && install -m 755 "${DIST_DIR}/bin/mipster-avx2"  "${STAGE}/usr/bin/"
[ -f "${DIST_DIR}/bin/mipster-neon" ]  && install -m 755 "${DIST_DIR}/bin/mipster-neon"  "${STAGE}/usr/bin/"

# Shared libraries (if present in dist)
[ -d "${DIST_DIR}/lib/generic" ] && cp -a "${DIST_DIR}/lib/generic/." "${STAGE}/usr/lib/mipster/generic/"
[ -d "${DIST_DIR}/lib/avx2" ]    && cp -a "${DIST_DIR}/lib/avx2/."    "${STAGE}/usr/lib/mipster/avx2/"
[ -d "${DIST_DIR}/lib/neon" ]    && cp -a "${DIST_DIR}/lib/neon/."    "${STAGE}/usr/lib/mipster/neon/"

# Headers
[ -d "${DIST_DIR}/include/mipster" ] && cp -a "${DIST_DIR}/include/mipster/." "${STAGE}/usr/include/mipster/"

# Documentation
[ -d "${DIST_DIR}/share/doc/mipster" ] && \
  cp -a "${DIST_DIR}/share/doc/mipster/." "${STAGE}/usr/share/doc/mipster/"
[ -f "${DIST_DIR}/share/man/man1/mipster.1.gz" ] && \
  cp "${DIST_DIR}/share/man/man1/mipster.1.gz" "${STAGE}/usr/share/man/man1/"

# Desktop files
[ -f "${DIST_DIR}/share/applications/mipster.desktop" ] && \
  cp "${DIST_DIR}/share/applications/mipster.desktop"      "${STAGE}/usr/share/applications/"
[ -f "${DIST_DIR}/share/applications/mipster-docs.desktop" ] && \
  cp "${DIST_DIR}/share/applications/mipster-docs.desktop" "${STAGE}/usr/share/applications/"

echo ""
echo "==> Staging tree:"
find "${STAGE}" -type f | sed "s|${STAGE}||" | sort

# ── Map arch names ─────────────────────────────────────────────────────────────
DEB_ARCH="${PKG_ARCH}"
RPM_ARCH="${PKG_ARCH}"
case "${PKG_ARCH}" in
  x86_64)  DEB_ARCH="amd64"  ; RPM_ARCH="x86_64"  ;;
  aarch64) DEB_ARCH="arm64"  ; RPM_ARCH="aarch64"  ;;
esac

OUT_DIR="${PWD}"

# ── Build DEB ─────────────────────────────────────────────────────────────────
if command -v fpm &>/dev/null; then
  echo ""
  echo "==> Building DEB package..."
  fpm \
    --input-type dir \
    --output-type deb \
    --name mipster \
    --version "${VERSION}" \
    --architecture "${DEB_ARCH}" \
    --description "MIPster — Mixed-Integer Programming Solver" \
    --url "https://github.com/h-g-s/mipster" \
    --maintainer "h-g-s" \
    --license "EPL-2.0" \
    --category math \
    --deb-no-default-config-files \
    --package "${OUT_DIR}/mipster_${VERSION}_${DEB_ARCH}.deb" \
    --chdir "${STAGE}" \
    .
  echo "    DEB: mipster_${VERSION}_${DEB_ARCH}.deb  ($(du -sh "${OUT_DIR}/mipster_${VERSION}_${DEB_ARCH}.deb" | cut -f1))"

  echo ""
  echo "==> Building RPM package..."
  fpm \
    --input-type dir \
    --output-type rpm \
    --name mipster \
    --version "${VERSION}" \
    --architecture "${RPM_ARCH}" \
    --description "MIPster — Mixed-Integer Programming Solver" \
    --url "https://github.com/h-g-s/mipster" \
    --maintainer "h-g-s" \
    --license "EPL-2.0" \
    --category "Applications/Engineering" \
    --rpm-summary "Mixed-Integer Programming Solver" \
    --package "${OUT_DIR}/mipster-${VERSION}-1.${RPM_ARCH}.rpm" \
    --chdir "${STAGE}" \
    .
  echo "    RPM: mipster-${VERSION}-1.${RPM_ARCH}.rpm  ($(du -sh "${OUT_DIR}/mipster-${VERSION}-1.${RPM_ARCH}.rpm" | cut -f1))"
else
  echo "WARNING: fpm not found — skipping DEB/RPM build. Install with: gem install fpm"
fi

# ── Build Arch pkg.tar.zst ─────────────────────────────────────────────────────
if command -v makepkg &>/dev/null; then
  echo ""
  echo "==> Building Arch package (pkg.tar.zst)..."

  ARCH_BUILD="${WORK_DIR}/arch"
  mkdir -p "${ARCH_BUILD}"

  # Write PKGBUILD
  cat > "${ARCH_BUILD}/PKGBUILD" << PKGBUILD_EOF
# Maintainer: h-g-s
pkgname=mipster
pkgver=${VERSION}
pkgrel=1
pkgdesc="MIPster — Mixed-Integer Programming Solver"
arch=('${PKG_ARCH}')
url="https://github.com/h-g-s/mipster"
license=('EPL-2.0')
provides=('mipster')

package() {
  cp -a "${STAGE}/." "\${pkgdir}/"
}
PKGBUILD_EOF

  cd "${ARCH_BUILD}"
  PKGDEST="${ARCH_BUILD}" makepkg --noconfirm --nodeps 2>&1
  find "${ARCH_BUILD}" -maxdepth 1 -name '*.pkg.tar.zst' -exec mv {} "${OUT_DIR}/" \;
  echo "    Arch pkg written to ${OUT_DIR}"
  cd "${OUT_DIR}"
else
  echo "NOTE: makepkg not found — skipping Arch package build (requires Arch Linux)"
fi

echo ""
echo "==> Done. Packages in: ${OUT_DIR}"
ls -lh "${OUT_DIR}"/mipster*.deb "${OUT_DIR}"/mipster*.rpm "${OUT_DIR}"/mipster*.pkg.tar.zst 2>/dev/null || true
