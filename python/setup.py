"""Build a binary wheel from a pre-built mipster tarball/zip.

The mipster CI already produces self-contained per-platform artifacts:

    dist/mipster-linux-x86_64.tar.gz      → manylinux_2_34_x86_64
    dist/mipster-linux-aarch64.tar.gz     → manylinux_2_34_aarch64
    dist/mipster-macos-x86_64.tar.gz      → macosx_10_13_x86_64
    dist/mipster-macos-arm64.tar.gz       → macosx_11_0_arm64
    dist/mipster-windows-x86_64.zip       → win_amd64

This setup.py extracts ONE such artifact (path supplied via the
``MIPSTER_TARBALL`` env var) into the ``mipster/`` package and produces
a platform-tagged wheel.  It does NOT rebuild from source — that work
is done by the existing per-platform build scripts.

Required env vars:
    MIPSTER_TARBALL   absolute path to the tarball (.tar.gz or .zip)
    MIPSTER_VERSION   version string for the wheel (e.g. 3.0.0)
    MIPSTER_PLAT_TAG  PEP 425 platform tag (e.g. manylinux_2_34_x86_64)

Optional:
    MIPSTER_VARIANT_LAYOUT  comma-separated mapping  src_subdir=variant
                            (default chosen automatically per platform)
"""

import os
import re
import shutil
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path

from setuptools import setup
from wheel.bdist_wheel import bdist_wheel as _bdist_wheel


HERE = Path(__file__).resolve().parent
PKG_DIR = HERE / "mipster"


# ── Wheel customisation ──────────────────────────────────────────────────────

class platform_bdist_wheel(_bdist_wheel):
    """Force py3-none-<plat> tag and mark the wheel as platform-specific."""

    def finalize_options(self):
        super().finalize_options()
        self.root_is_pure = False

    def get_tag(self):
        python, abi, plat = super().get_tag()
        plat_override = os.environ.get("MIPSTER_PLAT_TAG", "").strip()
        return ("py3", "none", plat_override or plat)


# ── Tarball extraction & staging ─────────────────────────────────────────────

def _read_version() -> str:
    v = os.environ.get("MIPSTER_VERSION", "").strip()
    if not v:
        sys.exit("MIPSTER_VERSION not set (e.g. MIPSTER_VERSION=3.0.0).")
    if not re.match(r"^[0-9]+(\.[0-9]+)*([abc]?[0-9]+|\.dev[0-9]+|\.rc[0-9]+|\.post[0-9]+)?$", v):
        sys.exit(f"MIPSTER_VERSION={v!r} is not a valid PEP 440 version.")
    return v


def _extract(tarball: Path, dest: Path) -> Path:
    """Extract *tarball* (.tar.gz or .zip) into *dest* and return the top dir."""
    dest.mkdir(parents=True, exist_ok=True)
    name = tarball.name
    if name.endswith(".tar.gz") or name.endswith(".tgz"):
        with tarfile.open(tarball, "r:gz") as tf:
            members = tf.getmembers()
            top = members[0].name.split("/", 1)[0] if members else ""
            tf.extractall(dest)
    elif name.endswith(".zip"):
        with zipfile.ZipFile(tarball) as zf:
            members = zf.namelist()
            top = members[0].split("/", 1)[0] if members else ""
            zf.extractall(dest)
    else:
        sys.exit(f"Unrecognised archive format: {tarball}")
    return dest / top


def _stage_layout(extracted: Path, plat_tag: str) -> None:
    """Copy ``bin/``, ``lib/<variant>/``, ``include/`` from *extracted* into
    ``mipster/mipster_dist*/`` directories that match the runtime dispatcher's
    expected layout (see ``mipster/__init__.py``).

    Layout produced inside the wheel:

        Linux x86_64:
            mipster/mipster_dist/lib/libmipster.so*       (from extracted lib/generic/)
            mipster/mipster_dist/bin/mipster              (from bin/mipster-generic)
            mipster/mipster_dist/include/mipster/...
            mipster/mipster_dist_avx2/lib/libmipster.so*  (from extracted lib/avx2/)
            mipster/mipster_dist_avx2/bin/mipster         (from bin/mipster-avx2)

        Linux aarch64:   like x86_64 but optimized variant is 'neon'
        macOS x86_64:    optimized variant is 'haswell'
        macOS arm64:     single 'mipster_dist' (extracted lib/libmipster.dylib)
        Windows x86_64:  variant DLLs are libmipster-N-generic.dll / -avx2.dll in bin/
    """
    # Wipe any prior staging.
    for entry in PKG_DIR.iterdir():
        if entry.name.startswith("mipster_dist"):
            shutil.rmtree(entry)

    is_linux_x86  = plat_tag.startswith("manylinux") and plat_tag.endswith("_x86_64")
    is_linux_arm  = plat_tag.startswith("manylinux") and plat_tag.endswith("_aarch64")
    is_macos_x86  = plat_tag.startswith("macosx") and plat_tag.endswith("_x86_64")
    is_macos_arm  = plat_tag.startswith("macosx") and plat_tag.endswith("_arm64")
    is_windows    = plat_tag.startswith("win")

    if is_linux_x86 or is_linux_arm:
        opt = "avx2" if is_linux_x86 else "neon"
        _stage_linux(extracted, "generic", "generic")
        _stage_linux(extracted, opt, opt)
    elif is_macos_x86:
        _stage_macos_x86(extracted, "generic", "generic")
        _stage_macos_x86(extracted, "haswell", "haswell")
    elif is_macos_arm:
        _stage_macos_arm(extracted)
    elif is_windows:
        _stage_windows(extracted, "generic", "generic")
        _stage_windows(extracted, "avx2", "avx2")
    else:
        sys.exit(f"Unsupported MIPSTER_PLAT_TAG={plat_tag!r}")

    # Headers go alongside the generic build (any variant works — the C API
    # is identical across variants).  python-mip etc. read them from
    # mipster.include_dir() which always points at the selected dist dir.
    src_inc = extracted / "include" / "mipster"
    if src_inc.is_dir():
        for tgt in PKG_DIR.glob("mipster_dist*"):
            tgt_inc = tgt / "include" / "mipster"
            tgt_inc.parent.mkdir(parents=True, exist_ok=True)
            if not tgt_inc.exists():
                shutil.copytree(src_inc, tgt_inc)


def _stage_linux(extracted: Path, src_variant: str, wheel_variant: str) -> None:
    """Linux: copy entire lib/<src_variant>/ + bin/mipster-<src_variant>.

    The CI tarball already bundles libgfortran.so.5, libquadmath.so.0,
    libgcc_s.so.1, libstdc++.so.6, and libbz2.so.1 alongside libmipster.so.0
    in lib/<variant>/. RUNPATH=$ORIGIN on libmipster makes those resolve
    relative to the wheel's lib dir at runtime, so the wheel must ship them
    all — not just libmipster.so*.
    """
    dist = PKG_DIR / ("mipster_dist" if wheel_variant == "generic"
                      else f"mipster_dist_{wheel_variant}")
    (dist / "lib").mkdir(parents=True, exist_ok=True)
    (dist / "bin").mkdir(parents=True, exist_ok=True)
    src_lib = extracted / "lib" / src_variant
    if not src_lib.is_dir():
        sys.exit(f"Expected {src_lib} in tarball but it is missing.")
    for so in src_lib.iterdir():
        # Copy every .so / .so.N (real files and symlinks) — preserves the
        # bundled non-mipster runtime libs that libmipster.so depends on.
        if so.is_symlink() or ".so" in so.name:
            shutil.copy2(so, dist / "lib" / so.name, follow_symlinks=False)
    src_bin = extracted / "bin" / f"mipster-{src_variant}"
    if src_bin.is_file():
        out_bin = dist / "bin" / "mipster"
        shutil.copy2(src_bin, out_bin)
        os.chmod(out_bin, 0o755)
        # The tarball binary's RUNPATH is "$ORIGIN/../lib/<variant>" (because
        # the tarball ships shared bin/ and per-variant lib/<variant>/). In
        # the wheel each variant has its own mipster_dist_<variant>/lib/, so
        # the binary needs RUNPATH=$ORIGIN/../lib.
        _patch_runpath(out_bin, "$ORIGIN/../lib")


def _patch_runpath(elf: Path, new_runpath: str) -> None:
    """Set RUNPATH on an ELF binary. patchelf is required (apt: patchelf)."""
    patchelf = shutil.which("patchelf")
    if not patchelf:
        sys.exit(
            "patchelf not found on PATH but is required for Linux wheel staging "
            "(install with `apt-get install patchelf` or `pip install patchelf`)."
        )
    subprocess.run(
        [patchelf, "--set-rpath", new_runpath, str(elf)],
        check=True,
    )


def _stage_macos_x86(extracted: Path, src_variant: str, wheel_variant: str) -> None:
    dist = PKG_DIR / ("mipster_dist" if wheel_variant == "generic"
                      else f"mipster_dist_{wheel_variant}")
    (dist / "lib").mkdir(parents=True, exist_ok=True)
    (dist / "bin").mkdir(parents=True, exist_ok=True)
    src_lib = extracted / "lib" / src_variant
    if not src_lib.is_dir():
        sys.exit(f"Expected {src_lib} in tarball but it is missing.")
    for dy in src_lib.iterdir():
        if dy.is_symlink() or dy.name.endswith(".dylib"):
            shutil.copy2(dy, dist / "lib" / dy.name, follow_symlinks=False)
    src_bin = extracted / "bin" / f"mipster-{src_variant}"
    if src_bin.is_file():
        shutil.copy2(src_bin, dist / "bin" / "mipster")
        os.chmod(dist / "bin" / "mipster", 0o755)


def _stage_macos_arm(extracted: Path) -> None:
    """Single-variant macOS arm64 build (no dispatch required)."""
    dist = PKG_DIR / "mipster_dist"
    (dist / "lib").mkdir(parents=True, exist_ok=True)
    (dist / "bin").mkdir(parents=True, exist_ok=True)
    src_lib = extracted / "lib"
    for dy in src_lib.iterdir():
        if dy.is_symlink() or dy.name.endswith(".dylib"):
            shutil.copy2(dy, dist / "lib" / dy.name, follow_symlinks=False)
    src_bin = extracted / "bin" / "mipster"
    if src_bin.is_file():
        shutil.copy2(src_bin, dist / "bin" / "mipster")
        os.chmod(dist / "bin" / "mipster", 0o755)


def _stage_windows(extracted: Path, src_variant: str, wheel_variant: str) -> None:
    """Windows: copy bin/<src_variant>/ wholesale into mipster_dist*/bin/.

    The CI tarball lays each variant out as bin/<variant>/{mipster.exe,
    libmipster-0.dll, libgcc_s_seh-1.dll, libstdc++-6.dll,
    libwinpthread-1.dll, zlib1.dll}.  All of those must travel together
    into the wheel — the variant EXE imports libmipster-0.dll, which in
    turn imports the four MinGW runtime DLLs, and Windows' DLL search
    looks in the directory of the loading EXE first.
    """
    dist = PKG_DIR / ("mipster_dist" if wheel_variant == "generic"
                      else f"mipster_dist_{wheel_variant}")
    (dist / "bin").mkdir(parents=True, exist_ok=True)
    src_variant_dir = extracted / "bin" / src_variant
    if not src_variant_dir.is_dir():
        sys.exit(
            f"Expected {src_variant_dir}/ in Windows tarball but it is missing."
        )
    for f in src_variant_dir.iterdir():
        if f.is_file() and (f.suffix.lower() in (".dll", ".exe")):
            shutil.copy2(f, dist / "bin" / f.name)
    # The launcher itself lives in bin/mipster.exe (one level up). It picks
    # the variant subdir at runtime — but inside the wheel each variant has
    # its own mipster_dist*/bin/mipster.exe (already copied above), so the
    # Python dispatcher in mipster/__init__.py picks the dir directly and
    # the launcher is unused by the wheel. We still copy it for parity with
    # the tarball / for users who poke around in the wheel.
    launcher = extracted / "bin" / "mipster.exe"
    if launcher.is_file() and src_variant == "generic":
        # Only need one copy; goes into the generic dist.
        pass  # the per-variant mipster.exe (above) already serves this role


# ── Driver ───────────────────────────────────────────────────────────────────

def _main() -> None:
    version = _read_version()
    plat_tag = os.environ.get("MIPSTER_PLAT_TAG", "").strip()
    tarball = os.environ.get("MIPSTER_TARBALL", "").strip()

    is_packaging = any(
        cmd in sys.argv for cmd in ("bdist_wheel", "build", "install")
    )

    if is_packaging:
        if not plat_tag:
            sys.exit("MIPSTER_PLAT_TAG not set (e.g. MIPSTER_PLAT_TAG=manylinux_2_34_x86_64).")
        if not tarball:
            sys.exit("MIPSTER_TARBALL not set (path to mipster-<plat>.tar.gz / .zip).")
        tarball_path = Path(tarball)
        if not tarball_path.is_file():
            sys.exit(f"MIPSTER_TARBALL={tarball_path} does not exist.")

        # Wipe setuptools' build/ — package_data globs (mipster_dist*/**/*) would
        # otherwise pick up stale staged files from a previous platform's wheel
        # build (e.g. Linux .so files leaking into a Windows wheel).
        build_dir = HERE / "build"
        if build_dir.exists():
            shutil.rmtree(build_dir, ignore_errors=True)

        extract_root = HERE / "_extract"
        if extract_root.exists():
            shutil.rmtree(extract_root)
        extracted = _extract(tarball_path, extract_root)
        _stage_layout(extracted, plat_tag)

    setup(
        version=version,
        cmdclass={"bdist_wheel": platform_bdist_wheel},
        packages=["mipster"],
        zip_safe=False,
        include_package_data=True,
        package_data={
            "mipster": [
                "mipster_dist/**/*",
                "mipster_dist_*/**/*",
            ],
        },
    )


_main()
