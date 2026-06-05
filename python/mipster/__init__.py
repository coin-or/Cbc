"""mipster — pre-built MIPster MIP solver shared library.

Resolves the appropriate libmipster shared library for the running CPU.

Two variants ship per x86_64 and Linux aarch64 wheel:

* Linux x86_64 / macOS x86_64 / Windows x86_64:
    generic + AVX2 (Haswell 2013+, AMD Zen2 2019+)
* Linux aarch64:
    generic (ARMv8-A) + neon (ARMv8.2-A — Cortex-A55/A75 2017+, Graviton2)
* macOS arm64:
    single build, ARMv8.5-A (Apple M1+)

By default, the most-optimised variant supported by the running CPU is
selected.  Override with the ``MIPSTER_BUILD`` environment variable:

    export MIPSTER_BUILD=generic   # force the baseline build
    export MIPSTER_BUILD=avx2      # force AVX2 (x86_64 only; raises if absent)
    export MIPSTER_BUILD=neon      # force ARMv8.2-A (Linux aarch64 only)
    export MIPSTER_BUILD=haswell   # alias for avx2 on macOS x86_64
    export MIPSTER_BUILD=debug     # debug+ASan build, when present

Set ``MIPSTER_VERBOSE=1`` to print the selected variant on first use.
"""

import os
import platform
import subprocess
from typing import Optional

__version__ = "0.0.0.dev0"

__all__ = [
    "dist_dir",
    "lib_dir",
    "bin_path",
    "include_dir",
    "lib_path",
    "selected_variant",
]


# ── CPU feature detection ────────────────────────────────────────────────────

def _has_avx2() -> bool:
    """Return True if the running CPU supports AVX2."""
    if platform.machine().lower() not in ("x86_64", "amd64"):
        return False
    try:
        if platform.system() == "Linux":
            with open("/proc/cpuinfo") as f:
                for line in f:
                    if line.startswith("flags"):
                        return "avx2" in line.split(":", 1)[1].split()
        elif platform.system() == "Darwin":
            r = subprocess.run(
                ["sysctl", "-n", "hw.optional.avx2_0"],
                capture_output=True, text=True,
            )
            return r.stdout.strip() == "1"
        elif os.name == "nt":
            import ctypes
            # PF_AVX2_INSTRUCTIONS_AVAILABLE = 40
            return bool(ctypes.windll.kernel32.IsProcessorFeaturePresent(40))
    except Exception:
        pass
    return False


def _has_armv8_2() -> bool:
    """Return True on aarch64 CPUs implementing ARMv8.2-A or later.

    Detected via the ``asimdrdm`` (FEAT_RDM, mandatory in ARMv8.1) and
    ``lrcpc`` (FEAT_LRCPC, mandatory in ARMv8.3) HWCAP bits exposed by
    the Linux kernel in /proc/cpuinfo.  We require ``asimdhp``
    (FP16 — mandatory in ARMv8.2) which is the most reliable single
    indicator that the CPU is at least ARMv8.2-A.
    """
    if platform.machine().lower() not in ("aarch64", "arm64"):
        return False
    if platform.system() != "Linux":
        return False
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("Features"):
                    feats = line.split(":", 1)[1].split()
                    return "asimdhp" in feats
    except Exception:
        pass
    return False


# ── Path resolution ──────────────────────────────────────────────────────────

_PKG_DIR = os.path.abspath(os.path.dirname(__file__))


def _candidate_dirs() -> dict:
    """Return {variant: absolute_path} for every dist directory present."""
    candidates = {}
    for entry in os.listdir(_PKG_DIR):
        if entry.startswith("mipster_dist"):
            full = os.path.join(_PKG_DIR, entry)
            if not os.path.isdir(full):
                continue
            if entry == "mipster_dist":
                candidates["generic"] = full
            else:
                # mipster_dist_avx2, _neon, _haswell, _debug, ...
                variant = entry[len("mipster_dist_"):]
                candidates[variant] = full
    return candidates


def _auto_select(variants: dict) -> str:
    """Pick the best variant the current CPU supports.

    Preference order on x86_64: avx2 / haswell > generic
    Preference order on aarch64 (Linux): neon > generic
    macOS arm64 / fallback: generic
    """
    machine = platform.machine().lower()
    if machine in ("x86_64", "amd64") and _has_avx2():
        for v in ("avx2", "haswell"):
            if v in variants:
                return v
    if machine in ("aarch64", "arm64") and _has_armv8_2():
        if "neon" in variants:
            return "neon"
    return "generic"


_selected_variant: Optional[str] = None
_selected_dir: Optional[str] = None


def _resolve() -> str:
    """Return the absolute path of the selected mipster_dist directory."""
    global _selected_variant, _selected_dir
    if _selected_dir is not None:
        return _selected_dir

    variants = _candidate_dirs()
    if not variants:
        raise RuntimeError(
            f"No mipster_dist* directory found under {_PKG_DIR}. "
            "The wheel appears to be incomplete; please reinstall."
        )

    override = os.environ.get("MIPSTER_BUILD", "").strip().lower()
    if override:
        if override not in variants:
            raise RuntimeError(
                f"MIPSTER_BUILD={override!r} requested but not present in this "
                f"installation. Available variants: {sorted(variants)}."
            )
        chosen = override
    elif "generic" in variants:
        chosen = _auto_select(variants)
        if chosen not in variants:
            chosen = "generic"
    else:
        # No 'generic' (e.g. macOS arm64 ships a single un-suffixed dist
        # directory called mipster_dist; this branch is defensive).
        chosen = next(iter(variants))

    _selected_variant = chosen
    _selected_dir = variants[chosen]

    if os.environ.get("MIPSTER_VERBOSE", "").strip() == "1":
        print(
            f"[mipster] variant={chosen}  dir={_selected_dir}",
            flush=True,
        )

    return _selected_dir


# ── Public API ───────────────────────────────────────────────────────────────

def dist_dir() -> str:
    """Absolute path of the selected ``mipster_dist*`` directory."""
    return _resolve()


def selected_variant() -> str:
    """Variant name actually selected (``generic``, ``avx2``, ``neon``, …)."""
    _resolve()
    assert _selected_variant is not None
    return _selected_variant


def bin_path() -> str:
    """Absolute path of the bundled ``mipster`` (or ``mipster.exe``) binary."""
    exe = "mipster.exe" if os.name == "nt" else "mipster"
    return os.path.join(dist_dir(), "bin", exe)


def include_dir() -> str:
    """Absolute path of the bundled C API headers (``mipster/Cbc_C_Interface.h``)."""
    return os.path.join(dist_dir(), "include", "mipster")


def lib_dir() -> str:
    """Absolute path of the directory holding the runtime shared library.

    Note: on Windows the DLL ships in ``bin/`` (libtool/MinGW convention),
    so this returns ``bin/`` there.  On Linux/macOS it returns ``lib/``.
    """
    if os.name == "nt":
        return os.path.join(dist_dir(), "bin")
    return os.path.join(dist_dir(), "lib")


def lib_path() -> str:
    """Absolute path of the libmipster shared library file itself."""
    d = lib_dir()
    if os.name == "nt":
        # Windows: MinGW libtool produces lib<name>-<N>.dll where N is the
        # current ABI version. Pick the first match.
        import glob
        matches = sorted(glob.glob(os.path.join(d, "libmipster-*.dll")))
        if not matches:
            raise FileNotFoundError(
                f"libmipster-*.dll not found in {d}; wheel may be corrupted."
            )
        return matches[0]
    if platform.system() == "Darwin":
        return os.path.join(d, "libmipster.dylib")
    return os.path.join(d, "libmipster.so")
