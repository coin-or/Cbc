# mipster (Python wheel)

This directory packages MIPster as a binary Python wheel. It does not rebuild
the C++ sources — it repackages the platform-specific tarballs produced by the
existing CI scripts under `.github/scripts/`.

## Layout inside the wheel

```
mipster/
├── __init__.py              # CPU-dispatch logic, lib_dir() / bin_path() / ...
├── __main__.py              # `python -m mipster ...` CLI launcher
├── mipster_dist/            # baseline (always present)
│   ├── bin/mipster          # CLI binary
│   ├── lib/libmipster.so*   # (or .dylib / libmipster-N.dll on Windows)
│   └── include/mipster/...  # public C API headers
└── mipster_dist_avx2/       # x86_64 only — optimised (Haswell 2013+)
    ├── bin/mipster
    └── lib/libmipster.so*
```

On Linux aarch64 the optimised dir is `mipster_dist_neon` (ARMv8.2-A).
On macOS x86_64 it is `mipster_dist_haswell`. macOS arm64 ships only the
baseline (M1+ ARMv8.5-A — no dispatch needed).

## Building a wheel locally

```sh
export MIPSTER_VERSION=3.0.0.dev0
export MIPSTER_PLAT_TAG=manylinux_2_34_x86_64
export MIPSTER_TARBALL=/path/to/mipster-linux-x86_64.tar.gz
python -m build --wheel
```

## CPU dispatch at runtime

```python
import mipster
print(mipster.selected_variant())   # 'avx2' | 'generic' | 'neon' | ...
print(mipster.lib_path())           # absolute path to libmipster.{so,dylib,dll}
```

Override with the `MIPSTER_BUILD` env var:

```sh
export MIPSTER_BUILD=generic   # force baseline
export MIPSTER_BUILD=avx2      # force AVX2 (errors out if not bundled)
```
