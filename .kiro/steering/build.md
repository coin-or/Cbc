---
inclusion: fileMatch
fileMatchPattern: "Cbc/src/Makefile.am,**/configure.ac,Cbc/scripts/configster,cbc_build.sh,cbc_clean.sh"
---

# MIPster — Build System

## Environment Variables (set in shell profile)

```sh
export MIPSTER_PREFIX="$HOME/prog/cbc"           # install prefix (opt build)
export MIPSTER_CBC_2X="$HOME/dev/cbc-stable"     # stable Cbc 2.x reference build
export MIPSTER_CBC_MASTER="$HOME/dev/cbc-master" # upstream Cbc master reference build
```

Debug builds install to `$MIPSTER_PREFIX-dbg` (ASan) and `$MIPSTER_PREFIX-tsan` (TSan).

## Preferred Workflow: `Cbc/scripts/configster`

```sh
./Cbc/scripts/configster --opt --install          # optimised build
./Cbc/scripts/configster --debug --sanitizer=asan --install   # ASan
./Cbc/scripts/configster --debug --sanitizer=tsan --install   # TSan
./Cbc/scripts/configster --debug --sanitizer=valgrind --install  # Valgrind-friendly
```

Key options: `--static`/`--shared`/`--both`, `--prefix=PATH`, `--jobs=N`, `--dry-run`.

## Incremental Rebuild

```sh
./cbc_build.sh --install   # rebuild + reinstall (safest, no reconfigure)
./cbc_build.sh             # rebuild only
cd Cbc && make V=1 -j$(nproc) && make install   # direct
```

## Adding New Source Files

1. Edit `Cbc/src/Makefile.am` — add `.cpp` to `lib*_la_SOURCES`, `.hpp` to `includecoin_HEADERS`
2. Regenerate: `cd Cbc && autoreconf --install -I BuildTools`
3. Rebuild: `./Cbc/scripts/configster --opt --install`

## Running Tests

```sh
cd Cbc/test && make -j$(nproc) test
# Individual binaries: Cbc/test/CInterfaceTest, Cbc/test/cbc_unittest
```

## Code Formatting

```sh
cd Cbc && ./format-all-sources.sh     # all sources
clang-format -i path/to/file.cpp      # single file
```

Style: 2-space indent, no tabs, `PointerAlignment: Right`, `BreakBeforeBraces: WebKit`.

## Cleaning

```sh
./cbc_clean.sh                  # clean + distclean
./cbc_clean.sh --distclean-only # distclean only (removes Makefiles)
```
