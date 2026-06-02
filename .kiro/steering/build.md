---
inclusion: fileMatch
fileMatchPattern: "src/Makefile.am,**/configure.ac,configster"
---

# MIPster — Build System

## Environment Variables (set in shell profile)

```sh
export MIPSTER_PREFIX="$HOME/prog/cbc"           # install prefix (opt build)
export MIPSTER_CBC_2X="$HOME/dev/cbc-stable"     # stable Cbc 2.x reference build
export MIPSTER_CBC_MASTER="$HOME/dev/cbc-master" # upstream Cbc master reference build
```

Debug builds install to `$MIPSTER_PREFIX-dbg` (ASan) and `$MIPSTER_PREFIX-tsan` (TSan).

## Preferred Workflow: `configster`

```sh
./configster --opt --install          # optimised build
./configster --debug --sanitizer=asan --install   # ASan
./configster --debug --sanitizer=tsan --install   # TSan
./configster --debug --sanitizer=valgrind --install  # Valgrind-friendly
```

Key options: `--static`/`--shared`/`--both`, `--prefix=PATH`, `--jobs=N`, `--dry-run`.

## Incremental Rebuild

```sh
make -j$(nproc) && make install   # rebuild + reinstall (no reconfigure)
make V=1 -j$(nproc)               # verbose rebuild only
```

## Adding New Source Files

1. Edit `src/Makefile.am` — add `.cpp` to `lib*_la_SOURCES`, `.hpp` to `includecoin_HEADERS`
2. Regenerate: `autoreconf --install -I BuildTools`
3. Rebuild: `./configster --opt --install`

## Running Tests

```sh
cd test && make -j$(nproc) test
# Individual binaries: test/CInterfaceTest, test/CbcSolverLpTest
```

## Code Formatting

```sh
./format-all-sources.sh               # all sources
clang-format -i path/to/file.cpp      # single file
```

Style: 2-space indent, no tabs, `PointerAlignment: Right`, `BreakBeforeBraces: WebKit`.

## Cleaning

```sh
make clean        # remove build objects (keeps Makefiles)
make distclean    # full clean (removes Makefiles; requires reconfigure after)
```
