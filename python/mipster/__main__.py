"""``python -m mipster`` — invoke the bundled MIPster CLI binary.

Forwards all arguments verbatim and propagates the binary's exit code.
"""

import os
import sys

from . import bin_path


def main() -> int:
    exe = bin_path()
    if os.name == "nt":
        # On Windows os.execv replaces the current process but the exit
        # code propagation is fragile; use subprocess instead.
        import subprocess
        return subprocess.call([exe, *sys.argv[1:]])
    os.execv(exe, [exe, *sys.argv[1:]])


if __name__ == "__main__":
    sys.exit(main())
