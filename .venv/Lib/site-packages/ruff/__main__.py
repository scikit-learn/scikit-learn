from __future__ import annotations

import os
import sys

from ruff import find_ruff_bin


def _run() -> None:
    ruff = find_ruff_bin()

    if sys.platform == "win32":
        import subprocess

        # Avoid emitting a traceback on interrupt
        try:
            completed_process = subprocess.run([ruff, *sys.argv[1:]])
        except KeyboardInterrupt:
            sys.exit(2)

        sys.exit(completed_process.returncode)
    else:
        os.execvp(ruff, [ruff, *sys.argv[1:]])


if __name__ == "__main__":
    _run()
