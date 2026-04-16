# SPDX-License-Identifier: Apache-2.0
# Copyright 2015-2016 The Meson development team

from __future__ import annotations

'''CD into dir given as first argument and execute
the command given in the rest of the arguments.'''

import os, subprocess, sys
import typing as T

def run(args: T.List[str]) -> int:
    dirname = args[0]
    command = args[1:]

    os.chdir(dirname)
    return subprocess.call(command)

if __name__ == '__main__':
    sys.exit(run(sys.argv[1:]))
