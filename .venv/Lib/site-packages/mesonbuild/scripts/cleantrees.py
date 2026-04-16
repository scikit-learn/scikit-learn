# SPDX-License-Identifier: Apache-2.0
# Copyright 2016 The Meson development team

from __future__ import annotations

import os
import sys
import shutil
import pickle
import typing as T

from ..utils.platform import path_has_root

def rmtrees(build_dir: str, trees: T.List[str]) -> None:
    for t in trees:
        # Never delete trees outside of the builddir
        if path_has_root(t):
            print(f'Cannot delete dir with absolute path {t!r}')
            continue
        bt = os.path.join(build_dir, t)
        # Skip if it doesn't exist, or if it is not a directory
        if os.path.isdir(bt):
            shutil.rmtree(bt, ignore_errors=True)

def run(args: T.List[str]) -> int:
    if len(args) != 1:
        print('Cleaner script for Meson. Do not run on your own please.')
        print('cleantrees.py <data-file>')
        return 1
    with open(args[0], 'rb') as f:
        data = pickle.load(f)
    rmtrees(data.build_dir, data.trees)
    # Never fail cleaning
    return 0

if __name__ == '__main__':
    run(sys.argv[1:])
