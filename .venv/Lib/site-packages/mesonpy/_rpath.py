# SPDX-FileCopyrightText: 2023 The meson-python developers
#
# SPDX-License-Identifier: MIT

from __future__ import annotations

import os
import subprocess
import sys
import typing


if typing.TYPE_CHECKING:
    from typing import List

    from mesonpy._compat import Iterable, Path


if sys.platform == 'win32' or sys.platform == 'cygwin':

    def fix_rpath(filepath: Path, libs_relative_path: str) -> None:
        pass

elif sys.platform == 'darwin':

    def _get_rpath(filepath: Path) -> List[str]:
        rpath = []
        r = subprocess.run(['otool', '-l', os.fspath(filepath)], capture_output=True, text=True)
        rpath_tag = False
        for line in [x.split() for x in r.stdout.split('\n')]:
            if line == ['cmd', 'LC_RPATH']:
                rpath_tag = True
            elif len(line) >= 2 and line[0] == 'path' and rpath_tag:
                rpath.append(line[1])
                rpath_tag = False
        return rpath

    def _replace_rpath(filepath: Path, old: str, new: str) -> None:
        subprocess.run(['install_name_tool', '-rpath', old, new, os.fspath(filepath)], check=True)

    def fix_rpath(filepath: Path, libs_relative_path: str) -> None:
        for path in _get_rpath(filepath):
            if path.startswith('@loader_path/'):
                _replace_rpath(filepath, path, '@loader_path/' + libs_relative_path)

elif sys.platform == 'sunos5':

    def _get_rpath(filepath: Path) -> List[str]:
        rpath = []
        r = subprocess.run(['/usr/bin/elfedit', '-r', '-e', 'dyn:rpath', os.fspath(filepath)],
            capture_output=True, check=True, text=True)
        for line in [x.split() for x in r.stdout.split('\n')]:
            if len(line) >= 4 and line[1] in ['RPATH', 'RUNPATH']:
                for path in line[3].split(':'):
                    if path not in rpath:
                        rpath.append(path)
        return rpath

    def _set_rpath(filepath: Path, rpath: Iterable[str]) -> None:
        subprocess.run(['/usr/bin/elfedit', '-e', 'dyn:rpath ' + ':'.join(rpath), os.fspath(filepath)], check=True)

    def fix_rpath(filepath: Path, libs_relative_path: str) -> None:
        old_rpath = _get_rpath(filepath)
        new_rpath = []
        for path in old_rpath:
            if path.startswith('$ORIGIN/'):
                path = '$ORIGIN/' + libs_relative_path
            new_rpath.append(path)
        if new_rpath != old_rpath:
            _set_rpath(filepath, new_rpath)

else:
    # Assume that any other platform uses ELF binaries.

    def _get_rpath(filepath: Path) -> List[str]:
        r = subprocess.run(['patchelf', '--print-rpath', os.fspath(filepath)], capture_output=True, text=True)
        return r.stdout.strip().split(':')

    def _set_rpath(filepath: Path, rpath: Iterable[str]) -> None:
        subprocess.run(['patchelf','--set-rpath', ':'.join(rpath), os.fspath(filepath)], check=True)

    def fix_rpath(filepath: Path, libs_relative_path: str) -> None:
        old_rpath = _get_rpath(filepath)
        new_rpath = []
        for path in old_rpath:
            if path.startswith('$ORIGIN/'):
                path = '$ORIGIN/' + libs_relative_path
            new_rpath.append(path)
        if new_rpath != old_rpath:
            _set_rpath(filepath, new_rpath)
