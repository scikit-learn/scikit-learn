# SPDX-License-Identifier: Apache-2.0
# Copyright 2018 The Meson development team

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from .run_tool import run_clang_tool, run_with_buffered_output
from ..tooldetect import detect_clangformat
from ..mesonlib import version_compare
from ..programs import ExternalProgram
import typing as T

async def run_clang_format(fname: Path, exelist: T.List[str], options: argparse.Namespace, cformat_ver: T.Optional[str]) -> int:
    clangformat_10 = False
    if options.check and cformat_ver:
        if version_compare(cformat_ver, '>=10'):
            clangformat_10 = True
            exelist = exelist + ['--dry-run', '--Werror']
            # The option is not documented but it exists in version 10
            if options.color == 'always' or options.color == 'auto' and sys.stdout.isatty():
                exelist += ['--color=1']
        else:
            original = fname.read_bytes()
    before = fname.stat().st_mtime
    ret = await run_with_buffered_output(exelist + ['-style=file', '-i', str(fname)])
    after = fname.stat().st_mtime
    if before != after:
        print('File reformatted: ', fname)
        if options.check and not clangformat_10:
            # Restore the original if only checking.
            fname.write_bytes(original)
            return 1
    return ret

def run(args: T.List[str]) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--color', default='always')
    parser.add_argument('sourcedir')
    parser.add_argument('builddir')
    options = parser.parse_args(args)

    srcdir = Path(options.sourcedir)
    builddir = Path(options.builddir)

    exelist = detect_clangformat()
    if not exelist:
        print('Could not execute clang-format "%s"' % ' '.join(exelist))
        return 1

    if options.check:
        cformat_ver = ExternalProgram('clang-format', exelist, silent=True).get_version()
    else:
        cformat_ver = None

    return run_clang_tool('clang-format', srcdir, builddir, run_clang_format, exelist, options, cformat_ver)
