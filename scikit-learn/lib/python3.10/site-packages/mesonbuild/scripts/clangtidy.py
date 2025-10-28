# SPDX-License-Identifier: Apache-2.0
# Copyright 2019 The Meson development team

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path
import tempfile
import os
import shutil
import sys

from .run_tool import run_with_buffered_output, run_clang_tool_on_sources
from ..environment import detect_clangtidy, detect_clangapply
import typing as T

async def run_clang_tidy(fname: Path, tidyexe: list, builddir: Path, fixesdir: T.Optional[Path]) -> int:
    args = []
    if fixesdir is not None:
        handle, name = tempfile.mkstemp(prefix=fname.name + '.', suffix='.yaml', dir=fixesdir)
        os.close(handle)
        args.extend(['-export-fixes', name])
    return await run_with_buffered_output(tidyexe + args + ['-quiet', '-p', str(builddir), str(fname)])

def run(args: T.List[str]) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument('--fix', action='store_true')
    parser.add_argument('--color', default='always')
    parser.add_argument('sourcedir')
    parser.add_argument('builddir')
    options = parser.parse_args(args)

    srcdir = Path(options.sourcedir)
    builddir = Path(options.builddir)

    tidyexe = detect_clangtidy()
    if not tidyexe:
        print(f'Could not execute clang-tidy "{" ".join(tidyexe)}"')
        return 1

    if options.color == 'always' or options.color == 'auto' and sys.stdout.isatty():
        tidyexe += ['--use-color']

    fixesdir: T.Optional[Path] = None
    if options.fix:
        applyexe = detect_clangapply()
        if not applyexe:
            print(f'Could not execute clang-apply-replacements "{" ".join(applyexe)}"')
            return 1

        fixesdir = builddir / 'meson-private' / 'clang-tidy-fix'
        if fixesdir.is_dir():
            shutil.rmtree(fixesdir)
        elif fixesdir.exists():
            fixesdir.unlink()
        fixesdir.mkdir(parents=True)

    tidyret = run_clang_tool_on_sources('clang-tidy', srcdir, builddir, run_clang_tidy, tidyexe, builddir, fixesdir)
    if fixesdir is not None:
        print('Applying fix-its...')
        applyret = subprocess.run(applyexe + ['-format', '-style=file', '-ignore-insert-conflict', fixesdir]).returncode

    if tidyret != 0:
        print('Errors encountered while running clang-tidy', file=sys.stderr)
        return tidyret
    if fixesdir is not None and applyret != 0:
        print('Errors encountered while running clang-apply-replacements', file=sys.stderr)
        return applyret
    return 0
