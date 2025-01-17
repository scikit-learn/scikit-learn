# SPDX-License-Identifier: Apache-2.0
# Copyright 2024 The Meson development team

from __future__ import annotations

import sys, os, subprocess, shutil
import pathlib
import typing as T

if T.TYPE_CHECKING:
    import argparse

from ..mesonlib import get_meson_command

# Note: when adding arguments, please also add them to the completion
# scripts in $MESONSRC/data/shell-completions/
def add_arguments(parser: 'argparse.ArgumentParser') -> None:
    parser.add_argument('--intermediaries',
                        default=False,
                        action='store_true',
                        help='Check intermediate files.')
    parser.add_argument('mesonargs', nargs='*',
                        help='Arguments to pass to "meson setup".')

IGNORE_PATTERNS = ('.ninja_log',
                   '.ninja_deps',
                   'meson-private',
                   'meson-logs',
                   'meson-info',
                   )

INTERMEDIATE_EXTENSIONS = ('.gch',
                           '.pch',
                           '.o',
                           '.obj',
                           '.class',
                           )

class ReproTester:
    def __init__(self, options: T.Any):
        self.args = options.mesonargs
        self.meson = get_meson_command()[:]
        self.builddir = pathlib.Path('buildrepro')
        self.storagedir = pathlib.Path('buildrepro.1st')
        self.issues: T.List[str] = []
        self.check_intermediaries = options.intermediaries

    def run(self) -> int:
        if not os.path.isfile('meson.build'):
            sys.exit('This command needs to be run at your project source root.')
        self.disable_ccache()
        self.cleanup()
        self.build()
        self.check_output()
        self.print_results()
        if not self.issues:
            self.cleanup()
        return len(self.issues)

    def disable_ccache(self) -> None:
        os.environ['CCACHE_DISABLE'] = '1'

    def cleanup(self) -> None:
        if self.builddir.exists():
            shutil.rmtree(self.builddir)
        if self.storagedir.exists():
            shutil.rmtree(self.storagedir)

    def build(self) -> None:
        setup_command: T.Sequence[str] = self.meson + ['setup', str(self.builddir)] + self.args
        build_command: T.Sequence[str] = self.meson + ['compile', '-C', str(self.builddir)]
        subprocess.check_call(setup_command)
        subprocess.check_call(build_command)
        self.builddir.rename(self.storagedir)
        subprocess.check_call(setup_command)
        subprocess.check_call(build_command)

    def ignore_file(self, fstr: str) -> bool:
        for p in IGNORE_PATTERNS:
            if p in fstr:
                return True
        if not self.check_intermediaries:
            if fstr.endswith(INTERMEDIATE_EXTENSIONS):
                return True
        return False

    def check_contents(self, fromdir: str, todir: str, check_contents: bool) -> None:
        import filecmp
        frompath = fromdir + '/'
        topath = todir + '/'
        for fromfile in pathlib.Path(fromdir).glob('**/*'):
            if not fromfile.is_file():
                continue
            fstr = fromfile.as_posix()
            if self.ignore_file(fstr):
                continue
            assert fstr.startswith(frompath)
            tofile = pathlib.Path(fstr.replace(frompath, topath, 1))
            if not tofile.exists():
                self.issues.append(f'Missing file: {tofile}')
            elif check_contents:
                if not filecmp.cmp(fromfile, tofile, shallow=False):
                    self.issues.append(f'File contents differ: {fromfile}')

    def print_results(self) -> None:
        if self.issues:
            print('Build differences detected')
            for i in self.issues:
                print(i)
        else:
            print('No differences detected.')

    def check_output(self) -> None:
        self.check_contents('buildrepro', 'buildrepro.1st', True)
        self.check_contents('buildrepro.1st', 'buildrepro', False)

def run(options: T.Any) -> None:
    rt = ReproTester(options)
    try:
        sys.exit(rt.run())
    except Exception as e:
        print(e)
        sys.exit(1)
