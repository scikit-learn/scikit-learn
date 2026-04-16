from __future__ import annotations

import os
import shutil
import subprocess

from . import destdir_join

import argparse
import typing as T

parser = argparse.ArgumentParser()
parser.add_argument('--install')
parser.add_argument('--extra-extension-path', action="append", default=[])
parser.add_argument('--name')
parser.add_argument('--builddir')
parser.add_argument('--project-version')
parser.add_argument('--docdir')


def run(argv: T.List[str]) -> int:
    options, args = parser.parse_known_args(argv)
    subenv = os.environ.copy()

    val = subenv.get('PYTHONPATH')
    paths = [val] if val else []
    subenv['PYTHONPATH'] = os.pathsep.join(paths + options.extra_extension_path)

    res = subprocess.call(args, cwd=options.builddir, env=subenv)
    if res != 0:
        return res

    if options.install:
        source_dir = os.path.join(options.builddir, options.install)
        destdir = os.environ.get('DESTDIR', '')
        installdir = destdir_join(destdir, options.docdir)

        shutil.rmtree(installdir, ignore_errors=True)
        shutil.copytree(source_dir, installdir)
    return 0
