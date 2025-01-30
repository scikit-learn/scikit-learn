# SPDX-License-Identifier: Apache-2.0
# Copyright 2019 The Meson development team

from __future__ import annotations

import os
import argparse
import multiprocessing
import subprocess
from pathlib import Path
import typing as T

from ..mesonlib import Popen_safe, split_args

class ExternalProject:
    def __init__(self, options: argparse.Namespace):
        self.name = options.name
        self.src_dir = options.srcdir
        self.build_dir = options.builddir
        self.install_dir = options.installdir
        self.log_dir = options.logdir
        self.verbose = options.verbose
        self.stampfile = options.stampfile
        self.depfile = options.depfile
        self.make = split_args(options.make)

    def write_depfile(self) -> None:
        with open(self.depfile, 'w', encoding='utf-8') as f:
            f.write(f'{self.stampfile}: \\\n')
            for dirpath, dirnames, filenames in os.walk(self.src_dir):
                dirnames[:] = [d for d in dirnames if not d.startswith('.')]
                for fname in filenames:
                    if fname.startswith('.'):
                        continue
                    path = Path(dirpath, fname)
                    f.write('  {} \\\n'.format(path.as_posix().replace(' ', '\\ ')))

    def write_stampfile(self) -> None:
        with open(self.stampfile, 'w', encoding='utf-8'):
            pass

    def supports_jobs_flag(self) -> bool:
        p, o, e = Popen_safe(self.make + ['--version'])
        if p.returncode == 0 and ('GNU Make' in o or 'waf' in o):
            return True
        return False

    def build(self) -> int:
        make_cmd = self.make.copy()
        if self.supports_jobs_flag():
            make_cmd.append(f'-j{multiprocessing.cpu_count()}')
        rc = self._run('build', make_cmd)
        if rc != 0:
            return rc

        install_cmd = self.make.copy()
        install_env = {}
        install_env['DESTDIR'] = self.install_dir
        install_cmd.append('install')
        rc = self._run('install', install_cmd, install_env)
        if rc != 0:
            return rc

        self.write_depfile()
        self.write_stampfile()

        return 0

    def _run(self, step: str, command: T.List[str], env: T.Optional[T.Dict[str, str]] = None) -> int:
        m = 'Running command ' + str(command) + ' in directory ' + str(self.build_dir) + '\n'
        log_filename = Path(self.log_dir, f'{self.name}-{step}.log')
        output = None
        if not self.verbose:
            output = open(log_filename, 'w', encoding='utf-8')
            output.write(m + '\n')
            output.flush()
        else:
            print(m)
        run_env = os.environ.copy()
        if env:
            run_env.update(env)
        p, o, e = Popen_safe(command, stderr=subprocess.STDOUT, stdout=output,
                             cwd=self.build_dir,
                             env=run_env)
        if p.returncode != 0:
            m = f'{step} step returned error code {p.returncode}.'
            if not self.verbose:
                m += '\nSee logs: ' + str(log_filename)
            print(m)
        return p.returncode

def run(args: T.List[str]) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument('--name')
    parser.add_argument('--srcdir')
    parser.add_argument('--builddir')
    parser.add_argument('--installdir')
    parser.add_argument('--logdir')
    parser.add_argument('--make')
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('stampfile')
    parser.add_argument('depfile')

    options = parser.parse_args(args)
    ep = ExternalProject(options)
    return ep.build()
