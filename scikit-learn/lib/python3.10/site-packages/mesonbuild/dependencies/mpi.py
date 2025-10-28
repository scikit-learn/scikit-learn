# SPDX-License-Identifier: Apache-2.0
# Copyright 2013-2019 The Meson development team

from __future__ import annotations

import functools
import typing as T
import os
import re

from ..environment import detect_cpu_family
from ..mesonlib import Popen_safe
from .base import DependencyException, DependencyMethods, detect_compiler, SystemDependency
from .configtool import ConfigToolDependency
from .detect import packages
from .factory import factory_methods
from .pkgconfig import PkgConfigDependency

if T.TYPE_CHECKING:
    from .factory import DependencyGenerator
    from ..environment import Environment
    from ..mesonlib import MachineChoice


@factory_methods({DependencyMethods.PKGCONFIG, DependencyMethods.CONFIG_TOOL, DependencyMethods.SYSTEM})
def mpi_factory(env: 'Environment',
                for_machine: 'MachineChoice',
                kwargs: T.Dict[str, T.Any],
                methods: T.List[DependencyMethods]) -> T.List['DependencyGenerator']:
    language = kwargs.get('language')
    if language is None:
        language = 'c'
    if language not in {'c', 'cpp', 'fortran'}:
        # OpenMPI doesn't work without any other languages
        return []

    candidates: T.List['DependencyGenerator'] = []
    compiler = detect_compiler('mpi', env, for_machine, language)
    if not compiler:
        return []
    compiler_is_intel = compiler.get_id() in {'intel', 'intel-cl'}

    if DependencyMethods.CONFIG_TOOL in methods:
        nwargs = kwargs.copy()

        # We try the environment variables for the tools first, but then
        # fall back to the hardcoded names

        if language == 'c':
            env_vars = ['MPICC']
        elif language == 'cpp':
            env_vars = ['MPICXX']
        elif language == 'fortran':
            env_vars = ['MPIFC', 'MPIF90', 'MPIF77']

        tool_names = [os.environ.get(env_name) for env_name in env_vars]
        tool_names = [t for t in tool_names if t]  # remove empty environment variables

        if compiler_is_intel:
            if env.machines[for_machine].is_windows():
                nwargs['returncode_value'] = 3

            if language == 'c':
                tool_names.append('mpiicc')
            elif language == 'cpp':
                tool_names.append('mpiicpc')
            elif language == 'fortran':
                tool_names.append('mpiifort')

        # even with intel compilers, mpicc has to be considered
        if language == 'c':
            tool_names.append('mpicc')
        elif language == 'cpp':
            tool_names.extend(['mpic++', 'mpicxx', 'mpiCC'])
        elif language == 'fortran':
            tool_names.extend(['mpifort', 'mpif90', 'mpif77'])

        nwargs['tools'] = tool_names
        candidates.append(functools.partial(
            MPIConfigToolDependency, tool_names[0], env, nwargs, language=language))

    if DependencyMethods.SYSTEM in methods and env.machines[for_machine].is_windows():
        candidates.append(functools.partial(
            MSMPIDependency, 'msmpi', env, kwargs, language=language))

    # Only OpenMPI has pkg-config, and it doesn't work with the intel compilers
    # for MPI, environment variables and commands like mpicc should have priority
    if DependencyMethods.PKGCONFIG in methods and not compiler_is_intel:
        pkg_name = None
        if language == 'c':
            pkg_name = 'ompi-c'
        elif language == 'cpp':
            pkg_name = 'ompi-cxx'
        elif language == 'fortran':
            pkg_name = 'ompi-fort'
        candidates.append(functools.partial(
            PkgConfigDependency, pkg_name, env, kwargs, language=language))

    return candidates

packages['mpi'] = mpi_factory


class MPIConfigToolDependency(ConfigToolDependency):
    """Wrapper around mpicc, Intel's mpiicc and friends."""

    def __init__(self, name: str, env: 'Environment', kwargs: T.Dict[str, T.Any],
                 language: T.Optional[str] = None):
        super().__init__(name, env, kwargs, language=language)
        if not self.is_found:
            return

        # --showme for OpenMPI, -compile_info/-link_info for MPICH and IntelMPI
        for comp, link in [('--showme:compile', '--showme:link'), ('-compile_info', '-link_info'), ('-show', None)]:
            try:
                c_args = self.get_config_value([comp], 'compile_args')
                l_args = self.get_config_value([link], 'link_args') if link is not None else c_args
            except DependencyException:
                continue
            else:
                break
        else:
            self.is_found = False
            return

        self.compile_args = self._filter_compile_args(c_args)
        self.link_args = self._filter_link_args(l_args)

    def _filter_compile_args(self, args: T.List[str]) -> T.List[str]:
        """
        MPI wrappers return a bunch of garbage args.
        Drop -O2 and everything that is not needed.
        """
        result = []
        multi_args: T.Tuple[str, ...] = ('-I', )
        if self.language == 'fortran':
            fc = self.env.coredata.compilers[self.for_machine]['fortran']
            multi_args += fc.get_module_incdir_args()

        include_next = False
        for f in args:
            if f.startswith(('-D', '-f') + multi_args) or f == '-pthread' \
                    or (f.startswith('-W') and f != '-Wall' and not f.startswith('-Werror')):
                result.append(f)
                if f in multi_args:
                    # Path is a separate argument.
                    include_next = True
            elif include_next:
                include_next = False
                result.append(f)
        return result

    def _filter_link_args(self, args: T.List[str]) -> T.List[str]:
        """
        MPI wrappers return a bunch of garbage args.
        Drop -O2 and everything that is not needed.
        """
        result = []
        include_next = False
        for f in args:
            if self._is_link_arg(f):
                result.append(f)
                if f in {'-L', '-Xlinker'}:
                    include_next = True
            elif include_next:
                include_next = False
                result.append(f)
        return result

    def _is_link_arg(self, f: str) -> bool:
        if self.clib_compiler.id == 'intel-cl':
            return f == '/link' or f.startswith('/LIBPATH') or f.endswith('.lib')   # always .lib whether static or dynamic
        else:
            return (f.startswith(('-L', '-l', '-Xlinker')) or
                    f == '-pthread' or
                    (f.startswith('-W') and f != '-Wall' and not f.startswith('-Werror')))

    def _check_and_get_version(self, tool: T.List[str], returncode: int) -> T.Tuple[bool, T.Union[str, None]]:
        p, out = Popen_safe(tool + ['--showme:version'])[:2]
        valid = p.returncode == returncode
        if valid:
            # OpenMPI
            v = re.search(r'\d+.\d+.\d+', out)
            if v:
                version = v.group(0)
            else:
                version = None
            return valid, version

        # --version is not the same as -v
        p, out = Popen_safe(tool + ['-v'])[:2]
        valid = p.returncode == returncode
        first_line = out.split('\n', maxsplit=1)[0]

        # cases like "mpicc for MPICH version 4.2.2"
        v = re.search(r'\d+.\d+.\d+', first_line)
        if v:
            return valid, v.group(0)

        # cases like "mpigcc for Intel(R) MPI library 2021.13"
        v = re.search(r'\d+.\d+', first_line)
        if v:
            return valid, v.group(0)

        # cases like "mpiifort for the Intel(R) MPI Library 2019 Update 9 for Linux*"
        v = re.search(r'(\d{4}) Update (\d)', first_line)
        if v:
            return valid, f'{v.group(1)}.{v.group(2)}'

        return valid, None


class MSMPIDependency(SystemDependency):

    """The Microsoft MPI."""

    def __init__(self, name: str, env: 'Environment', kwargs: T.Dict[str, T.Any],
                 language: T.Optional[str] = None):
        super().__init__(name, env, kwargs, language=language)
        # MSMPI only supports the C API
        if language not in {'c', 'fortran', None}:
            self.is_found = False
            return
        # MSMPI is only for windows, obviously
        if not self.env.machines[self.for_machine].is_windows():
            return

        incdir = os.environ.get('MSMPI_INC')
        arch = detect_cpu_family(self.env.coredata.compilers.host)
        libdir = None
        if arch == 'x86':
            libdir = os.environ.get('MSMPI_LIB32')
            post = 'x86'
        elif arch == 'x86_64':
            libdir = os.environ.get('MSMPI_LIB64')
            post = 'x64'

        if libdir is None or incdir is None:
            self.is_found = False
            return

        self.is_found = True
        self.link_args = ['-l' + os.path.join(libdir, 'msmpi')]
        self.compile_args = ['-I' + incdir, '-I' + os.path.join(incdir, post)]
        if self.language == 'fortran':
            self.link_args.append('-l' + os.path.join(libdir, 'msmpifec'))
