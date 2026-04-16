# SPDX-License-Identifier: Apache-2.0
# Copyright 2013-2019 The Meson development team

from __future__ import annotations

import typing as T
import os
import re

from ..envconfig import detect_cpu_family
from ..mesonlib import Popen_safe
from .base import DependencyCandidate, DependencyException, DependencyMethods, detect_compiler, SystemDependency
from .configtool import ConfigToolDependency
from .detect import packages
from .factory import factory_methods
from .pkgconfig import PkgConfigDependency

if T.TYPE_CHECKING:
    from .factory import DependencyGenerator
    from ..environment import Environment
    from .base import DependencyObjectKWs


@factory_methods({DependencyMethods.PKGCONFIG, DependencyMethods.CONFIG_TOOL, DependencyMethods.SYSTEM})
def mpi_factory(env: 'Environment',
                kwargs: DependencyObjectKWs,
                methods: T.List[DependencyMethods]) -> T.List['DependencyGenerator']:
    language = kwargs.get('language') or 'c'
    if language not in {'c', 'cpp', 'fortran'}:
        # OpenMPI doesn't work without any other languages
        return []

    for_machine = kwargs['native']

    candidates: T.List['DependencyGenerator'] = []
    compiler = detect_compiler('mpi', env, for_machine, language)
    if not compiler:
        return []
    compiler_is_intel = compiler.get_id().startswith('intel')

    if DependencyMethods.CONFIG_TOOL in methods and not env.machines[for_machine].is_windows():
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
            # The oneAPI compilers have different wrappers
            is_llvm_based = 'llvm' in compiler.id
            if env.machines[for_machine].is_windows():
                nwargs['returncode_value'] = 3

            if language == 'c':
                if is_llvm_based:
                    tool_names.append('mpiicx')
                else:
                    tool_names.append('mpiicc')
            elif language == 'cpp':
                if is_llvm_based:
                    tool_names.append('mpiicpx')
                else:
                    tool_names.append('mpiicpc')
            elif language == 'fortran':
                if is_llvm_based:
                    tool_names.append('mpiifx')
                else:
                    tool_names.append('mpiifort')

        # even with intel compilers, mpicc has to be considered
        if language == 'c':
            tool_names.append('mpicc')
        elif language == 'cpp':
            tool_names.extend(['mpic++', 'mpicxx', 'mpiCC'])
        elif language == 'fortran':
            tool_names.extend(['mpifort', 'mpif90', 'mpif77'])

        nwargs['tools'] = tool_names
        candidates.append(DependencyCandidate.from_dependency(
            tool_names[0], MPIConfigToolDependency, (env, nwargs)))

    if DependencyMethods.SYSTEM in methods and env.machines[for_machine].is_windows():
        candidates.append(DependencyCandidate.from_dependency(
            'msmpi', MSMPIDependency, (env, kwargs)))
        candidates.append(DependencyCandidate.from_dependency(
            'impi', IMPIDependency, (env, kwargs)))

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
        candidates.append(DependencyCandidate.from_dependency(
            pkg_name, PkgConfigDependency, (env, kwargs)))

    return candidates

packages['mpi'] = mpi_factory


class MPIConfigToolDependency(ConfigToolDependency):
    """Wrapper around mpicc, Intel's mpiicc and friends."""

    def __init__(self, name: str, env: 'Environment', kwargs: DependencyObjectKWs):
        super().__init__(name, env, kwargs)
        if not self.is_found:
            return

        for comp, link in [
            ('--showme:compile', '--showme:link'),  # for OpenMPI
            ('-show-compile-info', '-show-link-info'),  # for MPICH and Intel MPI
            ('-compile_info', '-link_info'),  # for older MPICH and Intel MPI
            ('-show', None),
        ]:
            try:
                # Set required=True to ensure that the next set of options is
                # tried when the current ones fail, even if the dependency is
                # not required
                c_args = self.get_config_value([comp], 'compile_args', required=True)
                l_args = self.get_config_value([link], 'link_args', required=True) if link is not None else c_args
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

    def __init__(self, name: str, env: 'Environment', kwargs: DependencyObjectKWs):
        super().__init__(name, env, kwargs)
        # MSMPI only supports the C API
        if self.language not in {'c', 'fortran', None}:
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


class IMPIDependency(SystemDependency):

    """Intel(R) MPI for Windows."""

    def __init__(self, name: str, env: Environment, kwargs: DependencyObjectKWs):
        super().__init__(name, env, kwargs)
        # only for windows
        if not self.env.machines[self.for_machine].is_windows():
            return
        # only for x86_64
        if self.env.machines[self.for_machine].cpu_family != 'x86_64':
            return

        rootdir = os.environ.get('I_MPI_ROOT')
        if rootdir is None:
            self.is_found = False
            return

        incdir = os.path.join(rootdir, 'include')
        libdir = os.path.join(rootdir, 'lib')

        debug = env.coredata.optstore.get_value_for('debug')
        assert isinstance(debug, bool)
        libdir_post = 'debug' if debug else 'release'
        for subdirs in (['mpi', libdir_post], [libdir_post]):
            libdir_buildtype = os.path.join(libdir, *subdirs)
            if os.path.isdir(libdir_buildtype):
                libdir = libdir_buildtype
                break

        found_header = os.path.isfile(os.path.join(incdir, 'mpi.h'))
        found_library = os.path.isfile(os.path.join(libdir, 'impi.lib'))
        if not found_header or not found_library:
            self.is_found = False
            return

        self.is_found = True
        self.compile_args = ['-I' + incdir]
        self.link_args = ['-l' + os.path.join(libdir, 'impi')]
        if self.language == 'cpp':
            # Some installations do not have the MPI C++ bindings library
            if not os.path.isfile(os.path.join(libdir, 'impicxx.lib')):
                self.is_found = False
                return
            self.link_args = ['-l' + os.path.join(libdir, 'impicxx')]
