# SPDX-License-Identifier: Apache-2.0
# Copyright 2013-2019 The Meson development team

# This file contains the detection logic for miscellaneous external dependencies.
from __future__ import annotations

import functools
import os
import re
from pathlib import Path

from ..mesonlib import OrderedSet, join_args
from .base import DependencyException, DependencyMethods
from .configtool import ConfigToolDependency
from .detect import packages
from .pkgconfig import PkgConfigDependency, PkgConfigInterface
from .factory import factory_methods
import typing as T

if T.TYPE_CHECKING:
    from .factory import DependencyGenerator
    from ..environment import Environment
    from ..mesonlib import MachineChoice


class HDF5PkgConfigDependency(PkgConfigDependency):

    """Handle brokenness in the HDF5 pkg-config files."""

    def __init__(self, name: str, environment: 'Environment', kwargs: T.Dict[str, T.Any], language: T.Optional[str] = None) -> None:
        language = language or 'c'
        if language not in {'c', 'cpp', 'fortran'}:
            raise DependencyException(f'Language {language} is not supported with HDF5.')

        super().__init__(name, environment, kwargs, language)
        if not self.is_found:
            return

        # some broken pkgconfig don't actually list the full path to the needed includes
        newinc: T.List[str] = []
        for arg in self.compile_args:
            if arg.startswith('-I'):
                stem = 'static' if self.static else 'shared'
                if (Path(arg[2:]) / stem).is_dir():
                    newinc.append('-I' + str(Path(arg[2:]) / stem))
        self.compile_args += newinc

        link_args: T.List[str] = []
        for larg in self.get_link_args():
            lpath = Path(larg)
            # some pkg-config hdf5.pc (e.g. Ubuntu) don't include the commonly-used HL HDF5 libraries,
            # so let's add them if they exist
            # additionally, some pkgconfig HDF5 HL files are malformed so let's be sure to find HL anyway
            if lpath.is_file():
                hl = []
                if language == 'cpp':
                    hl += ['_hl_cpp', '_cpp']
                elif language == 'fortran':
                    hl += ['_hl_fortran', 'hl_fortran', '_fortran']
                hl += ['_hl']  # C HL library, always needed

                suffix = '.' + lpath.name.split('.', 1)[1]  # in case of .dll.a
                for h in hl:
                    hlfn = lpath.parent / (lpath.name.split('.', 1)[0] + h + suffix)
                    if hlfn.is_file():
                        link_args.append(str(hlfn))
                # HDF5 C libs are required by other HDF5 languages
                link_args.append(larg)
            else:
                link_args.append(larg)

        self.link_args = link_args


class HDF5ConfigToolDependency(ConfigToolDependency):

    """Wrapper around hdf5 binary config tools."""

    version_arg = '-showconfig'

    def __init__(self, name: str, environment: 'Environment', kwargs: T.Dict[str, T.Any], language: T.Optional[str] = None) -> None:
        language = language or 'c'
        if language not in {'c', 'cpp', 'fortran'}:
            raise DependencyException(f'Language {language} is not supported with HDF5.')

        if language == 'c':
            cenv = 'CC'
            lenv = 'C'
            tools = ['h5cc', 'h5pcc']
        elif language == 'cpp':
            cenv = 'CXX'
            lenv = 'CXX'
            tools = ['h5c++', 'h5pc++']
        elif language == 'fortran':
            cenv = 'FC'
            lenv = 'F'
            tools = ['h5fc', 'h5pfc']
        else:
            raise DependencyException('How did you get here?')

        # We need this before we call super()
        for_machine = self.get_for_machine_from_kwargs(kwargs)

        nkwargs = kwargs.copy()
        nkwargs['tools'] = tools

        # Override the compiler that the config tools are going to use by
        # setting the environment variables that they use for the compiler and
        # linkers.
        compiler = environment.coredata.compilers[for_machine][language]
        try:
            os.environ[f'HDF5_{cenv}'] = join_args(compiler.get_exelist())
            os.environ[f'HDF5_{lenv}LINKER'] = join_args(compiler.get_linker_exelist())
            super().__init__(name, environment, nkwargs, language)
        finally:
            del os.environ[f'HDF5_{cenv}']
            del os.environ[f'HDF5_{lenv}LINKER']
        if not self.is_found:
            return

        # We first need to call the tool with -c to get the compile arguments
        # and then without -c to get the link arguments.
        args = self.get_config_value(['-show', '-c'], 'args')[1:]
        args += self.get_config_value(['-show', '-noshlib' if self.static else '-shlib'], 'args')[1:]
        found = False
        for arg in args:
            if arg.startswith(('-I', '-f', '-D')) or arg == '-pthread':
                self.compile_args.append(arg)
            elif arg.startswith(('-L', '-l', '-Wl')):
                self.link_args.append(arg)
                found = True
            elif Path(arg).is_file():
                self.link_args.append(arg)
                found = True

        # cmake h5cc is broken
        if not found:
            raise DependencyException('HDF5 was built with cmake instead of autotools, and h5cc is broken.')

    def _sanitize_version(self, ver: str) -> str:
        v = re.search(r'\s*HDF5 Version: (\d+\.\d+\.\d+)', ver)
        return v.group(1)


@factory_methods({DependencyMethods.PKGCONFIG, DependencyMethods.CONFIG_TOOL})
def hdf5_factory(env: 'Environment', for_machine: 'MachineChoice',
                 kwargs: T.Dict[str, T.Any], methods: T.List[DependencyMethods]) -> T.List['DependencyGenerator']:
    language = kwargs.get('language')
    candidates: T.List['DependencyGenerator'] = []

    if DependencyMethods.PKGCONFIG in methods:
        # Use an ordered set so that these remain the first tried pkg-config files
        pkgconfig_files = OrderedSet(['hdf5', 'hdf5-serial'])
        pkg = PkgConfigInterface.instance(env, for_machine, silent=False)
        if pkg:
            # some distros put hdf5-1.2.3.pc with version number in .pc filename.
            for mod in pkg.list_all():
                if mod.startswith('hdf5'):
                    pkgconfig_files.add(mod)
        for mod in pkgconfig_files:
            candidates.append(functools.partial(HDF5PkgConfigDependency, mod, env, kwargs, language))

    if DependencyMethods.CONFIG_TOOL in methods:
        candidates.append(functools.partial(HDF5ConfigToolDependency, 'hdf5', env, kwargs, language))

    return candidates

packages['hdf5'] = hdf5_factory
