# SPDX-License-Identifier: Apache-2.0
# Copyright 2012-2017 The Meson development team

from __future__ import annotations

import os.path
import typing as T

from .. import mlog
from .. import mesonlib
from ..mesonlib import version_compare, LibType
from ..options import OptionKey
from .compilers import CompileCheckMode, Compiler

if T.TYPE_CHECKING:
    from ..arglist import CompilerArgs
    from ..environment import Environment
    from ..mesonlib import MachineChoice
    from ..dependencies import Dependency
    from ..build import BuildTarget

class ValaCompiler(Compiler):

    language = 'vala'
    id = 'valac'

    def __init__(self, exelist: T.List[str], version: str, for_machine: MachineChoice,
                 environment: Environment):
        super().__init__([], exelist, version, for_machine, environment)
        self.version = version
        self.base_options = {OptionKey('b_colorout')}
        self.force_link = False
        self._has_color_support = version_compare(self.version, '>=0.37.1')
        self._has_posix_profile = version_compare(self.version, '>= 0.44')

    def needs_static_linker(self) -> bool:
        return False # Because compiles into C.

    def get_optimization_args(self, optimization_level: str) -> T.List[str]:
        return []

    def get_dependency_gen_args(self, outtarget: str, outfile: str) -> T.List[str]:
        if version_compare(self.version, '>=0.47.2'):
            return ['--depfile', outfile]
        return []

    def get_depfile_suffix(self) -> str:
        return 'depfile'

    def get_debug_args(self, is_debug: bool) -> T.List[str]:
        return ['--debug'] if is_debug else []

    def get_output_args(self, outputname: str) -> T.List[str]:
        return [] # Because compiles into C.

    def get_compile_only_args(self) -> T.List[str]:
        return [] # Because compiles into C.

    def get_compiler_args_for_mode(self, mode: CompileCheckMode) -> T.List[str]:
        args: T.List[str] = []
        if mode is CompileCheckMode.LINK and self.force_link:
            return args
        args += self.get_always_args()
        if mode is CompileCheckMode.COMPILE:
            args += self.get_compile_only_args()
        elif mode is CompileCheckMode.PREPROCESS:
            args += self.get_preprocess_only_args()
        return args

    def get_preprocess_only_args(self) -> T.List[str]:
        return []

    def get_pic_args(self) -> T.List[str]:
        return []

    def get_pie_args(self) -> T.List[str]:
        return []

    def get_pie_link_args(self) -> T.List[str]:
        return []

    def get_always_args(self) -> T.List[str]:
        return ['-C']

    def get_warn_args(self, level: str) -> T.List[str]:
        return []

    def get_werror_args(self) -> T.List[str]:
        return ['--fatal-warnings']

    def get_colorout_args(self, colortype: str) -> T.List[str]:
        if self._has_color_support:
            return ['--color=' + colortype]
        return []

    def compute_parameters_with_absolute_paths(self, parameter_list: T.List[str],
                                               build_dir: str) -> T.List[str]:
        for idx, i in enumerate(parameter_list):
            if i[:9] == '--girdir=':
                parameter_list[idx] = i[:9] + os.path.normpath(os.path.join(build_dir, i[9:]))
            if i[:10] == '--vapidir=':
                parameter_list[idx] = i[:10] + os.path.normpath(os.path.join(build_dir, i[10:]))
            if i[:13] == '--includedir=':
                parameter_list[idx] = i[:13] + os.path.normpath(os.path.join(build_dir, i[13:]))
            if i[:14] == '--metadatadir=':
                parameter_list[idx] = i[:14] + os.path.normpath(os.path.join(build_dir, i[14:]))

        return parameter_list

    def _sanity_check_source_code(self) -> str:
        return 'public static int main() { return 0; }'

    def _sanity_check_compile_args(self, sourcename: str, binname: str
                                   ) -> T.Tuple[T.List[str], T.List[str]]:
        args, largs = super()._sanity_check_compile_args(sourcename, binname)
        if self._has_posix_profile:
            # This removes the glib requirement. Posix and libc are equivalent,
            # but posix is available in older versions of valac
            args.append('--profile=posix')
        return args, largs

    def _transpiled_sanity_check_compile_args(
            self, compiler: Compiler, sourcename: str, binname: str
            ) -> T.Tuple[T.List[str], T.List[str]]:
        args, largs = super()._transpiled_sanity_check_compile_args(compiler, sourcename, binname)
        if self._has_posix_profile:
            return args, largs

        # If valac is too old for the posix profile then we need to find goobject-2.0 for linking.
        from ..dependencies import find_external_dependency
        with mlog.no_logging():
            dep = find_external_dependency('gobject-2.0', self.environment,
                                           {'required': False, 'native': self.for_machine})
        if not dep.found():
            raise mesonlib.EnvironmentException(
                'Valac < 0.44 requires gobject-2.0 for link testing, bit it could not be found.')

        args.extend(dep.get_all_compile_args())
        largs.extend(dep.get_all_link_args())
        return args, largs

    def _sanity_check_filenames(self) -> T.Tuple[str, T.Optional[str], str]:
        sourcename, _, binname = super()._sanity_check_filenames()
        return sourcename, f'{os.path.splitext(sourcename)[0]}.c', binname

    def find_library(self, libname: str, extra_dirs: T.List[str], libtype: LibType = LibType.PREFER_SHARED,
                     lib_prefix_warning: bool = True, ignore_system_dirs: bool = False,
                     skip_link_check: bool = False) -> T.Optional[T.List[str]]:
        if extra_dirs and isinstance(extra_dirs, str):
            extra_dirs = [extra_dirs]
        # Valac always looks in the default vapi dir, so only search there if
        # no extra dirs are specified.
        if not extra_dirs:
            code = 'class MesonFindLibrary : Object { }'
            args: T.List[str] = []
            args += self.environment.coredata.get_external_args(self.for_machine, self.language)
            vapi_args = ['--pkg', libname]
            args += vapi_args
            with self.cached_compile(code, extra_args=args, mode=CompileCheckMode.COMPILE) as p:
                if p.returncode == 0:
                    return vapi_args
        # Not found? Try to find the vapi file itself.
        for d in extra_dirs:
            vapi = os.path.join(d, libname + '.vapi')
            if os.path.isfile(vapi):
                return [vapi]
        mlog.debug(f'Searched {extra_dirs!r} and {libname!r} wasn\'t found')
        return None

    def thread_flags(self) -> T.List[str]:
        return []

    def thread_link_flags(self) -> T.List[str]:
        return []

    def get_option_link_args(self, target: 'BuildTarget', subproject: T.Optional[str] = None) -> T.List[str]:
        return []

    def build_wrapper_args(self,
                           extra_args: T.Union[None, CompilerArgs, T.List[str], T.Callable[[CompileCheckMode], T.List[str]]],
                           dependencies: T.Optional[T.List['Dependency']],
                           mode: CompileCheckMode = CompileCheckMode.COMPILE) -> CompilerArgs:
        if callable(extra_args):
            extra_args = extra_args(mode)
        if extra_args is None:
            extra_args = []
        if dependencies is None:
            dependencies = []

        # Collect compiler arguments
        args = self.compiler_args(self.get_compiler_check_args(mode))
        for d in dependencies:
            # Add compile flags needed by dependencies
            if mode is CompileCheckMode.LINK and self.force_link:
                # As we are passing the parameter to valac we don't need the dependent libraries.
                a = d.get_compile_args()
                if a:
                    p = a[0]
                    n = p[max(p.rfind('/'), p.rfind('\\'))+1:]
                    if not n == d.get_name():
                        args += ['--pkg=' + d.get_name()] # This is used by gio-2.0 among others.
                    else:
                        args += ['--pkg=' + n]
                else:
                    args += ['--Xcc=-l' + d.get_name()] # This is used by the maths library(-lm) among others.
            else:
                args += d.get_compile_args()
            if mode is CompileCheckMode.LINK:
                # Add link flags needed to find dependencies
                if not self.force_link: # There are no need for link dependencies when linking with valac.
                    args += d.get_link_args()

        if mode is CompileCheckMode.COMPILE:
            # Add DFLAGS from the env
            args += self.environment.coredata.get_external_args(self.for_machine, self.language)
        elif mode is CompileCheckMode.LINK:
            # Add LDFLAGS from the env
            args += self.environment.coredata.get_external_link_args(self.for_machine, self.language)
        # extra_args must override all other arguments, so we add them last
        args += extra_args
        return args

    def links(self, code: 'mesonlib.FileOrString', *,
              compiler: T.Optional['Compiler'] = None,
              extra_args: T.Union[None, T.List[str], CompilerArgs, T.Callable[[CompileCheckMode], T.List[str]]] = None,
              dependencies: T.Optional[T.List['Dependency']] = None,
              disable_cache: bool = False) -> T.Tuple[bool, bool]:
        self.force_link = True
        if compiler:
            with compiler._build_wrapper(code, dependencies=dependencies, want_output=True) as r:
                objfile = mesonlib.File.from_absolute_file(r.output_name)
                result = self.compiles(objfile, extra_args=extra_args,
                                       dependencies=dependencies, mode=CompileCheckMode.LINK, disable_cache=True)
                self.force_link = False
                return result
        result = self.compiles(code, extra_args=extra_args,
                               dependencies=dependencies, mode=CompileCheckMode.LINK, disable_cache=disable_cache)
        self.force_link = False
        return result
