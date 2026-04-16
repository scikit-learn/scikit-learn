# SPDX-License-Identifier: Apache-2.0
# Copyright 2012-2022 The Meson development team
# Copyright © 2023-2025 Intel Corporation

from __future__ import annotations

import argparse
import functools
import os.path
import textwrap
import re
import typing as T

from .. import options
from ..dependencies import InternalDependency
from ..mesonlib import EnvironmentException, MesonException, Popen_safe, Popen_safe_logged, version_compare
from ..linkers.linkers import VisualStudioLikeLinkerMixin
from ..options import OptionKey
from .compilers import Compiler, CompileCheckMode, clike_debug_args, is_library

if T.TYPE_CHECKING:
    from .. import build
    from ..options import MutableKeyedOptionDictType
    from ..environment import Environment  # noqa: F401
    from ..linkers.linkers import DynamicLinker
    from ..mesonlib import MachineChoice
    from ..dependencies import Dependency
    from ..build import BuildTarget

    from typing_extensions import Protocol

    class TargetParse(Protocol):
        target: T.Optional[str]


rust_optimization_args: T.Dict[str, T.List[str]] = {
    'plain': [],
    '0': [],
    'g': ['-C', 'opt-level=0'],
    '1': ['-C', 'opt-level=1'],
    '2': ['-C', 'opt-level=2'],
    '3': ['-C', 'opt-level=3'],
    's': ['-C', 'opt-level=s'],
}


class _TargetParser:

    """Helper for bindgen to look for --target in various command line arguments.

    Storing this as a helper class avoids the need to set up the ArgumentParser
    multiple times, and simplifies it's use as well as the typing.
    """

    def __init__(self) -> None:
        parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
        parser.add_argument('--target', action='store', default=None)
        self._parser = parser

    def parse(self, args: T.List[str]) -> T.Optional[str]:
        """Parse arguments looking for --target

        :param args: A list of arguments to search
        :return: the argument to --target if it exists, otherwise None
        """
        parsed = T.cast('TargetParse', self._parser.parse_known_args(args)[0])
        return parsed.target


parse_target = _TargetParser().parse

def get_rustup_run_and_args(exelist: T.List[str]) -> T.Optional[T.Tuple[T.List[str], T.List[str]]]:
    """Given the command for a rustc executable, check if it is invoked via
       "rustup run" and if so separate the "rustup [OPTIONS] run TOOLCHAIN"
       part from the arguments to rustc.  If the returned value is not None,
       other tools (for example clippy-driver or rustdoc) can be run by placing
       the name of the tool between the two elements of the tuple."""
    e = iter(exelist)
    try:
        if os.path.basename(next(e)) != 'rustup':
            return None
        # minimum three strings: "rustup run TOOLCHAIN"
        n = 3
        opt = next(e)

        # options come first
        while opt.startswith('-'):
            n += 1
            opt = next(e)

        # then "run TOOLCHAIN"
        if opt != 'run':
            return None

        next(e)
        next(e)
        return exelist[:n], list(e)
    except StopIteration:
        return None

def rustc_link_args(args: T.List[str]) -> T.List[str]:
    if not args:
        return args
    rustc_args: T.List[str] = []
    for arg in args:
        rustc_args.append('-C')
        rustc_args.append(f'link-arg={arg}')
    return rustc_args


class RustSystemDependency(InternalDependency):
    pass


class RustCompiler(Compiler):

    # rustc doesn't invoke the compiler itself, it doesn't need a LINKER_PREFIX
    language = 'rust'
    id = 'rustc'

    USED_FOR_SEPARATE_LINKING_STEP = False

    _WARNING_LEVELS: T.Dict[str, T.List[str]] = {
        '0': ['--cap-lints', 'allow'],
        '1': [],
        '2': [],
        '3': ['-W', 'warnings'],
        'everything': ['-W', 'warnings'],
    }

    allow_nightly: bool

    # libcore can be compiled with either static or dynamic CRT, so disable
    # both of them just in case.
    MSVCRT_ARGS: T.Mapping[str, T.List[str]] = {
        'none': [],
        'md': ['-Clink-arg=/nodefaultlib:libcmt', '-Clink-arg=/defaultlib:msvcrt'],
        'mdd': ['-Clink-arg=/nodefaultlib:libcmt', '-Clink-arg=/nodefaultlib:msvcrt', '-Clink-arg=/defaultlib:msvcrtd'],
        'mt': ['-Clink-arg=/defaultlib:libcmt', '-Clink-arg=/nodefaultlib:msvcrt'],
        'mtd': ['-Clink-arg=/nodefaultlib:libcmt', '-Clink-arg=/nodefaultlib:msvcrt', '-Clink-arg=/defaultlib:libcmtd'],
    }

    def __init__(self, exelist: T.List[str], version: str, for_machine: MachineChoice,
                 env: Environment, full_version: T.Optional[str] = None,
                 linker: T.Optional['DynamicLinker'] = None):
        super().__init__([], exelist, version, for_machine, env,
                         full_version=full_version, linker=linker)
        self.rustup_run_and_args: T.Optional[T.Tuple[T.List[str], T.List[str]]] = get_rustup_run_and_args(exelist)
        self.base_options.update({OptionKey(o) for o in ['b_colorout', 'b_coverage', 'b_ndebug', 'b_pgo']})
        if isinstance(self.linker, VisualStudioLikeLinkerMixin):
            self.base_options.add(OptionKey('b_vscrt'))
        self.native_static_libs: T.List[str] = []
        self.is_beta = '-beta' in full_version
        self.is_nightly = '-nightly' in full_version
        self.has_check_cfg = version_compare(version, '>=1.80.0')

    def init_from_options(self) -> None:
        nightly_opt = self.get_compileropt_value('nightly', None)
        if nightly_opt == 'enabled' and not self.is_nightly:
            raise EnvironmentException(f'Rust compiler {self.name_string()} is not a nightly compiler as required by the "nightly" option.')
        self.allow_nightly = nightly_opt != 'disabled' and self.is_nightly

    def needs_static_linker(self) -> bool:
        return False

    def _sanity_check_compile_args(self, sourcename: str, binname: str
                                   ) -> T.Tuple[T.List[str], T.List[str]]:
        cmdlist = self.exelist.copy()
        largs: T.List[str] = []
        assert self.linker is not None, 'for mypy'
        if self.info.kernel == 'none' and 'ld.' in self.linker.id:
            largs.extend(rustc_link_args(['-nostartfiles']))
        cmdlist.extend(self.get_output_args(binname))
        cmdlist.append(sourcename)
        return cmdlist, largs

    def _sanity_check_source_code(self) -> str:
        if self.info.kernel != 'none':
            return textwrap.dedent(
                '''fn main() {
                }
                ''')
        return textwrap.dedent(
            '''#![no_std]
            #![no_main]
            #[no_mangle]
            pub fn _start() {
            }
            #[panic_handler]
            fn panic(_info: &core::panic::PanicInfo) -> ! {
                loop {}
            }
            ''')

    def sanity_check(self, work_dir: str) -> None:
        super().sanity_check(work_dir)
        source_name = self._sanity_check_filenames()[0]
        self._native_static_libs(work_dir, source_name)

    def _native_static_libs(self, work_dir: str, source_name: str) -> None:
        # Get libraries needed to link with a Rust staticlib
        if self.native_static_libs:
            return

        cmdlist = self.exelist + ['--crate-type', 'staticlib', '--print', 'native-static-libs', source_name]
        p, stdo, stde = Popen_safe_logged(cmdlist, cwd=work_dir)
        if p.returncode != 0:
            raise EnvironmentException('Rust compiler cannot compile staticlib.')
        match = re.search('native-static-libs: (.*)$', stde, re.MULTILINE)
        if not match:
            if self.info.kernel == 'none':
                # no match and kernel == none (i.e. baremetal) is a valid use case.
                # return and let native_static_libs list empty
                return
            if self.info.system == 'emscripten':
                # no match and emscripten is valid after rustc 1.84
                return
            raise EnvironmentException('Failed to find native-static-libs in Rust compiler output.')
        # Exclude some well known libraries that we don't need because they
        # are always part of C/C++ linkers. Rustc probably should not print
        # them, pkg-config for example never specify them.
        # FIXME: https://github.com/rust-lang/rust/issues/55120
        exclude = {'-lc', '-lgcc_s', '-lkernel32', '-ladvapi32', '/defaultlib:msvcrt'}
        self.native_static_libs = [i for i in match.group(1).split() if i not in exclude]

    def get_dependency_gen_args(self, outtarget: str, outfile: str) -> T.List[str]:
        return ['--emit', f'dep-info={outfile}']

    def get_output_args(self, outputname: str) -> T.List[str]:
        return ['--emit', f'link={outputname}']

    @functools.lru_cache(maxsize=None)
    def get_sysroot(self) -> str:
        cmd = self.get_exelist(ccache=False) + ['--print', 'sysroot']
        p, stdo, stde = Popen_safe_logged(cmd)
        return stdo.split('\n', maxsplit=1)[0]

    @functools.lru_cache(maxsize=None)
    def get_target_libdir(self) -> str:
        cmd = self.get_exelist(ccache=False) + ['--print', 'target-libdir']
        p, stdo, stde = Popen_safe_logged(cmd)
        return stdo.split('\n', maxsplit=1)[0]

    @functools.lru_cache(maxsize=None)
    def get_cfgs(self) -> T.List[str]:
        cmd = self.get_exelist(ccache=False) + ['--print', 'cfg']
        p, stdo, stde = Popen_safe_logged(cmd)
        return stdo.splitlines()

    @functools.lru_cache(maxsize=None)
    def get_target_triple(self) -> str:
        # First check if --target is explicitly set in the compiler command
        target = parse_target(self.get_exe_args())
        if target:
            return target
        # Fall back to parsing the host triple from `rustc -vV`
        cmd = self.get_exelist(ccache=False) + ['-vV']
        p, stdo, stde = Popen_safe(cmd)
        for line in stdo.splitlines():
            if line.startswith('host:'):
                return line.split(':', 1)[1].strip()
        raise EnvironmentException('Could not determine Rust target triple')

    @functools.lru_cache(maxsize=None)
    def get_crt_static(self) -> bool:
        return 'target_feature="crt-static"' in self.get_cfgs()

    def get_nightly(self, target: T.Optional[BuildTarget]) -> bool:
        if not target:
            return self.allow_nightly
        key = self.form_compileropt_key('nightly')
        nightly_opt = self.environment.coredata.get_option_for_target(target, key)
        if nightly_opt == 'enabled' and not self.is_nightly:
            raise EnvironmentException(f'Rust compiler {self.name_string()} is not a nightly compiler as required by the "nightly" option.')
        return nightly_opt != 'disabled' and self.is_nightly

    def sanitizer_link_args(self, target: T.Optional[BuildTarget], value: T.List[str]) -> T.List[str]:
        # Sanitizers are not supported yet for Rust code.  Nightly supports that
        # with -Zsanitizer=, but procedural macros cannot use them.  But even if
        # Rust code cannot be instrumented, we can link in the sanitizer libraries
        # for the sake of C/C++ code
        return rustc_link_args(super().sanitizer_link_args(target, value))

    def get_soname_args(self, prefix: str, shlib_name: str, suffix: str, soversion: str,
                        darwin_versions: T.Tuple[str, str]) -> T.List[str]:
        return rustc_link_args(super().get_soname_args(prefix, shlib_name, suffix, soversion, darwin_versions))

    @functools.lru_cache(maxsize=None)
    def has_verbatim(self) -> bool:
        if version_compare(self.version, '< 1.67.0'):
            return False
        # GNU ld support '-l:PATH'
        if 'ld.' in self.linker.id and self.linker.id != 'ld.wasm':
            return True
        # -l:+verbatim does not work (yet?) with MSVC link or Apple ld64
        # (https://github.com/rust-lang/rust/pull/138753).  For ld64, it
        # works together with -l:+whole_archive because -force_load (the macOS
        # equivalent of --whole-archive), receives the full path to the library
        # being linked.  However, Meson uses "bundle", not "whole_archive".
        return False

    def lib_file_to_l_arg(self, libname: str) -> T.Optional[str]:
        """Undo the effects of -l on the filename, returning the
           argument that can be passed to -l, or None if the
           library name is not supported."""
        if not is_library(libname):
            return None
        libname, ext = os.path.splitext(libname)

        # On Windows, rustc's -lfoo searches either foo.lib or libfoo.a.
        # Elsewhere, it searches both static and shared libraries and always with
        # the "lib" prefix; for simplicity just skip .lib on non-Windows.
        if self.info.is_windows():
            if ext == '.lib':
                return libname
            if ext != '.a':
                return None
        else:
            if ext == '.lib':
                return None

        if not libname.startswith('lib'):
            return None
        libname = libname[3:]
        return libname

    def get_debug_args(self, is_debug: bool) -> T.List[str]:
        return clike_debug_args[is_debug]

    def get_optimization_args(self, optimization_level: str) -> T.List[str]:
        return rust_optimization_args[optimization_level]

    def build_rpath_args(self, build_dir: str, from_dir: str, target: BuildTarget,
                         extra_paths: T.Optional[T.List[str]] = None
                         ) -> T.Tuple[T.List[str], T.Set[bytes]]:
        # add rustc's sysroot to account for rustup installations
        args, to_remove = super().build_rpath_args(
            build_dir, from_dir, target, [self.get_target_libdir()])
        return rustc_link_args(args), to_remove

    def compute_parameters_with_absolute_paths(self, parameter_list: T.List[str],
                                               build_dir: str) -> T.List[str]:
        for idx, i in enumerate(parameter_list):
            if i[:2] == '-L':
                for j in ['dependency', 'crate', 'native', 'framework', 'all']:
                    combined_len = len(j) + 3
                    if i[:combined_len] == f'-L{j}=':
                        parameter_list[idx] = i[:combined_len] + os.path.normpath(os.path.join(build_dir, i[combined_len:]))
                        break

        return parameter_list

    @classmethod
    def use_linker_args(cls, linker: str, version: str) -> T.List[str]:
        return ['-C', f'linker={linker}']

    def get_options(self) -> MutableKeyedOptionDictType:
        opts = super().get_options()

        key = self.form_compileropt_key('std')
        opts[key] = options.UserComboOption(
            self.make_option_name(key),
            'Rust edition to use',
            'none',
            choices=['none', '2015', '2018', '2021', '2024'])

        key = self.form_compileropt_key('dynamic_std')
        opts[key] = options.UserBooleanOption(
            self.make_option_name(key),
            'Whether to link Rust build targets to a dynamic libstd',
            False)

        key = self.form_compileropt_key('nightly')
        opts[key] = options.UserFeatureOption(
            self.make_option_name(key),
            "Nightly Rust compiler (enabled=required, disabled=don't use nightly feature, auto=use nightly feature if available)",
            'auto')

        return opts

    def get_dependency_compile_args(self, dep: 'Dependency') -> T.List[str]:
        if isinstance(dep, RustSystemDependency):
            return dep.get_compile_args()
        # Rust doesn't have dependency compile arguments so simply return
        # nothing here. Dependencies are linked and all required metadata is
        # provided by the linker flags.
        return []

    def get_option_std_args(self, target: BuildTarget, subproject: T.Optional[str] = None) -> T.List[str]:
        args = []
        std = self.get_compileropt_value('std', target, subproject)
        assert isinstance(std, str)
        if std != 'none':
            args.append('--edition=' + std)
        return args

    def get_crt_compile_args(self, crt_val: str, env: Environment) -> T.List[str]:
        # Rust handles this for us, we don't need to do anything
        return []

    def get_crt_link_args(self, crt_val: str, env: Environment) -> T.List[str]:
        if not isinstance(self.linker, VisualStudioLikeLinkerMixin):
            return []
        # Rustc always use non-debug Windows runtime. Inject the one selected
        # by Meson options instead.
        # https://github.com/rust-lang/rust/issues/39016
        return self.MSVCRT_ARGS[self.get_crt_val(crt_val, env)]

    def get_colorout_args(self, colortype: str) -> T.List[str]:
        if colortype in {'always', 'never', 'auto'}:
            return [f'--color={colortype}']
        raise MesonException(f'Invalid color type for rust {colortype}')

    @functools.lru_cache(maxsize=None)
    def get_linker_always_args(self) -> T.List[str]:
        return rustc_link_args(super().get_linker_always_args()) + ['-Cdefault-linker-libraries']

    def get_embed_bitcode_args(self, bitcode: bool, lto: bool) -> T.List[str]:
        if bitcode:
            return ['-C', 'embed-bitcode=yes']
        elif lto:
            return []
        else:
            return ['-C', 'embed-bitcode=no']

    def get_lto_compile_args(self, *, target: T.Optional[BuildTarget] = None, threads: int = 0,
                             mode: str = 'default') -> T.List[str]:
        if target.rust_crate_type in {'dylib', 'proc-macro'}:
            return []

        # TODO: what about -Clinker-plugin-lto?
        rustc_lto = 'lto=thin' if mode == 'thin' else 'lto'
        return ['-C', rustc_lto]

    def get_lto_link_args(self, *, target: T.Optional[BuildTarget] = None, threads: int = 0,
                          mode: str = 'default', thinlto_cache_dir: T.Optional[str] = None) -> T.List[str]:
        # no need to specify anything because the rustc command line
        # includes the result of get_lto_compile_args()
        return []

    def get_lto_obj_cache_path(self, path: str) -> T.List[str]:
        return rustc_link_args(super().get_lto_obj_cache_path(path))

    def get_coverage_args(self) -> T.List[str]:
        return ['-C', 'instrument-coverage']

    def get_coverage_link_args(self) -> T.List[str]:
        return rustc_link_args(super().get_coverage_link_args())

    def gen_vs_module_defs_args(self, defsfile: str) -> T.List[str]:
        return rustc_link_args(super().gen_vs_module_defs_args(defsfile))

    def gen_export_dynamic_link_args(self) -> T.List[str]:
        return rustc_link_args(self.linker.export_dynamic_args())

    def get_profile_generate_args(self) -> T.List[str]:
        return ['-C', 'profile-generate']

    def get_profile_use_args(self) -> T.List[str]:
        return ['-C', 'profile-use']

    @functools.lru_cache(maxsize=None)
    def get_asneeded_args(self) -> T.List[str]:
        return rustc_link_args(super().get_asneeded_args())

    def bitcode_args(self) -> T.List[str]:
        return ['-C', 'embed-bitcode=yes']

    @functools.lru_cache(maxsize=None)
    def headerpad_args(self) -> T.List[str]:
        return rustc_link_args(super().headerpad_args())

    @functools.lru_cache(maxsize=None)
    def get_allow_undefined_link_args(self) -> T.List[str]:
        return rustc_link_args(super().get_allow_undefined_link_args())

    def get_build_link_args(self, target: BuildTarget, build: build.Build) -> T.List[str]:
        return rustc_link_args(super().get_build_link_args(target, build))

    def get_target_link_args(self, target: 'BuildTarget') -> T.List[str]:
        return rustc_link_args(super().get_target_link_args(target))

    def get_win_subsystem_args(self, value: str) -> T.List[str]:
        return rustc_link_args(super().get_win_subsystem_args(value))

    def get_werror_args(self) -> T.List[str]:
        # Use -D warnings, which makes every warning not explicitly allowed an
        # error
        return ['-D', 'warnings']

    def get_warn_args(self, level: str) -> T.List[str]:
        # TODO: I'm not really sure what to put here, Rustc doesn't have warning
        return self._WARNING_LEVELS[level]

    def get_pic_args(self) -> T.List[str]:
        # relocation-model=pic is rustc's default already.
        return []

    def get_std_link_args(self, env: Environment, is_thin: bool) -> T.List[str]:
        # Rust handles static library creation via --crate-type
        return []

    def get_std_shared_lib_link_args(self) -> T.List[str]:
        # Rust handles shared library creation via --crate-type
        return []

    def get_std_shared_module_link_args(self, target: BuildTarget) -> T.List[str]:
        # Rust handles shared module creation via --crate-type
        return []

    def get_pie_args(self) -> T.List[str]:
        # Rustc currently has no way to toggle this, it's controlled by whether
        # pic is on by rustc
        return []

    def get_compile_only_args(self) -> T.List[str]:
        return ['--crate-type', 'lib']

    def get_pie_link_args(self) -> T.List[str]:
        # Rustc currently has no way to toggle this, it's controlled by whether
        # pic is on by rustc
        return []

    def get_assert_args(self, disable: bool) -> T.List[str]:
        action = "no" if disable else "yes"
        return ['-C', f'debug-assertions={action}', '-C', 'overflow-checks=no']

    def get_rust_tool(self, name: str) -> T.List[str]:
        if self.rustup_run_and_args:
            rustup_exelist, args = self.rustup_run_and_args
            # do not use extend so that exelist is copied
            exelist = rustup_exelist + [name]
        else:
            exelist = [name]
            args = self.get_exe_args()

        from ..programs import find_external_program
        for prog in find_external_program(self.environment, self.for_machine, exelist[0], exelist[0],
                                          [exelist[0]], allow_default_for_cross=False):
            exelist[0] = prog.path
            break
        else:
            return []

        return exelist + args

    def has_multi_arguments(self, args: T.List[str]) -> T.Tuple[bool, bool]:
        return self.compiles('fn main() { std::process::exit(0) }\n', extra_args=args, mode=CompileCheckMode.COMPILE)

    def has_multi_link_arguments(self, args: T.List[str], to_host_args: bool = True) -> T.Tuple[bool, bool]:
        if to_host_args:
            args = rustc_link_args(args)
        args = rustc_link_args(self.linker.fatal_warnings()) + args
        return self.compiles('fn main() { std::process::exit(0) }\n', extra_args=args, mode=CompileCheckMode.LINK)

    @functools.lru_cache(maxsize=None)
    def get_rustdoc(self) -> T.Optional[RustdocTestCompiler]:
        exelist = self.get_rust_tool('rustdoc')
        if not exelist:
            return None

        return RustdocTestCompiler(exelist, self.version, self.for_machine,
                                   self.environment,
                                   full_version=self.full_version,
                                   linker=self.linker, rustc=self)

    def enable_env_set_args(self) -> T.Optional[T.List[str]]:
        '''Extra arguments to enable --env-set support in rustc.
        Returns None if not supported.
        '''
        if version_compare(self.version, '>= 1.76') and self.allow_nightly:
            return ['-Z', 'unstable-options']
        return None


class ClippyRustCompiler(RustCompiler):

    """Clippy is a linter that wraps Rustc.

    This just provides us a different id
    """

    id = 'clippy-driver rustc'


class RustdocTestCompiler(RustCompiler):

    """We invoke Rustdoc to run doctests.  Some of the flags
       are different from rustc and some (e.g. --emit link) are
       ignored."""

    id = 'rustdoc --test'

    def __init__(self, exelist: T.List[str], version: str, for_machine: MachineChoice,
                 env: Environment, full_version: T.Optional[str],
                 linker: T.Optional['DynamicLinker'], rustc: RustCompiler):
        super().__init__(exelist, version, for_machine,
                         env, full_version, linker)
        self.rustc = rustc

    @functools.lru_cache(maxsize=None)
    def get_sysroot(self) -> str:
        return self.rustc.get_sysroot()

    @functools.lru_cache(maxsize=None)
    def get_target_libdir(self) -> str:
        return self.rustc.get_target_libdir()

    @functools.lru_cache(maxsize=None)
    def get_cfgs(self) -> T.List[str]:
        return self.rustc.get_cfgs()

    def get_debug_args(self, is_debug: bool) -> T.List[str]:
        return []

    def get_dependency_gen_args(self, outtarget: str, outfile: str) -> T.List[str]:
        return []

    def get_output_args(self, outputname: str) -> T.List[str]:
        return []
