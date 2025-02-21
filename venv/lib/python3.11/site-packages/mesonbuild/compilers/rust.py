# SPDX-License-Identifier: Apache-2.0
# Copyright 2012-2022 The Meson development team

from __future__ import annotations

import functools
import subprocess, os.path
import textwrap
import re
import typing as T

from .. import options
from ..mesonlib import EnvironmentException, MesonException, Popen_safe_logged
from ..options import OptionKey
from .compilers import Compiler, clike_debug_args

if T.TYPE_CHECKING:
    from ..coredata import MutableKeyedOptionDictType, KeyedOptionDictType
    from ..envconfig import MachineInfo
    from ..environment import Environment  # noqa: F401
    from ..linkers.linkers import DynamicLinker
    from ..mesonlib import MachineChoice
    from ..dependencies import Dependency


rust_optimization_args: T.Dict[str, T.List[str]] = {
    'plain': [],
    '0': [],
    'g': ['-C', 'opt-level=0'],
    '1': ['-C', 'opt-level=1'],
    '2': ['-C', 'opt-level=2'],
    '3': ['-C', 'opt-level=3'],
    's': ['-C', 'opt-level=s'],
}

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

class RustCompiler(Compiler):

    # rustc doesn't invoke the compiler itself, it doesn't need a LINKER_PREFIX
    language = 'rust'
    id = 'rustc'

    _WARNING_LEVELS: T.Dict[str, T.List[str]] = {
        '0': ['--cap-lints', 'allow'],
        '1': [],
        '2': [],
        '3': ['-W', 'warnings'],
        'everything': ['-W', 'warnings'],
    }

    # Those are static libraries, but we use dylib= here as workaround to avoid
    # rust --tests to use /WHOLEARCHIVE.
    # https://github.com/rust-lang/rust/issues/116910
    MSVCRT_ARGS: T.Mapping[str, T.List[str]] = {
        'none': [],
        'md': [], # this is the default, no need to inject anything
        'mdd': ['-l', 'dylib=msvcrtd'],
        'mt': ['-l', 'dylib=libcmt'],
        'mtd': ['-l', 'dylib=libcmtd'],
    }

    def __init__(self, exelist: T.List[str], version: str, for_machine: MachineChoice,
                 is_cross: bool, info: 'MachineInfo',
                 full_version: T.Optional[str] = None,
                 linker: T.Optional['DynamicLinker'] = None):
        super().__init__([], exelist, version, for_machine, info,
                         is_cross=is_cross, full_version=full_version,
                         linker=linker)
        self.rustup_run_and_args: T.Optional[T.Tuple[T.List[str], T.List[str]]] = get_rustup_run_and_args(exelist)
        self.base_options.update({OptionKey(o) for o in ['b_colorout', 'b_ndebug']})
        if 'link' in self.linker.id:
            self.base_options.add(OptionKey('b_vscrt'))
        self.native_static_libs: T.List[str] = []

    def needs_static_linker(self) -> bool:
        return False

    def sanity_check(self, work_dir: str, environment: 'Environment') -> None:
        source_name = os.path.join(work_dir, 'sanity.rs')
        output_name = os.path.join(work_dir, 'rusttest')
        cmdlist = self.exelist.copy()

        with open(source_name, 'w', encoding='utf-8') as ofile:
            # If machine kernel is not `none`, try to compile a dummy program.
            # If 'none', this is likely a `no-std`(i.e. bare metal) project.
            if self.info.kernel != 'none':
                ofile.write(textwrap.dedent(
                    '''fn main() {
                    }
                    '''))
            else:
                # If rustc linker is gcc, add `-nostartfiles`
                if 'ld.' in self.linker.id:
                    cmdlist.extend(['-C', 'link-arg=-nostartfiles'])
                ofile.write(textwrap.dedent(
                    '''#![no_std]
                    #![no_main]
                    #[no_mangle]
                    pub fn _start() {
                    }
                    #[panic_handler]
                    fn panic(_info: &core::panic::PanicInfo) -> ! {
                        loop {}
                    }
                    '''))

        cmdlist.extend(['-o', output_name, source_name])
        pc, stdo, stde = Popen_safe_logged(cmdlist, cwd=work_dir)
        if pc.returncode != 0:
            raise EnvironmentException(f'Rust compiler {self.name_string()} cannot compile programs.')
        self._native_static_libs(work_dir, source_name)
        if self.is_cross:
            if not environment.has_exe_wrapper():
                # Can't check if the binaries run so we have to assume they do
                return
            cmdlist = environment.exe_wrapper.get_command() + [output_name]
        else:
            cmdlist = [output_name]
        pe = subprocess.Popen(cmdlist, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        pe.wait()
        if pe.returncode != 0:
            raise EnvironmentException(f'Executables created by Rust compiler {self.name_string()} are not runnable.')

    def _native_static_libs(self, work_dir: str, source_name: str) -> None:
        # Get libraries needed to link with a Rust staticlib
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
            raise EnvironmentException('Failed to find native-static-libs in Rust compiler output.')
        # Exclude some well known libraries that we don't need because they
        # are always part of C/C++ linkers. Rustc probably should not print
        # them, pkg-config for example never specify them.
        # FIXME: https://github.com/rust-lang/rust/issues/55120
        exclude = {'-lc', '-lgcc_s', '-lkernel32', '-ladvapi32'}
        self.native_static_libs = [i for i in match.group(1).split() if i not in exclude]

    def get_dependency_gen_args(self, outtarget: str, outfile: str) -> T.List[str]:
        return ['--dep-info', outfile]

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
    def get_crt_static(self) -> bool:
        cmd = self.get_exelist(ccache=False) + ['--print', 'cfg']
        p, stdo, stde = Popen_safe_logged(cmd)
        return bool(re.search('^target_feature="crt-static"$', stdo, re.MULTILINE))

    def get_debug_args(self, is_debug: bool) -> T.List[str]:
        return clike_debug_args[is_debug]

    def get_optimization_args(self, optimization_level: str) -> T.List[str]:
        return rust_optimization_args[optimization_level]

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

    def get_output_args(self, outputname: str) -> T.List[str]:
        return ['-o', outputname]

    @classmethod
    def use_linker_args(cls, linker: str, version: str) -> T.List[str]:
        return ['-C', f'linker={linker}']

    # Rust does not have a use_linker_args because it dispatches to a gcc-like
    # C compiler for dynamic linking, as such we invoke the C compiler's
    # use_linker_args method instead.

    def get_options(self) -> MutableKeyedOptionDictType:
        return dict((self.create_option(options.UserComboOption,
                                        self.form_compileropt_key('std'),
                                        'Rust edition to use',
                                        ['none', '2015', '2018', '2021', '2024'],
                                        'none'),))

    def get_dependency_compile_args(self, dep: 'Dependency') -> T.List[str]:
        # Rust doesn't have dependency compile arguments so simply return
        # nothing here. Dependencies are linked and all required metadata is
        # provided by the linker flags.
        return []

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        args = []
        key = self.form_compileropt_key('std')
        std = options.get_value(key)
        if std != 'none':
            args.append('--edition=' + std)
        return args

    def get_crt_compile_args(self, crt_val: str, buildtype: str) -> T.List[str]:
        # Rust handles this for us, we don't need to do anything
        return []

    def get_crt_link_args(self, crt_val: str, buildtype: str) -> T.List[str]:
        if self.linker.id not in {'link', 'lld-link'}:
            return []
        return self.MSVCRT_ARGS[self.get_crt_val(crt_val, buildtype)]

    def get_colorout_args(self, colortype: str) -> T.List[str]:
        if colortype in {'always', 'never', 'auto'}:
            return [f'--color={colortype}']
        raise MesonException(f'Invalid color type for rust {colortype}')

    def get_linker_always_args(self) -> T.List[str]:
        args: T.List[str] = []
        for a in super().get_linker_always_args():
            args.extend(['-C', f'link-arg={a}'])
        return args

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

    def get_pie_args(self) -> T.List[str]:
        # Rustc currently has no way to toggle this, it's controlled by whether
        # pic is on by rustc
        return []

    def get_assert_args(self, disable: bool, env: 'Environment') -> T.List[str]:
        action = "no" if disable else "yes"
        return ['-C', f'debug-assertions={action}', '-C', 'overflow-checks=no']

    def get_rust_tool(self, name: str, env: Environment) -> T.List[str]:
        if self.rustup_run_and_args:
            rustup_exelist, args = self.rustup_run_and_args
            # do not use extend so that exelist is copied
            exelist = rustup_exelist + [name]
        else:
            exelist = [name]
            args = self.exelist[1:]

        from ..programs import find_external_program
        for prog in find_external_program(env, self.for_machine, exelist[0], exelist[0],
                                          [exelist[0]], allow_default_for_cross=False):
            exelist[0] = prog.path
            break
        else:
            return []

        return exelist + args


class ClippyRustCompiler(RustCompiler):

    """Clippy is a linter that wraps Rustc.

    This just provides us a different id
    """

    id = 'clippy-driver rustc'
