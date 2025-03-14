# SPDX-License-Identifier: Apache-2.0
# Copyright 2012-2022 The Meson development team

from __future__ import annotations

from ..mesonlib import (
    MesonException, EnvironmentException, MachineChoice, join_args,
    search_version, is_windows, Popen_safe, Popen_safe_logged, windows_proof_rm,
)
from ..envconfig import BinaryTable
from .. import mlog

from ..linkers import guess_win_linker, guess_nix_linker

import subprocess
import platform
import re
import shutil
import tempfile
import os
import typing as T

if T.TYPE_CHECKING:
    from .compilers import Compiler
    from .c import CCompiler
    from .cpp import CPPCompiler
    from .fortran import FortranCompiler
    from .rust import RustCompiler
    from ..linkers.linkers import StaticLinker, DynamicLinker
    from ..environment import Environment


# Default compilers and linkers
# =============================

defaults: T.Dict[str, T.List[str]] = {}

# List of potential compilers.
if is_windows():
    # Intel C and C++ compiler is icl on Windows, but icc and icpc elsewhere.
    # Search for icl before cl, since Intel "helpfully" provides a
    # cl.exe that returns *exactly the same thing* that Microsoft's
    # cl.exe does, and if icl is present, it's almost certainly what
    # you want.
    defaults['c'] = ['icl', 'cl', 'cc', 'gcc', 'clang', 'clang-cl', 'pgcc']
    # There is currently no pgc++ for Windows, only for  Mac and Linux.
    defaults['cpp'] = ['icl', 'cl', 'c++', 'g++', 'clang++', 'clang-cl']
    # the binary flang-new will be renamed to flang in the foreseeable future
    defaults['fortran'] = ['ifort', 'gfortran', 'flang-new', 'flang', 'pgfortran', 'g95']
    defaults['objc'] = ['clang', 'clang-cl', 'gcc']
    defaults['objcpp'] = ['clang++', 'clang-cl', 'g++']
    defaults['cs'] = ['csc', 'mcs']
else:
    if platform.machine().lower() == 'e2k':
        defaults['c'] = ['cc', 'gcc', 'lcc', 'clang']
        defaults['cpp'] = ['c++', 'g++', 'l++', 'clang++']
        defaults['objc'] = ['clang']
        defaults['objcpp'] = ['clang++']
    else:
        defaults['c'] = ['cc', 'gcc', 'clang', 'nvc', 'pgcc', 'icc', 'icx']
        defaults['cpp'] = ['c++', 'g++', 'clang++', 'nvc++', 'pgc++', 'icpc', 'icpx']
        defaults['objc'] = ['clang', 'gcc']
        defaults['objcpp'] = ['clang++', 'g++']
    # the binary flang-new will be renamed to flang in the foreseeable future
    defaults['fortran'] = ['gfortran', 'flang-new', 'flang', 'nvfortran', 'pgfortran', 'ifort', 'ifx', 'g95']
    defaults['cs'] = ['mcs', 'csc']
defaults['d'] = ['ldc2', 'ldc', 'gdc', 'dmd']
defaults['java'] = ['javac']
defaults['cuda'] = ['nvcc']
defaults['rust'] = ['rustc']
defaults['swift'] = ['swiftc']
defaults['vala'] = ['valac']
defaults['cython'] = ['cython', 'cython3'] # Official name is cython, but Debian renamed it to cython3.
defaults['static_linker'] = ['ar', 'gar']
defaults['strip'] = ['strip']
defaults['vs_static_linker'] = ['lib']
defaults['clang_cl_static_linker'] = ['llvm-lib']
defaults['cuda_static_linker'] = ['nvlink']
defaults['gcc_static_linker'] = ['gcc-ar']
defaults['clang_static_linker'] = ['llvm-ar']
defaults['nasm'] = ['nasm', 'yasm']


def compiler_from_language(env: 'Environment', lang: str, for_machine: MachineChoice) -> T.Optional[Compiler]:
    lang_map: T.Dict[str, T.Callable[['Environment', MachineChoice], Compiler]] = {
        'c': detect_c_compiler,
        'cpp': detect_cpp_compiler,
        'objc': detect_objc_compiler,
        'cuda': detect_cuda_compiler,
        'objcpp': detect_objcpp_compiler,
        'java': detect_java_compiler,
        'cs': detect_cs_compiler,
        'vala': detect_vala_compiler,
        'd': detect_d_compiler,
        'rust': detect_rust_compiler,
        'fortran': detect_fortran_compiler,
        'swift': detect_swift_compiler,
        'cython': detect_cython_compiler,
        'nasm': detect_nasm_compiler,
        'masm': detect_masm_compiler,
        'linearasm': detect_linearasm_compiler,
    }
    return lang_map[lang](env, for_machine) if lang in lang_map else None

def detect_compiler_for(env: 'Environment', lang: str, for_machine: MachineChoice, skip_sanity_check: bool, subproject: str) -> T.Optional[Compiler]:
    comp = compiler_from_language(env, lang, for_machine)
    if comp is None:
        return comp
    assert comp.for_machine == for_machine
    env.coredata.process_compiler_options(lang, comp, env, subproject)
    if not skip_sanity_check:
        comp.sanity_check(env.get_scratch_dir(), env)
    env.coredata.compilers[comp.for_machine][lang] = comp
    return comp


# Helpers
# =======

def _get_compilers(env: 'Environment', lang: str, for_machine: MachineChoice,
                   allow_build_machine: bool = False) -> T.Tuple[T.List[T.List[str]], T.List[str]]:
    '''
    The list of compilers is detected in the exact same way for
    C, C++, ObjC, ObjC++, Fortran, CS so consolidate it here.
    '''
    value = env.lookup_binary_entry(for_machine, lang)
    if value is not None:
        comp, ccache = BinaryTable.parse_entry(value)
        # Return value has to be a list of compiler 'choices'
        compilers = [comp]
    else:
        if not env.machines.matches_build_machine(for_machine):
            if allow_build_machine:
                return _get_compilers(env, lang, MachineChoice.BUILD)
            raise EnvironmentException(f'{lang!r} compiler binary not defined in cross file [binaries] section')
        compilers = [[x] for x in defaults[lang]]
        ccache = BinaryTable.detect_compiler_cache()

    return compilers, ccache

def _handle_exceptions(
        exceptions: T.Mapping[str, T.Union[Exception, str]],
        binaries: T.List[T.List[str]],
        bintype: str = 'compiler') -> T.NoReturn:
    errmsg = f'Unknown {bintype}(s): {binaries}'
    if exceptions:
        errmsg += '\nThe following exception(s) were encountered:'
        for c, e in exceptions.items():
            errmsg += f'\nRunning `{c}` gave "{e}"'
    raise EnvironmentException(errmsg)


# Linker specific
# ===============

def detect_static_linker(env: 'Environment', compiler: Compiler) -> StaticLinker:
    from . import d
    from ..linkers import linkers
    linker = env.lookup_binary_entry(compiler.for_machine, 'ar')
    if linker is not None:
        trials = [linker]
    else:
        default_linkers = [[l] for l in defaults['static_linker']]
        if compiler.language == 'cuda':
            trials = [defaults['cuda_static_linker']] + default_linkers
        elif compiler.get_argument_syntax() == 'msvc':
            trials = [defaults['vs_static_linker'], defaults['clang_cl_static_linker']]
        elif compiler.id == 'gcc':
            # Use gcc-ar if available; needed for LTO
            trials = [defaults['gcc_static_linker']] + default_linkers
        elif compiler.id == 'clang':
            # Use llvm-ar if available; needed for LTO
            llvm_ar = defaults['clang_static_linker']
            # Extract the version major of the compiler to use as a suffix
            suffix = compiler.version.split('.')[0]
            # Prefer suffixed llvm-ar first, then unsuffixed then the defaults
            trials = [[f'{llvm_ar[0]}-{suffix}'], llvm_ar] + default_linkers
        elif compiler.language == 'd':
            # Prefer static linkers over linkers used by D compilers
            if is_windows():
                trials = [defaults['vs_static_linker'], defaults['clang_cl_static_linker'], compiler.get_linker_exelist()]
            else:
                trials = default_linkers
        elif compiler.id == 'intel-cl' and compiler.language == 'c': # why not cpp? Is this a bug?
            # Intel has its own linker that acts like Microsoft's lib
            trials = [['xilib']]
        elif is_windows() and compiler.id == 'pgi': # this handles cpp / nvidia HPC, in addition to just c/fortran
            trials = [['ar']]  # For PGI on Windows, "ar" is just a wrapper calling link/lib.
        elif is_windows() and compiler.id == 'nasm':
            # This may well be LINK.EXE if it's under a MSVC environment
            trials = [defaults['vs_static_linker'], defaults['clang_cl_static_linker']] + default_linkers
        else:
            trials = default_linkers
    popen_exceptions = {}
    for linker in trials:
        linker_name = os.path.basename(linker[0])

        if any(os.path.basename(x) in {'lib', 'lib.exe', 'llvm-lib', 'llvm-lib.exe', 'xilib', 'xilib.exe'} for x in linker):
            arg = '/?'
        elif linker_name in {'ar2000', 'ar2000.exe', 'ar430', 'ar430.exe', 'armar', 'armar.exe', 'ar6x', 'ar6x.exe'}:
            arg = '?'
        else:
            arg = '--version'
        try:
            p, out, err = Popen_safe_logged(linker + [arg], msg='Detecting archiver via')
        except OSError as e:
            popen_exceptions[join_args(linker + [arg])] = e
            continue
        if "xilib: executing 'lib'" in err:
            return linkers.IntelVisualStudioLinker(linker, getattr(compiler, 'machine', None))
        if '/OUT:' in out.upper() or '/OUT:' in err.upper():
            return linkers.VisualStudioLinker(linker, getattr(compiler, 'machine', None))
        if 'ar-Error-Unknown switch: --version' in err:
            return linkers.PGIStaticLinker(linker)
        if p.returncode == 0 and 'armar' in linker_name:
            return linkers.ArmarLinker(linker)
        if 'DMD32 D Compiler' in out or 'DMD64 D Compiler' in out:
            assert isinstance(compiler, d.DCompiler)
            return linkers.DLinker(linker, compiler.arch)
        if 'LDC - the LLVM D compiler' in out:
            assert isinstance(compiler, d.DCompiler)
            return linkers.DLinker(linker, compiler.arch, rsp_syntax=compiler.rsp_file_syntax())
        if 'GDC' in out and ' based on D ' in out:
            assert isinstance(compiler, d.DCompiler)
            return linkers.DLinker(linker, compiler.arch)
        if err.startswith('Renesas') and 'rlink' in linker_name:
            return linkers.CcrxLinker(linker)
        if out.startswith('GNU ar') and 'xc16-ar' in linker_name:
            return linkers.Xc16Linker(linker)
        if 'Texas Instruments Incorporated' in out:
            if 'ar2000' in linker_name:
                return linkers.C2000Linker(linker)
            elif 'ar6000' in linker_name:
                return linkers.C6000Linker(linker)
            else:
                return linkers.TILinker(linker)
        if out.startswith('The CompCert'):
            return linkers.CompCertLinker(linker)
        if out.strip().startswith('Metrowerks') or out.strip().startswith('Freescale'):
            if 'ARM' in out:
                return linkers.MetrowerksStaticLinkerARM(linker)
            else:
                return linkers.MetrowerksStaticLinkerEmbeddedPowerPC(linker)
        if 'TASKING VX-toolset' in err:
            return linkers.TaskingStaticLinker(linker)
        if p.returncode == 0:
            return linkers.ArLinker(compiler.for_machine, linker)
        if p.returncode == 1 and err.startswith('usage'): # OSX
            return linkers.AppleArLinker(compiler.for_machine, linker)
        if p.returncode == 1 and err.startswith('Usage'): # AIX
            return linkers.AIXArLinker(linker)
        if p.returncode == 1 and err.startswith('ar: bad option: --'): # Solaris
            return linkers.ArLinker(compiler.for_machine, linker)
    _handle_exceptions(popen_exceptions, trials, 'linker')
    raise EnvironmentException('Unreachable code (exception to make mypy happy)')


# Compilers
# =========


def _detect_c_or_cpp_compiler(env: 'Environment', lang: str, for_machine: MachineChoice, *, override_compiler: T.Optional[T.List[str]] = None) -> Compiler:
    """Shared implementation for finding the C or C++ compiler to use.

    the override_compiler option is provided to allow compilers which use
    the compiler (GCC or Clang usually) as their shared linker, to find
    the linker they need.
    """
    from . import c, cpp
    from ..linkers import linkers
    popen_exceptions: T.Dict[str, T.Union[Exception, str]] = {}
    compilers, ccache = _get_compilers(env, lang, for_machine)
    if override_compiler is not None:
        compilers = [override_compiler]
    is_cross = env.is_cross_build(for_machine)
    info = env.machines[for_machine]
    cls: T.Union[T.Type[CCompiler], T.Type[CPPCompiler]]
    lnk: T.Union[T.Type[StaticLinker], T.Type[DynamicLinker]]

    for compiler in compilers:
        if isinstance(compiler, str):
            compiler = [compiler]
        compiler_name = os.path.basename(compiler[0])

        if any(os.path.basename(x) in {'cl', 'cl.exe', 'clang-cl', 'clang-cl.exe'} for x in compiler):
            # Watcom C provides its own cl.exe clone that mimics an older
            # version of Microsoft's compiler. Since Watcom's cl.exe is
            # just a wrapper, we skip using it if we detect its presence
            # so as not to confuse Meson when configuring for MSVC.
            #
            # Additionally the help text of Watcom's cl.exe is paged, and
            # the binary will not exit without human intervention. In
            # practice, Meson will block waiting for Watcom's cl.exe to
            # exit, which requires user input and thus will never exit.
            if 'WATCOM' in os.environ:
                def sanitize(p: T.Optional[str]) -> T.Optional[str]:
                    return os.path.normcase(os.path.abspath(p)) if p else None

                watcom_cls = [sanitize(os.path.join(os.environ['WATCOM'], 'BINNT', 'cl')),
                              sanitize(os.path.join(os.environ['WATCOM'], 'BINNT', 'cl.exe')),
                              sanitize(os.path.join(os.environ['WATCOM'], 'BINNT64', 'cl')),
                              sanitize(os.path.join(os.environ['WATCOM'], 'BINNT64', 'cl.exe'))]
                found_cl = sanitize(shutil.which('cl'))
                if found_cl in watcom_cls:
                    mlog.debug('Skipping unsupported cl.exe clone at:', found_cl)
                    continue
            arg = '/?'
        elif 'armcc' in compiler_name:
            arg = '--vsn'
        elif 'ccrx' in compiler_name:
            arg = '-v'
        elif 'xc16' in compiler_name:
            arg = '--version'
        elif 'ccomp' in compiler_name:
            arg = '-version'
        elif compiler_name in {'cl2000', 'cl2000.exe', 'cl430', 'cl430.exe', 'armcl', 'armcl.exe', 'cl6x', 'cl6x.exe'}:
            # TI compiler
            arg = '-version'
        elif compiler_name in {'icl', 'icl.exe'}:
            # if you pass anything to icl you get stuck in a pager
            arg = ''
        else:
            arg = '--version'

        cmd = compiler + [arg]
        try:
            p, out, err = Popen_safe_logged(cmd, msg='Detecting compiler via')
        except OSError as e:
            popen_exceptions[join_args(cmd)] = e
            continue

        if 'ccrx' in compiler_name:
            out = err

        full_version = out.split('\n', 1)[0]
        version = search_version(out)

        guess_gcc_or_lcc: T.Optional[str] = None
        if 'Free Software Foundation' in out or out.startswith('xt-'):
            guess_gcc_or_lcc = 'gcc'
        if 'e2k' in out and 'lcc' in out:
            guess_gcc_or_lcc = 'lcc'
        if 'Microchip Technology' in out:
            # this output has "Free Software Foundation" in its version
            guess_gcc_or_lcc = None

        if guess_gcc_or_lcc:
            defines = _get_gnu_compiler_defines(compiler, lang)
            if not defines:
                popen_exceptions[join_args(compiler)] = 'no pre-processor defines'
                continue

            if guess_gcc_or_lcc == 'lcc':
                version = _get_lcc_version_from_defines(defines)
                cls = c.ElbrusCCompiler if lang == 'c' else cpp.ElbrusCPPCompiler
            else:
                version = _get_gnu_version_from_defines(defines)
                cls = c.GnuCCompiler if lang == 'c' else cpp.GnuCPPCompiler

            linker = guess_nix_linker(env, compiler, cls, version, for_machine)

            return cls(
                ccache, compiler, version, for_machine, is_cross,
                info, defines=defines, full_version=full_version,
                linker=linker)

        if 'Emscripten' in out:
            cls = c.EmscriptenCCompiler if lang == 'c' else cpp.EmscriptenCPPCompiler
            env.coredata.add_lang_args(cls.language, cls, for_machine, env)

            # emcc requires a file input in order to pass arguments to the
            # linker. It'll exit with an error code, but still print the
            # linker version.
            with tempfile.NamedTemporaryFile(suffix='.c') as f:
                cmd = compiler + [cls.LINKER_PREFIX + "--version", f.name]
                _, o, _ = Popen_safe(cmd)

            linker = linkers.WASMDynamicLinker(
                compiler, for_machine, cls.LINKER_PREFIX,
                [], version=search_version(o))
            return cls(
                ccache, compiler, version, for_machine, is_cross, info,
                linker=linker, full_version=full_version)

        if 'Arm C/C++/Fortran Compiler' in out:
            arm_ver_match = re.search(r'version (\d+)\.(\d+)\.?(\d+)? \(build number (\d+)\)', out)
            assert arm_ver_match is not None, 'for mypy'  # because mypy *should* be complaining that this could be None
            version = '.'.join([x for x in arm_ver_match.groups() if x is not None])
            if lang == 'c':
                cls = c.ArmLtdClangCCompiler
            elif lang == 'cpp':
                cls = cpp.ArmLtdClangCPPCompiler
            linker = guess_nix_linker(env, compiler, cls, version, for_machine)
            return cls(
                ccache, compiler, version, for_machine, is_cross, info,
                linker=linker)
        if 'armclang' in out:
            # The compiler version is not present in the first line of output,
            # instead it is present in second line, startswith 'Component:'.
            # So, searching for the 'Component' in out although we know it is
            # present in second line, as we are not sure about the
            # output format in future versions
            arm_ver_match = re.search('.*Component.*', out)
            if arm_ver_match is None:
                popen_exceptions[join_args(compiler)] = 'version string not found'
                continue
            arm_ver_str = arm_ver_match.group(0)
            # Override previous values
            version = search_version(arm_ver_str)
            full_version = arm_ver_str
            cls = c.ArmclangCCompiler if lang == 'c' else cpp.ArmclangCPPCompiler
            linker = linkers.ArmClangDynamicLinker(for_machine, version=version)
            env.coredata.add_lang_args(cls.language, cls, for_machine, env)
            return cls(
                ccache, compiler, version, for_machine, is_cross, info,
                full_version=full_version, linker=linker)
        if 'CL.EXE COMPATIBILITY' in out:
            # if this is clang-cl masquerading as cl, detect it as cl, not
            # clang
            arg = '--version'
            try:
                p, out, err = Popen_safe(compiler + [arg])
            except OSError as e:
                popen_exceptions[join_args(compiler + [arg])] = e
            version = search_version(out)
            match = re.search('^Target: (.*?)-', out, re.MULTILINE)
            if match:
                target = match.group(1)
            else:
                target = 'unknown target'
            cls = c.ClangClCCompiler if lang == 'c' else cpp.ClangClCPPCompiler
            linker = guess_win_linker(env, ['lld-link'], cls, version, for_machine)
            return cls(
                compiler, version, for_machine, is_cross, info, target,
                linker=linker)

        # must be detected here before clang because TI compilers contain 'clang' in their output and so that they can be detected as 'clang'
        ti_compilers = {
           'TMS320C2000 C/C++': (c.C2000CCompiler, cpp.C2000CPPCompiler, linkers.C2000DynamicLinker),
           'TMS320C6x C/C++': (c.C6000CCompiler, cpp.C6000CPPCompiler, linkers.C6000DynamicLinker),
           'TI ARM C/C++ Compiler': (c.TICCompiler, cpp.TICPPCompiler, linkers.TIDynamicLinker),
           'MSP430 C/C++': (c.TICCompiler, cpp.TICPPCompiler, linkers.TIDynamicLinker)
        }
        for identifier, compiler_classes in ti_compilers.items():
            if identifier in out:
                cls = compiler_classes[0] if lang == 'c' else compiler_classes[1]
                lnk = compiler_classes[2]
                env.coredata.add_lang_args(cls.language, cls, for_machine, env)
                linker = lnk(compiler, for_machine, version=version)
                return cls(
                    ccache, compiler, version, for_machine, is_cross, info,
                    full_version=full_version, linker=linker)

        if 'clang' in out or 'Clang' in out:
            linker = None

            defines = _get_clang_compiler_defines(compiler, lang)

            # Even if the for_machine is darwin, we could be using vanilla
            # clang.
            if 'Apple' in out:
                cls = c.AppleClangCCompiler if lang == 'c' else cpp.AppleClangCPPCompiler
            else:
                cls = c.ClangCCompiler if lang == 'c' else cpp.ClangCPPCompiler

            if 'windows' in out or env.machines[for_machine].is_windows():
                # If we're in a MINGW context this actually will use a gnu
                # style ld, but for clang on "real" windows we'll use
                # either link.exe or lld-link.exe
                try:
                    linker = guess_win_linker(env, compiler, cls, version, for_machine, invoked_directly=False)
                except MesonException:
                    pass
            if linker is None:
                linker = guess_nix_linker(env, compiler, cls, version, for_machine)

            return cls(
                ccache, compiler, version, for_machine, is_cross, info,
                defines=defines, full_version=full_version, linker=linker)

        if 'Intel(R) C++ Intel(R)' in err:
            version = search_version(err)
            target = 'x86' if 'IA-32' in err else 'x86_64'
            cls = c.IntelClCCompiler if lang == 'c' else cpp.IntelClCPPCompiler
            env.coredata.add_lang_args(cls.language, cls, for_machine, env)
            linker = linkers.XilinkDynamicLinker(for_machine, [], version=version)
            return cls(
                compiler, version, for_machine, is_cross, info, target,
                linker=linker)
        if 'Intel(R) oneAPI DPC++/C++ Compiler for applications' in err:
            version = search_version(err)
            target = 'x86' if 'IA-32' in err else 'x86_64'
            cls = c.IntelLLVMClCCompiler if lang == 'c' else cpp.IntelLLVMClCPPCompiler
            env.coredata.add_lang_args(cls.language, cls, for_machine, env)
            linker = linkers.XilinkDynamicLinker(for_machine, [], version=version)
            return cls(
                compiler, version, for_machine, is_cross, info, target,
                linker=linker)
        if 'Microsoft' in out or 'Microsoft' in err:
            # Latest versions of Visual Studio print version
            # number to stderr but earlier ones print version
            # on stdout.  Why? Lord only knows.
            # Check both outputs to figure out version.
            for lookat in [err, out]:
                version = search_version(lookat)
                if version != 'unknown version':
                    break
            else:
                raise EnvironmentException(f'Failed to detect MSVC compiler version: stderr was\n{err!r}')
            cl_signature = lookat.split('\n', maxsplit=1)[0]
            match = re.search(r'.*(x86|x64|ARM|ARM64)([^_A-Za-z0-9]|$)', cl_signature)
            if match:
                target = match.group(1)
            else:
                m = f'Failed to detect MSVC compiler target architecture: \'cl /?\' output is\n{cl_signature}'
                raise EnvironmentException(m)
            cls = c.VisualStudioCCompiler if lang == 'c' else cpp.VisualStudioCPPCompiler
            linker = guess_win_linker(env, ['link'], cls, version, for_machine)
            # As of this writing, CCache does not support MSVC but sccache does.
            if 'sccache' not in ccache:
                ccache = []
            return cls(
                ccache, compiler, version, for_machine, is_cross, info, target,
                full_version=cl_signature, linker=linker)
        if 'PGI Compilers' in out:
            cls = c.PGICCompiler if lang == 'c' else cpp.PGICPPCompiler
            env.coredata.add_lang_args(cls.language, cls, for_machine, env)
            linker = linkers.PGIDynamicLinker(compiler, for_machine, cls.LINKER_PREFIX, [], version=version)
            return cls(
                ccache, compiler, version, for_machine, is_cross,
                info, linker=linker)
        if 'NVIDIA Compilers and Tools' in out:
            cls = c.NvidiaHPC_CCompiler if lang == 'c' else cpp.NvidiaHPC_CPPCompiler
            env.coredata.add_lang_args(cls.language, cls, for_machine, env)
            linker = linkers.NvidiaHPC_DynamicLinker(compiler, for_machine, cls.LINKER_PREFIX, [], version=version)
            return cls(
                ccache, compiler, version, for_machine, is_cross,
                info, linker=linker)
        if '(ICC)' in out:
            cls = c.IntelCCompiler if lang == 'c' else cpp.IntelCPPCompiler
            l = guess_nix_linker(env, compiler, cls, version, for_machine)
            return cls(
                ccache, compiler, version, for_machine, is_cross, info,
                full_version=full_version, linker=l)
        if 'Intel(R) oneAPI' in out:
            cls = c.IntelLLVMCCompiler if lang == 'c' else cpp.IntelLLVMCPPCompiler
            l = guess_nix_linker(env, compiler, cls, version, for_machine)
            return cls(
                ccache, compiler, version, for_machine, is_cross, info,
                full_version=full_version, linker=l)
        if 'ARM' in out and not ('Metrowerks' in out or 'Freescale' in out):
            cls = c.ArmCCompiler if lang == 'c' else cpp.ArmCPPCompiler
            env.coredata.add_lang_args(cls.language, cls, for_machine, env)
            linker = linkers.ArmDynamicLinker(for_machine, version=version)
            return cls(
                ccache, compiler, version, for_machine, is_cross,
                info, full_version=full_version, linker=linker)
        if 'RX Family' in out:
            cls = c.CcrxCCompiler if lang == 'c' else cpp.CcrxCPPCompiler
            env.coredata.add_lang_args(cls.language, cls, for_machine, env)
            linker = linkers.CcrxDynamicLinker(for_machine, version=version)
            return cls(
                ccache, compiler, version, for_machine, is_cross, info,
                full_version=full_version, linker=linker)

        if 'Microchip Technology' in out:
            cls = c.Xc16CCompiler
            env.coredata.add_lang_args(cls.language, cls, for_machine, env)
            linker = linkers.Xc16DynamicLinker(for_machine, version=version)
            return cls(
                ccache, compiler, version, for_machine, is_cross, info,
                full_version=full_version, linker=linker)

        if 'CompCert' in out:
            cls = c.CompCertCCompiler
            env.coredata.add_lang_args(cls.language, cls, for_machine, env)
            linker = linkers.CompCertDynamicLinker(for_machine, version=version)
            return cls(
                ccache, compiler, version, for_machine, is_cross, info,
                full_version=full_version, linker=linker)

        if 'Metrowerks C/C++' in out or 'Freescale C/C++' in out:
            if 'ARM' in out:
                cls = c.MetrowerksCCompilerARM if lang == 'c' else cpp.MetrowerksCPPCompilerARM
                lnk = linkers.MetrowerksLinkerARM
            else:
                cls = c.MetrowerksCCompilerEmbeddedPowerPC if lang == 'c' else cpp.MetrowerksCPPCompilerEmbeddedPowerPC
                lnk = linkers.MetrowerksLinkerEmbeddedPowerPC

            mwcc_ver_match = re.search(r'Version (\d+)\.(\d+)\.?(\d+)? build (\d+)', out)
            assert mwcc_ver_match is not None, 'for mypy'  # because mypy *should* be complaining that this could be None
            compiler_version = '.'.join(x for x in mwcc_ver_match.groups() if x is not None)

            env.coredata.add_lang_args(cls.language, cls, for_machine, env)
            ld = env.lookup_binary_entry(for_machine, cls.language + '_ld')

            if ld is not None:
                _, o_ld, _ = Popen_safe(ld + ['--version'])

                mwld_ver_match = re.search(r'Version (\d+)\.(\d+)\.?(\d+)? build (\d+)', o_ld)
                assert mwld_ver_match is not None, 'for mypy'  # because mypy *should* be complaining that this could be None
                linker_version = '.'.join(x for x in mwld_ver_match.groups() if x is not None)

                linker = lnk(ld, for_machine, version=linker_version)
            else:
                raise EnvironmentException(f'Failed to detect linker for {cls.id!r} compiler. Please update your cross file(s).')

            return cls(
                ccache, compiler, compiler_version, for_machine, is_cross, info,
                full_version=full_version, linker=linker)
        if 'TASKING VX-toolset' in err:
            cls = c.TaskingCCompiler
            lnk = linkers.TaskingLinker

            tasking_ver_match = re.search(r'v([0-9]+)\.([0-9]+)r([0-9]+) Build ([0-9]+)', err)
            assert tasking_ver_match is not None, 'for mypy'
            tasking_version = '.'.join(x for x in tasking_ver_match.groups() if x is not None)

            env.coredata.add_lang_args(cls.language, cls, for_machine, env)
            ld = env.lookup_binary_entry(for_machine, cls.language + '_ld')
            if ld is None:
                raise MesonException(f'{cls.language}_ld was not properly defined in your cross file')

            linker = lnk(ld, for_machine, version=tasking_version)
            return cls(
                ccache, compiler, tasking_version, for_machine, is_cross, info,
                full_version=full_version, linker=linker)

    _handle_exceptions(popen_exceptions, compilers)
    raise EnvironmentException(f'Unknown compiler {compilers}')

def detect_c_compiler(env: 'Environment', for_machine: MachineChoice) -> Compiler:
    return _detect_c_or_cpp_compiler(env, 'c', for_machine)

def detect_cpp_compiler(env: 'Environment', for_machine: MachineChoice) -> Compiler:
    return _detect_c_or_cpp_compiler(env, 'cpp', for_machine)

def detect_cuda_compiler(env: 'Environment', for_machine: MachineChoice) -> Compiler:
    from .cuda import CudaCompiler, Phase
    from ..options import OptionKey
    from ..linkers.linkers import CudaLinker
    popen_exceptions = {}
    is_cross = env.is_cross_build(for_machine)
    compilers, ccache = _get_compilers(env, 'cuda', for_machine)
    info = env.machines[for_machine]
    for compiler in compilers:
        arg = '--version'
        try:
            p, out, err = Popen_safe_logged(compiler + [arg], msg='Detecting compiler via')
        except OSError as e:
            popen_exceptions[join_args(compiler + [arg])] = e
            continue
        # Example nvcc printout:
        #
        #     nvcc: NVIDIA (R) Cuda compiler driver
        #     Copyright (c) 2005-2018 NVIDIA Corporation
        #     Built on Sat_Aug_25_21:08:01_CDT_2018
        #     Cuda compilation tools, release 10.0, V10.0.130
        #
        # search_version() first finds the "10.0" after "release",
        # rather than the more precise "10.0.130" after "V".
        # The patch version number is occasionally important; For
        # instance, on Linux,
        #    - CUDA Toolkit 8.0.44 requires NVIDIA Driver 367.48
        #    - CUDA Toolkit 8.0.61 requires NVIDIA Driver 375.26
        # Luckily, the "V" also makes it very simple to extract
        # the full version:
        version = out.strip().rsplit('V', maxsplit=1)[-1]
        cpp_compiler = detect_cpp_compiler(env, for_machine)
        cls = CudaCompiler
        env.coredata.add_lang_args(cls.language, cls, for_machine, env)
        key = OptionKey('cuda_link_args', machine=for_machine)
        if key in env.options:
            # To fix LDFLAGS issue
            val = env.options[key]
            assert isinstance(val, list)
            env.coredata.set_options({key: cls.to_host_flags_base(val, Phase.LINKER)})
        linker = CudaLinker(compiler, for_machine, CudaCompiler.LINKER_PREFIX, [], version=CudaLinker.parse_version())
        return cls(ccache, compiler, version, for_machine, is_cross, host_compiler=cpp_compiler, info=info, linker=linker)
    raise EnvironmentException(f'Could not find suitable CUDA compiler: "{"; ".join([" ".join(c) for c in compilers])}"')

def detect_fortran_compiler(env: 'Environment', for_machine: MachineChoice) -> Compiler:
    from . import fortran
    from ..linkers import linkers
    popen_exceptions: T.Dict[str, T.Union[Exception, str]] = {}
    compilers, ccache = _get_compilers(env, 'fortran', for_machine)
    is_cross = env.is_cross_build(for_machine)
    info = env.machines[for_machine]
    cls: T.Type[FortranCompiler]
    for compiler in compilers:
        # capture help text for possible fallback
        try:
            _, help_out, _ = Popen_safe_logged(compiler + ['--help'], msg='Detecting compiler via')
        except OSError as e:
            popen_exceptions[join_args(compiler + ['--help'])] = e
            help_out = ''

        for arg in ['--version', '-V']:
            try:
                p, out, err = Popen_safe_logged(compiler + [arg], msg='Detecting compiler via')
            except OSError as e:
                popen_exceptions[join_args(compiler + [arg])] = e
                continue

            version = search_version(out)
            full_version = out.split('\n', 1)[0]

            guess_gcc_or_lcc: T.Optional[str] = None
            if 'GNU Fortran' in out:
                guess_gcc_or_lcc = 'gcc'
            if 'e2k' in out and 'lcc' in out:
                guess_gcc_or_lcc = 'lcc'

            if guess_gcc_or_lcc:
                defines = _get_gnu_compiler_defines(compiler, 'fortran')
                if not defines:
                    popen_exceptions[join_args(compiler)] = 'no pre-processor defines'
                    continue
                if guess_gcc_or_lcc == 'lcc':
                    version = _get_lcc_version_from_defines(defines)
                    cls = fortran.ElbrusFortranCompiler
                    linker = guess_nix_linker(env, compiler, cls, version, for_machine)
                    return cls(
                        compiler, version, for_machine, is_cross, info,
                        defines, full_version=full_version, linker=linker)
                else:
                    version = _get_gnu_version_from_defines(defines)
                    cls = fortran.GnuFortranCompiler
                    linker = guess_nix_linker(env, compiler, cls, version, for_machine)
                    return cls(
                        compiler, version, for_machine, is_cross, info,
                        defines, full_version=full_version, linker=linker)

            if 'Arm C/C++/Fortran Compiler' in out:
                cls = fortran.ArmLtdFlangFortranCompiler
                arm_ver_match = re.search(r'version (\d+)\.(\d+)\.?(\d+)? \(build number (\d+)\)', out)
                assert arm_ver_match is not None, 'for mypy'  # because mypy *should* be complaining that this could be None
                version = '.'.join([x for x in arm_ver_match.groups() if x is not None])
                linker = guess_nix_linker(env, compiler, cls, version, for_machine)
                return cls(
                    compiler, version, for_machine, is_cross, info,
                    linker=linker)
            if 'G95' in out:
                cls = fortran.G95FortranCompiler
                linker = guess_nix_linker(env, compiler, cls, version, for_machine)
                return cls(
                    compiler, version, for_machine, is_cross, info,
                    full_version=full_version, linker=linker)

            if 'Sun Fortran' in err:
                version = search_version(err)
                cls = fortran.SunFortranCompiler
                linker = guess_nix_linker(env, compiler, cls, version, for_machine)
                return cls(
                    compiler, version, for_machine, is_cross, info,
                    full_version=full_version, linker=linker)

            if 'Intel(R) Fortran Compiler for applications' in err:
                version = search_version(err)
                target = 'x86' if 'IA-32' in err else 'x86_64'
                cls = fortran.IntelLLVMClFortranCompiler
                env.coredata.add_lang_args(cls.language, cls, for_machine, env)
                linker = linkers.XilinkDynamicLinker(for_machine, [], version=version)
                return cls(
                    compiler, version, for_machine, is_cross, info,
                    target, linker=linker)

            if 'Intel(R) Visual Fortran' in err or 'Intel(R) Fortran' in err:
                version = search_version(err)
                target = 'x86' if 'IA-32' in err else 'x86_64'
                cls = fortran.IntelClFortranCompiler
                env.coredata.add_lang_args(cls.language, cls, for_machine, env)
                linker = linkers.XilinkDynamicLinker(for_machine, [], version=version)
                return cls(
                    compiler, version, for_machine, is_cross, info,
                    target, linker=linker)

            if 'ifort (IFORT)' in out:
                cls = fortran.IntelFortranCompiler
                linker = guess_nix_linker(env, compiler, cls, version, for_machine)
                return cls(
                    compiler, version, for_machine, is_cross, info,
                    full_version=full_version, linker=linker)

            if 'ifx (IFORT)' in out or 'ifx (IFX)' in out:
                cls = fortran.IntelLLVMFortranCompiler
                linker = guess_nix_linker(env, compiler, cls, version, for_machine)
                return cls(
                    compiler, version, for_machine, is_cross, info,
                    full_version=full_version, linker=linker)

            if 'PathScale EKOPath(tm)' in err:
                return fortran.PathScaleFortranCompiler(
                    compiler, version, for_machine, is_cross, info,
                    full_version=full_version)

            if 'PGI Compilers' in out:
                cls = fortran.PGIFortranCompiler
                env.coredata.add_lang_args(cls.language, cls, for_machine, env)
                linker = linkers.PGIDynamicLinker(compiler, for_machine,
                                                  cls.LINKER_PREFIX, [], version=version)
                return cls(
                    compiler, version, for_machine, is_cross, info,
                    full_version=full_version, linker=linker)

            if 'NVIDIA Compilers and Tools' in out:
                cls = fortran.NvidiaHPC_FortranCompiler
                env.coredata.add_lang_args(cls.language, cls, for_machine, env)
                linker = linkers.PGIDynamicLinker(compiler, for_machine,
                                                  cls.LINKER_PREFIX, [], version=version)
                return cls(
                    compiler, version, for_machine, is_cross, info,
                    full_version=full_version, linker=linker)

            def _get_linker_try_windows(cls: T.Type['Compiler']) -> T.Optional['DynamicLinker']:
                linker = None
                if 'windows' in out or env.machines[for_machine].is_windows():
                    # If we're in a MINGW context this actually will use a gnu
                    # style ld, but for flang on "real" windows we'll use
                    # either link.exe or lld-link.exe
                    try:
                        linker = guess_win_linker(
                            env, compiler, cls, version,
                            for_machine, invoked_directly=False
                        )
                    except MesonException:
                        pass
                if linker is None:
                    linker = guess_nix_linker(env, compiler, cls,
                                              version, for_machine)
                return linker

            if 'flang-new' in out or 'flang LLVM compiler' in help_out:
                cls = fortran.LlvmFlangFortranCompiler
                linker = _get_linker_try_windows(cls)
                return cls(
                    compiler, version, for_machine, is_cross, info,
                    full_version=full_version, linker=linker)

            if 'flang' in out or 'clang' in out:
                cls = fortran.ClassicFlangFortranCompiler
                linker = _get_linker_try_windows(cls)
                return cls(
                    compiler, version, for_machine, is_cross, info,
                    full_version=full_version, linker=linker)

            if 'Open64 Compiler Suite' in err:
                cls = fortran.Open64FortranCompiler
                linker = guess_nix_linker(env,
                                          compiler, cls, version, for_machine)
                return cls(
                    compiler, version, for_machine, is_cross, info,
                    full_version=full_version, linker=linker)

            if 'NAG Fortran' in err:
                full_version = err.split('\n', 1)[0]
                version = full_version.split()[-1]
                cls = fortran.NAGFortranCompiler
                env.coredata.add_lang_args(cls.language, cls, for_machine, env)
                linker = linkers.NAGDynamicLinker(
                    compiler, for_machine, cls.LINKER_PREFIX, [],
                    version=version)
                return cls(
                    compiler, version, for_machine, is_cross, info,
                    full_version=full_version, linker=linker)

    _handle_exceptions(popen_exceptions, compilers)
    raise EnvironmentException('Unreachable code (exception to make mypy happy)')

def detect_objc_compiler(env: 'Environment', for_machine: MachineChoice) -> 'Compiler':
    return _detect_objc_or_objcpp_compiler(env, 'objc', for_machine)

def detect_objcpp_compiler(env: 'Environment', for_machine: MachineChoice) -> 'Compiler':
    return _detect_objc_or_objcpp_compiler(env, 'objcpp', for_machine)

def _detect_objc_or_objcpp_compiler(env: 'Environment', lang: str, for_machine: MachineChoice) -> 'Compiler':
    from . import objc, objcpp
    popen_exceptions: T.Dict[str, T.Union[Exception, str]] = {}
    compilers, ccache = _get_compilers(env, lang, for_machine)
    is_cross = env.is_cross_build(for_machine)
    info = env.machines[for_machine]
    comp: T.Union[T.Type[objc.ObjCCompiler], T.Type[objcpp.ObjCPPCompiler]]

    for compiler in compilers:
        arg = ['--version']
        try:
            p, out, err = Popen_safe_logged(compiler + arg, msg='Detecting compiler via')
        except OSError as e:
            popen_exceptions[join_args(compiler + arg)] = e
            continue
        version = search_version(out)
        if 'Free Software Foundation' in out:
            defines = _get_gnu_compiler_defines(compiler, lang)
            if not defines:
                popen_exceptions[join_args(compiler)] = 'no pre-processor defines'
                continue
            version = _get_gnu_version_from_defines(defines)
            comp = objc.GnuObjCCompiler if lang == 'objc' else objcpp.GnuObjCPPCompiler
            linker = guess_nix_linker(env, compiler, comp, version, for_machine)
            return comp(
                ccache, compiler, version, for_machine, is_cross, info,
                defines, linker=linker)
        if 'clang' in out:
            linker = None
            defines = _get_clang_compiler_defines(compiler, lang)
            if not defines:
                popen_exceptions[join_args(compiler)] = 'no pre-processor defines'
                continue
            if 'Apple' in out:
                comp = objc.AppleClangObjCCompiler if lang == 'objc' else objcpp.AppleClangObjCPPCompiler
            else:
                comp = objc.ClangObjCCompiler if lang == 'objc' else objcpp.ClangObjCPPCompiler
            if 'windows' in out or env.machines[for_machine].is_windows():
                # If we're in a MINGW context this actually will use a gnu style ld
                try:
                    linker = guess_win_linker(env, compiler, comp, version, for_machine)
                except MesonException:
                    pass

            if not linker:
                linker = guess_nix_linker(env, compiler, comp, version, for_machine)
            return comp(
                ccache, compiler, version, for_machine,
                is_cross, info, linker=linker, defines=defines)
    _handle_exceptions(popen_exceptions, compilers)
    raise EnvironmentException('Unreachable code (exception to make mypy happy)')

def detect_java_compiler(env: 'Environment', for_machine: MachineChoice) -> Compiler:
    from .java import JavaCompiler
    exelist = env.lookup_binary_entry(for_machine, 'java')
    info = env.machines[for_machine]
    if exelist is None:
        # TODO support fallback
        exelist = [defaults['java'][0]]

    try:
        p, out, err = Popen_safe_logged(exelist + ['-version'], msg='Detecting compiler via')
    except OSError:
        raise EnvironmentException('Could not execute Java compiler: {}'.format(join_args(exelist)))
    if 'javac' in out or 'javac' in err:
        version = search_version(err if 'javac' in err else out)
        if not version or version == 'unknown version':
            parts = (err if 'javac' in err else out).split()
            if len(parts) > 1:
                version = parts[1]
        comp_class = JavaCompiler
        env.coredata.add_lang_args(comp_class.language, comp_class, for_machine, env)
        return comp_class(exelist, version, for_machine, info)
    raise EnvironmentException('Unknown compiler: ' + join_args(exelist))

def detect_cs_compiler(env: 'Environment', for_machine: MachineChoice) -> Compiler:
    from . import cs
    compilers, ccache = _get_compilers(env, 'cs', for_machine)
    popen_exceptions = {}
    info = env.machines[for_machine]
    for comp in compilers:
        try:
            p, out, err = Popen_safe_logged(comp + ['--version'], msg='Detecting compiler via')
        except OSError as e:
            popen_exceptions[join_args(comp + ['--version'])] = e
            continue

        version = search_version(out)
        cls: T.Type[cs.CsCompiler]
        if 'Mono' in out:
            cls = cs.MonoCompiler
        elif "Visual C#" in out:
            cls = cs.VisualStudioCsCompiler
        else:
            continue
        env.coredata.add_lang_args(cls.language, cls, for_machine, env)
        return cls(comp, version, for_machine, info)

    _handle_exceptions(popen_exceptions, compilers)
    raise EnvironmentException('Unreachable code (exception to make mypy happy)')

def detect_cython_compiler(env: 'Environment', for_machine: MachineChoice) -> Compiler:
    """Search for a cython compiler."""
    from .cython import CythonCompiler
    compilers, _ = _get_compilers(env, 'cython', MachineChoice.BUILD)
    is_cross = env.is_cross_build(for_machine)
    info = env.machines[for_machine]

    popen_exceptions: T.Dict[str, Exception] = {}
    for comp in compilers:
        try:
            _, out, err = Popen_safe_logged(comp + ['-V'], msg='Detecting compiler via')
        except OSError as e:
            popen_exceptions[join_args(comp + ['-V'])] = e
            continue

        version: T.Optional[str] = None
        # 3.0
        if 'Cython' in out:
            version = search_version(out)
        # older
        elif 'Cython' in err:
            version = search_version(err)
        if version is not None:
            comp_class = CythonCompiler
            env.coredata.add_lang_args(comp_class.language, comp_class, for_machine, env)
            return comp_class([], comp, version, for_machine, info, is_cross=is_cross)
    _handle_exceptions(popen_exceptions, compilers)
    raise EnvironmentException('Unreachable code (exception to make mypy happy)')

def detect_vala_compiler(env: 'Environment', for_machine: MachineChoice) -> Compiler:
    from .vala import ValaCompiler
    exelist = env.lookup_binary_entry(MachineChoice.BUILD, 'vala')
    is_cross = env.is_cross_build(for_machine)
    info = env.machines[for_machine]
    if exelist is None:
        # TODO support fallback
        exelist = [defaults['vala'][0]]

    try:
        p, out = Popen_safe_logged(exelist + ['--version'], msg='Detecting compiler via')[0:2]
    except OSError:
        raise EnvironmentException('Could not execute Vala compiler: {}'.format(join_args(exelist)))
    version = search_version(out)
    if 'Vala' in out:
        comp_class = ValaCompiler
        env.coredata.add_lang_args(comp_class.language, comp_class, for_machine, env)
        return comp_class(exelist, version, for_machine, is_cross, info)
    raise EnvironmentException('Unknown compiler: ' + join_args(exelist))

def detect_rust_compiler(env: 'Environment', for_machine: MachineChoice) -> RustCompiler:
    from . import rust
    from ..linkers import linkers
    popen_exceptions: T.Dict[str, Exception] = {}
    compilers, _ = _get_compilers(env, 'rust', for_machine)
    is_cross = env.is_cross_build(for_machine)
    info = env.machines[for_machine]

    cc = detect_c_compiler(env, for_machine)
    is_link_exe = isinstance(cc.linker, linkers.VisualStudioLikeLinkerMixin)
    override = env.lookup_binary_entry(for_machine, 'rust_ld')

    for compiler in compilers:
        arg = ['--version']
        try:
            out = Popen_safe_logged(compiler + arg, msg='Detecting compiler via')[1]
        except OSError as e:
            popen_exceptions[join_args(compiler + arg)] = e
            continue

        version = search_version(out)
        cls: T.Type[RustCompiler] = rust.RustCompiler

        # Clippy is a wrapper around rustc, but it doesn't have rustc in its
        # output. We can otherwise treat it as rustc.
        if 'clippy' in out:
            # clippy returns its own version and not the rustc version by
            # default so try harder here to get the correct version.
            # Also replace the whole output with the rustc output in
            # case this is later used for other purposes.
            arg = ['--rustc', '--version']
            try:
                out = Popen_safe(compiler + arg)[1]
            except OSError as e:
                popen_exceptions[join_args(compiler + arg)] = e
                continue
            version = search_version(out)

            cls = rust.ClippyRustCompiler
            mlog.deprecation(
                'clippy-driver is not intended as a general purpose compiler. '
                'You can use "ninja clippy" in order to run clippy on a '
                'meson project.')

        if 'rustc' in out:
            # On Linux and mac rustc will invoke gcc (clang for mac
            # presumably) and it can do this windows, for dynamic linking.
            # this means the easiest way to C compiler for dynamic linking.
            # figure out what linker to use is to just get the value of the
            # C compiler and use that as the basis of the rust linker.
            # However, there are two things we need to change, if CC is not
            # the default use that, and second add the necessary arguments
            # to rust to use -fuse-ld

            if any(a.startswith('linker=') for a in compiler):
                mlog.warning(
                    'Please do not put -C linker= in your compiler '
                    'command, set rust_ld=command in your cross file '
                    'or use the RUSTC_LD environment variable, otherwise meson '
                    'will override your selection.')

            compiler = compiler.copy()  # avoid mutating the original list

            if override is None:
                extra_args: T.Dict[str, T.Union[str, bool]] = {}
                always_args: T.List[str] = []
                if is_link_exe:
                    compiler.extend(cls.use_linker_args(cc.linker.exelist[0], ''))
                    extra_args['direct'] = True
                    extra_args['machine'] = cc.linker.machine
                else:
                    exelist = cc.linker.exelist + cc.linker.get_always_args()
                    if os.path.basename(exelist[0]) in {'ccache', 'sccache'}:
                        del exelist[0]
                    c = exelist.pop(0)
                    compiler.extend(cls.use_linker_args(c, ''))

                    # Also ensure that we pass any extra arguments to the linker
                    for l in exelist:
                        compiler.extend(['-C', f'link-arg={l}'])

                # This trickery with type() gets us the class of the linker
                # so we can initialize a new copy for the Rust Compiler
                # TODO rewrite this without type: ignore
                assert cc.linker is not None, 'for mypy'
                if is_link_exe:
                    linker = type(cc.linker)(for_machine, always_args, exelist=cc.linker.exelist,   # type: ignore
                                             version=cc.linker.version, **extra_args)               # type: ignore
                else:
                    linker = type(cc.linker)(compiler, for_machine, cc.LINKER_PREFIX,
                                             always_args=always_args, version=cc.linker.version,
                                             **extra_args)
            elif 'link' in override[0]:
                linker = guess_win_linker(env,
                                          override, cls, version, for_machine, use_linker_prefix=False)
                # rustc takes linker arguments without a prefix, and
                # inserts the correct prefix itself.
                assert isinstance(linker, linkers.VisualStudioLikeLinkerMixin)
                linker.direct = True
                compiler.extend(cls.use_linker_args(linker.exelist[0], ''))
            else:
                # On linux and macos rust will invoke the c compiler for
                # linking, on windows it will use lld-link or link.exe.
                # we will simply ask for the C compiler that corresponds to
                # it, and use that.
                cc = _detect_c_or_cpp_compiler(env, 'c', for_machine, override_compiler=override)
                linker = cc.linker

                # Of course, we're not going to use any of that, we just
                # need it to get the proper arguments to pass to rustc
                c = linker.exelist[1] if linker.exelist[0].endswith('ccache') else linker.exelist[0]
                compiler.extend(cls.use_linker_args(c, ''))

            env.coredata.add_lang_args(cls.language, cls, for_machine, env)
            return cls(
                compiler, version, for_machine, is_cross, info,
                linker=linker)

    _handle_exceptions(popen_exceptions, compilers)
    raise EnvironmentException('Unreachable code (exception to make mypy happy)')

def detect_d_compiler(env: 'Environment', for_machine: MachineChoice) -> Compiler:
    from . import c, d
    info = env.machines[for_machine]

    # Detect the target architecture, required for proper architecture handling on Windows.
    # MSVC compiler is required for correct platform detection.
    c_compiler = {'c': detect_c_compiler(env, for_machine)}
    is_msvc = isinstance(c_compiler['c'], c.VisualStudioCCompiler)
    if not is_msvc:
        c_compiler = {}

    # Import here to avoid circular imports
    from ..environment import detect_cpu_family
    arch = detect_cpu_family(c_compiler)
    if is_msvc and arch == 'x86':
        arch = 'x86_mscoff'

    popen_exceptions = {}
    is_cross = env.is_cross_build(for_machine)
    compilers, ccache = _get_compilers(env, 'd', for_machine)
    cls: T.Type[d.DCompiler]
    for exelist in compilers:
        # Search for a D compiler.
        # We prefer LDC over GDC unless overridden with the DC
        # environment variable because LDC has a much more
        # up to date language version at time (2016).
        if os.path.basename(exelist[-1]).startswith(('ldmd', 'gdmd')):
            raise EnvironmentException(
                f'Meson does not support {exelist[-1]} as it is only a DMD frontend for another compiler.'
                'Please provide a valid value for DC or unset it so that Meson can resolve the compiler by itself.')
        try:
            p, out = Popen_safe(exelist + ['--version'])[0:2]
        except OSError as e:
            popen_exceptions[join_args(exelist + ['--version'])] = e
            continue
        version = search_version(out)
        full_version = out.split('\n', 1)[0]

        if 'LLVM D compiler' in out:
            cls = d.LLVMDCompiler
            # LDC seems to require a file
            # We cannot use NamedTemporaryFile on windows, its documented
            # to not work for our uses. So, just use mkstemp and only have
            # one path for simplicity.
            o, f = tempfile.mkstemp('.d')
            os.close(o)

            try:
                if info.is_windows() or info.is_cygwin():
                    objfile = os.path.basename(f)[:-1] + 'obj'
                    extra_args = [f]
                    if is_cross:
                        extra_args.append(f'-mtriple={info.cpu}-windows')

                    linker = guess_win_linker(env,
                                              exelist,
                                              cls, full_version, for_machine,
                                              use_linker_prefix=True, invoked_directly=False,
                                              extra_args=extra_args)
                else:
                    # LDC writes an object file to the current working directory.
                    # Clean it up.
                    objfile = os.path.basename(f)[:-1] + 'o'
                    linker = guess_nix_linker(env,
                                              exelist, cls, full_version, for_machine,
                                              extra_args=[f])
            finally:
                windows_proof_rm(f)
                windows_proof_rm(objfile)

            return cls(
                exelist, version, for_machine, info, arch,
                full_version=full_version, linker=linker,
                is_cross=is_cross, version_output=out)
        elif 'gdc' in out:
            cls = d.GnuDCompiler
            linker = guess_nix_linker(env, exelist, cls, version, for_machine)
            return cls(
                exelist, version, for_machine, info, arch,
                is_cross=is_cross, full_version=full_version, linker=linker)
        elif 'The D Language Foundation' in out or 'Digital Mars' in out:
            cls = d.DmdDCompiler
            # DMD seems to require a file
            # We cannot use NamedTemporaryFile on windows, its documented
            # to not work for our uses. So, just use mkstemp and only have
            # one path for simplicity.
            o, f = tempfile.mkstemp('.d')
            os.close(o)

            # DMD as different detection logic for x86 and x86_64
            arch_arg = '-m64' if arch == 'x86_64' else '-m32'

            try:
                if info.is_windows() or info.is_cygwin():
                    objfile = os.path.basename(f)[:-1] + 'obj'
                    linker = guess_win_linker(env,
                                              exelist, cls, full_version, for_machine,
                                              invoked_directly=False, extra_args=[f, arch_arg])
                else:
                    objfile = os.path.basename(f)[:-1] + 'o'
                    linker = guess_nix_linker(env,
                                              exelist, cls, full_version, for_machine,
                                              extra_args=[f, arch_arg])
            finally:
                windows_proof_rm(f)
                windows_proof_rm(objfile)

            return cls(
                exelist, version, for_machine, info, arch,
                full_version=full_version, linker=linker)
        raise EnvironmentException('Unknown compiler: ' + join_args(exelist))

    _handle_exceptions(popen_exceptions, compilers)
    raise EnvironmentException('Unreachable code (exception to make mypy happy)')

def detect_swift_compiler(env: 'Environment', for_machine: MachineChoice) -> Compiler:
    from .swift import SwiftCompiler
    exelist = env.lookup_binary_entry(for_machine, 'swift')
    is_cross = env.is_cross_build(for_machine)
    info = env.machines[for_machine]
    if exelist is None:
        # TODO support fallback
        exelist = [defaults['swift'][0]]

    try:
        p, _, err = Popen_safe_logged(exelist + ['-v'], msg='Detecting compiler via')
    except OSError:
        raise EnvironmentException('Could not execute Swift compiler: {}'.format(join_args(exelist)))
    version = search_version(err)
    if 'Swift' in err:
        # As for 5.0.1 swiftc *requires* a file to check the linker:
        with tempfile.NamedTemporaryFile(suffix='.swift') as f:
            cls = SwiftCompiler
            linker = guess_nix_linker(env,
                                      exelist, cls, version, for_machine,
                                      extra_args=[f.name, '-o', '/dev/null'])
        return cls(
            exelist, version, for_machine, is_cross, info, linker=linker)

    raise EnvironmentException('Unknown compiler: ' + join_args(exelist))

def detect_nasm_compiler(env: 'Environment', for_machine: MachineChoice) -> Compiler:
    from .asm import NasmCompiler, YasmCompiler, MetrowerksAsmCompilerARM, MetrowerksAsmCompilerEmbeddedPowerPC
    is_cross = env.is_cross_build(for_machine)

    # When cross compiling and nasm is not defined in the cross file we can
    # fallback to the build machine nasm.
    compilers, _ = _get_compilers(env, 'nasm', for_machine, allow_build_machine=True)

    # We need a C compiler to properly detect the machine info and linker
    cc = detect_c_compiler(env, for_machine)
    if not is_cross:
        from ..environment import detect_machine_info
        info = detect_machine_info({'c': cc})
    else:
        info = env.machines[for_machine]

    popen_exceptions: T.Dict[str, Exception] = {}
    for comp in compilers:
        if comp == ['nasm'] and is_windows() and not shutil.which(comp[0]):
            # nasm is not in PATH on Windows by default
            default_path = os.path.join(os.environ['ProgramFiles'], 'NASM')
            comp[0] = shutil.which(comp[0], path=default_path) or comp[0]
        try:
            output = Popen_safe_logged(comp + ['--version'], msg='Detecting compiler via')[1]
        except OSError as e:
            popen_exceptions[' '.join(comp + ['--version'])] = e
            continue

        version = search_version(output)
        if 'NASM' in output:
            comp_class = NasmCompiler
            env.coredata.add_lang_args(comp_class.language, comp_class, for_machine, env)
            return comp_class([], comp, version, for_machine, info, cc.linker, is_cross=is_cross)
        elif 'yasm' in output:
            comp_class = YasmCompiler
            env.coredata.add_lang_args(comp_class.language, comp_class, for_machine, env)
            return comp_class([], comp, version, for_machine, info, cc.linker, is_cross=is_cross)
        elif 'Metrowerks' in output or 'Freescale' in output:
            if 'ARM' in output:
                comp_class_mwasmarm = MetrowerksAsmCompilerARM
                env.coredata.add_lang_args(comp_class_mwasmarm.language, comp_class_mwasmarm, for_machine, env)
                return comp_class_mwasmarm([], comp, version, for_machine, info, cc.linker, is_cross=is_cross)
            else:
                comp_class_mwasmeppc = MetrowerksAsmCompilerEmbeddedPowerPC
                env.coredata.add_lang_args(comp_class_mwasmeppc.language, comp_class_mwasmeppc, for_machine, env)
                return comp_class_mwasmeppc([], comp, version, for_machine, info, cc.linker, is_cross=is_cross)

    _handle_exceptions(popen_exceptions, compilers)
    raise EnvironmentException('Unreachable code (exception to make mypy happy)')

def detect_masm_compiler(env: 'Environment', for_machine: MachineChoice) -> Compiler:
    # We need a C compiler to properly detect the machine info and linker
    is_cross = env.is_cross_build(for_machine)
    cc = detect_c_compiler(env, for_machine)
    if not is_cross:
        from ..environment import detect_machine_info
        info = detect_machine_info({'c': cc})
    else:
        info = env.machines[for_machine]

    from .asm import MasmCompiler, MasmARMCompiler
    comp_class: T.Type[Compiler]
    if info.cpu_family == 'x86':
        comp = ['ml']
        comp_class = MasmCompiler
        arg = '/?'
    elif info.cpu_family == 'x86_64':
        comp = ['ml64']
        comp_class = MasmCompiler
        arg = '/?'
    elif info.cpu_family == 'arm':
        comp = ['armasm']
        comp_class = MasmARMCompiler
        arg = '-h'
    elif info.cpu_family == 'aarch64':
        comp = ['armasm64']
        comp_class = MasmARMCompiler
        arg = '-h'
    else:
        raise EnvironmentException(f'Platform {info.cpu_family} not supported by MASM')

    popen_exceptions: T.Dict[str, Exception] = {}
    try:
        output = Popen_safe(comp + [arg])[2]
        version = search_version(output)
        env.coredata.add_lang_args(comp_class.language, comp_class, for_machine, env)
        return comp_class([], comp, version, for_machine, info, cc.linker, is_cross=is_cross)
    except OSError as e:
        popen_exceptions[' '.join(comp + [arg])] = e
    _handle_exceptions(popen_exceptions, [comp])
    raise EnvironmentException('Unreachable code (exception to make mypy happy)')

def detect_linearasm_compiler(env: Environment, for_machine: MachineChoice) -> Compiler:
    from .asm import TILinearAsmCompiler
    comp = ['cl6x']
    comp_class: T.Type[Compiler] = TILinearAsmCompiler
    arg = '-h'
    info = env.machines[for_machine]
    cc = detect_c_compiler(env, for_machine)
    is_cross = env.is_cross_build(for_machine)

    popen_exceptions: T.Dict[str, Exception] = {}
    try:
        output = Popen_safe(comp + [arg])[2]
        version = search_version(output)
        env.coredata.add_lang_args(comp_class.language, comp_class, for_machine, env)
        return comp_class([], comp, version, for_machine, info, cc.linker, is_cross=is_cross)
    except OSError as e:
        popen_exceptions[' '.join(comp + [arg])] = e
    _handle_exceptions(popen_exceptions, [comp])
    raise EnvironmentException('Unreachable code (exception to make mypy happy)')

# GNU/Clang defines and version
# =============================

def _get_gnu_compiler_defines(compiler: T.List[str], lang: str) -> T.Dict[str, str]:
    """
    Get the list of GCC pre-processor defines
    """
    from .mixins.gnu import gnu_lang_map

    def _try_obtain_compiler_defines(args: T.List[str]) -> str:
        mlog.debug(f'Running command: {join_args(args)}')
        p, output, error = Popen_safe(compiler + args, write='', stdin=subprocess.PIPE)
        if p.returncode != 0:
            raise EnvironmentException('Unable to get gcc pre-processor defines:\n'
                                       f'Compiler stdout:\n{output}\n-----\n'
                                       f'Compiler stderr:\n{error}\n-----\n')
        return output

    # Arguments to output compiler pre-processor defines to stdout
    # gcc, g++, and gfortran all support these arguments
    baseline_test_args = ['-E', '-dM', '-']
    try:
        # We assume that when _get_gnu_compiler_defines is called, it's
        # close enough to a GCCish compiler so we reuse the _LANG_MAP
        # from the GCC mixin. This isn't a dangerous assumption because
        # we fallback if the detection fails anyway.

        # We might not have a match for Fortran, so fallback to detection
        # based on the driver.
        lang = gnu_lang_map[lang]

        # The compiler may not infer the target language based on the driver name
        # so first, try with '-cpp -x lang', then fallback without given it's less
        # portable. We try with '-cpp' as GCC needs it for Fortran at least, and
        # it seems to do no harm.
        output = _try_obtain_compiler_defines(['-cpp', '-x', lang] + baseline_test_args)
    except (EnvironmentException, KeyError):
        mlog.debug(f'pre-processor extraction using -cpp -x {lang} failed, falling back w/o lang')
        output = _try_obtain_compiler_defines(baseline_test_args)

    # Parse several lines of the type:
    # `#define ___SOME_DEF some_value`
    # and extract `___SOME_DEF`
    defines: T.Dict[str, str] = {}
    for line in output.split('\n'):
        if not line:
            continue
        d, *rest = line.split(' ', 2)
        if d != '#define':
            continue
        if len(rest) == 1:
            defines[rest[0]] = ''
        if len(rest) == 2:
            defines[rest[0]] = rest[1]
    return defines

def _get_clang_compiler_defines(compiler: T.List[str], lang: str) -> T.Dict[str, str]:
    """
    Get the list of Clang pre-processor defines
    """
    from .mixins.clang import clang_lang_map

    def _try_obtain_compiler_defines(args: T.List[str]) -> str:
        mlog.debug(f'Running command: {join_args(args)}')
        p, output, error = Popen_safe(compiler + args, write='', stdin=subprocess.PIPE)
        if p.returncode != 0:
            raise EnvironmentException('Unable to get clang pre-processor defines:\n'
                                       f'Compiler stdout:\n{output}\n-----\n'
                                       f'Compiler stderr:\n{error}\n-----\n')
        return output

    # Arguments to output compiler pre-processor defines to stdout
    baseline_test_args = ['-E', '-dM', '-']
    try:
        # We assume that when _get_clang_compiler_defines is called, it's
        # close enough to a Clangish compiler so we reuse the _LANG_MAP
        # from the Clang mixin. This isn't a dangerous assumption because
        # we fallback if the detection fails anyway.

        # We might not have a match for Fortran, so fallback to detection
        # based on the driver.
        lang = clang_lang_map[lang]

        # The compiler may not infer the target language based on the driver name.
        # Try first with '-x lang' to supported systemwide language level overrides,
        # then fallback to without since it's a more recent option.
        output = _try_obtain_compiler_defines(['-x', lang] + baseline_test_args)
    except (EnvironmentException, KeyError):
        mlog.debug(f'pre-processor extraction using -x {lang} failed, falling back w/o lang')
        output = _try_obtain_compiler_defines(baseline_test_args)

    defines: T.Dict[str, str] = {}
    for line in output.split('\n'):
        if not line:
            continue
        d, *rest = line.split(' ', 2)
        if d != '#define':
            continue
        if len(rest) == 1:
            defines[rest[0]] = ''
        if len(rest) == 2:
            defines[rest[0]] = rest[1]
    return defines

def _get_gnu_version_from_defines(defines: T.Dict[str, str]) -> str:
    dot = '.'
    major = defines.get('__GNUC__', '0')
    minor = defines.get('__GNUC_MINOR__', '0')
    patch = defines.get('__GNUC_PATCHLEVEL__', '0')
    return dot.join((major, minor, patch))

def _get_lcc_version_from_defines(defines: T.Dict[str, str]) -> str:
    dot = '.'
    generation_and_major = defines.get('__LCC__', '100')
    generation = generation_and_major[:1]
    major = generation_and_major[1:]
    minor = defines.get('__LCC_MINOR__', '0')
    return dot.join((generation, major, minor))
