# SPDX-License-Identifier: Apache-2.0
# Copyright 2012-2020 The Meson development team
# Copyright © 2023-2025 Intel Corporation

from __future__ import annotations

import os
import shutil
import typing as T

from . import coredata
from . import mesonlib
from . import mlog
from .mesonlib import MachineChoice, Popen_safe, search_version, quote_arg, split_args
from .programs import ExternalProgram


def detect_gcovr(gcovr_exe: str = 'gcovr', min_version: str = '3.3', log: bool = False) \
        -> T.Union[T.Tuple[None, None], T.Tuple[str, str]]:
    try:
        p, found = Popen_safe([gcovr_exe, '--version'])[0:2]
    except (FileNotFoundError, PermissionError):
        # Doesn't exist in PATH or isn't executable
        return None, None
    found = search_version(found)
    if p.returncode == 0 and mesonlib.version_compare(found, '>=' + min_version):
        if log:
            mlog.log('Found gcovr-{} at {}'.format(found, quote_arg(shutil.which(gcovr_exe))))
        return gcovr_exe, found
    return None, None

def detect_lcov(lcov_exe: str = 'lcov', log: bool = False) \
        -> T.Union[T.Tuple[None, None], T.Tuple[str, str]]:
    try:
        p, found = Popen_safe([lcov_exe, '--version'])[0:2]
    except (FileNotFoundError, PermissionError):
        # Doesn't exist in PATH or isn't executable
        return None, None
    found = search_version(found)
    if p.returncode == 0 and found:
        if log:
            mlog.log('Found lcov-{} at {}'.format(found, quote_arg(shutil.which(lcov_exe))))
        return lcov_exe, found
    return None, None

def detect_llvm_cov(suffix: T.Optional[str] = None) -> T.Optional[str]:
    # If there's a known suffix or forced lack of suffix, use that
    if suffix is not None:
        if suffix == '':
            tool = 'llvm-cov'
        else:
            tool = f'llvm-cov-{suffix}'
        if shutil.which(tool) is not None:
            return tool
    else:
        # Otherwise guess in the dark
        tools = get_llvm_tool_names('llvm-cov')
        for tool in tools:
            if shutil.which(tool):
                return tool
    return None

def compute_llvm_suffix(coredata: coredata.CoreData) -> T.Optional[str]:
    # Check to see if the user is trying to do coverage for either a C or C++ project
    compilers = coredata.compilers[MachineChoice.BUILD]
    cpp_compiler_is_clang = 'cpp' in compilers and compilers['cpp'].id == 'clang'
    c_compiler_is_clang = 'c' in compilers and compilers['c'].id == 'clang'
    # Extract first the C++ compiler if available. If it's a Clang of some kind, compute the suffix if possible
    if cpp_compiler_is_clang:
        suffix = compilers['cpp'].version.split('.')[0]
        return suffix

    # Then the C compiler, again checking if it's some kind of Clang and computing the suffix
    if c_compiler_is_clang:
        suffix = compilers['c'].version.split('.')[0]
        return suffix

    # Neither compiler is a Clang, or no compilers are for C or C++
    return None

def detect_lcov_genhtml(lcov_exe: str = 'lcov', genhtml_exe: str = 'genhtml') \
        -> T.Tuple[str, T.Optional[str], str]:
    lcov_exe, lcov_version = detect_lcov(lcov_exe)
    if shutil.which(genhtml_exe) is None:
        genhtml_exe = None

    return lcov_exe, lcov_version, genhtml_exe

def find_coverage_tools(coredata: coredata.CoreData) -> T.Tuple[T.Optional[str], T.Optional[str], T.Optional[str], T.Optional[str], T.Optional[str], T.Optional[str]]:
    gcovr_exe, gcovr_version = detect_gcovr()

    llvm_cov_exe = detect_llvm_cov(compute_llvm_suffix(coredata))
    # Some platforms may provide versioned clang but only non-versioned llvm utils
    if llvm_cov_exe is None:
        llvm_cov_exe = detect_llvm_cov('')

    lcov_exe, lcov_version, genhtml_exe = detect_lcov_genhtml()

    return gcovr_exe, gcovr_version, lcov_exe, lcov_version, genhtml_exe, llvm_cov_exe

def detect_ninja(version: str = '1.8.2', log: bool = False) -> T.Optional[T.List[str]]:
    r = detect_ninja_command_and_version(version, log)
    return r[0] if r else None

def detect_ninja_command_and_version(version: str = '1.8.2', log: bool = False) -> T.Optional[T.Tuple[T.List[str], str]]:
    env_ninja = os.environ.get('NINJA', None)
    for n in [env_ninja] if env_ninja else ['ninja', 'ninja-build', 'samu']:
        prog = ExternalProgram(n, silent=True)
        if not prog.found():
            continue
        try:
            p, found = Popen_safe(prog.command + ['--version'])[0:2]
        except (FileNotFoundError, PermissionError):
            # Doesn't exist in PATH or isn't executable
            continue
        found = found.strip()
        # Perhaps we should add a way for the caller to know the failure mode
        # (not found or too old)
        if p.returncode == 0 and mesonlib.version_compare(found, '>=' + version):
            if log:
                name = os.path.basename(n)
                if name.endswith('-' + found):
                    name = name[0:-1 - len(found)]
                if name == 'ninja-build':
                    name = 'ninja'
                if name == 'samu':
                    name = 'samurai'
                mlog.log('Found {}-{} at {}'.format(name, found,
                         ' '.join([quote_arg(x) for x in prog.command])))
            return (prog.command, found)
    return None

def get_llvm_tool_names(tool: str) -> T.List[str]:
    # Ordered list of possible suffixes of LLVM executables to try. Start with
    # base, then try newest back to oldest (3.5 is arbitrary), and finally the
    # devel version. Please note that the development snapshot in Debian does
    # not have a distinct name. Do not move it to the beginning of the list
    # unless it becomes a stable release.
    suffixes = [
        '', # base (no suffix)
        '-22.1', '22.1',
        '-22', '22',
        '-21.1', '21.1',
        '-21',  '21',
        '-20.1', '20.1',
        '-20',  '20',
        '-19.1', '19.1',
        '-19',  '19',
        '-18.1', '18.1',
        '-18',  '18',
        '-17',  '17',
        '-16',  '16',
        '-15',  '15',
        '-14',  '14',
        '-13',  '13',
        '-12',  '12',
        '-11',  '11',
        '-10',  '10',
        '-9',   '90',
        '-8',   '80',
        '-7',   '70',
        '-6.0', '60',
        '-5.0', '50',
        '-4.0', '40',
        '-3.9', '39',
        '-3.8', '38',
        '-3.7', '37',
        '-3.6', '36',
        '-3.5', '35',
        #'-20',    # Debian development snapshot
        '-devel', # FreeBSD development snapshot
    ]
    names: T.List[str] = []
    for suffix in suffixes:
        names.append(tool + suffix)
    return names

def detect_scanbuild() -> T.List[str]:
    """ Look for scan-build binary on build platform

    First, if a SCANBUILD env variable has been provided, give it precedence
    on all platforms.

    For most platforms, scan-build is found is the PATH contains a binary
    named "scan-build". However, some distribution's package manager (FreeBSD)
    don't. For those, loop through a list of candidates to see if one is
    available.

    Return: a single-element list of the found scan-build binary ready to be
        passed to Popen()
    """
    exelist: T.List[str] = []
    if 'SCANBUILD' in os.environ:
        exelist = split_args(os.environ['SCANBUILD'])

    else:
        tools = get_llvm_tool_names('scan-build')
        for tool in tools:
            which = shutil.which(tool)
            if which is not None:
                exelist = [which]
                break

    if exelist:
        tool = exelist[0]
        if os.path.isfile(tool) and os.access(tool, os.X_OK):
            return [tool]
    return []

def detect_clangformat() -> T.List[str]:
    """ Look for clang-format binary on build platform

    Do the same thing as detect_scanbuild to find clang-format except it
    currently does not check the environment variable.

    Return: a single-element list of the found clang-format binary ready to be
        passed to Popen()
    """
    tools = get_llvm_tool_names('clang-format')
    for tool in tools:
        path = shutil.which(tool)
        if path is not None:
            return [path]
    return []

def detect_clangtidy() -> T.List[str]:
    """ Look for clang-tidy binary on build platform

    Return: a single-element list of the found clang-tidy binary ready to be
        passed to Popen()
    """
    tools = get_llvm_tool_names('clang-tidy')
    for tool in tools:
        path = shutil.which(tool)
        if path is not None:
            return [path]
    return []

def detect_clangapply() -> T.List[str]:
    """ Look for clang-apply-replacements binary on build platform

    Return: a single-element list of the found clang-apply-replacements binary
        ready to be passed to Popen()
    """
    tools = get_llvm_tool_names('clang-apply-replacements')
    for tool in tools:
        path = shutil.which(tool)
        if path is not None:
            return [path]
    return []
