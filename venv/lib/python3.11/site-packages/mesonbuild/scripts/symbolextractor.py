# SPDX-License-Identifier: Apache-2.0
# Copyright 2013-2016 The Meson development team

# This script extracts the symbols of a given shared library
# into a file. If the symbols have not changed, the file is not
# touched. This information is used to skip link steps if the
# ABI has not changed.

# This file is basically a reimplementation of
# http://cgit.freedesktop.org/libreoffice/core/commit/?id=3213cd54b76bc80a6f0516aac75a48ff3b2ad67c
from __future__ import annotations

import typing as T
import os, sys
from .. import mesonlib
from .. import mlog
from ..mesonlib import Popen_safe
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--cross-host', default=None, dest='cross_host',
                    help='cross compilation host platform')
parser.add_argument('args', nargs='+')

TOOL_WARNING_FILE = None
RELINKING_WARNING = 'Relinking will always happen on source changes.'

def dummy_syms(outfilename: str) -> None:
    """Just touch it so relinking happens always."""
    with open(outfilename, 'w', encoding='utf-8'):
        pass

def write_if_changed(text: str, outfilename: str) -> None:
    try:
        with open(outfilename, encoding='utf-8') as f:
            oldtext = f.read()
        if text == oldtext:
            return
    except FileNotFoundError:
        pass
    with open(outfilename, 'w', encoding='utf-8') as f:
        f.write(text)

def print_tool_warning(tools: T.List[str], msg: str, stderr: T.Optional[str] = None) -> None:
    if os.path.exists(TOOL_WARNING_FILE):
        return
    m = f'{tools!r} {msg}. {RELINKING_WARNING}'
    if stderr:
        m += '\n' + stderr
    mlog.warning(m)
    # Write it out so we don't warn again
    with open(TOOL_WARNING_FILE, 'w', encoding='utf-8'):
        pass

def get_tool(name: str) -> T.List[str]:
    evar = name.upper()
    if evar in os.environ:
        import shlex
        return shlex.split(os.environ[evar])
    return [name]

def call_tool(name: str, args: T.List[str], **kwargs: T.Any) -> str:
    tool = get_tool(name)
    try:
        p, output, e = Popen_safe(tool + args, **kwargs)
    except FileNotFoundError:
        print_tool_warning(tool, 'not found')
        return None
    except PermissionError:
        print_tool_warning(tool, 'not usable')
        return None
    if p.returncode != 0:
        print_tool_warning(tool, 'does not work', e)
        return None
    return output

def call_tool_nowarn(tool: T.List[str], **kwargs: T.Any) -> T.Tuple[str, str]:
    try:
        p, output, e = Popen_safe(tool, **kwargs)
    except FileNotFoundError:
        return None, '{!r} not found\n'.format(tool[0])
    except PermissionError:
        return None, '{!r} not usable\n'.format(tool[0])
    if p.returncode != 0:
        return None, e
    return output, None

def gnu_syms(libfilename: str, outfilename: str) -> None:
    # Get the name of the library
    output = call_tool('readelf', ['-d', libfilename])
    if not output:
        dummy_syms(outfilename)
        return
    result = [x for x in output.split('\n') if 'SONAME' in x]
    assert len(result) <= 1
    # Get a list of all symbols exported
    output = call_tool('nm', ['--dynamic', '--extern-only', '--defined-only',
                              '--format=posix', libfilename])
    if not output:
        dummy_syms(outfilename)
        return
    for line in output.split('\n'):
        if not line:
            continue
        line_split = line.split()
        entry = line_split[0:2]
        # Store the size of symbols pointing to data objects so we relink
        # when those change, which is needed because of copy relocations
        # https://github.com/mesonbuild/meson/pull/7132#issuecomment-628353702
        if line_split[1].upper() in {'B', 'G', 'D'} and len(line_split) >= 4:
            entry += [line_split[3]]
        result += [' '.join(entry)]
    write_if_changed('\n'.join(result) + '\n', outfilename)

def solaris_syms(libfilename: str, outfilename: str) -> None:
    # gnu_syms() works with GNU nm & readelf, not Solaris nm & elfdump
    origpath = os.environ['PATH']
    try:
        os.environ['PATH'] = '/usr/gnu/bin:' + origpath
        gnu_syms(libfilename, outfilename)
    finally:
        os.environ['PATH'] = origpath

def osx_syms(libfilename: str, outfilename: str) -> None:
    # Get the name of the library
    output = call_tool('otool', ['-l', libfilename])
    if not output:
        dummy_syms(outfilename)
        return
    arr = output.split('\n')
    for (i, val) in enumerate(arr):
        if 'LC_ID_DYLIB' in val:
            match = i
            break
    result = [arr[match + 2], arr[match + 5]] # Libreoffice stores all 5 lines but the others seem irrelevant.
    # Get a list of all symbols exported.  `nm -g -U -P` is equivalent to, and more portable than,
    # `nm --extern-only --defined-only --format=posix`; cctools-port only understands the one-character form,
    # as does `nm` on very old macOS versions, (see meson#11131). `llvm-nm` understands both forms.
    output = call_tool('nm', ['-g', '-U', '-P', libfilename])
    if not output:
        dummy_syms(outfilename)
        return
    result += [' '.join(x.split()[0:2]) for x in output.split('\n')]
    write_if_changed('\n'.join(result) + '\n', outfilename)

def openbsd_syms(libfilename: str, outfilename: str) -> None:
    # Get the name of the library
    output = call_tool('readelf', ['-d', libfilename])
    if not output:
        dummy_syms(outfilename)
        return
    result = [x for x in output.split('\n') if 'SONAME' in x]
    assert len(result) <= 1
    # Get a list of all symbols exported
    output = call_tool('nm', ['-D', '-P', '-g', libfilename])
    if not output:
        dummy_syms(outfilename)
        return
    # U = undefined (cope with the lack of --defined-only option)
    result += [' '.join(x.split()[0:2]) for x in output.split('\n') if x and not x.endswith('U ')]
    write_if_changed('\n'.join(result) + '\n', outfilename)

def freebsd_syms(libfilename: str, outfilename: str) -> None:
    # Get the name of the library
    output = call_tool('readelf', ['-d', libfilename])
    if not output:
        dummy_syms(outfilename)
        return
    result = [x for x in output.split('\n') if 'SONAME' in x]
    assert len(result) <= 1
    # Get a list of all symbols exported
    output = call_tool('nm', ['--dynamic', '--extern-only', '--defined-only',
                              '--format=posix', libfilename])
    if not output:
        dummy_syms(outfilename)
        return

    result += [' '.join(x.split()[0:2]) for x in output.split('\n')]
    write_if_changed('\n'.join(result) + '\n', outfilename)

def cygwin_syms(impfilename: str, outfilename: str) -> None:
    # Get the name of the library
    output = call_tool('dlltool', ['-I', impfilename])
    if not output:
        dummy_syms(outfilename)
        return
    result = [output]
    # Get the list of all symbols exported
    output = call_tool('nm', ['--extern-only', '--defined-only',
                              '--format=posix', impfilename])
    if not output:
        dummy_syms(outfilename)
        return
    for line in output.split('\n'):
        if ' T ' not in line:
            continue
        result.append(line.split(maxsplit=1)[0])
    write_if_changed('\n'.join(result) + '\n', outfilename)

def _get_implib_dllname(impfilename: str) -> T.Tuple[T.List[str], str]:
    all_stderr = ''
    # First try lib.exe, which is provided by MSVC. Then llvm-lib.exe, by LLVM
    # for clang-cl.
    #
    # We cannot call get_tool on `lib` because it will look at the `LIB` env
    # var which is the list of library paths MSVC will search for import
    # libraries while linking.
    for lib in (['lib'], get_tool('llvm-lib')):
        output, e = call_tool_nowarn(lib + ['-list', impfilename])
        if output:
            # The output is a list of DLLs that each symbol exported by the import
            # library is available in. We only build import libraries that point to
            # a single DLL, so we can pick any of these. Pick the last one for
            # simplicity. Also skip the last line, which is empty.
            return output.split('\n')[-2:-1], None
        all_stderr += e
    # Next, try dlltool.exe which is provided by MinGW
    output, e = call_tool_nowarn(get_tool('dlltool') + ['-I', impfilename])
    if output:
        return [output], None
    all_stderr += e
    return ([], all_stderr)

def _get_implib_exports(impfilename: str) -> T.Tuple[T.List[str], str]:
    all_stderr = ''
    # Force dumpbin.exe to use en-US so we can parse its output
    env = os.environ.copy()
    env['VSLANG'] = '1033'
    output, e = call_tool_nowarn(get_tool('dumpbin') + ['-exports', impfilename], env=env)
    if output:
        lines = output.split('\n')
        start = lines.index('File Type: LIBRARY')
        end = lines.index('  Summary')
        return lines[start:end], None
    all_stderr += e
    # Next, try llvm-nm.exe provided by LLVM, then nm.exe provided by MinGW
    for nm in ('llvm-nm', 'nm'):
        output, e = call_tool_nowarn(get_tool(nm) + ['--extern-only', '--defined-only',
                                                     '--format=posix', impfilename])
        if output:
            result = []
            for line in output.split('\n'):
                if ' T ' not in line or line.startswith('.text'):
                    continue
                result.append(line.split(maxsplit=1)[0])
            return result, None
        all_stderr += e
    return ([], all_stderr)

def windows_syms(impfilename: str, outfilename: str) -> None:
    # Get the name of the library
    result, e = _get_implib_dllname(impfilename)
    if not result:
        print_tool_warning(['lib', 'llvm-lib', 'dlltool'], 'do not work or were not found', e)
        dummy_syms(outfilename)
        return
    # Get a list of all symbols exported
    symbols, e = _get_implib_exports(impfilename)
    if not symbols:
        print_tool_warning(['dumpbin', 'llvm-nm', 'nm'], 'do not work or were not found', e)
        dummy_syms(outfilename)
        return
    result += symbols
    write_if_changed('\n'.join(result) + '\n', outfilename)

def gen_symbols(libfilename: str, impfilename: str, outfilename: str, cross_host: str) -> None:
    if cross_host is not None:
        # In case of cross builds just always relink. In theory we could
        # determine the correct toolset, but we would need to use the correct
        # `nm`, `readelf`, etc, from the cross info which requires refactoring.
        dummy_syms(outfilename)
    elif mesonlib.is_linux() or mesonlib.is_hurd():
        gnu_syms(libfilename, outfilename)
    elif mesonlib.is_osx():
        osx_syms(libfilename, outfilename)
    elif mesonlib.is_openbsd():
        openbsd_syms(libfilename, outfilename)
    elif mesonlib.is_freebsd():
        freebsd_syms(libfilename, outfilename)
    elif mesonlib.is_netbsd():
        freebsd_syms(libfilename, outfilename)
    elif mesonlib.is_windows():
        if os.path.isfile(impfilename):
            windows_syms(impfilename, outfilename)
        else:
            # No import library. Not sure how the DLL is being used, so just
            # rebuild everything that links to it every time.
            dummy_syms(outfilename)
    elif mesonlib.is_cygwin():
        if os.path.isfile(impfilename):
            cygwin_syms(impfilename, outfilename)
        else:
            # No import library. Not sure how the DLL is being used, so just
            # rebuild everything that links to it every time.
            dummy_syms(outfilename)
    elif mesonlib.is_sunos():
        solaris_syms(libfilename, outfilename)
    else:
        if not os.path.exists(TOOL_WARNING_FILE):
            mlog.warning('Symbol extracting has not been implemented for this '
                         'platform. ' + RELINKING_WARNING)
            # Write it out so we don't warn again
            with open(TOOL_WARNING_FILE, 'w', encoding='utf-8'):
                pass
        dummy_syms(outfilename)

def run(args: T.List[str]) -> int:
    global TOOL_WARNING_FILE  # pylint: disable=global-statement
    options = parser.parse_args(args)
    if len(options.args) != 4:
        print('symbolextractor.py <shared library file> <import library> <output file>')
        sys.exit(1)
    privdir = os.path.join(options.args[0], 'meson-private')
    TOOL_WARNING_FILE = os.path.join(privdir, 'symbolextractor_tool_warning_printed')
    libfile = options.args[1]
    impfile = options.args[2] # Only used on Windows
    outfile = options.args[3]
    gen_symbols(libfile, impfile, outfile, options.cross_host)
    return 0

if __name__ == '__main__':
    sys.exit(run(sys.argv[1:]))
