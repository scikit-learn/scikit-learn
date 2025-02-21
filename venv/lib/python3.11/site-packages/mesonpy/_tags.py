# SPDX-FileCopyrightText: 2022 The meson-python developers
#
# SPDX-License-Identifier: MIT

from __future__ import annotations

import os
import platform
import struct
import sys
import sysconfig
import typing


if typing.TYPE_CHECKING:  # pragma: no cover
    from typing import Optional, Union


# https://peps.python.org/pep-0425/#python-tag
INTERPRETERS = {
    'python': 'py',
    'cpython': 'cp',
    'pypy': 'pp',
    'ironpython': 'ip',
    'jython': 'jy',
}


_32_BIT_INTERPRETER = struct.calcsize('P') == 4


def get_interpreter_tag() -> str:
    name = sys.implementation.name
    name = INTERPRETERS.get(name, name)
    version = sys.version_info
    return f'{name}{version[0]}{version[1]}'


def _get_config_var(name: str, default: Union[str, int, None] = None) -> Union[str, int, None]:
    value: Union[str, int, None] = sysconfig.get_config_var(name)
    if value is None:
        return default
    return value


def _get_cpython_abi() -> str:
    version = sys.version_info
    debug = pymalloc = ''
    if _get_config_var('Py_DEBUG', hasattr(sys, 'gettotalrefcount')):
        debug = 'd'
    if version < (3, 8) and _get_config_var('WITH_PYMALLOC', True):
        pymalloc = 'm'
    return f'cp{version[0]}{version[1]}{debug}{pymalloc}'


def get_abi_tag() -> str:
    # The best solution to obtain the Python ABI is to parse the
    # $SOABI or $EXT_SUFFIX sysconfig variables as defined in PEP-314.

    # PyPy reports a $SOABI that does not agree with $EXT_SUFFIX.
    # Using $EXT_SUFFIX will not break when PyPy will fix this.
    # See https://foss.heptapod.net/pypy/pypy/-/issues/3816 and
    # https://github.com/pypa/packaging/pull/607.
    try:
        empty, abi, ext = str(sysconfig.get_config_var('EXT_SUFFIX')).split('.')
    except ValueError as exc:
        # CPython <= 3.8.7 on Windows does not implement PEP3149 and
        # uses '.pyd' as $EXT_SUFFIX, which does not allow to extract
        # the interpreter ABI.  Check that the fallback is not hit for
        # any other Python implementation.
        if sys.implementation.name != 'cpython':
            raise NotImplementedError from exc
        return _get_cpython_abi()

    # The packaging module initially based his understanding of the
    # $SOABI variable on the inconsistent value reported by PyPy, and
    # did not strip architecture information from it.  Therefore the
    # ABI tag for later Python implementations (all the ones not
    # explicitly handled below) contains architecture information too.
    # Unfortunately, fixing this now would break compatibility.

    if abi.startswith('cpython'):
        abi = 'cp' + abi.split('-')[1]
    elif abi.startswith('cp'):
        abi = abi.split('-')[0]
    elif abi.startswith('pypy'):
        abi = '_'.join(abi.split('-')[:2])
    elif abi.startswith('graalpy'):
        abi = '_'.join(abi.split('-')[:3])

    return abi.replace('.', '_').replace('-', '_')


def _get_macosx_platform_tag() -> str:
    ver, _, arch = platform.mac_ver()

    # Override the architecture with the one provided in the
    # _PYTHON_HOST_PLATFORM environment variable.  This environment
    # variable affects the sysconfig.get_platform() return value and
    # is used to cross-compile python extensions on macOS for a
    # different architecture.  We base the platform tag computation on
    # platform.mac_ver() but respect the content of the environment
    # variable.
    try:
        arch = os.environ.get('_PYTHON_HOST_PLATFORM', '').split('-')[2]
    except IndexError:
        pass

    # Override the macOS version if one is provided via the
    # MACOSX_DEPLOYMENT_TARGET environment variable.
    try:
        version = tuple(map(int, os.environ.get('MACOSX_DEPLOYMENT_TARGET', '').split('.')))[:2]
    except ValueError:
        version = tuple(map(int, ver.split('.')))[:2]

    # Python built with older macOS SDK on macOS 11, reports an
    # unexising macOS 10.16 version instead of the real version.
    #
    # The packaging module introduced a workaround
    # https://github.com/pypa/packaging/commit/67c4a2820c549070bbfc4bfbf5e2a250075048da
    #
    # This results in packaging versions up to 21.3 generating
    # platform tags like "macosx_10_16_x86_64" and later versions
    # generating "macosx_11_0_x86_64".  Using the latter would be more
    # correct but prevents the resulting wheel from being installed on
    # systems using packaging 21.3 or earlier (pip 22.3 or earlier).
    #
    # Fortunately packaging versions carrying the workaround still
    # accepts "macosx_10_16_x86_64" as a compatible platform tag.  We
    # can therefore ignore the issue and generate the slightly
    # incorrect tag.

    # The minimum macOS ABI version on arm64 is 11.0.  The macOS SDK
    # on arm64 silently bumps any compatibility version specified via
    # the MACOSX_DEPLOYMENT_TARGET environment variable to 11.0.
    # Despite the platform ABI tag being intended to be a minimum
    # compatibility version, pip refuses to install wheels with a
    # platform tag specifying an ABI version lower than 11.0.  Use
    # 11.0 as minimum ABI version on arm64.
    if arch == 'arm64' and version < (11, 0):
        version = (11, 0)

    major, minor = version

    if major >= 11:
        # For macOS reelases up to 10.15, the major version number is
        # actually part of the OS name and the minor version is the
        # actual OS release.  Starting with macOS 11, the major
        # version number is the OS release and the minor version is
        # the patch level.  Reset the patch level to zero.
        minor = 0

    if _32_BIT_INTERPRETER:
        # 32-bit Python running on a 64-bit kernel.
        if arch == 'ppc64':
            arch = 'ppc'
        if arch == 'x86_64':
            arch = 'i386'

    return f'macosx_{major}_{minor}_{arch}'


def get_platform_tag() -> str:
    platform = sysconfig.get_platform()
    if platform.startswith('macosx'):
        return _get_macosx_platform_tag()
    if _32_BIT_INTERPRETER:
        # 32-bit Python running on a 64-bit kernel.
        if platform == 'linux-x86_64':
            return 'linux_i686'
        if platform == 'linux-aarch64':
            return 'linux_armv7l'
    return platform.replace('-', '_').replace('.', '_').lower()


class Tag:
    def __init__(self, interpreter: Optional[str] = None, abi: Optional[str] = None, platform: Optional[str] = None):
        self.interpreter = interpreter or get_interpreter_tag()
        self.abi = abi or get_abi_tag()
        self.platform = platform or get_platform_tag()

    def __str__(self) -> str:
        return f'{self.interpreter}-{self.abi}-{self.platform}'
