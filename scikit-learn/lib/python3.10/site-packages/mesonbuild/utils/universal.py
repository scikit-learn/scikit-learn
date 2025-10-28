# SPDX-License-Identifier: Apache-2.0
# Copyright 2012-2020 The Meson development team


"""A library of random helper functionality."""

from __future__ import annotations
from pathlib import Path
import argparse
import ast
import enum
import sys
import stat
import time
import abc
import multiprocessing
import platform, subprocess, operator, os, shlex, shutil, re
import collections
from functools import lru_cache, wraps
from itertools import tee
from tempfile import TemporaryDirectory, NamedTemporaryFile
import typing as T
import textwrap
import pickle
import errno
import json
import dataclasses

from mesonbuild import mlog
from .core import MesonException, HoldableObject

if T.TYPE_CHECKING:
    from typing_extensions import Literal, Protocol

    from .._typing import ImmutableListProtocol
    from ..build import ConfigurationData
    from ..coredata import StrOrBytesPath
    from ..environment import Environment
    from ..compilers.compilers import Compiler
    from ..interpreterbase.baseobjects import SubProject
    from .. import programs

    class _EnvPickleLoadable(Protocol):

        environment: Environment

    class _VerPickleLoadable(Protocol):

        version: str

    # A generic type for pickle_load. This allows any type that has either a
    # .version or a .environment to be passed.
    _PL = T.TypeVar('_PL', bound=T.Union[_EnvPickleLoadable, _VerPickleLoadable])

FileOrString = T.Union['File', str]

_T = T.TypeVar('_T')
_U = T.TypeVar('_U')

__all__ = [
    'GIT',
    'python_command',
    'NoProjectVersion',
    'project_meson_versions',
    'SecondLevelHolder',
    'File',
    'FileMode',
    'GitException',
    'LibType',
    'MachineChoice',
    'EnvironmentException',
    'FileOrString',
    'GitException',
    'dump_conf_header',
    'OrderedSet',
    'PerMachine',
    'PerMachineDefaultable',
    'PerThreeMachine',
    'PerThreeMachineDefaultable',
    'ProgressBar',
    'RealPathAction',
    'TemporaryDirectoryWinProof',
    'Version',
    'check_direntry_issues',
    'classify_unity_sources',
    'current_vs_supports_modules',
    'darwin_get_object_archs',
    'default_libdir',
    'default_libexecdir',
    'default_prefix',
    'default_datadir',
    'default_includedir',
    'default_infodir',
    'default_localedir',
    'default_mandir',
    'default_sbindir',
    'default_sysconfdir',
    'detect_subprojects',
    'detect_vcs',
    'determine_worker_count',
    'do_conf_file',
    'do_conf_str',
    'do_replacement',
    'expand_arguments',
    'extract_as_list',
    'first',
    'generate_list',
    'get_compiler_for_source',
    'get_filenames_templates_dict',
    'get_rsp_threshold',
    'get_variable_regex',
    'get_wine_shortpath',
    'git',
    'has_path_sep',
    'is_aix',
    'is_android',
    'is_ascii_string',
    'is_cygwin',
    'is_debianlike',
    'is_dragonflybsd',
    'is_freebsd',
    'is_haiku',
    'is_hurd',
    'is_irix',
    'is_linux',
    'is_netbsd',
    'is_openbsd',
    'is_osx',
    'is_parent_path',
    'is_qnx',
    'is_sunos',
    'is_windows',
    'is_wsl',
    'iter_regexin_iter',
    'join_args',
    'lazy_property',
    'listify',
    'listify_array_value',
    'partition',
    'path_is_in_root',
    'pickle_load',
    'Popen_safe',
    'Popen_safe_logged',
    'quiet_git',
    'quote_arg',
    'relative_to_if_possible',
    'relpath',
    'replace_if_different',
    'run_once',
    'get_meson_command',
    'set_meson_command',
    'split_args',
    'stringlistify',
    'substitute_values',
    'substring_is_in_list',
    'typeslistify',
    'verbose_git',
    'version_compare',
    'version_compare_condition_with_min',
    'version_compare_many',
    'search_version',
    'windows_detect_native_arch',
    'windows_proof_rm',
    'windows_proof_rmtree',
]


class NoProjectVersion:
    pass

# TODO: this is such a hack, this really should be either in coredata or in the
# interpreter
# {subproject: project_meson_version}
project_meson_versions: T.Dict[str, T.Union[str, NoProjectVersion]] = {}


from glob import glob

if getattr(sys, 'frozen', False):
    # Using e.g. a PyInstaller bundle, such as the MSI installed executable.
    # It is conventional for freeze programs to set this attribute to indicate
    # that the program is self hosted, and for example there is no associated
    # "python" executable.
    python_command = [sys.executable, 'runpython']
else:
    python_command = [sys.executable]
_meson_command: T.Optional['ImmutableListProtocol[str]'] = None


class EnvironmentException(MesonException):
    '''Exceptions thrown while processing and creating the build environment'''

class GitException(MesonException):
    def __init__(self, msg: str, output: T.Optional[str] = None):
        super().__init__(msg)
        self.output = output.strip() if output else ''

GIT = shutil.which('git')
def git(cmd: T.List[str], workingdir: StrOrBytesPath, check: bool = False, **kwargs: T.Any) -> T.Tuple[subprocess.Popen[str], str, str]:
    assert GIT is not None, 'Callers should make sure it exists'
    cmd = [GIT, *cmd]
    p, o, e = Popen_safe(cmd, cwd=workingdir, **kwargs)
    if check and p.returncode != 0:
        raise GitException('Git command failed: ' + str(cmd), e)
    return p, o, e

def quiet_git(cmd: T.List[str], workingdir: StrOrBytesPath, check: bool = False) -> T.Tuple[bool, str]:
    if not GIT:
        m = 'Git program not found.'
        if check:
            raise GitException(m)
        return False, m
    p, o, e = git(cmd, workingdir, check)
    if p.returncode != 0:
        return False, e
    return True, o

def verbose_git(cmd: T.List[str], workingdir: StrOrBytesPath, check: bool = False) -> bool:
    if not GIT:
        m = 'Git program not found.'
        if check:
            raise GitException(m)
        return False
    p, _, _ = git(cmd, workingdir, check, stdout=None, stderr=None)
    return p.returncode == 0

def set_meson_command(mainfile: str) -> None:
    global _meson_command  # pylint: disable=global-statement
    # On UNIX-like systems `meson` is a Python script
    # On Windows `meson` and `meson.exe` are wrapper exes
    if not mainfile.endswith('.py'):
        _meson_command = [mainfile]
    elif os.path.isabs(mainfile) and mainfile.endswith('mesonmain.py'):
        # Can't actually run meson with an absolute path to mesonmain.py, it must be run as -m mesonbuild.mesonmain
        _meson_command = python_command + ['-m', 'mesonbuild.mesonmain']
    else:
        # Either run uninstalled, or full path to meson-script.py
        _meson_command = python_command + [mainfile]
    # We print this value for unit tests.
    if 'MESON_COMMAND_TESTS' in os.environ:
        mlog.log(f'meson_command is {_meson_command!r}')


def get_meson_command() -> T.Optional['ImmutableListProtocol[str]']:
    return _meson_command


def is_ascii_string(astring: T.Union[str, bytes]) -> bool:
    try:
        if isinstance(astring, str):
            astring.encode('ascii')
        elif isinstance(astring, bytes):
            astring.decode('ascii')
    except UnicodeDecodeError:
        return False
    return True


def check_direntry_issues(direntry_array: T.Union[T.Iterable[T.Union[str, bytes]], str, bytes]) -> None:
    import locale
    # Warn if the locale is not UTF-8. This can cause various unfixable issues
    # such as os.stat not being able to decode filenames with unicode in them.
    # There is no way to reset both the preferred encoding and the filesystem
    # encoding, so we can just warn about it.
    e = locale.getpreferredencoding()
    if e.upper() != 'UTF-8' and not is_windows():
        if isinstance(direntry_array, (str, bytes)):
            direntry_array = [direntry_array]
        for de in direntry_array:
            if is_ascii_string(de):
                continue
            mlog.warning(textwrap.dedent(f'''
                You are using {e!r} which is not a Unicode-compatible
                locale but you are trying to access a file system entry called {de!r} which is
                not pure ASCII. This may cause problems.
                '''))

class SecondLevelHolder(HoldableObject, metaclass=abc.ABCMeta):
    ''' A second level object holder. The primary purpose
        of such objects is to hold multiple objects with one
        default option. '''

    @abc.abstractmethod
    def get_default_object(self) -> HoldableObject: ...

class FileMode:
    # The first triad is for owner permissions, the second for group permissions,
    # and the third for others (everyone else).
    # For the 1st character:
    #  'r' means can read
    #  '-' means not allowed
    # For the 2nd character:
    #  'w' means can write
    #  '-' means not allowed
    # For the 3rd character:
    #  'x' means can execute
    #  's' means can execute and setuid/setgid is set (owner/group triads only)
    #  'S' means cannot execute and setuid/setgid is set (owner/group triads only)
    #  't' means can execute and sticky bit is set ("others" triads only)
    #  'T' means cannot execute and sticky bit is set ("others" triads only)
    #  '-' means none of these are allowed
    #
    # The meanings of 'rwx' perms is not obvious for directories; see:
    # https://www.hackinglinuxexposed.com/articles/20030424.html
    #
    # For information on this notation such as setuid/setgid/sticky bits, see:
    # https://en.wikipedia.org/wiki/File_system_permissions#Symbolic_notation
    symbolic_perms_regex = re.compile('[r-][w-][xsS-]' # Owner perms
                                      '[r-][w-][xsS-]' # Group perms
                                      '[r-][w-][xtT-]') # Others perms

    def __init__(self, perms: T.Optional[str] = None, owner: T.Union[str, int, None] = None,
                 group: T.Union[str, int, None] = None):
        self.perms_s = perms
        self.perms = self.perms_s_to_bits(perms)
        self.owner = owner
        self.group = group

    def __repr__(self) -> str:
        ret = '<FileMode: {!r} owner={} group={}'
        return ret.format(self.perms_s, self.owner, self.group)

    @classmethod
    def perms_s_to_bits(cls, perms_s: T.Optional[str]) -> int:
        '''
        Does the opposite of stat.filemode(), converts strings of the form
        'rwxr-xr-x' to st_mode enums which can be passed to os.chmod()
        '''
        if perms_s is None:
            # No perms specified, we will not touch the permissions
            return -1
        eg = 'rwxr-xr-x'
        if not isinstance(perms_s, str):
            raise MesonException(f'Install perms must be a string. For example, {eg!r}')
        if len(perms_s) != 9 or not cls.symbolic_perms_regex.match(perms_s):
            raise MesonException(f'File perms {perms_s!r} must be exactly 9 chars. For example, {eg!r}')
        perms = 0
        # Owner perms
        if perms_s[0] == 'r':
            perms |= stat.S_IRUSR
        if perms_s[1] == 'w':
            perms |= stat.S_IWUSR
        if perms_s[2] == 'x':
            perms |= stat.S_IXUSR
        elif perms_s[2] == 'S':
            perms |= stat.S_ISUID
        elif perms_s[2] == 's':
            perms |= stat.S_IXUSR
            perms |= stat.S_ISUID
        # Group perms
        if perms_s[3] == 'r':
            perms |= stat.S_IRGRP
        if perms_s[4] == 'w':
            perms |= stat.S_IWGRP
        if perms_s[5] == 'x':
            perms |= stat.S_IXGRP
        elif perms_s[5] == 'S':
            perms |= stat.S_ISGID
        elif perms_s[5] == 's':
            perms |= stat.S_IXGRP
            perms |= stat.S_ISGID
        # Others perms
        if perms_s[6] == 'r':
            perms |= stat.S_IROTH
        if perms_s[7] == 'w':
            perms |= stat.S_IWOTH
        if perms_s[8] == 'x':
            perms |= stat.S_IXOTH
        elif perms_s[8] == 'T':
            perms |= stat.S_ISVTX
        elif perms_s[8] == 't':
            perms |= stat.S_IXOTH
            perms |= stat.S_ISVTX
        return perms

dot_C_dot_H_warning = """You are using .C or .H files in your project. This is deprecated.
         Currently, Meson treats this as C++ code, but they
            used to be treated as C code.
         Note that the situation is a bit more complex if you are using the
         Visual Studio compiler, as it treats .C files as C code, unless you add
         the /TP compiler flag, but this is unreliable.
         See https://github.com/mesonbuild/meson/pull/8747 for the discussions."""
class File(HoldableObject):
    def __init__(self, is_built: bool, subdir: str, fname: str):
        if fname.endswith(".C") or fname.endswith(".H"):
            mlog.warning(dot_C_dot_H_warning, once=True)
        self.is_built = is_built
        self.subdir = subdir
        self.fname = fname
        self.hash = hash((is_built, subdir, fname))

    def __str__(self) -> str:
        return self.relative_name()

    def __repr__(self) -> str:
        ret = '<File: {0}'
        if not self.is_built:
            ret += ' (not built)'
        ret += '>'
        return ret.format(self.relative_name())

    @staticmethod
    @lru_cache(maxsize=None)
    def from_source_file(source_root: str, subdir: str, fname: str) -> File:
        if not os.path.isfile(os.path.join(source_root, subdir, fname)):
            raise MesonException(f'File {fname} does not exist.')
        return File(False, subdir, fname)

    @staticmethod
    @lru_cache(maxsize=None)
    def from_built_file(subdir: str, fname: str) -> 'File':
        return File(True, subdir, fname)

    @staticmethod
    @lru_cache(maxsize=None)
    def from_built_relative(relative: str) -> 'File':
        dirpart, fnamepart = os.path.split(relative)
        return File.from_built_file(dirpart, fnamepart)

    @staticmethod
    def from_absolute_file(fname: str) -> 'File':
        return File(False, '', fname)

    @lru_cache(maxsize=None)
    def rel_to_builddir(self, build_to_src: str) -> str:
        if self.is_built:
            return self.relative_name()
        else:
            return os.path.join(build_to_src, self.subdir, self.fname)

    @lru_cache(maxsize=None)
    def absolute_path(self, srcdir: str, builddir: str) -> str:
        absdir = srcdir
        if self.is_built:
            absdir = builddir
        return os.path.normpath(os.path.join(absdir, self.relative_name()))

    @property
    def suffix(self) -> str:
        return os.path.splitext(self.fname)[1][1:].lower()

    def endswith(self, ending: T.Union[str, T.Tuple[str, ...]]) -> bool:
        return self.fname.endswith(ending)

    def split(self, s: str, maxsplit: int = -1) -> T.List[str]:
        return self.fname.split(s, maxsplit=maxsplit)

    def rsplit(self, s: str, maxsplit: int = -1) -> T.List[str]:
        return self.fname.rsplit(s, maxsplit=maxsplit)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, File):
            return NotImplemented
        if self.hash != other.hash:
            return False
        return (self.fname, self.subdir, self.is_built) == (other.fname, other.subdir, other.is_built)

    def __hash__(self) -> int:
        return self.hash

    @lru_cache(maxsize=None)
    def relative_name(self) -> str:
        return os.path.join(self.subdir, self.fname)


def get_compiler_for_source(compilers: T.Iterable['Compiler'], src: 'FileOrString') -> 'Compiler':
    """Given a set of compilers and a source, find the compiler for that source type."""
    for comp in compilers:
        if comp.can_compile(src):
            return comp
    raise MesonException(f'No specified compiler can handle file {src!s}')


def classify_unity_sources(compilers: T.Iterable['Compiler'], sources: T.Sequence['FileOrString']) -> T.Dict['Compiler', T.List['FileOrString']]:
    compsrclist: T.Dict['Compiler', T.List['FileOrString']] = {}
    for src in sources:
        comp = get_compiler_for_source(compilers, src)
        if comp not in compsrclist:
            compsrclist[comp] = [src]
        else:
            compsrclist[comp].append(src)
    return compsrclist


MACHINE_NAMES = ['build', 'host']
MACHINE_PREFIXES = ['build.', '']


class MachineChoice(enum.IntEnum):

    """Enum class representing one of the two abstract machine names used in
    most places: the build, and host, machines.
    """

    BUILD = 0
    HOST = 1

    def __str__(self) -> str:
        return f'{self.get_lower_case_name()} machine'

    def get_lower_case_name(self) -> str:
        return MACHINE_NAMES[self.value]

    def get_prefix(self) -> str:
        return MACHINE_PREFIXES[self.value]


@dataclasses.dataclass(eq=False, order=False)
class PerMachine(T.Generic[_T]):
    build: _T
    host: _T

    def __getitem__(self, machine: MachineChoice) -> _T:
        return [self.build, self.host][machine.value]

    def __setitem__(self, machine: MachineChoice, val: _T) -> None:
        setattr(self, machine.get_lower_case_name(), val)

    def miss_defaulting(self) -> PerMachineDefaultable[T.Optional[_T]]:
        """Unset definition duplicated from their previous to None

        This is the inverse of ''default_missing''. By removing defaulted
        machines, we can elaborate the original and then redefault them and thus
        avoid repeating the elaboration explicitly.
        """
        unfreeze: PerMachineDefaultable[T.Optional[_T]] = PerMachineDefaultable()
        unfreeze.build = self.build
        unfreeze.host = self.host
        if unfreeze.host == unfreeze.build:
            unfreeze.host = None
        return unfreeze

    def assign(self, build: _T, host: _T) -> None:
        self.build = build
        self.host = host


@dataclasses.dataclass(eq=False, order=False)
class PerThreeMachine(PerMachine[_T]):
    """Like `PerMachine` but includes `target` too.

    It turns out just one thing do we need track the target machine. There's no
    need to computer the `target` field so we don't bother overriding the
    `__getitem__`/`__setitem__` methods.
    """

    target: _T

    def miss_defaulting(self) -> "PerThreeMachineDefaultable[T.Optional[_T]]":
        """Unset definition duplicated from their previous to None

        This is the inverse of ''default_missing''. By removing defaulted
        machines, we can elaborate the original and then redefault them and thus
        avoid repeating the elaboration explicitly.
        """
        unfreeze: PerThreeMachineDefaultable[T.Optional[_T]] = PerThreeMachineDefaultable()
        unfreeze.build = self.build
        unfreeze.host = self.host
        unfreeze.target = self.target
        if unfreeze.target == unfreeze.host:
            unfreeze.target = None
        if unfreeze.host == unfreeze.build:
            unfreeze.host = None
        return unfreeze

    def matches_build_machine(self, machine: MachineChoice) -> bool:
        return self.build == self[machine]


@dataclasses.dataclass(eq=False, order=False)
class PerMachineDefaultable(PerMachine[T.Optional[_T]]):
    """Extends `PerMachine` with the ability to default from `None`s.
    """

    build: T.Optional[_T] = None
    host: T.Optional[_T] = None

    def default_missing(self) -> PerMachine[_T]:
        """Default host to build

        This allows just specifying nothing in the native case, and just host in the
        cross non-compiler case.
        """
        assert self.build is not None, 'Cannot fill in missing when all fields are empty'
        return PerMachine(self.build, self.host if self.host is not None else self.build)

    @classmethod
    def default(cls, is_cross: bool, build: _T, host: _T) -> PerMachine[_T]:
        """Easy way to get a defaulted value

        This allows simplifying the case where you can control whether host and
        build are separate or not with a boolean. If the is_cross value is set
        to true then the optional host value will be used, otherwise the host
        will be set to the build value.
        """
        m = cls(build)
        if is_cross:
            m.host = host
        return m.default_missing()


@dataclasses.dataclass(eq=False, order=False)
class PerThreeMachineDefaultable(PerMachineDefaultable[T.Optional[_T]], PerThreeMachine[T.Optional[_T]]):
    """Extends `PerThreeMachine` with the ability to default from `None`s.
    """

    target: T.Optional[_T] = None

    def default_missing(self) -> PerThreeMachine[_T]:
        """Default host to build and target to host.

        This allows just specifying nothing in the native case, just host in the
        cross non-compiler case, and just target in the native-built
        cross-compiler case.
        """
        assert self.build is not None, 'Cannot default a PerMachine when all values are None'
        host = self.host if self.host is not None else self.build
        target = self.target if self.target is not None else host
        return PerThreeMachine(self.build, host, target)


def is_sunos() -> bool:
    return platform.system().lower() == 'sunos'


def is_osx() -> bool:
    return platform.system().lower() == 'darwin'


def is_linux() -> bool:
    return platform.system().lower() == 'linux'


def is_android() -> bool:
    return platform.system().lower() == 'android'


def is_haiku() -> bool:
    return platform.system().lower() == 'haiku'


def is_openbsd() -> bool:
    return platform.system().lower() == 'openbsd'


def is_windows() -> bool:
    platname = platform.system().lower()
    return platname == 'windows'

def is_wsl() -> bool:
    return is_linux() and 'microsoft' in platform.release().lower()

def is_cygwin() -> bool:
    return sys.platform == 'cygwin'


def is_debianlike() -> bool:
    return os.path.isfile('/etc/debian_version')


def is_dragonflybsd() -> bool:
    return platform.system().lower() == 'dragonfly'


def is_netbsd() -> bool:
    return platform.system().lower() == 'netbsd'


def is_freebsd() -> bool:
    return platform.system().lower() == 'freebsd'

def is_irix() -> bool:
    return platform.system().startswith('irix')

def is_hurd() -> bool:
    return platform.system().lower() == 'gnu'

def is_qnx() -> bool:
    return platform.system().lower() == 'qnx'

def is_aix() -> bool:
    return platform.system().lower() == 'aix'

@lru_cache(maxsize=None)
def darwin_get_object_archs(objpath: str) -> 'ImmutableListProtocol[str]':
    '''
    For a specific object (executable, static library, dylib, etc), run `lipo`
    to fetch the list of archs supported by it. Supports both thin objects and
    'fat' objects.
    '''
    _, stdo, stderr = Popen_safe(['lipo', '-info', objpath])
    if not stdo:
        mlog.debug(f'lipo {objpath}: {stderr}')
        return None
    stdo = stdo.rsplit(': ', 1)[1]

    # Convert from lipo-style archs to meson-style CPUs
    map_arch = {
        'i386': 'x86',
        'arm64': 'aarch64',
        'arm64e': 'aarch64',
        'ppc7400': 'ppc',
        'ppc970': 'ppc',
    }
    lipo_archs = stdo.split()
    meson_archs = [map_arch.get(lipo_arch, lipo_arch) for lipo_arch in lipo_archs]

    # Add generic name for armv7 and armv7s
    if 'armv7' in stdo:
        meson_archs.append('arm')

    return meson_archs

def windows_detect_native_arch() -> str:
    """
    The architecture of Windows itself: x86, amd64 or arm64
    """
    if sys.platform != 'win32':
        return ''
    try:
        import ctypes
        process_arch = ctypes.c_ushort()
        native_arch = ctypes.c_ushort()
        kernel32 = ctypes.windll.kernel32
        process = ctypes.c_void_p(kernel32.GetCurrentProcess())
        # This is the only reliable way to detect an arm system if we are an x86/x64 process being emulated
        if kernel32.IsWow64Process2(process, ctypes.byref(process_arch), ctypes.byref(native_arch)):
            # https://docs.microsoft.com/en-us/windows/win32/sysinfo/image-file-machine-constants
            if native_arch.value == 0x8664:
                return 'amd64'
            elif native_arch.value == 0x014C:
                return 'x86'
            elif native_arch.value == 0xAA64:
                return 'arm64'
            elif native_arch.value == 0x01C4:
                return 'arm'
    except (OSError, AttributeError):
        pass
    # These env variables are always available. See:
    # https://msdn.microsoft.com/en-us/library/aa384274(VS.85).aspx
    # https://blogs.msdn.microsoft.com/david.wang/2006/03/27/howto-detect-process-bitness/
    arch = os.environ.get('PROCESSOR_ARCHITEW6432', '').lower()
    if not arch:
        try:
            # If this doesn't exist, something is messing with the environment
            arch = os.environ['PROCESSOR_ARCHITECTURE'].lower()
        except KeyError:
            raise EnvironmentException('Unable to detect native OS architecture')
    return arch

@dataclasses.dataclass
class VcsData:
    name: str
    cmd: str
    repo_dir: str
    get_rev: T.List[str]
    rev_regex: str
    dep: str
    wc_dir: T.Optional[str] = None
    repo_can_be_file: bool = False

    def repo_exists(self, curdir: Path) -> bool:
        if not shutil.which(self.cmd):
            return False

        repo = curdir / self.repo_dir
        if repo.is_dir():
            return True
        if repo.is_file() and self.repo_can_be_file:
            return True

        return False


def detect_vcs(source_dir: T.Union[str, Path]) -> T.Optional[VcsData]:
    vcs_systems = [
        VcsData(
            name = 'git',
            cmd = 'git',
            repo_dir = '.git',
            get_rev = ['git', 'describe', '--dirty=+', '--always'],
            rev_regex = '(.*)',
            dep = '.git/logs/HEAD',
            repo_can_be_file=True,
        ),
        VcsData(
            name = 'mercurial',
            cmd = 'hg',
            repo_dir = '.hg',
            get_rev = ['hg', 'id', '-i'],
            rev_regex = '(.*)',
            dep= '.hg/dirstate',
        ),
        VcsData(
            name = 'subversion',
            cmd = 'svn',
            repo_dir = '.svn',
            get_rev = ['svn', 'info'],
            rev_regex = 'Revision: (.*)',
            dep = '.svn/wc.db',
        ),
        VcsData(
            name = 'bazaar',
            cmd = 'bzr',
            repo_dir = '.bzr',
            get_rev = ['bzr', 'revno'],
            rev_regex = '(.*)',
            dep = '.bzr',
        ),
    ]
    if isinstance(source_dir, str):
        source_dir = Path(source_dir)

    parent_paths_and_self = collections.deque(source_dir.parents)
    # Prepend the source directory to the front so we can check it;
    # source_dir.parents doesn't include source_dir
    parent_paths_and_self.appendleft(source_dir)
    for curdir in parent_paths_and_self:
        for vcs in vcs_systems:
            if vcs.repo_exists(curdir):
                vcs.wc_dir = str(curdir)
                return vcs
    return None

def current_vs_supports_modules() -> bool:
    vsver = os.environ.get('VSCMD_VER', '')
    nums = vsver.split('.', 2)
    major = int(nums[0])
    if major >= 17:
        return True
    if major == 16 and int(nums[1]) >= 10:
        return True
    return vsver.startswith('16.9.0') and '-pre.' in vsver

_VERSION_TOK_RE = re.compile(r'(\d+)|([a-zA-Z]+)')

# a helper class which implements the same version ordering as RPM
class Version:
    def __init__(self, s: str) -> None:
        self._s = s

        # extract numeric and alphabetic sequences
        # numeric sequences are converted from strings to ints
        self._v = [
                int(m.group(1)) if m.group(1) else m.group(2)
                for m in _VERSION_TOK_RE.finditer(s)]

    def __str__(self) -> str:
        return '{} (V={})'.format(self._s, str(self._v))

    def __repr__(self) -> str:
        return f'<Version: {self._s}>'

    def __lt__(self, other: object) -> bool:
        if isinstance(other, Version):
            return self.__cmp(other, operator.lt)
        return NotImplemented

    def __gt__(self, other: object) -> bool:
        if isinstance(other, Version):
            return self.__cmp(other, operator.gt)
        return NotImplemented

    def __le__(self, other: object) -> bool:
        if isinstance(other, Version):
            return self.__cmp(other, operator.le)
        return NotImplemented

    def __ge__(self, other: object) -> bool:
        if isinstance(other, Version):
            return self.__cmp(other, operator.ge)
        return NotImplemented

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Version):
            return self._v == other._v
        return NotImplemented

    def __ne__(self, other: object) -> bool:
        if isinstance(other, Version):
            return self._v != other._v
        return NotImplemented

    def __cmp(self, other: 'Version', comparator: T.Callable[[T.Any, T.Any], bool]) -> bool:
        # compare each sequence in order
        for ours, theirs in zip(self._v, other._v):
            # sort a non-digit sequence before a digit sequence
            ours_is_int = isinstance(ours, int)
            theirs_is_int = isinstance(theirs, int)
            if ours_is_int != theirs_is_int:
                return comparator(ours_is_int, theirs_is_int)

            if ours != theirs:
                return comparator(ours, theirs)

        # if equal length, all components have matched, so equal
        # otherwise, the version with a suffix remaining is greater
        return comparator(len(self._v), len(other._v))


def _version_extract_cmpop(vstr2: str) -> T.Tuple[T.Callable[[T.Any, T.Any], bool], str]:
    if vstr2.startswith('>='):
        cmpop = operator.ge
        vstr2 = vstr2[2:]
    elif vstr2.startswith('<='):
        cmpop = operator.le
        vstr2 = vstr2[2:]
    elif vstr2.startswith('!='):
        cmpop = operator.ne
        vstr2 = vstr2[2:]
    elif vstr2.startswith('=='):
        cmpop = operator.eq
        vstr2 = vstr2[2:]
    elif vstr2.startswith('='):
        cmpop = operator.eq
        vstr2 = vstr2[1:]
    elif vstr2.startswith('>'):
        cmpop = operator.gt
        vstr2 = vstr2[1:]
    elif vstr2.startswith('<'):
        cmpop = operator.lt
        vstr2 = vstr2[1:]
    else:
        cmpop = operator.eq

    return (cmpop, vstr2)


def version_compare(vstr1: str, vstr2: str) -> bool:
    (cmpop, vstr2) = _version_extract_cmpop(vstr2)
    return cmpop(Version(vstr1), Version(vstr2))


def version_compare_many(vstr1: str, conditions: T.Union[str, T.Iterable[str]]) -> T.Tuple[bool, T.List[str], T.List[str]]:
    if isinstance(conditions, str):
        conditions = [conditions]
    found: T.List[str] = []
    not_found: T.List[str] = []
    for req in conditions:
        if not version_compare(vstr1, req):
            not_found.append(req)
        else:
            found.append(req)
    return not not_found, not_found, found


# determine if the minimum version satisfying the condition |condition| exceeds
# the minimum version for a feature |minimum|
def version_compare_condition_with_min(condition: str, minimum: str) -> bool:
    if condition.startswith('>='):
        cmpop = operator.le
        condition = condition[2:]
    elif condition.startswith('<='):
        return False
    elif condition.startswith('!='):
        return False
    elif condition.startswith('=='):
        cmpop = operator.le
        condition = condition[2:]
    elif condition.startswith('='):
        cmpop = operator.le
        condition = condition[1:]
    elif condition.startswith('>'):
        cmpop = operator.lt
        condition = condition[1:]
    elif condition.startswith('<'):
        return False
    else:
        cmpop = operator.le

    # Declaring a project(meson_version: '>=0.46') and then using features in
    # 0.46.0 is valid, because (knowing the meson versioning scheme) '0.46.0' is
    # the lowest version which satisfies the constraint '>=0.46'.
    #
    # But this will fail here, because the minimum version required by the
    # version constraint ('0.46') is strictly less (in our version comparison)
    # than the minimum version needed for the feature ('0.46.0').
    #
    # Map versions in the constraint of the form '0.46' to '0.46.0', to embed
    # this knowledge of the meson versioning scheme.
    condition = condition.strip()
    if re.match(r'^\d+.\d+$', condition):
        condition += '.0'

    return T.cast('bool', cmpop(Version(minimum), Version(condition)))

def search_version(text: str) -> str:
    # Usually of the type 4.1.4 but compiler output may contain
    # stuff like this:
    # (Sourcery CodeBench Lite 2014.05-29) 4.8.3 20140320 (prerelease)
    # Limiting major version number to two digits seems to work
    # thus far. When we get to GCC 100, this will break, but
    # if we are still relevant when that happens, it can be
    # considered an achievement in itself.
    #
    # This regex is reaching magic levels. If it ever needs
    # to be updated, do not complexify but convert to something
    # saner instead.
    # We'll demystify it a bit with a verbose definition.
    version_regex = re.compile(r"""
    (?<!                # Zero-width negative lookbehind assertion
        (
            \d          # One digit
            | \.        # Or one period
        )               # One occurrence
    )
    # Following pattern must not follow a digit or period
    (
        \d{1,2}         # One or two digits
        (
            \.\d+       # Period and one or more digits
        )+              # One or more occurrences
        (
            -[a-zA-Z0-9]+   # Hyphen and one or more alphanumeric
        )?              # Zero or one occurrence
    )                   # One occurrence
    """, re.VERBOSE)
    match = version_regex.search(text)
    if match:
        return match.group(0)

    # try a simpler regex that has like "blah 2020.01.100 foo" or "blah 2020.01 foo"
    version_regex = re.compile(r"(\d{1,4}\.\d{1,4}\.?\d{0,4})")
    match = version_regex.search(text)
    if match:
        return match.group(0)

    return 'unknown version'


def default_libdir() -> str:
    if is_debianlike():
        try:
            pc = subprocess.Popen(['dpkg-architecture', '-qDEB_HOST_MULTIARCH'],
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.DEVNULL)
            (stdo, _) = pc.communicate()
            if pc.returncode == 0:
                archpath = stdo.decode().strip()
                return 'lib/' + archpath
        except Exception:
            pass
    if is_freebsd() or is_irix():
        return 'lib'
    if os.path.isdir('/usr/lib64') and not os.path.islink('/usr/lib64'):
        return 'lib64'
    return 'lib'


def default_libexecdir() -> str:
    if is_haiku():
        return 'lib'
    # There is no way to auto-detect this, so it must be set at build time
    return 'libexec'


def default_prefix() -> str:
    if is_windows():
        return 'c:/'
    if is_haiku():
        return '/boot/system/non-packaged'
    return '/usr/local'


def default_datadir() -> str:
    if is_haiku():
        return 'data'
    return 'share'


def default_includedir() -> str:
    if is_haiku():
        return 'develop/headers'
    return 'include'


def default_infodir() -> str:
    if is_haiku():
        return 'documentation/info'
    return 'share/info'


def default_localedir() -> str:
    if is_haiku():
        return 'data/locale'
    return 'share/locale'


def default_mandir() -> str:
    if is_haiku():
        return 'documentation/man'
    return 'share/man'


def default_sbindir() -> str:
    if is_haiku():
        return 'bin'
    return 'sbin'


def default_sysconfdir() -> str:
    if is_haiku():
        return 'settings'
    return 'etc'


def determine_worker_count(varnames: T.Optional[T.List[str]] = None) -> int:
    num_workers = 0
    varnames = varnames or []
    # Add MESON_NUM_PROCESSES last, so it will prevail if more than one
    # variable is present.
    varnames.append('MESON_NUM_PROCESSES')
    for varname in varnames:
        if varname in os.environ:
            try:
                num_workers = int(os.environ[varname])
                if num_workers < 0:
                    raise ValueError
            except ValueError:
                print(f'Invalid value in {varname}, using 1 thread.')
                num_workers = 1

    if num_workers == 0:
        try:
            # Fails in some weird environments such as Debian
            # reproducible build.
            num_workers = multiprocessing.cpu_count()
        except Exception:
            num_workers = 1
    return num_workers

def is_parent_path(parent: str, trial: str) -> bool:
    '''Checks if @trial is a file under the directory @parent. Both @trial and @parent should be
       adequately normalized, though empty and '.' segments in @parent and @trial are accepted
       and discarded, matching the behavior of os.path.commonpath.  Either both or none should
       be absolute.'''
    assert os.path.isabs(parent) == os.path.isabs(trial)
    if is_windows():
        parent = parent.replace('\\', '/')
        trial = trial.replace('\\', '/')

    split_parent = parent.split('/')
    split_trial = trial.split('/')

    split_parent = [c for c in split_parent if c and c != '.']
    split_trial = [c for c in split_trial if c and c != '.']

    components = len(split_parent)
    return len(split_trial) >= components and split_trial[:components] == split_parent


def has_path_sep(name: str, sep: str = '/\\') -> bool:
    'Checks if any of the specified @sep path separators are in @name'
    for each in sep:
        if each in name:
            return True
    return False


if is_windows():
    # shlex.split is not suitable for splitting command line on Window (https://bugs.python.org/issue1724822);
    # shlex.quote is similarly problematic. Below are "proper" implementations of these functions according to
    # https://docs.microsoft.com/en-us/cpp/c-language/parsing-c-command-line-arguments and
    # https://blogs.msdn.microsoft.com/twistylittlepassagesallalike/2011/04/23/everyone-quotes-command-line-arguments-the-wrong-way/

    _whitespace = ' \t\n\r'
    _find_unsafe_char = re.compile(fr'[{_whitespace}"]').search

    def quote_arg(arg: str) -> str:
        if arg and not _find_unsafe_char(arg):
            return arg

        result = '"'
        num_backslashes = 0
        for c in arg:
            if c == '\\':
                num_backslashes += 1
            else:
                if c == '"':
                    # Escape all backslashes and the following double quotation mark
                    num_backslashes = num_backslashes * 2 + 1

                result += num_backslashes * '\\' + c
                num_backslashes = 0

        # Escape all backslashes, but let the terminating double quotation
        # mark we add below be interpreted as a metacharacter
        result += (num_backslashes * 2) * '\\' + '"'
        return result

    def split_args(cmd: str) -> T.List[str]:
        result: T.List[str] = []
        arg = ''
        num_backslashes = 0
        num_quotes = 0
        in_quotes = False
        for c in cmd:
            if c == '\\':
                num_backslashes += 1
            else:
                if c == '"' and not num_backslashes % 2:
                    # unescaped quote, eat it
                    arg += (num_backslashes // 2) * '\\'
                    num_quotes += 1
                    in_quotes = not in_quotes
                elif c in _whitespace and not in_quotes:
                    if arg or num_quotes:
                        # reached the end of the argument
                        result.append(arg)
                        arg = ''
                        num_quotes = 0
                else:
                    if c == '"':
                        # escaped quote
                        num_backslashes = (num_backslashes - 1) // 2

                    arg += num_backslashes * '\\' + c

                num_backslashes = 0

        if arg or num_quotes:
            result.append(arg)

        return result
else:
    def quote_arg(arg: str) -> str:
        return shlex.quote(arg)

    def split_args(cmd: str) -> T.List[str]:
        return shlex.split(cmd)


def join_args(args: T.Iterable[str]) -> str:
    return ' '.join([quote_arg(x) for x in args])


def do_replacement(regex: T.Pattern[str], line: str,
                   variable_format: Literal['meson', 'cmake', 'cmake@'],
                   confdata: T.Union[T.Dict[str, T.Tuple[str, T.Optional[str]]], 'ConfigurationData']) -> T.Tuple[str, T.Set[str]]:
    if variable_format == 'meson':
        return do_replacement_meson(regex, line, confdata)
    elif variable_format in {'cmake', 'cmake@'}:
        return do_replacement_cmake(line, variable_format == 'cmake@', confdata)
    else:
        raise MesonException('Invalid variable format')

def do_replacement_meson(regex: T.Pattern[str], line: str,
                         confdata: T.Union[T.Dict[str, T.Tuple[str, T.Optional[str]]], 'ConfigurationData']) -> T.Tuple[str, T.Set[str]]:
    missing_variables: T.Set[str] = set()

    def variable_replace(match: T.Match[str]) -> str:
        # Pairs of escape characters before '@', '\@', '${' or '\${'
        if match.group(0).endswith('\\'):
            num_escapes = match.end(0) - match.start(0)
            return '\\' * (num_escapes // 2)
        # \@escaped\@ variables
        elif match.groupdict().get('escaped') is not None:
            return match.group('escaped')[1:-2]+'@'
        else:
            # Template variable to be replaced
            varname = match.group('variable')
            var_str = ''
            if varname in confdata:
                var, _ = confdata.get(varname)
                if isinstance(var, str):
                    var_str = var
                elif isinstance(var, int):
                    var_str = str(var)
                else:
                    msg = f'Tried to replace variable {varname!r} value with ' \
                          f'something other than a string or int: {var!r}'
                    raise MesonException(msg)
            else:
                missing_variables.add(varname)
            return var_str
    return re.sub(regex, variable_replace, line), missing_variables

def do_replacement_cmake(line: str, at_only: bool,
                         confdata: T.Union[T.Dict[str, T.Tuple[str, T.Optional[str]]], 'ConfigurationData']) -> T.Tuple[str, T.Set[str]]:
    missing_variables: T.Set[str] = set()

    character_regex = re.compile(r'''
        [^a-zA-Z0-9_/.+\-]
    ''', re.VERBOSE)

    def variable_get(varname: str) -> str:
        var_str = ''
        if varname in confdata:
            var, _ = confdata.get(varname)
            if isinstance(var, str):
                var_str = var
            elif isinstance(var, bool):
                var_str = str(int(var))
            elif isinstance(var, int):
                var_str = str(var)
            else:
                msg = f'Tried to replace variable {varname!r} value with ' \
                      f'something other than a string or int: {var!r}'
                raise MesonException(msg)
        else:
            missing_variables.add(varname)
        return var_str

    def parse_line(line: str) -> str:
        index = 0
        while len(line) > index:
            if line[index] == '@':
                next_at = line.find("@", index+1)
                if next_at > index+1:
                    varname = line[index+1:next_at]
                    match = character_regex.search(varname)

                    # at substituion doesn't occur if they key isn't valid
                    # however it also doesn't raise an error
                    if not match:
                        value = variable_get(varname)
                        line = line[:index] + value + line[next_at+1:]

            elif not at_only and line[index:index+2] == '${':
                bracket_count = 1
                end_bracket = index + 2
                try:
                    while bracket_count > 0:
                        if line[end_bracket:end_bracket+2] == "${":
                            end_bracket += 2
                            bracket_count += 1
                        elif line[end_bracket] == "}":
                            end_bracket += 1
                            bracket_count -= 1
                        elif line[end_bracket] in {"@", "\n"}:
                            # these aren't valid variable characters
                            # but they are inconsequential at this point
                            end_bracket += 1
                        elif character_regex.search(line[end_bracket]):
                            invalid_character = line[end_bracket]
                            variable = line[index+2:end_bracket]
                            msg = f'Found invalid character {invalid_character!r}' \
                                  f' in variable {variable!r}'
                            raise MesonException(msg)
                        else:
                            end_bracket += 1
                except IndexError:
                    msg = f'Found incomplete variable {line[index:-1]!r}'
                    raise MesonException(msg)

                if bracket_count == 0:
                    varname = parse_line(line[index+2:end_bracket-1])
                    match = character_regex.search(varname)
                    if match:
                        invalid_character = line[end_bracket-2]
                        variable = line[index+2:end_bracket-3]
                        msg = f'Found invalid character {invalid_character!r}' \
                              f' in variable {variable!r}'
                        raise MesonException(msg)

                    value = variable_get(varname)
                    line = line[:index] + value + line[end_bracket:]

            index += 1

        return line

    return parse_line(line), missing_variables

def do_define_meson(regex: T.Pattern[str], line: str, confdata: 'ConfigurationData',
                    subproject: T.Optional[SubProject] = None) -> str:

    arr = line.split()
    if len(arr) != 2:
        raise MesonException('#mesondefine does not contain exactly two tokens: %s' % line.strip())

    varname = arr[1]
    try:
        v, _ = confdata.get(varname)
    except KeyError:
        return '/* #undef %s */\n' % varname

    if isinstance(v, str):
        result = f'#define {varname} {v}'.strip() + '\n'
        result, _ = do_replacement_meson(regex, result, confdata)
        return result
    elif isinstance(v, bool):
        if v:
            return '#define %s\n' % varname
        else:
            return '#undef %s\n' % varname
    elif isinstance(v, int):
        return '#define %s %d\n' % (varname, v)
    else:
        raise MesonException('#mesondefine argument "%s" is of unknown type.' % varname)

def do_define_cmake(line: str, confdata: 'ConfigurationData', at_only: bool,
                    subproject: T.Optional[SubProject] = None) -> str:
    cmake_bool_define = 'cmakedefine01' in line

    def get_cmake_define(line: str, confdata: 'ConfigurationData') -> str:
        arr = line[1:].split()

        if cmake_bool_define:
            (v, desc) = confdata.get(arr[1])
            return str(int(bool(v)))

        define_value: T.List[str] = []
        for token in arr[2:]:
            try:
                v, _ = confdata.get(token)
                define_value += [str(v)]
            except KeyError:
                define_value += [token]
        return ' '.join(define_value)

    arr = line[1:].split()

    if len(arr) != 2 and subproject is not None:
        from ..interpreterbase.decorators import FeatureNew
        FeatureNew.single_use('cmakedefine without exactly two tokens', '0.54.1', subproject)

    varname = arr[1]
    try:
        v, _ = confdata.get(varname)
    except KeyError:
        if cmake_bool_define:
            return '#define %s 0\n' % varname
        else:
            return '/* #undef %s */\n' % varname

    if not cmake_bool_define and not v:
        return '/* #undef %s */\n' % varname

    result = get_cmake_define(line, confdata)
    result = f'#define {varname} {result}'.strip() + '\n'
    result, _ = do_replacement_cmake(result, at_only, confdata)
    return result

def get_variable_regex(variable_format: Literal['meson', 'cmake', 'cmake@'] = 'meson') -> T.Pattern[str]:
    # Only allow (a-z, A-Z, 0-9, _, -) as valid characters for a define
    if variable_format == 'meson':
        # Also allow escaping pairs of '@' with '\@'
        regex = re.compile(r'''
            (?:\\\\)+(?=\\?@)  # Matches multiple backslashes followed by an @ symbol
            |                  # OR
            (?<!\\)@(?P<variable>[-a-zA-Z0-9_]+)@  # Match a variable enclosed in @ symbols and capture the variable name; no matches beginning with '\@'
            |                  # OR
            (?P<escaped>\\@[-a-zA-Z0-9_]+\\@)  # Match an escaped variable enclosed in @ symbols
        ''', re.VERBOSE)
    elif variable_format == 'cmake@':
        regex = re.compile(r'''
            (?<!\\)@(?P<variable>[-a-zA-Z0-9_]+)@  # Match a variable enclosed in @ symbols and capture the variable name; no matches beginning with '\@'
        ''', re.VERBOSE)
    elif variable_format == "cmake":
        regex = re.compile(r'''
            \${(?P<variable>[-a-zA-Z0-9_]*)}  # Match a variable enclosed in curly braces and capture the variable name
        ''', re.VERBOSE)
    return regex

def do_conf_str(src: str, data: T.List[str], confdata: 'ConfigurationData',
                variable_format: Literal['meson', 'cmake', 'cmake@'],
                subproject: T.Optional[SubProject] = None) -> T.Tuple[T.List[str], T.Set[str], bool]:
    if variable_format == 'meson':
        return do_conf_str_meson(src, data, confdata, subproject)
    elif variable_format in {'cmake', 'cmake@'}:
        return do_conf_str_cmake(src, data, confdata, variable_format == 'cmake@', subproject)
    else:
        raise MesonException('Invalid variable format')

def do_conf_str_meson(src: str, data: T.List[str], confdata: 'ConfigurationData',
                      subproject: T.Optional[SubProject] = None) -> T.Tuple[T.List[str], T.Set[str], bool]:

    regex = get_variable_regex('meson')

    search_token = '#mesondefine'

    result: T.List[str] = []
    missing_variables: T.Set[str] = set()
    # Detect when the configuration data is empty and no tokens were found
    # during substitution so we can warn the user to use the `copy:` kwarg.
    confdata_useless = not confdata.keys()
    for line in data:
        if line.lstrip().startswith(search_token):
            confdata_useless = False
            line = do_define_meson(regex, line, confdata, subproject)
        else:
            if '#cmakedefine' in line:
                raise MesonException(f'Format error in {src}: saw "{line.strip()}" when format set to "meson"')
            line, missing = do_replacement_meson(regex, line, confdata)
            missing_variables.update(missing)
            if missing:
                confdata_useless = False
        result.append(line)

    return result, missing_variables, confdata_useless

def do_conf_str_cmake(src: str, data: T.List[str], confdata: 'ConfigurationData', at_only: bool,
                      subproject: T.Optional[SubProject] = None) -> T.Tuple[T.List[str], T.Set[str], bool]:

    variable_format: Literal['cmake', 'cmake@'] = 'cmake'
    if at_only:
        variable_format = 'cmake@'

    search_token = 'cmakedefine'

    result: T.List[str] = []
    missing_variables: T.Set[str] = set()
    # Detect when the configuration data is empty and no tokens were found
    # during substitution so we can warn the user to use the `copy:` kwarg.
    confdata_useless = not confdata.keys()
    for line in data:
        stripped_line = line.lstrip()
        if len(stripped_line) >= 2 and stripped_line[0] == '#' and stripped_line[1:].lstrip().startswith(search_token):
            confdata_useless = False

            line = do_define_cmake(line, confdata, at_only, subproject)
        else:
            if '#mesondefine' in line:
                raise MesonException(f'Format error in {src}: saw "{line.strip()}" when format set to "{variable_format}"')
            line, missing = do_replacement_cmake(line, at_only, confdata)
            missing_variables.update(missing)
            if missing:
                confdata_useless = False
        result.append(line)

    return result, missing_variables, confdata_useless

def do_conf_file(src: str, dst: str, confdata: 'ConfigurationData',
                 variable_format: Literal['meson', 'cmake', 'cmake@'],
                 encoding: str = 'utf-8', subproject: T.Optional[SubProject] = None) -> T.Tuple[T.Set[str], bool]:
    try:
        with open(src, encoding=encoding, newline='') as f:
            data = f.readlines()
    except Exception as e:
        raise MesonException(f'Could not read input file {src}: {e!s}')

    (result, missing_variables, confdata_useless) = do_conf_str(src, data, confdata, variable_format, subproject)
    dst_tmp = dst + '~'
    try:
        with open(dst_tmp, 'w', encoding=encoding, newline='') as f:
            f.writelines(result)
    except Exception as e:
        raise MesonException(f'Could not write output file {dst}: {e!s}')
    shutil.copymode(src, dst_tmp)
    replace_if_different(dst, dst_tmp)
    return missing_variables, confdata_useless

CONF_C_PRELUDE = '''/*
 * Autogenerated by the Meson build system.
 * Do not edit, your changes will be lost.
 */

{}

'''

CONF_NASM_PRELUDE = '''; Autogenerated by the Meson build system.
; Do not edit, your changes will be lost.

'''

def _dump_c_header(ofile: T.TextIO,
                   cdata: ConfigurationData,
                   output_format: Literal['c', 'nasm'],
                   macro_name: T.Optional[str]) -> None:
    format_desc: T.Callable[[str], str]
    if output_format == 'c':
        if macro_name:
            prelude = CONF_C_PRELUDE.format('#ifndef {0}\n#define {0}'.format(macro_name))
        else:
            prelude = CONF_C_PRELUDE.format('#pragma once')
        prefix = '#'
        format_desc = lambda desc: f'/* {desc} */\n'
    else:  # nasm
        prelude = CONF_NASM_PRELUDE
        prefix = '%'
        format_desc = lambda desc: '; ' + '\n; '.join(desc.splitlines()) + '\n'

    ofile.write(prelude)
    for k in sorted(cdata.keys()):
        (v, desc) = cdata.get(k)
        if desc:
            ofile.write(format_desc(desc))
        if isinstance(v, bool):
            if v:
                ofile.write(f'{prefix}define {k}\n\n')
            else:
                ofile.write(f'{prefix}undef {k}\n\n')
        elif isinstance(v, (int, str)):
            ofile.write(f'{prefix}define {k} {v}\n\n')
        else:
            raise MesonException('Unknown data type in configuration file entry: ' + k)
    if output_format == 'c' and macro_name:
        ofile.write('#endif\n')


def dump_conf_header(ofilename: str, cdata: ConfigurationData,
                     output_format: Literal['c', 'nasm', 'json'],
                     macro_name: T.Optional[str]) -> None:
    ofilename_tmp = ofilename + '~'
    with open(ofilename_tmp, 'w', encoding='utf-8') as ofile:
        if output_format == 'json':
            data = {k: v[0] for k, v in cdata.values.items()}
            json.dump(data, ofile, sort_keys=True)
        else:  # c, nasm
            _dump_c_header(ofile, cdata, output_format, macro_name)

    replace_if_different(ofilename, ofilename_tmp)


def replace_if_different(dst: str, dst_tmp: str) -> None:
    # If contents are identical, don't touch the file to prevent
    # unnecessary rebuilds.
    different = True
    try:
        with open(dst, 'rb') as f1, open(dst_tmp, 'rb') as f2:
            if f1.read() == f2.read():
                different = False
    except FileNotFoundError:
        pass
    if different:
        os.replace(dst_tmp, dst)
    else:
        os.unlink(dst_tmp)


def listify(item: T.Any, flatten: bool = True) -> T.List[T.Any]:
    '''
    Returns a list with all args embedded in a list if they are not a list.
    This function preserves order.
    @flatten: Convert lists of lists to a flat list
    '''
    if not isinstance(item, list):
        return [item]
    result: T.List[T.Any] = []
    for i in item:
        if flatten and isinstance(i, list):
            result += listify(i, flatten=True)
        else:
            result.append(i)
    return result

def listify_array_value(value: object, shlex_split_args: bool = False) -> T.List[str]:
    if isinstance(value, str):
        if value.startswith('['):
            try:
                newvalue = ast.literal_eval(value)
            except ValueError:
                raise MesonException(f'malformed value {value}')
        elif value == '':
            newvalue = []
        else:
            if shlex_split_args:
                newvalue = split_args(value)
            else:
                newvalue = [v.strip() for v in value.split(',')]
    elif isinstance(value, list):
        newvalue = value
    else:
        raise MesonException(f'"{value}" should be a string array, but it is not')
    assert isinstance(newvalue, list)
    return newvalue

def extract_as_list(dict_object: T.Dict[_T, _U], key: _T, pop: bool = False) -> T.List[_U]:
    '''
    Extracts all values from given dict_object and listifies them.
    '''
    fetch: T.Callable[[_T], _U] = dict_object.get
    if pop:
        fetch = dict_object.pop
    # If there's only one key, we don't return a list with one element
    return listify(fetch(key) or [], flatten=True)


def typeslistify(item: 'T.Union[_T, T.Sequence[_T]]',
                 types: 'T.Union[T.Type[_T], T.Tuple[T.Type[_T]]]') -> T.List[_T]:
    '''
    Ensure that type(@item) is one of @types or a
    list of items all of which are of type @types
    '''
    if isinstance(item, types):
        item = T.cast('T.List[_T]', [item])
    if not isinstance(item, list):
        raise MesonException('Item must be a list or one of {!r}, not {!r}'.format(types, type(item)))
    for i in item:
        if i is not None and not isinstance(i, types):
            raise MesonException('List item must be one of {!r}, not {!r}'.format(types, type(i)))
    return item


def stringlistify(item: T.Union[T.Any, T.Sequence[T.Any]]) -> T.List[str]:
    return typeslistify(item, str)


def expand_arguments(args: T.Iterable[str]) -> T.Optional[T.List[str]]:
    expended_args: T.List[str] = []
    for arg in args:
        if not arg.startswith('@'):
            expended_args.append(arg)
            continue

        args_file = arg[1:]
        try:
            with open(args_file, encoding='utf-8') as f:
                extended_args = f.read().split()
            expended_args += extended_args
        except Exception as e:
            mlog.error('Expanding command line arguments:',  args_file, 'not found')
            mlog.exception(e)
            return None
    return expended_args


def partition(pred: T.Callable[[_T], object], iterable: T.Iterable[_T]) -> T.Tuple[T.Iterator[_T], T.Iterator[_T]]:
    """Use a predicate to partition entries into false entries and true
    entries.

    >>> x, y = partition(is_odd, range(10))
    >>> (list(x), list(y))
    ([0, 2, 4, 6, 8], [1, 3, 5, 7, 9])
    """
    t1, t2 = tee(iterable)
    return (t for t in t1 if not pred(t)), (t for t in t2 if pred(t))


def Popen_safe(args: T.List[str], write: T.Optional[str] = None,
               stdin: T.Union[None, T.TextIO, T.BinaryIO, int] = subprocess.DEVNULL,
               stdout: T.Union[None, T.TextIO, T.BinaryIO, int] = subprocess.PIPE,
               stderr: T.Union[None, T.TextIO, T.BinaryIO, int] = subprocess.PIPE,
               **kwargs: T.Any) -> T.Tuple['subprocess.Popen[str]', str, str]:
    import locale
    encoding = locale.getpreferredencoding()
    # Stdin defaults to DEVNULL otherwise the command run by us here might mess
    # up the console and ANSI colors will stop working on Windows.
    # If write is not None, set stdin to PIPE so data can be sent.
    if write is not None:
        stdin = subprocess.PIPE

    try:
        if not sys.stdout.encoding or encoding.upper() != 'UTF-8':
            p, o, e = Popen_safe_legacy(args, write=write, stdin=stdin, stdout=stdout, stderr=stderr, **kwargs)
        else:
            p = subprocess.Popen(args, universal_newlines=True, encoding=encoding, close_fds=False,
                                 stdin=stdin, stdout=stdout, stderr=stderr, **kwargs)
            o, e = p.communicate(write)
    except OSError as oserr:
        if oserr.errno == errno.ENOEXEC:
            raise MesonException(f'Failed running {args[0]!r}, binary or interpreter not executable.\n'
                                 'Possibly wrong architecture or the executable bit is not set.')
        raise
    # Sometimes the command that we run will call another command which will be
    # without the above stdin workaround, so set the console mode again just in
    # case.
    mlog.setup_console()
    return p, o, e


def Popen_safe_legacy(args: T.List[str], write: T.Optional[str] = None,
                      stdin: T.Union[None, T.TextIO, T.BinaryIO, int] = subprocess.DEVNULL,
                      stdout: T.Union[None, T.TextIO, T.BinaryIO, int] = subprocess.PIPE,
                      stderr: T.Union[None, T.TextIO, T.BinaryIO, int] = subprocess.PIPE,
                      **kwargs: T.Any) -> T.Tuple['subprocess.Popen[str]', str, str]:
    p = subprocess.Popen(args, universal_newlines=False, close_fds=False,
                         stdin=stdin, stdout=stdout, stderr=stderr, **kwargs)
    input_: T.Optional[bytes] = None
    if write is not None:
        input_ = write.encode('utf-8')
    o, e = p.communicate(input_)
    if o is not None:
        if sys.stdout.encoding is not None:
            o = o.decode(encoding=sys.stdout.encoding, errors='replace').replace('\r\n', '\n')
        else:
            o = o.decode(errors='replace').replace('\r\n', '\n')
    if e is not None:
        if sys.stderr is not None and sys.stderr.encoding:
            e = e.decode(encoding=sys.stderr.encoding, errors='replace').replace('\r\n', '\n')
        else:
            e = e.decode(errors='replace').replace('\r\n', '\n')
    return p, o, e


def Popen_safe_logged(args: T.List[str], msg: str = 'Called', **kwargs: T.Any) -> T.Tuple['subprocess.Popen[str]', str, str]:
    '''
    Wrapper around Popen_safe that assumes standard piped o/e and logs this to the meson log.
    '''
    try:
        p, o, e = Popen_safe(args, **kwargs)
    except Exception as excp:
        mlog.debug('-----------')
        mlog.debug(f'{msg}: `{join_args(args)}` -> {excp}')
        raise

    rc, out, err = p.returncode, o.strip() if o else None, e.strip() if e else None
    mlog.debug('-----------')
    mlog.debug(f'{msg}: `{join_args(args)}` -> {rc}')
    if out:
        mlog.debug(f'stdout:\n{out}\n-----------')
    if err:
        mlog.debug(f'stderr:\n{err}\n-----------')
    return p, o, e


def iter_regexin_iter(regexiter: T.Iterable[str], initer: T.Iterable[str | programs.ExternalProgram]) -> T.Optional[str]:
    '''
    Takes each regular expression in @regexiter and tries to search for it in
    every item in @initer. If there is a match, returns that match.
    Else returns False.
    '''
    for regex in regexiter:
        for ii in initer:
            if not isinstance(ii, str):
                continue
            match = re.search(regex, ii)
            if match:
                return match.group()
    return None


def _substitute_values_check_errors(command: T.List[str | programs.ExternalProgram], values: T.Dict[str, T.Union[str, T.List[str]]]) -> None:
    # Error checking
    inregex: T.List[str] = ['@INPUT([0-9]+)?@', '@PLAINNAME@', '@BASENAME@']
    outregex: T.List[str] = ['@OUTPUT([0-9]+)?@', '@OUTDIR@']
    if '@INPUT@' not in values:
        # Error out if any input-derived templates are present in the command
        match = iter_regexin_iter(inregex, command)
        if match:
            raise MesonException(f'Command cannot have {match!r}, since no input files were specified')
    else:
        if len(values['@INPUT@']) > 1:
            # Error out if @PLAINNAME@ or @BASENAME@ is present in the command
            match = iter_regexin_iter(inregex[1:], command)
            if match:
                raise MesonException(f'Command cannot have {match!r} when there is '
                                     'more than one input file')
        # Error out if an invalid @INPUTnn@ template was specified
        for each in command:
            if not isinstance(each, str):
                continue
            match2 = re.search(inregex[0], each)
            if match2 and match2.group() not in values:
                m = 'Command cannot have {!r} since there are only {!r} inputs'
                raise MesonException(m.format(match2.group(), len(values['@INPUT@'])))
    if '@OUTPUT@' not in values:
        # Error out if any output-derived templates are present in the command
        match = iter_regexin_iter(outregex, command)
        if match:
            raise MesonException(f'Command cannot have {match!r} since there are no outputs')
    else:
        # Error out if an invalid @OUTPUTnn@ template was specified
        for each in command:
            if not isinstance(each, str):
                continue
            match2 = re.search(outregex[0], each)
            if match2 and match2.group() not in values:
                m = 'Command cannot have {!r} since there are only {!r} outputs'
                raise MesonException(m.format(match2.group(), len(values['@OUTPUT@'])))


def substitute_values(command: T.List[str | programs.ExternalProgram], values: T.Dict[str, T.Union[str, T.List[str]]]) -> T.List[str | programs.ExternalProgram]:
    '''
    Substitute the template strings in the @values dict into the list of
    strings @command and return a new list. For a full list of the templates,
    see get_filenames_templates_dict()

    If multiple inputs/outputs are given in the @values dictionary, we
    substitute @INPUT@ and @OUTPUT@ only if they are the entire string, not
    just a part of it, and in that case we substitute *all* of them.

    The typing of this function is difficult, as only @OUTPUT@ and @INPUT@ can
    be lists, everything else is a string. However, TypeDict cannot represent
    this, as you can have optional keys, but not extra keys. We end up just
    having to us asserts to convince type checkers that this is okay.

    https://github.com/python/mypy/issues/4617
    '''

    def replace(m: T.Match[str]) -> str:
        v = values[m.group(0)]
        assert isinstance(v, str), 'for mypy'
        return v

    # Error checking
    _substitute_values_check_errors(command, values)

    # Substitution
    outcmd: T.List[str | programs.ExternalProgram] = []
    rx_keys = [re.escape(key) for key in values if key not in ('@INPUT@', '@OUTPUT@')]
    value_rx = re.compile('|'.join(rx_keys)) if rx_keys else None
    for vv in command:
        more: T.Optional[str] = None
        if not isinstance(vv, str):
            outcmd.append(vv)
        elif '@INPUT@' in vv:
            inputs = values['@INPUT@']
            if vv == '@INPUT@':
                outcmd += inputs
            elif len(inputs) == 1:
                outcmd.append(vv.replace('@INPUT@', inputs[0]))
            else:
                raise MesonException("Command has '@INPUT@' as part of a "
                                     "string and more than one input file")
        elif '@OUTPUT@' in vv:
            outputs = values['@OUTPUT@']
            if vv == '@OUTPUT@':
                outcmd += outputs
            elif len(outputs) == 1:
                outcmd.append(vv.replace('@OUTPUT@', outputs[0]))
            else:
                raise MesonException("Command has '@OUTPUT@' as part of a "
                                     "string and more than one output file")

        # Append values that are exactly a template string.
        # This is faster than a string replace.
        elif vv in values:
            o = values[vv]
            assert isinstance(o, str), 'for mypy'
            more = o
        # Substitute everything else with replacement
        elif value_rx:
            more = value_rx.sub(replace, vv)
        else:
            more = vv

        if more is not None:
            outcmd.append(more)

    return outcmd


def get_filenames_templates_dict(inputs: T.List[str], outputs: T.List[str]) -> T.Dict[str, T.Union[str, T.List[str]]]:
    '''
    Create a dictionary with template strings as keys and values as values for
    the following templates:

    @INPUT@  - the full path to one or more input files, from @inputs
    @OUTPUT@ - the full path to one or more output files, from @outputs
    @OUTDIR@ - the full path to the directory containing the output files

    If there is only one input file, the following keys are also created:

    @PLAINNAME@ - the filename of the input file
    @BASENAME@ - the filename of the input file with the extension removed

    If there is more than one input file, the following keys are also created:

    @INPUT0@, @INPUT1@, ... one for each input file
    @PLAINNAME0@, @PLAINNAME1@, ... one for each input file
    @BASENAME0@, @BASENAME1@, ... one for each input file

    If there is more than one output file, the following keys are also created:

    @OUTPUT0@, @OUTPUT1@, ... one for each output file
    '''
    values: T.Dict[str, T.Union[str, T.List[str]]] = {}
    # Gather values derived from the input
    if inputs:
        # We want to substitute all the inputs.
        values['@INPUT@'] = inputs
        for (ii, vv) in enumerate(inputs):
            # Write out @INPUT0@, @INPUT1@, ...
            values[f'@INPUT{ii}@'] = vv
            plain = os.path.basename(vv)
            values[f'@PLAINNAME{ii}@'] = plain
            values[f'@BASENAME{ii}@'] = os.path.splitext(plain)[0]
        if len(inputs) == 1:
            # Just one value, substitute @PLAINNAME@ and @BASENAME@
            values['@PLAINNAME@'] = plain = os.path.basename(inputs[0])
            values['@BASENAME@'] = os.path.splitext(plain)[0]
    if outputs:
        # Gather values derived from the outputs, similar to above.
        values['@OUTPUT@'] = outputs
        for (ii, vv) in enumerate(outputs):
            values[f'@OUTPUT{ii}@'] = vv
        # Outdir should be the same for all outputs
        values['@OUTDIR@'] = os.path.dirname(outputs[0])
        # Many external programs fail on empty arguments.
        if values['@OUTDIR@'] == '':
            values['@OUTDIR@'] = '.'
    return values


def _make_tree_writable(topdir: T.Union[str, Path]) -> None:
    # Ensure all files and directories under topdir are writable
    # (and readable) by owner.
    for d, _, files in os.walk(topdir):
        os.chmod(d, os.stat(d).st_mode | stat.S_IWRITE | stat.S_IREAD)
        for fname in files:
            fpath = os.path.join(d, fname)
            if not os.path.islink(fpath) and os.path.isfile(fpath):
                os.chmod(fpath, os.stat(fpath).st_mode | stat.S_IWRITE | stat.S_IREAD)


def windows_proof_rmtree(f:  T.Union[str, Path]) -> None:
    # On Windows if anyone is holding a file open you can't
    # delete it. As an example an anti virus scanner might
    # be scanning files you are trying to delete. The only
    # way to fix this is to try again and again.
    delays = [0.1, 0.1, 0.2, 0.2, 0.2, 0.5, 0.5, 1, 1, 1, 1, 2]
    writable = False
    for d in delays:
        try:
            # Start by making the tree writable.
            if not writable:
                _make_tree_writable(f)
                writable = True
        except PermissionError:
            time.sleep(d)
            continue
        try:
            shutil.rmtree(f)
            return
        except FileNotFoundError:
            return
        except OSError:
            time.sleep(d)
    # Try one last time and throw if it fails.
    shutil.rmtree(f)


def windows_proof_rm(fpath: T.Union[str, Path]) -> None:
    """Like windows_proof_rmtree, but for a single file."""
    if os.path.isfile(fpath):
        os.chmod(fpath, os.stat(fpath).st_mode | stat.S_IWRITE | stat.S_IREAD)
    delays = [0.1, 0.1, 0.2, 0.2, 0.2, 0.5, 0.5, 1, 1, 1, 1, 2]
    for d in delays:
        try:
            os.unlink(fpath)
            return
        except FileNotFoundError:
            return
        except OSError:
            time.sleep(d)
    os.unlink(fpath)


class TemporaryDirectoryWinProof(TemporaryDirectory):
    """
    Like TemporaryDirectory, but cleans things up using
    windows_proof_rmtree()
    """

    def __exit__(self, exc: T.Any, value: T.Any, tb: T.Any) -> None:
        try:
            super().__exit__(exc, value, tb)
        except OSError:
            windows_proof_rmtree(self.name)

    def cleanup(self) -> None:
        try:
            super().cleanup()
        except OSError:
            windows_proof_rmtree(self.name)


def detect_subprojects(spdir_name: str, current_dir: str = '',
                       result: T.Optional[T.Dict[str, T.List[str]]] = None) -> T.Dict[str, T.List[str]]:
    if result is None:
        result = {}
    spdir = os.path.join(current_dir, spdir_name)
    if not os.path.exists(spdir):
        return result
    for trial in glob(os.path.join(spdir, '*')):
        basename = os.path.basename(trial)
        if trial == 'packagecache':
            continue
        append_this = True
        if os.path.isdir(trial):
            detect_subprojects(spdir_name, trial, result)
        elif trial.endswith('.wrap') and os.path.isfile(trial):
            basename = os.path.splitext(basename)[0]
        else:
            append_this = False
        if append_this:
            if basename in result:
                result[basename].append(trial)
            else:
                result[basename] = [trial]
    return result


def substring_is_in_list(substr: str, strlist: T.List[str]) -> bool:
    for s in strlist:
        if substr in s:
            return True
    return False


class OrderedSet(T.MutableSet[_T]):
    """A set that preserves the order in which items are added, by first
    insertion.
    """
    def __init__(self, iterable: T.Optional[T.Iterable[_T]] = None):
        self.__container: T.OrderedDict[_T, None] = collections.OrderedDict()
        if iterable:
            self.update(iterable)

    def __contains__(self, value: object) -> bool:
        return value in self.__container

    def __iter__(self) -> T.Iterator[_T]:
        return iter(self.__container.keys())

    def __len__(self) -> int:
        return len(self.__container)

    def __repr__(self) -> str:
        # Don't print 'OrderedSet("")' for an empty set.
        if self.__container:
            return 'OrderedSet([{}])'.format(
                ', '.join(repr(e) for e in self.__container.keys()))
        return 'OrderedSet()'

    def __reversed__(self) -> T.Iterator[_T]:
        return reversed(self.__container.keys())

    def add(self, value: _T) -> None:
        self.__container[value] = None

    def discard(self, value: _T) -> None:
        if value in self.__container:
            del self.__container[value]

    def move_to_end(self, value: _T, last: bool = True) -> None:
        self.__container.move_to_end(value, last)

    def pop(self, last: bool = True) -> _T:
        item, _ = self.__container.popitem(last)
        return item

    def update(self, iterable: T.Iterable[_T]) -> None:
        for item in iterable:
            self.__container[item] = None

    def difference(self, set_: T.Iterable[_T]) -> 'OrderedSet[_T]':
        return type(self)(e for e in self if e not in set_)

    def difference_update(self, iterable: T.Iterable[_T]) -> None:
        for item in iterable:
            self.discard(item)

def relpath(path: T.Union[str, Path], start: T.Union[str, Path]) -> str:
    # On Windows a relative path can't be evaluated for paths on two different
    # drives (i.e. c:\foo and f:\bar).  The only thing left to do is to use the
    # original absolute path.
    try:
        return os.path.relpath(path, start)
    except (TypeError, ValueError):
        return str(path)

def path_is_in_root(path: Path, root: Path, resolve: bool = False) -> bool:
    # Check whether a path is within the root directory root
    try:
        if resolve:
            path.resolve().relative_to(root.resolve())
        else:
            path.relative_to(root)
    except ValueError:
        return False
    return True

def relative_to_if_possible(path: Path, root: Path, resolve: bool = False) -> Path:
    try:
        if resolve:
            return path.resolve().relative_to(root.resolve())
        else:
            return path.relative_to(root)
    except ValueError:
        return path

class LibType(enum.IntEnum):

    """Enumeration for library types."""

    SHARED = 0
    STATIC = 1
    PREFER_SHARED = 2
    PREFER_STATIC = 3


class ProgressBarFallback:  # lgtm [py/iter-returns-non-self]
    '''
    Fallback progress bar implementation when tqdm is not foundclass OptionType(enum.IntEnum):

    """Enum used to specify what kind of argument a thing is."""

    BUILTIN = 0
    BACKEND = 1
    BASE = 2
    COMPILER = 3
    PROJECT = 4

# This is copied from coredata. There is no way to share this, because this
# is used in the OptionKey constructor, and the coredata lists are
# OptionKeys...
_BUILTIN_NAMES = {
    'prefix',
    'bindir',
    'datadir',
    'includedir',
    'infodir',
    'libdir',
    'licensedir',
    'libexecdir',
    'localedir',
    'localstatedir',
    'mandir',
    'sbindir',
    'sharedstatedir',
    'sysconfdir',
    'auto_features',
    'backend',
    'buildtype',
    'debug',
    'default_library',
    'errorlogs',
    'genvslite',
    'install_umask',
    'layout',
    'optimization',
    'prefer_static',
    'stdsplit',
    'strip',
    'unity',
    'unity_size',
    'warning_level',
    'werror',
    'wrap_mode',
    'force_fallback_for',
    'pkg_config_path',
    'cmake_prefix_path',
    'vsenv',
}


    Since this class is not an actual iterator, but only provides a minimal
    fallback, it is safe to ignore the 'Iterator does not return self from
    __iter__ method' warning.
    '''
    def __init__(self, iterable: T.Optional[T.Iterable[str]] = None, total: T.Optional[int] = None,
                 bar_type: T.Optional[str] = None, desc: T.Optional[str] = None,
                 disable: T.Optional[bool] = None):
        if iterable is not None:
            self.iterable = iter(iterable)
            return
        self.total = total
        self.done = 0
        self.printed_dots = 0
        self.disable = not mlog.colorize_console() if disable is None else disable
        if not self.disable:
            if self.total and bar_type == 'download':
                print('Download size:', self.total)
            if desc:
                print(f'{desc}: ', end='')

    # Pretend to be an iterator when called as one and don't print any
    # progress
    def __iter__(self) -> T.Iterator[str]:
        return self.iterable

    def __next__(self) -> str:
        return next(self.iterable)

    def print_dot(self) -> None:
        if not self.disable:
            print('.', end='')
            sys.stdout.flush()
        self.printed_dots += 1

    def update(self, progress: int) -> None:
        self.done += progress
        if not self.total:
            # Just print one dot per call if we don't have a total length
            self.print_dot()
            return
        ratio = int(self.done / self.total * 10)
        while self.printed_dots < ratio:
            self.print_dot()

    def close(self) -> None:
        if not self.disable:
            print()

try:
    from tqdm import tqdm
except ImportError:
    # ideally we would use a typing.Protocol here, but it's part of typing_extensions until 3.8
    ProgressBar: T.Union[T.Type[ProgressBarFallback], T.Type[ProgressBarTqdm]] = ProgressBarFallback
else:
    class ProgressBarTqdm(tqdm):
        def __init__(self, *args: T.Any, bar_type: T.Optional[str] = None, **kwargs: T.Any) -> None:
            if bar_type == 'download':
                kwargs.update({'unit': 'B',
                               'unit_scale': True,
                               'unit_divisor': 1024,
                               'leave': True,
                               'bar_format': '{l_bar}{bar}| {n_fmt}/{total_fmt} {rate_fmt} eta {remaining}',
                               })

            else:
                kwargs.update({'leave': False,
                               'bar_format': '{l_bar}{bar}| {n_fmt}/{total_fmt} eta {remaining}',
                               })
            super().__init__(*args, **kwargs)

    ProgressBar = ProgressBarTqdm


class RealPathAction(argparse.Action):
    def __init__(self, option_strings: T.List[str], dest: str, default: str = '.', **kwargs: T.Any):
        default = os.path.abspath(os.path.realpath(default))
        super().__init__(option_strings, dest, nargs=None, default=default, **kwargs)

    def __call__(self, parser: argparse.ArgumentParser, namespace: argparse.Namespace,
                 values: T.Union[str, T.Sequence[T.Any], None], option_string: T.Optional[str] = None) -> None:
        assert isinstance(values, str)
        setattr(namespace, self.dest, os.path.abspath(os.path.realpath(values)))


def get_wine_shortpath(winecmd: T.List[str], wine_paths: T.List[str],
                       workdir: T.Optional[str] = None) -> str:
    '''
    WINEPATH size is limited to 1024 bytes which can easily be exceeded when
    adding the path to every dll inside build directory. See
    https://bugs.winehq.org/show_bug.cgi?id=45810.

    To shorten it as much as possible we use path relative to `workdir`
    where possible and convert absolute paths to Windows shortpath (e.g.
    "/usr/x86_64-w64-mingw32/lib" to "Z:\\usr\\X86_~FWL\\lib").

    This limitation reportedly has been fixed with wine >= 6.4
    '''

    # Remove duplicates
    wine_paths = list(OrderedSet(wine_paths))

    # Check if it's already short enough
    wine_path = ';'.join(wine_paths)
    if len(wine_path) <= 1024:
        return wine_path

    # Check if we have wine >= 6.4
    from ..programs import ExternalProgram
    wine = ExternalProgram('wine', winecmd, silent=True)
    if version_compare(wine.get_version(), '>=6.4'):
        return wine_path

    # Check paths that can be reduced by making them relative to workdir.
    rel_paths: T.List[str] = []
    if workdir:
        abs_paths: T.List[str] = []
        for p in wine_paths:
            try:
                rel = Path(p).relative_to(workdir)
                rel_paths.append(str(rel))
            except ValueError:
                abs_paths.append(p)
        wine_paths = abs_paths

    if wine_paths:
        # BAT script that takes a list of paths in argv and prints semi-colon separated shortpaths
        with NamedTemporaryFile('w', suffix='.bat', encoding='utf-8', delete=False) as bat_file:
            bat_file.write('''
            @ECHO OFF
            for %%x in (%*) do (
                echo|set /p=;%~sx
            )
            ''')
        try:
            stdout = subprocess.check_output(winecmd + ['cmd', '/C', bat_file.name] + wine_paths,
                                             encoding='utf-8', stderr=subprocess.DEVNULL)
            stdout = stdout.strip(';')
            if stdout:
                wine_paths = stdout.split(';')
            else:
                mlog.warning('Could not shorten WINEPATH: empty stdout')
        except subprocess.CalledProcessError as e:
            mlog.warning(f'Could not shorten WINEPATH: {str(e)}')
        finally:
            os.unlink(bat_file.name)
    wine_path = ';'.join(rel_paths + wine_paths)
    if len(wine_path) > 1024:
        mlog.warning('WINEPATH exceeds 1024 characters which could cause issues')
    return wine_path


def run_once(func: T.Callable[..., _T]) -> T.Callable[..., _T]:
    ret: T.List[_T] = []

    @wraps(func)
    def wrapper(*args: T.Any, **kwargs: T.Any) -> _T:
        if ret:
            return ret[0]

        val = func(*args, **kwargs)
        ret.append(val)
        return val

    return wrapper


def generate_list(func: T.Callable[..., T.Generator[_T, None, None]]) -> T.Callable[..., T.List[_T]]:
    @wraps(func)
    def wrapper(*args: T.Any, **kwargs: T.Any) -> T.List[_T]:
        return list(func(*args, **kwargs))

    return wrapper


def pickle_load(filename: str, object_name: str, object_type: T.Type[_PL], suggest_reconfigure: bool = True) -> _PL:
    load_fail_msg = f'{object_name} file {filename!r} is corrupted.'
    extra_msg = ' Consider reconfiguring the directory with "meson setup --reconfigure".' if suggest_reconfigure else ''
    try:
        with open(filename, 'rb') as f:
            obj = pickle.load(f)
    except (pickle.UnpicklingError, EOFError):
        raise MesonException(load_fail_msg + extra_msg)
    except (TypeError, ModuleNotFoundError, AttributeError):
        raise MesonException(
            f"{object_name} file {filename!r} references functions or classes that don't "
            "exist. This probably means that it was generated with an old "
            "version of meson." + extra_msg)

    if not isinstance(obj, object_type):
        raise MesonException(load_fail_msg + extra_msg)

    # Because these Protocols are not available at runtime (and cannot be made
    # available at runtime until we drop support for Python < 3.8), we have to
    # do a bit of hackery so that mypy understands what's going on here
    version: str
    if hasattr(obj, 'version'):
        version = T.cast('_VerPickleLoadable', obj).version
    else:
        version = T.cast('_EnvPickleLoadable', obj).environment.coredata.version

    from ..coredata import version as coredata_version
    from ..coredata import major_versions_differ, MesonVersionMismatchException
    if major_versions_differ(version, coredata_version):
        raise MesonVersionMismatchException(version, coredata_version, extra_msg)
    return obj


def first(iter: T.Iterable[_T], predicate: T.Callable[[_T], bool]) -> T.Optional[_T]:
    """Find the first entry in an iterable where the given predicate is true

    :param iter: The iterable to search
    :param predicate: A finding function that takes an element from the iterable
        and returns True if found, otherwise False
    :return: The first found element, or None if it is not found
    """
    for i in iter:
        if predicate(i):
            return i
    return None


def get_rsp_threshold() -> int:
    '''Return a conservative estimate of the commandline size in bytes
    above which a response file should be used.  May be overridden for
    debugging by setting environment variable MESON_RSP_THRESHOLD.'''

    if is_windows():
        # Usually 32k, but some projects might use cmd.exe,
        # and that has a limit of 8k.
        limit = 8192
    else:
        # Unix-like OSes usually have very large command line limits, (On Linux,
        # for example, this is limited by the kernel's MAX_ARG_STRLEN). However,
        # some programs place much lower limits, notably Wine which enforces a
        # 32k limit like Windows. Therefore, we limit the command line to 32k.
        limit = 32768

    # Be conservative
    limit = limit // 2
    return int(os.environ.get('MESON_RSP_THRESHOLD', limit))


class lazy_property(T.Generic[_T]):
    """Descriptor that replaces the function it wraps with the value generated.

    This property will only be calculated the first time it's queried, and will
    be cached and the cached value used for subsequent calls.

    This works by shadowing itself with the calculated value, in the instance.
    Due to Python's MRO that means that the calculated value will be found
    before this property, speeding up subsequent lookups.
    """
    def __init__(self, func: T.Callable[[T.Any], _T]) -> None:
        self.__name: T.Optional[str] = None
        self.__func = func

    def __set_name__(self, owner: T.Any, name: str) -> None:
        if self.__name is None:
            self.__name = name
        else:
            assert name == self.__name

    def __get__(self, instance: object, cls: T.Type) -> _T:
        value = self.__func(instance)
        setattr(instance, self.__name, value)
        return value
