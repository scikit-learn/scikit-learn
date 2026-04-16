# SPDX-License-Identifier: Apache-2.0
# Copyright 2012-2016 The Meson development team

from __future__ import annotations

from dataclasses import dataclass
import typing as T
from enum import Enum
import os
import platform
import sys

from . import mesonlib
from .mesonlib import EnvironmentException, HoldableObject, lazy_property, Popen_safe
from .programs import ExternalProgram
from . import mlog
from pathlib import Path, PurePath, PurePosixPath, PureWindowsPath

if T.TYPE_CHECKING:
    from .options import ElementaryOptionValues
    from .compilers.compilers import CompilerDict
    from .compilers.mixins.visualstudio import VisualStudioLikeCompiler
    from ._typing import ImmutableListProtocol


# These classes contains all the data pulled from configuration files (native
# and cross file currently), and also assists with the reading environment
# variables.
#
# At this time there isn't an ironclad difference between this and other sources
# of state like `coredata`. But one rough guide is much what is in `coredata` is
# the *output* of the configuration process: the final decisions after tests.
# This, on the other hand has *inputs*. The config files are parsed, but
# otherwise minimally transformed. When more complex fallbacks (environment
# detection) exist, they are defined elsewhere as functions that construct
# instances of these classes.


known_cpu_families = (
    'aarch64',
    'alpha',
    'arc',
    'arm',
    'avr',
    'c2000',
    'c6000',
    'csky',
    'dspic',
    'e2k',
    'ft32',
    'ia64',
    'loongarch64',
    'm68k',
    'microblaze',
    'mips',
    'mips64',
    'msp430',
    'parisc',
    'pic24',
    'pic32',
    'ppc',
    'ppc64',
    'riscv32',
    'riscv64',
    'rl78',
    'rx',
    's390',
    's390x',
    'sh4',
    'sparc',
    'sparc64',
    'sw_64',
    'wasm32',
    'wasm64',
    'x86',
    'x86_64',
    'tricore'
)

# It would feel more natural to call this "64_BIT_CPU_FAMILIES", but
# python identifiers cannot start with numbers
CPU_FAMILIES_64_BIT = [
    'aarch64',
    'alpha',
    'ia64',
    'loongarch64',
    'mips64',
    'ppc64',
    'riscv64',
    's390x',
    'sparc64',
    'sw_64',
    'wasm64',
    'x86_64',
]

# Map from language identifiers to environment variables.
ENV_VAR_COMPILER_MAP: T.Mapping[str, ImmutableListProtocol[str]] = {
    # Compilers
    'c': ['CC'],
    'cpp': ['CXX'],
    'cs': ['CSC'],
    'cython': ['CYTHON'],
    'd': ['DC'],
    'fortran': ['FC'],
    'objc': ['OBJC'],
    'objcpp': ['OBJCXX'],
    'rust': ['RUSTC'],
    'vala': ['VALAC'],
    'nasm': ['NASM'],

    # Linkers
    'c_ld': ['CC_LD'],
    'cpp_ld': ['CXX_LD'],
    'd_ld': ['DC_LD'],
    'fortran_ld': ['FC_LD'],
    'objc_ld': ['OBJC_LD'],
    'objcpp_ld': ['OBJCXX_LD'],
    'rust_ld': ['RUSTC_LD'],
}

# Map from utility names to environment variables.
ENV_VAR_TOOL_MAP: T.Mapping[str, ImmutableListProtocol[str]] = {
    # Binutils
    'ar': ['AR'],
    'as': ['AS'],
    'ld': ['LD'],
    'nm': ['NM'],
    'objcopy': ['OBJCOPY'],
    'objdump': ['OBJDUMP'],
    'ranlib': ['RANLIB'],
    'readelf': ['READELF'],
    'size': ['SIZE'],
    'strings': ['STRINGS'],
    'strip': ['STRIP'],
    'windres': ['RC', 'WINDRES'],

    # Other tools
    'cmake': ['CMAKE'],
    'qmake': ['QMAKE'],
    'pkg-config': ['PKG_CONFIG'],
    'make': ['MAKE'],
    'vapigen': ['VAPIGEN'],
    'llvm-config': ['LLVM_CONFIG'],
}

ENV_VAR_PROG_MAP = {**ENV_VAR_COMPILER_MAP, **ENV_VAR_TOOL_MAP}

# Deprecated environment variables mapped from the new variable to the old one
# Deprecated in 0.54.0
DEPRECATED_ENV_PROG_MAP: T.Mapping[str, ImmutableListProtocol[str]] = {
    'd_ld': ['D_LD'],
    'fortran_ld': ['F_LD'],
    'rust_ld': ['RUST_LD'],
    'objcpp_ld': ['OBJCPP_LD'],
}

class CMakeSkipCompilerTest(Enum):
    ALWAYS = 'always'
    NEVER = 'never'
    DEP_ONLY = 'dep_only'

class Properties:
    def __init__(
            self,
            properties: T.Optional[T.Dict[str, ElementaryOptionValues]] = None,
    ):
        self.properties = properties or {}

    def has_stdlib(self, language: str) -> bool:
        return language + '_stdlib' in self.properties

    # Some of get_stdlib, get_root, get_sys_root are wider than is actually
    # true, but without heterogeneous dict annotations it's not practical to
    # narrow them
    def get_stdlib(self, language: str) -> T.Union[str, T.List[str]]:
        stdlib = self.properties[language + '_stdlib']
        if isinstance(stdlib, str):
            return stdlib
        assert isinstance(stdlib, list)
        for i in stdlib:
            assert isinstance(i, str)
        return stdlib

    def get_root(self) -> T.Optional[str]:
        root = self.properties.get('root', None)
        assert root is None or isinstance(root, str)
        return root

    def get_sys_root(self) -> T.Optional[str]:
        sys_root = self.properties.get('sys_root', None)
        assert sys_root is None or isinstance(sys_root, str)
        return sys_root

    def get_pkg_config_libdir(self) -> T.Optional[T.List[str]]:
        p = self.properties.get('pkg_config_libdir', None)
        if p is None:
            return p
        res = mesonlib.listify(p)
        for i in res:
            assert isinstance(i, str)
        return res

    def get_cmake_defaults(self) -> bool:
        if 'cmake_defaults' not in self.properties:
            return True
        res = self.properties['cmake_defaults']
        assert isinstance(res, bool)
        return res

    def get_cmake_toolchain_file(self) -> T.Optional[Path]:
        if 'cmake_toolchain_file' not in self.properties:
            return None
        raw = self.properties['cmake_toolchain_file']
        assert isinstance(raw, str)
        cmake_toolchain_file = Path(raw)
        if not cmake_toolchain_file.is_absolute():
            raise EnvironmentException(f'cmake_toolchain_file ({raw}) is not absolute')
        return cmake_toolchain_file

    def get_cmake_skip_compiler_test(self) -> CMakeSkipCompilerTest:
        if 'cmake_skip_compiler_test' not in self.properties:
            return CMakeSkipCompilerTest.DEP_ONLY
        raw = self.properties['cmake_skip_compiler_test']
        assert isinstance(raw, str)
        try:
            return CMakeSkipCompilerTest(raw)
        except ValueError:
            raise EnvironmentException(
                '"{}" is not a valid value for cmake_skip_compiler_test. Supported values are {}'
                .format(raw, [e.value for e in CMakeSkipCompilerTest]))

    def get_cmake_use_exe_wrapper(self) -> bool:
        if 'cmake_use_exe_wrapper' not in self.properties:
            return True
        res = self.properties['cmake_use_exe_wrapper']
        assert isinstance(res, bool)
        return res

    def get_java_home(self) -> T.Optional[Path]:
        value = T.cast('T.Optional[str]', self.properties.get('java_home'))
        return Path(value) if value else None

    def get_bindgen_clang_args(self) -> T.List[str]:
        value = mesonlib.listify(self.properties.get('bindgen_clang_arguments', []))
        if not all(isinstance(v, str) for v in value):
            raise EnvironmentException('bindgen_clang_arguments must be a string or an array of strings')
        return T.cast('T.List[str]', value)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, type(self)):
            return self.properties == other.properties
        return NotImplemented

    # TODO consider removing so Properties is less freeform
    def __getitem__(self, key: str) -> T.Optional[T.Union[str, bool, int, T.List[str]]]:
        return self.properties[key]

    # TODO consider removing so Properties is less freeform
    def __contains__(self, item: T.Union[str, bool, int, T.List[str]]) -> bool:
        return item in self.properties

    # TODO consider removing, for same reasons as above
    def get(self, key: str, default: T.Optional[T.Union[str, bool, int, T.List[str]]] = None) -> T.Optional[T.Union[str, bool, int, T.List[str]]]:
        return self.properties.get(key, default)

@dataclass(unsafe_hash=True)
class MachineInfo(HoldableObject):
    system: str
    cpu_family: str
    cpu: str
    endian: str
    kernel: T.Optional[str]
    subsystem: T.Optional[str]

    def __post_init__(self) -> None:
        self.is_64_bit: bool = self.cpu_family in CPU_FAMILIES_64_BIT

    def __repr__(self) -> str:
        return f'<MachineInfo: {self.system} {self.cpu_family} ({self.cpu})>'

    @classmethod
    def from_literal(cls, raw: T.Dict[str, ElementaryOptionValues]) -> 'MachineInfo':
        # We don't have enough type information to be sure of what we loaded
        # So we need to accept that this might have ElementaryOptionValues, but
        # then ensure that it's actually strings, since that's what the
        # [*_machine] section should have.
        assert all(isinstance(v, str) for v in raw.values()), 'for mypy'
        literal = T.cast('T.Dict[str, str]', raw)
        minimum_literal = {'cpu', 'cpu_family', 'endian', 'system'}
        if minimum_literal - set(literal):
            raise EnvironmentException(
                f'Machine info is currently {literal}\n' +
                'but is missing {}.'.format(minimum_literal - set(literal)))

        cpu_family = literal['cpu_family']
        if cpu_family not in known_cpu_families:
            mlog.warning(f'Unknown CPU family {cpu_family}, please report this at https://github.com/mesonbuild/meson/issues/new')

        endian = literal['endian']
        if endian not in ('little', 'big'):
            mlog.warning(f'Unknown endian {endian}')

        system = literal['system']
        kernel = literal.get('kernel', None)
        subsystem = literal.get('subsystem', None)

        return cls(system, cpu_family, literal['cpu'], endian, kernel, subsystem)

    def is_windows(self) -> bool:
        """
        Machine is windows?
        """
        return self.system == 'windows'

    def is_cygwin(self) -> bool:
        """
        Machine is cygwin?
        """
        return self.system == 'cygwin'

    @lazy_property
    def pure_path_class(self) -> T.Type[PurePath]:
        """Get the appropriate PurePath class for this machine."""
        if self.is_windows():
            return PureWindowsPath
        return PurePosixPath

    def is_linux(self) -> bool:
        """
        Machine is linux?
        """
        return self.system == 'linux'

    def is_darwin(self) -> bool:
        """
        Machine is Darwin (macOS/iOS/tvOS/visionOS/watchOS)?
        """
        return self.system in {'darwin', 'ios', 'tvos', 'visionos', 'watchos'}

    def is_android(self) -> bool:
        """
        Machine is Android?
        """
        return self.system == 'android'

    def is_haiku(self) -> bool:
        """
        Machine is Haiku?
        """
        return self.system == 'haiku'

    def is_netbsd(self) -> bool:
        """
        Machine is NetBSD?
        """
        return self.system == 'netbsd'

    def is_openbsd(self) -> bool:
        """
        Machine is OpenBSD?
        """
        return self.system == 'openbsd'

    def is_dragonflybsd(self) -> bool:
        """Machine is DragonflyBSD?"""
        return self.system == 'dragonfly'

    def is_freebsd(self) -> bool:
        """Machine is FreeBSD?"""
        return self.system == 'freebsd'

    def is_sunos(self) -> bool:
        """Machine is illumos or Solaris?"""
        return self.system == 'sunos'

    def is_hurd(self) -> bool:
        """
        Machine is GNU/Hurd?
        """
        return self.system == 'gnu'

    def is_aix(self) -> bool:
        """
        Machine is aix?
        """
        return self.system == 'aix'

    def is_irix(self) -> bool:
        """Machine is IRIX?"""
        return self.system.startswith('irix')

    def is_os2(self) -> bool:
        """
        Machine is OS/2?
        """
        return self.system == 'os/2'

    # Various prefixes and suffixes for import libraries, shared libraries,
    # static libraries, and executables.
    # Versioning is added to these names in the backends as-needed.
    def get_exe_suffix(self) -> str:
        if self.is_windows() or self.is_cygwin() or self.is_os2():
            return 'exe'
        else:
            return ''

    def get_object_suffix(self) -> str:
        if self.is_windows():
            return 'obj'
        else:
            return 'o'

    def libdir_layout_is_win(self) -> bool:
        return self.is_windows() or self.is_cygwin()

class BinaryTable:

    def __init__(
            self,
            binaries: T.Optional[T.Mapping[str, ElementaryOptionValues]] = None,
    ):
        self.binaries: T.Dict[str, T.List[str]] = {}
        if binaries:
            for name, command in binaries.items():
                if not isinstance(command, (list, str)):
                    raise mesonlib.MesonException(
                        f'Invalid type {command!r} for entry {name!r} in cross file')
                self.binaries[name] = mesonlib.listify(command)
            if 'pkgconfig' in self.binaries:
                if 'pkg-config' not in self.binaries:
                    mlog.deprecation('"pkgconfig" entry is deprecated and should be replaced by "pkg-config"', fatal=False)
                    self.binaries['pkg-config'] = self.binaries['pkgconfig']
                elif self.binaries['pkgconfig'] != self.binaries['pkg-config']:
                    raise mesonlib.MesonException('Mismatched pkgconfig and pkg-config binaries in the machine file.')
                else:
                    # Both are defined with the same value, this is allowed
                    # for backward compatibility.
                    # FIXME: We should still print deprecation warning if the
                    # project targets Meson >= 1.3.0, but we have no way to know
                    # that here.
                    pass
                del self.binaries['pkgconfig']

    @staticmethod
    def detect_ccache() -> ExternalProgram:
        return ExternalProgram('ccache', silent=True)

    @staticmethod
    def detect_sccache() -> ExternalProgram:
        return ExternalProgram('sccache', silent=True)

    @staticmethod
    def detect_compiler_cache() -> ExternalProgram:
        # Sccache is "newer" so it is assumed that people would prefer it by default.
        cache = BinaryTable.detect_sccache()
        if cache.found():
            return cache
        return BinaryTable.detect_ccache()

    @classmethod
    def parse_entry(cls, entry: T.Union[str, T.List[str]]) -> T.Tuple[T.List[str], T.Union[None, ExternalProgram]]:
        parts = mesonlib.stringlistify(entry)
        # Ensure ccache exists and remove it if it doesn't
        if parts[0] == 'ccache':
            compiler = parts[1:]
            ccache = cls.detect_ccache()
        elif parts[0] == 'sccache':
            compiler = parts[1:]
            ccache = cls.detect_sccache()
        else:
            compiler = parts
            ccache = None
        if not compiler:
            raise EnvironmentException(f'Compiler cache specified without compiler: {parts[0]}')
        # Return value has to be a list of compiler 'choices'
        return compiler, ccache

    def lookup_entry(self, name: str) -> T.Optional[T.List[str]]:
        """Lookup binary in cross/native file and fallback to environment.

        Returns command with args as list if found, Returns `None` if nothing is
        found.
        """
        command = self.binaries.get(name)
        if not command:
            return None
        elif not command[0].strip():
            return None
        return command

class CMakeVariables:
    def __init__(self, variables: T.Optional[T.Dict[str, T.Any]] = None) -> None:
        variables = variables or {}
        self.variables: T.Dict[str, T.List[str]] = {}

        for key, value in variables.items():
            value = mesonlib.listify(value)
            for i in value:
                if not isinstance(i, str):
                    raise EnvironmentException(f"Value '{i}' of CMake variable '{key}' defined in a machine file is a {type(i).__name__} and not a str")
            self.variables[key] = value

    def get_variables(self) -> T.Dict[str, T.List[str]]:
        return self.variables


# Machine and platform detection functions
# ========================================

KERNEL_MAPPINGS: T.Mapping[str, str] = {'freebsd': 'freebsd',
                                        'openbsd': 'openbsd',
                                        'netbsd': 'netbsd',
                                        'windows': 'nt',
                                        'android': 'linux',
                                        'linux': 'linux',
                                        'cygwin': 'nt',
                                        'darwin': 'xnu',
                                        'ios': 'xnu',
                                        'tvos': 'xnu',
                                        'visionos': 'xnu',
                                        'watchos': 'xnu',
                                        'dragonfly': 'dragonfly',
                                        'haiku': 'haiku',
                                        'gnu': 'gnu',
                                        }

def detect_windows_arch(compilers: CompilerDict) -> str:
    """
    Detecting the 'native' architecture of Windows is not a trivial task. We
    cannot trust that the architecture that Python is built for is the 'native'
    one because you can run 32-bit apps on 64-bit Windows using WOW64 and
    people sometimes install 32-bit Python on 64-bit Windows.

    We also can't rely on the architecture of the OS itself, since it's
    perfectly normal to compile and run 32-bit applications on Windows as if
    they were native applications. It's a terrible experience to require the
    user to supply a cross-info file to compile 32-bit applications on 64-bit
    Windows. Thankfully, the only way to compile things with Visual Studio on
    Windows is by entering the 'msvc toolchain' environment, which can be
    easily detected.

    In the end, the sanest method is as follows:
    1. Check environment variables that are set by Windows and WOW64 to find out
       if this is x86 (possibly in WOW64), if so use that as our 'native'
       architecture.
    2. If the compiler toolchain target architecture is x86, use that as our
      'native' architecture.
    3. Otherwise, use the actual Windows architecture

    """
    os_arch = mesonlib.windows_detect_native_arch()
    if os_arch == 'x86':
        return os_arch
    # If we're on 64-bit Windows, 32-bit apps can be compiled without
    # cross-compilation. So if we're doing that, just set the native arch as
    # 32-bit and pretend like we're running under WOW64. Else, return the
    # actual Windows architecture that we deduced above.
    for compiler in compilers.values():
        compiler = T.cast('VisualStudioLikeCompiler', compiler)
        if compiler.id == 'msvc' and (compiler.target in {'x86', '80x86'}):
            return 'x86'
        if compiler.id == 'clang-cl' and (compiler.target in {'x86', 'i686'}):
            return 'x86'
        if compiler.id == 'gcc' and compiler.has_builtin_define('__i386__'):
            return 'x86'
    return os_arch

def any_compiler_has_define(compilers: CompilerDict, define: str) -> bool:
    for c in compilers.values():
        try:
            if c.has_builtin_define(define):
                return True
        except mesonlib.MesonException:
            # Ignore compilers that do not support has_builtin_define.
            pass
    return False

def detect_cpu_family(compilers: CompilerDict) -> str:
    """
    Python is inconsistent in its platform module.
    It returns different values for the same cpu.
    For x86 it might return 'x86', 'i686' or some such.
    Do some canonicalization.
    """
    if mesonlib.is_windows():
        trial = detect_windows_arch(compilers)
    elif mesonlib.is_freebsd() or mesonlib.is_netbsd() or mesonlib.is_openbsd() or mesonlib.is_qnx() or mesonlib.is_aix():
        trial = platform.processor().lower()
    else:
        trial = platform.machine().lower()
    if trial.startswith('i') and trial.endswith('86'):
        trial = 'x86'
    elif trial == 'bepc':
        trial = 'x86'
    elif trial == 'arm64':
        trial = 'aarch64'
    elif trial.startswith('aarch64'):
        # This can be `aarch64_be`
        trial = 'aarch64'
    elif trial.startswith('arm') or trial.startswith('earm'):
        trial = 'arm'
    elif trial.startswith(('powerpc64', 'ppc64')):
        trial = 'ppc64'
    elif trial.startswith(('powerpc', 'ppc')) or trial in {'macppc', 'power macintosh'}:
        trial = 'ppc'
    elif trial in {'amd64', 'x64', 'i86pc'}:
        trial = 'x86_64'
    elif trial in {'sun4u', 'sun4v'}:
        trial = 'sparc64'
    elif trial.startswith('mips'):
        if '64' not in trial:
            trial = 'mips'
        else:
            trial = 'mips64'
    elif trial in {'ip30', 'ip35'}:
        trial = 'mips64'

    # On Linux (and maybe others) there can be any mixture of 32/64 bit code in
    # the kernel, Python, system, 32-bit chroot on 64-bit host, etc. The only
    # reliable way to know is to check the compiler defines.
    if trial == 'x86_64':
        if any_compiler_has_define(compilers, '__i386__'):
            trial = 'x86'
    elif trial == 'aarch64':
        if any_compiler_has_define(compilers, '__arm__'):
            trial = 'arm'
    # Add more quirks here as bugs are reported. Keep in sync with detect_cpu()
    # below.
    elif trial == 'parisc64':
        # ATM there is no 64 bit userland for PA-RISC. Thus always
        # report it as 32 bit for simplicity.
        trial = 'parisc'
    elif trial == 'ppc':
        # AIX always returns powerpc, check here for 64-bit
        if any_compiler_has_define(compilers, '__64BIT__'):
            trial = 'ppc64'
    # MIPS64 is able to run MIPS32 code natively, so there is a chance that
    # such mixture mentioned above exists.
    elif trial == 'mips64':
        if compilers and not any_compiler_has_define(compilers, '__mips64'):
            trial = 'mips'

    if trial not in known_cpu_families:
        mlog.warning(f'Unknown CPU family {trial!r}, please report this at '
                     'https://github.com/mesonbuild/meson/issues/new with the '
                     'output of `uname -a` and `cat /proc/cpuinfo`')

    return trial

def detect_cpu(compilers: CompilerDict) -> str:
    if mesonlib.is_windows():
        trial = detect_windows_arch(compilers)
    elif mesonlib.is_freebsd() or mesonlib.is_netbsd() or mesonlib.is_openbsd() or mesonlib.is_aix():
        trial = platform.processor().lower()
    else:
        trial = platform.machine().lower()

    if trial in {'amd64', 'x64', 'i86pc'}:
        trial = 'x86_64'
    if trial == 'x86_64':
        # Same check as above for cpu_family
        if any_compiler_has_define(compilers, '__i386__'):
            trial = 'i686' # All 64 bit cpus have at least this level of x86 support.
    elif trial.startswith('aarch64') or trial.startswith('arm64'):
        # Same check as above for cpu_family
        if any_compiler_has_define(compilers, '__arm__'):
            trial = 'arm'
        else:
            # for aarch64_be
            trial = 'aarch64'
    elif trial.startswith('earm'):
        trial = 'arm'
    elif trial == 'e2k':
        # Make more precise CPU detection for Elbrus platform.
        trial = platform.processor().lower()
    elif trial.startswith('mips'):
        if '64' not in trial:
            trial = 'mips'
        else:
            if compilers and not any_compiler_has_define(compilers, '__mips64'):
                trial = 'mips'
            else:
                trial = 'mips64'
    elif trial == 'ppc':
        # AIX always returns powerpc, check here for 64-bit
        if any_compiler_has_define(compilers, '__64BIT__'):
            trial = 'ppc64'

    # Add more quirks here as bugs are reported. Keep in sync with
    # detect_cpu_family() above.
    return trial

def detect_kernel(system: str) -> T.Optional[str]:
    if system == 'sunos':
        # Solaris 5.10 uname doesn't support the -o switch, and illumos started
        # with version 5.11 so shortcut the logic to report 'solaris' in such
        # cases where the version is 5.10 or below.
        if mesonlib.version_compare(platform.uname().release, '<=5.10'):
            return 'solaris'
        # This needs to be /usr/bin/uname because gnu-uname could be installed and
        # won't provide the necessary information
        p, out, _ = Popen_safe(['/usr/bin/uname', '-o'])
        if p.returncode != 0:
            raise mesonlib.MesonException('Failed to run "/usr/bin/uname -o"')
        out = out.lower().strip()
        if out not in {'illumos', 'solaris'}:
            mlog.warning(f'Got an unexpected value for kernel on a SunOS derived platform, expected either "illumos" or "solaris", but got "{out}".'
                         "Please open a Meson issue with the OS you're running and the value detected for your kernel.")
            return None
        return out
    return KERNEL_MAPPINGS.get(system, None)

def detect_subsystem(system: str) -> T.Optional[str]:
    if system == 'darwin':
        return 'macos'
    return system

def detect_system() -> str:
    if sys.platform == 'cygwin':
        return 'cygwin'
    return platform.system().lower()

def detect_msys2_arch() -> T.Optional[str]:
    return os.environ.get('MSYSTEM_CARCH', None)

def detect_machine_info(compilers: T.Optional[CompilerDict] = None) -> MachineInfo:
    """Detect the machine we're running on

    If compilers are not provided, we cannot know as much. None out those
    fields to avoid accidentally depending on partial knowledge. The
    underlying ''detect_*'' method can be called to explicitly use the
    partial information.
    """
    system = detect_system()
    return MachineInfo(
        system,
        detect_cpu_family(compilers) if compilers is not None else None,
        detect_cpu(compilers) if compilers is not None else None,
        sys.byteorder,
        detect_kernel(system),
        detect_subsystem(system))

# TODO make this compare two `MachineInfo`s purely. How important is the
# `detect_cpu_family({})` distinction? It is the one impediment to that.
def machine_info_can_run(machine_info: MachineInfo) -> bool:
    """Whether we can run binaries for this machine on the current machine.

    Can almost always run 32-bit binaries on 64-bit natively if the host
    and build systems are the same. We don't pass any compilers to
    detect_cpu_family() here because we always want to know the OS
    architecture, not what the compiler environment tells us.
    """
    system = detect_system()
    if machine_info.system != system:
        return False
    if machine_info.subsystem != detect_subsystem(system):
        return False
    true_build_cpu_family = detect_cpu_family({})
    assert machine_info.cpu_family is not None, 'called on incomplete machine_info'
    return \
        (machine_info.cpu_family == true_build_cpu_family) or \
        ((true_build_cpu_family == 'x86_64') and (machine_info.cpu_family == 'x86')) or \
        ((true_build_cpu_family == 'mips64') and (machine_info.cpu_family == 'mips'))
