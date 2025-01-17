# SPDX-License-Identifier: Apache-2.0
# Copyright 2012-2022 The Meson development team
# Copyright © 2023-2024 Intel Corporation

from __future__ import annotations

import abc
import contextlib, os.path, re
import enum
import itertools
import typing as T
from dataclasses import dataclass, field
from functools import lru_cache

from .. import mlog
from .. import mesonlib
from .. import options
from ..mesonlib import (
    HoldableObject,
    EnvironmentException, MesonException,
    Popen_safe_logged, LibType, TemporaryDirectoryWinProof,
)

from ..options import OptionKey

from ..arglist import CompilerArgs

if T.TYPE_CHECKING:
    from .. import coredata
    from ..build import BuildTarget, DFeatures
    from ..coredata import MutableKeyedOptionDictType, KeyedOptionDictType
    from ..envconfig import MachineInfo
    from ..environment import Environment
    from ..linkers import RSPFileSyntax
    from ..linkers.linkers import DynamicLinker
    from ..mesonlib import MachineChoice
    from ..dependencies import Dependency

    CompilerType = T.TypeVar('CompilerType', bound='Compiler')
    _T = T.TypeVar('_T')
    UserOptionType = T.TypeVar('UserOptionType', bound=options.UserOption)

"""This file contains the data files of all compilers Meson knows
about. To support a new compiler, add its information below.
Also add corresponding autodetection code in detect.py."""

header_suffixes = {'h', 'hh', 'hpp', 'hxx', 'H', 'ipp', 'moc', 'vapi', 'di'}
obj_suffixes = {'o', 'obj', 'res'}
# To the emscripten compiler, .js files are libraries
lib_suffixes = {'a', 'lib', 'dll', 'dll.a', 'dylib', 'so', 'js'}
# Mapping of language to suffixes of files that should always be in that language
# This means we can't include .h headers here since they could be C, C++, ObjC, etc.
# First suffix is the language's default.
lang_suffixes: T.Mapping[str, T.Tuple[str, ...]] = {
    'c': ('c',),
    'cpp': ('cpp', 'cppm', 'cc', 'cp', 'cxx', 'c++', 'hh', 'hp', 'hpp', 'ipp', 'hxx', 'h++', 'ino', 'ixx', 'CPP', 'C', 'HPP', 'H'),
    'cuda': ('cu',),
    # f90, f95, f03, f08 are for free-form fortran ('f90' recommended)
    # f, for, ftn, fpp are for fixed-form fortran ('f' or 'for' recommended)
    'fortran': ('f90', 'f95', 'f03', 'f08', 'f', 'for', 'ftn', 'fpp'),
    'd': ('d', 'di'),
    'objc': ('m',),
    'objcpp': ('mm',),
    'rust': ('rs',),
    'vala': ('vala', 'vapi', 'gs'),
    'cs': ('cs',),
    'swift': ('swift',),
    'java': ('java',),
    'cython': ('pyx', ),
    'nasm': ('asm', 'nasm',),
    'masm': ('masm',),
}
all_languages = lang_suffixes.keys()
c_cpp_suffixes = {'h'}
cpp_suffixes = set(lang_suffixes['cpp']) | c_cpp_suffixes
c_suffixes = set(lang_suffixes['c']) | c_cpp_suffixes
assembler_suffixes = {'s', 'S', 'sx', 'asm', 'masm'}
llvm_ir_suffixes = {'ll'}
all_suffixes = set(itertools.chain(*lang_suffixes.values(), assembler_suffixes, llvm_ir_suffixes, c_cpp_suffixes))
source_suffixes = all_suffixes - header_suffixes
# List of languages that by default consume and output libraries following the
# C ABI; these can generally be used interchangeably
# This must be sorted, see sort_clink().
clib_langs = ('objcpp', 'cpp', 'objc', 'c', 'nasm', 'fortran')
# List of languages that can be linked with C code directly by the linker
# used in build.py:process_compilers() and build.py:get_dynamic_linker()
# This must be sorted, see sort_clink().
clink_langs = ('d', 'cuda') + clib_langs

SUFFIX_TO_LANG = dict(itertools.chain(*(
    [(suffix, lang) for suffix in v] for lang, v in lang_suffixes.items())))

# Languages that should use LDFLAGS arguments when linking.
LANGUAGES_USING_LDFLAGS = {'objcpp', 'cpp', 'objc', 'c', 'fortran', 'd', 'cuda'}
# Languages that should use CPPFLAGS arguments when linking.
LANGUAGES_USING_CPPFLAGS = {'c', 'cpp', 'objc', 'objcpp'}
soregex = re.compile(r'.*\.so(\.[0-9]+)?(\.[0-9]+)?(\.[0-9]+)?$')

# Environment variables that each lang uses.
CFLAGS_MAPPING: T.Mapping[str, str] = {
    'c': 'CFLAGS',
    'cpp': 'CXXFLAGS',
    'cuda': 'CUFLAGS',
    'objc': 'OBJCFLAGS',
    'objcpp': 'OBJCXXFLAGS',
    'fortran': 'FFLAGS',
    'd': 'DFLAGS',
    'vala': 'VALAFLAGS',
    'rust': 'RUSTFLAGS',
    'cython': 'CYTHONFLAGS',
    'cs': 'CSFLAGS', # This one might not be standard.
}

# All these are only for C-linkable languages; see `clink_langs` above.

def sort_clink(lang: str) -> int:
    '''
    Sorting function to sort the list of languages according to
    reversed(compilers.clink_langs) and append the unknown langs in the end.
    The purpose is to prefer C over C++ for files that can be compiled by
    both such as assembly, C, etc. Also applies to ObjC, ObjC++, etc.
    '''
    if lang not in clink_langs:
        return 1
    return -clink_langs.index(lang)

def is_header(fname: 'mesonlib.FileOrString') -> bool:
    if isinstance(fname, mesonlib.File):
        fname = fname.fname
    suffix = fname.split('.')[-1]
    return suffix in header_suffixes

def is_source_suffix(suffix: str) -> bool:
    return suffix in source_suffixes

def is_source(fname: 'mesonlib.FileOrString') -> bool:
    if isinstance(fname, mesonlib.File):
        fname = fname.fname
    suffix = fname.split('.')[-1].lower()
    return is_source_suffix(suffix)

def is_assembly(fname: 'mesonlib.FileOrString') -> bool:
    if isinstance(fname, mesonlib.File):
        fname = fname.fname
    suffix = fname.split('.')[-1]
    return suffix in assembler_suffixes

def is_llvm_ir(fname: 'mesonlib.FileOrString') -> bool:
    if isinstance(fname, mesonlib.File):
        fname = fname.fname
    suffix = fname.split('.')[-1]
    return suffix in llvm_ir_suffixes

@lru_cache(maxsize=None)
def cached_by_name(fname: 'mesonlib.FileOrString') -> bool:
    suffix = fname.split('.')[-1]
    return suffix in obj_suffixes

def is_object(fname: 'mesonlib.FileOrString') -> bool:
    if isinstance(fname, mesonlib.File):
        fname = fname.fname
    return cached_by_name(fname)

def is_library(fname: 'mesonlib.FileOrString') -> bool:
    if isinstance(fname, mesonlib.File):
        fname = fname.fname

    if soregex.match(fname):
        return True

    suffix = fname.split('.')[-1]
    return suffix in lib_suffixes

def is_known_suffix(fname: 'mesonlib.FileOrString') -> bool:
    if isinstance(fname, mesonlib.File):
        fname = fname.fname
    suffix = fname.split('.')[-1]

    return suffix in all_suffixes


class CompileCheckMode(enum.Enum):

    PREPROCESS = 'preprocess'
    COMPILE = 'compile'
    LINK = 'link'


gnu_winlibs = ['-lkernel32', '-luser32', '-lgdi32', '-lwinspool', '-lshell32',
               '-lole32', '-loleaut32', '-luuid', '-lcomdlg32', '-ladvapi32']

msvc_winlibs = ['kernel32.lib', 'user32.lib', 'gdi32.lib',
                'winspool.lib', 'shell32.lib', 'ole32.lib', 'oleaut32.lib',
                'uuid.lib', 'comdlg32.lib', 'advapi32.lib']

clike_optimization_args: T.Dict[str, T.List[str]] = {
    'plain': [],
    '0': [],
    'g': [],
    '1': ['-O1'],
    '2': ['-O2'],
    '3': ['-O3'],
    's': ['-Os'],
}

clike_debug_args: T.Dict[bool, T.List[str]] = {
    False: [],
    True: ['-g']
}


MSCRT_VALS = ['none', 'md', 'mdd', 'mt', 'mtd']

@dataclass
class BaseOption(T.Generic[options._T, options._U]):
    opt_type: T.Type[options._U]
    description: str
    default: T.Any = None
    choices: T.Any = None

    def init_option(self, name: OptionKey) -> options._U:
        keywords = {'value': self.default}
        if self.choices:
            keywords['choices'] = self.choices
        return self.opt_type(name.name, self.description, **keywords)

BASE_OPTIONS: T.Mapping[OptionKey, BaseOption] = {
    OptionKey('b_pch'): BaseOption(options.UserBooleanOption, 'Use precompiled headers', True),
    OptionKey('b_lto'): BaseOption(options.UserBooleanOption, 'Use link time optimization', False),
    OptionKey('b_lto_threads'): BaseOption(options.UserIntegerOption, 'Use multiple threads for Link Time Optimization', (None, None, 0)),
    OptionKey('b_lto_mode'): BaseOption(options.UserComboOption, 'Select between different LTO modes.', 'default',
                                        choices=['default', 'thin']),
    OptionKey('b_thinlto_cache'): BaseOption(options.UserBooleanOption, 'Use LLVM ThinLTO caching for faster incremental builds', False),
    OptionKey('b_thinlto_cache_dir'): BaseOption(options.UserStringOption, 'Directory to store ThinLTO cache objects', ''),
    OptionKey('b_sanitize'): BaseOption(options.UserComboOption, 'Code sanitizer to use', 'none',
                                        choices=['none', 'address', 'thread', 'undefined', 'memory', 'leak', 'address,undefined']),
    OptionKey('b_lundef'): BaseOption(options.UserBooleanOption, 'Use -Wl,--no-undefined when linking', True),
    OptionKey('b_asneeded'): BaseOption(options.UserBooleanOption, 'Use -Wl,--as-needed when linking', True),
    OptionKey('b_pgo'): BaseOption(options.UserComboOption, 'Use profile guided optimization', 'off',
                                   choices=['off', 'generate', 'use']),
    OptionKey('b_coverage'): BaseOption(options.UserBooleanOption, 'Enable coverage tracking.', False),
    OptionKey('b_colorout'): BaseOption(options.UserComboOption, 'Use colored output', 'always',
                                        choices=['auto', 'always', 'never']),
    OptionKey('b_ndebug'): BaseOption(options.UserComboOption, 'Disable asserts', 'false', choices=['true', 'false', 'if-release']),
    OptionKey('b_staticpic'): BaseOption(options.UserBooleanOption, 'Build static libraries as position independent', True),
    OptionKey('b_pie'): BaseOption(options.UserBooleanOption, 'Build executables as position independent', False),
    OptionKey('b_bitcode'): BaseOption(options.UserBooleanOption, 'Generate and embed bitcode (only macOS/iOS/tvOS)', False),
    OptionKey('b_vscrt'): BaseOption(options.UserComboOption, 'VS run-time library type to use.', 'from_buildtype',
                                     choices=MSCRT_VALS + ['from_buildtype', 'static_from_buildtype']),
}

base_options = {key: base_opt.init_option(key) for key, base_opt in BASE_OPTIONS.items()}

def option_enabled(boptions: T.Set[OptionKey], options: 'KeyedOptionDictType',
                   option: OptionKey) -> bool:
    try:
        if option not in boptions:
            return False
        ret = options.get_value(option)
        assert isinstance(ret, bool), 'must return bool'  # could also be str
        return ret
    except KeyError:
        return False


def get_option_value(options: 'KeyedOptionDictType', opt: OptionKey, fallback: '_T') -> '_T':
    """Get the value of an option, or the fallback value."""
    try:
        v: '_T' = options.get_value(opt)
    except (KeyError, AttributeError):
        return fallback

    assert isinstance(v, type(fallback)), f'Should have {type(fallback)!r} but was {type(v)!r}'
    # Mypy doesn't understand that the above assert ensures that v is type _T
    return v


def are_asserts_disabled(options: KeyedOptionDictType) -> bool:
    """Should debug assertions be disabled

    :param options: OptionDictionary
    :return: whether to disable assertions or not
    """
    return (options.get_value('b_ndebug') == 'true' or
            (options.get_value('b_ndebug') == 'if-release' and
             options.get_value('buildtype') in {'release', 'plain'}))


def get_base_compile_args(options: 'KeyedOptionDictType', compiler: 'Compiler', env: 'Environment') -> T.List[str]:
    args: T.List[str] = []
    try:
        if options.get_value(OptionKey('b_lto')):
            args.extend(compiler.get_lto_compile_args(
                threads=get_option_value(options, OptionKey('b_lto_threads'), 0),
                mode=get_option_value(options, OptionKey('b_lto_mode'), 'default')))
    except (KeyError, AttributeError):
        pass
    try:
        args += compiler.get_colorout_args(options.get_value(OptionKey('b_colorout')))
    except (KeyError, AttributeError):
        pass
    try:
        args += compiler.sanitizer_compile_args(options.get_value(OptionKey('b_sanitize')))
    except (KeyError, AttributeError):
        pass
    try:
        pgo_val = options.get_value(OptionKey('b_pgo'))
        if pgo_val == 'generate':
            args.extend(compiler.get_profile_generate_args())
        elif pgo_val == 'use':
            args.extend(compiler.get_profile_use_args())
    except (KeyError, AttributeError):
        pass
    try:
        if options.get_value(OptionKey('b_coverage')):
            args += compiler.get_coverage_args()
    except (KeyError, AttributeError):
        pass
    try:
        args += compiler.get_assert_args(are_asserts_disabled(options), env)
    except (KeyError, AttributeError):
        pass
    # This does not need a try...except
    if option_enabled(compiler.base_options, options, OptionKey('b_bitcode')):
        args.append('-fembed-bitcode')
    try:
        try:
            crt_val = options.get_value(OptionKey('b_vscrt'))
            buildtype = options.get_value(OptionKey('buildtype'))
            args += compiler.get_crt_compile_args(crt_val, buildtype)
        except AttributeError:
            pass
    except KeyError:
        pass
    return args

def get_base_link_args(options: 'KeyedOptionDictType', linker: 'Compiler',
                       is_shared_module: bool, build_dir: str) -> T.List[str]:
    args: T.List[str] = []
    try:
        if options.get_value('b_lto'):
            if options.get_value('werror'):
                args.extend(linker.get_werror_args())

            thinlto_cache_dir = None
            if get_option_value(options, OptionKey('b_thinlto_cache'), False):
                thinlto_cache_dir = get_option_value(options, OptionKey('b_thinlto_cache_dir'), '')
                if thinlto_cache_dir == '':
                    thinlto_cache_dir = os.path.join(build_dir, 'meson-private', 'thinlto-cache')
            args.extend(linker.get_lto_link_args(
                threads=get_option_value(options, OptionKey('b_lto_threads'), 0),
                mode=get_option_value(options, OptionKey('b_lto_mode'), 'default'),
                thinlto_cache_dir=thinlto_cache_dir))
    except (KeyError, AttributeError):
        pass
    try:
        args += linker.sanitizer_link_args(options.get_value('b_sanitize'))
    except (KeyError, AttributeError):
        pass
    try:
        pgo_val = options.get_value('b_pgo')
        if pgo_val == 'generate':
            args.extend(linker.get_profile_generate_args())
        elif pgo_val == 'use':
            args.extend(linker.get_profile_use_args())
    except (KeyError, AttributeError):
        pass
    try:
        if options.get_value('b_coverage'):
            args += linker.get_coverage_link_args()
    except (KeyError, AttributeError):
        pass

    as_needed = option_enabled(linker.base_options, options, OptionKey('b_asneeded'))
    bitcode = option_enabled(linker.base_options, options, OptionKey('b_bitcode'))
    # Shared modules cannot be built with bitcode_bundle because
    # -bitcode_bundle is incompatible with -undefined and -bundle
    if bitcode and not is_shared_module:
        args.extend(linker.bitcode_args())
    elif as_needed:
        # -Wl,-dead_strip_dylibs is incompatible with bitcode
        args.extend(linker.get_asneeded_args())

    # Apple's ld (the only one that supports bitcode) does not like -undefined
    # arguments or -headerpad_max_install_names when bitcode is enabled
    if not bitcode:
        args.extend(linker.headerpad_args())
        if (not is_shared_module and
                option_enabled(linker.base_options, options, OptionKey('b_lundef'))):
            args.extend(linker.no_undefined_link_args())
        else:
            args.extend(linker.get_allow_undefined_link_args())

    try:
        try:
            crt_val = options.get_value(OptionKey('b_vscrt'))
            buildtype = options.get_value(OptionKey('buildtype'))
            args += linker.get_crt_link_args(crt_val, buildtype)
        except AttributeError:
            pass
    except KeyError:
        pass
    return args


class CrossNoRunException(MesonException):
    pass

@dataclass
class RunResult(HoldableObject):
    compiled: bool
    returncode: int = 999
    stdout: str = 'UNDEFINED'
    stderr: str = 'UNDEFINED'
    cached: bool = False


@dataclass
class CompileResult(HoldableObject):

    """The result of Compiler.compiles (and friends)."""

    stdout: str
    stderr: str
    command: T.List[str]
    returncode: int
    input_name: str
    output_name: T.Optional[str] = field(default=None, init=False)
    cached: bool = field(default=False, init=False)


class Compiler(HoldableObject, metaclass=abc.ABCMeta):
    # Libraries to ignore in find_library() since they are provided by the
    # compiler or the C library. Currently only used for MSVC.
    ignore_libs: T.List[str] = []
    # Libraries that are internal compiler implementations, and must not be
    # manually searched.
    internal_libs: T.List[str] = []

    LINKER_PREFIX: T.Union[None, str, T.List[str]] = None
    INVOKES_LINKER = True

    language: str
    id: str
    warn_args: T.Dict[str, T.List[str]]
    mode = 'COMPILER'

    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str,
                 for_machine: MachineChoice, info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None, is_cross: bool = False):
        self.exelist = ccache + exelist
        self.exelist_no_ccache = exelist
        self.file_suffixes = lang_suffixes[self.language]
        self.can_compile_suffixes = set(self.file_suffixes)
        self.default_suffix = self.file_suffixes[0]
        self.version = version
        self.full_version = full_version
        self.for_machine = for_machine
        self.base_options: T.Set[OptionKey] = set()
        self.linker = linker
        self.info = info
        self.is_cross = is_cross
        self.modes: T.List[Compiler] = []

    def __repr__(self) -> str:
        repr_str = "<{0}: v{1} `{2}`>"
        return repr_str.format(self.__class__.__name__, self.version,
                               ' '.join(self.exelist))

    @lru_cache(maxsize=None)
    def can_compile(self, src: 'mesonlib.FileOrString') -> bool:
        if isinstance(src, mesonlib.File):
            src = src.fname
        suffix = os.path.splitext(src)[1]
        if suffix != '.C':
            suffix = suffix.lower()
        return bool(suffix) and suffix[1:] in self.can_compile_suffixes

    def get_id(self) -> str:
        return self.id

    def get_modes(self) -> T.List[Compiler]:
        return self.modes

    def get_linker_id(self) -> str:
        # There is not guarantee that we have a dynamic linker instance, as
        # some languages don't have separate linkers and compilers. In those
        # cases return the compiler id
        try:
            return self.linker.id
        except AttributeError:
            return self.id

    def get_version_string(self) -> str:
        details = [self.id, self.version]
        if self.full_version:
            details += ['"%s"' % (self.full_version)]
        return '(%s)' % (' '.join(details))

    def get_language(self) -> str:
        return self.language

    @classmethod
    def get_display_language(cls) -> str:
        return cls.language.capitalize()

    def get_default_suffix(self) -> str:
        return self.default_suffix

    def get_define(self, dname: str, prefix: str, env: 'Environment',
                   extra_args: T.Union[T.List[str], T.Callable[[CompileCheckMode], T.List[str]]],
                   dependencies: T.List['Dependency'],
                   disable_cache: bool = False) -> T.Tuple[str, bool]:
        raise EnvironmentException('%s does not support get_define ' % self.get_id())

    def compute_int(self, expression: str, low: T.Optional[int], high: T.Optional[int],
                    guess: T.Optional[int], prefix: str, env: 'Environment', *,
                    extra_args: T.Union[None, T.List[str], T.Callable[[CompileCheckMode], T.List[str]]],
                    dependencies: T.Optional[T.List['Dependency']]) -> int:
        raise EnvironmentException('%s does not support compute_int ' % self.get_id())

    def compute_parameters_with_absolute_paths(self, parameter_list: T.List[str],
                                               build_dir: str) -> T.List[str]:
        raise EnvironmentException('%s does not support compute_parameters_with_absolute_paths ' % self.get_id())

    def has_members(self, typename: str, membernames: T.List[str],
                    prefix: str, env: 'Environment', *,
                    extra_args: T.Union[None, T.List[str], T.Callable[[CompileCheckMode], T.List[str]]] = None,
                    dependencies: T.Optional[T.List['Dependency']] = None) -> T.Tuple[bool, bool]:
        raise EnvironmentException('%s does not support has_member(s) ' % self.get_id())

    def has_type(self, typename: str, prefix: str, env: 'Environment',
                 extra_args: T.Union[T.List[str], T.Callable[[CompileCheckMode], T.List[str]]], *,
                 dependencies: T.Optional[T.List['Dependency']] = None) -> T.Tuple[bool, bool]:
        raise EnvironmentException('%s does not support has_type ' % self.get_id())

    def symbols_have_underscore_prefix(self, env: 'Environment') -> bool:
        raise EnvironmentException('%s does not support symbols_have_underscore_prefix ' % self.get_id())

    def get_exelist(self, ccache: bool = True) -> T.List[str]:
        return self.exelist.copy() if ccache else self.exelist_no_ccache.copy()

    def get_linker_exelist(self) -> T.List[str]:
        return self.linker.get_exelist() if self.linker else self.get_exelist()

    @abc.abstractmethod
    def get_output_args(self, outputname: str) -> T.List[str]:
        pass

    def get_linker_output_args(self, outputname: str) -> T.List[str]:
        return self.linker.get_output_args(outputname)

    def get_linker_search_args(self, dirname: str) -> T.List[str]:
        return self.linker.get_search_args(dirname)

    def get_builtin_define(self, define: str) -> T.Optional[str]:
        raise EnvironmentException('%s does not support get_builtin_define.' % self.id)

    def has_builtin_define(self, define: str) -> bool:
        raise EnvironmentException('%s does not support has_builtin_define.' % self.id)

    def get_always_args(self) -> T.List[str]:
        return []

    def can_linker_accept_rsp(self) -> bool:
        """
        Determines whether the linker can accept arguments using the @rsp syntax.
        """
        return self.linker.get_accepts_rsp()

    def get_linker_always_args(self) -> T.List[str]:
        return self.linker.get_always_args()

    def get_linker_lib_prefix(self) -> str:
        return self.linker.get_lib_prefix()

    def gen_import_library_args(self, implibname: str) -> T.List[str]:
        """
        Used only on Windows for libraries that need an import library.
        This currently means C, C++, Fortran.
        """
        return []

    def create_option(self, option_type: T.Type[UserOptionType], option_key: OptionKey, *args: T.Any, **kwargs: T.Any) -> T.Tuple[OptionKey, UserOptionType]:
        return option_key, option_type(f'{self.language}_{option_key.name}', *args, **kwargs)

    @staticmethod
    def update_options(options: MutableKeyedOptionDictType, *args: T.Tuple[OptionKey, UserOptionType]) -> MutableKeyedOptionDictType:
        options.update(args)
        return options

    def get_options(self) -> 'MutableKeyedOptionDictType':
        return {}

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        return []

    def get_option_link_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        return self.linker.get_option_args(options)

    def check_header(self, hname: str, prefix: str, env: 'Environment', *,
                     extra_args: T.Union[None, T.List[str], T.Callable[[CompileCheckMode], T.List[str]]] = None,
                     dependencies: T.Optional[T.List['Dependency']] = None) -> T.Tuple[bool, bool]:
        """Check that header is usable.

        Returns a two item tuple of bools. The first bool is whether the
        check succeeded, the second is whether the result was cached (True)
        or run fresh (False).
        """
        raise EnvironmentException('Language %s does not support header checks.' % self.get_display_language())

    def has_header(self, hname: str, prefix: str, env: 'Environment', *,
                   extra_args: T.Union[None, T.List[str], T.Callable[[CompileCheckMode], T.List[str]]] = None,
                   dependencies: T.Optional[T.List['Dependency']] = None,
                   disable_cache: bool = False) -> T.Tuple[bool, bool]:
        """Check that header is exists.

        This check will return true if the file exists, even if it contains:

        ```c
        # error "You thought you could use this, LOLZ!"
        ```

        Use check_header if your header only works in some cases.

        Returns a two item tuple of bools. The first bool is whether the
        check succeeded, the second is whether the result was cached (True)
        or run fresh (False).
        """
        raise EnvironmentException('Language %s does not support header checks.' % self.get_display_language())

    def has_header_symbol(self, hname: str, symbol: str, prefix: str,
                          env: 'Environment', *,
                          extra_args: T.Union[None, T.List[str], T.Callable[[CompileCheckMode], T.List[str]]] = None,
                          dependencies: T.Optional[T.List['Dependency']] = None) -> T.Tuple[bool, bool]:
        raise EnvironmentException('Language %s does not support header symbol checks.' % self.get_display_language())

    def run(self, code: 'mesonlib.FileOrString', env: 'Environment',
            extra_args: T.Union[T.List[str], T.Callable[[CompileCheckMode], T.List[str]], None] = None,
            dependencies: T.Optional[T.List['Dependency']] = None,
            run_env: T.Optional[T.Dict[str, str]] = None,
            run_cwd: T.Optional[str] = None) -> RunResult:
        need_exe_wrapper = env.need_exe_wrapper(self.for_machine)
        if need_exe_wrapper and not env.has_exe_wrapper():
            raise CrossNoRunException('Can not run test applications in this cross environment.')
        with self._build_wrapper(code, env, extra_args, dependencies, mode=CompileCheckMode.LINK, want_output=True) as p:
            if p.returncode != 0:
                mlog.debug(f'Could not compile test file {p.input_name}: {p.returncode}\n')
                return RunResult(False)
            if need_exe_wrapper:
                cmdlist = env.exe_wrapper.get_command() + [p.output_name]
            else:
                cmdlist = [p.output_name]
            try:
                pe, so, se = mesonlib.Popen_safe(cmdlist, env=run_env, cwd=run_cwd)
            except Exception as e:
                mlog.debug(f'Could not run: {cmdlist} (error: {e})\n')
                return RunResult(False)

        mlog.debug('Program stdout:\n')
        mlog.debug(so)
        mlog.debug('Program stderr:\n')
        mlog.debug(se)
        return RunResult(True, pe.returncode, so, se)

    # Caching run() in general seems too risky (no way to know what the program
    # depends on), but some callers know more about the programs they intend to
    # run.
    # For now we just accept code as a string, as that's what internal callers
    # need anyway. If we wanted to accept files, the cache key would need to
    # include mtime.
    def cached_run(self, code: str, env: 'Environment', *,
                   extra_args: T.Union[T.List[str], T.Callable[[CompileCheckMode], T.List[str]], None] = None,
                   dependencies: T.Optional[T.List['Dependency']] = None) -> RunResult:
        run_check_cache = env.coredata.run_check_cache
        args = self.build_wrapper_args(env, extra_args, dependencies, CompileCheckMode('link'))
        key = (code, tuple(args))
        if key in run_check_cache:
            p = run_check_cache[key]
            p.cached = True
            mlog.debug('Using cached run result:')
            mlog.debug('Code:\n', code)
            mlog.debug('Args:\n', extra_args)
            mlog.debug('Cached run returncode:\n', p.returncode)
            mlog.debug('Cached run stdout:\n', p.stdout)
            mlog.debug('Cached run stderr:\n', p.stderr)
        else:
            p = self.run(code, env, extra_args=extra_args, dependencies=dependencies)
            run_check_cache[key] = p
        return p

    def sizeof(self, typename: str, prefix: str, env: 'Environment', *,
               extra_args: T.Union[None, T.List[str], T.Callable[[CompileCheckMode], T.List[str]]] = None,
               dependencies: T.Optional[T.List['Dependency']] = None) -> T.Tuple[int, bool]:
        raise EnvironmentException('Language %s does not support sizeof checks.' % self.get_display_language())

    def alignment(self, typename: str, prefix: str, env: 'Environment', *,
                  extra_args: T.Optional[T.List[str]] = None,
                  dependencies: T.Optional[T.List['Dependency']] = None) -> T.Tuple[int, bool]:
        raise EnvironmentException('Language %s does not support alignment checks.' % self.get_display_language())

    def has_function(self, funcname: str, prefix: str, env: 'Environment', *,
                     extra_args: T.Optional[T.List[str]] = None,
                     dependencies: T.Optional[T.List['Dependency']] = None) -> T.Tuple[bool, bool]:
        """See if a function exists.

        Returns a two item tuple of bools. The first bool is whether the
        check succeeded, the second is whether the result was cached (True)
        or run fresh (False).
        """
        raise EnvironmentException('Language %s does not support function checks.' % self.get_display_language())

    @classmethod
    def _unix_args_to_native(cls, args: T.List[str], info: MachineInfo) -> T.List[str]:
        "Always returns a copy that can be independently mutated"
        return args.copy()

    def unix_args_to_native(self, args: T.List[str]) -> T.List[str]:
        return self._unix_args_to_native(args, self.info)

    @classmethod
    def native_args_to_unix(cls, args: T.List[str]) -> T.List[str]:
        "Always returns a copy that can be independently mutated"
        return args.copy()

    def find_library(self, libname: str, env: 'Environment', extra_dirs: T.List[str],
                     libtype: LibType = LibType.PREFER_SHARED, lib_prefix_warning: bool = True) -> T.Optional[T.List[str]]:
        raise EnvironmentException(f'Language {self.get_display_language()} does not support library finding.')

    def get_library_naming(self, env: 'Environment', libtype: LibType,
                           strict: bool = False) -> T.Optional[T.Tuple[str, ...]]:
        raise EnvironmentException(
            'Language {} does not support get_library_naming.'.format(
                self.get_display_language()))

    def get_program_dirs(self, env: 'Environment') -> T.List[str]:
        return []

    def has_multi_arguments(self, args: T.List[str], env: 'Environment') -> T.Tuple[bool, bool]:
        """Checks if the compiler has all of the arguments.

        :returns:
            A tuple of (bool, bool). The first value is whether the check
            succeeded, and the second is whether it was retrieved from a cache
        """
        raise EnvironmentException(
            'Language {} does not support has_multi_arguments.'.format(
                self.get_display_language()))

    def has_multi_link_arguments(self, args: T.List[str], env: 'Environment') -> T.Tuple[bool, bool]:
        """Checks if the linker has all of the arguments.

        :returns:
            A tuple of (bool, bool). The first value is whether the check
            succeeded, and the second is whether it was retrieved from a cache
        """
        return self.linker.has_multi_arguments(args, env)

    def _get_compile_output(self, dirname: str, mode: CompileCheckMode) -> str:
        assert mode != CompileCheckMode.PREPROCESS, 'In pre-processor mode, the output is sent to stdout and discarded'
        # Extension only matters if running results; '.exe' is
        # guaranteed to be executable on every platform.
        if mode == CompileCheckMode.LINK:
            suffix = 'exe'
        else:
            suffix = 'obj'
        return os.path.join(dirname, 'output.' + suffix)

    def get_compiler_args_for_mode(self, mode: CompileCheckMode) -> T.List[str]:
        args: T.List[str] = []
        args += self.get_always_args()
        if mode is CompileCheckMode.COMPILE:
            args += self.get_compile_only_args()
        elif mode is CompileCheckMode.PREPROCESS:
            args += self.get_preprocess_only_args()
        else:
            assert mode is CompileCheckMode.LINK
        return args

    def compiler_args(self, args: T.Optional[T.Iterable[str]] = None) -> CompilerArgs:
        """Return an appropriate CompilerArgs instance for this class."""
        return CompilerArgs(self, args)

    @contextlib.contextmanager
    def compile(self, code: 'mesonlib.FileOrString',
                extra_args: T.Union[None, CompilerArgs, T.List[str]] = None,
                *, mode: CompileCheckMode = CompileCheckMode.LINK, want_output: bool = False,
                temp_dir: T.Optional[str] = None) -> T.Iterator[CompileResult]:
        # TODO: there isn't really any reason for this to be a contextmanager

        if mode == CompileCheckMode.PREPROCESS:
            assert not want_output, 'In pre-processor mode, the output is sent to stdout and discarded'

        if extra_args is None:
            extra_args = []

        with TemporaryDirectoryWinProof(dir=temp_dir) as tmpdirname:
            no_ccache = False
            if isinstance(code, str):
                srcname = os.path.join(tmpdirname,
                                       'testfile.' + self.default_suffix)
                with open(srcname, 'w', encoding='utf-8') as ofile:
                    ofile.write(code)
                    if not code.endswith('\n'):
                        ofile.write('\n')
                # ccache would result in a cache miss
                no_ccache = True
                code_debug = f'Code:\n{code}'
            else:
                srcname = code.fname
                code_debug = f'Source file: {srcname}'

            # Construct the compiler command-line
            commands = self.compiler_args()
            commands.append(srcname)

            # Preprocess mode outputs to stdout, so no output args
            if mode != CompileCheckMode.PREPROCESS:
                output = self._get_compile_output(tmpdirname, mode)
                commands += self.get_output_args(output)
            commands.extend(self.get_compiler_args_for_mode(CompileCheckMode(mode)))

            # extra_args must be last because it could contain '/link' to
            # pass args to VisualStudio's linker. In that case everything
            # in the command line after '/link' is given to the linker.
            if extra_args:
                commands += extra_args
            # Generate full command-line with the exelist
            command_list = self.get_exelist(ccache=not no_ccache) + commands.to_native()
            mlog.debug('Running compile:')
            mlog.debug('Working directory: ', tmpdirname)
            mlog.debug(code_debug)
            os_env = os.environ.copy()
            os_env['LC_ALL'] = 'C'
            if no_ccache:
                os_env['CCACHE_DISABLE'] = '1'
            p, stdo, stde = Popen_safe_logged(command_list, msg='Command line', cwd=tmpdirname, env=os_env)

            result = CompileResult(stdo, stde, command_list, p.returncode, input_name=srcname)
            if want_output:
                result.output_name = output
            yield result

    @contextlib.contextmanager
    def cached_compile(self, code: 'mesonlib.FileOrString', cdata: coredata.CoreData, *,
                       extra_args: T.Union[None, T.List[str], CompilerArgs] = None,
                       mode: CompileCheckMode = CompileCheckMode.LINK,
                       temp_dir: T.Optional[str] = None) -> T.Iterator[CompileResult]:
        # TODO: There's isn't really any reason for this to be a context manager

        # Calculate the key
        textra_args: T.Tuple[str, ...] = tuple(extra_args) if extra_args is not None else tuple()
        key: coredata.CompilerCheckCacheKey = (tuple(self.exelist), self.version, code, textra_args, mode)

        # Check if not cached, and generate, otherwise get from the cache
        if key in cdata.compiler_check_cache:
            p = cdata.compiler_check_cache[key]
            p.cached = True
            mlog.debug('Using cached compile:')
            mlog.debug('Cached command line: ', ' '.join(p.command), '\n')
            mlog.debug('Code:\n', code)
            mlog.debug('Cached compiler stdout:\n', p.stdout)
            mlog.debug('Cached compiler stderr:\n', p.stderr)
            yield p
        else:
            with self.compile(code, extra_args=extra_args, mode=mode, want_output=False, temp_dir=temp_dir) as p:
                cdata.compiler_check_cache[key] = p
                yield p

    def get_colorout_args(self, colortype: str) -> T.List[str]:
        # TODO: colortype can probably be an emum
        return []

    # Some compilers (msvc) write debug info to a separate file.
    # These args specify where it should be written.
    def get_compile_debugfile_args(self, rel_obj: str, pch: bool = False) -> T.List[str]:
        return []

    def should_link_pch_object(self) -> bool:
        return False

    def get_link_debugfile_name(self, targetfile: str) -> T.Optional[str]:
        return self.linker.get_debugfile_name(targetfile)

    def get_link_debugfile_args(self, targetfile: str) -> T.List[str]:
        return self.linker.get_debugfile_args(targetfile)

    def get_std_shared_lib_link_args(self) -> T.List[str]:
        return self.linker.get_std_shared_lib_args()

    def get_std_shared_module_link_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        return self.linker.get_std_shared_module_args(options)

    def get_link_whole_for(self, args: T.List[str]) -> T.List[str]:
        return self.linker.get_link_whole_for(args)

    def get_allow_undefined_link_args(self) -> T.List[str]:
        return self.linker.get_allow_undefined_args()

    def no_undefined_link_args(self) -> T.List[str]:
        return self.linker.no_undefined_args()

    def get_instruction_set_args(self, instruction_set: str) -> T.Optional[T.List[str]]:
        """Compiler arguments needed to enable the given instruction set.

        Return type ay be an empty list meaning nothing needed or None
        meaning the given set is not supported.
        """
        return None

    def build_rpath_args(self, env: 'Environment', build_dir: str, from_dir: str,
                         rpath_paths: T.Tuple[str, ...], build_rpath: str,
                         install_rpath: str) -> T.Tuple[T.List[str], T.Set[bytes]]:
        return self.linker.build_rpath_args(
            env, build_dir, from_dir, rpath_paths, build_rpath, install_rpath)

    def get_archive_name(self, filename: str) -> str:
        return self.linker.get_archive_name(filename)

    def get_command_to_archive_shlib(self) -> T.List[str]:
        if not self.linker:
            return []
        return self.linker.get_command_to_archive_shlib()

    def thread_flags(self, env: 'Environment') -> T.List[str]:
        return []

    def thread_link_flags(self, env: 'Environment') -> T.List[str]:
        return self.linker.thread_flags(env)

    def openmp_flags(self, env: Environment) -> T.List[str]:
        raise EnvironmentException('Language %s does not support OpenMP flags.' % self.get_display_language())

    def openmp_link_flags(self, env: Environment) -> T.List[str]:
        return self.openmp_flags(env)

    def language_stdlib_only_link_flags(self, env: 'Environment') -> T.List[str]:
        return []

    def gnu_symbol_visibility_args(self, vistype: str) -> T.List[str]:
        return []

    def get_win_subsystem_args(self, value: str) -> T.List[str]:
        # By default the dynamic linker is going to return an empty
        # array in case it either doesn't support Windows subsystems
        # or does not target Windows
        return self.linker.get_win_subsystem_args(value)

    def has_func_attribute(self, name: str, env: 'Environment') -> T.Tuple[bool, bool]:
        raise EnvironmentException(
            f'Language {self.get_display_language()} does not support function attributes.')

    def get_pic_args(self) -> T.List[str]:
        m = 'Language {} does not support position-independent code'
        raise EnvironmentException(m.format(self.get_display_language()))

    def get_pie_args(self) -> T.List[str]:
        m = 'Language {} does not support position-independent executable'
        raise EnvironmentException(m.format(self.get_display_language()))

    def get_pie_link_args(self) -> T.List[str]:
        return self.linker.get_pie_args()

    @staticmethod
    def get_argument_syntax() -> str:
        """Returns the argument family type.

        Compilers fall into families if they try to emulate the command line
        interface of another compiler. For example, clang is in the GCC family
        since it accepts most of the same arguments as GCC. ICL (ICC on
        windows) is in the MSVC family since it accepts most of the same
        arguments as MSVC.
        """
        return 'other'

    def get_profile_generate_args(self) -> T.List[str]:
        raise EnvironmentException(
            '%s does not support get_profile_generate_args ' % self.get_id())

    def get_profile_use_args(self) -> T.List[str]:
        raise EnvironmentException(
            '%s does not support get_profile_use_args ' % self.get_id())

    def remove_linkerlike_args(self, args: T.List[str]) -> T.List[str]:
        rm_exact = ('-headerpad_max_install_names',)
        rm_prefixes = ('-Wl,', '-L',)
        rm_next = ('-L', '-framework',)
        ret: T.List[str] = []
        iargs = iter(args)
        for arg in iargs:
            # Remove this argument
            if arg in rm_exact:
                continue
            # If the argument starts with this, but is not *exactly* this
            # f.ex., '-L' should match ['-Lfoo'] but not ['-L', 'foo']
            if arg.startswith(rm_prefixes) and arg not in rm_prefixes:
                continue
            # Ignore this argument and the one after it
            if arg in rm_next:
                next(iargs)
                continue
            ret.append(arg)
        return ret

    def get_lto_compile_args(self, *, threads: int = 0, mode: str = 'default') -> T.List[str]:
        return []

    def get_lto_link_args(self, *, threads: int = 0, mode: str = 'default',
                          thinlto_cache_dir: T.Optional[str] = None) -> T.List[str]:
        return self.linker.get_lto_args()

    def sanitizer_compile_args(self, value: str) -> T.List[str]:
        return []

    def sanitizer_link_args(self, value: str) -> T.List[str]:
        return self.linker.sanitizer_args(value)

    def get_asneeded_args(self) -> T.List[str]:
        return self.linker.get_asneeded_args()

    def headerpad_args(self) -> T.List[str]:
        return self.linker.headerpad_args()

    def bitcode_args(self) -> T.List[str]:
        return self.linker.bitcode_args()

    def get_optimization_link_args(self, optimization_level: str) -> T.List[str]:
        return self.linker.get_optimization_link_args(optimization_level)

    def get_soname_args(self, env: 'Environment', prefix: str, shlib_name: str,
                        suffix: str, soversion: str,
                        darwin_versions: T.Tuple[str, str]) -> T.List[str]:
        return self.linker.get_soname_args(
            env, prefix, shlib_name, suffix, soversion,
            darwin_versions)

    def get_target_link_args(self, target: 'BuildTarget') -> T.List[str]:
        return target.link_args

    def get_dependency_compile_args(self, dep: 'Dependency') -> T.List[str]:
        return dep.get_compile_args()

    def get_dependency_link_args(self, dep: 'Dependency') -> T.List[str]:
        return dep.get_link_args(self.get_language())

    @classmethod
    def use_linker_args(cls, linker: str, version: str) -> T.List[str]:
        """Get a list of arguments to pass to the compiler to set the linker.
        """
        return []

    def get_coverage_args(self) -> T.List[str]:
        return []

    def get_coverage_link_args(self) -> T.List[str]:
        return self.linker.get_coverage_args()

    def get_assert_args(self, disable: bool, env: 'Environment') -> T.List[str]:
        """Get arguments to enable or disable assertion.

        :param disable: Whether to disable assertions
        :return: A list of string arguments for this compiler
        """
        return []

    def get_crt_val(self, crt_val: str, buildtype: str) -> str:
        if crt_val in MSCRT_VALS:
            return crt_val
        assert crt_val in {'from_buildtype', 'static_from_buildtype'}

        dbg = 'mdd'
        rel = 'md'
        if crt_val == 'static_from_buildtype':
            dbg = 'mtd'
            rel = 'mt'

        # Match what build type flags used to do.
        if buildtype == 'plain':
            return 'none'
        elif buildtype == 'debug':
            return dbg
        elif buildtype in {'debugoptimized', 'release', 'minsize'}:
            return rel
        else:
            assert buildtype == 'custom'
            raise EnvironmentException('Requested C runtime based on buildtype, but buildtype is "custom".')

    def get_crt_compile_args(self, crt_val: str, buildtype: str) -> T.List[str]:
        raise EnvironmentException('This compiler does not support Windows CRT selection')

    def get_crt_link_args(self, crt_val: str, buildtype: str) -> T.List[str]:
        raise EnvironmentException('This compiler does not support Windows CRT selection')

    def get_compile_only_args(self) -> T.List[str]:
        return []

    def get_preprocess_only_args(self) -> T.List[str]:
        raise EnvironmentException('This compiler does not have a preprocessor')

    def get_preprocess_to_file_args(self) -> T.List[str]:
        return self.get_preprocess_only_args()

    def get_default_include_dirs(self) -> T.List[str]:
        # TODO: This is a candidate for returning an immutable list
        return []

    def get_largefile_args(self) -> T.List[str]:
        '''Enable transparent large-file-support for 32-bit UNIX systems'''
        if not (self.get_argument_syntax() == 'msvc' or self.info.is_darwin()):
            # Enable large-file support unconditionally on all platforms other
            # than macOS and MSVC. macOS is now 64-bit-only so it doesn't
            # need anything special, and MSVC doesn't have automatic LFS.
            # You must use the 64-bit counterparts explicitly.
            # glibc, musl, and uclibc, and all BSD libcs support this. On Android,
            # support for transparent LFS is available depending on the version of
            # Bionic: https://github.com/android/platform_bionic#32-bit-abi-bugs
            # https://code.google.com/p/android/issues/detail?id=64613
            #
            # If this breaks your code, fix it! It's been 20+ years!
            return ['-D_FILE_OFFSET_BITS=64']
            # We don't enable -D_LARGEFILE64_SOURCE since that enables
            # transitionary features and must be enabled by programs that use
            # those features explicitly.
        return []

    def get_library_dirs(self, env: 'Environment',
                         elf_class: T.Optional[int] = None) -> T.List[str]:
        return []

    def get_return_value(self,
                         fname: str,
                         rtype: str,
                         prefix: str,
                         env: 'Environment',
                         extra_args: T.Optional[T.List[str]],
                         dependencies: T.Optional[T.List['Dependency']]) -> T.Union[str, int]:
        raise EnvironmentException(f'{self.id} does not support get_return_value')

    def find_framework(self,
                       name: str,
                       env: 'Environment',
                       extra_dirs: T.List[str],
                       allow_system: bool = True) -> T.Optional[T.List[str]]:
        raise EnvironmentException(f'{self.id} does not support find_framework')

    def find_framework_paths(self, env: 'Environment') -> T.List[str]:
        raise EnvironmentException(f'{self.id} does not support find_framework_paths')

    def attribute_check_func(self, name: str) -> str:
        raise EnvironmentException(f'{self.id} does not support attribute checks')

    def get_pch_suffix(self) -> str:
        raise EnvironmentException(f'{self.id} does not support pre compiled headers')

    def get_pch_name(self, name: str) -> str:
        raise EnvironmentException(f'{self.id} does not support pre compiled headers')

    def get_pch_use_args(self, pch_dir: str, header: str) -> T.List[str]:
        raise EnvironmentException(f'{self.id} does not support pre compiled headers')

    def get_has_func_attribute_extra_args(self, name: str) -> T.List[str]:
        raise EnvironmentException(f'{self.id} does not support function attributes')

    def name_string(self) -> str:
        return ' '.join(self.exelist)

    @abc.abstractmethod
    def sanity_check(self, work_dir: str, environment: 'Environment') -> None:
        """Check that this compiler actually works.

        This should provide a simple compile/link test. Something as simple as:
        ```python
        main(): return 0
        ```
        is good enough here.
        """

    def split_shlib_to_parts(self, fname: str) -> T.Tuple[T.Optional[str], str]:
        return None, fname

    def get_dependency_gen_args(self, outtarget: str, outfile: str) -> T.List[str]:
        return []

    def get_std_exe_link_args(self) -> T.List[str]:
        # TODO: is this a linker property?
        return []

    def get_include_args(self, path: str, is_system: bool) -> T.List[str]:
        return []

    def depfile_for_object(self, objfile: str) -> T.Optional[str]:
        return objfile + '.' + self.get_depfile_suffix()

    def get_depfile_suffix(self) -> str:
        raise EnvironmentException(f'{self.id} does not implement get_depfile_suffix')

    def get_no_stdinc_args(self) -> T.List[str]:
        """Arguments to turn off default inclusion of standard libraries."""
        return []

    def get_warn_args(self, level: str) -> T.List[str]:
        return []

    def get_werror_args(self) -> T.List[str]:
        return []

    @abc.abstractmethod
    def get_optimization_args(self, optimization_level: str) -> T.List[str]:
        pass

    def get_module_incdir_args(self) -> T.Tuple[str, ...]:
        raise EnvironmentException(f'{self.id} does not implement get_module_incdir_args')

    def get_module_outdir_args(self, path: str) -> T.List[str]:
        raise EnvironmentException(f'{self.id} does not implement get_module_outdir_args')

    def module_name_to_filename(self, module_name: str) -> str:
        raise EnvironmentException(f'{self.id} does not implement module_name_to_filename')

    def get_compiler_check_args(self, mode: CompileCheckMode) -> T.List[str]:
        """Arguments to pass the compiler and/or linker for checks.

        The default implementation turns off optimizations.

        Examples of things that go here:
          - extra arguments for error checking
          - Arguments required to make the compiler exit with a non-zero status
            when something is wrong.
        """
        return self.get_no_optimization_args()

    def get_no_optimization_args(self) -> T.List[str]:
        """Arguments to the compiler to turn off all optimizations."""
        return []

    def build_wrapper_args(self, env: 'Environment',
                           extra_args: T.Union[None, CompilerArgs, T.List[str], T.Callable[[CompileCheckMode], T.List[str]]],
                           dependencies: T.Optional[T.List['Dependency']],
                           mode: CompileCheckMode = CompileCheckMode.COMPILE) -> CompilerArgs:
        """Arguments to pass the build_wrapper helper.

        This generally needs to be set on a per-language basis. It provides
        a hook for languages to handle dependencies and extra args. The base
        implementation handles the most common cases, namely adding the
        check_arguments, unwrapping dependencies, and appending extra args.
        """
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
            args += d.get_compile_args()
            if mode is CompileCheckMode.LINK:
                # Add link flags needed to find dependencies
                args += d.get_link_args()

        if mode is CompileCheckMode.COMPILE:
            # Add DFLAGS from the env
            args += env.coredata.get_external_args(self.for_machine, self.language)
        elif mode is CompileCheckMode.LINK:
            # Add LDFLAGS from the env
            args += env.coredata.get_external_link_args(self.for_machine, self.language)
        # extra_args must override all other arguments, so we add them last
        args += extra_args
        return args

    @contextlib.contextmanager
    def _build_wrapper(self, code: 'mesonlib.FileOrString', env: 'Environment',
                       extra_args: T.Union[None, CompilerArgs, T.List[str], T.Callable[[CompileCheckMode], T.List[str]]] = None,
                       dependencies: T.Optional[T.List['Dependency']] = None,
                       mode: CompileCheckMode = CompileCheckMode.COMPILE, want_output: bool = False,
                       disable_cache: bool = False) -> T.Iterator[CompileResult]:
        """Helper for getting a cached value when possible.

        This method isn't meant to be called externally, it's mean to be
        wrapped by other methods like compiles() and links().
        """
        args = self.build_wrapper_args(env, extra_args, dependencies, mode)
        if disable_cache or want_output:
            with self.compile(code, extra_args=args, mode=mode, want_output=want_output, temp_dir=env.scratch_dir) as r:
                yield r
        else:
            with self.cached_compile(code, env.coredata, extra_args=args, mode=mode, temp_dir=env.scratch_dir) as r:
                yield r

    def compiles(self, code: 'mesonlib.FileOrString', env: 'Environment', *,
                 extra_args: T.Union[None, T.List[str], CompilerArgs, T.Callable[[CompileCheckMode], T.List[str]]] = None,
                 dependencies: T.Optional[T.List['Dependency']] = None,
                 mode: CompileCheckMode = CompileCheckMode.COMPILE,
                 disable_cache: bool = False) -> T.Tuple[bool, bool]:
        """Run a compilation or link test to see if code can be compiled/linked.

        :returns:
            A tuple of (bool, bool). The first value is whether the check
            succeeded, and the second is whether it was retrieved from a cache
        """
        with self._build_wrapper(code, env, extra_args, dependencies, mode, disable_cache=disable_cache) as p:
            return p.returncode == 0, p.cached

    def links(self, code: 'mesonlib.FileOrString', env: 'Environment', *,
              compiler: T.Optional['Compiler'] = None,
              extra_args: T.Union[None, T.List[str], CompilerArgs, T.Callable[[CompileCheckMode], T.List[str]]] = None,
              dependencies: T.Optional[T.List['Dependency']] = None,
              disable_cache: bool = False) -> T.Tuple[bool, bool]:
        if compiler:
            with compiler._build_wrapper(code, env, dependencies=dependencies, want_output=True) as r:
                objfile = mesonlib.File.from_absolute_file(r.output_name)
                return self.compiles(objfile, env, extra_args=extra_args,
                                     dependencies=dependencies, mode=CompileCheckMode.LINK, disable_cache=True)

        return self.compiles(code, env, extra_args=extra_args,
                             dependencies=dependencies, mode=CompileCheckMode.LINK, disable_cache=disable_cache)

    def get_feature_args(self, kwargs: DFeatures, build_to_src: str) -> T.List[str]:
        """Used by D for extra language features."""
        # TODO: using a TypeDict here would improve this
        raise EnvironmentException(f'{self.id} does not implement get_feature_args')

    def get_prelink_args(self, prelink_name: str, obj_list: T.List[str]) -> T.List[str]:
        raise EnvironmentException(f'{self.id} does not know how to do prelinking.')

    def rsp_file_syntax(self) -> 'RSPFileSyntax':
        """The format of the RSP file that this compiler supports.

        If `self.can_linker_accept_rsp()` returns True, then this needs to
        be implemented
        """
        return self.linker.rsp_file_syntax()

    def get_debug_args(self, is_debug: bool) -> T.List[str]:
        """Arguments required for a debug build."""
        return []

    def needs_static_linker(self) -> bool:
        raise NotImplementedError(f'There is no static linker for {self.language}')

    def get_preprocessor(self) -> Compiler:
        """Get compiler's preprocessor.
        """
        raise EnvironmentException(f'{self.get_id()} does not support preprocessor')

    def form_compileropt_key(self, basename: str) -> OptionKey:
        return OptionKey(f'{self.language}_{basename}', machine=self.for_machine)

def get_global_options(lang: str,
                       comp: T.Type[Compiler],
                       for_machine: MachineChoice,
                       env: 'Environment') -> dict[OptionKey, options.UserOption[T.Any]]:
    """Retrieve options that apply to all compilers for a given language."""
    description = f'Extra arguments passed to the {lang}'
    argkey = OptionKey(f'{lang}_args', machine=for_machine)
    largkey = argkey.evolve(f'{lang}_link_args')
    envkey = argkey.evolve(f'{lang}_env_args')

    comp_key = argkey if argkey in env.options else envkey

    comp_options = env.options.get(comp_key, [])
    link_options = env.options.get(largkey, [])

    cargs = options.UserArrayOption(
        f'{lang}_{argkey.name}',
        description + ' compiler',
        comp_options, split_args=True, allow_dups=True)

    largs = options.UserArrayOption(
        f'{lang}_{largkey.name}',
        description + ' linker',
        link_options, split_args=True, allow_dups=True)

    if comp.INVOKES_LINKER and comp_key == envkey:
        # If the compiler acts as a linker driver, and we're using the
        # environment variable flags for both the compiler and linker
        # arguments, then put the compiler flags in the linker flags as well.
        # This is how autotools works, and the env vars feature is for
        # autotools compatibility.
        largs.extend_value(comp_options)

    opts: dict[OptionKey, options.UserOption[T.Any]] = {argkey: cargs, largkey: largs}

    return opts
