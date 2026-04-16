# SPDX-License-Identifier: Apache-2.0
# Copyright 2012-2017 The Meson development team

from __future__ import annotations
from collections import defaultdict, deque, OrderedDict
from dataclasses import dataclass, field
from functools import lru_cache
import abc
import copy
import hashlib
import itertools, pathlib
import os
import pickle
import re
import textwrap
import typing as T

from . import coredata
from . import dependencies
from . import mlog
from . import programs
from .mesonlib import (
    HoldableObject, SecondLevelHolder,
    File, MesonException, MachineChoice, PerMachine, OrderedSet, listify,
    classify_unity_sources,
    get_filenames_templates_dict, substitute_values, has_path_sep,
    is_parent_path, relpath, PerMachineDefaultable,
    MesonBugException, EnvironmentVariables, pickle_load, lazy_property,
)
from .options import OptionKey

from .compilers import (
    is_header, is_object, is_source, clink_langs, sort_clink, all_languages,
    is_known_suffix, detect_static_linker, LANGUAGES_USING_LDFLAGS
)
from .interpreterbase import FeatureNew, FeatureDeprecated, SubProject

if T.TYPE_CHECKING:
    from typing_extensions import Literal, TypeAlias, TypedDict

    from .environment import Environment
    from ._typing import ImmutableListProtocol
    from .backend.backends import Backend
    from .compilers.compilers import Compiler, CompilerDict, Language
    from .interpreter.interpreter import SourceOutputs, Interpreter
    from .interpreter.interpreterobjects import Test, Doctest
    from .linkers.linkers import StaticLinker
    from .mesonlib import ExecutableSerialisation, FileMode, FileOrString
    from .mparser import BaseNode

    GeneratedTypes: TypeAlias = T.Union['CustomTarget', 'CustomTargetIndex', 'GeneratedList']
    LibTypes: TypeAlias = T.Union['SharedLibrary', 'StaticLibrary', 'CustomTarget', 'CustomTargetIndex']
    BuildTargetTypes: TypeAlias = T.Union['BuildTarget', 'CustomTarget', 'CustomTargetIndex']
    StaticTargetTypes: TypeAlias = T.Union['StaticLibrary', 'CustomTarget', 'CustomTargetIndex']
    ObjectTypes: TypeAlias = T.Union[str, 'File', 'ExtractedObjects', 'GeneratedTypes']
    AnyTargetType: TypeAlias = T.Union['Target', 'CustomTargetIndex']
    RustCrateType: TypeAlias = Literal['bin', 'lib', 'rlib', 'dylib', 'cdylib', 'staticlib', 'proc-macro']

    class DFeatures(TypedDict):

        unittest: bool
        debug: T.List[T.Union[str, int]]
        import_dirs: T.List[IncludeDirs]
        versions: T.List[T.Union[str, int]]

    class BuildTargetKeywordArguments(TypedDict, total=False):

        # The use of Sequence[str] here is intentional to get the generally
        # undesirable behavior of Sequence[str]

        build_by_default: bool
        build_rpath: str
        build_subdir: T.Optional[str]
        c_pch: T.Optional[T.Tuple[str, T.Optional[str]]]
        cpp_pch: T.Optional[T.Tuple[str, T.Optional[str]]]
        d_debug: T.List[T.Union[str, int]]
        d_import_dirs: T.List[IncludeDirs]
        d_module_versions: T.List[T.Union[str, int]]
        d_unittest: bool
        dependencies: T.List[dependencies.Dependency]
        extra_files: T.List[File]
        gnu_symbol_visibility: Literal['default', 'internal', 'hidden', 'protected', 'inlineshidden', '']
        implicit_include_directories: bool
        include_directories: T.List[IncludeDirs]
        install: bool
        install_dir: T.List[T.Union[str, bool]]
        install_mode: FileMode
        install_rpath: str
        install_tag: T.List[T.Optional[str]]
        language_args: T.DefaultDict[Language, T.List[str]]
        link_args: T.List[str]
        link_early_args: T.List[str]
        link_depends: T.List[T.Union[str, File, CustomTarget, CustomTargetIndex]]
        link_language: Language
        link_whole: T.List[StaticTargetTypes]
        link_with: T.List[BuildTargetTypes]
        name_prefix: T.Optional[str]
        name_suffix: T.Optional[str]
        native: MachineChoice
        override_options: T.Dict[OptionKey, str]
        resources: T.List[str]
        swift_interoperability_mode: Literal['c', 'cpp']
        swift_module_name: str
        rust_crate_type: RustCrateType
        rust_dependency_map: T.Dict[str, str]
        vala_gir: T.Optional[str]
        install_vala_gir: bool
        install_vala_gir_dir: T.Optional[str]
        vala_header: T.Optional[str]
        install_vala_header: bool
        install_vala_header_dir: T.Optional[str]
        vala_vapi: T.Optional[str]
        install_vala_vapi: bool
        install_vala_vapi_dir: T.Optional[str]
        win_subsystem: str

        _allow_no_sources: bool

    class ExecutableKeywordArguments(BuildTargetKeywordArguments, total=False):

        android_exe_type: T.Optional[Literal['application', 'executable']]
        implib: T.Optional[str]
        export_dynamic: bool
        pie: bool
        vs_module_defs: T.Union[str, File, CustomTarget, CustomTargetIndex]

    class SharedModuleKeywordArguments(BuildTargetKeywordArguments, total=False):

        vs_module_defs: T.Union[str, File, CustomTarget, CustomTargetIndex]

    class SharedLibraryKeywordArguments(SharedModuleKeywordArguments, total=False):

        version: str
        soversion: str
        darwin_versions: T.Tuple[str, str]
        shortname: str

    class StaticLibraryKeywordArguments(BuildTargetKeywordArguments, total=False):

        pic: bool
        prelink: bool

_T = T.TypeVar('_T')

DEFAULT_STATIC_LIBRARY_NAMES: T.Mapping[str, T.Tuple[str, str]] = {
    'unix': ('lib', 'a'),
    'windows': ('', 'lib'),
    'darwin': ('lib', 'a'),
    'cygwin': ('lib', 'a'),
}

DEFAULT_SHARED_LIBRARY_NAMES: T.Mapping[str, T.Tuple[str, str, str]] = {
    'unix': ('lib', 'so', ''),
    'windows': ('', 'dll', 'dll.lib'),
    'darwin': ('lib', 'dylib', ''),
    'cygwin': ('cyg', 'dll', 'dll.a'),
}

pch_kwargs = {'c_pch', 'cpp_pch'}

lang_arg_kwargs = {f'{lang}_args' for lang in all_languages}
lang_arg_kwargs |= {
    'd_import_dirs',
    'd_unittest',
    'd_module_versions',
    'd_debug',
}

vala_kwargs = {'vala_header', 'vala_gir', 'vala_vapi', 'install_vala_header', 'install_vala_gir', 'install_vala_vapi'}
rust_kwargs = {'rust_crate_type', 'rust_dependency_map'}
cs_kwargs = {'resources', 'cs_args'}
swift_kwargs = {'swift_interoperability_mode', 'swift_module_name'}

buildtarget_kwargs = {
    'build_by_default',
    'build_rpath',
    'build_subdir',
    'dependencies',
    'extra_files',
    'gui_app',
    'link_with',
    'link_whole',
    'link_args',
    'link_early_args',
    'link_depends',
    'implicit_include_directories',
    'include_directories',
    'install',
    'install_rpath',
    'install_dir',
    'install_mode',
    'install_tag',
    'name_prefix',
    'name_suffix',
    'native',
    'objects',
    'override_options',
    'sources',
    'gnu_symbol_visibility',
    'link_language',
    'win_subsystem',
}

known_build_target_kwargs = (
    buildtarget_kwargs |
    lang_arg_kwargs |
    pch_kwargs |
    vala_kwargs |
    rust_kwargs |
    cs_kwargs |
    swift_kwargs)

known_exe_kwargs = known_build_target_kwargs | {'implib', 'export_dynamic', 'pie', 'vs_module_defs', 'android_exe_type'}
known_shlib_kwargs = known_build_target_kwargs | {'version', 'soversion', 'vs_module_defs', 'darwin_versions', 'rust_abi', 'shortname'}
known_shmod_kwargs = known_build_target_kwargs | {'vs_module_defs', 'rust_abi'}
known_stlib_kwargs = known_build_target_kwargs | {'pic', 'prelink', 'rust_abi'}
known_jar_kwargs = known_exe_kwargs | {'main_class', 'java_resources'}

def _process_install_tag(install_tag: T.Optional[T.List[T.Optional[str]]],
                         num_outputs: int) -> T.List[T.Optional[str]]:
    _install_tag: T.List[T.Optional[str]]
    if not install_tag:
        _install_tag = [None] * num_outputs
    elif len(install_tag) == 1:
        _install_tag = install_tag * num_outputs
    else:
        _install_tag = install_tag
    return _install_tag


@lru_cache(maxsize=None)
def get_target_macos_dylib_install_name(ld: SharedLibrary) -> str:
    name = ['@rpath/', ld.prefix, ld.name]
    if ld.soversion is not None:
        name.append('.' + ld.soversion)
    name.append('.dylib')
    return ''.join(name)


class InvalidArguments(MesonException):
    pass

@dataclass(eq=False)
class DependencyOverride(HoldableObject):
    dep: dependencies.Dependency
    node: 'BaseNode'
    explicit: bool = True

@dataclass(eq=False)
class Headers(HoldableObject):
    sources: T.List[File]
    install_subdir: T.Optional[str]
    custom_install_dir: T.Optional[str]
    custom_install_mode: 'FileMode'
    subproject: str
    follow_symlinks: T.Optional[bool] = None
    install_tag: T.Optional[str] = None

    def __post_init__(self) -> None:
        if self.install_tag is None:
            self.install_tag = 'devel'

    # TODO: we really don't need any of these methods, but they're preserved to
    # keep APIs relying on them working.

    def set_install_subdir(self, subdir: str) -> None:
        self.install_subdir = subdir

    def get_install_subdir(self) -> T.Optional[str]:
        return self.install_subdir

    def get_sources(self) -> T.List[File]:
        return self.sources

    def get_custom_install_dir(self) -> T.Optional[str]:
        return self.custom_install_dir

    def get_custom_install_mode(self) -> 'FileMode':
        return self.custom_install_mode


@dataclass(eq=False)
class Man(HoldableObject):
    sources: T.List[File]
    custom_install_dir: T.Optional[str]
    custom_install_mode: 'FileMode'
    subproject: str
    locale: T.Optional[str]
    install_tag: T.Optional[str] = None

    def get_custom_install_dir(self) -> T.Optional[str]:
        return self.custom_install_dir

    def get_custom_install_mode(self) -> 'FileMode':
        return self.custom_install_mode

    def get_sources(self) -> T.List['File']:
        return self.sources


@dataclass(eq=False)
class EmptyDir(HoldableObject):
    path: str
    install_mode: 'FileMode'
    subproject: str
    install_tag: T.Optional[str] = None


@dataclass(eq=False)
class InstallDir(HoldableObject):
    source_subdir: str
    installable_subdir: str
    install_dir: str
    install_dir_name: str
    install_mode: 'FileMode'
    exclude: T.Tuple[T.Set[str], T.Set[str]]
    strip_directory: bool
    subproject: str
    from_source_dir: bool = True
    install_tag: T.Optional[str] = None
    follow_symlinks: T.Optional[bool] = None

@dataclass(eq=False)
class DepManifest:
    version: str
    license: T.List[str]
    license_files: T.List[T.Tuple[str, File]]
    subproject: str

    def license_mapping(self) -> T.List[T.Tuple[str, str]]:
        ret = []
        for ifilename, name in self.license_files:
            fname = os.path.join(*(x for x in pathlib.PurePath(os.path.normpath(name.fname)).parts if x != '..'))
            ret.append((ifilename, os.path.join(name.subdir, fname)))
        return ret

    def to_json(self) -> T.Dict[str, T.Union[str, T.List[str]]]:
        return {
            'version': self.version,
            'license': self.license,
            'license_files': [l[1] for l in self.license_mapping()],
        }


@dataclass(eq=False)
class BuildProject:
    name: str
    version: str
    project_args: PerMachine[T.Dict[str, T.List[str]]] = field(default_factory=lambda: PerMachine({}, {}))
    project_link_args: PerMachine[T.Dict[str, T.List[str]]] = field(default_factory=lambda: PerMachine({}, {}))


# literally everything isn't dataclass stuff
class Build:
    """A class that holds the status of one build including
    all dependencies and so on.
    """

    def __init__(self, environment: Environment):
        self.version = coredata.version
        self._def_files: T.Optional[T.List[str]] = None
        self.project_name = 'name of master project'
        self.project_version: T.Optional[str] = None
        self.environment = environment
        self.projects: T.Dict[SubProject, BuildProject] = {}
        self.targets: 'T.OrderedDict[str, T.Union[CustomTarget, BuildTarget]]' = OrderedDict()
        self.targetnames: T.Set[T.Tuple[str, str]] = set() # Set of executable names and their subdir
        self.global_args: PerMachine[T.Dict[str, T.List[str]]] = PerMachine({}, {})
        self.global_link_args: PerMachine[T.Dict[str, T.List[str]]] = PerMachine({}, {})
        self.tests: T.List['Test'] = []
        self.benchmarks: T.List['Test'] = []
        self.headers: T.List[Headers] = []
        self.man: T.List[Man] = []
        self.emptydir: T.List[EmptyDir] = []
        self.data: T.List[Data] = []
        self.symlinks: T.List[SymlinkData] = []
        self.static_linker: PerMachine[T.Optional[StaticLinker]] = PerMachine(None, None)
        self.subproject_dir = ''
        self.install_scripts: T.List['ExecutableSerialisation'] = []
        self.postconf_scripts: T.List['ExecutableSerialisation'] = []
        self.dist_scripts: T.List['ExecutableSerialisation'] = []
        self.install_dirs: T.List[InstallDir] = []
        self.dep_manifest_name: T.Optional[str] = None
        self.dep_manifest: T.Dict[str, DepManifest] = {}
        self.test_setups: T.Dict[str, TestSetup] = {}
        self.test_setup_default_name = None
        self.searched_programs: T.Set[str] = set() # The list of all programs that have been searched for.

        # If we are doing a cross build we need two caches, if we're doing a
        # build == host compilation the both caches should point to the same place.
        self.stdlibs: PerMachine[T.Dict[str, dependencies.Dependency]] = PerMachineDefaultable.default(
            environment.is_cross_build(), {}, {})
        self.dependency_overrides: PerMachine[T.Dict[T.Tuple, DependencyOverride]] = PerMachineDefaultable.default(
            environment.is_cross_build(), {}, {})

        self.devenv: T.List[EnvironmentVariables] = []
        self.modules: T.Set[str] = set()
        """Used to track which modules are enabled in all subprojects.

        Needed for tracking whether a modules options needs to be exposed to the user.
        """

    @property
    def def_files(self) -> T.List[str]:
        if self._def_files is None:
            raise MesonBugException('build.def_files has not been set yet')
        return self._def_files

    @def_files.setter
    def def_files(self, value: T.List[str]) -> None:
        if self._def_files is not None:
            raise MesonBugException('build.def_files already set')
        self._def_files = value

    def get_build_targets(self) -> OrderedDict[str, BuildTarget]:
        build_targets = OrderedDict()
        for name, t in self.targets.items():
            if isinstance(t, BuildTarget):
                build_targets[name] = t
        return build_targets

    def get_custom_targets(self) -> OrderedDict[str, CustomTarget]:
        custom_targets = OrderedDict()
        for name, t in self.targets.items():
            if isinstance(t, CustomTarget):
                custom_targets[name] = t
        return custom_targets

    def copy(self) -> Build:
        def copy_value(v: _T) -> _T:
            if isinstance(v, PerMachine):
                build = copy_value(v.build)
                host = build if v.build is v.host else copy_value(v.host)
                return PerMachine(build=build, host=host)
            if isinstance(v, (list, dict, set, OrderedDict)):
                return v.copy()
            return v

        other = Build.__new__(Build)
        for k, v in self.__dict__.items():
            other.__dict__[k] = copy_value(v)
        return other

    def merge(self, other: Build) -> None:
        for k, v in other.__dict__.items():
            # These are not modified in subprojects
            if k in {'global_args', 'global_link_args'}:
                continue

            if isinstance(v, PerMachine):
                dest = self.__dict__[k]
                dest.build = v.build
                dest.host = v.host
            else:
                self.__dict__[k] = v

    def ensure_static_linker(self, compiler: Compiler) -> None:
        for_machine = compiler.for_machine if self.environment.is_cross_build() else MachineChoice.HOST
        if self.static_linker[for_machine] is None and compiler.needs_static_linker():
            self.static_linker[for_machine] = detect_static_linker(self.environment, compiler)
            if not self.environment.is_cross_build():
                self.static_linker[MachineChoice.BUILD] = self.static_linker[MachineChoice.HOST]

    def get_project(self) -> str:
        return self.projects[SubProject('')].name

    def get_subproject_dir(self) -> str:
        return self.subproject_dir

    def get_targets(self) -> 'T.OrderedDict[str, T.Union[CustomTarget, BuildTarget]]':
        return self.targets

    def get_tests(self) -> T.List['Test']:
        return self.tests

    def get_benchmarks(self) -> T.List['Test']:
        return self.benchmarks

    def get_headers(self) -> T.List['Headers']:
        return self.headers

    def get_man(self) -> T.List['Man']:
        return self.man

    def get_data(self) -> T.List['Data']:
        return self.data

    def get_symlinks(self) -> T.List['SymlinkData']:
        return self.symlinks

    def get_emptydir(self) -> T.List['EmptyDir']:
        return self.emptydir

    def get_install_subdirs(self) -> T.List['InstallDir']:
        return self.install_dirs

    def get_global_args(self, compiler: 'Compiler', for_machine: 'MachineChoice') -> T.List[str]:
        d = self.global_args[for_machine]
        return d.get(compiler.get_language(), [])

    def get_project_args(self, compiler: 'Compiler', target: BuildTarget) -> T.List[str]:
        args = self.projects[target.subproject].project_args[target.for_machine]
        if not args:
            return []
        return args.get(compiler.get_language(), [])

    def get_global_link_args(self, compiler: 'Compiler', for_machine: 'MachineChoice') -> T.List[str]:
        d = self.global_link_args[for_machine]
        return d.get(compiler.get_language(), [])

    def get_project_link_args(self, compiler: 'Compiler', target: BuildTarget) -> T.List[str]:
        link_args = self.projects[target.subproject].project_link_args[target.for_machine]
        if not link_args:
            return []

        return link_args.get(compiler.get_language(), [])

@dataclass(eq=False)
class IncludeDirs(HoldableObject):

    """Internal representation of an include_directories call.

    :param curdir: The directory from which this IncludeDirs is declared
    :param incdirs: The paths relative to curdir the include
    :param is_system: Should these be treated as `-isystem` (true) or `-I
    :param extra_build_dirs: Extra build directories to include
    """

    curdir: str
    incdirs: T.List[str]
    is_system: bool
    extra_build_dirs: T.List[str] = field(default_factory=list)

    def abs_string_list(self, sourcedir: str, builddir: str) -> T.List[str]:
        """Convert IncludeDirs object to a list of absolute string paths.

        :param sourcedir: The absolute source directory
        :param builddir: The absolute build directory, option, build dir will not
            be added if this is unset
        :returns: A list of strings (without compiler argument)
        """
        strlist: T.List[str] = []
        for idir in self.incdirs:
            strlist.append(os.path.join(sourcedir, self.curdir, idir))
            strlist.append(os.path.join(builddir, self.curdir, idir))
        for idir in self.extra_build_dirs:
            strlist.append(os.path.join(builddir, self.curdir, idir))
        return strlist

    def rel_string_list(self, build_to_src: str, build_root: T.Optional[str] = None) -> T.List[str]:
        """Convert IncludeDirs object to a list of relative string paths.

        :param build_to_src: The relative path from the build dir to source dir
        :param build_root: The absolute build root. When provided, build
            directories that have not been created will be ignored. Default: None.
        :return: A list if strings (without compiler argument)
        """
        strlist: T.List[str] = []
        for idirs, add_src in [(self.incdirs, True), (self.extra_build_dirs, False)]:
            for idir in idirs:
                bld_dir = os.path.normpath(os.path.join(self.curdir, idir))
                if idir not in {'', '.'}:
                    expdir = bld_dir
                else:
                    expdir = self.curdir
                if build_root is None or os.path.isdir(os.path.join(build_root, expdir)):
                    strlist.append(bld_dir)
                if add_src:
                    strlist.append(os.path.normpath(os.path.join(build_to_src, expdir)))

        return strlist


@dataclass(eq=False)
class ExtractedObjects(HoldableObject):
    '''
    Holds a list of sources for which the objects must be extracted
    '''
    target: 'BuildTarget'
    srclist: T.List[File] = field(default_factory=list)
    genlist: T.List['GeneratedTypes'] = field(default_factory=list)
    objlist: T.List[T.Union[str, 'File', 'ExtractedObjects']] = field(default_factory=list)
    recursive: bool = True
    pch: bool = False

    def __repr__(self) -> str:
        r = '<{0} {1!r}: {2}>'
        return r.format(self.__class__.__name__, self.target.name, self.srclist)

    @staticmethod
    def get_sources(sources: T.Sequence['FileOrString'], generated_sources: T.Sequence['GeneratedTypes']) -> T.List['FileOrString']:
        # Merge sources and generated sources
        sources = list(sources)
        for gensrc in generated_sources:
            for s in gensrc.get_outputs():
                # We cannot know the path where this source will be generated,
                # but all we need here is the file extension to determine the
                # compiler.
                sources.append(s)

        # Filter out headers and all non-source files
        return [s for s in sources if is_source(s)]

    def classify_all_sources(self, sources: T.List[FileOrString], generated_sources: T.Sequence['GeneratedTypes']) -> T.Dict['Compiler', T.List['FileOrString']]:
        sources_ = self.get_sources(sources, generated_sources)
        return classify_unity_sources(self.target.compilers.values(), sources_)

    def check_unity_compatible(self) -> None:
        # Figure out if the extracted object list is compatible with a Unity
        # build. When we're doing a Unified build, we go through the sources,
        # and create a single source file from each subset of the sources that
        # can be compiled with a specific compiler. Then we create one object
        # from each unified source file. So for each compiler we can either
        # extra all its sources or none.
        cmpsrcs = self.classify_all_sources(self.target.sources, self.target.generated)
        extracted_cmpsrcs = self.classify_all_sources(self.srclist, self.genlist)

        for comp, srcs in extracted_cmpsrcs.items():
            if set(srcs) != set(cmpsrcs[comp]):
                raise MesonException('Single object files cannot be extracted '
                                     'in Unity builds. You can only extract all '
                                     'the object files for each compiler at once.')


@dataclass(eq=False, order=False)
class StructuredSources(HoldableObject):

    """A container for sources in languages that use filesystem hierarchy.

    Languages like Rust and Cython rely on the layout of files in the filesystem
    as part of the compiler implementation. This structure allows us to
    represent the required filesystem layout.
    """

    sources: T.DefaultDict[str, T.List[T.Union[File, GeneratedTypes]]] = field(
        default_factory=lambda: defaultdict(list))

    def __add__(self, other: StructuredSources) -> StructuredSources:
        sources = self.sources.copy()
        for k, v in other.sources.items():
            sources[k].extend(v)
        return StructuredSources(sources)

    def __bool__(self) -> bool:
        return bool(self.sources)

    def as_list(self) -> T.List[T.Union[File, GeneratedTypes]]:
        return list(itertools.chain.from_iterable(self.sources.values()))

    def needs_copy(self) -> bool:
        """Do we need to create a structure in the build directory.

        This allows us to avoid making copies if the structures exists in the
        source dir. Which could happen in situations where a generated source
        only exists in some configurations
        """
        for files in self.sources.values():
            for f in files:
                if isinstance(f, File):
                    if f.is_built:
                        return True
                else:
                    return True
        return False


@dataclass(eq=False)
class Target(HoldableObject, metaclass=abc.ABCMeta):

    name: str
    subdir: str
    subproject: 'SubProject'
    build_by_default: bool
    for_machine: MachineChoice
    environment: Environment
    install: bool = False
    build_always_stale: bool = False
    extra_files: T.List[File] = field(default_factory=list)
    build_subdir: str = ''

    @abc.abstractproperty
    def typename(self) -> str:
        pass

    @abc.abstractmethod
    def type_suffix(self) -> str:
        pass

    def __post_init__(self) -> None:
        # XXX: this should happen in the interpreter
        if has_path_sep(self.name):
            # Fix failing test 53 when this becomes an error.
            mlog.warning(textwrap.dedent(f'''\
                Target "{self.name}" has a path separator in its name.
                This is not supported, it can cause unexpected failures and will become
                a hard error in the future.'''))
        self.builddir = self.subdir
        if self.build_subdir:
            self.builddir = os.path.join(self.subdir, self.build_subdir)

    # dataclass comparators?
    def __lt__(self, other: object) -> bool:
        if not isinstance(other, Target):
            return NotImplemented
        return self.get_id() < other.get_id()

    def __le__(self, other: object) -> bool:
        if not isinstance(other, Target):
            return NotImplemented
        return self.get_id() <= other.get_id()

    def __gt__(self, other: object) -> bool:
        if not isinstance(other, Target):
            return NotImplemented
        return self.get_id() > other.get_id()

    def __ge__(self, other: object) -> bool:
        if not isinstance(other, Target):
            return NotImplemented
        return self.get_id() >= other.get_id()

    def get_custom_install_dir(self) -> T.List[T.Union[str, Literal[False]]]:
        raise NotImplementedError

    def get_basename(self) -> str:
        return self.name

    def get_subdir(self) -> str:
        return self.subdir

    def get_typename(self) -> str:
        return self.typename

    def get_build_subdir(self) -> str:
        return self.build_subdir

    def get_builddir(self) -> str:
        return self.builddir

    @staticmethod
    def _get_id_hash(target_id: str) -> str:
        # We don't really need cryptographic security here.
        # Small-digest hash function with unlikely collision is good enough.
        h = hashlib.sha256()
        h.update(target_id.encode(encoding='utf-8', errors='replace'))
        # This ID should be case-insensitive and should work in Visual Studio,
        # e.g. it should not start with leading '-'.
        return h.hexdigest()[:7]

    @staticmethod
    def construct_id_from_path(subdir: str, name: str, type_suffix: str) -> str:
        """Construct target ID from subdir, name and type suffix.

        This helper function is made public mostly for tests."""
        # This ID must also be a valid file name on all OSs.
        # It should also avoid shell metacharacters for obvious
        # reasons. '@' is not used as often as '_' in source code names.
        # In case of collisions consider using checksums.
        # FIXME replace with assert when slash in names is prohibited
        name_part = name.replace('/', '@').replace('\\', '@')
        assert not has_path_sep(type_suffix)
        my_id = name_part + type_suffix
        if subdir:
            subdir_part = Target._get_id_hash(subdir)
            # preserve myid for better debuggability
            return subdir_part + '@@' + my_id
        return my_id

    @lazy_property
    def id(self) -> str:
        return self.construct_id_from_path(
            self.builddir, self.name, self.type_suffix())

    def get_id(self) -> str:
        return self.id

    def get_override(self, name: str) -> T.Optional[str]:
        return self.raw_overrides.get(name, None)

    def is_linkable_target(self) -> bool:
        return False

    def get_outputs(self) -> T.List[str]:
        return []

    def should_install(self) -> bool:
        return False

class BuildTarget(Target):
    known_kwargs = known_build_target_kwargs

    rust_crate_type: RustCrateType

    # This set contains all the languages a linker can link natively
    # without extra flags. For instance, nvcc (cuda) can link C++
    # without injecting -lc++/-lstdc++, see
    #   https://github.com/mesonbuild/meson/issues/10570
    _MASK_LANGS: T.FrozenSet[T.Tuple[str, str]] = frozenset([
        # (language, linker)
        ('cpp', 'cuda'),
    ])

    def __init__(
            self,
            name: str,
            subdir: str,
            subproject: SubProject,
            for_machine: MachineChoice,
            sources: T.List['SourceOutputs'],
            structured_sources: T.Optional[StructuredSources],
            objects: T.List[ObjectTypes],
            environment: Environment,
            compilers: CompilerDict,
            kwargs: BuildTargetKeywordArguments):
        super().__init__(name, subdir, subproject, True, for_machine, environment, install=kwargs.get('install', False), build_subdir=kwargs.get('build_subdir', ''))
        # all_compilers is a reference to Interpreter.compilers, as such we
        # cannot mutate it inside build. Use a Mapping to get help from the
        # static type checker
        self.all_compilers: T.Mapping[Language, Compiler] = compilers
        self.compilers: CompilerDict = {}
        self.objects: T.List[ObjectTypes] = []
        self.structured_sources = structured_sources
        self.external_deps: T.List[dependencies.Dependency] = []
        self.include_dirs: T.List['IncludeDirs'] = []
        self.link_language: T.Optional[Language] = kwargs.get('link_language')
        self.link_targets: T.List[BuildTargetTypes] = []
        self.link_whole_targets: T.List[StaticTargetTypes] = []
        self.depend_files: T.List[File] = []
        self.link_depends: T.List[T.Union[File, BuildTargetTypes]] = []
        self.added_deps: T.Set[dependencies.Dependency] = set()
        self.name_prefix_set = False
        self.name_suffix_set = False
        self.filename = 'no_name'
        self.doctests: T.Optional[Doctest] = None
        # The debugging information file this target will generate
        self.debug_filename = None
        # The list of all files outputted by this target. Useful in cases such
        # as Vala which generates .vapi and .h besides the compiled output.
        self.outputs = [self.filename]
        self.pch: T.Dict[Language, T.Optional[T.Tuple[str, T.Optional[str]]]] = {}
        self.extra_args: T.DefaultDict[Language, T.List[str]] = kwargs.get('language_args', defaultdict(list))
        self.sources: T.List[File] = []
        # If the same source is defined multiple times, use it only once.
        self.seen_sources: T.Set[File] = set()
        self.generated: T.List['GeneratedTypes'] = []
        self.extra_files: T.List[File] = []
        self.d_features: DFeatures = {
            'debug': kwargs.get('d_debug', []),
            'import_dirs': kwargs.get('d_import_dirs', []),
            'versions': kwargs.get('d_module_versions', []),
            'unittest': kwargs.get('d_unittest', False),
        }
        self.pic = False
        self.pie = False
        self.both_lib: T.Optional[T.Union[StaticLibrary, SharedLibrary]] = None
        # Track build_rpath entries so we can remove them at install time
        self.rpath_dirs_to_remove: T.Set[bytes] = set()
        self.vala_header: T.Optional[str] = None
        self.vala_vapi: T.Optional[str] = None
        self.vala_gir: T.Optional[str] = None
        self.process_sourcelist(sources)
        # Objects can be:
        # 1. Preexisting objects provided by the user with the `objects:` kwarg
        # 2. Compiled objects created by and extracted from another target
        self.process_objectlist(objects)
        self.process_kwargs(kwargs)
        self.missing_languages = self.process_compilers()

        # self.link_targets and self.link_whole_targets contains libraries from
        # dependencies (see add_deps()). They have not been processed yet because
        # we have to call process_compilers() first and we need to process libraries
        # from link_with and link_whole first.
        # See https://github.com/mesonbuild/meson/pull/11957#issuecomment-1629243208.
        link_targets = self.extract_targets_as_list(kwargs, 'link_with')
        link_whole_targets = self.extract_targets_as_list(kwargs, 'link_whole')
        self.link_targets.clear()
        self.link_whole_targets.clear()
        self.link(link_targets)
        self.link_whole(link_whole_targets)
        self._set_vala_args(kwargs)

        if not any([[src for src in self.sources if not is_header(src)], self.generated, self.objects,
                    self.link_whole_targets, self.structured_sources, kwargs.pop('_allow_no_sources', False)]):
            mlog.warning(f'Build target {name} has no sources. '
                         'This was never supposed to be allowed but did because of a bug, '
                         'support will be removed in a future release of Meson')
        self.check_unknown_kwargs(kwargs)
        self.validate_install()
        self.check_module_linking()

    @lazy_property
    def id(self) -> str:
        name = self.name
        if self.name_suffix_set:
            name += '.' + self.suffix
        return self.construct_id_from_path(
            self.builddir, name, self.type_suffix())

    def _set_vala_args(self, kwargs: BuildTargetKeywordArguments) -> None:
        if self.uses_vala():
            self.vala_header = kwargs.get('vala_header') or self.name + '.h'
            inst = kwargs.get('install_vala_header', False)
            if inst is True:
                dir_ = kwargs.get('install_vala_header_dir') or self.environment.get_includedir()
                self.install_dir.append(dir_)
            else:
                self.install_dir.append(False)

            self.vala_vapi = kwargs.get('vala_vapi') or self.name + '.vapi'
            inst = kwargs.get('install_vala_vapi', False)
            if inst is True:
                dir_ = kwargs.get('install_vala_vapi_dir') or os.path.join(self.environment.get_datadir(), 'vala', 'vapi')
                self.install_dir.append(dir_)
            else:
                self.install_dir.append(False)

            self.vala_gir = kwargs.get('vala_gir')
            if self.vala_gir:
                inst = kwargs.get('install_vala_gir', False)
                if inst is True:
                    assert self.vala_gir is not None, 'Installing Vala GIR without generating one?'
                    dir_ = kwargs.get('install_vala_gir_dir') or os.path.join(self.environment.get_datadir(), 'gir-1.0')
                    self.install_dir.append(dir_)
                elif self.vala_gir is not None:
                    self.install_dir.append(False)

    def post_init(self) -> None:
        ''' Initialisations and checks requiring the final list of compilers to be known
        '''
        self.validate_sources()
        if self.uses_rust():
            if self.link_language and self.link_language != 'rust':
                raise MesonException('cannot build Rust sources with a different link_language')
            if self.structured_sources:
                # TODO: the interpreter should be able to generate a better error message?
                if any((s.endswith('.rs') for s in self.sources)) or \
                       any(any((s.endswith('.rs') for s in g.get_outputs())) for g in self.generated):
                    raise MesonException('cannot mix Rust structured sources and unstructured sources')

            # relocation-model=pic is rustc's default and Meson does not
            # currently have a way to disable PIC.
            self.pic = True
            self.pie = True
        else:
            if self.structured_sources:
                raise MesonException('structured sources are only supported in Rust targets')

        if self.is_linkable_target():
            if self.vala_header is not None:
                self.outputs.append(self.vala_header)
                self.install_tag.append('devel')
            if self.vala_vapi is not None:
                self.outputs.append(self.vala_vapi)
                self.install_tag.append('devel')
            if self.vala_gir is not None:
                self.outputs.append(self.vala_gir)
                self.install_tag.append('devel')

    def __repr__(self) -> str:
        repr_str = "<{0} {1}: {2}>"
        return repr_str.format(self.__class__.__name__, self.get_id(), self.filename)

    def __str__(self) -> str:
        return f"{self.name}"

    def validate_install(self) -> None:
        if self.for_machine is MachineChoice.BUILD and self.install:
            if self.environment.is_cross_build():
                raise InvalidArguments('Tried to install a target for the build machine in a cross build.')
            else:
                mlog.warning('Installing target build for the build machine. This will fail in a cross build.')

    def check_unknown_kwargs(self, kwargs: BuildTargetKeywordArguments) -> None:
        # Override this method in derived classes that have more
        # keywords.
        self.check_unknown_kwargs_int(kwargs, self.known_kwargs)

    def check_unknown_kwargs_int(self, kwargs: BuildTargetKeywordArguments, known_kwargs: T.Set[str]) -> None:
        unknowns = []
        for k in kwargs:
            if k in {'language_args', 'install_vala_header', 'install_vala_vapi',
                     'install_vala_gir', 'install_vala_header_dir',
                     'install_vala_vapi_dir', 'install_vala_gir_dir'}:
                continue
            if k not in known_kwargs:
                unknowns.append(k)
        if len(unknowns) > 0:
            mlog.warning('Unknown keyword argument(s) in target {}: {}.'.format(self.name, ', '.join(unknowns)))

    def process_objectlist(self, objects: T.List[ObjectTypes]) -> None:
        assert isinstance(objects, list)
        deprecated_non_objects = []
        for s in objects:
            if isinstance(s, (str, File, ExtractedObjects)):
                self.objects.append(s)
                if not isinstance(s, ExtractedObjects) and not is_object(s):
                    deprecated_non_objects.append(s)
            elif isinstance(s, (CustomTarget, CustomTargetIndex, GeneratedList)):
                non_objects = [o for o in s.get_outputs() if not is_object(o)]
                if non_objects:
                    raise InvalidArguments(f'Generated file {non_objects[0]} in the \'objects\' kwarg is not an object.')
                self.generated.append(s)
            else:
                raise InvalidArguments(f'Bad object of type {type(s).__name__!r} in target {self.name!r}.')
        if deprecated_non_objects:
            FeatureDeprecated.single_use(f'Source file {deprecated_non_objects[0]} in the \'objects\' kwarg is not an object.',
                                         '1.3.0', self.subproject)

    def process_sourcelist(self, sources: T.List['SourceOutputs']) -> None:
        """Split sources into generated and static sources.

        Sources can be:
        1. Preexisting source files in the source tree (static)
        2. Preexisting sources generated by configure_file in the build tree.
           (static as they are only regenerated if meson itself is regenerated)
        3. Sources files generated by another target or a Generator (generated)
        """
        for s in sources:
            if isinstance(s, File):
                if s not in self.seen_sources:
                    self.sources.append(s)
                    self.seen_sources.add(s)
            elif isinstance(s, (CustomTarget, CustomTargetIndex, GeneratedList)):
                self.generated.append(s)

    @staticmethod
    def can_compile_remove_sources(compiler: 'Compiler', sources: T.List['FileOrString']) -> bool:
        removed = False
        for s in sources[:]:
            if compiler.can_compile(s):
                sources.remove(s)
                removed = True
        return removed

    def process_compilers_late(self) -> None:
        """Processes additional compilers after kwargs have been evaluated.

        This can add extra compilers that might be required by keyword
        arguments, such as link_with or dependencies. It will also try to guess
        which compiler to use if one hasn't been selected already.
        """
        for lang in self.missing_languages:
            try:
                self.compilers[lang] = self.all_compilers[lang]
            except KeyError:
                # It is intentional not to update self.all_compilers here, since
                # that is a reference to Interpreter.compilers, which should
                # only contain the compilers enabled explicilty by the project

                is_error: bool = True

                # In the case of cython it's possible that we have an
                # implementation detail language
                if self.uses_cython() and lang == self.environment.coredata.get_option_for_target(self, 'cython_language'):
                    self.compilers[lang] = self.environment.coredata.compilers[self.for_machine][lang]
                    is_error = False

                # C++ is an implementation detail of Cuda, and may not be
                # explicitly enabled.
                if self.uses_cuda() and lang == 'cpp':
                    self.compilers[lang] = self.environment.coredata.compilers[self.for_machine][lang]
                    is_error = False

                if is_error:
                    raise

        # did user override clink_langs for this target?
        link_langs = (self.link_language, ) if self.link_language else clink_langs

        if self.link_language:
            if self.link_language not in self.all_compilers:
                m = f'Target {self.name} requires {self.link_language} compiler not part of the project'
                raise MesonException(m)

        # If this library is linked against another library we need to consider
        # the languages of those libraries as well.
        if self.link_targets or self.link_whole_targets:
            for t in itertools.chain(self.link_targets, self.link_whole_targets):
                if isinstance(t, (CustomTarget, CustomTargetIndex)):
                    continue # We can't know anything about these.
                target_langs = [t.link_language] if t.link_language else list(t.compilers)
                for lang in target_langs:
                    # Rust is always linked through a C-ABI target, so do not add
                    # the compiler here
                    if lang != 'rust' and lang in link_langs:
                        if lang not in self.compilers:
                            try:
                                self.compilers[lang] = t.all_compilers[lang]
                            except KeyError:
                                # The language needed might be provided by a
                                # subproject, and will not be in this project's
                                # compiler set
                                self.compilers[lang] = self.environment.coredata.compilers[self.for_machine][lang]

        if not self.compilers:
            # No source files or parent targets, target consists of only object
            # files of unknown origin. Just add the first clink compiler
            # that we have and hope that it can link these objects

            # TODO: pyright understands this, mypy doesn't
            for lang in T.cast('T.Iterable[Language]', reversed(link_langs)):
                if lang in self.all_compilers:
                    self.compilers[lang] = self.all_compilers[lang]
                    break

        # Now that we have the final list of compilers we can sort it according
        # to clink_langs and do sanity checks.
        self.compilers = OrderedDict(sorted(self.compilers.items(),
                                            key=lambda t: sort_clink(t[0])))
        self.post_init()

    def process_compilers(self) -> T.List[Language]:
        '''
        Populate self.compilers, which is the list of compilers that this
        target will use for compiling all its sources.
        We also add compilers that were used by extracted objects to simplify
        dynamic linker determination.
        Returns a list of missing languages that we can add implicitly, such as
        C/C++ compiler for cython.
        '''
        missing_languages: T.List[Language] = []
        if not any([self.sources, self.generated, self.objects, self.structured_sources]):
            return missing_languages
        # Preexisting sources
        sources: T.List['FileOrString'] = list(self.sources)
        generated = self.generated.copy()

        if self.structured_sources:
            for v in self.structured_sources.sources.values():
                for src in v:
                    if isinstance(src, (str, File)):
                        sources.append(src)
                    else:
                        generated.append(src)

        # All generated sources
        for gensrc in generated:
            for s in gensrc.get_outputs():
                # Generated objects can't be compiled, so don't use them for
                # compiler detection. If our target only has generated objects,
                # we will fall back to using the first c-like compiler we find,
                # which is what we need.
                if not is_object(s):
                    sources.append(s)
        for d in self.external_deps:
            for s in d.sources:
                if isinstance(s, (str, File)):
                    sources.append(s)

        # Sources that were used to create our extracted objects
        for o in self.objects:
            if not isinstance(o, ExtractedObjects):
                continue
            compsrcs = o.classify_all_sources(o.srclist, [])
            for comp in compsrcs:
                # Don't add Vala sources since that will pull in the Vala
                # compiler even though we will never use it since we are
                # dealing with compiled C code.
                if comp.language == 'vala':
                    continue
                if comp.language not in self.compilers:
                    self.compilers[comp.language] = comp
        if sources:
            # For each source, try to add one compiler that can compile it.
            #
            # If it has a suffix that belongs to a known language, we must have
            # a compiler for that language.
            #
            # Otherwise, it's ok if no compilers can compile it, because users
            # are expected to be able to add arbitrary non-source files to the
            # sources list
            for s in sources:
                for lang, compiler in self.all_compilers.items():
                    if compiler.can_compile(s):
                        if lang not in self.compilers:
                            self.compilers[lang] = compiler
                        break
                else:
                    if is_known_suffix(s) and not is_header(s):
                        path = pathlib.Path(str(s)).as_posix()
                        m = f'No {self.for_machine.get_lower_case_name()} machine compiler for {path!r}'
                        raise MesonException(m)

        # If all our sources are Vala, our target also needs the C compiler but
        # it won't get added above.
        if 'vala' in self.compilers and 'c' not in self.compilers:
            self.compilers['c'] = self.all_compilers['c']
        if 'cython' in self.compilers:
            _value = self.environment.coredata.get_option_for_target(self, 'cython_language')
            assert isinstance(_value, str), 'for mypy'
            value = T.cast('Language', _value)
            try:
                self.compilers[value] = self.all_compilers[value]
            except KeyError:
                # This might be an implementation detail (C or C++)
                #
                # It is intentional not to update self.all_compilers here, since
                # that is a reference to Interpreter.compilers, which should
                # only contain the compilers enabled explicilty by the project
                try:
                    self.compilers[value] = self.environment.coredata.compilers[self.for_machine][value]
                except KeyError:
                    missing_languages.append(value)

        return missing_languages

    def validate_sources(self) -> None:
        if len(self.compilers) > 1 and any(lang in self.compilers for lang in ['cs', 'java']):
            langs = ', '.join(self.compilers.keys())
            raise InvalidArguments(f'Cannot mix those languages into a target: {langs}')

    def process_link_depends(self, sources: T.Iterable[T.Union[str, File, BuildTargetTypes]]) -> None:
        """Process the link_depends keyword argument.

        This is designed to handle strings, Files, and the output of Custom
        Targets. Notably it doesn't handle generator() returned objects, since
        adding them as a link depends would inherently cause them to be
        generated twice, since the output needs to be passed to the ld_args and
        link_depends.
        """
        for s in sources:
            if isinstance(s, File):
                self.link_depends.append(s)
            elif isinstance(s, str):
                self.link_depends.append(
                    File.from_source_file(self.environment.source_dir, self.get_subdir(), s))
            else:
                self.link_depends.append(s)

    def extract_objects(self, srclist: T.List[T.Union['FileOrString', 'GeneratedTypes']], is_unity: bool) -> ExtractedObjects:
        sources_set = set(self.sources)
        generated_set = set(self.generated)

        obj_src: T.List['File'] = []
        obj_gen: T.List['GeneratedTypes'] = []
        for src in srclist:
            if isinstance(src, (str, File)):
                if isinstance(src, str):
                    src = File(False, self.subdir, src)
                else:
                    FeatureNew.single_use('File argument for extract_objects', '0.50.0', self.subproject)
                if src not in sources_set:
                    raise MesonException(f'Tried to extract unknown source {src}.')
                obj_src.append(src)
            elif isinstance(src, (CustomTarget, CustomTargetIndex, GeneratedList)):
                FeatureNew.single_use('Generated sources for extract_objects', '0.61.0', self.subproject)
                target = src.target if isinstance(src, CustomTargetIndex) else src
                if src not in generated_set and target not in generated_set:
                    raise MesonException(f'Tried to extract unknown source {target.get_basename()}.')
                obj_gen.append(src)
            else:
                raise MesonException(f'Object extraction arguments must be strings, Files or targets (got {type(src).__name__}).')
        eobjs = ExtractedObjects(self, obj_src, obj_gen)
        if is_unity:
            eobjs.check_unity_compatible()
        return eobjs

    def extract_all_objects(self, recursive: bool = True) -> ExtractedObjects:
        return ExtractedObjects(self, self.sources, self.generated, self.objects,
                                recursive, pch=True)

    @lru_cache(maxsize=None)
    def get_all_link_deps(self) -> ImmutableListProtocol[BuildTargetTypes]:
        """ Get all shared libraries dependencies
        This returns all shared libraries in the entire dependency tree. Those
        are libraries needed at runtime which is different from the set needed
        at link time, see get_dependencies() for that.
        """
        result: OrderedSet[BuildTargetTypes] = OrderedSet()
        nonresults: T.Set[BuildTargetTypes] = set()
        stack: T.Deque[BuildTargetTypes] = deque()
        stack.appendleft(self)
        while stack:
            t = stack.pop()
            if t in result or t in nonresults:
                continue
            if isinstance(t, CustomTargetIndex):
                stack.appendleft(t.target)
                continue
            if isinstance(t, SharedLibrary):
                result.add(t)
            else:
                nonresults.add(t)
            if isinstance(t, BuildTarget):
                stack.extendleft((t2 for t2 in t.link_targets if t2 not in nonresults))
                stack.extendleft((t2 for t2 in t.link_whole_targets if t2 not in nonresults))
        return list(result)

    @lru_cache(maxsize=None)
    def get_all_linked_targets(self) -> ImmutableListProtocol[BuildTargetTypes]:
        """Get all targets that have been linked with this one.

        This is useful for cases where we need to analyze these links, such as
        for module information.

        This includes static libraries and static libraries linked with static
        libraries. This differs from :method:`get_all_link_deps` in that it does
        add static libs, and differs from `:method:`get_dependencies`, which
        does not look for targets that are not directly linked, such as those
        that are added with `link_whole`.

        :returns: An immutable list of BuildTargets
        """
        result: OrderedSet[BuildTargetTypes] = OrderedSet()
        stack: T.Deque[BuildTargetTypes] = deque()
        stack.extendleft(self.link_targets)
        stack.extendleft(self.link_whole_targets)
        while stack:
            t = stack.pop()
            if t in result:
                continue
            if isinstance(t, CustomTargetIndex):
                stack.appendleft(t.target)
                continue
            if isinstance(t, BuildTarget):
                result.add(t)
                stack.extendleft(t.link_targets)
                stack.extendleft(t.link_whole_targets)
        assert self not in result, 'should not have self'
        return list(result)

    def get_link_deps_mapping(self, prefix: str) -> T.Mapping[str, str]:
        return self.get_transitive_link_deps_mapping(prefix)

    @lru_cache(maxsize=None)
    def get_transitive_link_deps_mapping(self, prefix: str) -> T.Mapping[str, str]:
        result: T.Dict[str, str] = {}
        for i in self.link_targets:
            mapping = i.get_link_deps_mapping(prefix)
            #we are merging two dictionaries, while keeping the earlier one dominant
            result_tmp = mapping.copy()
            result_tmp.update(result)
            result = result_tmp
        return result

    @lru_cache(maxsize=None)
    def get_link_dep_subdirs(self) -> T.AbstractSet[str]:
        result: OrderedSet[str] = OrderedSet()
        for i in self.link_targets:
            if not isinstance(i, StaticLibrary):
                result.add(i.get_builddir())
            result.update(i.get_link_dep_subdirs())
        return result

    def get_default_install_dir(self) -> T.Tuple[str, str]:
        return self.environment.get_libdir(), '{libdir}'

    def get_custom_install_dir(self) -> T.List[T.Union[str, Literal[False]]]:
        return self.install_dir

    def get_custom_install_mode(self) -> T.Optional['FileMode']:
        return self.install_mode

    def process_kwargs(self, kwargs: BuildTargetKeywordArguments) -> None:
        self.original_kwargs = kwargs

        if 'build_by_default' in kwargs:
            self.build_by_default = kwargs['build_by_default']

        if not self.build_by_default and kwargs.get('install', False):
            # For backward compatibility, if build_by_default is not explicitly
            # set, use the value of 'install' if it's enabled.
            self.build_by_default = True

        self.raw_overrides = kwargs.get('override_options', {})

        self.pch['c'] = kwargs.get('c_pch')
        self.pch['cpp'] = kwargs.get('cpp_pch')

        self.link_args = kwargs.get('link_args', [])
        for l in self.link_args:
            if '-Wl,-rpath' in l or l.startswith('-rpath'):
                mlog.warning(textwrap.dedent('''\
                    Please do not define rpath with a linker argument, use install_rpath
                    or build_rpath properties instead.
                    This will become a hard error in a future Meson release.
                '''))
        self.link_early_args = kwargs.get('link_early_args', [])
        for i in self.link_early_args:
            if not isinstance(i, str):
                raise InvalidArguments('link_early_args values must be strings.')
        self.process_link_depends(kwargs.get('link_depends', []))
        # Target-specific include dirs must be added BEFORE include dirs from
        # internal deps (added inside self.add_deps()) to override them.
        self.add_include_dirs(kwargs.get('include_directories', []))
        # Add dependencies (which also have include_directories)
        self.add_deps(kwargs.get('dependencies', []))

        self.has_custom_install_dir = False
        i = kwargs.get('install_dir', [])
        install_dir = i[0] if i else True
        default_install_dir = self.get_default_install_dir()[0]
        if install_dir is True:
            install_dir = default_install_dir
        elif install_dir != default_install_dir:
            self.has_custom_install_dir = True
        self.install_dir: T.List[T.Union[str, T.Literal[False]]] = [install_dir]

        self.install_mode = kwargs.get('install_mode', None)
        self.install_tag: T.List[T.Optional[str]] = kwargs.get('install_tag') or [None]
        self.extra_files = kwargs.get('extra_files', [])
        self.install_rpath: str = kwargs.get('install_rpath', '')
        self.build_rpath = kwargs.get('build_rpath', '')
        self.resources = kwargs.get('resources', [])
        name_prefix = kwargs.get('name_prefix')
        if name_prefix is not None:
            self.prefix = name_prefix
            self.name_prefix_set = True
        name_suffix = kwargs.get('name_suffix')
        if name_suffix is not None:
            self.suffix = name_suffix
            self.name_suffix_set = True
        self.implicit_include_directories = kwargs.get('implicit_include_directories', True)
        self.gnu_symbol_visibility = kwargs.get('gnu_symbol_visibility', '')
        self.rust_dependency_map = kwargs.get('rust_dependency_map', {})

        self.swift_interoperability_mode = kwargs.get('swift_interoperability_mode', 'c')
        self.swift_module_name = kwargs.get('swift_module_name') or self.name

    @T.overload
    def _extract_pic_pie(self, kwargs: StaticLibraryKeywordArguments, arg: Literal['pic'],
                         option: Literal['b_staticpic']) -> bool: ...

    @T.overload
    def _extract_pic_pie(self, kwargs: T.Union[StaticLibraryKeywordArguments, ExecutableKeywordArguments],
                         arg: Literal['pie'], option: Literal['b_pie']) -> bool: ...

    def _extract_pic_pie(self, kwargs: T.Union[StaticLibraryKeywordArguments, ExecutableKeywordArguments],
                         arg: Literal['pic', 'pie'], option: Literal['b_staticpic', 'b_pie']) -> bool:
        # You can't disable PIC on OS X. The compiler ignores -fno-PIC.
        # PIC is always on for Windows (all code is position-independent
        # since library loading is done differently)
        m = self.environment.machines[self.for_machine]
        assert m is not None, 'for mypy'
        if arg == 'pic' and (m.is_darwin() or m.is_windows()):
            return True

        # Executables must be PIE on Android
        if arg == 'pie' and m.is_android():
            return True

        # Check if we have -fPIC, -fpic, -fPIE, or -fpie in cflags
        all_flags = self.extra_args['c'] + self.extra_args['cpp']
        if '-f' + arg.lower() in all_flags or '-f' + arg.upper() in all_flags:
            mlog.warning(f"Use the '{arg}' kwarg instead of passing '-f{arg}' manually to {self.name!r}")
            return True

        k = OptionKey(option)
        if kwargs.get(arg) is not None:
            return kwargs[arg]
        elif k in self.environment.coredata.optstore:
            val = self.environment.coredata.get_option_for_target(self, k)
            assert isinstance(val, bool), 'for mypy'
            return val
        return False

    def install_dir_names(self) -> T.List[T.Optional[str]]:
        install_dir_names: T.List[T.Optional[str]]
        if self.has_custom_install_dir:
            install_dir_names = [getattr(i, 'optname', None) for i in self.install_dir]
        else:
            default = self.get_default_install_dir()[1]
            install_dir_names = T.cast('T.List[T.Optional[str]]', [default]) * len(self.install_dir)

        return install_dir_names

    def get_filename(self) -> str:
        return self.filename

    def get_debug_filename(self) -> T.Optional[str]:
        """
        The name of debuginfo file that will be created by the compiler

        Returns None if the build won't create any debuginfo file
        """
        return self.debug_filename

    def get_outputs(self) -> T.List[str]:
        return self.outputs

    def get_extra_args(self, language: Language) -> T.List[str]:
        return self.extra_args[language]

    @lru_cache(maxsize=None)
    def get_dependencies(self) -> OrderedSet[BuildTargetTypes]:
        # Get all targets needed for linking. This includes all link_with and
        # link_whole targets, and also all dependencies of static libraries
        # recursively. The algorithm here is closely related to what we do in
        # get_internal_static_libraries(): Installed static libraries include
        # objects from all their dependencies already.
        result: OrderedSet[BuildTargetTypes] = OrderedSet()
        visited: T.Set[T.Tuple[BuildTargetTypes, bool, bool]] = set()
        for t in itertools.chain(self.link_targets, self.link_whole_targets):
            if t not in result:
                result.add(t)
                if isinstance(t, StaticLibrary):
                    t.get_dependencies_recurse(result, visited, handled_by_rustc=self.uses_rust())
        return result

    def get_dependencies_recurse(self, result: OrderedSet[BuildTargetTypes],
                                 visited: T.Set[T.Tuple[BuildTargetTypes, bool, bool]],
                                 include_internals: bool = True, handled_by_rustc: bool = False) -> None:
        # self is always a static library because we don't need to pull dependencies
        # of shared libraries. If self is installed (not internal) it already
        # include objects extracted from all its internal dependencies so we can
        # skip them.
        include_internals = include_internals and self.is_internal()
        key = (self, include_internals, handled_by_rustc)
        if key in visited:
            return
        visited.add(key)

        for t in self.link_targets:
            uses_rust_abi = isinstance(t, BuildTarget) and t.uses_rust_abi()
            if not handled_by_rustc and uses_rust_abi:
                # Rules for including libraries via Rust rlibs and staticlibs are complex:
                # - proc-macro crates should be skipped completely when the build product
                #   is not a Rust program, because only rustc knows that they are
                #   special build-machine shared libraries
                # - rlibs must always be returned for Rust programs, because even though
                #   the -l flag is implicitly added, the -L flag is not.  ninjabackend.py
                #   handles leaving out the -l flag
                # - rlibs are bundled into staticlibs and need not be in the command line of
                #   non-Rust programs; this is the case when t is not added to result.
                # - C-ABI libraries included in link_with use -lstatic:-bundle, so even for
                #   staticlibs we do need to recurse into rlibs and collect these non-bundled
                #   libraries.  So don't return, unlike for procedural macros
                if t.rust_crate_type == 'proc-macro':
                    continue

            elif not include_internals and t.is_internal():
                pass

            else:
                result.add(t)
            if isinstance(t, StaticLibrary):
                t.get_dependencies_recurse(result, visited, include_internals, handled_by_rustc and uses_rust_abi)
        for t in self.link_whole_targets:
            uses_rust_abi = isinstance(t, BuildTarget) and t.uses_rust_abi()
            t.get_dependencies_recurse(result, visited, include_internals, handled_by_rustc and uses_rust_abi)

    def get_sources(self) -> T.List[File]:
        return self.sources

    def get_objects(self) -> T.List[T.Union[str, 'File', 'ExtractedObjects']]:
        return self.objects

    def get_generated_sources(self) -> T.List['GeneratedTypes']:
        return self.generated

    def should_install(self) -> bool:
        return self.install

    def has_pch(self) -> bool:
        return any(x is not None for x in self.pch.values())

    def get_include_dirs(self) -> T.List['IncludeDirs']:
        return self.include_dirs

    def add_deps(self, deps: T.List[dependencies.Dependency]) -> None:
        for dep in deps:
            if dep in self.added_deps:
                # Prefer to add dependencies to added_deps which have a name
                if dep.is_named():
                    self.added_deps.remove(dep)
                    self.added_deps.add(dep)
                continue

            if isinstance(dep, dependencies.InternalDependency):
                # Those parts that are internal.
                self.process_sourcelist(dep.sources)
                self.extra_files.extend(f for f in dep.extra_files if f not in self.extra_files)
                self.add_include_dirs(dep.include_directories, dep.get_include_type())
                self.objects.extend(dep.objects)
                self.link_targets.extend(dep.libraries)
                self.link_whole_targets.extend(dep.whole_libraries)
                if dep.get_compile_args() or dep.get_link_args():
                    # Those parts that are external.
                    extpart = type(dep)(dep.version,
                                        compile_args=dep.get_compile_args(),
                                        link_args=dep.get_link_args(),
                                        name=dep.name)
                    self.external_deps.append(extpart)
                # Deps of deps.
                self.add_deps(dep.ext_deps)
            else:
                if dep not in self.external_deps:
                    self.external_deps.append(dep)
                    self.process_sourcelist(dep.get_sources())
                self.add_deps(dep.ext_deps)

            dep_d_features = dep.d_features

            for feature in ('versions', 'import_dirs'):
                if feature in dep_d_features:
                    self.d_features[feature].extend(dep_d_features[feature])

            self.added_deps.add(dep)

    def get_external_deps(self) -> T.List[dependencies.Dependency]:
        return self.external_deps

    def is_internal(self) -> bool:
        return False

    def link(self, targets: T.List[BuildTargetTypes]) -> None:
        for t in targets:
            self.check_can_link_together(t)
            self.link_targets.append(t)

    def link_whole(
            self,
            targets: T.List[StaticTargetTypes],
            promoted: bool = False) -> None:
        for t in targets:
            self.check_can_link_together(t)
            self.link_whole_targets.append(t)

    @lru_cache(maxsize=None)
    def get_internal_static_libraries(self) -> OrderedSet[StaticTargetTypes]:
        result: OrderedSet[StaticTargetTypes] = OrderedSet()
        self.get_internal_static_libraries_recurse(result)
        return result

    def get_internal_static_libraries_recurse(self, result: OrderedSet[StaticTargetTypes]) -> None:
        for t in self.link_targets:
            if t.is_internal() and t not in result:
                result.add(t)
                t.get_internal_static_libraries_recurse(result)
        for t in self.link_whole_targets:
            if t.is_internal():
                t.get_internal_static_libraries_recurse(result)

    def check_can_link_together(self, t: BuildTargetTypes) -> None:
        links_with_rust_abi = isinstance(t, BuildTarget) and t.uses_rust_abi()
        if not self.uses_rust() and links_with_rust_abi:
            raise InvalidArguments(f'Try to link Rust ABI library {t.name!r} with a non-Rust target {self.name!r}')
        if self.for_machine is not t.for_machine and (not links_with_rust_abi or t.rust_crate_type != 'proc-macro'):
            msg = f'Tried to mix a {t.for_machine} library ("{t.name}") with a {self.for_machine} target "{self.name}"'
            if self.environment.is_cross_build():
                raise InvalidArguments(msg + ' This is not possible in a cross build.')
            else:
                mlog.warning(msg + ' This will fail in cross build.')

    def add_include_dirs(self, args: T.Sequence['IncludeDirs'], set_is_system: str = 'preserve') -> None:
        if set_is_system != 'preserve':
            is_system = set_is_system == 'system'
            self.include_dirs.extend([IncludeDirs(x.curdir, x.incdirs, is_system, x.extra_build_dirs) for x in args])
        else:
            self.include_dirs.extend(args)

    def get_aliases(self) -> T.List[T.Tuple[str, str, str]]:
        return []

    def get_langs_used_by_deps(self) -> T.List[Language]:
        '''
        Sometimes you want to link to a C++ library that exports C API, which
        means the linker must link in the C++ stdlib, and we must use a C++
        compiler for linking. The same is also applicable for objc/objc++, etc,
        so we can keep using clink_langs for the priority order.

        See: https://github.com/mesonbuild/meson/issues/1653
        '''
        langs: T.List[Language] = []

        # Check if any of the external libraries were written in this language
        for dep in self.external_deps:
            if dep.language is None:
                continue
            if dep.language not in langs:
                langs.append(dep.language)
        # Check if any of the internal libraries this target links to were
        # written in this language
        for t in itertools.chain(self.link_targets, self.link_whole_targets):
            if isinstance(t, (CustomTarget, CustomTargetIndex)):
                continue
            target_langs = [t.link_language] if t.link_language else list(t.compilers)
            for language in target_langs:
                if language == 'rust' and not t.uses_rust_abi():
                    # All Rust dependencies must go through a C-ABI dependency, so ignore it
                    continue
                if language not in langs:
                    langs.append(language)

        return langs

    def get_prelinker(self) -> Compiler:
        if self.link_language:
            comp = self.all_compilers[self.link_language]
            return comp
        for l in T.cast('T.Tuple[Language, ...]', clink_langs):
            if l in self.compilers:
                try:
                    prelinker = self.all_compilers[l]
                except KeyError:
                    raise MesonException(
                        f'Could not get a prelinker linker for build target {self.name!r}. '
                        f'Requires a compiler for language "{l}", but that is not '
                        'a project language.')
                return prelinker
        raise MesonException(f'Could not determine prelinker for {self.name!r}.')

    def get_clink_dynamic_linker_and_stdlibs(self) -> T.Tuple['Compiler', T.List[str]]:
        '''
        We use the order of languages in `clink_langs` to determine which
        linker to use in case the target has sources compiled with multiple
        compilers. All languages other than those in this list have their own
        linker.
        Note that Vala outputs C code, so Vala sources can use any linker
        that can link compiled C. We don't actually need to add an exception
        for Vala here because of that.
        '''
        # If the user set the link_language, just return that.
        if self.link_language:
            comp = self.all_compilers[self.link_language]
            return comp, comp.language_stdlib_only_link_flags()

        # Since dependencies could come from subprojects, they could have
        # languages we don't have in self.all_compilers. Use the global list of
        # all compilers here.
        all_compilers = self.environment.coredata.compilers[self.for_machine]

        # Languages used by dependencies
        dep_langs = self.get_langs_used_by_deps()

        # Pick a compiler based on the language priority-order
        for l in T.cast('T.Tuple[Language, ...]', clink_langs):
            if l in self.compilers or l in dep_langs:
                try:
                    linker = all_compilers[l]
                except KeyError:
                    raise MesonException(
                        f'Could not get a dynamic linker for build target {self.name!r}. '
                        f'Requires a linker for language "{l}", but that is not '
                        'a project language.')
                stdlib_args: T.List[str] = self.get_used_stdlib_args(linker.language)
                # Type of var 'linker' is Compiler.
                # Pretty hard to fix because the return value is passed everywhere
                return linker, stdlib_args

        # None of our compilers can do clink, this happens for example if the
        # target only has ASM sources. Pick the first capable compiler.
        for l in T.cast('T.Tuple[Language, ...]', clink_langs):
            try:
                comp = self.all_compilers[l]
                return comp, comp.language_stdlib_only_link_flags()
            except KeyError:
                pass

        raise AssertionError(f'Could not get a dynamic linker for build target {self.name!r}')

    def get_used_stdlib_args(self, link_language: Language) -> T.List[str]:
        all_compilers = self.environment.coredata.compilers[self.for_machine]
        all_langs = set(self.compilers).union(self.get_langs_used_by_deps())
        stdlib_args: T.List[str] = []
        for dl in all_langs:
            if dl != link_language and (dl, link_language) not in self._MASK_LANGS:
                # We need to use all_compilers here because
                # get_langs_used_by_deps could return a language from a
                # subproject
                stdlib_args.extend(all_compilers[dl].language_stdlib_only_link_flags())
        return stdlib_args

    def uses_rust(self) -> bool:
        return 'rust' in self.compilers

    def uses_rust_abi(self) -> bool:
        return self.uses_rust() and self.rust_crate_type in {'dylib', 'rlib', 'proc-macro'}

    def uses_fortran(self) -> bool:
        return 'fortran' in self.compilers

    def uses_cython(self) -> bool:
        return 'cython' in self.compilers

    def uses_cuda(self) -> bool:
        return 'cuda' in self.compilers

    def uses_vala(self) -> bool:
        return 'vala' in self.compilers

    def uses_swift_cpp_interop(self) -> bool:
        return self.swift_interoperability_mode == 'cpp' and 'swift' in self.compilers

    def get_using_msvc(self) -> bool:
        '''
        Check if the dynamic linker is MSVC. Used by Executable, StaticLibrary,
        and SharedLibrary for deciding when to use MSVC-specific file naming
        and debug filenames.

        If at least some code is built with MSVC and the final library is
        linked with MSVC, we can be sure that some debug info will be
        generated. We only check the dynamic linker here because the static
        linker is guaranteed to be of the same type.

        Interesting cases:
        1. The Vala compiler outputs C code to be compiled by whatever
           C compiler we're using, so all objects will still be created by the
           MSVC compiler.
        2. If the target contains only objects, process_compilers guesses and
           picks the first compiler that smells right.
        '''
        # Rustc can use msvc style linkers
        if self.uses_rust():
            compiler = self.all_compilers['rust']
        else:
            compiler, _ = self.get_clink_dynamic_linker_and_stdlibs()
        # Mixing many languages with MSVC is not supported yet so ignore stdlibs.
        return compiler and compiler.get_linker_id() in {'link', 'lld-link', 'xilink', 'optlink'}

    def check_module_linking(self) -> None:
        '''
        Warn if shared modules are linked with target: (link_with) #2865
        '''
        for link_target in self.link_targets:
            if isinstance(link_target, SharedModule) and not link_target.force_soname:
                if self.environment.machines[self.for_machine].is_darwin():
                    raise MesonException(
                        f'target {self.name} links against shared module {link_target.name}. This is not permitted on OSX')
                elif self.environment.machines[self.for_machine].is_android() and isinstance(self, SharedModule):
                    # Android requires shared modules that use symbols from other shared modules to
                    # be linked before they can be dlopen()ed in the correct order. Not doing so
                    # leads to a missing symbol error: https://github.com/android/ndk/issues/201
                    link_target.force_soname = True
                else:
                    mlog.deprecation(f'target {self.name} links against shared module {link_target.name}, which is incorrect.'
                                     '\n             '
                                     f'This will be an error in meson 2.0, so please use shared_library() for {link_target.name} instead.'
                                     '\n             '
                                     f'If shared_module() was used for {link_target.name} because it has references to undefined symbols,'
                                     '\n             '
                                     'use shared_library() with `override_options: [\'b_lundef=false\']` instead.')
                    link_target.force_soname = True

    def process_vs_module_defs_kw(self, kwargs: ExecutableKeywordArguments) -> None:
        path = kwargs.get('vs_module_defs')
        if path is None:
            return

        if isinstance(path, str):
            if os.path.isabs(path):
                self.vs_module_defs = File.from_absolute_file(path)
            else:
                self.vs_module_defs = File.from_source_file(self.environment.source_dir, self.subdir, path)
        elif isinstance(path, File):
            # When passing a generated file.
            self.vs_module_defs = path
        else:
            # When passing output of a Custom Target
            self.vs_module_defs = File.from_built_file(path.get_builddir(), path.get_filename())
        self.process_link_depends([path])

    def extract_targets_as_list(self, kwargs: BuildTargetKeywordArguments, key: T.Literal['link_with', 'link_whole']) -> T.List[LibTypes]:
        bl_type = self.environment.coredata.optstore.get_value_for(OptionKey('default_both_libraries'))
        if bl_type == 'auto':
            if isinstance(self, StaticLibrary):
                bl_type = 'static'
            elif isinstance(self, SharedLibrary):
                bl_type = 'shared'

        self_libs: T.List[LibTypes] = self.link_targets if key == 'link_with' else self.link_whole_targets

        lib_list = []
        for lib in listify(kwargs.get(key, [])) + self_libs:
            if isinstance(lib, (Target, BothLibraries)):
                lib_list.append(lib.get(bl_type))
            else:
                lib_list.append(lib)
        return lib_list

    def get(self, lib_type: T.Literal['static', 'shared']) -> LibTypes:
        """Base case used by BothLibraries"""
        return self

    def determine_rpath_dirs(self) -> T.Tuple[str, ...]:
        result: OrderedSet[str]
        if self.environment.coredata.optstore.get_value_for(OptionKey('layout')) == 'mirror':
            # Need a copy here
            result = OrderedSet(self.get_link_dep_subdirs())
        else:
            result = OrderedSet()
            result.add('meson-out')
        result.update(self.rpaths_for_non_system_absolute_shared_libraries())
        self.rpath_dirs_to_remove.update([d.encode('utf-8') for d in result])
        return tuple(result)

    @lru_cache(maxsize=None)
    def rpaths_for_non_system_absolute_shared_libraries(self, exclude_system: bool = True) -> ImmutableListProtocol[str]:
        paths: OrderedSet[str] = OrderedSet()
        srcdir = self.environment.get_source_dir()

        system_dirs = set()
        if exclude_system:
            for cc in self.compilers.values():
                system_dirs.update(cc.get_library_dirs())

        external_rpaths = self.get_external_rpath_dirs()
        build_to_src = relpath(self.environment.get_source_dir(),
                               self.environment.get_build_dir())

        for dep in self.external_deps:
            if dep.type_name not in {'library', 'pkgconfig', 'cmake'}:
                continue
            for libpath in dep.link_args:
                if libpath.startswith('-'):
                    continue
                # For all link args that are absolute paths to a library file, add RPATH args
                if not os.path.isabs(libpath):
                    continue
                libdir, libname = os.path.split(libpath)
                # Windows doesn't support rpaths, but we use this function to
                # emulate rpaths by setting PATH
                # .dll is there for mingw gcc
                # .so's may be extended with version information, e.g. libxyz.so.1.2.3
                if not (
                    libname.endswith(('.dll', '.lib', '.so', '.dylib'))
                    or '.so.' in libname
                ):
                    continue

                # Don't remove rpaths specified in LDFLAGS.
                if libdir in external_rpaths:
                    continue
                if system_dirs and os.path.normpath(libdir) in system_dirs:
                    # No point in adding system paths.
                    continue

                if is_parent_path(srcdir, libdir):
                    rel_to_src = libdir[len(srcdir) + 1:]
                    assert not os.path.isabs(rel_to_src), f'rel_to_src: {rel_to_src} is absolute'
                    paths.add(os.path.join(build_to_src, rel_to_src))
                else:
                    paths.add(libdir)
            # Don't remove rpaths specified by the dependency
            paths.difference_update(self.get_rpath_dirs_from_link_args(dep.link_args))
        for i in itertools.chain(self.link_targets, self.link_whole_targets):
            if isinstance(i, BuildTarget):
                paths.update(i.rpaths_for_non_system_absolute_shared_libraries(exclude_system))
        return list(paths)

    def get_external_rpath_dirs(self) -> T.Set[str]:
        args: T.List[str] = []
        for lang in LANGUAGES_USING_LDFLAGS:
            try:
                args += self.environment.coredata.get_external_link_args(self.for_machine, lang)
            except KeyError:
                pass
        return self.get_rpath_dirs_from_link_args(args)

    # Match rpath formats:
    # -Wl,-rpath=
    # -Wl,-rpath,
    _rpath_regex = re.compile(r'-Wl,-rpath[=,]([^,]+)')
    # Match solaris style compat runpath formats:
    # -Wl,-R
    # -Wl,-R,
    _runpath_regex = re.compile(r'-Wl,-R[,]?([^,]+)')
    # Match symbols formats:
    # -Wl,--just-symbols=
    # -Wl,--just-symbols,
    _symbols_regex = re.compile(r'-Wl,--just-symbols[=,]([^,]+)')

    @classmethod
    def get_rpath_dirs_from_link_args(cls, args: T.List[str]) -> T.Set[str]:
        dirs: T.Set[str] = set()

        for arg in args:
            if not arg.startswith('-Wl,'):
                continue

            rpath_match = cls._rpath_regex.match(arg)
            if rpath_match:
                for dir in rpath_match.group(1).split(':'):
                    dirs.add(dir)
            runpath_match = cls._runpath_regex.match(arg)
            if runpath_match:
                for dir in runpath_match.group(1).split(':'):
                    # The symbols arg is an rpath if the path is a directory
                    if os.path.isdir(dir):
                        dirs.add(dir)
            symbols_match = cls._symbols_regex.match(arg)
            if symbols_match:
                for dir in symbols_match.group(1).split(':'):
                    # Prevent usage of --just-symbols to specify rpath
                    if os.path.isdir(dir):
                        raise MesonException(f'Invalid arg for --just-symbols, {dir} is a directory.')
        return dirs

    def get_platform_scheme_name(self) -> str:
        m = self.environment.machines[self.for_machine]
        if m.is_cygwin():
            return 'cygwin'
        elif m.is_windows():
            return 'windows'
        elif m.is_darwin():
            return 'darwin'
        else:
            return 'unix'


class FileInTargetPrivateDir:
    """Represents a file with the path '/path/to/build/target_private_dir/fname'.
       target_private_dir is the return value of get_target_private_dir which is e.g. 'subdir/target.p'.
    """

    def __init__(self, fname: str):
        self.fname = fname

    def __str__(self) -> str:
        return self.fname

class FileMaybeInTargetPrivateDir:
    """Union between 'File' and 'FileInTargetPrivateDir'"""

    def __init__(self, inner: T.Union[File, FileInTargetPrivateDir]):
        self.inner = inner

    @property
    def fname(self) -> str:
        return self.inner.fname

    def rel_to_builddir(self, build_to_src: str, target_private_dir: str) -> str:
        if isinstance(self.inner, FileInTargetPrivateDir):
            return os.path.join(target_private_dir, self.inner.fname)
        return self.inner.rel_to_builddir(build_to_src)

    def absolute_path(self, srcdir: str, builddir: str) -> str:
        if isinstance(self.inner, FileInTargetPrivateDir):
            raise RuntimeError('Unreachable code')
        return self.inner.absolute_path(srcdir, builddir)

    def __str__(self) -> str:
        return self.fname

class Generator(HoldableObject):
    def __init__(self, env: Environment,
                 exe: programs.Program,
                 arguments: T.List[str],
                 output: T.List[str],
                 # how2dataclass
                 *,
                 depfile: T.Optional[str] = None,
                 capture: bool = False,
                 depends: T.Optional[T.List[BuildTargetTypes]] = None,
                 name: str = 'Generator'):
        self.environment = env
        self.exe = exe
        self.depfile = depfile
        self.capture = capture
        self.depends: T.List[BuildTargetTypes] = depends or []
        self.arglist = arguments
        self.outputs = output
        self.name = name

    def __repr__(self) -> str:
        repr_str = "<{0}: {1}>"
        return repr_str.format(self.__class__.__name__, self.exe)

    def get_exe(self) -> programs.Program:
        return self.exe

    def get_base_outnames(self, inname: str) -> T.List[str]:
        plainname = os.path.basename(inname)
        basename = os.path.splitext(plainname)[0]
        bases = [x.replace('@BASENAME@', basename).replace('@PLAINNAME@', plainname) for x in self.outputs]
        return bases

    def get_dep_outname(self, inname: str) -> T.List[str]:
        if self.depfile is None:
            raise InvalidArguments('Tried to get dep name for rule that does not have dependency file defined.')
        plainname = os.path.basename(inname)
        basename = os.path.splitext(plainname)[0]
        return self.depfile.replace('@BASENAME@', basename).replace('@PLAINNAME@', plainname)

    def get_arglist(self, inname: str) -> T.List[str]:
        plainname = os.path.basename(inname)
        basename = os.path.splitext(plainname)[0]
        return [x.replace('@BASENAME@', basename).replace('@PLAINNAME@', plainname) for x in self.arglist]

    def process_files(self, files: T.Iterable[T.Union[str, File, GeneratedTypes]],
                      subdir: str = '',
                      preserve_path_from: T.Optional[str] = None,
                      extra_args: T.Optional[T.List[str]] = None,
                      env: T.Optional[EnvironmentVariables] = None) -> 'GeneratedList':
        output = GeneratedList(
            self,
            subdir,
            preserve_path_from,
            extra_args=extra_args if extra_args is not None else [],
            env=env if env is not None else EnvironmentVariables())

        for e in files:
            if isinstance(e, (CustomTarget, CustomTargetIndex)):
                output.depends.add(e)
                fs = [File.from_built_file(e.get_builddir(), f) for f in e.get_outputs()]
            elif isinstance(e, GeneratedList):
                if preserve_path_from:
                    raise InvalidArguments("generator.process: 'preserve_path_from' is not allowed if one input is a 'generated_list'.")
                output.depends.add(e)
                fs = [FileInTargetPrivateDir(f) for f in e.get_outputs()]
            elif isinstance(e, str):
                fs = [File.from_source_file(self.environment.source_dir, subdir, e)]
            else:
                fs = [e]

            for f in fs:
                if preserve_path_from:
                    abs_f = f.absolute_path(self.environment.source_dir, self.environment.build_dir)
                    if not is_parent_path(preserve_path_from, abs_f):
                        raise InvalidArguments('generator.process: When using preserve_path_from, all input files must be in a subdirectory of the given dir.')
                f = FileMaybeInTargetPrivateDir(f)
                output.add_file(f, self.environment)
        return output


@dataclass(eq=False)
class GeneratedList(HoldableObject):

    """The output of generator.process."""

    generator: Generator
    subdir: str
    preserve_path_from: T.Optional[str]
    extra_args: T.List[str]
    env: T.Optional[EnvironmentVariables]

    def __post_init__(self) -> None:
        self.name = self.generator.exe
        self.depends: T.Set[GeneratedTypes] = set()
        self.infilelist: T.List[FileMaybeInTargetPrivateDir] = []
        self.outfilelist: T.List[str] = []
        self.outmap: T.Dict[FileMaybeInTargetPrivateDir, T.List[str]] = {}
        self.extra_depends: T.List[BuildTargetTypes] = []
        self.depend_files: T.List[File] = []

        if self.extra_args is None:
            self.extra_args: T.List[str] = []

        if self.env is None:
            self.env: EnvironmentVariables = EnvironmentVariables()

        if isinstance(self.generator.exe, programs.Program):
            if not self.generator.exe.found():
                raise InvalidArguments('Tried to use not-found external program as generator')
        if isinstance(self.generator.exe, LocalProgram):
            self.extra_depends.append(self.generator.exe.program)
        else:
            path = self.generator.exe.get_path()
            if os.path.isabs(path):
                # Can only add a dependency on an external program which we
                # know the absolute path of
                self.depend_files.append(File.from_absolute_file(path))

    def add_preserved_path_segment(self, infile: FileMaybeInTargetPrivateDir, outfiles: T.List[str], environment: Environment) -> T.List[str]:
        result: T.List[str] = []
        in_abs = infile.absolute_path(environment.source_dir, environment.build_dir)
        assert os.path.isabs(self.preserve_path_from)
        rel = os.path.relpath(in_abs, self.preserve_path_from)
        path_segment = os.path.dirname(rel)
        for of in outfiles:
            result.append(os.path.join(path_segment, of))
        return result

    def add_file(self, newfile: FileMaybeInTargetPrivateDir, environment: Environment) -> None:
        self.infilelist.append(newfile)
        outfiles = self.generator.get_base_outnames(newfile.fname)
        if self.preserve_path_from:
            outfiles = self.add_preserved_path_segment(newfile, outfiles, environment)
        self.outfilelist += outfiles
        self.outmap[newfile] = outfiles

    def get_inputs(self) -> T.List[FileMaybeInTargetPrivateDir]:
        return self.infilelist

    def get_outputs(self) -> T.List[str]:
        return self.outfilelist

    def get_outputs_for(self, filename: FileMaybeInTargetPrivateDir) -> T.List[str]:
        return self.outmap[filename]

    def get_generator(self) -> 'Generator':
        return self.generator

    def get_extra_args(self) -> T.List[str]:
        return self.extra_args

    def get_subdir(self) -> str:
        return self.subdir


class Executable(BuildTarget):
    known_kwargs = known_exe_kwargs

    typename = 'executable'

    def __init__(
            self,
            name: str,
            subdir: str,
            subproject: SubProject,
            for_machine: MachineChoice,
            sources: T.List['SourceOutputs'],
            structured_sources: T.Optional[StructuredSources],
            objects: T.List[ObjectTypes],
            environment: Environment,
            compilers: CompilerDict,
            kwargs: ExecutableKeywordArguments):
        self.export_dynamic = kwargs.get('export_dynamic', False)
        self.rust_crate_type = kwargs.get('rust_crate_type', 'bin')
        super().__init__(name, subdir, subproject, for_machine, sources, structured_sources, objects,
                         environment, compilers, kwargs)
        self.win_subsystem = kwargs.get('win_subsystem') or 'console'
        self.pie = self._extract_pic_pie(kwargs, 'pie', 'b_pie')
        # Check for export_dynamic
        self.implib_name = kwargs.get('implib')
        # Only linkwithable if using export_dynamic
        self.is_linkwithable = self.export_dynamic

        self.vs_module_defs: T.Optional[File] = None
        self.process_vs_module_defs_kw(kwargs)

    def _set_vala_args(self, kwargs: BuildTargetKeywordArguments) -> None:
        # These don't get generated if the executable doesn't have
        # export_dynamic set to true.
        if self.export_dynamic:
            super()._set_vala_args(kwargs)

    def post_init(self) -> None:
        super().post_init()
        machine = self.environment.machines[self.for_machine]
        # Unless overridden, executables have no suffix or prefix. Except on
        # Windows and with C#/Mono executables where the suffix is 'exe'
        if not hasattr(self, 'prefix'):
            self.prefix = ''
        if not hasattr(self, 'suffix'):
            # Executable for Windows or C#/Mono
            if machine.is_windows() or machine.is_cygwin() or 'cs' in self.compilers:
                self.suffix = 'exe'
            elif machine.system.startswith('wasm') or machine.system == 'emscripten':
                self.suffix = 'js'
            elif ('c' in self.compilers and self.compilers['c'].get_id().startswith('armclang') or
                  'cpp' in self.compilers and self.compilers['cpp'].get_id().startswith('armclang')):
                self.suffix = 'axf'
            elif ('c' in self.compilers and self.compilers['c'].get_id().startswith('ccrx') or
                  'cpp' in self.compilers and self.compilers['cpp'].get_id().startswith('ccrx')):
                self.suffix = 'abs'
            elif ('c' in self.compilers and self.compilers['c'].get_id().startswith('xc16')):
                self.suffix = 'elf'
            elif ('c' in self.compilers and self.compilers['c'].get_id() in {'ti', 'c2000', 'c6000'} or
                  'cpp' in self.compilers and self.compilers['cpp'].get_id() in {'ti', 'c2000', 'c6000'}):
                self.suffix = 'out'
            elif ('c' in self.compilers and self.compilers['c'].get_id() in {'mwccarm', 'mwcceppc'} or
                  'cpp' in self.compilers and self.compilers['cpp'].get_id() in {'mwccarm', 'mwcceppc'}):
                self.suffix = 'nef'
            elif ('c' in self.compilers and self.compilers['c'].get_id() == 'tasking'):
                self.suffix = 'elf'
            else:
                self.suffix = machine.get_exe_suffix()
        self.filename = self.name
        if self.prefix:
            self.filename = self.prefix + self.filename
        if self.suffix:
            self.filename += '.' + self.suffix
        self.outputs[0] = self.filename

        # The import library this target will generate
        self.import_filename = None
        # The debugging information file this target will generate
        self.debug_filename = None

        # If using export_dynamic, set the import library name
        if self.export_dynamic:
            implib_basename = self.implib_name or self.name + '.exe'
            if machine.is_windows() or machine.is_cygwin():
                if self.get_using_msvc():
                    self.import_filename = f'{implib_basename}.lib'
                else:
                    self.import_filename = f'lib{implib_basename}.a'

        create_debug_file = (
            machine.is_windows()
            and ('cs' in self.compilers or self.uses_rust() or self.get_using_msvc())
            # .pdb file is created only when debug symbols are enabled
            and self.environment.coredata.optstore.get_value_for(OptionKey("debug"))
        )
        if create_debug_file:
            # If the target is has a standard exe extension (i.e. 'foo.exe'),
            # then the pdb name simply becomes 'foo.pdb'. If the extension is
            # something exotic, then include that in the name for uniqueness
            # reasons (e.g. 'foo_com.pdb').
            name = self.name
            if getattr(self, 'suffix', 'exe') != 'exe':
                name += '_' + self.suffix
            self.debug_filename = name + '.pdb'

    def get_default_install_dir(self) -> T.Tuple[str, str]:
        return self.environment.get_bindir(), '{bindir}'

    def type_suffix(self) -> str:
        return "@exe"

    def get_import_filename(self) -> T.Optional[str]:
        """
        The name of the import library that will be outputted by the compiler

        Returns None if there is no import library required for this platform
        """
        return self.import_filename

    def get_debug_filename(self) -> T.Optional[str]:
        """
        The name of debuginfo file that will be created by the compiler

        Returns None if the build won't create any debuginfo file
        """
        return self.debug_filename

    def is_linkable_target(self) -> bool:
        return self.is_linkwithable


class StaticLibrary(BuildTarget):
    known_kwargs = known_stlib_kwargs

    typename = 'static library'

    def __init__(
            self,
            name: str,
            subdir: str,
            subproject: SubProject,
            for_machine: MachineChoice,
            sources: T.List['SourceOutputs'],
            structured_sources: T.Optional[StructuredSources],
            objects: T.List[ObjectTypes],
            environment: Environment,
            compilers: CompilerDict,
            kwargs: StaticLibraryKeywordArguments):
        self.prelink = kwargs.get('prelink', False)
        self.rust_crate_type = kwargs.get('rust_crate_type', 'rlib')
        super().__init__(name, subdir, subproject, for_machine, sources, structured_sources, objects,
                         environment, compilers, kwargs)
        self.pic = self._extract_pic_pie(kwargs, 'pic', 'b_staticpic')
        if not self.pic:
            self.pie = self._extract_pic_pie(kwargs, 'pie', 'b_pie')

    def post_init(self) -> None:
        super().post_init()
        if 'cs' in self.compilers:
            raise InvalidArguments('Static libraries not supported for C#.')
        if self.uses_rust():
            # See https://github.com/rust-lang/rust/issues/110460
            if self.rust_crate_type == 'rlib' and any(c in self.name for c in ['-', ' ', '.']):
                raise InvalidArguments(f'Rust crate {self.name} type {self.rust_crate_type} does not allow spaces, '
                                       'periods or dashes in the library name due to a limitation of rustc. '
                                       'Replace them with underscores, for example')
            if self.rust_crate_type == 'staticlib':
                # FIXME: In the case of no-std we should not add those libraries,
                # but we have no way to know currently.

                # XXX:
                #  In the case of no-std, we are likely in a bare metal case
                #  and thus, machine_info kernel should be set to 'none'.
                #  In that case, native_static_libs list is empty.
                rustc = self.compilers['rust']
                if rustc.get_crt_static():
                    # musl targets need self-contained for libc.a, libunwind.a etc.
                    link_args = ['-L' + rustc.get_target_libdir() + '/self-contained']
                else:
                    # Avoid embedding self-contained libc.a into shared libraries —
                    # creates a second musl instance with uninitialized __libc state.
                    link_args = []
                link_args += rustc.native_static_libs
                d = dependencies.InternalDependency(version='undefined', link_args=link_args,
                                                    name='_rust_native_static_libs')
                self.external_deps.append(d)

        default_prefix, default_suffix = self.determine_default_prefix_and_suffix()
        if not self.name_prefix_set:
            self.prefix = default_prefix
        if not self.name_suffix_set:
            self.suffix = default_suffix
        self.filename = self.prefix + self.name + '.' + self.suffix
        self.outputs[0] = self.filename

    def determine_default_prefix_and_suffix(self) -> T.Tuple[str, str]:
        scheme = self.environment.coredata.get_option_for_target(self, 'namingscheme')
        assert isinstance(scheme, str), 'for mypy'
        if scheme == 'platform':
            schemename = self.get_platform_scheme_name()
            prefix, suffix = DEFAULT_STATIC_LIBRARY_NAMES[schemename]
        else:
            prefix = ''
            suffix = ''
        # By default a static library is named libfoo.a even on Windows because
        # MSVC does not have a consistent convention for what static libraries
        # are called. The MSVC CRT uses libfoo.lib syntax but nothing else uses
        # it and GCC only looks for static libraries called foo.lib and
        # libfoo.a. However, we cannot use foo.lib because that's the same as
        # the import library. Using libfoo.a is ok because people using MSVC
        # always pass the library filename while linking anyway.
        #
        # See our FAQ for more detailed rationale:
        # https://mesonbuild.com/FAQ.html#why-does-building-my-project-with-msvc-output-static-libraries-called-libfooa
        if not hasattr(self, 'prefix'):
            prefix = 'lib'
        if not hasattr(self, 'suffix'):
            if self.uses_rust():
                if self.rust_crate_type == 'rlib':
                    # default Rust static library suffix
                    suffix = 'rlib'
                elif self.rust_crate_type == 'staticlib':
                    suffix = 'a'
            elif self.environment.machines[self.for_machine].is_os2() and self.environment.coredata.optstore.get_value_for(OptionKey('os2_emxomf')):
                suffix = 'lib'
            else:
                suffix = 'a'
                if 'c' in self.compilers and self.compilers['c'].get_id() == 'tasking' and not self.prelink:
                    key = OptionKey('b_lto', self.subproject, self.for_machine)
                    v = self.environment.coredata.get_option_for_target(self, key)
                    assert isinstance(v, bool), 'for mypy'
                    if v:
                        suffix = 'ma'
        return (prefix, suffix)

    def get_link_deps_mapping(self, prefix: str) -> T.Mapping[str, str]:
        return {}

    def get_default_install_dir(self) -> T.Tuple[str, str]:
        return self.environment.get_static_lib_dir(), '{libdir_static}'

    def type_suffix(self) -> str:
        return "@rlib" if self.uses_rust_abi() else "@sta"

    def is_linkable_target(self) -> bool:
        return True

    def is_internal(self) -> bool:
        return not self.install

    def set_shared(self, shared_library: SharedLibrary) -> None:
        self.both_lib = copy.copy(shared_library)
        self.both_lib.both_lib = None

    def get(self, lib_type: T.Literal['static', 'shared'], recursive: bool = False) -> LibTypes:
        result = self
        if lib_type == 'shared':
            result = self.both_lib or self
        if recursive:
            result.link_targets = [t.get(lib_type, True) for t in self.link_targets]
        return result

    def link_whole(
            self,
            targets: T.List[StaticTargetTypes],
            promoted: bool = False) -> None:
        for t in targets:
            self.check_can_link_together(t)

            # When we're a static library and we link_whole: to another static
            # library, we need to add that target's objects to ourselves.
            self._bundle_static_library(t, promoted)

            # If we install this static library we also need to include objects
            # from all uninstalled static libraries it depends on.
            if self.install:
                for lib in t.get_internal_static_libraries():
                    self._bundle_static_library(lib, True)
            self.link_whole_targets.append(t)

    def link(self, targets: T.List[BuildTargetTypes]) -> None:
        for t in targets:
            if self.install and t.is_internal():
                # When we're a static library and we link_with to an
                # internal/convenience library, promote to link_whole.
                self.link_whole([t], promoted=True)
                continue
            self.check_can_link_together(t)
            self.link_targets.append(t)

    def _bundle_static_library(self, t: StaticTargetTypes, promoted: bool = False) -> None:
        if self.uses_rust():
            # Rustc can bundle static libraries, no need to extract objects.
            self.link_whole_targets.append(t)
        elif isinstance(t, (CustomTarget, CustomTargetIndex)) or t.uses_rust():
            # To extract objects from a custom target we would have to extract
            # the archive, WIP implementation can be found in
            # https://github.com/mesonbuild/meson/pull/9218.
            # For Rust C ABI we could in theory have access to objects, but we
            # don't currently build them in such a way that this is possible:
            # https://github.com/mesonbuild/meson/issues/10724
            m = (f'Cannot link_whole a custom or Rust target {t.name!r} into a static library {self.name!r}. '
                 'Instead, pass individual object files with the "objects:" keyword argument if possible.')
            if promoted:
                m += (f' Meson had to promote link to link_whole because {self.name!r} is installed but not {t.name!r},'
                      f' and thus has to include objects from {t.name!r} to be usable.')
            raise InvalidArguments(m)
        else:
            self.objects.append(t.extract_all_objects())

class SharedLibrary(BuildTarget):
    known_kwargs = known_shlib_kwargs

    typename = 'shared library'

    # Used by AIX to decide whether to archive shared library or not.
    aix_so_archive = True

    def __init__(
            self,
            name: str,
            subdir: str,
            subproject: SubProject,
            for_machine: MachineChoice,
            sources: T.List['SourceOutputs'],
            structured_sources: T.Optional[StructuredSources],
            objects: T.List[ObjectTypes],
            environment: Environment,
            compilers: CompilerDict,
            kwargs):
        self.soversion: T.Optional[str] = None
        self.ltversion: T.Optional[str] = None
        # Max length 2, first element is compatibility_version, second is current_version
        self.darwin_versions: T.Optional[T.Tuple[str, str]] = None
        self.vs_module_defs = None
        self.shortname: T.Optional[str] = None
        # The import library this target will generate
        self.import_filename = None
        # The debugging information file this target will generate
        self.debug_filename = None
        # Use by the pkgconfig module
        self.shared_library_only = False
        self.rust_crate_type = kwargs.get('rust_crate_type', 'dylib')
        super().__init__(name, subdir, subproject, for_machine, sources, structured_sources, objects,
                         environment, compilers, kwargs)

    def post_init(self) -> None:
        super().post_init()
        if self.uses_rust():
            # See https://github.com/rust-lang/rust/issues/110460
            if self.rust_crate_type != 'cdylib' and any(c in self.name for c in ['-', ' ', '.']):
                raise InvalidArguments(f'Rust crate {self.name} type {self.rust_crate_type} does not allow spaces, '
                                       'periods or dashes in the library name due to a limitation of rustc. '
                                       'Replace them with underscores, for example')

        if not hasattr(self, 'prefix'):
            self.prefix = None
        if not hasattr(self, 'suffix'):
            self.suffix = None
        self.basic_filename_tpl = '{0.prefix}{0.name}.{0.suffix}'
        self.determine_filenames()

    def get_link_deps_mapping(self, prefix: str) -> T.Mapping[str, str]:
        result: T.Dict[str, str] = {}
        mappings = self.get_transitive_link_deps_mapping(prefix)
        old = get_target_macos_dylib_install_name(self)
        if old not in mappings:
            fname = self.get_filename()
            install_dir = self.install_dir[0]
            assert install_dir is not False, 'This is a bug'
            new = os.path.join(prefix, install_dir, fname)
            result.update({old: new})
        mappings.update(result)
        return mappings

    def get_default_install_dir(self) -> T.Tuple[str, str]:
        return self.environment.get_shared_lib_dir(), '{libdir_shared}'

    def determine_naming_info(self) -> T.Tuple[str, str, str, str, bool]:
        scheme = self.environment.coredata.get_option_for_target(self, 'namingscheme')
        assert isinstance(scheme, str), 'for mypy'
        if scheme == 'platform':
            schemename = self.get_platform_scheme_name()
            prefix, suffix, import_suffix = DEFAULT_SHARED_LIBRARY_NAMES[schemename]
        else:
            prefix = None
            suffix = None
            import_suffix = None
        filename_tpl = self.basic_filename_tpl
        create_debug_file = False
        create_debug_file = False
        self.filename_tpl = self.basic_filename_tpl
        import_filename_tpl = None
        # NOTE: manual prefix/suffix override is currently only tested for C/C++
        # C# and Mono
        if 'cs' in self.compilers:
            prefix = ''
            suffix = 'dll'
            filename_tpl = '{0.prefix}{0.name}.{0.suffix}'
            create_debug_file = True
        # C, C++, Swift, Vala
        # Only Windows and OS/2 uses a separate import library for linking
        # For all other targets/platforms import_filename stays None
        elif self.environment.machines[self.for_machine].is_windows():
            suffix = suffix if suffix is not None else 'dll'
            if self.uses_rust():
                # Shared library is of the form foo.dll
                prefix = prefix if prefix is not None else ''
                # Import library is called foo.dll.lib
                import_filename_tpl = '{0.prefix}{0.name}.dll.lib'
                # .pdb file is only created when debug symbols are enabled
                create_debug_file = self.environment.coredata.optstore.get_value_for(OptionKey("debug"))
            elif self.get_using_msvc():
                # Shared library is of the form foo.dll
                prefix = prefix if prefix is not None else ''
                # Import library is called foo.lib
                import_suffix = import_suffix if import_suffix is not None else 'lib'
                import_filename_tpl = '{0.prefix}{0.name}.' + import_suffix
                # .pdb file is only created when debug symbols are enabled
                create_debug_file = self.environment.coredata.optstore.get_value_for(OptionKey("debug"))
            # Assume GCC-compatible naming
            else:
                # Shared library is of the form libfoo.dll
                prefix = prefix if prefix is not None else 'lib'
                # Import library is called libfoo.dll.a
                import_suffix = import_suffix if import_suffix is not None else '.dll.a'
                import_filename_tpl = '{0.prefix}{0.name}' + import_suffix
            # Shared library has the soversion if it is defined
            if self.soversion:
                filename_tpl = '{0.prefix}{0.name}-{0.soversion}.{0.suffix}'
            else:
                filename_tpl = '{0.prefix}{0.name}.{0.suffix}'
        elif self.environment.machines[self.for_machine].is_cygwin():
            suffix = 'dll'
            # Shared library is of the form cygfoo.dll
            # (ld --dll-search-prefix=cyg is the default)
            prefix = 'cyg'
            # Import library is called libfoo.dll.a
            import_suffix = import_suffix if import_suffix is not None else '.dll.a'
            import_prefix = self.prefix if self.prefix is not None else 'lib'
            import_filename_tpl = import_prefix + '{0.name}' + import_suffix
            if self.soversion:
                filename_tpl = '{0.prefix}{0.name}-{0.soversion}.{0.suffix}'
            else:
                filename_tpl = '{0.prefix}{0.name}.{0.suffix}'
        elif self.environment.machines[self.for_machine].is_darwin():
            prefix = prefix if prefix is not None else 'lib'
            suffix = suffix if suffix is not None else 'dylib'
            # On macOS, the filename can only contain the major version
            if self.soversion:
                # libfoo.X.dylib
                filename_tpl = '{0.prefix}{0.name}.{0.soversion}.{0.suffix}'
            else:
                # libfoo.dylib
                filename_tpl = '{0.prefix}{0.name}.{0.suffix}'
        elif self.environment.machines[self.for_machine].is_android():
            prefix = prefix if prefix is not None else 'lib'
            suffix = suffix if suffix is not None else 'so'
            # Android doesn't support shared_library versioning
            filename_tpl = '{0.prefix}{0.name}.{0.suffix}'
        elif self.environment.machines[self.for_machine].is_os2():
            # Shared library is of the form foo.dll
            prefix = prefix if prefix is not None else ''
            suffix = suffix if suffix is not None else 'dll'
            # Import library is called foo_dll.a or foo_dll.lib
            if import_suffix is None:
                import_suffix = '_dll.lib' if self.environment.coredata.optstore.get_value_for(OptionKey('os2_emxomf')) else '_dll.a'
            import_filename_tpl = '{0.prefix}{0.name}' + import_suffix
            filename_tpl = '{0.shortname}' if self.shortname else '{0.prefix}{0.name}'
            if self.soversion:
                # fooX.dll
                filename_tpl += '{0.soversion}'
            filename_tpl += '.{0.suffix}'
        else:
            prefix = prefix if prefix is not None else 'lib'
            suffix = suffix if suffix is not None else 'so'
            if self.ltversion:
                # libfoo.so.X[.Y[.Z]] (.Y and .Z are optional)
                filename_tpl = '{0.prefix}{0.name}.{0.suffix}.{0.ltversion}'
            elif self.soversion:
                # libfoo.so.X
                filename_tpl = '{0.prefix}{0.name}.{0.suffix}.{0.soversion}'
            else:
                # No versioning, libfoo.so
                filename_tpl = '{0.prefix}{0.name}.{0.suffix}'
        return (prefix, suffix, filename_tpl, import_filename_tpl, create_debug_file)

    def determine_filenames(self) -> None:
        """
        See https://github.com/mesonbuild/meson/pull/417 for details.

        First we determine the filename template (self.filename_tpl), then we
        set the output filename (self.filename).

        The template is needed while creating aliases (self.get_aliases),
        which are needed while generating .so shared libraries for Linux.

        Besides this, there's also the import library name (self.import_filename),
        which is only used on Windows and OS/2 since on that platform the linker uses a
        separate library called the "import library" during linking instead of
        the shared library (DLL).
        """
        prefix, suffix, filename_tpl, import_filename_tpl, create_debug_file = self.determine_naming_info()
        if self.prefix is None:
            self.prefix = prefix
        if self.suffix is None:
            self.suffix = suffix
        self.filename_tpl = filename_tpl
        self.filename = self.filename_tpl.format(self)
        if self.environment.machines[self.for_machine].is_os2():
            # OS/2 does not allow a longer DLL name than 8 chars
            name = os.path.splitext(self.filename)[0]
            if len(name) > 8:
                name = name[:8]
                if self.soversion:
                    name = name[:-len(self.soversion)] + self.soversion
            self.filename = '{}.{}'.format(name, self.suffix)
        if import_filename_tpl:
            self.import_filename = import_filename_tpl.format(self)
        # There may have been more outputs added by the time we get here, so
        # only replace the first entry
        self.outputs[0] = self.filename
        if create_debug_file:
            self.debug_filename = os.path.splitext(self.filename)[0] + '.pdb'

    def process_kwargs(self, kwargs: SharedLibraryKeywordArguments) -> None:
        super().process_kwargs(kwargs)

        if not self.environment.machines[self.for_machine].is_android():
            # Shared library version
            self.ltversion = kwargs.get('version')
            self.soversion = kwargs.get('soversion')
            if self.soversion is None and self.ltversion is not None:
                # library version is defined, get the soversion from that
                # We replicate what Autotools does here and take the first
                # number of the version by default.
                self.soversion = self.ltversion.split('.')[0]
            # macOS, iOS and tvOS dylib compatibility_version and current_version
            self.darwin_versions = kwargs.get('darwin_versions')
            if self.darwin_versions is None and self.soversion is not None:
                # If unspecified, pick the soversion
                self.darwin_versions = (self.soversion, self.soversion)

        # Visual Studio module-definitions file
        self.process_vs_module_defs_kw(kwargs)

        # OS/2 uses a 8.3 name for a DLL
        self.shortname = kwargs.get('shortname')

    def get_import_filename(self) -> T.Optional[str]:
        """
        The name of the import library that will be outputted by the compiler

        Returns None if there is no import library required for this platform
        """
        return self.import_filename

    def get_debug_filename(self) -> T.Optional[str]:
        """
        The name of debuginfo file that will be created by the compiler

        Returns None if the build won't create any debuginfo file
        """
        return self.debug_filename

    def get_aliases(self) -> T.List[T.Tuple[str, str, str]]:
        """
        If the versioned library name is libfoo.so.0.100.0, aliases are:
        * libfoo.so.0 (soversion) -> libfoo.so.0.100.0
        * libfoo.so (unversioned; for linking) -> libfoo.so.0
        Same for dylib:
        * libfoo.dylib (unversioned; for linking) -> libfoo.0.dylib
        """
        aliases: T.List[T.Tuple[str, str, str]] = []
        # Aliases are only useful with .so and .dylib libraries. Also if
        # there's no self.soversion (no versioning), we don't need aliases.
        if self.suffix not in ('so', 'dylib') or not self.soversion:
            return aliases
        # With .so libraries, the minor and micro versions are also in the
        # filename. If ltversion != soversion we create an soversion alias:
        # libfoo.so.0 -> libfoo.so.0.100.0
        # Where libfoo.so.0.100.0 is the actual library
        if self.suffix == 'so' and self.ltversion and self.ltversion != self.soversion:
            alias_tpl = self.filename_tpl.replace('ltversion', 'soversion')
            ltversion_filename = alias_tpl.format(self)
            tag = self.install_tag[0] or 'runtime'
            aliases.append((ltversion_filename, self.filename, tag))
        # libfoo.so.0/libfoo.0.dylib is the actual library
        else:
            ltversion_filename = self.filename
        # Unversioned alias:
        #  libfoo.so -> libfoo.so.0
        #  libfoo.dylib -> libfoo.0.dylib
        tag = self.install_tag[0] or 'devel'
        aliases.append((self.basic_filename_tpl.format(self), ltversion_filename, tag))
        return aliases

    def type_suffix(self) -> str:
        return "@sha"

    def is_linkable_target(self) -> bool:
        return True

    def set_static(self, static_library: StaticLibrary) -> None:
        self.both_lib = copy.copy(static_library)
        self.both_lib.both_lib = None

    def get(self, lib_type: T.Literal['static', 'shared'], recursive: bool = False) -> LibTypes:
        result = self
        if lib_type == 'static':
            result = self.both_lib or self
        if recursive:
            result.link_targets = [t.get(lib_type, True) for t in self.link_targets]
        return result

    def link_whole(
            self,
            targets: T.List[StaticTargetTypes],
            promoted: bool = False) -> None:
        for t in targets:
            self.check_can_link_together(t)
            if not getattr(t, 'pic', True):
                msg = f"Can't link non-PIC static library {t.name!r} into shared library {self.name!r}. "
                msg += "Use the 'pic' option to static_library to build with PIC."
                raise InvalidArguments(msg)
            self.link_whole_targets.append(t)

    def link(self, targets: T.List[BuildTargetTypes]) -> None:
        for t in targets:
            if isinstance(t, StaticLibrary) and not t.pic:
                msg = f"Can't link non-PIC static library {t.name!r} into shared library {self.name!r}. "
                msg += "Use the 'pic' option to static_library to build with PIC."
                raise InvalidArguments(msg)
            self.check_can_link_together(t)
            self.link_targets.append(t)

# A shared library that is meant to be used with dlopen rather than linking
# into something else.
class SharedModule(SharedLibrary):
    known_kwargs = known_shmod_kwargs

    typename = 'shared module'

    # Used by AIX to not archive shared library for dlopen mechanism
    aix_so_archive = False

    def __init__(
            self,
            name: str,
            subdir: str,
            subproject: SubProject,
            for_machine: MachineChoice,
            sources: T.List['SourceOutputs'],
            structured_sources: T.Optional[StructuredSources],
            objects: T.List[ObjectTypes],
            environment: Environment,
            compilers: CompilerDict,
            kwargs):
        if 'version' in kwargs:
            raise MesonException('Shared modules must not specify the version kwarg.')
        if 'soversion' in kwargs:
            raise MesonException('Shared modules must not specify the soversion kwarg.')
        super().__init__(name, subdir, subproject, for_machine, sources,
                         structured_sources, objects, environment, compilers, kwargs)
        # We need to set the soname in cases where build files link the module
        # to build targets, see: https://github.com/mesonbuild/meson/issues/9492
        self.force_soname = False

    def get_default_install_dir(self) -> T.Tuple[str, str]:
        return self.environment.get_shared_module_dir(), '{moduledir_shared}'

class BothLibraries(SecondLevelHolder):
    def __init__(self, shared: SharedLibrary, static: StaticLibrary, preferred_library: Literal['shared', 'static']) -> None:
        self._preferred_library = preferred_library
        self.shared = shared
        self.static = static
        self.subproject = self.shared.subproject

    def __repr__(self) -> str:
        return f'<BothLibraries: static={repr(self.static)}; shared={repr(self.shared)}>'

    def get(self, lib_type: T.Literal['static', 'shared']) -> T.Union[StaticLibrary, SharedLibrary]:
        if lib_type == 'static':
            return self.static
        if lib_type == 'shared':
            return self.shared
        return self.get_default_object()

    def get_default_object(self) -> T.Union[StaticLibrary, SharedLibrary]:
        if self._preferred_library == 'shared':
            return self.shared
        elif self._preferred_library == 'static':
            return self.static
        raise MesonBugException(f'self._preferred_library == "{self._preferred_library}" is neither "shared" nor "static".')

    def get_id(self) -> str:
        return self.get_default_object().get_id()

    def is_linkable_target(self) -> bool:
        # For polymorphism with build targets
        return True

class CommandBase:

    depend_files: T.List[File]
    dependencies: T.List[T.Union[BuildTarget, 'CustomTarget']]
    subproject: str

    def flatten_command(self, cmd: T.Sequence[T.Union[str, File, programs.Program, BuildTargetTypes]]) -> \
            T.List[T.Union[str, File, BuildTarget, CustomTarget, programs.Program]]:
        cmd = listify(cmd)
        final_cmd: T.List[T.Union[str, File, BuildTarget, 'CustomTarget']] = []
        for c in cmd:
            if isinstance(c, LocalProgram):
                c = c.program
            if isinstance(c, str):
                final_cmd.append(c)
            elif isinstance(c, File):
                self.depend_files.append(c)
                final_cmd.append(c)
            elif isinstance(c, programs.Program):
                if not c.found():
                    raise InvalidArguments('Tried to use not-found external program in "command"')
                path = c.get_path()
                if os.path.isabs(path):
                    # Can only add a dependency on an external program which we
                    # know the absolute path of
                    self.depend_files.append(File.from_absolute_file(path))
                # Do NOT flatten -- it is needed for later parsing
                final_cmd.append(c)
            elif isinstance(c, (BuildTarget, CustomTarget)):
                self.dependencies.append(c)
                final_cmd.append(c)
            elif isinstance(c, CustomTargetIndex):
                FeatureNew.single_use('CustomTargetIndex for command argument', '0.60', self.subproject)
                self.dependencies.append(c.target)
                final_cmd += self.flatten_command(File.from_built_file(c.get_subdir(), c.get_filename()))
            elif isinstance(c, list):
                final_cmd += self.flatten_command(c)
            else:
                raise InvalidArguments(f'Argument {c!r} in "command" is invalid')
        return final_cmd

class CustomTargetBase:
    ''' Base class for CustomTarget and CustomTargetIndex

    This base class can be used to provide a dummy implementation of some
    private methods to avoid repeating `isinstance(t, BuildTarget)` when dealing
    with custom targets.
    '''

    rust_crate_type = ''

    def get_dependencies_recurse(self, result: OrderedSet[BuildTargetTypes], visited: T.Set[BuildTargetTypes],
                                 include_internals: bool = True) -> None:
        pass

    def get_internal_static_libraries(self) -> OrderedSet[StaticTargetTypes]:
        return OrderedSet()

    def get_internal_static_libraries_recurse(self, result: OrderedSet[StaticTargetTypes]) -> None:
        pass

    def get_all_linked_targets(self) -> ImmutableListProtocol[BuildTargetTypes]:
        return []

    def get(self, lib_type: T.Literal['static', 'shared'], recursive: bool = False) -> LibTypes:
        """Base case used by BothLibraries"""
        return self

class CustomTarget(Target, CustomTargetBase, CommandBase):

    typename = 'custom'

    def __init__(self,
                 name: T.Optional[str],
                 subdir: str,
                 subproject: str,
                 environment: Environment,
                 command: T.Sequence[T.Union[
                     str, BuildTargetTypes, GeneratedList,
                     programs.Program, File]],
                 sources: T.Sequence[T.Union[
                     str, File, BuildTargetTypes, ExtractedObjects,
                     GeneratedList, programs.Program]],
                 outputs: T.List[str],
                 *,
                 build_always_stale: bool = False,
                 build_by_default: T.Optional[bool] = None,
                 capture: bool = False,
                 console: bool = False,
                 depend_files: T.Optional[T.Sequence[FileOrString]] = None,
                 extra_depends: T.Optional[T.Sequence[T.Union[str, SourceOutputs]]] = None,
                 depfile: T.Optional[str] = None,
                 depfile_type: T.Optional[Literal['gcc', 'msvc']] = None,
                 env: T.Optional[EnvironmentVariables] = None,
                 feed: bool = False,
                 install: bool = False,
                 install_dir: T.Optional[T.List[T.Union[str, Literal[False]]]] = None,
                 install_mode: T.Optional[FileMode] = None,
                 install_tag: T.Optional[T.List[T.Optional[str]]] = None,
                 rspable: bool = False,
                 absolute_paths: bool = False,
                 backend: T.Optional['Backend'] = None,
                 description: str = 'Generating {} with a custom command',
                 build_subdir: str = '',
                 ):
        # TODO expose keyword arg to make MachineChoice.HOST configurable
        super().__init__(name, subdir, subproject, False, MachineChoice.HOST, environment,
                         install, build_always_stale, build_subdir = build_subdir)
        self.sources = list(sources)
        self.outputs = substitute_values(
            outputs, get_filenames_templates_dict(
                get_sources_string_names(sources, backend),
                []))
        self.build_by_default = build_by_default if build_by_default is not None else install
        self.capture = capture
        self.console = console
        self.depend_files = list(depend_files or [])
        self.dependencies: T.List[T.Union[CustomTarget, BuildTarget]] = []
        # must be after depend_files and dependencies
        self.command = self.flatten_command(command)
        self.depfile = depfile
        self.depfile_type = 'gcc' if depfile else depfile_type
        self.env = env or EnvironmentVariables()
        self.extra_depends = list(extra_depends or [])
        self.feed = feed
        self.install_dir = list(install_dir or [])
        self.has_custom_install_dir = bool(self.install_dir)
        self.install_mode = install_mode
        self.install_tag = _process_install_tag(install_tag, len(self.outputs))
        self.name = name if name else self.outputs[0]
        self.description = description

        # Whether to use absolute paths for all files on the commandline
        self.absolute_paths = absolute_paths

        # Whether to enable using response files for the underlying tool
        self.rspable = rspable

    def install_dir_names(self) -> T.List[T.Optional[str]]:
        install_dir_names: T.List[T.Optional[str]] = []
        if self.has_custom_install_dir:
            install_dir_names = [getattr(i, 'optname', None) for i in self.install_dir]

        return install_dir_names

    def __repr__(self):
        repr_str = "<{0} {1}: {2}>"
        return repr_str.format(self.__class__.__name__, self.get_id(), self.command)

    def get_target_dependencies(self) -> T.List[T.Union[SourceOutputs, str]]:
        deps: T.List[T.Union[SourceOutputs, str]] = []
        deps.extend(self.dependencies)
        deps.extend(self.extra_depends)
        for c in self.sources:
            if isinstance(c, CustomTargetIndex):
                deps.append(c.target)
            elif not isinstance(c, programs.ExternalProgram):
                deps.append(c)
        return deps

    def get_transitive_build_target_deps(self) -> T.Set[T.Union[BuildTarget, 'CustomTarget']]:
        '''
        Recursively fetch the build targets that this custom target depends on,
        whether through `command:`, `depends:`, or `sources:` The recursion is
        only performed on custom targets.
        This is useful for setting PATH on Windows for finding required DLLs.
        F.ex, if you have a python script that loads a C module that links to
        other DLLs in your project.
        '''
        bdeps: T.Set[T.Union[BuildTarget, 'CustomTarget']] = set()
        deps = self.get_target_dependencies()
        for d in deps:
            if isinstance(d, BuildTarget):
                bdeps.add(d)
            elif isinstance(d, CustomTarget):
                bdeps.update(d.get_transitive_build_target_deps())
        return bdeps

    def get_dependencies(self) -> T.List[T.Union[CustomTarget, coredata.BuildTarget]]:
        return self.dependencies

    def should_install(self) -> bool:
        return self.install

    def get_custom_install_dir(self) -> T.List[T.Union[str, Literal[False]]]:
        return self.install_dir

    def get_custom_install_mode(self) -> T.Optional['FileMode']:
        return self.install_mode

    def get_outputs(self) -> T.List[str]:
        return self.outputs

    def get_filename(self) -> str:
        return self.outputs[0]

    def get_sources(self) -> T.List[T.Union[str, File, BuildTarget, GeneratedTypes, ExtractedObjects, programs.Program]]:
        return self.sources

    def get_generated_lists(self) -> T.List[GeneratedList]:
        genlists: T.List[GeneratedList] = []
        for c in self.sources:
            if isinstance(c, GeneratedList):
                genlists.append(c)
        return genlists

    def get_generated_sources(self) -> T.List[GeneratedList]:
        return self.get_generated_lists()

    def get_dep_outname(self, infilenames):
        if self.depfile is None:
            raise InvalidArguments('Tried to get depfile name for custom_target that does not have depfile defined.')
        if infilenames:
            plainname = os.path.basename(infilenames[0])
            basename = os.path.splitext(plainname)[0]
            return self.depfile.replace('@BASENAME@', basename).replace('@PLAINNAME@', plainname)
        else:
            if '@BASENAME@' in self.depfile or '@PLAINNAME@' in self.depfile:
                raise InvalidArguments('Substitution in depfile for custom_target that does not have an input file.')
            return self.depfile

    def is_linkable_output(self, output: str) -> bool:
        if output.endswith(('.a', '.dll', '.lib', '.so', '.dylib')):
            return True
        # libfoo.so.X soname
        if re.search(r'\.so(\.\d+)*$', output):
            return True
        return False

    def is_linkable_target(self) -> bool:
        if len(self.outputs) != 1:
            return False
        return self.is_linkable_output(self.outputs[0])

    def links_dynamically(self) -> bool:
        """Whether this target links dynamically or statically

        Does not assert the target is linkable, just that it is not shared

        :return: True if is dynamically linked, otherwise False
        """
        suf = os.path.splitext(self.outputs[0])[-1]
        return suf not in {'.a', '.lib'}

    def get_link_deps_mapping(self, prefix: str) -> T.Mapping[str, str]:
        return {}

    def get_link_dep_subdirs(self) -> T.AbstractSet[str]:
        return OrderedSet()

    def get_all_link_deps(self) -> ImmutableListProtocol[BuildTargetTypes]:
        return []

    def is_internal(self) -> bool:
        '''
        Returns True if this is a not installed static library.
        '''
        if len(self.outputs) != 1:
            return False
        return CustomTargetIndex(self, self.outputs[0]).is_internal()

    def extract_all_objects(self) -> T.List[T.Union[str, 'ExtractedObjects']]:
        return self.get_outputs()

    def type_suffix(self) -> str:
        return "@cus"

    def __getitem__(self, index: int) -> 'CustomTargetIndex':
        return CustomTargetIndex(self, self.outputs[index])

    def __setitem__(self, index, value):
        raise NotImplementedError

    def __delitem__(self, index):
        raise NotImplementedError

    def __iter__(self):
        for i in self.outputs:
            yield CustomTargetIndex(self, i)

    def __len__(self) -> int:
        return len(self.outputs)

class CompileTarget(BuildTarget):
    '''
    Target that only compile sources without linking them together.
    It can be used as preprocessor, or transpiler.
    '''

    typename = 'compile'

    def __init__(self,
                 name: str,
                 subdir: str,
                 subproject: str,
                 environment: Environment,
                 sources: T.List['SourceOutputs'],
                 output_templ: str,
                 compiler: Compiler,
                 backend: Backend,
                 compile_args: T.List[str],
                 include_directories: T.List[IncludeDirs],
                 dependencies: T.List[dependencies.Dependency],
                 depends: T.List[BuildTargetTypes]):
        compilers = {compiler.get_language(): compiler}
        kwargs = {
            'build_by_default': False,
            'language_args': {compiler.language: compile_args},
            'include_directories': include_directories,
            'dependencies': dependencies,
        }
        super().__init__(name, subdir, subproject, compiler.for_machine,
                         sources, None, [], environment, compilers, kwargs)
        self.filename = name
        self.compiler = compiler
        self.output_templ = output_templ
        self.outputs = []
        self.sources_map: T.Dict[File, str] = {}
        self.depends = list(depends or [])
        for f in self.sources:
            self._add_output(f)
        for gensrc in self.generated:
            for s in gensrc.get_outputs():
                rel_src = backend.get_target_generated_dir(self, gensrc, s)
                self._add_output(File.from_built_relative(rel_src))

    def type_suffix(self) -> str:
        return "@compile"

    def _add_output(self, f: File) -> None:
        plainname = os.path.basename(f.fname)
        basename = os.path.splitext(plainname)[0]
        o = self.output_templ.replace('@BASENAME@', basename).replace('@PLAINNAME@', plainname)
        self.outputs.append(o)
        self.sources_map[f] = o

    def get_generated_headers(self) -> T.List[File]:
        gen_headers: T.List[File] = []
        for dep in self.depends:
            gen_headers += [File(True, dep.subdir, o) for o in dep.get_outputs()]
        return gen_headers

class RunTarget(Target, CommandBase):

    typename = 'run'

    def __init__(self, name: str,
                 command: T.Sequence[T.Union[str, File, BuildTargetTypes, programs.Program]],
                 dependencies: T.Sequence[AnyTargetType],
                 subdir: str,
                 subproject: str,
                 environment: Environment,
                 env: T.Optional[EnvironmentVariables] = None,
                 default_env: bool = True):
        # These don't produce output artifacts
        super().__init__(name, subdir, subproject, False, MachineChoice.BUILD, environment)
        self.dependencies = dependencies
        self.depend_files = []
        self.command = self.flatten_command(command)
        self.absolute_paths = False
        self.env = env
        self.default_env = default_env

    def __repr__(self) -> str:
        repr_str = "<{0} {1}: {2}>"
        return repr_str.format(self.__class__.__name__, self.get_id(), self.command[0])

    def get_dependencies(self) -> T.List[BuildTargetTypes]:
        return self.dependencies

    def get_generated_sources(self) -> T.List[GeneratedTypes]:
        return []

    def get_sources(self) -> T.List[File]:
        return []

    def should_install(self) -> bool:
        return False

    def get_filename(self) -> str:
        return self.name

    def get_outputs(self) -> T.List[str]:
        if isinstance(self.name, str):
            return [self.name]
        elif isinstance(self.name, list):
            return self.name
        else:
            raise RuntimeError('RunTarget: self.name is neither a list nor a string. This is a bug')

    def type_suffix(self) -> str:
        return "@run"

class AliasTarget(RunTarget):

    typename = 'alias'

    def __init__(self, name: str, dependencies: T.Sequence[Target],
                 subdir: str, subproject: str, environment: Environment):
        super().__init__(name, [], dependencies, subdir, subproject, environment)

    def __repr__(self) -> str:
        repr_str = "<{0} {1}>"
        return repr_str.format(self.__class__.__name__, self.get_id())

class Jar(BuildTarget):
    known_kwargs = known_jar_kwargs

    typename = 'jar'
    rust_crate_type = ''  # type: ignore[assignment]

    def __init__(self, name: str, subdir: str, subproject: str, for_machine: MachineChoice,
                 sources: T.List[SourceOutputs], structured_sources: T.Optional['StructuredSources'],
                 objects, environment: Environment, compilers: CompilerDict,
                 kwargs):
        super().__init__(name, subdir, subproject, for_machine, sources, structured_sources, objects,
                         environment, compilers, kwargs)
        for s in self.sources:
            if not s.endswith('.java'):
                raise InvalidArguments(f'Jar source {s} is not a java file.')
        for t in self.link_targets:
            if not isinstance(t, Jar):
                raise InvalidArguments(f'Link target {t} is not a jar target.')
        if self.structured_sources:
            raise InvalidArguments('structured sources are not supported in Java targets.')
        self.filename = self.name + '.jar'
        self.outputs = [self.filename]
        self.java_args = self.extra_args['java']
        self.main_class = kwargs.get('main_class', '')
        self.java_resources: T.Optional[StructuredSources] = kwargs.get('java_resources', None)

    def get_main_class(self) -> str:
        return self.main_class

    def type_suffix(self) -> str:
        return "@jar"

    def get_java_args(self) -> T.List[str]:
        return self.java_args

    def get_java_resources(self) -> T.Optional[StructuredSources]:
        return self.java_resources

    def validate_install(self) -> None:
        # All jar targets are installable.
        pass

    def is_linkable_target(self) -> bool:
        return True

    def get_classpath_args(self) -> T.List[str]:
        cp_paths = [os.path.join(l.get_subdir(), l.get_filename()) for l in self.link_targets]
        cp_string = os.pathsep.join(cp_paths)
        if cp_string:
            return ['-cp', os.pathsep.join(cp_paths)]
        return []

    def get_default_install_dir(self) -> T.Tuple[str, str]:
        return self.environment.get_jar_dir(), '{jardir}'

@dataclass(eq=False)
class CustomTargetIndex(CustomTargetBase, HoldableObject):

    """A special opaque object returned by indexing a CustomTarget. This object
    exists in Meson, but acts as a proxy in the backends, making targets depend
    on the CustomTarget it's derived from, but only adding one source file to
    the sources.
    """

    typename: T.ClassVar[str] = 'custom'

    target: T.Union[CustomTarget, CompileTarget]
    output: str

    def __post_init__(self) -> None:
        self.for_machine = self.target.for_machine

    @lazy_property
    def __index(self) -> int:
        # Must be detected lazily becuase preporcessed sources don't end up in the
        # outputs, but can still be in a CustomTargetIndex.
        return self.target.outputs.index(self.output)

    @property
    def name(self) -> str:
        return f'{self.target.name}[{self.output}]'

    @property
    def depend_files(self) -> T.List[File]:
        return self.target.depend_files

    @property
    def subdir(self) -> str:
        return self.target.subdir

    @property
    def has_custom_install_dir(self) -> bool:
        return self.target.has_custom_install_dir

    @property
    def install_dir(self) -> T.List[T.Union[str, Literal[False]]]:
        return [self.target.install_dir[self.__index]]

    def __repr__(self):
        return '<CustomTargetIndex: {!r}[{}]>'.format(self.target, self.output)

    def get_outputs(self) -> T.List[str]:
        return [self.output]

    def get_subdir(self) -> str:
        return self.target.get_subdir()

    def get_build_subdir(self) -> str:
        return self.target.get_build_subdir()

    def get_builddir(self) -> str:
        return self.target.get_builddir()

    def get_filename(self) -> str:
        return self.output

    def get_id(self) -> str:
        return self.target.get_id()

    def get_all_link_deps(self) -> ImmutableListProtocol[BuildTargetTypes]:
        return self.target.get_all_link_deps()

    def get_link_deps_mapping(self, prefix: str) -> T.Mapping[str, str]:
        return self.target.get_link_deps_mapping(prefix)

    def get_link_dep_subdirs(self) -> T.AbstractSet[str]:
        return self.target.get_link_dep_subdirs()

    def is_linkable_target(self) -> bool:
        return self.target.is_linkable_output(self.output)

    def links_dynamically(self) -> bool:
        """Whether this target links dynamically or statically

        Does not assert the target is linkable, just that it is not shared

        :return: True if is dynamically linked, otherwise False
        """
        suf = os.path.splitext(self.output)[-1]
        return suf not in {'.a', '.lib'}

    def should_install(self) -> bool:
        return self.target.should_install()

    def is_internal(self) -> bool:
        '''
        Returns True if this is a not installed static library
        '''
        suf = os.path.splitext(self.output)[-1]
        return suf in {'.a', '.lib'} and not self.should_install()

    def extract_all_objects(self) -> T.List[T.Union[str, 'ExtractedObjects']]:
        return self.target.extract_all_objects()

    def get_custom_install_dir(self) -> T.List[T.Union[str, Literal[False]]]:
        return self.target.get_custom_install_dir()

    def get_basename(self) -> str:
        return self.target.get_basename()

    def install_dir_names(self) -> T.List[T.Optional[str]]:
        install_dir_names = self.target.install_dir_names()
        return [install_dir_names[self.__index]]


class ConfigurationData(HoldableObject):
    def __init__(self, initial_values: T.Optional[T.Union[
                T.Dict[str, T.Tuple[T.Union[str, int, bool], T.Optional[str]]],
                T.Dict[str, T.Union[str, int, bool]]]
            ] = None):
        super().__init__()
        self.values: T.Dict[str, T.Tuple[T.Union[str, int, bool], T.Optional[str]]] = \
            {k: v if isinstance(v, tuple) else (v, None) for k, v in initial_values.items()} if initial_values else {}
        self.used: bool = False

    def __repr__(self) -> str:
        return repr(self.values)

    def __contains__(self, value: str) -> bool:
        return value in self.values

    def __bool__(self) -> bool:
        return bool(self.values)

    def get(self, name: str) -> T.Tuple[T.Union[str, int, bool], T.Optional[str]]:
        return self.values[name] # (val, desc)

    def keys(self) -> T.Iterator[str]:
        return self.values.keys()

class LocalProgram(programs.Program):
    def __init__(self, program: T.Union[programs.ExternalProgram, Executable, CustomTarget, CustomTargetIndex], version: str) -> None:
        super().__init__()
        if isinstance(program, CustomTarget):
            if len(program.outputs) != 1:
                raise InvalidArguments('CustomTarget used as LocalProgram must have exactly one output.')
        self.name = program.name
        self.for_machine = program.for_machine
        self.program = program
        self.version = version

    def found(self) -> bool:
        return True

    def get_version(self, interpreter: T.Optional[Interpreter] = None) -> str:
        return self.version

    def get_command(self) -> T.List[str]:
        if isinstance(self.program, programs.ExternalProgram):
            return self.program.get_command()
        # Only the backend knows the actual path to the build program.
        raise MesonBugException('Cannot call get_command() on program that is a build target.')

    def get_path(self) -> str:
        if isinstance(self.program, programs.ExternalProgram):
            return self.program.get_path()
        # Only the backend knows the actual path to the build program.
        raise MesonBugException('Cannot call get_path() on program that is a build target.')

    def description(self) -> str:
        if isinstance(self.program, programs.ExternalProgram):
            return self.program.description()
        if isinstance(self.program, Executable):
            return self.program.name
        return self.program.get_filename()

    def runnable(self) -> bool:
        return isinstance(self.program, programs.ExternalProgram)

# A bit poorly named, but this represents plain data files to copy
# during install.
@dataclass(eq=False)
class Data(HoldableObject):
    sources: T.List[File]
    install_dir: str
    install_dir_name: str
    install_mode: 'FileMode'
    subproject: str
    rename: T.List[str] = None
    install_tag: T.Optional[str] = None
    data_type: str = None
    follow_symlinks: T.Optional[bool] = None

    def __post_init__(self) -> None:
        if self.rename is None:
            self.rename = [os.path.basename(f.fname) for f in self.sources]

@dataclass(eq=False)
class SymlinkData(HoldableObject):
    target: str
    name: str
    install_dir: str
    subproject: str
    install_tag: T.Optional[str] = None

    def __post_init__(self) -> None:
        if self.name != os.path.basename(self.name):
            raise InvalidArguments(f'Link name is "{self.name}", but link names cannot contain path separators. '
                                   'The dir part should be in install_dir.')

@dataclass(eq=False)
class TestSetup:
    exe_wrapper: T.List[str]
    gdb: bool
    timeout_multiplier: int
    env: EnvironmentVariables
    exclude_suites: T.List[str]

def get_sources_string_names(sources, backend):
    '''
    For the specified list of @sources which can be strings, Files, or targets,
    get all the output basenames.
    '''
    names = []
    for s in sources:
        if isinstance(s, str):
            names.append(s)
        elif isinstance(s, (BuildTarget, CustomTarget, CustomTargetIndex, GeneratedList)):
            names += s.get_outputs()
        elif isinstance(s, ExtractedObjects):
            names += backend.determine_ext_objs(s)
        elif isinstance(s, File):
            names.append(s.fname)
        else:
            raise AssertionError(f'Unknown source type: {s!r}')
    return names

def load(build_dir: str) -> Build:
    filename = os.path.join(build_dir, 'meson-private', 'build.dat')
    try:
        b = pickle_load(filename, 'Build data', Build)
        # We excluded coredata when saving Build object, load it separately
        b.environment.coredata = coredata.load(build_dir)
        return b
    except FileNotFoundError:
        raise MesonException(f'No such build data file as {filename!r}.')


def save(obj: Build, filename: str) -> None:
    # Exclude coredata because we pickle it separately already
    cdata = obj.environment.coredata
    obj.environment.coredata = None
    try:
        with open(filename, 'wb') as f:
            pickle.dump(obj, f)
    finally:
        obj.environment.coredata = cdata
