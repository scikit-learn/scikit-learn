# SPDX-License-Identifier: Apache-2.0
# Copyright 2013-2018 The Meson development team
# Copyright © 2024 Intel Corporation

# This file contains the detection logic for external dependencies.
# Custom logic for several other packages are in separate files.

from __future__ import annotations
import copy
import os
import collections
import itertools
import typing as T
import uuid
from enum import Enum

from .. import mlog, mesonlib
from ..compilers import clib_langs
from ..mesonlib import LibType, MachineChoice, MesonException, HoldableObject, version_compare_many
from ..options import OptionKey
#from ..interpreterbase import FeatureDeprecated, FeatureNew

if T.TYPE_CHECKING:
    from ..compilers.compilers import Compiler
    from ..environment import Environment
    from ..interpreterbase import FeatureCheckBase
    from ..build import (
        CustomTarget, IncludeDirs, CustomTargetIndex, LibTypes,
        StaticLibrary, StructuredSources, ExtractedObjects, GeneratedTypes
    )
    from ..interpreter.type_checking import PkgConfigDefineType

    _MissingCompilerBase = Compiler
else:
    _MissingCompilerBase = object


class DependencyException(MesonException):
    '''Exceptions raised while trying to find dependencies'''


class MissingCompiler(_MissingCompilerBase):
    """Represent a None Compiler - when no tool chain is found.
    replacing AttributeError with DependencyException"""

    # These are needed in type checking mode to avoid errors, but we don't want
    # the extra overhead at runtime
    if T.TYPE_CHECKING:
        def __init__(self) -> None:
            pass

        def get_optimization_args(self, optimization_level: str) -> T.List[str]:
            return []

        def get_output_args(self, outputname: str) -> T.List[str]:
            return []

        def sanity_check(self, work_dir: str, environment: 'Environment') -> None:
            return None

    def __getattr__(self, item: str) -> T.Any:
        if item.startswith('__'):
            raise AttributeError()
        raise DependencyException('no toolchain found')

    def __bool__(self) -> bool:
        return False


class DependencyMethods(Enum):
    # Auto means to use whatever dependency checking mechanisms in whatever order meson thinks is best.
    AUTO = 'auto'
    PKGCONFIG = 'pkg-config'
    CMAKE = 'cmake'
    # The dependency is provided by the standard library and does not need to be linked
    BUILTIN = 'builtin'
    # Just specify the standard link arguments, assuming the operating system provides the library.
    SYSTEM = 'system'
    # This is only supported on OSX - search the frameworks directory by name.
    EXTRAFRAMEWORK = 'extraframework'
    # Detect using the sysconfig module.
    SYSCONFIG = 'sysconfig'
    # Specify using a "program"-config style tool
    CONFIG_TOOL = 'config-tool'
    # For backwards compatibility
    SDLCONFIG = 'sdlconfig'
    CUPSCONFIG = 'cups-config'
    PCAPCONFIG = 'pcap-config'
    LIBWMFCONFIG = 'libwmf-config'
    QMAKE = 'qmake'
    # Misc
    DUB = 'dub'


DependencyTypeName = T.NewType('DependencyTypeName', str)


class Dependency(HoldableObject):

    @classmethod
    def _process_include_type_kw(cls, kwargs: T.Dict[str, T.Any]) -> str:
        if 'include_type' not in kwargs:
            return 'preserve'
        if not isinstance(kwargs['include_type'], str):
            raise DependencyException('The include_type kwarg must be a string type')
        if kwargs['include_type'] not in ['preserve', 'system', 'non-system']:
            raise DependencyException("include_type may only be one of ['preserve', 'system', 'non-system']")
        return kwargs['include_type']

    def __init__(self, type_name: DependencyTypeName, kwargs: T.Dict[str, T.Any]) -> None:
        # This allows two Dependencies to be compared even after being copied.
        # The purpose is to allow the name to be changed, but still have a proper comparison
        self._id = uuid.uuid4().int
        self.name = f'dep{id(self)}'
        self.version:  T.Optional[str] = None
        self.language: T.Optional[str] = None # None means C-like
        self.is_found = False
        self.type_name = type_name
        self.compile_args: T.List[str] = []
        self.link_args:    T.List[str] = []
        # Raw -L and -l arguments without manual library searching
        # If None, self.link_args will be used
        self.raw_link_args: T.Optional[T.List[str]] = None
        self.sources: T.List[T.Union[mesonlib.File, GeneratedTypes, 'StructuredSources']] = []
        self.extra_files: T.List[mesonlib.File] = []
        self.include_type = self._process_include_type_kw(kwargs)
        self.ext_deps: T.List[Dependency] = []
        self.d_features: T.DefaultDict[str, T.List[T.Any]] = collections.defaultdict(list)
        self.featurechecks: T.List['FeatureCheckBase'] = []
        self.feature_since: T.Optional[T.Tuple[str, str]] = None

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Dependency):
            return NotImplemented
        return self._id == other._id

    def __hash__(self) -> int:
        return self._id

    def __repr__(self) -> str:
        return f'<{self.__class__.__name__} {self.name}: {self.is_found}>'

    def is_built(self) -> bool:
        return False

    def summary_value(self) -> T.Union[str, mlog.AnsiDecorator, mlog.AnsiText]:
        if not self.found():
            return mlog.red('NO')
        if not self.version:
            return mlog.green('YES')
        return mlog.AnsiText(mlog.green('YES'), ' ', mlog.cyan(self.version))

    def get_compile_args(self) -> T.List[str]:
        if self.include_type == 'system':
            converted = []
            for i in self.compile_args:
                if i.startswith('-I') or i.startswith('/I'):
                    converted += ['-isystem' + i[2:]]
                else:
                    converted += [i]
            return converted
        if self.include_type == 'non-system':
            converted = []
            for i in self.compile_args:
                if i.startswith('-isystem'):
                    converted += ['-I' + i[8:]]
                else:
                    converted += [i]
            return converted
        return self.compile_args

    def get_all_compile_args(self) -> T.List[str]:
        """Get the compile arguments from this dependency and its sub dependencies."""
        return list(itertools.chain(self.get_compile_args(),
                                    *(d.get_all_compile_args() for d in self.ext_deps)))

    def get_link_args(self, language: T.Optional[str] = None, raw: bool = False) -> T.List[str]:
        if raw and self.raw_link_args is not None:
            return self.raw_link_args
        return self.link_args

    def get_all_link_args(self) -> T.List[str]:
        """Get the link arguments from this dependency and its sub dependencies."""
        return list(itertools.chain(self.get_link_args(),
                                    *(d.get_all_link_args() for d in self.ext_deps)))

    def found(self) -> bool:
        return self.is_found

    def get_sources(self) -> T.List[T.Union[mesonlib.File, GeneratedTypes, 'StructuredSources']]:
        """Source files that need to be added to the target.
        As an example, gtest-all.cc when using GTest."""
        return self.sources

    def get_extra_files(self) -> T.List[mesonlib.File]:
        """Mostly for introspection and IDEs"""
        return self.extra_files

    def get_name(self) -> str:
        return self.name

    def get_version(self) -> str:
        if self.version:
            return self.version
        else:
            return 'unknown'

    def get_include_dirs(self) -> T.List['IncludeDirs']:
        return []

    def get_include_type(self) -> str:
        return self.include_type

    def get_exe_args(self, compiler: 'Compiler') -> T.List[str]:
        return []

    def get_partial_dependency(self, *, compile_args: bool = False,
                               link_args: bool = False, links: bool = False,
                               includes: bool = False, sources: bool = False) -> 'Dependency':
        """Create a new dependency that contains part of the parent dependency.

        The following options can be inherited:
            links -- all link_with arguments
            includes -- all include_directory and -I/-isystem calls
            sources -- any source, header, or generated sources
            compile_args -- any compile args
            link_args -- any link args

        Additionally the new dependency will have the version parameter of its
        parent (if any) and the requested values of any dependencies will be
        added as well.
        """
        raise RuntimeError('Unreachable code in partial_dependency called')

    def _add_sub_dependency(self, deplist: T.Iterable[T.Callable[[], 'Dependency']]) -> bool:
        """Add an internal dependency from a list of possible dependencies.

        This method is intended to make it easier to add additional
        dependencies to another dependency internally.

        Returns true if the dependency was successfully added, false
        otherwise.
        """
        for d in deplist:
            dep = d()
            if dep.is_found:
                self.ext_deps.append(dep)
                return True
        return False

    def get_variable(self, *, cmake: T.Optional[str] = None, pkgconfig: T.Optional[str] = None,
                     configtool: T.Optional[str] = None, internal: T.Optional[str] = None,
                     system: T.Optional[str] = None, default_value: T.Optional[str] = None,
                     pkgconfig_define: PkgConfigDefineType = None) -> str:
        if default_value is not None:
            return default_value
        raise DependencyException(f'No default provided for dependency {self!r}, which is not pkg-config, cmake, or config-tool based.')

    def generate_system_dependency(self, include_type: str) -> 'Dependency':
        new_dep = copy.deepcopy(self)
        new_dep.include_type = self._process_include_type_kw({'include_type': include_type})
        return new_dep

    def get_as_static(self, recursive: bool) -> Dependency:
        """Used as base case for internal_dependency"""
        return self

    def get_as_shared(self, recursive: bool) -> Dependency:
        """Used as base case for internal_dependency"""
        return self

class InternalDependency(Dependency):
    def __init__(self, version: str, incdirs: T.List['IncludeDirs'], compile_args: T.List[str],
                 link_args: T.List[str],
                 libraries: T.List[LibTypes],
                 whole_libraries: T.List[T.Union[StaticLibrary, CustomTarget, CustomTargetIndex]],
                 sources: T.Sequence[T.Union[mesonlib.File, GeneratedTypes, StructuredSources]],
                 extra_files: T.Sequence[mesonlib.File],
                 ext_deps: T.List[Dependency], variables: T.Dict[str, str],
                 d_module_versions: T.List[T.Union[str, int]], d_import_dirs: T.List['IncludeDirs'],
                 objects: T.List['ExtractedObjects']):
        super().__init__(DependencyTypeName('internal'), {})
        self.version = version
        self.is_found = True
        self.include_directories = incdirs
        self.compile_args = compile_args
        self.link_args = link_args
        self.libraries = libraries
        self.whole_libraries = whole_libraries
        self.sources = list(sources)
        self.extra_files = list(extra_files)
        self.ext_deps = ext_deps
        self.variables = variables
        self.objects = objects
        if d_module_versions:
            self.d_features['versions'] = d_module_versions
        if d_import_dirs:
            self.d_features['import_dirs'] = d_import_dirs

    def __deepcopy__(self, memo: T.Dict[int, 'InternalDependency']) -> 'InternalDependency':
        result = self.__class__.__new__(self.__class__)
        assert isinstance(result, InternalDependency)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            if k in {'libraries', 'whole_libraries'}:
                setattr(result, k, copy.copy(v))
            else:
                setattr(result, k, copy.deepcopy(v, memo))
        return result

    def summary_value(self) -> mlog.AnsiDecorator:
        # Omit the version.  Most of the time it will be just the project
        # version, which is uninteresting in the summary.
        return mlog.green('YES')

    def is_built(self) -> bool:
        if self.sources or self.libraries or self.whole_libraries:
            return True
        return any(d.is_built() for d in self.ext_deps)

    def get_partial_dependency(self, *, compile_args: bool = False,
                               link_args: bool = False, links: bool = False,
                               includes: bool = False, sources: bool = False,
                               extra_files: bool = False) -> InternalDependency:
        final_compile_args = self.compile_args.copy() if compile_args else []
        final_link_args = self.link_args.copy() if link_args else []
        final_libraries = self.libraries.copy() if links else []
        final_whole_libraries = self.whole_libraries.copy() if links else []
        final_sources = self.sources.copy() if sources else []
        final_extra_files = self.extra_files.copy() if extra_files else []
        final_includes = self.include_directories.copy() if includes else []
        final_deps = [d.get_partial_dependency(
            compile_args=compile_args, link_args=link_args, links=links,
            includes=includes, sources=sources) for d in self.ext_deps]
        return InternalDependency(
            self.version, final_includes, final_compile_args,
            final_link_args, final_libraries, final_whole_libraries,
            final_sources, final_extra_files, final_deps, self.variables, [], [], [])

    def get_include_dirs(self) -> T.List['IncludeDirs']:
        return self.include_directories

    def get_variable(self, *, cmake: T.Optional[str] = None, pkgconfig: T.Optional[str] = None,
                     configtool: T.Optional[str] = None, internal: T.Optional[str] = None,
                     system: T.Optional[str] = None, default_value: T.Optional[str] = None,
                     pkgconfig_define: PkgConfigDefineType = None) -> str:
        val = self.variables.get(internal, default_value)
        if val is not None:
            return val
        raise DependencyException(f'Could not get an internal variable and no default provided for {self!r}')

    def generate_link_whole_dependency(self) -> Dependency:
        from ..build import SharedLibrary, CustomTarget, CustomTargetIndex
        new_dep = copy.deepcopy(self)
        for x in new_dep.libraries:
            if isinstance(x, SharedLibrary):
                raise MesonException('Cannot convert a dependency to link_whole when it contains a '
                                     'SharedLibrary')
            elif isinstance(x, (CustomTarget, CustomTargetIndex)) and x.links_dynamically():
                raise MesonException('Cannot convert a dependency to link_whole when it contains a '
                                     'CustomTarget or CustomTargetIndex which is a shared library')

        # Mypy doesn't understand that the above is a TypeGuard
        new_dep.whole_libraries += T.cast('T.List[T.Union[StaticLibrary, CustomTarget, CustomTargetIndex]]',
                                          new_dep.libraries)
        new_dep.libraries = []
        return new_dep

    def get_as_static(self, recursive: bool) -> InternalDependency:
        new_dep = copy.copy(self)
        new_dep.libraries = [lib.get('static') for lib in self.libraries]
        if recursive:
            new_dep.ext_deps = [dep.get_as_static(True) for dep in self.ext_deps]
        return new_dep

    def get_as_shared(self, recursive: bool) -> InternalDependency:
        new_dep = copy.copy(self)
        new_dep.libraries = [lib.get('shared') for lib in self.libraries]
        if recursive:
            new_dep.ext_deps = [dep.get_as_shared(True) for dep in self.ext_deps]
        return new_dep

class HasNativeKwarg:
    def __init__(self, kwargs: T.Dict[str, T.Any]):
        self.for_machine = self.get_for_machine_from_kwargs(kwargs)

    def get_for_machine_from_kwargs(self, kwargs: T.Dict[str, T.Any]) -> MachineChoice:
        return MachineChoice.BUILD if kwargs.get('native', False) else MachineChoice.HOST

class ExternalDependency(Dependency, HasNativeKwarg):
    def __init__(self, type_name: DependencyTypeName, environment: 'Environment', kwargs: T.Dict[str, T.Any], language: T.Optional[str] = None):
        Dependency.__init__(self, type_name, kwargs)
        self.env = environment
        self.name = type_name # default
        self.is_found = False
        self.language = language
        version_reqs = kwargs.get('version', None)
        if isinstance(version_reqs, str):
            version_reqs = [version_reqs]
        self.version_reqs: T.Optional[T.List[str]] = version_reqs
        self.required = kwargs.get('required', True)
        self.silent = kwargs.get('silent', False)
        self.static = kwargs.get('static', self.env.coredata.get_option(OptionKey('prefer_static')))
        self.libtype = LibType.STATIC if self.static else LibType.PREFER_SHARED
        if not isinstance(self.static, bool):
            raise DependencyException('Static keyword must be boolean')
        # Is this dependency to be run on the build platform?
        HasNativeKwarg.__init__(self, kwargs)
        self.clib_compiler = detect_compiler(self.name, environment, self.for_machine, self.language)

    def get_compiler(self) -> T.Union['MissingCompiler', 'Compiler']:
        return self.clib_compiler

    def get_partial_dependency(self, *, compile_args: bool = False,
                               link_args: bool = False, links: bool = False,
                               includes: bool = False, sources: bool = False) -> Dependency:
        new = copy.copy(self)
        new._id = uuid.uuid4().int
        if not compile_args:
            new.compile_args = []
        if not link_args:
            new.link_args = []
        if not sources:
            new.sources = []
        if not includes:
            pass # TODO maybe filter compile_args?
        if not sources:
            new.sources = []

        return new

    def log_details(self) -> str:
        return ''

    def log_info(self) -> str:
        return ''

    @staticmethod
    def log_tried() -> str:
        return ''

    # Check if dependency version meets the requirements
    def _check_version(self) -> None:
        if not self.is_found:
            return

        if self.version_reqs:
            # an unknown version can never satisfy any requirement
            if not self.version:
                self.is_found = False
                found_msg: mlog.TV_LoggableList = []
                found_msg += ['Dependency', mlog.bold(self.name), 'found:']
                found_msg += [str(mlog.red('NO')) + '.', 'Unknown version, but need:', self.version_reqs]
                mlog.log(*found_msg)

                if self.required:
                    m = f'Unknown version, but need {self.version_reqs!r}.'
                    raise DependencyException(m)

            else:
                (self.is_found, not_found, found) = \
                    version_compare_many(self.version, self.version_reqs)
                if not self.is_found:
                    found_msg = ['Dependency', mlog.bold(self.name), 'found:']
                    found_msg += [str(mlog.red('NO')) + '.',
                                  'Found', mlog.normal_cyan(self.version), 'but need:',
                                  mlog.bold(', '.join([f"'{e}'" for e in not_found]))]
                    if found:
                        found_msg += ['; matched:',
                                      ', '.join([f"'{e}'" for e in found])]
                    mlog.log(*found_msg)

                    if self.required:
                        m = 'Invalid version, need {!r} {!r} found {!r}.'
                        raise DependencyException(m.format(self.name, not_found, self.version))
                    return


class NotFoundDependency(Dependency):
    def __init__(self, name: str, environment: 'Environment') -> None:
        super().__init__(DependencyTypeName('not-found'), {})
        self.env = environment
        self.name = name
        self.is_found = False

    def get_partial_dependency(self, *, compile_args: bool = False,
                               link_args: bool = False, links: bool = False,
                               includes: bool = False, sources: bool = False) -> 'NotFoundDependency':
        new = copy.copy(self)
        new._id = uuid.uuid4().int
        return new


class ExternalLibrary(ExternalDependency):
    def __init__(self, name: str, link_args: T.List[str], environment: 'Environment',
                 language: str, silent: bool = False) -> None:
        super().__init__(DependencyTypeName('library'), environment, {}, language=language)
        self.name = name
        self.language = language
        self.is_found = False
        if link_args:
            self.is_found = True
            self.link_args = link_args
        if not silent:
            if self.is_found:
                mlog.log('Library', mlog.bold(name), 'found:', mlog.green('YES'))
            else:
                mlog.log('Library', mlog.bold(name), 'found:', mlog.red('NO'))

    def get_link_args(self, language: T.Optional[str] = None, raw: bool = False) -> T.List[str]:
        '''
        External libraries detected using a compiler must only be used with
        compatible code. For instance, Vala libraries (.vapi files) cannot be
        used with C code, and not all Rust library types can be linked with
        C-like code. Note that C++ libraries *can* be linked with C code with
        a C++ linker (and vice-versa).
        '''
        # Using a vala library in a non-vala target, or a non-vala library in a vala target
        # XXX: This should be extended to other non-C linkers such as Rust
        if (self.language == 'vala' and language != 'vala') or \
           (language == 'vala' and self.language != 'vala'):
            return []
        return super().get_link_args(language=language, raw=raw)

    def get_partial_dependency(self, *, compile_args: bool = False,
                               link_args: bool = False, links: bool = False,
                               includes: bool = False, sources: bool = False) -> 'ExternalLibrary':
        # External library only has link_args, so ignore the rest of the
        # interface.
        new = copy.copy(self)
        new._id = uuid.uuid4().int
        if not link_args:
            new.link_args = []
        return new


def get_leaf_external_dependencies(deps: T.List[Dependency]) -> T.List[Dependency]:
    if not deps:
        # Ensure that we always return a new instance
        return deps.copy()
    final_deps = []
    while deps:
        next_deps = []
        for d in mesonlib.listify(deps):
            if not isinstance(d, Dependency) or d.is_built():
                raise DependencyException('Dependencies must be external dependencies')
            final_deps.append(d)
            next_deps.extend(d.ext_deps)
        deps = next_deps
    return final_deps


def sort_libpaths(libpaths: T.List[str], refpaths: T.List[str]) -> T.List[str]:
    """Sort <libpaths> according to <refpaths>

    It is intended to be used to sort -L flags returned by pkg-config.
    Pkg-config returns flags in random order which cannot be relied on.
    """
    if len(refpaths) == 0:
        return list(libpaths)

    def key_func(libpath: str) -> T.Tuple[int, int]:
        common_lengths: T.List[int] = []
        for refpath in refpaths:
            try:
                common_path: str = os.path.commonpath([libpath, refpath])
            except ValueError:
                common_path = ''
            common_lengths.append(len(common_path))
        max_length = max(common_lengths)
        max_index = common_lengths.index(max_length)
        reversed_max_length = len(refpaths[max_index]) - max_length
        return (max_index, reversed_max_length)
    return sorted(libpaths, key=key_func)

def strip_system_libdirs(environment: 'Environment', for_machine: MachineChoice, link_args: T.List[str]) -> T.List[str]:
    """Remove -L<system path> arguments.

    leaving these in will break builds where a user has a version of a library
    in the system path, and a different version not in the system path if they
    want to link against the non-system path version.
    """
    exclude = {f'-L{p}' for p in environment.get_compiler_system_lib_dirs(for_machine)}
    return [l for l in link_args if l not in exclude]

def strip_system_includedirs(environment: 'Environment', for_machine: MachineChoice, include_args: T.List[str]) -> T.List[str]:
    """Remove -I<system path> arguments.

    leaving these in will break builds where user want dependencies with system
    include-type used in rust.bindgen targets as if will cause system headers
    to not be found.
    """

    exclude = {f'-I{p}' for p in environment.get_compiler_system_include_dirs(for_machine)}
    return [i for i in include_args if i not in exclude]

def process_method_kw(possible: T.Iterable[DependencyMethods], kwargs: T.Dict[str, T.Any]) -> T.List[DependencyMethods]:
    method: T.Union[DependencyMethods, str] = kwargs.get('method', 'auto')
    if isinstance(method, DependencyMethods):
        return [method]
    # TODO: try/except?
    if method not in [e.value for e in DependencyMethods]:
        raise DependencyException(f'method {method!r} is invalid')
    method = DependencyMethods(method)

    # Raise FeatureNew where appropriate
    if method is DependencyMethods.CONFIG_TOOL:
        # FIXME: needs to get a handle on the subproject
        # FeatureNew.single_use('Configuration method "config-tool"', '0.44.0')
        pass
    # This sets per-tool config methods which are deprecated to to the new
    # generic CONFIG_TOOL value.
    if method in [DependencyMethods.SDLCONFIG, DependencyMethods.CUPSCONFIG,
                  DependencyMethods.PCAPCONFIG, DependencyMethods.LIBWMFCONFIG]:
        # FIXME: needs to get a handle on the subproject
        #FeatureDeprecated.single_use(f'Configuration method {method.value}', '0.44', 'Use "config-tool" instead.')
        method = DependencyMethods.CONFIG_TOOL
    if method is DependencyMethods.QMAKE:
        # FIXME: needs to get a handle on the subproject
        # FeatureDeprecated.single_use('Configuration method "qmake"', '0.58', 'Use "config-tool" instead.')
        method = DependencyMethods.CONFIG_TOOL

    # Set the detection method. If the method is set to auto, use any available method.
    # If method is set to a specific string, allow only that detection method.
    if method == DependencyMethods.AUTO:
        methods = list(possible)
    elif method in possible:
        methods = [method]
    else:
        raise DependencyException(
            'Unsupported detection method: {}, allowed methods are {}'.format(
                method.value,
                mlog.format_list([x.value for x in [DependencyMethods.AUTO] + list(possible)])))

    return methods

def detect_compiler(name: str, env: 'Environment', for_machine: MachineChoice,
                    language: T.Optional[str]) -> T.Union['MissingCompiler', 'Compiler']:
    """Given a language and environment find the compiler used."""
    compilers = env.coredata.compilers[for_machine]

    # Set the compiler for this dependency if a language is specified,
    # else try to pick something that looks usable.
    if language:
        if language not in compilers:
            m = name.capitalize() + ' requires a {0} compiler, but ' \
                '{0} is not in the list of project languages'
            raise DependencyException(m.format(language.capitalize()))
        return compilers[language]
    else:
        for lang in clib_langs:
            try:
                return compilers[lang]
            except KeyError:
                continue
    return MissingCompiler()


class SystemDependency(ExternalDependency):

    """Dependency base for System type dependencies."""

    def __init__(self, name: str, env: 'Environment', kwargs: T.Dict[str, T.Any],
                 language: T.Optional[str] = None) -> None:
        super().__init__(DependencyTypeName('system'), env, kwargs, language=language)
        self.name = name

    @staticmethod
    def log_tried() -> str:
        return 'system'


class BuiltinDependency(ExternalDependency):

    """Dependency base for Builtin type dependencies."""

    def __init__(self, name: str, env: 'Environment', kwargs: T.Dict[str, T.Any],
                 language: T.Optional[str] = None) -> None:
        super().__init__(DependencyTypeName('builtin'), env, kwargs, language=language)
        self.name = name

    @staticmethod
    def log_tried() -> str:
        return 'builtin'
