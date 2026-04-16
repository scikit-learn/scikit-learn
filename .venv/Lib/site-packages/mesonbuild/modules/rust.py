# SPDX-License-Identifier: Apache-2.0
# Copyright © 2020-2025 Intel Corporation

from __future__ import annotations
import itertools
import os
import re
import textwrap
import typing as T

from mesonbuild.interpreterbase.decorators import FeatureNew

from . import ExtensionModule, ModuleReturnValue, ModuleInfo, ModuleObject
from .. import mesonlib, mlog
from ..build import (BothLibraries, BuildTarget, CustomTargetIndex, Executable, ExtractedObjects, GeneratedList,
                     CustomTarget, InvalidArguments, Jar, StructuredSources, SharedLibrary, StaticLibrary,
                     SharedModule)
from ..compilers.compilers import are_asserts_disabled_for_subproject, lang_suffixes
from ..compilers.rust import parse_target, RustSystemDependency
from ..dependencies import Dependency
from ..interpreter.type_checking import (
    DEPENDENCIES_KW, LINK_WITH_KW, LINK_WHOLE_KW, SHARED_LIB_KWS, TEST_KWS, TEST_KWS_NO_ARGS,
    OUTPUT_KW, INCLUDE_DIRECTORIES, SOURCES_VARARGS, NATIVE_KW, NoneType, in_set_validator,
    EXECUTABLE_KWS, LIBRARY_KWS, SHARED_MOD_KWS, _BASE_LANG_KW
)
from ..interpreterbase import ContainerTypeInfo, InterpreterException, KwargInfo, typed_kwargs, typed_pos_args, noKwargs, noPosargs, permittedKwargs
from ..interpreter.interpreterobjects import Doctest
from ..mesonlib import (is_parent_path, File, MachineChoice, MesonException, PerMachine)
from ..programs import ExternalProgram, NonExistingExternalProgram

if T.TYPE_CHECKING:
    from . import ModuleState
    from .. import cargo
    from ..build import BuildTargetTypes, ExecutableKeywordArguments, IncludeDirs, LibTypes
    from ..cargo.interpreter import RUST_ABI
    from ..compilers.compilers import Language
    from ..compilers.rust import RustCompiler
    from ..dependencies import ExternalLibrary
    from ..interpreter import Interpreter
    from ..interpreter import kwargs as _kwargs
    from ..interpreter.interpreter import SourceInputs, SourceOutputs
    from ..interpreter.interpreterobjects import Test
    from ..interpreterbase import TYPE_kwargs
    from ..programs import Program
    from ..interpreter.type_checking import SourcesVarargsType

    from typing_extensions import TypedDict, Literal

    ArgsType = T.TypeVar('ArgsType')

    class FuncRustTest(_kwargs.BaseTest, T.Generic[ArgsType]):
        args: T.List[ArgsType]
        dependencies: T.List[T.Union[Dependency, ExternalLibrary]]
        is_parallel: bool
        link_with: T.List[LibTypes]
        link_whole: T.List[T.Union[StaticLibrary, CustomTarget, CustomTargetIndex]]
        rust_args: T.List[str]

    FuncTest = FuncRustTest[_kwargs.TestArgs]
    FuncDoctest = FuncRustTest[str]

    class FuncBindgen(TypedDict):

        args: T.List[str]
        c_args: T.List[str]
        include_directories: T.List[IncludeDirs]
        input: T.List[SourceInputs]
        output: str
        output_inline_wrapper: str
        dependencies: T.List[T.Union[Dependency, ExternalLibrary]]
        language: T.Optional[Literal['c', 'cpp']]
        bindgen_version: T.List[str]

    class FuncWorkspace(TypedDict):
        default_features: T.Optional[bool]
        features: T.List[str]

    class FuncDependency(TypedDict):
        rust_abi: T.Optional[RUST_ABI]

    class RustPackageDependencies(TypedDict):
        dependencies: bool
        dev_dependencies: bool
        system_dependencies: bool

    class RustPackageExecutable(_kwargs.Executable):
        pass

    class RustPackageLibrary(_kwargs.Library):
        pass

RUST_TEST_KWS: T.List[KwargInfo] = [
     KwargInfo(
         'rust_args',
         ContainerTypeInfo(list, str),
         listify=True,
         default=[],
         since='1.2.0',
     ),
     KwargInfo('is_parallel', bool, default=False),
]

# The native argument should be passed to the "package" method
EXECUTABLE_KWS_NO_NATIVE = [kwi for kwi in EXECUTABLE_KWS if kwi.name != 'native']
LIBRARY_KWS_NO_NATIVE = [kwi for kwi in LIBRARY_KWS if kwi.name != 'native']
SHARED_LIB_KWS_NO_NATIVE = [kwi for kwi in SHARED_LIB_KWS if kwi.name != 'native']
SHARED_MOD_KWS_NO_NATIVE = [kwi for kwi in SHARED_MOD_KWS if kwi.name != 'native']


def no_spaces_validator(arg: T.Optional[T.Union[str, T.List]]) -> T.Optional[str]:
    if any(bool(re.search(r'\s', x)) for x in arg):
        return 'must not contain spaces due to limitations of rustdoc'
    return None

def dep_to_system_dependency(dep: Dependency, depname: str) -> Dependency:
    if not dep.found():
        return dep
    if not depname:
        if not dep.name:
            raise MesonException("rust.to_system_dependency() called with an unnamed dependency and no explicit name")
        depname = dep.name
    depname = re.sub(r'[^a-zA-Z0-9]', '_', depname)
    rust_args = ['--cfg', f'system_deps_have_{depname}']
    return RustSystemDependency(dep.version, compile_args=rust_args, ext_deps=[dep], name=dep.name)

class RustWorkspace(ModuleObject):
    """Represents a Rust workspace, controlling the build of packages
       recorded in a Cargo.lock file."""

    def __init__(self, interpreter: Interpreter, ws: cargo.WorkspaceState) -> None:
        super().__init__()
        self.interpreter = interpreter
        self.ws = ws
        self.methods.update({
            'packages': self.packages_method,
            'package': self.package_method,
            'subproject': self.subproject_method,
        })

    @property
    def subdir(self) -> str:
        return self.ws.subdir

    @noPosargs
    @noKwargs
    def packages_method(self, state: ModuleState, args: T.List, kwargs: TYPE_kwargs) -> T.List[str]:
        """Returns list of package names in workspace."""
        package_names = [pkg.manifest.package.name for pkg in self.ws.packages.values()]
        return sorted(package_names)

    @typed_pos_args('workspace.package', optargs=[str])
    def package_method(self, state: 'ModuleState', args: T.List, kwargs: TYPE_kwargs) -> RustPackage:
        """Returns a package object."""
        package_name = args[0] if args else None
        return RustPackage(self, self.interpreter.cargo.load_package(self.ws, package_name))

    def _do_subproject(self, pkg: cargo.PackageState) -> None:
        kw: _kwargs.DoSubproject = {
            'required': True,
            'version': None,
            'options': None,
            'cmake_options': [],
            'default_options': {},
        }
        subp_name = pkg.get_subproject_name()
        self.interpreter.do_subproject(subp_name, kw, force_method='cargo')

    @typed_pos_args('workspace.subproject', str, optargs=[str])
    @noKwargs
    def subproject_method(self, state: ModuleState, args: T.Tuple[str, T.Optional[str]], kwargs: TYPE_kwargs) -> RustSubproject:
        """Returns a package object for a subproject package."""
        package_name = args[0]
        pkg = self.interpreter.cargo.resolve_package(package_name, args[1] or '')
        if pkg is None:
            if args[1]:
                raise MesonException(f'No version of cargo package "{package_name}" provides API {args[1]}')
            else:
                raise MesonException(f'Cargo package "{package_name}" not available')

        self._do_subproject(pkg)
        return RustSubproject(self, pkg)


class RustCrate(ModuleObject):
    """Abstract base class for Rust crate representations."""

    def __init__(self, rust_ws: RustWorkspace, package: cargo.PackageState) -> None:
        super().__init__()
        self.rust_ws = rust_ws
        self.package = package
        self.methods.update({
            'all_features': self.all_features_method,
            'api': self.api_method,
            'features': self.features_method,
            'name': self.name_method,
            'version': self.version_method,
            'rust_args': self.rust_args_method,
            'env': self.env_method, # type: ignore[dict-item]
            'rust_dependency_map': self.rust_dependency_map_method, # type: ignore[dict-item]
        })

    @noPosargs
    @noKwargs
    def name_method(self, state: ModuleState, args: T.List, kwargs: TYPE_kwargs) -> str:
        """Returns the name of the package."""
        return self.package.manifest.package.name

    @noPosargs
    @noKwargs
    def api_method(self, state: ModuleState, args: T.List, kwargs: TYPE_kwargs) -> str:
        """Returns the API version of the package."""
        return self.package.manifest.package.api

    @noPosargs
    @noKwargs
    def version_method(self, state: ModuleState, args: T.List, kwargs: TYPE_kwargs) -> str:
        """Returns the version of the package."""
        return self.package.manifest.package.version

    @noPosargs
    @noKwargs
    def all_features_method(self, state: ModuleState, args: T.List, kwargs: TYPE_kwargs) -> T.List[str]:
        """Returns all features for specific package."""
        return sorted(list(self.package.manifest.features.keys()))

    @noPosargs
    @noKwargs
    def features_method(self, state: ModuleState, args: T.List, kwargs: TYPE_kwargs) -> T.List[str]:
        """Returns chosen features for specific package."""
        return sorted(list(self.package.cfg.features))

    @noPosargs
    @noKwargs
    def rust_args_method(self, state: ModuleState, args: T.List, kwargs: TYPE_kwargs) -> T.List[str]:
        """Returns rustc arguments for this package."""
        return self.package.get_rustc_args(state.environment, state.subdir, mesonlib.MachineChoice.HOST)

    @noPosargs
    @noKwargs
    def env_method(self, state: ModuleState, args: T.List, kwargs: TYPE_kwargs) -> T.Dict[str, str]:
        """Returns environment variables for this package."""
        return self.package.get_env_dict(state.environment, state.subdir)

    @noPosargs
    @noKwargs
    def rust_dependency_map_method(self, state: ModuleState, args: T.List, kwargs: TYPE_kwargs) -> T.Dict[str, str]:
        """Returns rust dependency mapping for this package."""
        return self.package.cfg.get_dependency_map(self.package.manifest)


class RustPackage(RustCrate):
    """Represents a Rust package within a workspace."""

    def __init__(self, rust_ws: RustWorkspace, package: cargo.PackageState) -> None:
        super().__init__(rust_ws, package)
        self.methods.update({
            'dependencies': self.dependencies_method,
            'library': self.library_method,
            'proc_macro': self.proc_macro_method,
            'shared_module': self.shared_module_method,
            'executable': self.executable_method,
            'override_dependency': self.override_dependency_method,
        })

    def _dependencies_method(self, state: ModuleState, kwargs: RustPackageDependencies,
                             for_machine: MachineChoice) -> T.List[Dependency]:
        dependencies: T.List[Dependency] = []
        cfg = self.package.cfg

        if kwargs['dependencies']:
            for dep_key, dep_pkg in cfg.dep_packages.items():
                if dep_pkg.manifest.lib:
                    if dep_pkg.ws_subdir != self.rust_ws.subdir or \
                        is_parent_path(os.path.join(self.rust_ws.subdir, state.subproject_dir),
                                       dep_pkg.path):
                        self.rust_ws._do_subproject(dep_pkg)
                    # Get the dependency name for this package (rust or proc-macro ABI)
                    depname = dep_pkg.get_rust_dependency_name()
                    dependency = state.overridden_dependency(depname, for_machine)
                    dependencies.append(dependency)

        if kwargs['dev_dependencies']:
            raise MesonException('dev_dependencies is not implemented yet')

        if kwargs['system_dependencies']:
            for name, sys_dep in self.package.manifest.system_dependencies.items():
                if sys_dep.enabled(cfg.features):
                    # System dependencies use the original dependency name from Cargo.toml
                    dependency = state.dependency(sys_dep.name, native=(for_machine == MachineChoice.BUILD),
                                                  required=not sys_dep.optional,
                                                  wanted=sys_dep.meson_version)
                    dependencies.append(dep_to_system_dependency(dependency, name))

        return dependencies

    @noPosargs
    @typed_kwargs('package.dependencies',
                  KwargInfo('dependencies', bool, default=True),
                  KwargInfo('dev_dependencies', bool, default=False),
                  KwargInfo('system_dependencies', bool, default=True))
    def dependencies_method(self, state: ModuleState, args: T.List, kwargs: RustPackageDependencies) -> T.List[Dependency]:
        """Returns the dependencies for this package."""
        return self._dependencies_method(state, kwargs, MachineChoice.HOST)

    @staticmethod
    def validate_pos_args(name: str, args: T.Tuple[
            T.Optional[T.Union[str, StructuredSources]],
            T.Optional[StructuredSources]]) -> T.Tuple[T.Optional[str], T.Optional[StructuredSources]]:
        if isinstance(args[0], str):
            return args[0], args[1]
        if args[1] is not None:
            raise MesonException(f"{name} only accepts one StructuredSources parameter")
        return None, args[0]

    def merge_kw_args(self, state: ModuleState, kwargs: T.Union[RustPackageExecutable, RustPackageLibrary]) -> None:
        kwargs.setdefault('native', MachineChoice.HOST)

        deps = kwargs['dependencies']
        kwargs['dependencies'] = self._dependencies_method(state, {
            'dependencies': True,
            'dev_dependencies': False,
            'system_dependencies': True,
        }, kwargs['native'])
        kwargs['dependencies'].extend(deps)

        depmap = kwargs['rust_dependency_map']
        kwargs['rust_dependency_map'] = \
            self.package.cfg.get_dependency_map(self.package.manifest)
        kwargs['rust_dependency_map'].update(depmap)

        rust_args = kwargs['rust_args']
        kwargs['rust_args'] = \
            self.package.get_rustc_args(state.environment, state.subdir, kwargs['native'])
        kwargs['rust_args'].extend(rust_args)

        kwargs['override_options'].setdefault('rust_std', self.package.manifest.package.edition)

    def _library_method(self, state: ModuleState, args: T.Tuple[
            T.Optional[T.Union[str, StructuredSources]],
            T.Optional[StructuredSources]], kwargs: RustPackageLibrary,
            static: bool, shared: bool,
            shared_mod: bool = False) -> T.Union[BothLibraries, SharedLibrary, StaticLibrary]:
        tgt_args = self.validate_pos_args('package.library', args)
        if not self.package.manifest.lib:
            raise MesonException("no [lib] section in Cargo package")

        sources: T.Union[StructuredSources, str]
        tgt_name, sources = tgt_args
        if not tgt_name:
            rust_abi: RUST_ABI
            if kwargs['rust_crate_type'] is not None:
                rust_abi = 'rust' if kwargs['rust_crate_type'] in {'lib', 'rlib', 'dylib', 'proc-macro'} else 'c'
            else:
                rust_abi = kwargs['rust_abi']
            tgt_name = self.package.library_name(rust_abi)
        if not sources:
            sources = os.path.relpath(os.path.join(self.package.path, self.package.manifest.lib.path),
                                      state.subdir)

        lib_args: T.Tuple[str, SourcesVarargsType] = (tgt_name, [sources])
        self.merge_kw_args(state, kwargs)

        if shared_mod:
            return state._interpreter.build_target(state.current_node, lib_args,
                                                   T.cast('_kwargs.SharedModule', kwargs),
                                                   SharedModule)

        if static and shared:
            return state._interpreter.build_both_libraries(state.current_node, lib_args, kwargs)
        elif shared:
            return state._interpreter.build_target(state.current_node, lib_args,
                                                   T.cast('_kwargs.SharedLibrary', kwargs),
                                                   SharedLibrary)
        else:
            return state._interpreter.build_target(state.current_node, lib_args,
                                                   T.cast('_kwargs.StaticLibrary', kwargs),
                                                   StaticLibrary)

    def _proc_macro_method(self, state: 'ModuleState', args: T.Tuple[
            T.Optional[T.Union[str, StructuredSources]],
            T.Optional[StructuredSources]], kwargs: RustPackageLibrary) -> SharedLibrary:
        kwargs['native'] = MachineChoice.BUILD
        kwargs['rust_abi'] = None
        kwargs['rust_crate_type'] = 'proc-macro'
        kwargs['rust_args'] = kwargs['rust_args'] + ['--extern', 'proc_macro']
        result = self._library_method(state, args, kwargs, shared=True, static=False)
        return T.cast('SharedLibrary', result)

    @typed_pos_args('package.override_dependency', Dependency)
    @typed_kwargs('package.override_dependency',
                  KwargInfo('rust_abi', (str, NoneType), default=None, validator=in_set_validator({'rust', 'c', 'proc-macro'})))
    def override_dependency_method(self, state: ModuleState, args: T.Tuple[Dependency], kwargs: FuncDependency) -> None:
        dep = args[0]
        rust_abi = self.package.abi_resolve_default(kwargs['rust_abi'])
        depname = self.package.get_dependency_name(rust_abi)

        if rust_abi == 'proc-macro':
            state.override_dependency(depname, dep, for_machine=MachineChoice.HOST)
            state.override_dependency(depname, dep, static=False, for_machine=MachineChoice.HOST)
            if state.environment.is_cross_build():
                state.override_dependency(depname, dep, for_machine=MachineChoice.BUILD)
                state.override_dependency(depname, dep, static=False, for_machine=MachineChoice.BUILD)
            return

        state.override_dependency(depname, dep)
        if self.package.abi_has_static(rust_abi):
            state.override_dependency(depname, dep, static=True)
        if self.package.abi_has_shared(rust_abi):
            state.override_dependency(depname, dep, static=False)

    @typed_pos_args('package.library', optargs=[(str, StructuredSources), StructuredSources])
    @typed_kwargs(
        'package.library',
        *LIBRARY_KWS_NO_NATIVE,
        DEPENDENCIES_KW,
        LINK_WITH_KW,
        LINK_WHOLE_KW,
        _BASE_LANG_KW.evolve(name='rust_args'),
    )
    def library_method(self, state: ModuleState, args: T.Tuple[
            T.Optional[T.Union[str, StructuredSources]],
            T.Optional[StructuredSources]], kwargs: RustPackageLibrary) -> T.Union[BothLibraries, SharedLibrary, StaticLibrary]:
        if not self.package.manifest.lib:
            raise MesonException("no [lib] section in Cargo package")
        if kwargs['rust_crate_type'] is not None:
            static = kwargs['rust_crate_type'] in {'lib', 'rlib', 'staticlib'}
            shared = kwargs['rust_crate_type'] in {'dylib', 'cdylib', 'proc-macro'}
        else:
            rust_abi = self.package.abi_resolve_default(kwargs['rust_abi'])
            static = self.package.abi_has_static(rust_abi)
            shared = self.package.abi_has_shared(rust_abi)
            if rust_abi == 'proc-macro':
                kwargs['rust_crate_type'] = 'proc-macro'
                kwargs['rust_abi'] = None
            else:
                kwargs['rust_abi'] = rust_abi
        return self._library_method(state, args, kwargs, static=static, shared=shared)

    @typed_pos_args('package.proc_macro', optargs=[(str, StructuredSources), StructuredSources])
    @typed_kwargs(
        'package.proc_macro',
        *SHARED_LIB_KWS_NO_NATIVE,
        DEPENDENCIES_KW,
        LINK_WITH_KW,
        LINK_WHOLE_KW,
        _BASE_LANG_KW.evolve(name='rust_args'),
    )
    def proc_macro_method(self, state: 'ModuleState', args: T.Tuple[
            T.Optional[T.Union[str, StructuredSources]],
            T.Optional[StructuredSources]], kwargs: RustPackageLibrary) -> SharedLibrary:
        if not self.package.manifest.lib:
            raise MesonException("no [lib] section in Cargo package")
        if 'proc-macro' not in self.package.manifest.lib.crate_type:
            raise MesonException("not a procedural macro crate")
        return self._proc_macro_method(state, args, kwargs)

    @typed_pos_args('package.shared_module', optargs=[(str, StructuredSources), StructuredSources])
    @typed_kwargs(
        'package.shared_module',
        *SHARED_MOD_KWS_NO_NATIVE,
        DEPENDENCIES_KW,
        LINK_WITH_KW,
        LINK_WHOLE_KW,
        _BASE_LANG_KW.evolve(name='rust_args'),
    )
    def shared_module_method(self, state: 'ModuleState', args: T.Tuple[
            T.Optional[T.Union[str, StructuredSources]],
            T.Optional[StructuredSources]], kwargs: RustPackageLibrary) -> SharedModule:
        if not self.package.manifest.lib:
            raise MesonException("no [lib] section in Cargo package")
        if 'cdylib' not in self.package.manifest.lib.crate_type:
            raise MesonException("not a cdylib crate")

        kwargs['rust_abi'] = None
        kwargs['rust_crate_type'] = 'cdylib'
        result = self._library_method(state, args, kwargs, shared=True, static=False, shared_mod=True)
        return T.cast('SharedModule', result)

    @typed_pos_args('package.executable', optargs=[(str, StructuredSources), StructuredSources])
    @typed_kwargs(
        'package.executable',
        *EXECUTABLE_KWS_NO_NATIVE,
        DEPENDENCIES_KW,
        LINK_WITH_KW,
        LINK_WHOLE_KW,
        _BASE_LANG_KW.evolve(name='rust_args'),
    )
    def executable_method(self, state: 'ModuleState', args: T.Tuple[
            T.Optional[T.Union[str, StructuredSources]],
            T.Optional[StructuredSources]], kwargs: RustPackageExecutable) -> Executable:
        """Builds executable targets from workspace bins."""
        tgt_args = self.validate_pos_args('package.executable', args)
        if not self.package.manifest.bin:
            raise MesonException("no [[bin]] section in Cargo package")

        sources: T.Union[StructuredSources, str]
        tgt_name, sources = tgt_args
        # If there's more than one binary, the first argument must be specified
        # and must be one of the keys in pkg.bin
        if not tgt_name:
            if len(self.package.manifest.bin) > 1:
                raise MesonException("Package has multiple binaries, you must specify which one to build as the first argument")
            # Single binary, use it
            tgt_name = next(iter(self.package.manifest.bin.keys()))
        else:
            if tgt_name not in self.package.manifest.bin:
                raise MesonException(f"Binary '{tgt_name}' not found.")

        if not sources:
            sources = os.path.relpath(os.path.join(self.package.path, self.package.manifest.bin[tgt_name].path),
                                      state.subdir)

        exe_args: T.Tuple[str, SourcesVarargsType] = (tgt_name, [sources])
        self.merge_kw_args(state, kwargs)
        return state._interpreter.build_target(state.current_node, exe_args, kwargs, Executable)


class RustSubproject(RustCrate):
    """Represents a Cargo subproject."""

    def __init__(self, rust_ws: RustWorkspace, package: cargo.PackageState) -> None:
        super().__init__(rust_ws, package)
        self.methods.update({
            'dependency': self.dependency_method,
        })

    @noPosargs
    @typed_kwargs('package.dependency',
                  KwargInfo('rust_abi', (str, NoneType), default=None, validator=in_set_validator({'rust', 'c', 'proc-macro'})))
    def dependency_method(self, state: ModuleState, args: T.List, kwargs: FuncDependency) -> Dependency:
        """Returns dependency for the package with the given ABI."""
        depname = self.package.get_dependency_name(kwargs['rust_abi'])
        return state.overridden_dependency(depname)


class RustModule(ExtensionModule):

    """A module that holds helper functions for rust."""

    INFO = ModuleInfo('rust', '0.57.0', stabilized='1.0.0')
    _bindgen_rust_target: T.Optional[str]
    rustdoc: PerMachine[T.Optional[ExternalProgram]] = PerMachine(None, None)

    def __init__(self, interpreter: Interpreter) -> None:
        super().__init__(interpreter)
        self._bindgen_bin: T.Optional[Program] = None
        if 'rust' in interpreter.compilers.host:
            rustc = T.cast('RustCompiler', interpreter.compilers.host['rust'])
            self._bindgen_rust_target = 'nightly' if rustc.is_nightly else rustc.version
        else:
            self._bindgen_rust_target = None
        self._bindgen_set_std = False
        self.methods.update({
            'test': self.test,
            'doctest': self.doctest,
            'bindgen': self.bindgen,
            'compiler_target': self.compiler_target,
            'proc_macro': self.proc_macro,
            'to_system_dependency': self.to_system_dependency,
            'workspace': self.workspace,
        })

    def test_common(self, funcname: str, state: ModuleState, args: T.Tuple[str, BuildTarget], kwargs: FuncRustTest) -> T.Tuple[Executable, _kwargs.FuncTest]:
        """Generate a rust test target from a given rust target.

        Rust puts its unitests inside its main source files, unlike most
        languages that put them in external files. This means that normally
        you have to define two separate targets with basically the same
        arguments to get tests:

        ```meson
        rust_lib_sources = [...]
        rust_lib = static_library(
            'rust_lib',
            rust_lib_sources,
        )

        rust_lib_test = executable(
            'rust_lib_test',
            rust_lib_sources,
            rust_args : ['--test'],
        )

        test(
            'rust_lib_test',
            rust_lib_test,
            protocol : 'rust',
        )
        ```

        This is all fine, but not very DRY. This method makes it much easier
        to define rust tests:

        ```meson
        rust = import('unstable-rust')

        rust_lib = static_library(
            'rust_lib',
            [sources],
        )

        rust.test('rust_lib_test', rust_lib)
        ```
        """
        if any(isinstance(t, Jar) for t in kwargs.get('link_with', [])):
            raise InvalidArguments('Rust tests cannot link with Jar targets')
        if any(isinstance(t, Jar) for t in kwargs.get('link_whole', [])):
            raise InvalidArguments('Rust tests cannot link with Jar targets')

        name = args[0]
        base_target: BuildTarget = args[1]
        if not base_target.uses_rust():
            raise InterpreterException(f'Second positional argument to rustmod.{funcname}() must be a rust based target')
        extra_args = kwargs['args']

        # Delete any arguments we don't want passed
        if '--test' in extra_args:
            mlog.warning(f'Do not add --test to rustmod.{funcname}() arguments')
            extra_args.remove('--test')
        if '--format' in extra_args:
            mlog.warning(f'Do not add --format to rustmod.{funcname}() arguments')
            i = extra_args.index('--format')
            # Also delete the argument to --format
            del extra_args[i + 1]
            del extra_args[i]
        for i, a in enumerate(extra_args):
            if isinstance(a, str) and a.startswith('--format='):
                del extra_args[i]
                break

        # We need to cast here, as currently these don't have protocol in them, but test itself does.
        tkwargs = T.cast('_kwargs.FuncTest', kwargs.copy())

        tkwargs['args'] = extra_args + ['--test', '--format', 'pretty']
        tkwargs['protocol'] = 'rust'

        new_target_kwargs = T.cast('ExecutableKeywordArguments', base_target.original_kwargs.copy())
        del new_target_kwargs['rust_crate_type']
        for kw in ('pic', 'prelink', 'rust_abi', 'version', 'soversion', 'darwin_versions', 'shortname'):
            if kw in new_target_kwargs:
                del new_target_kwargs[kw]  # type: ignore[misc]

        new_target_kwargs['install'] = False
        new_target_kwargs['dependencies'] = new_target_kwargs.get('dependencies', []) + kwargs['dependencies']
        new_target_kwargs['link_with'] = new_target_kwargs.get('link_with', []) + T.cast('T.List[BuildTargetTypes]', kwargs['link_with'])
        new_target_kwargs['link_whole'] = new_target_kwargs.get('link_whole', []) + kwargs['link_whole']

        lang_args = base_target.extra_args.copy()
        lang_args['rust'] = base_target.extra_args['rust'] + kwargs['rust_args'] + ['--test']
        new_target_kwargs['language_args'] = lang_args

        sources = T.cast('T.List[SourceOutputs]', base_target.sources.copy())
        sources.extend(base_target.generated)

        new_target = Executable(
            name, base_target.subdir, state.subproject, base_target.for_machine,
            sources, base_target.structured_sources,
            base_target.objects, base_target.environment, base_target.compilers,
            new_target_kwargs)
        return new_target, tkwargs

    @typed_pos_args('rust.test', str, BuildTarget)
    @typed_kwargs(
        'rust.test',
        *TEST_KWS,
        DEPENDENCIES_KW,
        LINK_WITH_KW.evolve(since='1.2.0'),
        LINK_WHOLE_KW.evolve(since='1.8.0'),
        *RUST_TEST_KWS,
    )
    def test(self, state: ModuleState, args: T.Tuple[str, BuildTarget], kwargs: FuncTest) -> ModuleReturnValue:
        name, _ = args
        new_target, tkwargs = self.test_common('test', state, args, kwargs)
        test: Test = self.interpreter.make_test(
            self.interpreter.current_node, (name, new_target), tkwargs)

        return ModuleReturnValue(None, [new_target, test])

    @FeatureNew('rust.doctest', '1.8.0')
    @typed_pos_args('rust.doctest', str, BuildTarget)
    @typed_kwargs(
        'rust.doctest',
        *TEST_KWS_NO_ARGS,
        DEPENDENCIES_KW,
        LINK_WITH_KW,
        LINK_WHOLE_KW,
        *RUST_TEST_KWS,
        KwargInfo(
            'args',
            ContainerTypeInfo(list, str),
            listify=True,
            default=[],
            validator=no_spaces_validator,
        ),
    )
    def doctest(self, state: ModuleState, args: T.Tuple[str, T.Union[SharedLibrary, StaticLibrary]], kwargs: FuncDoctest) -> ModuleReturnValue:
        name, base_target = args

        if not base_target.uses_rust():
            raise MesonException('doc tests are only supported for Rust targets')
        if not base_target.uses_rust_abi():
            raise MesonException("doc tests are not supported for rust_abi: 'c'")
        if state.environment.is_cross_build() and state.environment.need_exe_wrapper(base_target.for_machine):
            mlog.notice('skipping Rust doctests due to cross compilation', once=True)
            return ModuleReturnValue(None, [])

        # Link the base target's crate into the tests
        kwargs['link_with'].append(base_target)
        kwargs['depends'].append(base_target)
        workdir = kwargs['workdir']
        kwargs['workdir'] = None
        new_target, tkwargs = self.test_common('doctest', state, args, kwargs)

        # added automatically by rustdoc; keep things simple
        tkwargs['args'].remove('--test')

        # --test-args= is "parsed" simply via the Rust function split_whitespace().
        # This means no quoting nightmares (pfew) but it also means no spaces.
        # Unfortunately it's pretty hard at this point to accept e.g. CustomTarget,
        # because their paths may not be known.  This is not a big deal because the
        # user does not control the test harness, so make things easy and allow
        # strings only.
        if tkwargs['args']:
            tkwargs['args'] = ['--test-args=' + ' '.join(T.cast('T.Sequence[str]', tkwargs['args']))]
        if workdir:
            tkwargs['args'].append('--test-run-directory=' + workdir)

        if self.rustdoc[base_target.for_machine] is None:
            rustc = T.cast('RustCompiler', base_target.compilers['rust'])
            rustdoc = rustc.get_rustdoc()
            if rustdoc:
                self.rustdoc[base_target.for_machine] = ExternalProgram(rustdoc.get_exe())
            else:
                self.rustdoc[base_target.for_machine] = NonExistingExternalProgram()

        rustdoc_prog = self.rustdoc[base_target.for_machine]
        if not rustdoc_prog.found():
            raise MesonException(f'could not find rustdoc for {base_target.for_machine} machine')

        doctests: Doctest = self.interpreter.make_test(
            self.interpreter.current_node, (name, rustdoc_prog), tkwargs, Doctest)

        # Note that the new_target is intentionally not returned, as it
        # is only reached via the base_target and never built by "ninja",
        # so we need to complete its initialization here
        new_target.process_compilers_late()
        doctests.target = new_target
        base_target.doctests = doctests
        return ModuleReturnValue(None, [doctests])

    @noPosargs
    @typed_kwargs(
        'rust.bindgen',
        KwargInfo('c_args', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo('args', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo(
            'input',
            ContainerTypeInfo(list, (File, GeneratedList, BuildTarget, BothLibraries, ExtractedObjects, CustomTargetIndex, CustomTarget, str), allow_empty=False),
            default=[],
            listify=True,
            required=True,
        ),
        KwargInfo('language', (str, NoneType), since='1.4.0', validator=in_set_validator({'c', 'cpp'})),
        KwargInfo('bindgen_version', ContainerTypeInfo(list, str), default=[], listify=True, since='1.4.0'),
        INCLUDE_DIRECTORIES.evolve(since_values={ContainerTypeInfo(list, str): '1.0.0'}),
        OUTPUT_KW,
        KwargInfo(
            'output_inline_wrapper',
            str,
            default='',
            since='1.4.0',
        ),
        DEPENDENCIES_KW.evolve(since='1.0.0'),
    )
    def bindgen(self, state: ModuleState, args: T.List, kwargs: FuncBindgen) -> ModuleReturnValue:
        """Wrapper around bindgen to simplify its use.

        The main thing this simplifies is the use of `include_directory`
        objects, instead of having to pass a plethora of `-I` arguments.
        """
        header, *_deps = self.interpreter.source_strings_to_files(kwargs['input'])

        # Split File and Target dependencies to add pass to CustomTarget
        depends: T.List[SourceOutputs] = []
        depend_files: T.List[File] = []
        for d in _deps:
            if isinstance(d, File):
                depend_files.append(d)
            else:
                depends.append(d)

        # Copy to avoid subsequent calls mutating the original
        # TODO: if we want this to be per-machine we'll need a native kwarg
        clang_args = state.environment.properties.host.get_bindgen_clang_args().copy()

        # Look for --target in the rust command itself if there isn't one passed in clang_args
        target_arg = parse_target(clang_args)
        if not target_arg and 'rust' in state._interpreter.compilers.host:
            rust_args = state._interpreter.compilers.host['rust'].get_exe_args()
            target_arg = parse_target(rust_args)
            if target_arg:
                clang_args.append(f'--target={target_arg}')

        # Find the first C'ish compiler to fetch the default compiler flags
        # from. Append those to the bindgen flags to ensure we use a compatible
        # environment.
        comp = mesonlib.first(
            # The cast here shouldn't be necessary
            [state.environment.coredata.compilers.host.get(T.cast('Language', l))
             for l in ('c', 'cpp', 'objc', 'objcpp')],
            lambda x: x is not None,
        )
        if comp:
            clang_args.extend(comp.get_always_args())
        else:
            mlog.warning(textwrap.dedent('''\
                Using `rust.bindgen` without configuring C (or a C-like)
                language in Meson will skip compiler detection and can cause
                ABI incompatibilities due to missing crucial compiler flags.
                Consider calling `add_languages('c')` in your Meson build
                files.
            '''))

        for i in state.process_include_dirs(kwargs['include_directories']):
            # bindgen always uses clang, so it's safe to hardcode -I here
            clang_args.extend([f'-I{x}' for x in i.abs_string_list(
                state.environment.get_source_dir(), state.environment.get_build_dir())])
        if are_asserts_disabled_for_subproject(state.subproject, state.environment):
            clang_args.append('-DNDEBUG')

        for de in kwargs['dependencies']:
            for i in de.get_include_dirs():
                clang_args.extend([f'-I{x}' for x in i.abs_string_list(
                    state.environment.get_source_dir(), state.environment.get_build_dir())])
            clang_args.extend(de.get_all_compile_args())
            for s in de.get_sources():
                if isinstance(s, File):
                    depend_files.append(s)
                elif isinstance(s, CustomTarget):
                    depends.append(s)

        if self._bindgen_bin is None:
            self._bindgen_bin = state.find_program('bindgen', wanted=kwargs['bindgen_version'])
            if self._bindgen_rust_target is not None:
                _, _, err = mesonlib.Popen_safe(self._bindgen_bin.get_command() + ['--rust-target', self._bindgen_rust_target])
                # < 0.71: Sometimes this is "invalid Rust target" and
                # sometimes "invalid # rust target"
                # >= 0.71: error: invalid value '...' for '--rust-target <RUST_TARGET>': "..." is not a valid Rust target, accepted values are of the form ...
                # It's also much harder to hit this in 0.71 than in previous versions
                if 'Got an invalid' in err or 'is not a valid Rust target' in err:
                    self._bindgen_rust_target = None

            self._bindgen_set_std = mesonlib.version_compare(self._bindgen_bin.get_version(), '>= 0.71')

        name: str
        if isinstance(header, File):
            name = header.fname
        elif isinstance(header, (BuildTarget, BothLibraries, ExtractedObjects, StructuredSources)):
            raise InterpreterException('bindgen source file must be a C header, not an object or build target')
        else:
            name = header.get_outputs()[0]

        # bindgen assumes that C++ headers will be called .hpp. We want to
        # ensure that anything Meson considers a C++ header is treated as one.
        language = kwargs['language']
        if language is None:
            ext = os.path.splitext(name)[1][1:]
            if ext in lang_suffixes['cpp']:
                language = 'cpp'
            elif ext == 'h':
                language = 'c'
            else:
                raise InterpreterException(f'Unknown file type extension for: {name}')

        # We only want include directories and defines, other things may not be valid
        cargs = state.get_option(f'{language}_args', state.subproject)
        assert isinstance(cargs, list), 'for mypy'
        for a in itertools.chain(state.global_args.get(language, []), state.project_args.get(language, []), cargs):
            if a.startswith(('-I', '/I', '-D', '/D', '-U', '/U')):
                clang_args.append(a)

        if language == 'cpp':
            clang_args.extend(['-x', 'c++'])

        # Add the C++ standard to the clang arguments. Attempt to translate VS
        # extension versions into the nearest standard version
        std = state.get_option(f'{language}_std')
        assert isinstance(std, str), 'for mypy'
        if std.startswith('vc++'):
            if std.endswith('latest'):
                mlog.warning('Attempting to translate vc++latest into a clang compatible version.',
                             'Currently this is hardcoded for c++20', once=True, fatal=False)
                std = 'c++20'
            else:
                mlog.debug('The current C++ standard is a Visual Studio extension version.',
                           'bindgen will use a the nearest C++ standard instead')
                std = std[1:]

        if std != 'none':
            clang_args.append(f'-std={std}')

        inline_wrapper_args: T.List[str] = []
        outputs = [kwargs['output']]
        if kwargs['output_inline_wrapper']:
            # Todo drop this isinstance once Executable supports version_compare
            if isinstance(self._bindgen_bin, ExternalProgram):
                if mesonlib.version_compare(self._bindgen_bin.get_version(), '< 0.65'):
                    raise InterpreterException('\'output_inline_wrapper\' parameter of rust.bindgen requires bindgen-0.65 or newer')

            outputs.append(kwargs['output_inline_wrapper'])
            inline_wrapper_args = [
                '--experimental', '--wrap-static-fns',
                '--wrap-static-fns-path', os.path.join(state.environment.build_dir, '@OUTPUT1@')
            ]

        cmd = self._bindgen_bin.get_command() + \
            [
                '@INPUT@', '--output',
                os.path.join(state.environment.build_dir, '@OUTPUT0@')
            ] + \
            kwargs['args'] + inline_wrapper_args
        if self._bindgen_rust_target and '--rust-target' not in cmd:
            cmd.extend(['--rust-target', self._bindgen_rust_target])
        if self._bindgen_set_std and '--rust-edition' not in cmd:
            try:
                rust_std = state.environment.coredata.optstore.get_value_for('rust_std')
            except KeyError:
                rust_std = 'none'
            assert isinstance(rust_std, str), 'for mypy'
            if rust_std != 'none':
                cmd.extend(['--rust-edition', rust_std])
        cmd.append('--')
        cmd.extend(kwargs['c_args'])
        cmd.extend(clang_args)
        cmd.extend(['-MD', '-MQ', '@INPUT@', '-MF', '@DEPFILE@'])

        target = CustomTarget(
            f'rustmod-bindgen-{name}'.replace('/', '_'),
            state.subdir,
            state.subproject,
            state.environment,
            cmd,
            [header],
            outputs,
            depfile='@PLAINNAME@.d',
            extra_depends=depends,
            depend_files=depend_files,
            backend=state.backend,
            description='Generating bindings for Rust {}',
        )

        return ModuleReturnValue(target, [target])

    @FeatureNew('rust.compiler_target', '1.11.0')
    @noPosargs
    @typed_kwargs('rust.compiler_target', NATIVE_KW)
    def compiler_target(self, state: ModuleState, args: T.List, kwargs: '_kwargs.NativeKW') -> str:
        """Returns the Rust target triple for the specified machine's Rust compiler."""
        for_machine = kwargs['native']
        compilers = state._interpreter.coredata.compilers[for_machine]
        if 'rust' in compilers:
            rustc = T.cast('RustCompiler', compilers['rust'])
            return rustc.get_target_triple()
        else:
            raise MesonException(f'No Rust compiler was requested for the {for_machine} machine')

    # Allow a limited set of kwargs, but still use the full set of typed_kwargs()
    # because it could be setting required default values.
    @FeatureNew('rust.proc_macro', '1.3.0')
    @permittedKwargs({'rust_args', 'rust_dependency_map', 'sources', 'dependencies', 'extra_files',
                      'link_args', 'link_depends', 'link_with', 'override_options'})
    @typed_pos_args('rust.proc_macro', str, varargs=SOURCES_VARARGS)
    @typed_kwargs('rust.proc_macro', *SHARED_LIB_KWS)
    def proc_macro(self, state: ModuleState, args: T.Tuple[str, SourcesVarargsType], kwargs: _kwargs.SharedLibrary) -> SharedLibrary:
        kwargs['native'] = MachineChoice.BUILD
        kwargs['rust_crate_type'] = 'proc-macro'
        kwargs['rust_args'] = kwargs['rust_args'] + ['--extern', 'proc_macro']
        target = state._interpreter.build_target(state.current_node, args, kwargs, SharedLibrary)
        return target

    @FeatureNew('rust.to_system_dependency', '1.11.0')
    @typed_pos_args('rust.to_system_dependency', Dependency, optargs=[str])
    @noKwargs
    def to_system_dependency(self, state: ModuleState, args: T.Tuple[Dependency, T.Optional[str]], kwargs: TYPE_kwargs) -> Dependency:
        dep, depname = args
        return dep_to_system_dependency(dep, depname)

    @FeatureNew('rust.workspace', '1.11.0')
    @noPosargs
    @typed_kwargs(
        'rust.workspace',
        KwargInfo('default_features', (bool, NoneType), default=None),
        KwargInfo(
            'features',
            (ContainerTypeInfo(list, str), NoneType),
            default=None,
            listify=True,
        ),
    )
    def workspace(self, state: ModuleState, args: T.List, kwargs: FuncWorkspace) -> RustWorkspace:
        """Creates a Rust workspace object, controlling the build of
           all the packages in a Cargo.lock file."""
        if self.interpreter.cargo is None:
            raise MesonException("rust.workspace() requires a Cargo project (Cargo.toml and Cargo.lock)")

        self.interpreter.add_languages(['rust'], True, MachineChoice.HOST)
        self.interpreter.add_languages(['rust'], True, MachineChoice.BUILD)

        default_features = kwargs['default_features']
        features = kwargs['features']
        if default_features is not None or features is not None:
            # If custom features are provided, default_features = None should be treated as True
            if default_features is None:
                default_features = True

            cargo_features = ['default'] if default_features else []
            if features is not None:
                cargo_features.extend(features)
            self.interpreter.cargo.features = cargo_features

        ws = self.interpreter.cargo.load_workspace(state.root_subdir)

        # Cargo projects may not have a subprojects directory, because
        # dependencies are declared in Cargo.toml rather than .wrap files.
        # Create it now so that the wrap resolver can download the crates there
        os.makedirs(self.interpreter.environment.wrap_resolver.subdir_root, exist_ok=True)

        return RustWorkspace(self.interpreter, ws)


def initialize(interp: Interpreter) -> RustModule:
    return RustModule(interp)
