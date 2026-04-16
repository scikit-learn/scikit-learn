# SPDX-License-Identifier: Apache-2.0
# Copyright © 2022-2024 Intel Corporation

"""Interpreter for converting Cargo Toml definitions to Meson AST

There are some notable limits here. We don't even try to convert something with
a build.rs: there's so few limits on what Cargo allows a build.rs (basically
none), and no good way for us to convert them. In that case, an actual meson
port will be required.
"""

from __future__ import annotations
import dataclasses
import functools
import itertools
import os
import pathlib
import collections
import urllib.parse
import typing as T
from pathlib import PurePath

from . import builder, version
from .cfg import eval_cfg
from .toml import load_toml
from .manifest import Manifest, CargoLock, CargoLockPackage, Workspace, fixup_meson_varname
from ..interpreterbase import SubProject
from ..mesonlib import (
    is_parent_path, lazy_property, MesonException, MachineChoice,
    unique_list, version_compare)
from .. import coredata, mlog
from ..wrap.wrap import PackageDefinition

if T.TYPE_CHECKING:
    from . import raw
    from .. import mparser
    from typing_extensions import Literal

    from .manifest import Dependency
    from ..environment import Environment
    from ..compilers.rust import RustCompiler

    RUST_ABI = Literal['rust', 'c', 'proc-macro']

def _dependency_name(package_name: str, api: str, suffix: str = '-rs') -> str:
    basename = package_name[:-len(suffix)] if suffix and package_name.endswith(suffix) else package_name
    return f'{basename}-{api}{suffix}'


def _extra_args_varname() -> str:
    return 'extra_args'


def _extra_deps_varname() -> str:
    return 'extra_deps'


@dataclasses.dataclass
class PackageConfiguration:
    """Configuration for a package during dependency resolution."""
    features: T.Set[str] = dataclasses.field(default_factory=set)
    required_deps: T.Set[str] = dataclasses.field(default_factory=set)
    optional_deps_features: T.Dict[str, T.Set[str]] = dataclasses.field(default_factory=lambda: collections.defaultdict(set))
    # Cache of resolved dependency packages
    dep_packages: T.Dict[PackageKey, PackageState] = dataclasses.field(default_factory=dict)

    def get_features_args(self) -> T.List[str]:
        """Get feature configuration arguments."""
        args: T.List[str] = []
        for feature in sorted(self.features):
            args.extend(['--cfg', f'feature="{feature}"'])
        return args

    def get_dependency_map(self, manifest: Manifest) -> T.Dict[str, str]:
        """Get the rust dependency mapping for this package configuration."""
        dependency_map: T.Dict[str, str] = {}
        for name in sorted(self.required_deps):
            dep = manifest.dependencies[name]
            dep_key = PackageKey(dep.package, dep.api)
            dep_pkg = self.dep_packages[dep_key]
            dep_lib_name = dep_pkg.library_name()
            dep_crate_name = name if name != dep.package else dep_pkg.manifest.lib.name
            dependency_map[dep_lib_name] = dep_crate_name
        return dependency_map


@dataclasses.dataclass
class PackageState:
    manifest: Manifest
    downloaded: bool = False
    # If this package is member of a workspace.
    ws_subdir: T.Optional[str] = None
    ws_member: T.Optional[str] = None
    # Package configuration state
    cfg: T.Optional[PackageConfiguration] = None
    # Subproject name as known to the wrap resolver (may differ from the
    # meson dep name for git sources, where the wrap is named after the
    # git directory rather than the crate name + api version).
    subproject_name: T.Optional[str] = None

    @lazy_property
    def path(self) -> T.Optional[str]:
        if not self.ws_subdir:
            return None
        return os.path.normpath(os.path.join(self.ws_subdir, self.ws_member))

    def library_name(self, lib_type: RUST_ABI = 'rust') -> str:
        # Add the API version to the library name to avoid conflicts when multiple
        # versions of the same crate are used. The Ninja backend removed everything
        # after the + to form the crate name.
        name = fixup_meson_varname(self.manifest.package.name)
        if lib_type == 'c':
            return name
        return f'{name}+{self.manifest.package.api.replace(".", "_")}'

    def get_env_dict(self, environment: Environment, subdir: str) -> T.Dict[str, str]:
        """Get environment variables for this package."""
        # Common variables for build.rs and crates
        # https://doc.rust-lang.org/cargo/reference/environment-variables.html
        # OUT_DIR is the directory where build.rs generate files. In our case,
        # it's the directory where meson/meson.build places generated files.
        out_dir = os.path.join(environment.build_dir, subdir, 'meson')
        os.makedirs(out_dir, exist_ok=True)
        version_arr = self.manifest.package.version.split('.')
        version_arr += [''] * (4 - len(version_arr))

        return {
            'OUT_DIR': out_dir,
            'CARGO_MANIFEST_DIR': os.path.join(environment.source_dir, subdir),
            'CARGO_MANIFEST_PATH': os.path.join(environment.source_dir, subdir, 'Cargo.toml'),
            'CARGO_PKG_VERSION': self.manifest.package.version,
            'CARGO_PKG_VERSION_MAJOR': version_arr[0],
            'CARGO_PKG_VERSION_MINOR': version_arr[1],
            'CARGO_PKG_VERSION_PATCH': version_arr[2],
            'CARGO_PKG_VERSION_PRE': version_arr[3],
            'CARGO_PKG_AUTHORS': ','.join(self.manifest.package.authors),
            'CARGO_PKG_NAME': self.manifest.package.name,
            # FIXME: description can contain newlines which breaks ninja.
            #'CARGO_PKG_DESCRIPTION': self.manifest.package.description or '',
            'CARGO_PKG_HOMEPAGE': self.manifest.package.homepage or '',
            'CARGO_PKG_REPOSITORY': self.manifest.package.repository or '',
            'CARGO_PKG_LICENSE': self.manifest.package.license or '',
            'CARGO_PKG_LICENSE_FILE': self.manifest.package.license_file or '',
            'CARGO_PKG_RUST_VERSION': self.manifest.package.rust_version or '',
            'CARGO_PKG_README': self.manifest.package.readme or '',
            'CARGO_CRATE_NAME': fixup_meson_varname(self.manifest.package.name),
        }

    def get_lint_args(self, rustc: RustCompiler) -> T.List[str]:
        """Get lint arguments for this package."""
        args: T.List[str] = []
        has_check_cfg = rustc.has_check_cfg

        for lint in self.manifest.lints:
            args.extend(lint.to_arguments(has_check_cfg))

        if has_check_cfg:
            args.append('--check-cfg')
            args.append('cfg(test)')
            for feature in self.manifest.features:
                if feature != 'default':
                    args.append('--check-cfg')
                    args.append(f'cfg(feature,values("{feature}"))')
            for name in self.manifest.system_dependencies:
                args.append('--check-cfg')
                args.append(f'cfg(system_deps_have_{fixup_meson_varname(name)})')

        return args

    def get_env_args(self, rustc: RustCompiler, environment: Environment, subdir: str) -> T.List[str]:
        """Get environment variable arguments for rustc."""
        enable_env_set_args = rustc.enable_env_set_args()
        if enable_env_set_args is None:
            return []

        env_dict = self.get_env_dict(environment, subdir)
        env_args = list(enable_env_set_args)
        for k, v in env_dict.items():
            env_args.extend(['--env-set', f'{k}={v}'])
        return env_args

    def get_rustc_args(self, environment: Environment, subdir: str, machine: MachineChoice) -> T.List[str]:
        """Get rustc arguments for this package."""
        if not environment.is_cross_build():
            machine = MachineChoice.HOST

        rustc = T.cast('RustCompiler', environment.coredata.compilers[machine]['rust'])

        cfg = self.cfg

        args: T.List[str] = []
        args.extend(self.get_lint_args(rustc))
        args.extend(cfg.get_features_args())
        args.extend(self.get_env_args(rustc, environment, subdir))
        return args

    def supported_abis(self) -> T.Set[RUST_ABI]:
        """Return which ABIs are exposed by the package's crate_types."""
        crate_types = self.manifest.lib.crate_type
        abis: T.Set[RUST_ABI] = set()
        if any(ct in {'lib', 'rlib', 'dylib'} for ct in crate_types):
            abis.add('rust')
        if any(ct in {'staticlib', 'cdylib'} for ct in crate_types):
            abis.add('c')
        if 'proc-macro' in crate_types:
            abis.add('proc-macro')
        return abis

    def get_subproject_name(self) -> SubProject:
        if self.subproject_name is not None:
            return SubProject(self.subproject_name)
        dep = _dependency_name(self.manifest.package.name, self.manifest.package.api)
        return SubProject(dep)

    def abi_resolve_default(self, rust_abi: T.Optional[RUST_ABI]) -> RUST_ABI:
        supported_abis = self.supported_abis()
        if rust_abi is None:
            if len(supported_abis) > 1:
                raise MesonException(f'Package {self.manifest.package.name} support more than one ABI')
            return next(iter(supported_abis))
        else:
            if rust_abi not in supported_abis:
                raise MesonException(f'Package {self.manifest.package.name} does not support ABI {rust_abi}')
            return rust_abi

    def abi_has_shared(self, rust_abi: RUST_ABI) -> bool:
        if rust_abi == 'proc-macro':
            return True
        return ('cdylib' if rust_abi == 'c' else 'dylib') in self.manifest.lib.crate_type

    def abi_has_static(self, rust_abi: RUST_ABI) -> bool:
        if rust_abi == 'proc-macro':
            return False
        crate_type = self.manifest.lib.crate_type
        if rust_abi == 'c':
            return 'staticlib' in crate_type
        return 'lib' in crate_type or 'rlib' in crate_type

    def get_dependency_name(self, rust_abi: T.Optional[RUST_ABI]) -> str:
        """Get the dependency name for a package with the given ABI."""
        rust_abi = self.abi_resolve_default(rust_abi)
        package_name = self.manifest.package.name
        api = self.manifest.package.api

        if rust_abi in {'rust', 'proc-macro'}:
            return _dependency_name(package_name, api)
        elif rust_abi == 'c':
            return _dependency_name(package_name, api, '')
        else:
            raise MesonException(f'Unknown rust_abi: {rust_abi}')

    def get_rust_dependency_name(self) -> str:
        """Get the dependency name for a package with the rust or proc-macro ABI."""
        supported_abis = self.supported_abis()
        package_name = self.manifest.package.name
        if 'rust' in supported_abis or 'proc-macro' in supported_abis:
            return _dependency_name(package_name, self.manifest.package.api)
        raise MesonException(f'Package {package_name} does not support rust or proc-macro ABI')

@dataclasses.dataclass(frozen=True)
class PackageKey:
    package_name: str
    api: str


@dataclasses.dataclass
class WorkspaceState:
    workspace: Workspace
    subdir: str
    downloaded: bool = False
    # member path -> PackageState, for all members of this workspace
    packages: T.Dict[str, PackageState] = dataclasses.field(default_factory=dict)
    # package name to member path, for all members of this workspace
    packages_to_member: T.Dict[str, str] = dataclasses.field(default_factory=dict)
    # member paths that are required to be built
    required_members: T.List[str] = dataclasses.field(default_factory=list)


class Interpreter:
    _features: T.Optional[T.List[str]] = None

    def __init__(self, env: Environment, subdir: str, subprojects_dir: str) -> None:
        self.environment = env
        self.subprojects_dir = subprojects_dir
        # Map Cargo.toml's subdir to loaded manifest.
        self.manifests: T.Dict[str, T.Union[Manifest, Workspace]] = {}
        # Map of cargo package (name + api) to its state
        self.packages: T.Dict[PackageKey, PackageState] = {}
        # Map subdir to workspace
        self.workspaces: T.Dict[str, WorkspaceState] = {}
        # Files that should trigger a reconfigure if modified
        self.build_def_files: T.List[str] = []
        # Cargo packages
        filename = os.path.join(self.environment.get_source_dir(), subdir, 'Cargo.lock')
        subprojects_dir = os.path.join(self.environment.get_source_dir(), subdir, subprojects_dir)
        self.cargolock = load_cargo_lock(filename, subprojects_dir)
        if self.cargolock:
            self.environment.wrap_resolver.merge_wraps(self.cargolock.wraps)
            self.build_def_files.append(filename)

    @property
    def features(self) -> T.List[str]:
        """Get the features list. Once read, it cannot be modified."""
        if self._features is None:
            self._features = ['default']
        return self._features

    @features.setter
    def features(self, value: T.List[str]) -> None:
        """Set the features list. Can only be set before first read."""
        value_unique = sorted(unique_list(value))
        if self._features is not None and value_unique != self._features:
            raise MesonException("Cannot modify features after they have been selected or used")
        self._features = value_unique

    def get_build_def_files(self) -> T.List[str]:
        return self.build_def_files

    def load_workspace(self, subdir: str) -> WorkspaceState:
        """Load the root Cargo.toml package and prepare it with features and dependencies."""
        subdir = os.path.normpath(subdir)
        manifest, cached = self._load_manifest(subdir)
        ws = self._get_workspace(manifest, subdir, False)
        if not cached:
            self._prepare_entry_point(ws)
        return ws

    def _prepare_entry_point(self, ws: WorkspaceState) -> None:
        pkgs = [self._require_workspace_member(ws, m) for m in ws.workspace.default_members]
        for pkg in pkgs:
            self._prepare_package(pkg)
            for feature in self.features:
                self._enable_feature(pkg, feature)

    def load_package(self, ws: WorkspaceState, package_name: T.Optional[str]) -> PackageState:
        if package_name is None:
            if not ws.workspace.root_package:
                raise MesonException('no root package in workspace')
            path = '.'
        else:
            try:
                path = ws.packages_to_member[package_name]
            except KeyError:
                raise MesonException(f'workspace member "{package_name}" not found')

        if is_parent_path(self.subprojects_dir, path):
            raise MesonException('argument to package() cannot be a subproject')
        return ws.packages[path]

    def interpret(self, subdir: str, project_root: T.Optional[str] = None) -> mparser.CodeBlockNode:
        filename = os.path.join(self.environment.source_dir, subdir, 'Cargo.toml')
        build = builder.Builder(filename)
        if project_root:
            # this is a subdir()
            manifest, _ = self._load_manifest(subdir)
            assert isinstance(manifest, Manifest)
            return self.interpret_package(manifest, build, subdir, project_root)
        else:
            ws = self.load_workspace(subdir)
            return self.interpret_workspace(ws, build, subdir)

    def interpret_package(self, manifest: Manifest, build: builder.Builder, subdir: str, project_root: str) -> mparser.CodeBlockNode:
        # Build an AST for this package
        ws = self.workspaces[project_root]
        member = ws.packages_to_member[manifest.package.name]
        pkg = ws.packages[member]
        ast = self._create_package(pkg, build, subdir)
        return build.block(ast)

    def _create_package(self, pkg: PackageState, build: builder.Builder, subdir: str) -> T.List[mparser.BaseNode]:
        ast: T.List[mparser.BaseNode] = [
            build.assign(build.method('package', build.identifier('cargo_ws'),
                                      [build.string(pkg.manifest.package.name)]), 'pkg_obj'),
            build.assign(build.method('features', build.identifier('pkg_obj')), 'features'),
            build.function('message', [
                build.string('Enabled features:'),
                build.identifier('features'),
            ]),
        ]
        ast += self._create_feature_checks(pkg, build)
        ast += self._create_meson_subdir(build)

        if pkg.manifest.lib:
            crate_type = pkg.manifest.lib.crate_type
            if 'dylib' in crate_type and 'cdylib' in crate_type:
                raise MesonException('Cannot build both dylib and cdylib due to file name conflict')
            for abi in pkg.supported_abis():
                ast.extend(self._create_lib(pkg, build, subdir, abi))

        return ast

    def interpret_workspace(self, ws: WorkspaceState, build: builder.Builder, subdir: str) -> mparser.CodeBlockNode:
        name = os.path.basename(subdir)
        subprojects_dir = os.path.join(subdir, 'subprojects')
        self.environment.wrap_resolver.load_and_merge(subprojects_dir, SubProject(name))
        ast: T.List[mparser.BaseNode] = []

        # Call subdir() for each required member of the workspace. The order is
        # important, if a member depends on another member, that member must be
        # processed first.
        processed_members: T.Dict[str, PackageState] = {}

        def _process_member(member: str) -> None:
            if member in processed_members:
                return
            pkg = ws.packages[member]
            cfg = pkg.cfg
            if not cfg:
                raise MesonException(f'Package {pkg.manifest.package.name!r} is not enabled for this build '
                                     'configuration. Maybe you forgot to enable a Cargo feature, or to check '
                                     'a Meson option?')
            for depname in cfg.required_deps:
                dep = pkg.manifest.dependencies[depname]
                if dep.path:
                    dep_member = os.path.normpath(os.path.join(pkg.ws_member, dep.path))
                    _process_member(dep_member)
            if member == '.':
                ast.extend(self._create_package(pkg, build, subdir))
            elif is_parent_path(self.subprojects_dir, member):
                depname = _dependency_name(pkg.manifest.package.name, pkg.manifest.package.api)
                ast.append(build.function('subproject', [build.string(depname)]))
            else:
                ast.append(build.function('subdir', [build.string(member)]))
            processed_members[member] = pkg

        for member in ws.required_members:
            _process_member(member)
        ast = self._create_project(name, processed_members.get('.'), build) + ast
        return build.block(ast)

    def _load_workspace_member(self, ws: WorkspaceState, m: str) -> None:
        m = os.path.normpath(m)
        if m in ws.packages:
            return
        # Load member's manifest
        m_subdir = os.path.join(ws.subdir, m)
        manifest_, _ = self._load_manifest(m_subdir, ws.workspace, m)
        assert isinstance(manifest_, Manifest)
        self._add_workspace_member(manifest_, ws, m)

    def _add_workspace_member(self, manifest_: Manifest, ws: WorkspaceState, m: str) -> None:
        key = PackageKey(manifest_.package.name, manifest_.package.api)
        ws.packages_to_member[manifest_.package.name] = m
        if key in self.packages:
            ws.packages[m] = self.packages[key]
            self._require_workspace_member(ws, m)
        else:
            ws.packages[m] = PackageState(manifest_, ws_subdir=ws.subdir, ws_member=m, downloaded=ws.downloaded)

    def _get_workspace(self, manifest: T.Union[Workspace, Manifest], subdir: str, downloaded: bool) -> WorkspaceState:
        ws = self.workspaces.get(subdir)
        if ws:
            return ws
        workspace = manifest if isinstance(manifest, Workspace) else \
            Workspace(root_package=manifest, members=['.'], default_members=['.'])
        ws = WorkspaceState(workspace, subdir, downloaded=downloaded)
        if workspace.root_package:
            self._add_workspace_member(workspace.root_package, ws, '.')
        for m in workspace.members:
            self._load_workspace_member(ws, m)
        self.workspaces[subdir] = ws
        return ws

    def _record_package(self, pkg: PackageState) -> None:
        key = PackageKey(pkg.manifest.package.name, pkg.manifest.package.api)
        if key not in self.packages:
            self.packages[key] = pkg

    def _require_workspace_member(self, ws: WorkspaceState, member: str) -> PackageState:
        member = os.path.normpath(member)
        pkg = ws.packages[member]
        if member not in ws.required_members:
            self._record_package(pkg)
            ws.required_members.append(member)
        return pkg

    def _fetch_package(self, package_name: str, api: str) -> PackageState:
        key = PackageKey(package_name, api)
        pkg = self.packages.get(key)
        if pkg:
            return pkg
        return self._fetch_package_from_provider(package_name, api)

    def _resolve_package(self, package_name: str, version_constraints: T.List[str]) -> T.Optional[CargoLockPackage]:
        """From all available versions from Cargo.lock, pick the most recent
           satisfying the constraints and return it."""
        if self.cargolock:
            cargo_lock_pkgs = self.cargolock.named(package_name)
        else:
            cargo_lock_pkgs = []
        for cargo_pkg in cargo_lock_pkgs:
            if all(version_compare(cargo_pkg.version, v) for v in version_constraints):
                return cargo_pkg

        if not version_constraints:
            raise MesonException(f'Cannot determine version of cargo package {package_name}')
        return None

    def resolve_package(self, package_name: str, api: str) -> T.Optional[PackageState]:
        cargo_pkg = self._resolve_package(package_name, version.convert(api))
        if not cargo_pkg:
            return None
        api = version.api(cargo_pkg.version)
        return self._fetch_package(package_name, api)

    def _fetch_package_from_provider(self, package_name: str, api: str) -> PackageState:
        meson_depname = _dependency_name(package_name, api)
        subp_name, _ = self.environment.wrap_resolver.find_dep_provider(meson_depname)
        if subp_name is None:
            if self.cargolock is None:
                raise MesonException(f'Dependency {meson_depname!r} not found in any wrap files.')
            # If Cargo.lock has a different version, this could be a resolution
            # bug, but maybe also a version mismatch?  I am not sure yet...
            similar_deps = [pkg.subproject
                            for pkg in self.cargolock.named(package_name)]
            if similar_deps:
                similar_msg = f'Cargo.lock provides: {", ".join(similar_deps)}.'
            else:
                similar_msg = 'Cargo.lock does not contain this crate name.'
            raise MesonException(f'Dependency {meson_depname!r} not found in any wrap files or Cargo.lock; {similar_msg} This could be a Meson bug, please report it.')

        return self._fetch_package_from_subproject(package_name, subp_name)

    def _fetch_package_from_subproject(self, package_name: str, subp_name: str) -> PackageState:
        subdir, _ = self.environment.wrap_resolver.resolve(subp_name)
        subprojects_dir = os.path.join(subdir, 'subprojects')
        self.environment.wrap_resolver.load_and_merge(subprojects_dir, SubProject(subp_name))
        manifest, _ = self._load_manifest(subdir)
        downloaded = \
            subp_name in self.environment.wrap_resolver.wraps and \
            self.environment.wrap_resolver.wraps[subp_name].type is not None

        ws = self._get_workspace(manifest, subdir, downloaded=downloaded)
        member = ws.packages_to_member[package_name]
        pkg = self._require_workspace_member(ws, member)
        pkg.subproject_name = subp_name
        return pkg

    def _prepare_package(self, pkg: PackageState) -> None:
        key = PackageKey(pkg.manifest.package.name, pkg.manifest.package.api)
        assert key in self.packages
        if pkg.cfg:
            return

        pkg.cfg = PackageConfiguration()
        # Merge target specific dependencies that are enabled
        cfgs = self._get_cfgs(MachineChoice.HOST)
        for condition, dependencies in pkg.manifest.target.items():
            if eval_cfg(condition, cfgs):
                pkg.manifest.dependencies.update(dependencies)

        # If you specify the optional dependency with the dep: prefix anywhere in the [features]
        # table, that disables the implicit feature.
        deps = set(feature[4:]
                   for feature in itertools.chain.from_iterable(pkg.manifest.features.values())
                   if feature.startswith('dep:'))
        for name, dep in itertools.chain(pkg.manifest.dependencies.items(),
                                         pkg.manifest.dev_dependencies.items(),
                                         pkg.manifest.build_dependencies.items()):
            if dep.optional and name not in deps:
                pkg.manifest.features.setdefault(name, [])
                pkg.manifest.features[name].append(f'dep:{name}')
                deps.add(name)

        # Fetch required dependencies recursively.
        for depname, dep in pkg.manifest.dependencies.items():
            if not dep.optional:
                self._add_dependency(pkg, depname)

    def _dep_package(self, pkg: PackageState, dep: Dependency) -> PackageState:
        if dep.path:
            ws = self.workspaces[pkg.ws_subdir]
            dep_member = os.path.normpath(os.path.join(pkg.ws_member, dep.path))
            if is_parent_path(self.subprojects_dir, dep_member):
                if len(pathlib.PurePath(dep_member).parts) != 2:
                    raise MesonException('found "{self.subprojects_dir}" in path but it is not a valid subproject path')
            self._load_workspace_member(ws, dep_member)
            dep_pkg = self._require_workspace_member(ws, dep_member)
        elif dep.git:
            _, _, directory = _parse_git_url(dep.git, dep.branch)
            dep_pkg = self._fetch_package_from_subproject(dep.package, directory)
        else:
            cargo_pkg = self._resolve_package(dep.package, dep.meson_version)
            if cargo_pkg:
                dep.update_version(f'={cargo_pkg.version}')
            dep_pkg = self._fetch_package(dep.package, dep.api)

        if not dep.version:
            dep.update_version(f'={dep_pkg.manifest.package.version}')

        dep_key = PackageKey(dep.package, dep.api)
        pkg.cfg.dep_packages.setdefault(dep_key, dep_pkg)
        assert pkg.cfg.dep_packages[dep_key] == dep_pkg
        return dep_pkg

    def _load_manifest(self, subdir: str, workspace: T.Optional[Workspace] = None, member_path: str = '') -> T.Tuple[T.Union[Manifest, Workspace], bool]:
        manifest_ = self.manifests.get(subdir)
        if manifest_:
            return manifest_, True
        path = os.path.join(self.environment.source_dir, subdir)
        filename = os.path.join(path, 'Cargo.toml')
        try:
            raw_manifest = T.cast('raw.Manifest', load_toml(filename))
        except OSError as e:
            raise MesonException(f'could not load {subdir}/Cargo.toml: {e}')

        self.build_def_files.append(filename)
        if 'workspace' in raw_manifest:
            manifest_ = Workspace.from_raw(raw_manifest, path)
        elif 'package' in raw_manifest:
            manifest_ = Manifest.from_raw(raw_manifest, path, workspace, member_path)
        else:
            raise MesonException(f'{subdir}/Cargo.toml does not have [package] or [workspace] section')
        self.manifests[subdir] = manifest_
        return manifest_, False

    def _add_dependency(self, pkg: PackageState, depname: str) -> None:
        cfg = pkg.cfg
        if depname in cfg.required_deps:
            return
        dep = pkg.manifest.dependencies.get(depname)
        if not dep:
            # It could be build/dev/target dependency. Just ignore it.
            return
        cfg.required_deps.add(depname)
        dep_pkg = self._dep_package(pkg, dep)
        self._prepare_package(dep_pkg)
        if dep.default_features:
            self._enable_feature(dep_pkg, 'default')
        for f in dep.features:
            self._enable_feature(dep_pkg, f)
        for f in cfg.optional_deps_features[depname]:
            self._enable_feature(dep_pkg, f)

    def _enable_feature(self, pkg: PackageState, feature: str) -> None:
        cfg = pkg.cfg
        if feature in cfg.features:
            return
        cfg.features.add(feature)
        # Recurse on extra features and dependencies this feature pulls.
        # https://doc.rust-lang.org/cargo/reference/features.html#the-features-section
        for f in pkg.manifest.features.get(feature, []):
            if '/' in f:
                depname, dep_f = f.split('/', 1)
                if depname[-1] == '?':
                    depname = depname[:-1]
                else:
                    self._add_dependency(pkg, depname)
                if depname in cfg.required_deps:
                    dep = pkg.manifest.dependencies[depname]
                    dep_pkg = self._dep_package(pkg, dep)
                    self._enable_feature(dep_pkg, dep_f)
                else:
                    # This feature will be enabled only if that dependency
                    # is later added.
                    cfg.optional_deps_features[depname].add(dep_f)
            elif f.startswith('dep:'):
                self._add_dependency(pkg, f[4:])
            else:
                self._enable_feature(pkg, f)

    def has_check_cfg(self, machine: MachineChoice) -> bool:
        if not self.environment.is_cross_build():
            machine = MachineChoice.HOST
        rustc = T.cast('RustCompiler', self.environment.coredata.compilers[machine]['rust'])
        return rustc.has_check_cfg

    @functools.lru_cache(maxsize=None)
    def _get_cfgs(self, machine: MachineChoice) -> T.Dict[str, str]:
        if not self.environment.is_cross_build():
            machine = MachineChoice.HOST
        rustc = T.cast('RustCompiler', self.environment.coredata.compilers[machine]['rust'])
        cfgs = rustc.get_cfgs().copy()
        rustflags = self.environment.coredata.get_external_args(machine, 'rust')
        rustflags_i = iter(rustflags)
        for i in rustflags_i:
            if i == '--cfg':
                cfgs.append(next(rustflags_i))
        return dict(self._split_cfg(i) for i in cfgs)

    @staticmethod
    def _split_cfg(cfg: str) -> T.Tuple[str, str]:
        pair = cfg.split('=', maxsplit=1)
        value = pair[1] if len(pair) > 1 else ''
        if value and value[0] == '"':
            value = value[1:-1]
        return pair[0], value

    def _create_project(self, name: str, pkg: T.Optional[PackageState], build: builder.Builder) -> T.List[mparser.BaseNode]:
        """Create the project() function call

        :param pkg: The package to generate from
        :param build: The AST builder
        :return: a list nodes
        """
        args: T.List[mparser.BaseNode] = [
            build.string(name),
            build.string('rust'),
        ]
        kwargs: T.Dict[str, mparser.BaseNode] = {
            # Always assume that the generated meson is using the latest features
            # This will warn when when we generate deprecated code, which is helpful
            # for the upkeep of the module
            'meson_version': build.string(f'>= {coredata.stable_version}'),
        }
        if pkg:
            default_options: T.Dict[str, mparser.BaseNode] = {}
            if pkg.downloaded:
                default_options['warning_level'] = build.string('0')

            kwargs.update({
                'version': build.string(pkg.manifest.package.version),
                'default_options': build.dict({build.string(k): v for k, v in default_options.items()}),
            })
            if pkg.manifest.package.license:
                kwargs['license'] = build.string(pkg.manifest.package.license)
            elif pkg.manifest.package.license_file:
                kwargs['license_files'] = build.string(pkg.manifest.package.license_file)

        # project(...)
        # rust = import('rust')
        # cargo_ws = rust.workspace()
        return [
            build.function('project', args, kwargs),
            build.assign(build.function('import', [build.string('rust')]),
                         'rust'),
            build.assign(build.method('workspace', build.identifier('rust'), []),
                         'cargo_ws')
        ]

    def _create_feature_checks(self, pkg: PackageState, build: builder.Builder) -> T.List[mparser.BaseNode]:
        cfg = pkg.cfg
        ast: T.List[mparser.BaseNode] = []
        for depname in cfg.required_deps:
            dep = pkg.manifest.dependencies[depname]
            dep_pkg = self._dep_package(pkg, dep)
            if dep_pkg.manifest.lib:
                ast += self._create_feature_check(dep_pkg, dep, build)
        return ast

    def _create_feature_check(self, pkg: PackageState, dep: Dependency, build: builder.Builder) -> T.List[mparser.BaseNode]:
        cfg = pkg.cfg
        feat_obj: mparser.BaseNode
        feat_pkg = self.cargolock and self.resolve_package(dep.package, dep.api)
        if feat_pkg:
            if feat_pkg.ws_subdir == pkg.ws_subdir:
                return []

            feat_obj = build.method(
                'features',
                build.method(
                    'subproject',
                    build.identifier('cargo_ws'),
                    [build.string(dep.package), build.string(dep.api)]))
        else:
            version_ = dep.meson_version or [pkg.manifest.package.version]
            kw = {
                'version': build.array([build.string(s) for s in version_]),
            }
            # actual_features = dependency(...).get_variable('features', default_value : '').split(',')
            dep_obj = build.function(
                 'dependency',
                 [build.string(_dependency_name(dep.package, dep.api))],
                 kw)
            feat_obj = build.method(
                'split',
                build.method(
                    'get_variable',
                    dep_obj,
                    [build.string('features')],
                    {'default_value': build.string('')}
                ),
                [build.string(',')],
            )

        # However, this subproject could have been previously configured with a
        # different set of features. Cargo collects the set of features globally
        # but Meson can only use features enabled by the first call that triggered
        # the configuration of that subproject.
        #
        # Verify all features that we need are actually enabled for that dependency,
        # otherwise abort with an error message. The user has to set the corresponding
        # option manually with -Dxxx-rs:feature-yyy=true, or the main project can do
        # that in its project(..., default_options: ['xxx-rs:feature-yyy=true']).
        return [
            # actual_features = dependency(...).get_variable('features', default_value : '').split(',')
            build.assign(
                feat_obj,
                'actual_features'
            ),
            # needed_features = [f1, f2, ...]
            # foreach f : needed_features
            #   if f not in actual_features
            #     error()
            #   endif
            # endforeach
            build.assign(build.array([build.string(f) for f in cfg.features]), 'needed_features'),
            build.foreach(['f'], build.identifier('needed_features'), build.block([
                build.if_(build.not_in(build.identifier('f'), build.identifier('actual_features')), build.block([
                    build.function('error', [
                        build.string('Dependency'),
                        build.string(_dependency_name(dep.package, dep.api)),
                        build.string('previously configured with features'),
                        build.identifier('actual_features'),
                        build.string('but need'),
                        build.identifier('needed_features'),
                    ])
                ]))
            ])),
        ]

    def _create_meson_subdir(self, build: builder.Builder) -> T.List[mparser.BaseNode]:
        # Allow Cargo subprojects to add extra Rust args in meson/meson.build file.
        # This is used to replace build.rs logic.

        # extra_args = []
        # extra_deps = []
        # fs = import('fs')
        # if fs.is_dir('meson')
        #  subdir('meson')
        # endif
        return [
            build.assign(build.array([]), _extra_args_varname()),
            build.assign(build.array([]), _extra_deps_varname()),
            build.assign(build.function('import', [build.string('fs')]), 'fs'),
            build.if_(build.method('is_dir', build.identifier('fs'), [build.string('meson')]),
                      build.block([build.function('subdir', [build.string('meson')])]))
        ]

    def _create_lib(self, pkg: PackageState, build: builder.Builder, subdir: str,
                    lib_type: RUST_ABI) -> T.List[mparser.BaseNode]:
        posargs: T.List[mparser.BaseNode] = [
            build.string(pkg.library_name(lib_type)),
        ]

        kwargs: T.Dict[str, mparser.BaseNode] = {
            'dependencies': build.identifier(_extra_deps_varname()),
            'rust_args': build.identifier(_extra_args_varname()),
        }

        if lib_type == 'proc-macro':
            lib = build.method('proc_macro', build.identifier('pkg_obj'), posargs, kwargs)
        else:
            kwargs['rust_abi'] = build.string(lib_type)
            lib = build.method('library', build.identifier('pkg_obj'), posargs, kwargs)

        # lib = xxx_library()
        # dep = declare_dependency()
        # meson.override_dependency()
        return [
            build.assign(lib, 'lib'),
            build.assign(
                build.function(
                    'declare_dependency',
                    kw={
                        'link_with': build.identifier('lib'),
                        'variables': build.dict({
                            build.string('features'): build.method('join', build.string(','),
                                                                   [build.identifier('features')]),
                        }),
                        'version': build.method('version', build.identifier('pkg_obj')),
                    },
                ),
                'dep'
            ),
            build.method(
                'override_dependency',
                build.identifier('pkg_obj'),
                [build.identifier('dep')],
                {'rust_abi': build.string(lib_type)}
            ),
        ]


def _parse_git_url(url: str, branch: T.Optional[str] = None) -> T.Tuple[str, str, str]:
    if url.startswith('git+'):
        url = url[4:]
    parts = urllib.parse.urlparse(url)
    query = urllib.parse.parse_qs(parts.query)
    query_branch = query['branch'][0] if 'branch' in query else ''
    branch = branch or query_branch
    revision = parts.fragment or branch
    directory = PurePath(parts.path).name
    if directory.endswith('.git'):
        directory = directory[:-4]
    if branch:
        directory += f'-{branch}'
    url = urllib.parse.urlunparse(parts._replace(params='', query='', fragment=''))
    return url, revision, directory


def load_cargo_lock(filename: str, subproject_dir: str) -> T.Optional[CargoLock]:
    """ Convert Cargo.lock into a list of wraps """

    # Map directory -> PackageDefinition, to avoid duplicates. Multiple packages
    # can have the same source URL, in that case we have a single wrap that
    # provides multiple dependency names.
    if os.path.exists(filename):
        toml = load_toml(filename)
        raw_cargolock = T.cast('raw.CargoLock', toml)
        cargolock = CargoLock.from_raw(raw_cargolock)
        packagefiles_dir = os.path.join(subproject_dir, 'packagefiles')
        wraps: T.Dict[str, PackageDefinition] = {}
        for package in cargolock.package:
            meson_depname = _dependency_name(package.name, version.api(package.version))
            if package.source is None:
                # This is project's package, or one of its workspace members.
                continue
            elif package.source == 'registry+https://github.com/rust-lang/crates.io-index':
                checksum = package.checksum
                if checksum is None:
                    checksum = cargolock.metadata[f'checksum {package.name} {package.version} ({package.source})']
                url = f'https://crates.io/api/v1/crates/{package.name}/{package.version}/download'
                directory = f'{package.name}-{package.version}'
                name = meson_depname
                wrap_type = 'file'
                cfg = {
                    'directory': directory,
                    'source_url': url,
                    'source_filename': f'{directory}.tar.gz',
                    'source_hash': checksum,
                    'method': 'cargo',
                }
            elif package.source.startswith('git+'):
                url, revision, directory = _parse_git_url(package.source)
                name = directory
                wrap_type = 'git'
                cfg = {
                    'url': url,
                    'revision': revision,
                    'method': 'cargo',
                }
            else:
                mlog.warning(f'Unsupported source URL in {filename}: {package.source}')
                continue
            if os.path.isdir(os.path.join(packagefiles_dir, name)):
                cfg['patch_directory'] = name
            if directory not in wraps:
                wraps[directory] = PackageDefinition.from_values(name, subproject_dir, wrap_type, cfg)
            wraps[directory].add_provided_dep(meson_depname)
        cargolock.wraps = {w.name: w for w in wraps.values()}
        return cargolock
    return None
