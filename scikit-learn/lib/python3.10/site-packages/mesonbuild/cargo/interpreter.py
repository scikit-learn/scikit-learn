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
import os
import collections
import urllib.parse
import itertools
import typing as T

from . import builder, version, cfg
from .toml import load_toml, TomlImplementationMissing
from .manifest import Manifest, CargoLock, fixup_meson_varname
from ..mesonlib import MesonException, MachineChoice
from .. import coredata, mlog
from ..wrap.wrap import PackageDefinition

if T.TYPE_CHECKING:
    from . import raw
    from .. import mparser
    from .manifest import Dependency, SystemDependency
    from ..environment import Environment
    from ..interpreterbase import SubProject
    from ..compilers.rust import RustCompiler

def _dependency_name(package_name: str, api: str, suffix: str = '-rs') -> str:
    basename = package_name[:-len(suffix)] if package_name.endswith(suffix) else package_name
    return f'{basename}-{api}{suffix}'


def _dependency_varname(package_name: str) -> str:
    return f'{fixup_meson_varname(package_name)}_dep'


def _extra_args_varname() -> str:
    return 'extra_args'


def _extra_deps_varname() -> str:
    return 'extra_deps'


@dataclasses.dataclass
class PackageState:
    manifest: Manifest
    downloaded: bool = False
    features: T.Set[str] = dataclasses.field(default_factory=set)
    required_deps: T.Set[str] = dataclasses.field(default_factory=set)
    optional_deps_features: T.Dict[str, T.Set[str]] = dataclasses.field(default_factory=lambda: collections.defaultdict(set))


@dataclasses.dataclass(frozen=True)
class PackageKey:
    package_name: str
    api: str


class Interpreter:
    def __init__(self, env: Environment) -> None:
        self.environment = env
        self.host_rustc = T.cast('RustCompiler', self.environment.coredata.compilers[MachineChoice.HOST]['rust'])
        # Map Cargo.toml's subdir to loaded manifest.
        self.manifests: T.Dict[str, Manifest] = {}
        # Map of cargo package (name + api) to its state
        self.packages: T.Dict[PackageKey, PackageState] = {}
        # Rustc's config
        self.cfgs = self._get_cfgs()

    def get_build_def_files(self) -> T.List[str]:
        return [os.path.join(subdir, 'Cargo.toml') for subdir in self.manifests]

    def interpret(self, subdir: str) -> mparser.CodeBlockNode:
        manifest = self._load_manifest(subdir)
        pkg, cached = self._fetch_package(manifest.package.name, manifest.package.api)
        if not cached:
            # This is an entry point, always enable the 'default' feature.
            # FIXME: We should have a Meson option similar to `cargo build --no-default-features`
            self._enable_feature(pkg, 'default')

        # Build an AST for this package
        filename = os.path.join(self.environment.source_dir, subdir, 'Cargo.toml')
        build = builder.Builder(filename)
        ast = self._create_project(pkg, build)
        ast += [
            build.assign(build.function('import', [build.string('rust')]), 'rust'),
            build.function('message', [
                build.string('Enabled features:'),
                build.array([build.string(f) for f in pkg.features]),
            ]),
        ]
        ast += self._create_dependencies(pkg, build)
        ast += self._create_meson_subdir(build)

        if pkg.manifest.lib:
            for crate_type in pkg.manifest.lib.crate_type:
                ast.extend(self._create_lib(pkg, build, crate_type))

        return build.block(ast)

    def _fetch_package(self, package_name: str, api: str) -> T.Tuple[PackageState, bool]:
        key = PackageKey(package_name, api)
        pkg = self.packages.get(key)
        if pkg:
            return pkg, True
        meson_depname = _dependency_name(package_name, api)
        subdir, _ = self.environment.wrap_resolver.resolve(meson_depname)
        subprojects_dir = os.path.join(subdir, 'subprojects')
        self.environment.wrap_resolver.load_and_merge(subprojects_dir, T.cast('SubProject', meson_depname))
        manifest = self._load_manifest(subdir)
        downloaded = \
            meson_depname in self.environment.wrap_resolver.wraps and \
            self.environment.wrap_resolver.wraps[meson_depname].type is not None
        pkg = PackageState(manifest, downloaded)
        self.packages[key] = pkg
        # Merge target specific dependencies that are enabled
        for condition, dependencies in manifest.target.items():
            if cfg.eval_cfg(condition, self.cfgs):
                manifest.dependencies.update(dependencies)
        # Fetch required dependencies recursively.
        for depname, dep in manifest.dependencies.items():
            if not dep.optional:
                self._add_dependency(pkg, depname)
        return pkg, False

    def _dep_package(self, dep: Dependency) -> PackageState:
        return self.packages[PackageKey(dep.package, dep.api)]

    def _load_manifest(self, subdir: str) -> Manifest:
        manifest_ = self.manifests.get(subdir)
        if not manifest_:
            path = os.path.join(self.environment.source_dir, subdir)
            filename = os.path.join(path, 'Cargo.toml')
            toml = load_toml(filename)
            if 'package' in toml:
                raw_manifest = T.cast('raw.Manifest', toml)
                manifest_ = Manifest.from_raw(raw_manifest, path)
                self.manifests[subdir] = manifest_
            else:
                raise MesonException(f'{subdir}/Cargo.toml does not have [package] section')
        return manifest_

    def _add_dependency(self, pkg: PackageState, depname: str) -> None:
        if depname in pkg.required_deps:
            return
        dep = pkg.manifest.dependencies.get(depname)
        if not dep:
            if depname in itertools.chain(pkg.manifest.dev_dependencies, pkg.manifest.build_dependencies):
                # FIXME: Not supported yet
                return
            raise MesonException(f'Dependency {depname} not defined in {pkg.manifest.package.name} manifest')
        pkg.required_deps.add(depname)
        dep_pkg, _ = self._fetch_package(dep.package, dep.api)
        if dep.default_features:
            self._enable_feature(dep_pkg, 'default')
        for f in dep.features:
            self._enable_feature(dep_pkg, f)
        for f in pkg.optional_deps_features[depname]:
            self._enable_feature(dep_pkg, f)

    def _enable_feature(self, pkg: PackageState, feature: str) -> None:
        if feature in pkg.features:
            return
        pkg.features.add(feature)
        # A feature can also be a dependency.
        if feature in pkg.manifest.dependencies:
            self._add_dependency(pkg, feature)
        # Recurse on extra features and dependencies this feature pulls.
        # https://doc.rust-lang.org/cargo/reference/features.html#the-features-section
        for f in pkg.manifest.features.get(feature, []):
            if '/' in f:
                depname, dep_f = f.split('/', 1)
                if depname[-1] == '?':
                    depname = depname[:-1]
                    if depname in pkg.required_deps:
                        dep = pkg.manifest.dependencies[depname]
                        dep_pkg = self._dep_package(dep)
                        self._enable_feature(dep_pkg, dep_f)
                    else:
                        # This feature will be enabled only if that dependency
                        # is later added.
                        pkg.optional_deps_features[depname].add(dep_f)
                else:
                    self._add_dependency(pkg, depname)
                    dep = pkg.manifest.dependencies.get(depname)
                    if dep:
                        dep_pkg = self._dep_package(dep)
                        self._enable_feature(dep_pkg, dep_f)
            elif f.startswith('dep:'):
                self._add_dependency(pkg, f[4:])
            else:
                self._enable_feature(pkg, f)

    def _get_cfgs(self) -> T.Dict[str, str]:
        cfgs = self.host_rustc.get_cfgs().copy()
        rustflags = self.environment.coredata.get_external_args(MachineChoice.HOST, 'rust')
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

    def _create_project(self, pkg: PackageState, build: builder.Builder) -> T.List[mparser.BaseNode]:
        """Create the project() function call

        :param pkg: The package to generate from
        :param build: The AST builder
        :return: a list nodes
        """
        default_options: T.List[mparser.BaseNode] = []
        default_options.append(build.string(f'rust_std={pkg.manifest.package.edition}'))
        default_options.append(build.string(f'build.rust_std={pkg.manifest.package.edition}'))
        if pkg.downloaded:
            default_options.append(build.string('warning_level=0'))

        args: T.List[mparser.BaseNode] = []
        args.extend([
            build.string(pkg.manifest.package.name),
            build.string('rust'),
        ])
        kwargs: T.Dict[str, mparser.BaseNode] = {
            'version': build.string(pkg.manifest.package.version),
            # Always assume that the generated meson is using the latest features
            # This will warn when when we generate deprecated code, which is helpful
            # for the upkeep of the module
            'meson_version': build.string(f'>= {coredata.stable_version}'),
            'default_options': build.array(default_options),
        }
        if pkg.manifest.package.license:
            kwargs['license'] = build.string(pkg.manifest.package.license)
        elif pkg.manifest.package.license_file:
            kwargs['license_files'] = build.string(pkg.manifest.package.license_file)

        return [build.function('project', args, kwargs)]

    def _create_dependencies(self, pkg: PackageState, build: builder.Builder) -> T.List[mparser.BaseNode]:
        ast: T.List[mparser.BaseNode] = []
        for depname in pkg.required_deps:
            dep = pkg.manifest.dependencies[depname]
            ast += self._create_dependency(dep, build)
        ast.append(build.assign(build.array([]), 'system_deps_args'))
        for name, sys_dep in pkg.manifest.system_dependencies.items():
            if sys_dep.enabled(pkg.features):
                ast += self._create_system_dependency(name, sys_dep, build)
        return ast

    def _create_system_dependency(self, name: str, dep: SystemDependency, build: builder.Builder) -> T.List[mparser.BaseNode]:
        # TODO: handle feature_overrides
        kw = {
            'version': build.array([build.string(s) for s in dep.meson_version]),
            'required': build.bool(not dep.optional),
        }
        varname = f'{fixup_meson_varname(name)}_system_dep'
        cfg = f'system_deps_have_{fixup_meson_varname(name)}'
        return [
            build.assign(
                build.function(
                    'dependency',
                    [build.string(dep.name)],
                    kw,
                ),
                varname,
            ),
            build.if_(
                build.method('found', build.identifier(varname)), build.block([
                    build.plusassign(
                        build.array([build.string('--cfg'), build.string(cfg)]),
                        'system_deps_args'
                    ),
                ])
            ),
        ]

    def _create_dependency(self, dep: Dependency, build: builder.Builder) -> T.List[mparser.BaseNode]:
        pkg = self._dep_package(dep)
        kw = {
            'version': build.array([build.string(s) for s in dep.meson_version]),
        }
        # Lookup for this dependency with the features we want in default_options kwarg.
        #
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
            # xxx_dep = dependency('xxx', version : ...)
            build.assign(
                build.function(
                    'dependency',
                    [build.string(_dependency_name(dep.package, dep.api))],
                    kw,
                ),
                _dependency_varname(dep.package),
            ),
            # actual_features = xxx_dep.get_variable('features', default_value : '').split(',')
            build.assign(
                build.method(
                    'split',
                    build.method(
                        'get_variable',
                        build.identifier(_dependency_varname(dep.package)),
                        [build.string('features')],
                        {'default_value': build.string('')}
                    ),
                    [build.string(',')],
                ),
                'actual_features'
            ),
            # needed_features = [f1, f2, ...]
            # foreach f : needed_features
            #   if f not in actual_features
            #     error()
            #   endif
            # endforeach
            build.assign(build.array([build.string(f) for f in pkg.features]), 'needed_features'),
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

    def _create_lib(self, pkg: PackageState, build: builder.Builder, crate_type: raw.CRATE_TYPE) -> T.List[mparser.BaseNode]:
        dependencies: T.List[mparser.BaseNode] = []
        dependency_map: T.Dict[mparser.BaseNode, mparser.BaseNode] = {}
        for name in pkg.required_deps:
            dep = pkg.manifest.dependencies[name]
            dependencies.append(build.identifier(_dependency_varname(dep.package)))
            if name != dep.package:
                dep_pkg = self._dep_package(dep)
                dep_lib_name = dep_pkg.manifest.lib.name
                dependency_map[build.string(fixup_meson_varname(dep_lib_name))] = build.string(name)
        for name, sys_dep in pkg.manifest.system_dependencies.items():
            if sys_dep.enabled(pkg.features):
                dependencies.append(build.identifier(f'{fixup_meson_varname(name)}_system_dep'))

        rust_args: T.List[mparser.BaseNode] = [
            build.identifier('features_args'),
            build.identifier(_extra_args_varname()),
            build.identifier('system_deps_args'),
        ]

        dependencies.append(build.identifier(_extra_deps_varname()))

        posargs: T.List[mparser.BaseNode] = [
            build.string(fixup_meson_varname(pkg.manifest.lib.name)),
            build.string(pkg.manifest.lib.path),
        ]

        kwargs: T.Dict[str, mparser.BaseNode] = {
            'dependencies': build.array(dependencies),
            'rust_dependency_map': build.dict(dependency_map),
            'rust_args': build.array(rust_args),
        }

        depname_suffix = '-rs' if crate_type in {'lib', 'rlib', 'proc-macro'} else f'-{crate_type}'
        depname = _dependency_name(pkg.manifest.package.name, pkg.manifest.package.api, depname_suffix)

        lib: mparser.BaseNode
        if pkg.manifest.lib.proc_macro or crate_type == 'proc-macro':
            lib = build.method('proc_macro', build.identifier('rust'), posargs, kwargs)
        else:
            if crate_type in {'lib', 'rlib', 'staticlib'}:
                target_type = 'static_library'
            elif crate_type in {'dylib', 'cdylib'}:
                target_type = 'shared_library'
            else:
                raise MesonException(f'Unsupported crate type {crate_type}')
            if crate_type in {'staticlib', 'cdylib'}:
                kwargs['rust_abi'] = build.string('c')
            lib = build.function(target_type, posargs, kwargs)

        features_args: T.List[mparser.BaseNode] = []
        for f in pkg.features:
            features_args += [build.string('--cfg'), build.string(f'feature="{f}"')]

        # features_args = ['--cfg', 'feature="f1"', ...]
        # lib = xxx_library()
        # dep = declare_dependency()
        # meson.override_dependency()
        return [
            build.assign(build.array(features_args), 'features_args'),
            build.assign(lib, 'lib'),
            build.assign(
                build.function(
                    'declare_dependency',
                    kw={
                        'link_with': build.identifier('lib'),
                        'variables': build.dict({
                            build.string('features'): build.string(','.join(pkg.features)),
                        }),
                        'version': build.string(pkg.manifest.package.version),
                    },
                ),
                'dep'
            ),
            build.method(
                'override_dependency',
                build.identifier('meson'),
                [
                    build.string(depname),
                    build.identifier('dep'),
                ],
            ),
        ]


def load_wraps(source_dir: str, subproject_dir: str) -> T.List[PackageDefinition]:
    """ Convert Cargo.lock into a list of wraps """

    wraps: T.List[PackageDefinition] = []
    filename = os.path.join(source_dir, 'Cargo.lock')
    if os.path.exists(filename):
        try:
            toml = load_toml(filename)
        except TomlImplementationMissing as e:
            mlog.warning('Failed to load Cargo.lock:', str(e), fatal=False)
            return wraps
        raw_cargolock = T.cast('raw.CargoLock', toml)
        cargolock = CargoLock.from_raw(raw_cargolock)
        for package in cargolock.package:
            subp_name = _dependency_name(package.name, version.api(package.version))
            if package.source is None:
                # This is project's package, or one of its workspace members.
                pass
            elif package.source == 'registry+https://github.com/rust-lang/crates.io-index':
                checksum = package.checksum
                if checksum is None:
                    checksum = cargolock.metadata[f'checksum {package.name} {package.version} ({package.source})']
                url = f'https://crates.io/api/v1/crates/{package.name}/{package.version}/download'
                directory = f'{package.name}-{package.version}'
                wraps.append(PackageDefinition.from_values(subp_name, subproject_dir, 'file', {
                    'directory': directory,
                    'source_url': url,
                    'source_filename': f'{directory}.tar.gz',
                    'source_hash': checksum,
                    'method': 'cargo',
                }))
            elif package.source.startswith('git+'):
                parts = urllib.parse.urlparse(package.source[4:])
                query = urllib.parse.parse_qs(parts.query)
                branch = query['branch'][0] if 'branch' in query else ''
                revision = parts.fragment or branch
                url = urllib.parse.urlunparse(parts._replace(params='', query='', fragment=''))
                wraps.append(PackageDefinition.from_values(subp_name, subproject_dir, 'git', {
                    'directory': package.name,
                    'url': url,
                    'revision': revision,
                    'method': 'cargo',
                }))
            else:
                mlog.warning(f'Unsupported source URL in {filename}: {package.source}')
    return wraps
