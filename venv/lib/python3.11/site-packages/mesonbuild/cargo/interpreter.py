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
import importlib
import json
import os
import shutil
import collections
import urllib.parse
import itertools
import typing as T

from . import builder
from . import version
from ..mesonlib import MesonException, Popen_safe
from .. import coredata, mlog
from ..wrap.wrap import PackageDefinition

if T.TYPE_CHECKING:
    from types import ModuleType

    from typing_extensions import Protocol, Self

    from . import manifest
    from .. import mparser
    from ..environment import Environment
    from ..interpreterbase import SubProject

    # Copied from typeshed. Blarg that they don't expose this
    class DataclassInstance(Protocol):
        __dataclass_fields__: T.ClassVar[dict[str, dataclasses.Field[T.Any]]]

    _UnknownKeysT = T.TypeVar('_UnknownKeysT', manifest.FixedPackage,
                              manifest.FixedDependency, manifest.FixedLibTarget,
                              manifest.FixedBuildTarget)


# tomllib is present in python 3.11, before that it is a pypi module called tomli,
# we try to import tomllib, then tomli,
# TODO: add a fallback to toml2json?
tomllib: T.Optional[ModuleType] = None
toml2json: T.Optional[str] = None
for t in ['tomllib', 'tomli']:
    try:
        tomllib = importlib.import_module(t)
        break
    except ImportError:
        pass
else:
    # TODO: it would be better to use an Executable here, which could be looked
    # up in the cross file or provided by a wrap. However, that will have to be
    # passed in externally, since we don't have (and I don't think we should),
    # have access to the `Environment` for that in this module.
    toml2json = shutil.which('toml2json')


_EXTRA_KEYS_WARNING = (
    "This may (unlikely) be an error in the cargo manifest, or may be a missing "
    "implementation in Meson. If this issue can be reproduced with the latest "
    "version of Meson, please help us by opening an issue at "
    "https://github.com/mesonbuild/meson/issues. Please include the crate and "
    "version that is generating this warning if possible."
)

class TomlImplementationMissing(MesonException):
    pass


def load_toml(filename: str) -> T.Dict[object, object]:
    if tomllib:
        with open(filename, 'rb') as f:
            raw = tomllib.load(f)
    else:
        if toml2json is None:
            raise TomlImplementationMissing('Could not find an implementation of tomllib, nor toml2json')

        p, out, err = Popen_safe([toml2json, filename])
        if p.returncode != 0:
            raise MesonException('toml2json failed to decode output\n', err)

        raw = json.loads(out)

    if not isinstance(raw, dict):
        raise MesonException("Cargo.toml isn't a dictionary? How did that happen?")

    return raw


def fixup_meson_varname(name: str) -> str:
    """Fixup a meson variable name

    :param name: The name to fix
    :return: the fixed name
    """
    return name.replace('-', '_')


# Pylance can figure out that these do not, in fact, overlap, but mypy can't
@T.overload
def _fixup_raw_mappings(d: manifest.BuildTarget) -> manifest.FixedBuildTarget: ...  # type: ignore

@T.overload
def _fixup_raw_mappings(d: manifest.LibTarget) -> manifest.FixedLibTarget: ...  # type: ignore

@T.overload
def _fixup_raw_mappings(d: manifest.Dependency) -> manifest.FixedDependency: ...

def _fixup_raw_mappings(d: T.Union[manifest.BuildTarget, manifest.LibTarget, manifest.Dependency]
                        ) -> T.Union[manifest.FixedBuildTarget, manifest.FixedLibTarget,
                                     manifest.FixedDependency]:
    """Fixup raw cargo mappings to ones more suitable for python to consume.

    This does the following:
    * replaces any `-` with `_`, cargo likes the former, but python dicts make
      keys with `-` in them awkward to work with
    * Convert Dependency versions from the cargo format to something meson
      understands

    :param d: The mapping to fix
    :return: the fixed string
    """
    raw = {fixup_meson_varname(k): v for k, v in d.items()}
    if 'version' in raw:
        assert isinstance(raw['version'], str), 'for mypy'
        raw['version'] = version.convert(raw['version'])
    return T.cast('T.Union[manifest.FixedBuildTarget, manifest.FixedLibTarget, manifest.FixedDependency]', raw)


def _handle_unknown_keys(data: _UnknownKeysT, cls: T.Union[DataclassInstance, T.Type[DataclassInstance]],
                         msg: str) -> _UnknownKeysT:
    """Remove and warn on keys that are coming from cargo, but are unknown to
    our representations.

    This is intended to give users the possibility of things proceeding when a
    new key is added to Cargo.toml that we don't yet handle, but to still warn
    them that things might not work.

    :param data: The raw data to look at
    :param cls: The Dataclass derived type that will be created
    :param msg: the header for the error message. Usually something like "In N structure".
    :return: The original data structure, but with all unknown keys removed.
    """
    unexpected = set(data) - {x.name for x in dataclasses.fields(cls)}
    if unexpected:
        mlog.warning(msg, 'has unexpected keys', '"{}".'.format(', '.join(sorted(unexpected))),
                     _EXTRA_KEYS_WARNING)
        for k in unexpected:
            # Mypy and Pyright can't prove that this is okay
            del data[k]  # type: ignore[misc]
    return data


@dataclasses.dataclass
class Package:

    """Representation of a Cargo Package entry, with defaults filled in."""

    name: str
    version: str
    description: T.Optional[str] = None
    resolver: T.Optional[str] = None
    authors: T.List[str] = dataclasses.field(default_factory=list)
    edition: manifest.EDITION = '2015'
    rust_version: T.Optional[str] = None
    documentation: T.Optional[str] = None
    readme: T.Optional[str] = None
    homepage: T.Optional[str] = None
    repository: T.Optional[str] = None
    license: T.Optional[str] = None
    license_file: T.Optional[str] = None
    keywords: T.List[str] = dataclasses.field(default_factory=list)
    categories: T.List[str] = dataclasses.field(default_factory=list)
    workspace: T.Optional[str] = None
    build: T.Optional[str] = None
    links: T.Optional[str] = None
    exclude: T.List[str] = dataclasses.field(default_factory=list)
    include: T.List[str] = dataclasses.field(default_factory=list)
    publish: bool = True
    metadata: T.Dict[str, T.Any] = dataclasses.field(default_factory=dict)
    default_run: T.Optional[str] = None
    autolib: bool = True
    autobins: bool = True
    autoexamples: bool = True
    autotests: bool = True
    autobenches: bool = True
    api: str = dataclasses.field(init=False)

    def __post_init__(self) -> None:
        self.api = _version_to_api(self.version)

    @classmethod
    def from_raw(cls, raw: manifest.Package) -> Self:
        pkg = T.cast('manifest.FixedPackage',
                     {fixup_meson_varname(k): v for k, v in raw.items()})
        pkg = _handle_unknown_keys(pkg, cls, f'Package entry {pkg["name"]}')
        return cls(**pkg)

@dataclasses.dataclass
class SystemDependency:

    """ Representation of a Cargo system-deps entry
        https://docs.rs/system-deps/latest/system_deps
    """

    name: str
    version: T.List[str]
    optional: bool = False
    feature: T.Optional[str] = None
    feature_overrides: T.Dict[str, T.Dict[str, str]] = dataclasses.field(default_factory=dict)

    @classmethod
    def from_raw(cls, name: str, raw: T.Any) -> SystemDependency:
        if isinstance(raw, str):
            return cls(name, SystemDependency.convert_version(raw))
        name = raw.get('name', name)
        version = SystemDependency.convert_version(raw.get('version'))
        optional = raw.get('optional', False)
        feature = raw.get('feature')
        # Everything else are overrides when certain features are enabled.
        feature_overrides = {k: v for k, v in raw.items() if k not in {'name', 'version', 'optional', 'feature'}}
        return cls(name, version, optional, feature, feature_overrides)

    @staticmethod
    def convert_version(version: T.Optional[str]) -> T.List[str]:
        vers = version.split(',') if version is not None else []
        result: T.List[str] = []
        for v in vers:
            v = v.strip()
            if v[0] not in '><=':
                v = f'>={v}'
            result.append(v)
        return result

    def enabled(self, features: T.Set[str]) -> bool:
        return self.feature is None or self.feature in features

@dataclasses.dataclass
class Dependency:

    """Representation of a Cargo Dependency Entry."""

    name: dataclasses.InitVar[str]
    version: T.List[str]
    registry: T.Optional[str] = None
    git: T.Optional[str] = None
    branch: T.Optional[str] = None
    rev: T.Optional[str] = None
    path: T.Optional[str] = None
    optional: bool = False
    package: str = ''
    default_features: bool = True
    features: T.List[str] = dataclasses.field(default_factory=list)
    api: str = dataclasses.field(init=False)

    def __post_init__(self, name: str) -> None:
        self.package = self.package or name
        # Extract wanted API version from version constraints.
        api = set()
        for v in self.version:
            if v.startswith(('>=', '==')):
                api.add(_version_to_api(v[2:].strip()))
            elif v.startswith('='):
                api.add(_version_to_api(v[1:].strip()))
        if not api:
            self.api = '0'
        elif len(api) == 1:
            self.api = api.pop()
        else:
            raise MesonException(f'Cannot determine minimum API version from {self.version}.')

    @classmethod
    def from_raw(cls, name: str, raw: manifest.DependencyV) -> Dependency:
        """Create a dependency from a raw cargo dictionary"""
        if isinstance(raw, str):
            return cls(name, version.convert(raw))
        fixed = _handle_unknown_keys(_fixup_raw_mappings(raw), cls, f'Dependency entry {name}')
        return cls(name, **fixed)


@dataclasses.dataclass
class BuildTarget:

    name: str
    crate_type: T.List[manifest.CRATE_TYPE] = dataclasses.field(default_factory=lambda: ['lib'])
    path: dataclasses.InitVar[T.Optional[str]] = None

    # https://doc.rust-lang.org/cargo/reference/cargo-targets.html#the-test-field
    # True for lib, bin, test
    test: bool = True

    # https://doc.rust-lang.org/cargo/reference/cargo-targets.html#the-doctest-field
    # True for lib
    doctest: bool = False

    # https://doc.rust-lang.org/cargo/reference/cargo-targets.html#the-bench-field
    # True for lib, bin, benchmark
    bench: bool = True

    # https://doc.rust-lang.org/cargo/reference/cargo-targets.html#the-doc-field
    # True for libraries and binaries
    doc: bool = False

    harness: bool = True
    edition: manifest.EDITION = '2015'
    required_features: T.List[str] = dataclasses.field(default_factory=list)
    plugin: bool = False

    @classmethod
    def from_raw(cls, raw: manifest.BuildTarget) -> Self:
        name = raw.get('name', '<anonymous>')
        build = _handle_unknown_keys(_fixup_raw_mappings(raw), cls, f'Binary entry {name}')
        return cls(**build)

@dataclasses.dataclass
class Library(BuildTarget):

    """Representation of a Cargo Library Entry."""

    doctest: bool = True
    doc: bool = True
    path: str = os.path.join('src', 'lib.rs')
    proc_macro: bool = False
    crate_type: T.List[manifest.CRATE_TYPE] = dataclasses.field(default_factory=lambda: ['lib'])
    doc_scrape_examples: bool = True

    @classmethod
    def from_raw(cls, raw: manifest.LibTarget, fallback_name: str) -> Self:  # type: ignore[override]
        fixed = _fixup_raw_mappings(raw)

        # We need to set the name field if it's not set manually, including if
        # other fields are set in the lib section
        if 'name' not in fixed:
            fixed['name'] = fallback_name
        fixed = _handle_unknown_keys(fixed, cls, f'Library entry {fixed["name"]}')

        return cls(**fixed)


@dataclasses.dataclass
class Binary(BuildTarget):

    """Representation of a Cargo Bin Entry."""

    doc: bool = True


@dataclasses.dataclass
class Test(BuildTarget):

    """Representation of a Cargo Test Entry."""

    bench: bool = True


@dataclasses.dataclass
class Benchmark(BuildTarget):

    """Representation of a Cargo Benchmark Entry."""

    test: bool = True


@dataclasses.dataclass
class Example(BuildTarget):

    """Representation of a Cargo Example Entry."""

    crate_type: T.List[manifest.CRATE_TYPE] = dataclasses.field(default_factory=lambda: ['bin'])


@dataclasses.dataclass
class Manifest:

    """Cargo Manifest definition.

    Most of these values map up to the Cargo Manifest, but with default values
    if not provided.

    Cargo subprojects can contain what Meson wants to treat as multiple,
    interdependent, subprojects.

    :param path: the path within the cargo subproject.
    """

    package: Package
    dependencies: T.Dict[str, Dependency]
    dev_dependencies: T.Dict[str, Dependency]
    build_dependencies: T.Dict[str, Dependency]
    system_dependencies: T.Dict[str, SystemDependency] = dataclasses.field(init=False)
    lib: Library
    bin: T.List[Binary]
    test: T.List[Test]
    bench: T.List[Benchmark]
    example: T.List[Example]
    features: T.Dict[str, T.List[str]]
    target: T.Dict[str, T.Dict[str, Dependency]]
    path: str = ''

    def __post_init__(self) -> None:
        self.features.setdefault('default', [])
        self.system_dependencies = {k: SystemDependency.from_raw(k, v) for k, v in self.package.metadata.get('system-deps', {}).items()}


def _convert_manifest(raw_manifest: manifest.Manifest, subdir: str, path: str = '') -> Manifest:
    return Manifest(
        Package.from_raw(raw_manifest['package']),
        {k: Dependency.from_raw(k, v) for k, v in raw_manifest.get('dependencies', {}).items()},
        {k: Dependency.from_raw(k, v) for k, v in raw_manifest.get('dev-dependencies', {}).items()},
        {k: Dependency.from_raw(k, v) for k, v in raw_manifest.get('build-dependencies', {}).items()},
        Library.from_raw(raw_manifest.get('lib', {}), raw_manifest['package']['name']),
        [Binary.from_raw(b) for b in raw_manifest.get('bin', {})],
        [Test.from_raw(b) for b in raw_manifest.get('test', {})],
        [Benchmark.from_raw(b) for b in raw_manifest.get('bench', {})],
        [Example.from_raw(b) for b in raw_manifest.get('example', {})],
        raw_manifest.get('features', {}),
        {k: {k2: Dependency.from_raw(k2, v2) for k2, v2 in v.get('dependencies', {}).items()}
         for k, v in raw_manifest.get('target', {}).items()},
        path,
    )


def _version_to_api(version: str) -> str:
    # x.y.z -> x
    # 0.x.y -> 0.x
    # 0.0.x -> 0
    vers = version.split('.')
    if int(vers[0]) != 0:
        return vers[0]
    elif len(vers) >= 2 and int(vers[1]) != 0:
        return f'0.{vers[1]}'
    return '0'


def _dependency_name(package_name: str, api: str) -> str:
    basename = package_name[:-3] if package_name.endswith('-rs') else package_name
    return f'{basename}-{api}-rs'


def _dependency_varname(package_name: str) -> str:
    return f'{fixup_meson_varname(package_name)}_dep'


def _extra_args_varname() -> str:
    return 'extra_args'


def _extra_deps_varname() -> str:
    return 'extra_deps'


class PackageState:
    def __init__(self, manifest: Manifest, downloaded: bool) -> None:
        self.manifest = manifest
        self.downloaded = downloaded
        self.features: T.Set[str] = set()
        self.required_deps: T.Set[str] = set()
        self.optional_deps_features: T.Dict[str, T.Set[str]] = collections.defaultdict(set)


@dataclasses.dataclass(frozen=True)
class PackageKey:
    package_name: str
    api: str


class Interpreter:
    def __init__(self, env: Environment) -> None:
        self.environment = env
        # Map Cargo.toml's subdir to loaded manifest.
        self.manifests: T.Dict[str, Manifest] = {}
        # Map of cargo package (name + api) to its state
        self.packages: T.Dict[PackageKey, PackageState] = {}

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

        # Libs are always auto-discovered and there's no other way to handle them,
        # which is unfortunate for reproducability
        if os.path.exists(os.path.join(self.environment.source_dir, subdir, pkg.manifest.path, pkg.manifest.lib.path)):
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
            filename = os.path.join(self.environment.source_dir, subdir, 'Cargo.toml')
            raw = load_toml(filename)
            if 'package' in raw:
                raw_manifest = T.cast('manifest.Manifest', raw)
                manifest_ = _convert_manifest(raw_manifest, subdir)
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

    def _create_project(self, pkg: PackageState, build: builder.Builder) -> T.List[mparser.BaseNode]:
        """Create the project() function call

        :param pkg: The package to generate from
        :param build: The AST builder
        :return: a list nodes
        """
        default_options: T.List[mparser.BaseNode] = []
        default_options.append(build.string(f'rust_std={pkg.manifest.package.edition}'))
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
        kw = {
            'version': build.array([build.string(s) for s in dep.version]),
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
            'version': build.array([build.string(s) for s in dep.version]),
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

    def _create_lib(self, pkg: PackageState, build: builder.Builder, crate_type: manifest.CRATE_TYPE) -> T.List[mparser.BaseNode]:
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
                        })
                    },
                ),
                'dep'
            ),
            build.method(
                'override_dependency',
                build.identifier('meson'),
                [
                    build.string(_dependency_name(pkg.manifest.package.name, pkg.manifest.package.api)),
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
            cargolock = T.cast('manifest.CargoLock', load_toml(filename))
        except TomlImplementationMissing as e:
            mlog.warning('Failed to load Cargo.lock:', str(e), fatal=False)
            return wraps
        for package in cargolock['package']:
            name = package['name']
            version = package['version']
            subp_name = _dependency_name(name, _version_to_api(version))
            source = package.get('source')
            if source is None:
                # This is project's package, or one of its workspace members.
                pass
            elif source == 'registry+https://github.com/rust-lang/crates.io-index':
                checksum = package.get('checksum')
                if checksum is None:
                    checksum = cargolock['metadata'][f'checksum {name} {version} ({source})']
                url = f'https://crates.io/api/v1/crates/{name}/{version}/download'
                directory = f'{name}-{version}'
                wraps.append(PackageDefinition.from_values(subp_name, subproject_dir, 'file', {
                    'directory': directory,
                    'source_url': url,
                    'source_filename': f'{directory}.tar.gz',
                    'source_hash': checksum,
                    'method': 'cargo',
                }))
            elif source.startswith('git+'):
                parts = urllib.parse.urlparse(source[4:])
                query = urllib.parse.parse_qs(parts.query)
                branch = query['branch'][0] if 'branch' in query else ''
                revision = parts.fragment or branch
                url = urllib.parse.urlunparse(parts._replace(params='', query='', fragment=''))
                wraps.append(PackageDefinition.from_values(subp_name, subproject_dir, 'git', {
                    'directory': name,
                    'url': url,
                    'revision': revision,
                    'method': 'cargo',
                }))
            else:
                mlog.warning(f'Unsupported source URL in {filename}: {source}')
    return wraps
