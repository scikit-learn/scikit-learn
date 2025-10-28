# SPDX-License-Identifier: Apache-2.0
# Copyright © 2022-2024 Intel Corporation

"""Type definitions for cargo manifest files."""

from __future__ import annotations

import dataclasses
import os
import typing as T

from . import version
from ..mesonlib import MesonException, lazy_property
from .. import mlog

if T.TYPE_CHECKING:
    from typing_extensions import Protocol, Self

    from . import raw
    from .raw import EDITION, CRATE_TYPE

    # Copied from typeshed. Blarg that they don't expose this
    class DataclassInstance(Protocol):
        __dataclass_fields__: T.ClassVar[dict[str, dataclasses.Field[T.Any]]]

_DI = T.TypeVar('_DI', bound='DataclassInstance')
_R = T.TypeVar('_R', bound='raw._BaseBuildTarget')

_EXTRA_KEYS_WARNING = (
    "This may (unlikely) be an error in the cargo manifest, or may be a missing "
    "implementation in Meson. If this issue can be reproduced with the latest "
    "version of Meson, please help us by opening an issue at "
    "https://github.com/mesonbuild/meson/issues. Please include the crate and "
    "version that is generating this warning if possible."
)


def fixup_meson_varname(name: str) -> str:
    """Fixup a meson variable name

    :param name: The name to fix
    :return: the fixed name
    """
    return name.replace('-', '_')


@T.overload
def _depv_to_dep(depv: raw.FromWorkspace) -> raw.FromWorkspace: ...

@T.overload
def _depv_to_dep(depv: raw.DependencyV) -> raw.Dependency: ...

def _depv_to_dep(depv: T.Union[raw.FromWorkspace, raw.DependencyV]) -> T.Union[raw.FromWorkspace, raw.Dependency]:
    return {'version': depv} if isinstance(depv, str) else depv


def _raw_to_dataclass(raw: T.Mapping[str, object], cls: T.Type[_DI],
                      msg: str, **kwargs: T.Callable[[T.Any], object]) -> _DI:
    """Fixup raw cargo mappings to ones more suitable for python to consume as dataclass.

    * Replaces any `-` with `_` in the keys.
    * Optionally pass values through the functions in kwargs, in order to do
      recursive conversions.
    * Remove and warn on keys that are coming from cargo, but are unknown to
      our representations.

    This is intended to give users the possibility of things proceeding when a
    new key is added to Cargo.toml that we don't yet handle, but to still warn
    them that things might not work.

    :param data: The raw data to look at
    :param cls: The Dataclass derived type that will be created
    :param msg: the header for the error message. Usually something like "In N structure".
    :return: The original data structure, but with all unknown keys removed.
    """
    new_dict = {}
    unexpected = set()
    fields = {x.name for x in dataclasses.fields(cls)}
    for orig_k, v in raw.items():
        k = fixup_meson_varname(orig_k)
        if k not in fields:
            unexpected.add(orig_k)
            continue
        if k in kwargs:
            v = kwargs[k](v)
        new_dict[k] = v

    if unexpected:
        mlog.warning(msg, 'has unexpected keys', '"{}".'.format(', '.join(sorted(unexpected))),
                     _EXTRA_KEYS_WARNING)
    return cls(**new_dict)


@T.overload
def _inherit_from_workspace(raw: raw.Package,
                            raw_from_workspace: T.Optional[T.Mapping[str, object]],
                            msg: str,
                            **kwargs: T.Callable[[T.Any, T.Any], object]) -> raw.Package: ...

@T.overload
def _inherit_from_workspace(raw: T.Union[raw.FromWorkspace, raw.Dependency],
                            raw_from_workspace: T.Optional[T.Mapping[str, object]],
                            msg: str,
                            **kwargs: T.Callable[[T.Any, T.Any], object]) -> raw.Dependency: ...

def _inherit_from_workspace(raw_: T.Union[raw.FromWorkspace, raw.Package, raw.Dependency], # type: ignore[misc]
                            raw_from_workspace: T.Optional[T.Mapping[str, object]],
                            msg: str,
                            **kwargs: T.Callable[[T.Any, T.Any], object]) -> T.Mapping[str, object]:
    # allow accesses by non-literal key below
    raw = T.cast('T.Mapping[str, object]', raw_)

    if not raw_from_workspace:
        if raw.get('workspace', False) or \
                any(isinstance(v, dict) and v.get('workspace', False) for v in raw):
            raise MesonException(f'Cargo.toml file requests {msg} from workspace')

        return raw

    result = {k: v for k, v in raw.items() if k != 'workspace'}
    for k, v in raw.items():
        if isinstance(v, dict) and v.get('workspace', False):
            if k in raw_from_workspace:
                result[k] = raw_from_workspace[k]
                if k in kwargs:
                    result[k] = kwargs[k](v, result[k])
            else:
                del result[k]

    if raw.get('workspace', False):
        for k, v in raw_from_workspace.items():
            if k not in result or k in kwargs:
                if k in kwargs:
                    v = kwargs[k](raw.get(k), v)
                result[k] = v
    return result


@dataclasses.dataclass
class Package:

    """Representation of a Cargo Package entry, with defaults filled in."""

    name: str
    version: str
    description: T.Optional[str] = None
    resolver: T.Optional[str] = None
    authors: T.List[str] = dataclasses.field(default_factory=list)
    edition: EDITION = '2015'
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

    @lazy_property
    def api(self) -> str:
        return version.api(self.version)

    @classmethod
    def from_raw(cls, raw_pkg: raw.Package, workspace: T.Optional[Workspace] = None) -> Self:
        raw_ws_pkg = None
        if workspace is not None:
            raw_ws_pkg = workspace.package

        raw_pkg = _inherit_from_workspace(raw_pkg, raw_ws_pkg, f'Package entry {raw_pkg["name"]}')
        return _raw_to_dataclass(raw_pkg, cls, f'Package entry {raw_pkg["name"]}')

@dataclasses.dataclass
class SystemDependency:

    """ Representation of a Cargo system-deps entry
        https://docs.rs/system-deps/latest/system_deps
    """

    name: str
    version: str = ''
    optional: bool = False
    feature: T.Optional[str] = None
    # TODO: convert values to dataclass
    feature_overrides: T.Dict[str, T.Dict[str, str]] = dataclasses.field(default_factory=dict)

    @classmethod
    def from_raw(cls, name: str, raw: T.Union[T.Dict[str, T.Any], str]) -> SystemDependency:
        if isinstance(raw, str):
            raw = {'version': raw}
        name = raw.get('name', name)
        version = raw.get('version', '')
        optional = raw.get('optional', False)
        feature = raw.get('feature')
        # Everything else are overrides when certain features are enabled.
        feature_overrides = {k: v for k, v in raw.items() if k not in {'name', 'version', 'optional', 'feature'}}
        return cls(name, version, optional, feature, feature_overrides)

    @lazy_property
    def meson_version(self) -> T.List[str]:
        vers = self.version.split(',') if self.version else []
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

    package: str
    version: str = ''
    registry: T.Optional[str] = None
    git: T.Optional[str] = None
    branch: T.Optional[str] = None
    rev: T.Optional[str] = None
    path: T.Optional[str] = None
    optional: bool = False
    default_features: bool = True
    features: T.List[str] = dataclasses.field(default_factory=list)

    @lazy_property
    def meson_version(self) -> T.List[str]:
        return version.convert(self.version)

    @lazy_property
    def api(self) -> str:
        # Extract wanted API version from version constraints.
        api = set()
        for v in self.meson_version:
            if v.startswith(('>=', '==')):
                api.add(version.api(v[2:].strip()))
            elif v.startswith('='):
                api.add(version.api(v[1:].strip()))
        if not api:
            return '0'
        elif len(api) == 1:
            return api.pop()
        else:
            raise MesonException(f'Cannot determine minimum API version from {self.version}.')

    @classmethod
    def from_raw_dict(cls, name: str, raw_dep: T.Union[raw.FromWorkspace, raw.Dependency], member_path: str = '', raw_ws_dep: T.Optional[raw.Dependency] = None) -> Dependency:
        raw_dep = _inherit_from_workspace(raw_dep, raw_ws_dep,
                                          f'Dependency entry {name}',
                                          path=lambda pkg_path, ws_path: os.path.relpath(ws_path, member_path),
                                          features=lambda pkg_path, ws_path: (pkg_path or []) + (ws_path or []))
        raw_dep.setdefault('package', name)
        return _raw_to_dataclass(raw_dep, cls, f'Dependency entry {name}')

    @classmethod
    def from_raw(cls, name: str, raw_depv: T.Union[raw.FromWorkspace, raw.DependencyV], member_path: str = '', workspace: T.Optional[Workspace] = None) -> Dependency:
        """Create a dependency from a raw cargo dictionary or string"""
        raw_ws_dep: T.Optional[raw.Dependency] = None
        if workspace is not None:
            raw_ws_depv = workspace.dependencies.get(name, {})
            raw_ws_dep = _depv_to_dep(raw_ws_depv)

        raw_dep = _depv_to_dep(raw_depv)
        return cls.from_raw_dict(name, raw_dep, member_path, raw_ws_dep)


@dataclasses.dataclass
class BuildTarget(T.Generic[_R]):

    name: str
    path: str
    crate_type: T.List[CRATE_TYPE]

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
    edition: EDITION = '2015'
    required_features: T.List[str] = dataclasses.field(default_factory=list)
    plugin: bool = False

    @classmethod
    def from_raw(cls, raw: _R) -> Self:
        name = raw.get('name', '<anonymous>')
        return _raw_to_dataclass(raw, cls, f'Binary entry {name}')

@dataclasses.dataclass
class Library(BuildTarget['raw.LibTarget']):

    """Representation of a Cargo Library Entry."""

    doctest: bool = True
    doc: bool = True
    path: str = os.path.join('src', 'lib.rs')
    proc_macro: bool = False
    crate_type: T.List[CRATE_TYPE] = dataclasses.field(default_factory=lambda: ['lib'])
    doc_scrape_examples: bool = True

    @classmethod
    def from_raw(cls, raw: raw.LibTarget, fallback_name: str) -> Self:  # type: ignore[override]
        # We need to set the name field if it's not set manually, including if
        # other fields are set in the lib section
        raw.setdefault('name', fallback_name)
        return _raw_to_dataclass(raw, cls, f'Library entry {raw["name"]}')


@dataclasses.dataclass
class Binary(BuildTarget['raw.BuildTarget']):

    """Representation of a Cargo Bin Entry."""

    doc: bool = True
    crate_type: T.List[CRATE_TYPE] = dataclasses.field(default_factory=lambda: ['bin'])

    @classmethod
    def from_raw(cls, raw: raw.BuildTarget) -> Self:
        if 'path' not in raw:
            raw['path'] = os.path.join('bin', raw['name'] + '.rs')
        return super().from_raw(raw)


@dataclasses.dataclass
class Test(BuildTarget['raw.BuildTarget']):

    """Representation of a Cargo Test Entry."""

    bench: bool = True
    crate_type: T.List[CRATE_TYPE] = dataclasses.field(default_factory=lambda: ['bin'])

    @classmethod
    def from_raw(cls, raw: raw.BuildTarget) -> Self:
        if 'path' not in raw:
            raw['path'] = os.path.join('tests', raw['name'] + '.rs')
        return super().from_raw(raw)

@dataclasses.dataclass
class Benchmark(BuildTarget['raw.BuildTarget']):

    """Representation of a Cargo Benchmark Entry."""

    test: bool = True
    crate_type: T.List[CRATE_TYPE] = dataclasses.field(default_factory=lambda: ['bin'])

    @classmethod
    def from_raw(cls, raw: raw.BuildTarget) -> Self:
        if 'path' not in raw:
            raw['path'] = os.path.join('benches', raw['name'] + '.rs')
        return super().from_raw(raw)


@dataclasses.dataclass
class Example(BuildTarget['raw.BuildTarget']):

    """Representation of a Cargo Example Entry."""

    crate_type: T.List[CRATE_TYPE] = dataclasses.field(default_factory=lambda: ['bin'])

    @classmethod
    def from_raw(cls, raw: raw.BuildTarget) -> Self:
        if 'path' not in raw:
            raw['path'] = os.path.join('examples', raw['name'] + '.rs')
        return super().from_raw(raw)


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
    dependencies: T.Dict[str, Dependency] = dataclasses.field(default_factory=dict)
    dev_dependencies: T.Dict[str, Dependency] = dataclasses.field(default_factory=dict)
    build_dependencies: T.Dict[str, Dependency] = dataclasses.field(default_factory=dict)
    lib: T.Optional[Library] = None
    bin: T.List[Binary] = dataclasses.field(default_factory=list)
    test: T.List[Test] = dataclasses.field(default_factory=list)
    bench: T.List[Benchmark] = dataclasses.field(default_factory=list)
    example: T.List[Example] = dataclasses.field(default_factory=list)
    features: T.Dict[str, T.List[str]] = dataclasses.field(default_factory=dict)
    target: T.Dict[str, T.Dict[str, Dependency]] = dataclasses.field(default_factory=dict)

    path: str = ''

    def __post_init__(self) -> None:
        self.features.setdefault('default', [])

    @lazy_property
    def system_dependencies(self) -> T.Dict[str, SystemDependency]:
        return {k: SystemDependency.from_raw(k, v) for k, v in self.package.metadata.get('system-deps', {}).items()}

    @classmethod
    def from_raw(cls, raw: raw.Manifest, path: str = '', workspace: T.Optional[Workspace] = None, member_path: str = '') -> Self:
        # Libs are always auto-discovered and there's no other way to handle them,
        # which is unfortunate for reproducability
        pkg = Package.from_raw(raw['package'], workspace)
        if pkg.autolib and 'lib' not in raw and \
                os.path.exists(os.path.join(path, 'src/lib.rs')):
            raw['lib'] = {}
        fixed = _raw_to_dataclass(raw, cls, f'Cargo.toml package {raw["package"]["name"]}',
                                  package=lambda x: pkg,
                                  dependencies=lambda x: {k: Dependency.from_raw(k, v, member_path, workspace) for k, v in x.items()},
                                  dev_dependencies=lambda x: {k: Dependency.from_raw(k, v, member_path, workspace) for k, v in x.items()},
                                  build_dependencies=lambda x: {k: Dependency.from_raw(k, v, member_path, workspace) for k, v in x.items()},
                                  lib=lambda x: Library.from_raw(x, raw['package']['name']),
                                  bin=lambda x: [Binary.from_raw(b) for b in x],
                                  test=lambda x: [Test.from_raw(b) for b in x],
                                  bench=lambda x: [Benchmark.from_raw(b) for b in x],
                                  example=lambda x: [Example.from_raw(b) for b in x],
                                  target=lambda x: {k: {k2: Dependency.from_raw(k2, v2, member_path, workspace) for k2, v2 in v.get('dependencies', {}).items()}
                                                    for k, v in x.items()})
        fixed.path = path
        return fixed


@dataclasses.dataclass
class Workspace:

    """Cargo Workspace definition.
    """

    resolver: str = dataclasses.field(default_factory=lambda: '2')
    members: T.List[str] = dataclasses.field(default_factory=list)
    exclude: T.List[str] = dataclasses.field(default_factory=list)
    default_members: T.List[str] = dataclasses.field(default_factory=list)

    # inheritable settings are kept in raw format, for use with _inherit_from_workspace
    package: T.Optional[raw.Package] = None
    dependencies: T.Dict[str, raw.Dependency] = dataclasses.field(default_factory=dict)
    lints: T.Dict[str, T.Any] = dataclasses.field(default_factory=dict)
    metadata: T.Dict[str, T.Any] = dataclasses.field(default_factory=dict)

    # A workspace can also have a root package.
    root_package: T.Optional[Manifest] = dataclasses.field(init=False)

    @classmethod
    def from_raw(cls, raw: raw.VirtualManifest) -> Workspace:
        ws_raw = raw['workspace']
        fixed = _raw_to_dataclass(ws_raw, cls, 'Workspace')
        return fixed


@dataclasses.dataclass
class CargoLockPackage:

    """A description of a package in the Cargo.lock file format."""

    name: str
    version: str
    source: T.Optional[str] = None
    checksum: T.Optional[str] = None
    dependencies: T.List[str] = dataclasses.field(default_factory=list)

    @classmethod
    def from_raw(cls, raw: raw.CargoLockPackage) -> CargoLockPackage:
        return _raw_to_dataclass(raw, cls, 'Cargo.lock package')


@dataclasses.dataclass
class CargoLock:

    """A description of the Cargo.lock file format."""

    version: int = 1
    package: T.List[CargoLockPackage] = dataclasses.field(default_factory=list)
    metadata: T.Dict[str, str] = dataclasses.field(default_factory=dict)

    @classmethod
    def from_raw(cls, raw: raw.CargoLock) -> CargoLock:
        return _raw_to_dataclass(raw, cls, 'Cargo.lock',
                                 package=lambda x: [CargoLockPackage.from_raw(p) for p in x])
