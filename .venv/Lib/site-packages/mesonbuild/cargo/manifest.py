# SPDX-License-Identifier: Apache-2.0
# Copyright © 2022-2024 Intel Corporation

"""Type definitions for cargo manifest files."""

from __future__ import annotations

import collections
import dataclasses
import os
import typing as T


from . import version
from ..mesonlib import MesonException, lazy_property, Version
from .. import mlog

if T.TYPE_CHECKING:
    from typing_extensions import Protocol, Self

    from . import raw
    from .raw import EDITION, CRATE_TYPE, LINT_LEVEL
    from ..wrap.wrap import PackageDefinition

    # Copied from typeshed. Blarg that they don't expose this
    class DataclassInstance(Protocol):
        __dataclass_fields__: T.ClassVar[dict[str, dataclasses.Field[T.Any]]]

_DI = T.TypeVar('_DI', bound='DataclassInstance')

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


class DefaultValue:
    """Base class to converts a raw value from cargo manifest to a meson value

    It returns the value from current manifest, or fallback to the
    workspace value. If both are None, its default value is used. Subclasses can
    override the convert() method to implement custom conversion logic.
    """

    def __init__(self, default: object = None) -> None:
        self.default = default

    def convert(self, v: T.Any, ws_v: T.Any) -> object:
        return v if v is not None else ws_v


class MergeValue(DefaultValue):
    def __init__(self, func: T.Callable[[T.Any, T.Any], object], default: object = None) -> None:
        super().__init__(default)
        self.func = func

    def convert(self, v: T.Any, ws_v: T.Any) -> object:
        return self.func(v, ws_v)


class ConvertValue(DefaultValue):
    def __init__(self, func: T.Callable[[T.Any], object], default: object = None) -> None:
        super().__init__(default)
        self.func = func

    def convert(self, v: T.Any, ws_v: T.Any) -> object:
        return self.func(v if v is not None else ws_v)


class DictMergeValue(ConvertValue):
    """Merge the incoming array of tables with a dictionary;
       a user-provided function maps each table to one of the
       entries of the dictionary."""

    def __init__(self, func: T.Callable[[T.Any], T.List[object]],
                 merge_key: T.Callable[[T.Any], str],
                 out_key: T.Callable[[T.Any], str],
                 base: T.Mapping[str, object] = None) -> None:
        super().__init__(func, base)
        self.merge_key = merge_key
        self.out_key = out_key

    def convert(self, v: T.Any, ws_v: T.Any) -> object:
        out = self.func(v if v is not None else ws_v)
        assert isinstance(out, list) # for mypy
        assert isinstance(self.default, dict) # for mypy

        explicit: T.Set[str] = set(self.merge_key(x) for x in out)
        out_d: T.Dict[str, object] = {self.out_key(x): x for x in out}
        for k, v in self.default.items():
            if self.merge_key(v) not in explicit:
                out_d[self.out_key(v)] = v
        return out_d


def _raw_to_dataclass(raw: T.Mapping[str, object], cls: T.Type[_DI], msg: str,
                      raw_from_workspace: T.Optional[T.Mapping[str, object]] = None,
                      ignored_fields: T.Optional[T.List[str]] = None,
                      **kwargs: DefaultValue) -> _DI:
    """Fixup raw cargo mappings to a dataclass.

    * Inherit values from the workspace.
    * Replaces any `-` with `_` in the keys.
    * Optionally pass values through the functions in kwargs, in order to do
      recursive conversions.
    * Remove and warn on keys that are coming from cargo, but are unknown to
      our representations.

    This is intended to give users the possibility of things proceeding when a
    new key is added to Cargo.toml that we don't yet handle, but to still warn
    them that things might not work.

    :param raw: The raw data to look at
    :param cls: The Dataclass derived type that will be created
    :param msg: the header for the error message. Usually something like "In N structure".
    :param raw_from_workspace: If inheriting from a workspace, the raw data from the workspace.
    :param kwargs: DefaultValue instances to convert values.
    :return: A @cls instance.
    """
    new_dict = {}
    unexpected = set()
    fields = {x.name for x in dataclasses.fields(cls)}
    raw_from_workspace = raw_from_workspace or {}
    ignored_fields = ignored_fields or []
    inherit = raw.get('workspace', False)

    for orig_k, v in raw.items():
        if orig_k == 'workspace':
            continue
        ws_v = None
        if isinstance(v, dict) and v.get('workspace', False):
            # foo.workspace = true, take value from workspace.
            ws_v = raw_from_workspace[orig_k]
            v = None
        elif inherit:
            # foo = {}, give the workspace value, if any, to the converter
            # function in the case it wants to merge values.
            ws_v = raw_from_workspace.get(orig_k)
        k = fixup_meson_varname(orig_k)
        if k not in fields:
            if orig_k not in ignored_fields:
                unexpected.add(orig_k)
            continue
        if k in kwargs:
            new_dict[k] = kwargs[k].convert(v, ws_v)
        else:
            new_dict[k] = v if v is not None else ws_v

    if inherit:
        # Inherit any keys from the workspace that we don't have yet.
        for orig_k, ws_v in raw_from_workspace.items():
            k = fixup_meson_varname(orig_k)
            if k not in fields:
                if orig_k not in ignored_fields:
                    unexpected.add(orig_k)
                continue
            if k in new_dict:
                continue
            if k in kwargs:
                new_dict[k] = kwargs[k].convert(None, ws_v)
            else:
                new_dict[k] = ws_v

    # Finally, set default values.
    for k, convertor in kwargs.items():
        if k not in new_dict and convertor.default is not None:
            new_dict[k] = convertor.default

    if unexpected:
        mlog.warning(msg, 'has unexpected keys', '"{}".'.format(', '.join(sorted(unexpected))),
                     _EXTRA_KEYS_WARNING)

    return cls(**new_dict)


@dataclasses.dataclass
class Package:

    """Representation of a Cargo Package entry, with defaults filled in."""

    name: str
    version: str = "0"
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
        raw_ws_pkg = workspace.package if workspace else None
        return _raw_to_dataclass(raw_pkg, cls, f'Package entry {raw_pkg["name"]}', raw_ws_pkg)


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

    @classmethod
    def from_raw(cls, name: str, raw: T.Union[T.Dict[str, T.Any], str]) -> Self:
        if isinstance(raw, str):
            raw = {'version': raw}
        name = raw.get('name', name)
        version = raw.get('version', '')
        optional = raw.get('optional', False)
        feature = raw.get('feature')
        # Everything else are overrides when certain features are enabled.
        feature_overrides = {k: v for k, v in raw.items() if k not in {'name', 'version', 'optional', 'feature'}}
        return cls(name, version, optional, feature, feature_overrides)


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
            return ''
        elif len(api) == 1:
            return api.pop()
        else:
            raise MesonException(f'Cannot determine minimum API version from {self.version}.')

    def update_version(self, v: str) -> None:
        self.version = v
        try:
            delattr(self, 'api')
        except AttributeError:
            pass
        try:
            delattr(self, 'meson_version')
        except AttributeError:
            pass

    @T.overload
    @staticmethod
    def _depv_to_dep(depv: raw.FromWorkspace) -> raw.FromWorkspace: ...

    @T.overload
    @staticmethod
    def _depv_to_dep(depv: raw.DependencyV) -> raw.Dependency: ...

    @staticmethod
    def _depv_to_dep(depv: T.Union[raw.FromWorkspace, raw.DependencyV]) -> T.Union[raw.FromWorkspace, raw.Dependency]:
        return {'version': depv} if isinstance(depv, str) else depv

    @classmethod
    def from_raw(cls, name: str, raw_depv: T.Union[raw.FromWorkspace, raw.DependencyV], member_path: str = '', workspace: T.Optional[Workspace] = None) -> Self:
        """Create a dependency from a raw cargo dictionary or string"""
        raw_ws_dep = workspace.dependencies.get(name) if workspace else None
        raw_ws_dep = cls._depv_to_dep(raw_ws_dep or {})
        raw_dep = cls._depv_to_dep(raw_depv)

        def path_convertor(path: T.Optional[str], ws_path: T.Optional[str]) -> T.Optional[str]:
            if path:
                return path
            if ws_path:
                return os.path.relpath(ws_path, member_path)
            return None

        return _raw_to_dataclass(raw_dep, cls, f'Dependency entry {name}', raw_ws_dep,
                                 package=DefaultValue(name),
                                 path=MergeValue(path_convertor),
                                 features=MergeValue(lambda features, ws_features: (features or []) + (ws_features or [])))


@dataclasses.dataclass
class BuildTarget:

    # https://doc.rust-lang.org/cargo/reference/cargo-targets.html
    # Some default values are overridden in subclasses
    name: str
    path: str
    edition: EDITION
    test: bool = True
    doctest: bool = True
    bench: bool = True
    doc: bool = True
    harness: bool = True
    crate_type: T.List[CRATE_TYPE] = dataclasses.field(default_factory=lambda: ['bin'])
    required_features: T.List[str] = dataclasses.field(default_factory=list)
    plugin: bool = False


@dataclasses.dataclass
class Library(BuildTarget):

    """Representation of a Cargo Library Entry."""

    @classmethod
    def from_raw(cls, raw: raw.LibTarget, pkg: Package) -> Self:
        name = raw.get('name', fixup_meson_varname(pkg.name))
        # If proc_macro is True, it takes precedence and sets crate_type to proc-macro
        proc_macro = raw.get('proc-macro', False)
        return _raw_to_dataclass(raw, cls, f'Library entry {name}',
                                 ignored_fields=['proc-macro'],
                                 name=DefaultValue(name),
                                 path=DefaultValue('src/lib.rs'),
                                 edition=DefaultValue(pkg.edition),
                                 crate_type=ConvertValue(lambda x: ['proc-macro'] if proc_macro else x,
                                                         ['proc-macro'] if proc_macro else ['lib']))


@dataclasses.dataclass
class Binary(BuildTarget):

    """Representation of a Cargo Bin Entry."""

    @classmethod
    def from_raw(cls, raw: raw.BuildTarget, pkg: Package) -> Self:
        name = raw["name"]
        return _raw_to_dataclass(raw, cls, f'Binary entry {name}',
                                 path=DefaultValue('src/main.rs'),
                                 edition=DefaultValue(pkg.edition))


@dataclasses.dataclass
class Test(BuildTarget):

    """Representation of a Cargo Test Entry."""

    @classmethod
    def from_raw(cls, raw: raw.BuildTarget, pkg: Package) -> Self:
        name = raw["name"]
        return _raw_to_dataclass(raw, cls, f'Test entry {name}',
                                 path=DefaultValue(f'tests/{name}.rs'),
                                 edition=DefaultValue(pkg.edition),
                                 bench=DefaultValue(False),
                                 doc=DefaultValue(False))


@dataclasses.dataclass
class Benchmark(BuildTarget):

    """Representation of a Cargo Benchmark Entry."""

    @classmethod
    def from_raw(cls, raw: raw.BuildTarget, pkg: Package) -> Self:
        name = raw["name"]
        return _raw_to_dataclass(raw, cls, f'Benchmark entry {name}',
                                 path=DefaultValue(f'benches/{name}.rs'),
                                 edition=DefaultValue(pkg.edition),
                                 test=DefaultValue(False),
                                 doc=DefaultValue(False))


@dataclasses.dataclass
class Example(BuildTarget):

    """Representation of a Cargo Example Entry."""

    @classmethod
    def from_raw(cls, raw: raw.BuildTarget, pkg: Package) -> Self:
        name = raw["name"]
        return _raw_to_dataclass(raw, cls, f'Example entry {name}',
                                 path=DefaultValue(f'examples/{name}.rs'),
                                 edition=DefaultValue(pkg.edition),
                                 test=DefaultValue(False),
                                 bench=DefaultValue(False),
                                 doc=DefaultValue(False))


@dataclasses.dataclass
class Lint:

    """Cargo Lint definition.
    """

    name: str
    level: LINT_LEVEL
    priority: int
    check_cfg: T.Optional[T.List[str]]

    @classmethod
    def from_raw(cls, r: T.Union[raw.FromWorkspace, T.Dict[str, T.Dict[str, raw.LintV]]]) -> T.List[Lint]:
        r = T.cast('T.Dict[str, T.Dict[str, raw.LintV]]', r)
        lints: T.Dict[str, Lint] = {}
        for tool, raw_lints in r.items():
            prefix = '' if tool == 'rust' else f'{tool}::'
            for name, settings in raw_lints.items():
                name = prefix + name
                if isinstance(settings, str):
                    settings = T.cast('raw.Lint', {'level': settings})
                check_cfg = None
                if name == 'unexpected_cfgs':
                    check_cfg = settings.get('check-cfg', [])
                lints[name] = Lint(name=name,
                                   level=settings['level'],
                                   priority=settings.get('priority', 0),
                                   check_cfg=check_cfg)

        lints_final = list(lints.values())
        lints_final.sort(key=lambda x: x.priority)
        return lints_final

    def to_arguments(self, check_cfg: bool) -> T.List[str]:
        if self.level == "deny":
            flag = "-D"
        elif self.level == "allow":
            flag = "-A"
        elif self.level == "warn":
            flag = "-W"
        elif self.level == "forbid":
            flag = "-F"
        else:
            raise MesonException(f"invalid level {self.level!r} for {self.name}")
        args = [flag, self.name]
        if check_cfg and self.check_cfg:
            for arg in self.check_cfg:
                args.append('--check-cfg')
                args.append(arg)
        return args


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
    bin: T.Dict[str, Binary] = dataclasses.field(default_factory=dict)
    test: T.List[Test] = dataclasses.field(default_factory=list)
    bench: T.List[Benchmark] = dataclasses.field(default_factory=list)
    example: T.List[Example] = dataclasses.field(default_factory=list)
    features: T.Dict[str, T.List[str]] = dataclasses.field(default_factory=dict)
    target: T.Dict[str, T.Dict[str, Dependency]] = dataclasses.field(default_factory=dict)
    lints: T.List[Lint] = dataclasses.field(default_factory=list)

    # missing: profile

    def __post_init__(self) -> None:
        self.features.setdefault('default', [])

    @lazy_property
    def system_dependencies(self) -> T.Dict[str, SystemDependency]:
        return {k: SystemDependency.from_raw(k, v) for k, v in self.package.metadata.get('system-deps', {}).items()}

    @classmethod
    def from_raw(cls, raw: raw.Manifest, path: str, workspace: T.Optional[Workspace] = None, member_path: str = '') -> Self:
        pkg = Package.from_raw(raw['package'], workspace)

        autolib = None
        if pkg.autolib and os.path.exists(os.path.join(path, 'src/lib.rs')):
            autolib = Library.from_raw({}, pkg)

        def _discover_targets(subdir: str) -> T.Generator[T.Tuple[str, str], None, None]:
            """Discover .rs files in a subdirectory and yield (name, path) tuples."""
            target_dir = os.path.join(path, subdir)
            if os.path.isdir(target_dir):
                for entry in os.listdir(target_dir):
                    if entry.endswith('.rs'):
                        target_name = entry[:-3]  # Remove .rs extension
                        yield target_name, f'{subdir}/{entry}'

        autobins: T.Dict[str, Binary] = {}
        if pkg.autobins:
            # Check for default binary (src/main.rs)
            if os.path.exists(os.path.join(path, 'src/main.rs')):
                autobins[pkg.name] = Binary.from_raw({'name': pkg.name, 'path': 'src/main.rs'}, pkg)
            # Add additional binaries from src/bin/
            for bin_name, bin_path in _discover_targets('src/bin'):
                if bin_name in autobins:
                    raise MesonException(f'Binary target {bin_name!r} is defined more than once '
                                         f'({autobins[bin_name].path} and {bin_path})')
                autobins[bin_name] = Binary.from_raw({'name': bin_name, 'path': bin_path}, pkg)

        def dependencies_from_raw(x: T.Dict[str, T.Any]) -> T.Dict[str, Dependency]:
            return {k: Dependency.from_raw(k, v, member_path, workspace) for k, v in x.items()}

        return _raw_to_dataclass(raw, cls, f'Cargo.toml package {pkg.name}',
                                 raw_from_workspace=workspace.inheritable if workspace else None,
                                 ignored_fields=['badges', 'workspace'],
                                 package=ConvertValue(lambda _: pkg),
                                 dependencies=ConvertValue(dependencies_from_raw),
                                 dev_dependencies=ConvertValue(dependencies_from_raw),
                                 build_dependencies=ConvertValue(dependencies_from_raw),
                                 lints=ConvertValue(Lint.from_raw),
                                 lib=ConvertValue(lambda x: Library.from_raw(x, pkg), default=autolib),
                                 bin=DictMergeValue(lambda x: [Binary.from_raw(b, pkg) for b in x],
                                                    merge_key=lambda x: x.path,
                                                    out_key=lambda x: x.name,
                                                    base=autobins),
                                 test=ConvertValue(lambda x: [Test.from_raw(b, pkg) for b in x]),
                                 bench=ConvertValue(lambda x: [Benchmark.from_raw(b, pkg) for b in x]),
                                 example=ConvertValue(lambda x: [Example.from_raw(b, pkg) for b in x]),
                                 target=ConvertValue(lambda x: {k: dependencies_from_raw(v.get('dependencies', {})) for k, v in x.items()}))


@dataclasses.dataclass
class Workspace:

    """Cargo Workspace definition.
    """

    resolver: str = dataclasses.field(default_factory=lambda: '2')
    members: T.List[str] = dataclasses.field(default_factory=list)
    exclude: T.List[str] = dataclasses.field(default_factory=list)
    default_members: T.List[str] = dataclasses.field(default_factory=list)

    # inheritable settings are kept in raw format, for use with _raw_to_dataclass
    package: T.Optional[raw.Package] = None
    dependencies: T.Dict[str, raw.Dependency] = dataclasses.field(default_factory=dict)
    lints: T.Dict[str, T.Dict[str, raw.LintV]] = dataclasses.field(default_factory=dict)
    metadata: T.Dict[str, T.Any] = dataclasses.field(default_factory=dict)

    # A workspace can also have a root package.
    root_package: T.Optional[Manifest] = None

    @lazy_property
    def inheritable(self) -> T.Dict[str, object]:
        # the whole lints table is inherited.  Do not add package, dependencies
        # etc. because they can only be inherited a field at a time.
        return {
            'lints': self.lints,
        }

    @classmethod
    def from_raw(cls, raw: raw.Manifest, path: str) -> Self:
        ws = _raw_to_dataclass(raw['workspace'], cls, 'Workspace')
        if 'package' in raw:
            ws.root_package = Manifest.from_raw(raw, path, ws, '.')
        if not ws.default_members:
            ws.default_members = ['.'] if ws.root_package else ws.members
        return ws


@dataclasses.dataclass
class CargoLockPackage:

    """A description of a package in the Cargo.lock file format."""

    name: str
    version: str
    source: T.Optional[str] = None
    checksum: T.Optional[str] = None
    dependencies: T.List[str] = dataclasses.field(default_factory=list)

    @lazy_property
    def api(self) -> str:
        return version.api(self.version)

    @lazy_property
    def subproject(self) -> str:
        return f'{self.name}-{self.api}-rs'

    @classmethod
    def from_raw(cls, raw: raw.CargoLockPackage) -> Self:
        return _raw_to_dataclass(raw, cls, 'Cargo.lock package')


@dataclasses.dataclass
class CargoLock:

    """A description of the Cargo.lock file format."""

    version: int = 1
    package: T.List[CargoLockPackage] = dataclasses.field(default_factory=list)
    metadata: T.Dict[str, str] = dataclasses.field(default_factory=dict)
    wraps: T.Dict[str, PackageDefinition] = dataclasses.field(default_factory=dict)

    def named(self, name: str) -> T.Sequence[CargoLockPackage]:
        return self._versions[name]

    @lazy_property
    def _versions(self) -> T.Dict[str, T.List[CargoLockPackage]]:
        versions = collections.defaultdict(list)
        for pkg in self.package:
            versions[pkg.name].append(pkg)
        for pkg_versions in versions.values():
            pkg_versions.sort(reverse=True, key=lambda pkg: Version(pkg.version))
        return versions

    @classmethod
    def from_raw(cls, raw: raw.CargoLock) -> Self:
        return _raw_to_dataclass(raw, cls, 'Cargo.lock',
                                 package=ConvertValue(lambda x: [CargoLockPackage.from_raw(p) for p in x]))
