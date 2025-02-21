# SPDX-FileCopyrightText: 2021 Filipe Laíns <lains@riseup.net>
# SPDX-FileCopyrightText: 2021 Quansight, LLC
# SPDX-FileCopyrightText: 2022 The meson-python developers
#
# SPDX-License-Identifier: MIT

"""Meson Python build backend

Implements PEP 517 hooks.
"""

from __future__ import annotations

import argparse
import collections
import contextlib
import difflib
import functools
import importlib.machinery
import io
import itertools
import json
import os
import pathlib
import platform
import re
import shutil
import subprocess
import sys
import sysconfig
import tarfile
import tempfile
import textwrap
import typing
import warnings


if sys.version_info < (3, 11):
    import tomli as tomllib
else:
    import tomllib

import packaging.utils
import packaging.version
import pyproject_metadata

import mesonpy._compat
import mesonpy._rpath
import mesonpy._tags
import mesonpy._util
import mesonpy._wheelfile

from mesonpy._compat import cached_property, read_binary


if typing.TYPE_CHECKING:  # pragma: no cover
    from typing import Any, Callable, DefaultDict, Dict, List, Literal, Optional, Sequence, TextIO, Tuple, Type, TypeVar, Union

    from mesonpy._compat import Collection, Iterator, Mapping, ParamSpec, Path, Self

    P = ParamSpec('P')
    T = TypeVar('T')

    MesonArgsKeys = Literal['dist', 'setup', 'compile', 'install']
    MesonArgs = Mapping[MesonArgsKeys, List[str]]


__version__ = '0.17.1'


_NINJA_REQUIRED_VERSION = '1.8.2'
_MESON_REQUIRED_VERSION = '0.63.3' # keep in sync with the version requirement in pyproject.toml

_MESON_ARGS_KEYS = ['dist', 'setup', 'compile', 'install']

_SUFFIXES = importlib.machinery.all_suffixes()
_EXTENSION_SUFFIX_REGEX = re.compile(r'^[^.]+\.(?:(?P<abi>[^.]+)\.)?(?:so|pyd|dll)$')
assert all(re.match(_EXTENSION_SUFFIX_REGEX, f'foo{x}') for x in importlib.machinery.EXTENSION_SUFFIXES)

# Map Meson installation path placeholders to wheel installation paths.
# See https://docs.python.org/3/library/sysconfig.html#installation-paths
_INSTALLATION_PATH_MAP = {
    '{bindir}': 'scripts',
    '{py_purelib}': 'purelib',
    '{py_platlib}': 'platlib',
    '{moduledir_shared}': 'platlib',
    '{includedir}': 'headers',
    '{datadir}': 'data',
    # custom location
    '{libdir}': 'mesonpy-libs',
    '{libdir_shared}': 'mesonpy-libs',
}


def _map_to_wheel(sources: Dict[str, Dict[str, Any]]) -> DefaultDict[str, List[Tuple[pathlib.Path, str]]]:
    """Map files to the wheel, organized by wheel installation directory."""
    wheel_files: DefaultDict[str, List[Tuple[pathlib.Path, str]]] = collections.defaultdict(list)
    packages: Dict[str, str] = {}

    for key, group in sources.items():
        for src, target in group.items():
            destination = pathlib.Path(target['destination'])
            anchor = destination.parts[0]
            dst = pathlib.Path(*destination.parts[1:])

            path = _INSTALLATION_PATH_MAP.get(anchor)
            if path is None:
                raise BuildError(f'Could not map installation path to an equivalent wheel directory: {str(destination)!r}')

            if path == 'purelib' or path == 'platlib':
                package = destination.parts[1]
                other = packages.setdefault(package, path)
                if other != path:
                    this = os.fspath(pathlib.Path(path, *destination.parts[1:]))
                    that = os.fspath(other / next(d for d, s in wheel_files[other] if d.parts[0] == destination.parts[1]))
                    raise BuildError(
                        f'The {package} package is split between {path} and {other}: '
                        f'{this!r} and {that!r}, a "pure: false" argument may be missing in meson.build. '
                        f'It is recommended to set it in "import(\'python\').find_installation()"')

            if key == 'install_subdirs' or key == 'targets' and os.path.isdir(src):
                exclude_files = {os.path.normpath(x) for x in target.get('exclude_files', [])}
                exclude_dirs = {os.path.normpath(x) for x in target.get('exclude_dirs', [])}
                for root, dirnames, filenames in os.walk(src):
                    for name in dirnames.copy():
                        dirsrc = os.path.join(root, name)
                        relpath = os.path.relpath(dirsrc, src)
                        if relpath in exclude_dirs:
                            dirnames.remove(name)
                    # sort to process directories determninistically
                    dirnames.sort()
                    for name in sorted(filenames):
                        filesrc = os.path.join(root, name)
                        relpath = os.path.relpath(filesrc, src)
                        if relpath in exclude_files:
                            continue
                        filedst = dst / relpath
                        wheel_files[path].append((filedst, filesrc))
            else:
                wheel_files[path].append((dst, src))

    return wheel_files


class style:
    ERROR = '\33[31m'  # red
    WARNING = '\33[93m'  # bright yellow
    INFO = '\33[36m\33[1m'  # cyan, bold
    RESET = '\33[0m'

    @staticmethod
    def strip(string: str) -> str:
        """Strip ANSI escape sequences from string."""
        return re.sub(r'\033\[[;?0-9]*[a-zA-Z]', '', string)


@functools.lru_cache()
def _use_ansi_escapes() -> bool:
    """Determine whether logging should use ANSI escapes."""

    # We print log messages and error messages that may contain file
    # names containing characters that cannot be represented in the
    # stdout encoding. Use replacement markers for those instead than
    # raising UnicodeEncodeError.
    sys.stdout.reconfigure(errors='replace')  # type: ignore[attr-defined]

    if 'NO_COLOR' in os.environ:
        return False

    if 'FORCE_COLOR' in os.environ or sys.stdout.isatty() and os.environ.get('TERM') != 'dumb':
        if sys.platform == 'win32' and not os.environ.get('ANSICON'):
            return mesonpy._util.setup_windows_console()
        return True

    return False


def _log(string: str , **kwargs: Any) -> None:
    if not _use_ansi_escapes():
        string = style.strip(string)
    print(string, **kwargs)


def _showwarning(
    message: Union[Warning, str],
    category: Type[Warning],
    filename: str,
    lineno: int,
    file: Optional[TextIO] = None,
    line: Optional[str] = None,
) -> None:  # pragma: no cover
    """Callable to override the default warning handler, to have colored output."""
    _log(f'{style.WARNING}meson-python: warning:{style.RESET} {message}')


class _clicounter:
    def __init__(self, total: int) -> None:
        self._total = total
        self._count = itertools.count(start=1)

    def __enter__(self) -> Self:
        return self

    def update(self, description: str) -> None:
        line = f'[{next(self._count)}/{self._total}] {description}'
        if _use_ansi_escapes():
            print('\r', line, sep='', end='\33[0K', flush=True)
        else:
            print(line)

    def __exit__(self, exc_type: Any, exc_value: Any, exc_tb: Any) -> None:
        if _use_ansi_escapes():
            print()


class Error(RuntimeError):
    def __str__(self) -> str:
        return str(self.args[0])


class ConfigError(Error):
    """Error in the backend configuration."""


class BuildError(Error):
    """Error when building the wheel."""


class MesonBuilderError(Error):
    """Error when building the Meson package."""


class Metadata(pyproject_metadata.StandardMetadata):
    def __init__(self, name: str, *args: Any, **kwargs: Any):
        super().__init__(name, *args, **kwargs)
        # Local fix for https://github.com/FFY00/python-pyproject-metadata/issues/60
        self.name = self._validate_name(name)

    @staticmethod
    def _validate_name(name: str) -> str:
        # See https://packaging.python.org/en/latest/specifications/core-metadata/#name
        if not re.match(r'^([A-Z0-9]|[A-Z0-9][A-Z0-9._-]*[A-Z0-9])$', name, re.IGNORECASE):
            raise pyproject_metadata.ConfigurationError(
                f'Invalid project name "{name}". A valid name consists only of ASCII letters and '
                f'numbers, period, underscore and hyphen. It must start and end with a letter or number')
        return name

    @classmethod
    def from_pyproject(  # type: ignore[override]
        cls,
        data: Mapping[str, Any],
        project_dir: Path = os.path.curdir,
        metadata_version: Optional[str] = None
    ) -> Self:
        metadata = super().from_pyproject(data, project_dir)

        # Check for missing version field.
        if not metadata.version and 'version' not in metadata.dynamic:
            raise pyproject_metadata.ConfigurationError(
                'Required "project.version" field is missing and not declared as dynamic')

        # Check for unsupported dynamic fields.
        unsupported_dynamic = set(metadata.dynamic) - {'version', }
        if unsupported_dynamic:
            fields = ', '.join(f'"{x}"' for x in unsupported_dynamic)
            raise pyproject_metadata.ConfigurationError(f'Unsupported dynamic fields: {fields}')

        return metadata

    # Local fix for a bug in pyproject-metadata. See
    # https://github.com/mesonbuild/meson-python/issues/454
    def _update_dynamic(self, value: Any) -> None:
        if value and 'version' in self.dynamic:
            self.dynamic.remove('version')

    @property
    def canonical_name(self) -> str:
        # See https://packaging.python.org/en/latest/specifications/name-normalization/#normalization
        return packaging.utils.canonicalize_name(self.name)

    @property
    def distribution_name(self) -> str:
        """Name to be used in wheel and sdist file names."""
        return self.canonical_name.replace('-', '_')


def _is_native(file: Path) -> bool:
    """Check if file is a native file."""

    with open(file, 'rb') as f:
        if sys.platform == 'darwin':
            return f.read(4) in (
                b'\xfe\xed\xfa\xce',  # 32-bit
                b'\xfe\xed\xfa\xcf',  # 64-bit
                b'\xcf\xfa\xed\xfe',  # arm64
                b'\xca\xfe\xba\xbe',  # universal / fat (same as java class so beware!)
            )
        elif sys.platform == 'win32' or sys.platform == 'cygwin':
            return f.read(2) == b'MZ'
        else:
            # Assume that any other platform uses ELF binaries.
            return f.read(4) == b'\x7fELF'  # ELF


class _WheelBuilder():
    """Helper class to build wheels from projects."""

    def __init__(
        self,
        metadata: Metadata,
        manifest: Dict[str, List[Tuple[pathlib.Path, str]]],
        limited_api: bool,
    ) -> None:
        self._metadata = metadata
        self._manifest = manifest
        self._limited_api = limited_api

    @property
    def _has_internal_libs(self) -> bool:
        return bool(self._manifest.get('mesonpy-libs'))

    @property
    def _has_extension_modules(self) -> bool:
        # Assume that all code installed in {platlib} is Python ABI dependent.
        return bool(self._manifest.get('platlib'))

    @cached_property
    def _pure(self) -> bool:
        """Whether the wheel is architecture independent"""
        if self._manifest['platlib'] or self._manifest['mesonpy-libs']:
            return False
        for _, file in self._manifest['scripts']:
            if _is_native(file):
                return False
        return True

    @property
    def tag(self) -> mesonpy._tags.Tag:
        """Wheel tags."""
        if self._pure:
            return mesonpy._tags.Tag('py3', 'none', 'any')
        if not self._has_extension_modules:
            # The wheel has platform dependent code (is not pure) but
            # does not contain any extension module (does not
            # distribute any file in {platlib}) thus use generic
            # implementation and ABI tags.
            return mesonpy._tags.Tag('py3', 'none', None)
        return mesonpy._tags.Tag(None, self._stable_abi, None)

    @property
    def name(self) -> str:
        """Wheel name, this includes the basename and tag."""
        return f'{self._metadata.distribution_name}-{self._metadata.version}-{self.tag}'

    @property
    def _distinfo_dir(self) -> str:
        return f'{self._metadata.distribution_name}-{self._metadata.version}.dist-info'

    @property
    def _data_dir(self) -> str:
        return f'{self._metadata.distribution_name}-{self._metadata.version}.data'

    @property
    def _libs_dir(self) -> str:
        return f'.{self._metadata.distribution_name}.mesonpy.libs'

    @property
    def _license_file(self) -> Optional[pathlib.Path]:
        license_ = self._metadata.license
        if license_ and isinstance(license_, pyproject_metadata.License):
            return license_.file
        return None

    @property
    def wheel(self) -> bytes:
        """Return WHEEL file for dist-info."""
        return textwrap.dedent('''
            Wheel-Version: 1.0
            Generator: meson
            Root-Is-Purelib: {is_purelib}
            Tag: {tag}
        ''').strip().format(
            is_purelib='true' if self._pure else 'false',
            tag=self.tag,
        ).encode()

    @property
    def entrypoints_txt(self) -> bytes:
        """dist-info entry_points.txt."""
        data = self._metadata.entrypoints.copy()
        data.update({
            'console_scripts': self._metadata.scripts,
            'gui_scripts': self._metadata.gui_scripts,
        })

        text = ''
        for entrypoint in data:
            if data[entrypoint]:
                text += f'[{entrypoint}]\n'
                for name, target in data[entrypoint].items():
                    text += f'{name} = {target}\n'
                text += '\n'

        return text.encode()

    @cached_property
    def _stable_abi(self) -> Optional[str]:
        # PyPy supports the limited API but does not provide a stable
        # ABI, therefore extension modules using the limited API do
        # not use the stable ABI filename suffix and wheels should not
        # be tagged with the abi3 tag.
        if self._limited_api and '__pypy__' not in sys.builtin_module_names:
            # Verify stable ABI compatibility: examine files installed
            # in {platlib} that look like extension modules, and raise
            # an exception if any of them has a Python version
            # specific extension filename suffix ABI tag.
            for path, _ in self._manifest['platlib']:
                match = _EXTENSION_SUFFIX_REGEX.match(path.name)
                if match:
                    abi = match.group('abi')
                    if abi is not None and abi != 'abi3':
                        raise BuildError(
                            f'The package declares compatibility with Python limited API but extension '
                            f'module {os.fspath(path)!r} is tagged for a specific Python version.')
            return 'abi3'
        return None

    def _install_path(self, wheel_file: mesonpy._wheelfile.WheelFile, origin: Path, destination: pathlib.Path) -> None:
        """Add a file to the wheel."""

        if self._has_internal_libs:
            if _is_native(origin):
                # When an executable, libray, or Python extension module is
                # dynamically linked to a library built as part of the project,
                # Meson adds a library load path to it pointing to the build
                # directory, in the form of a relative RPATH entry. meson-python
                # relocates the shared libraries to the $project.mesonpy.libs
                # folder. Rewrite the RPATH to point to that folder instead.
                libspath = os.path.relpath(self._libs_dir, destination.parent)
                mesonpy._rpath.fix_rpath(origin, libspath)

        try:
            wheel_file.write(origin, destination.as_posix())
        except FileNotFoundError:
            # work around for Meson bug, see https://github.com/mesonbuild/meson/pull/11655
            if not os.fspath(origin).endswith('.pdb'):
                raise

    def _wheel_write_metadata(self, whl: mesonpy._wheelfile.WheelFile) -> None:
        # add metadata
        whl.writestr(f'{self._distinfo_dir}/METADATA', bytes(self._metadata.as_rfc822()))
        whl.writestr(f'{self._distinfo_dir}/WHEEL', self.wheel)
        if self.entrypoints_txt:
            whl.writestr(f'{self._distinfo_dir}/entry_points.txt', self.entrypoints_txt)

        # add license (see https://github.com/mesonbuild/meson-python/issues/88)
        if self._license_file:
            whl.write(self._license_file, f'{self._distinfo_dir}/{os.path.basename(self._license_file)}')

    def build(self, directory: Path) -> pathlib.Path:
        wheel_file = pathlib.Path(directory, f'{self.name}.whl')
        with mesonpy._wheelfile.WheelFile(wheel_file, 'w') as whl:
            self._wheel_write_metadata(whl)

            with _clicounter(sum(len(x) for x in self._manifest.values())) as counter:

                root = 'purelib' if self._pure else 'platlib'

                for path, entries in self._manifest.items():
                    for dst, src in entries:
                        counter.update(src)

                        if path == root:
                            pass
                        elif path == 'mesonpy-libs':
                            # custom installation path for bundled libraries
                            dst = pathlib.Path(self._libs_dir, dst)
                        else:
                            dst = pathlib.Path(self._data_dir, path, dst)

                        self._install_path(whl, src, dst)

        return wheel_file


class _EditableWheelBuilder(_WheelBuilder):

    @property
    def _top_level_modules(self) -> Collection[str]:
        modules = set()
        for type_ in self._manifest:
            for path, _ in self._manifest[type_]:
                name, dot, ext = path.parts[0].partition('.')
                if dot:
                    # module
                    suffix = dot + ext
                    if suffix in _SUFFIXES:
                        modules.add(name)
                else:
                    # package
                    modules.add(name)
        return modules

    def build(self, directory: Path, source_dir: pathlib.Path, build_dir: pathlib.Path,  # type: ignore[override]
              build_command: List[str], verbose: bool = False) -> pathlib.Path:

        wheel_file = pathlib.Path(directory, f'{self.name}.whl')
        with mesonpy._wheelfile.WheelFile(wheel_file, 'w') as whl:
            self._wheel_write_metadata(whl)
            whl.writestr(
                f'{self._distinfo_dir}/direct_url.json',
                source_dir.as_uri().encode('utf-8'))

            # install loader module
            loader_module_name = f'_{self._metadata.distribution_name}_editable_loader'
            whl.writestr(
                f'{loader_module_name}.py',
                read_binary('mesonpy', '_editable.py') + textwrap.dedent(f'''
                   install(
                       {self._metadata.name!r},
                       {self._top_level_modules!r},
                       {os.fspath(build_dir)!r},
                       {build_command!r},
                       {verbose!r},
                   )''').encode('utf-8'))

            # install .pth file
            whl.writestr(
                f'{self._metadata.canonical_name}-editable.pth',
                f'import {loader_module_name}'.encode('utf-8'))

        return wheel_file


def _validate_pyproject_config(pyproject: Dict[str, Any]) -> Dict[str, Any]:

    def _table(scheme: Dict[str, Callable[[Any, str], Any]]) -> Callable[[Any, str], Dict[str, Any]]:
        def func(value: Any, name: str) -> Dict[str, Any]:
            if not isinstance(value, dict):
                raise ConfigError(f'Configuration entry "{name}" must be a table')
            table = {}
            for key, val in value.items():
                check = scheme.get(key)
                if check is None:
                    raise ConfigError(f'Unknown configuration entry "{name}.{key}"')
                table[key] = check(val, f'{name}.{key}')
            return table
        return func

    def _strings(value: Any, name: str) -> List[str]:
        if not isinstance(value, list) or not all(isinstance(x, str) for x in value):
            raise ConfigError(f'Configuration entry "{name}" must be a list of strings')
        return value

    def _bool(value: Any, name: str) -> bool:
        if not isinstance(value, bool):
            raise ConfigError(f'Configuration entry "{name}" must be a boolean')
        return value

    def _string_or_path(value: Any, name: str) -> str:
        if not isinstance(value, str):
            raise ConfigError(f'Configuration entry "{name}" must be a string')
        if os.path.isfile(value):
            value = os.path.abspath(value)
        return value

    scheme = _table({
        'meson': _string_or_path,
        'limited-api': _bool,
        'args': _table({
            name: _strings for name in _MESON_ARGS_KEYS
        }),
    })

    table = pyproject.get('tool', {}).get('meson-python', {})
    return scheme(table, 'tool.meson-python')


def _validate_config_settings(config_settings: Dict[str, Any]) -> Dict[str, Any]:
    """Validate options received from build frontend."""

    def _string(value: Any, name: str) -> str:
        if not isinstance(value, str):
            raise ConfigError(f'Only one value for "{name}" can be specified')
        return value

    def _bool(value: Any, name: str) -> bool:
        return True

    def _string_or_strings(value: Any, name: str) -> List[str]:
        return list([value,] if isinstance(value, str) else value)

    options = {
        'builddir': _string,
        'build-dir': _string,
        'editable-verbose': _bool,
        'dist-args': _string_or_strings,
        'setup-args': _string_or_strings,
        'compile-args': _string_or_strings,
        'install-args': _string_or_strings,
    }
    assert all(f'{name}-args' in options for name in _MESON_ARGS_KEYS)

    config = {}
    for key, value in config_settings.items():
        parser = options.get(key)
        if parser is None:
            matches = difflib.get_close_matches(key, options.keys(), n=2)
            if matches:
                alternatives = ' or '.join(f'"{match}"' for match in matches)
                raise ConfigError(f'Unknown option "{key}". Did you mean {alternatives}?')
            else:
                raise ConfigError(f'Unknown option "{key}"')
        config[key] = parser(value, key)

    # Check backward compatibility aliases.
    aliases = {
        'build-dir': 'builddir',
    }
    for key, alt in aliases.items():
        if key in config and alt in config:
            raise ConfigError(f'Option "{alt}" is a backward compatibility alias for "{key}". Only one can be used')
        if alt in config:
            config[key] = config[alt]

    return config


class Project():
    """Meson project wrapper to generate Python artifacts."""

    def __init__(
        self,
        source_dir: Path,
        build_dir: Path,
        meson_args: Optional[MesonArgs] = None,
        editable_verbose: bool = False,
    ) -> None:
        self._source_dir = pathlib.Path(source_dir).absolute()
        self._build_dir = pathlib.Path(build_dir).absolute()
        self._editable_verbose = editable_verbose
        self._meson_native_file = self._build_dir / 'meson-python-native-file.ini'
        self._meson_cross_file = self._build_dir / 'meson-python-cross-file.ini'
        self._meson_args: MesonArgs = collections.defaultdict(list)
        self._limited_api = False

        # load pyproject.toml
        pyproject = tomllib.loads(self._source_dir.joinpath('pyproject.toml').read_text(encoding='utf-8'))

        # load meson args from pyproject.toml
        pyproject_config = _validate_pyproject_config(pyproject)
        for key, value in pyproject_config.get('args', {}).items():
            self._meson_args[key].extend(value)

        # meson arguments from the command line take precedence over
        # arguments from the configuration file thus are added later
        if meson_args:
            for key, value in meson_args.items():
                self._meson_args[key].extend(value)

        # determine command to invoke meson
        self._meson = _get_meson_command(pyproject_config.get('meson'))

        self._ninja = _env_ninja_command()
        if self._ninja is None:
            raise ConfigError(f'Could not find ninja version {_NINJA_REQUIRED_VERSION} or newer.')
        os.environ.setdefault('NINJA', self._ninja)

        # make sure the build dir exists
        self._build_dir.mkdir(exist_ok=True, parents=True)

        # if the build dir is empty, add .gitignore and .hgignore files
        if not any(self._build_dir.iterdir()):
            _add_ignore_files(self._build_dir)

        # setuptools-like ARCHFLAGS environment variable support
        if sysconfig.get_platform().startswith('macosx-'):
            archflags = os.environ.get('ARCHFLAGS', '').strip()
            if archflags:

                # parse the ARCHFLAGS environment variable
                parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
                parser.add_argument('-arch', action='append')
                args, unknown = parser.parse_known_args(archflags.split())
                if unknown:
                    raise ConfigError(f'Unknown flag specified in $ARCHFLAGS={archflags!r}')
                arch, *other = set(args.arch)
                if other:
                    raise ConfigError(f'Multi-architecture builds are not supported but $ARCHFLAGS={archflags!r}')

                macver, _, nativearch = platform.mac_ver()
                if arch != nativearch:
                    x = os.environ.setdefault('_PYTHON_HOST_PLATFORM', f'macosx-{macver}-{arch}')
                    if not x.endswith(arch):
                        raise ConfigError(f'$ARCHFLAGS={archflags!r} and $_PYTHON_HOST_PLATFORM={x!r} do not agree')
                    family = 'aarch64' if arch == 'arm64' else arch
                    cross_file_data = textwrap.dedent(f'''
                        [binaries]
                        c = ['cc', '-arch', {arch!r}]
                        cpp = ['c++', '-arch', {arch!r}]
                        objc = ['cc', '-arch', {arch!r}]
                        objcpp = ['c++', '-arch', {arch!r}]
                        [host_machine]
                        system = 'darwin'
                        cpu = {arch!r}
                        cpu_family = {family!r}
                        endian = 'little'
                    ''')
                    self._meson_cross_file.write_text(cross_file_data, encoding='utf-8')
                    self._meson_args['setup'].extend(('--cross-file', os.fspath(self._meson_cross_file)))

        # write the native file
        native_file_data = textwrap.dedent(f'''
            [binaries]
            python = '{sys.executable}'
        ''')
        self._meson_native_file.write_text(native_file_data, encoding='utf-8')

        # reconfigure if we have a valid Meson build directory. Meson
        # uses the presence of the 'meson-private/coredata.dat' file
        # in the build directory as indication that the build
        # directory has already been configured and arranges this file
        # to be created as late as possible or deleted if something
        # goes wrong during setup.
        reconfigure = self._build_dir.joinpath('meson-private/coredata.dat').is_file()

        # run meson setup
        self._configure(reconfigure=reconfigure)

        # package metadata
        if 'project' in pyproject:
            self._metadata = Metadata.from_pyproject(pyproject, self._source_dir)
            # set version from meson.build if version is declared as dynamic
            if 'version' in self._metadata.dynamic:
                version = self._meson_version
                if version == 'undefined':
                    raise pyproject_metadata.ConfigurationError(
                        'Field "version" declared as dynamic but version is not defined in meson.build')
                self._metadata.version = packaging.version.Version(version)
        else:
            # if project section is missing, use minimal metdata from meson.build
            name, version = self._meson_name, self._meson_version
            if version == 'undefined':
                raise pyproject_metadata.ConfigurationError(
                    'Section "project" missing in pyproject.toml and version is not defined in meson.build')
            self._metadata = Metadata(name=name, version=packaging.version.Version(version))

        # verify that we are running on a supported interpreter
        if self._metadata.requires_python:
            self._metadata.requires_python.prereleases = True
            if platform.python_version().rstrip('+') not in self._metadata.requires_python:
                raise MesonBuilderError(
                    f'The package requires Python version {self._metadata.requires_python}, '
                    f'running on {platform.python_version()}')

        # limited API
        self._limited_api = pyproject_config.get('limited-api', False)
        if self._limited_api:
            # check whether limited API is disabled for the Meson project
            options = self._info('intro-buildoptions')
            value = next((option['value'] for option in options if option['name'] == 'python.allow_limited_api'), None)
            if not value:
                self._limited_api = False

        if self._limited_api and bool(sysconfig.get_config_var('Py_GIL_DISABLED')):
            raise BuildError(
                'The package targets Python\'s Limited API, which is not supported by free-threaded CPython. '
                'The "python.allow_limited_api" Meson build option may be used to override the package default.')

    def _run(self, cmd: Sequence[str]) -> None:
        """Invoke a subprocess."""
        # Flush the line to ensure that the log line with the executed
        # command line appears before the command output. Without it,
        # the lines appear in the wrong order in pip output.
        _log('{style.INFO}+ {cmd}{style.RESET}'.format(style=style, cmd=' '.join(cmd)), flush=True)
        r = subprocess.run(cmd, cwd=self._build_dir)
        if r.returncode != 0:
            raise SystemExit(r.returncode)

    def _configure(self, reconfigure: bool = False) -> None:
        """Configure Meson project."""
        setup_args = [
            os.fspath(self._source_dir),
            os.fspath(self._build_dir),
            # default build options
            '-Dbuildtype=release',
            '-Db_ndebug=if-release',
            '-Db_vscrt=md',
            # user build options
            *self._meson_args['setup'],
            # pass native file last to have it override the python
            # interpreter path that may have been specified in user
            # provided native files
            f'--native-file={os.fspath(self._meson_native_file)}',
        ]
        if reconfigure:
            setup_args.insert(0, '--reconfigure')
        self._run(self._meson + ['setup', *setup_args])

    @property
    def _build_command(self) -> List[str]:
        assert self._ninja is not None  # help mypy out
        if sys.platform == 'win32':
            # On Windows use 'meson compile' to setup the MSVC compiler
            # environment. Using the --ninja-args option allows to
            # provide the exact same semantics for the compile arguments
            # provided by the users.
            cmd = self._meson + ['compile']
            args = list(self._meson_args['compile'])
            if args:
                cmd.append(f'--ninja-args={args!r}')
            return cmd
        return [self._ninja, *self._meson_args['compile']]

    @functools.lru_cache(maxsize=None)
    def build(self) -> None:
        """Build the Meson project."""
        self._run(self._build_command)

    @functools.lru_cache()
    def _info(self, name: str) -> Any:
        """Read info from meson-info directory."""
        info = self._build_dir.joinpath('meson-info', f'{name}.json')
        return json.loads(info.read_text(encoding='utf-8'))

    @property
    def _manifest(self) -> DefaultDict[str, List[Tuple[pathlib.Path, str]]]:
        """The files to be added to the wheel, organized by wheel path."""

        # Obtain the list of files Meson would install.
        install_plan = self._info('intro-install_plan')

        # Parse the 'meson install' args to extract --tags and --skip-subprojects
        parser = argparse.ArgumentParser(add_help=False)
        parser.add_argument('--tags')
        parser.add_argument('--skip-subprojects', nargs='?', const='*', default='')
        args, _ = parser.parse_known_args(self._meson_args['install'])
        install_tags = {t.strip() for t in args.tags.split(',')} if args.tags else None
        skip_subprojects = {p for p in (p.strip() for p in args.skip_subprojects.split(',')) if p}

        # Filter the install plan accordingly.
        sources: DefaultDict[str, Dict[str, Dict[str, str]]] = collections.defaultdict(dict)
        for key, targets in install_plan.items():
            for target, details in targets.items():
                if install_tags is not None and details['tag'] not in install_tags:
                    continue
                subproject = details.get('subproject')
                if subproject is not None and (subproject in skip_subprojects or '*' in skip_subprojects):
                    continue
                sources[key][target] = details

        # Map Meson installation locations to wheel paths.
        return _map_to_wheel(sources)

    @property
    def _meson_name(self) -> str:
        """Name in meson.build."""
        name = self._info('intro-projectinfo')['descriptive_name']
        assert isinstance(name, str)
        return name

    @property
    def _meson_version(self) -> str:
        """Version in meson.build."""
        name = self._info('intro-projectinfo')['version']
        assert isinstance(name, str)
        return name

    def sdist(self, directory: Path) -> pathlib.Path:
        """Generates a sdist (source distribution) in the specified directory."""
        # Generate meson dist file.
        self._run(self._meson + ['dist', '--allow-dirty', '--no-tests', '--formats', 'gztar', *self._meson_args['dist']])

        dist_name = f'{self._metadata.distribution_name}-{self._metadata.version}'
        meson_dist_name = f'{self._meson_name}-{self._meson_version}'
        meson_dist_path = pathlib.Path(self._build_dir, 'meson-dist', f'{meson_dist_name}.tar.gz')
        sdist_path = pathlib.Path(directory, f'{dist_name}.tar.gz')
        pyproject_toml_mtime = 0

        with tarfile.open(meson_dist_path, 'r:gz') as meson_dist, mesonpy._util.create_targz(sdist_path) as sdist:
            for member in meson_dist.getmembers():
                if member.isfile():
                    file = meson_dist.extractfile(member.name)

                    # Reset pax extended header.  The tar archive member may be
                    # using pax headers to store some file metadata.  The pax
                    # headers are not reset when the metadata is modified and
                    # they take precedence when the member is deserialized.
                    # This is relevant because when rewriting the member name,
                    # the length of the path may shrink from being more than
                    # 100 characters (requiring the path to be stored in the
                    # pax headers) to being less than 100 characters. When this
                    # happens, the tar archive member is serialized with the
                    # shorter name in the regular header and the longer one in
                    # the extended pax header.  The archives handled here are
                    # not expected to use extended pax headers other than for
                    # the ones required to encode file metadata.  The easiest
                    # solution is to reset the pax extended headers.
                    member.pax_headers = {}

                    # Rewrite the path to match the sdist distribution name.
                    stem = member.name.split('/', 1)[1]
                    member.name = '/'.join((dist_name, stem))

                    if stem == 'pyproject.toml':
                        pyproject_toml_mtime = member.mtime

                    # Reset owner and group to root:root.  This mimics what
                    # 'git archive' does and makes the sdist reproducible upon
                    # being built by different users.
                    member.uname = member.gname = 'root'
                    member.uid = member.gid = 0

                    sdist.addfile(member, file)

            # Add 'PKG-INFO'.
            member = tarfile.TarInfo(f'{dist_name}/PKG-INFO')
            member.uid = member.gid = 0
            member.uname = member.gname = 'root'

            # Set the 'PKG-INFO' modification time to the modification time of
            # 'pyproject.toml' in the archive generated by 'meson dist'.  In
            # turn this is the last commit time, unless touched by a dist
            # script.  This makes the sdist reproducible upon being built at
            # different times, when dist scripts are not used, which should be
            # the majority of cases.
            #
            # Note that support for dynamic version in project metadata allows
            # the version to depend on the build time.  Therefore, setting the
            # 'PKG-INFO' modification time to the 'pyproject.toml'
            # modification time can be seen as not strictly correct.  However,
            # the sdist standard does not dictate which modification time to
            # use for 'PKG-INFO'.  This choice allows to make the sdist
            # byte-for-byte reproducible in the most common case.
            member.mtime = pyproject_toml_mtime

            metadata = bytes(self._metadata.as_rfc822())
            member.size = len(metadata)
            sdist.addfile(member, io.BytesIO(metadata))

        return sdist_path

    def wheel(self, directory: Path) -> pathlib.Path:
        """Generates a wheel in the specified directory."""
        self.build()
        builder = _WheelBuilder(self._metadata, self._manifest, self._limited_api)
        return builder.build(directory)

    def editable(self, directory: Path) -> pathlib.Path:
        """Generates an editable wheel in the specified directory."""
        self.build()
        builder = _EditableWheelBuilder(self._metadata, self._manifest, self._limited_api)
        return builder.build(directory, self._source_dir, self._build_dir, self._build_command, self._editable_verbose)


@contextlib.contextmanager
def _project(config_settings: Optional[Dict[Any, Any]] = None) -> Iterator[Project]:
    """Create the project given the given config settings."""

    settings = _validate_config_settings(config_settings or {})
    meson_args = typing.cast('MesonArgs', {name: settings.get(f'{name}-args', []) for name in _MESON_ARGS_KEYS})
    source_dir = os.path.curdir
    build_dir = settings.get('build-dir')
    editable_verbose = bool(settings.get('editable-verbose'))

    with contextlib.ExitStack() as ctx:
        if build_dir is None:
            build_dir = ctx.enter_context(tempfile.TemporaryDirectory(prefix='.mesonpy-', dir=source_dir))
        yield Project(source_dir, build_dir, meson_args, editable_verbose)


def _parse_version_string(string: str) -> Tuple[int, ...]:
    """Parse version string."""
    try:
        return tuple(map(int, string.split('.')[:3]))
    except ValueError:
        return (0, )


def _get_meson_command(
        meson: Optional[str] = None, *, version: str = _MESON_REQUIRED_VERSION
    ) -> List[str]:
    """Return the command to invoke meson."""

    # The MESON env var, if set, overrides the config value from pyproject.toml.
    # The config value, if given, is an absolute path or the name of an executable.
    meson = os.environ.get('MESON', meson or 'meson')

    # If the specified Meson string ends in `.py`, we run it with the current
    # Python executable. This avoids problems for users on Windows, where
    # making a script executable isn't enough to get it to run when invoked
    # directly. For packages that vendor a forked Meson, the `meson.py` in the
    # root of the Meson repo can be used this way.
    if meson.endswith('.py'):
        if not os.path.exists(meson):
            raise ConfigError(f'Could not find the specified meson: "{meson}"')
        cmd = [sys.executable, meson]
    else:
        cmd = [meson]

    # The meson Python package is a dependency of the meson-python Python
    # package, however, it may occur that the meson Python package is installed
    # but the corresponding meson command is not available in $PATH. Implement
    # a runtime check to verify that the build environment is setup correcly.
    try:
        r = subprocess.run(cmd + ['--version'], text=True, capture_output=True)
    except FileNotFoundError as err:
        raise ConfigError(f'meson executable "{meson}" not found') from err
    if r.returncode != 0:
        raise ConfigError(f'Could not execute meson: {r.stderr.strip()}')
    meson_version = r.stdout.strip()

    if _parse_version_string(meson_version) < _parse_version_string(version):
        raise ConfigError(f'Could not find meson version {version} or newer, found {meson_version}.')

    return cmd


def _env_ninja_command(*, version: str = _NINJA_REQUIRED_VERSION) -> Optional[str]:
    """Returns the path to ninja, or None if no ninja found."""
    required_version = _parse_version_string(version)
    env_ninja = os.environ.get('NINJA')
    ninja_candidates = [env_ninja] if env_ninja else ['ninja', 'ninja-build', 'samu']
    for ninja in ninja_candidates:
        ninja_path = shutil.which(ninja)
        if ninja_path is not None:
            version = subprocess.run([ninja_path, '--version'], check=False, text=True, capture_output=True).stdout
            if _parse_version_string(version) >= required_version:
                return ninja_path
    return None


def _add_ignore_files(directory: pathlib.Path) -> None:
    directory.joinpath('.gitignore').write_text(textwrap.dedent('''
        # This file is generated by meson-python. It will not be recreated if deleted or modified.
        *
    '''), encoding='utf-8')
    directory.joinpath('.hgignore').write_text(textwrap.dedent('''
        # This file is generated by meson-python. It will not be recreated if deleted or modified.
        syntax: glob
        **/*
    '''), encoding='utf-8')


def _pyproject_hook(func: Callable[P, T]) -> Callable[P, T]:
    @functools.wraps(func)
    def wrapper(*args: P.args, **kwargs: P.kwargs) -> T:
        warnings.showwarning = _showwarning
        try:
            return func(*args, **kwargs)
        except (Error, pyproject_metadata.ConfigurationError) as exc:
            prefix = f'{style.ERROR}meson-python: error:{style.RESET} '
            _log('\n' + textwrap.indent(str(exc), prefix))
            raise SystemExit(1) from exc
    return wrapper


@_pyproject_hook
def get_requires_for_build_sdist(config_settings: Optional[Dict[str, str]] = None) -> List[str]:
    dependencies = []

    if os.environ.get('NINJA') is None and _env_ninja_command() is None:
        dependencies.append(f'ninja >= {_NINJA_REQUIRED_VERSION}')

    return dependencies


@_pyproject_hook
def get_requires_for_build_wheel(config_settings: Optional[Dict[str, str]] = None) -> List[str]:
    dependencies = []

    if os.environ.get('NINJA') is None and _env_ninja_command() is None:
        dependencies.append(f'ninja >= {_NINJA_REQUIRED_VERSION}')

    if sys.platform.startswith('linux') and not shutil.which('patchelf'):
        dependencies.append('patchelf >= 0.11.0')

    return dependencies


get_requires_for_build_editable = get_requires_for_build_wheel


@_pyproject_hook
def build_sdist(
    sdist_directory: str,
    config_settings: Optional[Dict[Any, Any]] = None,
) -> str:

    out = pathlib.Path(sdist_directory)
    with _project(config_settings) as project:
        return project.sdist(out).name


@_pyproject_hook
def build_wheel(
    wheel_directory: str, config_settings:
    Optional[Dict[Any, Any]] = None,
    metadata_directory: Optional[str] = None,
) -> str:

    out = pathlib.Path(wheel_directory)
    with _project(config_settings) as project:
        return project.wheel(out).name


@_pyproject_hook
def build_editable(
    wheel_directory: str,
    config_settings: Optional[Dict[Any, Any]] = None,
    metadata_directory: Optional[str] = None,
) -> str:

    # Force set a permanent build directory.
    if not config_settings:
        config_settings = {}
    if 'build-dir' not in config_settings and 'builddir' not in config_settings:
        config_settings['build-dir'] = 'build/' + mesonpy._tags.get_abi_tag()

    out = pathlib.Path(wheel_directory)
    with _project(config_settings) as project:
        return project.editable(out).name
