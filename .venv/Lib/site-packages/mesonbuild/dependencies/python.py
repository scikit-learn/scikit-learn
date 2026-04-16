# SPDX-License-Identifier: Apache-2.
# Copyright 2022 The Meson development team

from __future__ import annotations

import functools, json, operator, os, textwrap
from pathlib import Path
import typing as T

from .. import mesonlib, mlog
from .base import process_method_kw, DependencyCandidate, DependencyException, DependencyMethods, ExternalDependency, SystemDependency
from .configtool import ConfigToolDependency
from .detect import packages
from .factory import DependencyFactory
from .framework import ExtraFrameworkDependency
from .pkgconfig import PkgConfigDependency
from ..envconfig import detect_cpu_family
from ..mesonlib import MachineChoice, path_is_in_root
from ..programs import ExternalProgram
from ..options import OptionKey
from ..scripts import destdir_join

if T.TYPE_CHECKING:
    from typing_extensions import Final, TypedDict

    from .factory import DependencyGenerator
    from ..environment import Environment
    from .base import DependencyObjectKWs

    class PythonIntrospectionDict(TypedDict):

        install_paths: T.Dict[str, str]
        is_pypy: bool
        is_venv: bool
        is_freethreaded: bool
        link_libpython: bool
        sysconfig_paths: T.Dict[str, str]
        paths: T.Dict[str, str]
        platform: str
        suffix: str
        limited_api_suffix: str
        variables: T.Dict[str, str]
        version: str

    _Base = ExternalDependency
else:
    _Base = object


class Pybind11ConfigToolDependency(ConfigToolDependency):

    tools = ['pybind11-config']

    # any version of the tool is valid, since this is header-only
    allow_default_for_cross = True

    # pybind11 in 2.10.4 added --version, sanity-check another flag unique to it
    # in the meantime
    skip_version = '--pkgconfigdir'

    def __init__(self, name: str, environment: Environment, kwargs: DependencyObjectKWs):
        super().__init__(name, environment, kwargs)
        if not self.is_found:
            return
        self.compile_args = self.get_config_value(['--includes'], 'compile_args')


class NumPyConfigToolDependency(ConfigToolDependency):

    tools = ['numpy-config']

    def __init__(self, name: str, environment: Environment, kwargs: DependencyObjectKWs):
        super().__init__(name, environment, kwargs)
        if not self.is_found:
            return
        self.compile_args = self.get_config_value(['--cflags'], 'compile_args')


class PythonBuildConfig:
    """PEP 739 build-details.json config file."""

    IMPLEMENTED_VERSION: Final[str] = '1.0'
    """Schema version currently implemented."""
    _PATH_KEYS = (
        'base_interpreter',
        'libpython.dynamic',
        'libpython.dynamic_stableabi',
        'libpython.static',
        'c_api.headers',
        'c_api.pkgconfig_path',
    )
    """Path keys — may be relative, need to be expanded."""

    def __init__(self, path: str) -> None:
        self._path = Path(path)

        try:
            self._data = json.loads(self._path.read_text(encoding='utf8'))
        except OSError as e:
            raise DependencyException(f'Failed to read python.build_config: {e}') from e

        self._validate_data()
        self._expand_paths()

    def __getitem__(self, key: str) -> T.Any:
        return functools.reduce(operator.getitem, key.split('.'), self._data)

    def __contains__(self, key: str) -> bool:
        try:
            self[key]
        except KeyError:
            return False
        else:
            return True

    def get(self, key: str, default: T.Any = None) -> T.Any:
        try:
            return self[key]
        except KeyError:
            return default

    def _validate_data(self) -> None:
        schema_version = self._data['schema_version']
        if mesonlib.version_compare(schema_version, '< 1.0'):
            raise DependencyException(f'Invalid schema_version in python.build_config: {schema_version}')
        if mesonlib.version_compare(schema_version, '>= 2.0'):
            raise DependencyException(
                f'Unsupported schema_version {schema_version!r} in python.build_config, '
                f'but we only implement support for {self.IMPLEMENTED_VERSION!r}'
            )
        # Schema version that we currently understand
        if mesonlib.version_compare(schema_version, f'> {self.IMPLEMENTED_VERSION}'):
            mlog.log(
                f'python.build_config has schema_version {schema_version!r}, '
                f'but we only implement support for {self.IMPLEMENTED_VERSION!r}, '
                'new functionality might be missing'
            )

    def _expand_paths(self) -> None:
        """Expand relative path (they're relative to base_prefix)."""
        for key in self._PATH_KEYS:
            if key not in self:
                continue
            parent, _, child = key.rpartition('.')
            container = self[parent] if parent else self._data
            path = Path(container[child])
            if not path.is_absolute():
                container[child] = os.fspath(self.base_prefix / path)

    @property
    def config_path(self) -> Path:
        return self._path

    @mesonlib.lazy_property
    def base_prefix(self) -> Path:
        path = Path(self._data['base_prefix'])
        if path.is_absolute():
            return path
        # Non-absolute paths are relative to the build config directory
        return self.config_path.parent / path


class BasicPythonExternalProgram(ExternalProgram):
    def __init__(self, name: str, command: T.Optional[T.List[str]] = None,
                 ext_prog: T.Optional[ExternalProgram] = None,
                 build_config_path: T.Optional[str] = None):
        if ext_prog is None:
            super().__init__(name, command=command, silent=True)
        else:
            self.name = name
            self.command = ext_prog.command
            self.path = ext_prog.path
            self.cached_version = None
            self.version_arg = '--version'

        self.build_config = PythonBuildConfig(build_config_path) if build_config_path else None

        # We want strong key values, so we always populate this with bogus data.
        # Otherwise to make the type checkers happy we'd have to do .get() for
        # everycall, even though we know that the introspection data will be
        # complete
        self.info: 'PythonIntrospectionDict' = {
            'install_paths': {},
            'is_pypy': False,
            'is_venv': False,
            'is_freethreaded': False,
            'link_libpython': False,
            'sysconfig_paths': {},
            'paths': {},
            'platform': 'sentinel',
            'suffix': 'sentinel',
            'limited_api_suffix': 'sentinel',
            'variables': {},
            'version': '0.0',
        }
        self.pure: bool = True

    @property
    def version(self) -> str:
        if self.build_config:
            value = self.build_config['language']['version']
        else:
            value = self.info['variables'].get('LDVERSION') or self.info['version']
        assert isinstance(value, str)
        return value

    def _check_version(self, version: str) -> bool:
        if self.name == 'python2':
            return mesonlib.version_compare(version, '< 3.0')
        elif self.name == 'python3':
            return mesonlib.version_compare(version, '>= 3.0')
        return True

    def sanity(self) -> bool:
        # Sanity check, we expect to have something that at least quacks in tune

        if self.build_config:
            if not self.build_config['libpython']:
                mlog.debug('This Python installation does not provide a libpython')
                return False
            if not self.build_config['c_api']:
                mlog.debug('This Python installation does support the C API')
                return False

        import importlib.resources

        with importlib.resources.path('mesonbuild.scripts', 'python_info.py') as f:
            cmd = self.get_command() + [str(f)]
            env = os.environ.copy()
            env['SETUPTOOLS_USE_DISTUTILS'] = 'stdlib'
            p, stdout, stderr = mesonlib.Popen_safe(cmd, env=env)

        try:
            info = json.loads(stdout)
        except json.JSONDecodeError:
            info = None
            mlog.debug('Could not introspect Python (%s): exit code %d' % (str(p.args), p.returncode))
            mlog.debug('Program stdout:\n')
            mlog.debug(stdout)
            mlog.debug('Program stderr:\n')
            mlog.debug(stderr)

        if info is not None and self._check_version(info['version']):
            self.info = T.cast('PythonIntrospectionDict', info)
            return True
        else:
            return False


class _PythonDependencyBase(_Base):

    for_machine: MachineChoice

    def __init__(self, python_holder: 'BasicPythonExternalProgram', embed: bool):
        self.embed = embed
        self.build_config = python_holder.build_config

        if self.build_config:
            self.version = self.build_config['language']['version']
            self.platform = self.build_config['platform']
            self.is_freethreaded = 't' in self.build_config['abi']['flags']
            self.link_libpython = self.build_config['libpython']['link_extensions']
            # TODO: figure out how to deal with frameworks
            # see the logic at the bottom of PythonPkgConfigDependency.__init__()
            if self.env.machines.host.is_darwin():
                raise DependencyException('--python.build-config is not supported on Darwin')
        else:
            self.version = python_holder.info['version']
            self.platform = python_holder.info['platform']
            self.is_freethreaded = python_holder.info['is_freethreaded']
            self.link_libpython = python_holder.info['link_libpython']
            # This data shouldn't be needed when build_config is set
            self.is_pypy = python_holder.info['is_pypy']
            self.variables = python_holder.info['variables']

        self.paths = python_holder.info['paths']

        # The "-embed" version of python.pc / python-config was introduced in 3.8,
        # and distutils extension linking was changed to be considered a non embed
        # usage. Before then, this dependency always uses the embed=True handling
        # because that is the only one that exists.
        #
        # On macOS and some Linux distros (Debian) distutils doesn't link extensions
        # against libpython, even on 3.7 and below. We call into distutils and
        # mirror its behavior. See https://github.com/mesonbuild/meson/issues/4117
        if not self.link_libpython:
            self.link_libpython = embed

        self.info: T.Optional[T.Dict[str, str]] = None
        if mesonlib.version_compare(self.version, '>= 3.0'):
            self.major_version = 3
        else:
            self.major_version = 2

        # pyconfig.h is shared between regular and free-threaded builds in the
        # Windows installer from python.org, and hence does not define
        # Py_GIL_DISABLED correctly. So do it here:
        if mesonlib.is_windows() and self.is_freethreaded:
            self.compile_args += ['-DPy_GIL_DISABLED']

    def find_libpy(self, environment: 'Environment') -> None:
        if self.build_config:
            path = self.build_config['libpython'].get('dynamic')
            if not path:
                raise DependencyException('Python does not provide a dynamic libpython library')
            sysroot = environment.properties[self.for_machine].get_sys_root()
            if sysroot and not path_is_in_root(Path(path), Path(sysroot)):
                path = destdir_join(sysroot, path)
            if not os.path.isfile(path):
                raise DependencyException('Python dynamic library does not exist or is not a file')
            self.link_args = [path]
            self.is_found = True
            return

        if self.is_pypy:
            if self.major_version == 3:
                libname = 'pypy3-c'
            else:
                libname = 'pypy-c'
            libdir = os.path.join(self.variables.get('base'), 'bin')
            libdirs = [libdir]
        else:
            libname = f'python{self.version}'
            if 'DEBUG_EXT' in self.variables:
                libname += self.variables['DEBUG_EXT']
            if 'ABIFLAGS' in self.variables:
                libname += self.variables['ABIFLAGS']
            libdirs = []

        largs = self.clib_compiler.find_library(libname, libdirs)
        if largs is not None:
            self.link_args = largs
            self.is_found = True

    def get_windows_python_arch(self) -> str:
        if self.platform.startswith('mingw'):
            if 'x86_64' in self.platform:
                return 'x86_64'
            elif 'i686' in self.platform:
                return 'x86'
            elif 'aarch64' in self.platform:
                return 'aarch64'
            else:
                raise DependencyException(f'MinGW Python built with unknown platform {self.platform!r}, please file a bug')
        elif self.platform == 'win32':
            return 'x86'
        elif self.platform in {'win64', 'win-amd64'}:
            return 'x86_64'
        elif self.platform in {'win-arm64'}:
            return 'aarch64'
        raise DependencyException('Unknown Windows Python platform {self.platform!r}')

    def get_windows_link_args(self, limited_api: bool, environment: 'Environment') -> T.Optional[T.List[str]]:
        if self.build_config:
            if self.static:
                key = 'static'
            elif limited_api:
                key = 'dynamic-stableabi'
            else:
                key = 'dynamic'
            sysroot = environment.properties[self.for_machine].get_sys_root()
            path = self.build_config['libpython'][key]
            if sysroot and not path_is_in_root(Path(path), Path(sysroot)):
                path = destdir_join(sysroot, path)
            return [path]

        if self.platform.startswith('win'):
            vernum = self.variables.get('py_version_nodot')
            verdot = self.variables.get('py_version_short')
            imp_lower = self.variables.get('implementation_lower', 'python')
            if self.static:
                libpath = Path('libs') / f'libpython{vernum}.a'
            else:
                if limited_api:
                    vernum = vernum[0]
                comp = self.get_compiler()
                if comp.id == "gcc":
                    if imp_lower == 'pypy' and verdot == '3.8':
                        # The naming changed between 3.8 and 3.9
                        libpath = Path('libpypy3-c.dll')
                    elif imp_lower == 'pypy':
                        libpath = Path(f'libpypy{verdot}-c.dll')
                    else:
                        if self.is_freethreaded:
                            libpath = Path(f'python{vernum}t.dll')
                        else:
                            libpath = Path(f'python{vernum}.dll')
                else:
                    if self.is_freethreaded:
                        libpath = Path('libs') / f'python{vernum}t.lib'
                    else:
                        libpath = Path('libs') / f'python{vernum}.lib'
                    # For a debug build, pyconfig.h may force linking with
                    # pythonX_d.lib (see meson#10776). This cannot be avoided
                    # and won't work unless we also have a debug build of
                    # Python itself (except with pybind11, which has an ugly
                    # hack to work around this) - so emit a warning to explain
                    # the cause of the expected link error.
                    buildtype = self.env.coredata.optstore.get_value_for(OptionKey('buildtype'))
                    assert isinstance(buildtype, str)
                    debug = self.env.coredata.optstore.get_value_for(OptionKey('debug'))
                    # `debugoptimized` buildtype may not set debug=True currently, see gh-11645
                    is_debug_build = debug or buildtype == 'debug'
                    vscrt_debug = False
                    if OptionKey('b_vscrt') in self.env.coredata.optstore:
                        vscrt = self.env.coredata.optstore.get_value_for('b_vscrt')
                        if vscrt in {'mdd', 'mtd', 'from_buildtype', 'static_from_buildtype'}:
                            vscrt_debug = True
                    if is_debug_build and vscrt_debug and not self.variables.get('Py_DEBUG'):
                        mlog.warning(textwrap.dedent('''\
                            Using a debug build type with MSVC or an MSVC-compatible compiler
                            when the Python interpreter is not also a debug build will almost
                            certainly result in a failed build. Prefer using a release build
                            type or a debug Python interpreter.
                            '''))
            # base_prefix to allow for virtualenvs.
            lib = Path(self.variables.get('base_prefix')) / libpath
        elif self.platform.startswith('mingw'):
            if self.static:
                if limited_api:
                    libname = self.variables.get('ABI3DLLLIBRARY')
                else:
                    libname = self.variables.get('LIBRARY')
            else:
                if limited_api:
                    libname = self.variables.get('ABI3LDLIBRARY')
                else:
                    libname = self.variables.get('LDLIBRARY')
            lib = Path(self.variables.get('LIBDIR')) / libname
        else:
            raise mesonlib.MesonBugException(
                'On a Windows path, but the OS doesn\'t appear to be Windows or MinGW.')
        if not lib.exists():
            mlog.log('Could not find Python3 library {!r}'.format(str(lib)))
            return None
        return [str(lib)]

    def find_libpy_windows(self, env: 'Environment', limited_api: bool = False) -> None:
        '''
        Find python3 libraries on Windows and also verify that the arch matches
        what we are building for.
        '''
        try:
            pyarch = self.get_windows_python_arch()
        except DependencyException as e:
            mlog.log(str(e))
            self.is_found = False
            return
        arch = detect_cpu_family(env.coredata.compilers.host)
        if arch != pyarch:
            mlog.log('Need', mlog.bold(self.name), f'for {arch}, but found {pyarch}')
            self.is_found = False
            return
        # This can fail if the library is not found
        largs = self.get_windows_link_args(limited_api, env)
        if largs is None:
            self.is_found = False
            return
        self.link_args = largs
        self.is_found = True


class PythonPkgConfigDependency(PkgConfigDependency, _PythonDependencyBase):

    # name is needed for polymorphism
    def __init__(self, name: str, environment: Environment, kwargs: DependencyObjectKWs,
                 installation: 'BasicPythonExternalProgram'):
        embed = kwargs.get('embed', False)
        pkg_embed = '-embed' if embed and mesonlib.version_compare(installation.info['version'], '>=3.8') else ''
        pkg_name = f'python-{installation.version}{pkg_embed}'

        if installation.build_config:
            pkg_libdir = installation.build_config.get('c_api.pkgconfig_path')
            pkg_libdir_origin = 'c_api.pkgconfig_path from the Python build config'
        else:
            pkg_libdir = installation.info['variables'].get('LIBPC')
            pkg_libdir_origin = 'LIBPC'
        if pkg_libdir is None:
            # we do not fall back to system directories, since this could lead
            # to using pkg-config of another Python installation, for example
            # we could end up using CPython .pc file for PyPy
            mlog.debug(f'Skipping pkgconfig lookup, {pkg_libdir_origin} is unset')
            self.is_found = False
            return

        for_machine = kwargs['native']
        sysroot = environment.properties[for_machine].get_sys_root()
        if sysroot and not path_is_in_root(Path(pkg_libdir), Path(sysroot)):
            pkg_libdir = destdir_join(sysroot, pkg_libdir)

        mlog.debug(f'Searching for {pkg_libdir!r} via pkgconfig lookup in {pkg_libdir_origin}')
        pkgconfig_paths = [pkg_libdir] if pkg_libdir else []

        PkgConfigDependency.__init__(self, pkg_name, environment, kwargs, extra_paths=pkgconfig_paths)
        _PythonDependencyBase.__init__(self, installation, embed)

        if pkg_libdir and not self.is_found:
            mlog.debug(f'{pkg_name!r} could not be found in {pkg_libdir_origin}, '
                       'this is likely due to a relocated python installation')
            return

        # pkg-config files are usually accurate starting with python 3.8
        if not self.link_libpython and mesonlib.version_compare(self.version, '< 3.8'):
            self.link_args = []

        # But not Apple, because it's a framework
        if self.env.machines.host.is_darwin() and 'PYTHONFRAMEWORKPREFIX' in self.variables:
            framework_prefix = self.variables['PYTHONFRAMEWORKPREFIX']
            # Add rpath, will be de-duplicated if necessary
            if framework_prefix.startswith('/Applications/Xcode.app/'):
                self.link_args += ['-Wl,-rpath,' + framework_prefix]
                if self.raw_link_args is not None:
                    # When None, self.link_args is used
                    self.raw_link_args += ['-Wl,-rpath,' + framework_prefix]


class PythonFrameworkDependency(ExtraFrameworkDependency, _PythonDependencyBase):

    def __init__(self, name: str, environment: 'Environment',
                 kwargs: DependencyObjectKWs, installation: 'BasicPythonExternalProgram'):
        ExtraFrameworkDependency.__init__(self, name, environment, kwargs)
        _PythonDependencyBase.__init__(self, installation, kwargs.get('embed', False))


class PythonSystemDependency(SystemDependency, _PythonDependencyBase):

    def __init__(self, name: str, environment: 'Environment',
                 kwargs: DependencyObjectKWs, installation: BasicPythonExternalProgram):
        SystemDependency.__init__(self, name, environment, kwargs)
        _PythonDependencyBase.__init__(self, installation, kwargs.get('embed', False))

        # For most platforms, match pkg-config behavior. iOS is a special case;
        # check for that first, so that check takes priority over
        # `link_libpython` (which *shouldn't* be set, but just in case)
        if self.platform.startswith('ios-'):
            # iOS doesn't use link_libpython - it links with the *framework*.
            self.link_args = ['-framework', 'Python', '-F', self.variables.get('base_prefix')]
            self.is_found = True
        elif self.link_libpython:
            # link args
            if mesonlib.is_windows():
                self.find_libpy_windows(environment, limited_api=False)
            else:
                self.find_libpy(environment)
        else:
            self.is_found = True

        # compile args
        if self.build_config:
            sysroot = environment.properties[self.for_machine].get_sys_root()
            path = self.build_config['c_api']['headers']
            if sysroot and not path_is_in_root(Path(path), Path(sysroot)):
                path = destdir_join(sysroot, path)
            inc_paths = mesonlib.OrderedSet([path])
        else:
            inc_paths = mesonlib.OrderedSet([
                self.variables.get('INCLUDEPY'),
                self.paths.get('include'),
                self.paths.get('platinclude')])

        self.compile_args += ['-I' + path for path in inc_paths if path]

        # https://sourceforge.net/p/mingw-w64/mailman/message/30504611/
        # https://github.com/python/cpython/pull/100137
        if mesonlib.is_windows() and self.get_windows_python_arch().endswith('64') and mesonlib.version_compare(self.version, '<3.12'):
            self.compile_args += ['-DMS_WIN64=']

        if not self.clib_compiler.has_header('Python.h', '', extra_args=self.compile_args)[0]:
            self.is_found = False

def python_factory(env: Environment, kwargs: DependencyObjectKWs,
                   installation: T.Optional['BasicPythonExternalProgram'] = None) -> T.List['DependencyGenerator']:
    # We can't use the factory_methods decorator here, as we need to pass the
    # extra installation argument
    methods = process_method_kw({DependencyMethods.PKGCONFIG, DependencyMethods.SYSTEM}, kwargs)
    candidates: T.List['DependencyGenerator'] = []
    from_installation = installation is not None
    # When not invoked through the python module, default installation.
    if installation is None:
        installation = BasicPythonExternalProgram('python3', mesonlib.python_command)
        installation.sanity()

    if DependencyMethods.PKGCONFIG in methods:
        if from_installation:
            candidates.append(DependencyCandidate(
                functools.partial(PythonPkgConfigDependency, installation=installation),
                'python3', PythonPkgConfigDependency.type_name, arguments=(env, kwargs)))
        else:
            candidates.append(DependencyCandidate.from_dependency(
                'python3', PkgConfigDependency, (env, kwargs)))

    if DependencyMethods.SYSTEM in methods:
        # This is a unique log-tried.
        candidates.append(DependencyCandidate(
            functools.partial(PythonSystemDependency, installation=installation),
            'python', 'sysconfig', arguments=(env, kwargs)))

    if DependencyMethods.EXTRAFRAMEWORK in methods:
        nkwargs = kwargs.copy()
        if mesonlib.version_compare(installation.version, '>= 3'):
            # There is a python in /System/Library/Frameworks, but that's python 2.x,
            # Python 3 will always be in /Library
            nkwargs['paths'] = ['/Library/Frameworks']
        candidates.append(DependencyCandidate(
            functools.partial(PythonFrameworkDependency, installation=installation),
            'python', PythonPkgConfigDependency.type_name, arguments=(env, nkwargs)))

    return candidates

packages['python3'] = python_factory

packages['pybind11'] = pybind11_factory = DependencyFactory(
    'pybind11',
    [DependencyMethods.PKGCONFIG, DependencyMethods.CONFIG_TOOL, DependencyMethods.CMAKE],
    configtool=Pybind11ConfigToolDependency,
)

packages['numpy'] = numpy_factory = DependencyFactory(
    'numpy',
    [DependencyMethods.PKGCONFIG, DependencyMethods.CONFIG_TOOL],
    configtool=NumPyConfigToolDependency,
)
