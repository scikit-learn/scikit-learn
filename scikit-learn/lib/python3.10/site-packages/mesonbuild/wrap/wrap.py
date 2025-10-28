# SPDX-License-Identifier: Apache-2.0
# Copyright 2015 The Meson development team

from __future__ import annotations

from .. import mlog
import contextlib
from dataclasses import dataclass
import urllib.request
import urllib.error
import urllib.parse
import os
import hashlib
import shutil
import tempfile
import stat
import subprocess
import sys
import configparser
import time
import typing as T
import textwrap
import json
import gzip

from base64 import b64encode
from netrc import netrc
from pathlib import Path, PurePath
from functools import lru_cache

from . import WrapMode
from .. import coredata
from ..mesonlib import (
    DirectoryLock, DirectoryLockAction, quiet_git, GIT, ProgressBar, MesonException,
    windows_proof_rmtree, Popen_safe
)
from ..interpreterbase import FeatureNew
from ..interpreterbase import SubProject
from .. import mesonlib

if T.TYPE_CHECKING:
    import http.client
    from typing_extensions import Literal

    Method = Literal['meson', 'cmake', 'cargo']

try:
    # Importing is just done to check if SSL exists, so all warnings
    # regarding 'imported but unused' can be safely ignored
    import ssl  # noqa
    has_ssl = True
except ImportError:
    has_ssl = False

REQ_TIMEOUT = 30.0
WHITELIST_SUBDOMAIN = 'wrapdb.mesonbuild.com'

ALL_TYPES = ['file', 'git', 'hg', 'svn', 'redirect']

if mesonlib.is_windows():
    from ..programs import ExternalProgram
    from ..mesonlib import version_compare
    _exclude_paths: T.List[str] = []
    while True:
        _patch = ExternalProgram('patch', silent=True, exclude_paths=_exclude_paths)
        if not _patch.found():
            break
        if version_compare(_patch.get_version(), '>=2.6.1'):
            break
        _exclude_paths.append(os.path.dirname(_patch.get_path()))
    PATCH = _patch.get_path() if _patch.found() else None
else:
    PATCH = shutil.which('patch')


def whitelist_wrapdb(urlstr: str) -> urllib.parse.ParseResult:
    """ raises WrapException if not whitelisted subdomain """
    url = urllib.parse.urlparse(urlstr)
    if not url.hostname:
        raise WrapException(f'{urlstr} is not a valid URL')
    if not url.hostname.endswith(WHITELIST_SUBDOMAIN):
        raise WrapException(f'{urlstr} is not a whitelisted WrapDB URL')
    if has_ssl and not url.scheme == 'https':
        raise WrapException(f'WrapDB did not have expected SSL https url, instead got {urlstr}')
    return url

def open_wrapdburl(urlstring: str, allow_insecure: bool = False, have_opt: bool = False, allow_compression: bool = False) -> http.client.HTTPResponse:
    if have_opt:
        insecure_msg = '\n\n    To allow connecting anyway, pass `--allow-insecure`.'
    else:
        insecure_msg = ''

    def do_urlopen(url: urllib.parse.ParseResult) -> http.client.HTTPResponse:
        headers = {}
        if allow_compression:
            headers['Accept-Encoding'] = 'gzip'
        req = urllib.request.Request(urllib.parse.urlunparse(url), headers=headers)
        return T.cast('http.client.HTTPResponse', urllib.request.urlopen(req, timeout=REQ_TIMEOUT))

    url = whitelist_wrapdb(urlstring)
    if has_ssl:
        try:
            return do_urlopen(url)
        except OSError as excp:
            msg = f'WrapDB connection failed to {urlstring} with error {excp}.'
            if isinstance(excp, urllib.error.URLError) and isinstance(excp.reason, ssl.SSLCertVerificationError):
                if allow_insecure:
                    mlog.warning(f'{msg}\n\n    Proceeding without authentication.')
                else:
                    raise WrapException(f'{msg}{insecure_msg}')
            else:
                raise WrapException(msg)
    elif not allow_insecure:
        raise WrapException(f'SSL module not available in {sys.executable}: Cannot contact the WrapDB.{insecure_msg}')
    else:
        # following code is only for those without Python SSL
        mlog.warning(f'SSL module not available in {sys.executable}: WrapDB traffic not authenticated.', once=True)

    # If we got this far, allow_insecure was manually passed
    try:
        return do_urlopen(url._replace(scheme='http'))
    except OSError as excp:
        raise WrapException(f'WrapDB connection failed to {urlstring} with error {excp}')

def read_and_decompress(resp: http.client.HTTPResponse) -> bytes:
    data = resp.read()
    encoding = resp.headers['Content-Encoding']
    if encoding == 'gzip':
        return gzip.decompress(data)
    elif encoding:
        raise WrapException(f'Unexpected Content-Encoding for {resp.url}: {encoding}')
    else:
        return data

def get_releases_data(allow_insecure: bool) -> bytes:
    url = open_wrapdburl('https://wrapdb.mesonbuild.com/v2/releases.json', allow_insecure, True, True)
    return read_and_decompress(url)

@lru_cache(maxsize=None)
def get_releases(allow_insecure: bool) -> T.Dict[str, T.Any]:
    data = get_releases_data(allow_insecure)
    return T.cast('T.Dict[str, T.Any]', json.loads(data.decode()))

def update_wrap_file(wrapfile: str, name: str, new_version: str, new_revision: str, allow_insecure: bool) -> None:
    url = open_wrapdburl(f'https://wrapdb.mesonbuild.com/v2/{name}_{new_version}-{new_revision}/{name}.wrap',
                         allow_insecure, True, True)
    with open(wrapfile, 'wb') as f:
        f.write(read_and_decompress(url))

def parse_patch_url(patch_url: str) -> T.Tuple[str, str]:
    u = urllib.parse.urlparse(patch_url)
    if u.netloc != 'wrapdb.mesonbuild.com':
        raise WrapException(f'URL {patch_url} does not seems to be a wrapdb patch')
    arr = u.path.strip('/').split('/')
    if arr[0] == 'v1':
        # e.g. https://wrapdb.mesonbuild.com/v1/projects/zlib/1.2.11/5/get_zip
        return arr[-3], arr[-2]
    elif arr[0] == 'v2':
        # e.g. https://wrapdb.mesonbuild.com/v2/zlib_1.2.11-5/get_patch
        tag = arr[-2]
        _, version = tag.rsplit('_', 1)
        version, revision = version.rsplit('-', 1)
        return version, revision
    else:
        raise WrapException(f'Invalid wrapdb URL {patch_url}')

class WrapException(MesonException):
    pass

class WrapNotFoundException(WrapException):
    pass

class PackageDefinition:
    def __init__(self, name: str, subprojects_dir: str, type_: T.Optional[str] = None, values: T.Optional[T.Dict[str, str]] = None):
        self.name = name
        self.subprojects_dir = subprojects_dir
        self.type = type_
        self.values = values or {}
        self.provided_deps: T.Dict[str, T.Optional[str]] = {}
        self.provided_programs: T.List[str] = []
        self.diff_files: T.List[Path] = []
        self.wrapfile_hash: T.Optional[str] = None
        self.original_filename: T.Optional[str] = None
        self.redirected: bool = False
        self.filesdir = os.path.join(self.subprojects_dir, 'packagefiles')
        self.directory = self.values.get('directory', self.name)
        if os.path.dirname(self.directory):
            raise WrapException('Directory key must be a name and not a path')
        if 'diff_files' in self.values:
            for s in self.values['diff_files'].split(','):
                path = Path(s.strip())
                if path.is_absolute():
                    raise WrapException('diff_files paths cannot be absolute')
                if '..' in path.parts:
                    raise WrapException('diff_files paths cannot contain ".."')
                self.diff_files.append(path)
        # must be lowercase for consistency with dep=variable assignment
        self.provided_deps[self.name.lower()] = None

    @staticmethod
    def from_values(name: str, subprojects_dir: str, type_: str, values: T.Dict[str, str]) -> PackageDefinition:
        return PackageDefinition(name, subprojects_dir, type_, values)

    @staticmethod
    def from_directory(filename: str) -> PackageDefinition:
        name = os.path.basename(filename)
        subprojects_dir = os.path.dirname(filename)
        return PackageDefinition(name, subprojects_dir)

    @staticmethod
    def from_wrap_file(filename: str, subproject: SubProject = SubProject('')) -> PackageDefinition:
        config, type_, values = PackageDefinition._parse_wrap(filename)
        if 'diff_files' in values:
            FeatureNew('Wrap files with diff_files', '0.63.0').use(subproject)
        if 'patch_directory' in values:
            FeatureNew('Wrap files with patch_directory', '0.55.0').use(subproject)
        for what in ['patch', 'source']:
            if f'{what}_filename' in values and f'{what}_url' not in values:
                FeatureNew(f'Local wrap patch files without {what}_url', '0.55.0').use(subproject)

        subprojects_dir = os.path.dirname(filename)

        if type_ == 'redirect':
            # [wrap-redirect] have a `filename` value pointing to the real wrap
            # file we should parse instead. It must be relative to the current
            # wrap file location and must be in the form foo/subprojects/bar.wrap.
            fname = Path(values['filename'])
            for i, p in enumerate(fname.parts):
                if i % 2 == 0:
                    if p == '..':
                        raise WrapException('wrap-redirect filename cannot contain ".."')
                else:
                    if p != 'subprojects':
                        raise WrapException('wrap-redirect filename must be in the form foo/subprojects/bar.wrap')
            if fname.suffix != '.wrap':
                raise WrapException('wrap-redirect filename must be a .wrap file')
            fname = Path(subprojects_dir, fname)
            if not fname.is_file():
                raise WrapException(f'wrap-redirect {fname} filename does not exist')
            wrap = PackageDefinition.from_wrap_file(str(fname), subproject)
            wrap.original_filename = filename
            wrap.redirected = True
            return wrap

        name = os.path.basename(filename)[:-5]
        wrap = PackageDefinition.from_values(name, subprojects_dir, type_, values)
        wrap.original_filename = filename
        wrap.parse_provide_section(config)

        patch_url = values.get('patch_url')
        if patch_url and patch_url.startswith('https://wrapdb.mesonbuild.com/v1'):
            if name == 'sqlite':
                mlog.deprecation('sqlite wrap has been renamed to sqlite3, update using `meson wrap install sqlite3`')
            elif name == 'libjpeg':
                mlog.deprecation('libjpeg wrap has been renamed to libjpeg-turbo, update using `meson wrap install libjpeg-turbo`')
            else:
                mlog.deprecation(f'WrapDB v1 is deprecated, updated using `meson wrap update {name}`')

        with open(filename, 'r', encoding='utf-8') as file:
            wrap.wrapfile_hash = hashlib.sha256(file.read().encode('utf-8')).hexdigest()

        return wrap

    @staticmethod
    def _parse_wrap(filename: str) -> T.Tuple[configparser.ConfigParser, str, T.Dict[str, str]]:
        try:
            config = configparser.ConfigParser(interpolation=None)
            config.read(filename, encoding='utf-8')
        except configparser.Error as e:
            raise WrapException(f'Failed to parse {filename}: {e!s}')
        if len(config.sections()) < 1:
            raise WrapException(f'Missing sections in {filename}')
        wrap_section = config.sections()[0]
        if not wrap_section.startswith('wrap-'):
            raise WrapException(f'{wrap_section!r} is not a valid first section in {filename}')
        type_ = wrap_section[5:]
        if type_ not in ALL_TYPES:
            raise WrapException(f'Unknown wrap type {type_!r}')
        values = dict(config[wrap_section])
        return config, type_, values

    def parse_provide_section(self, config: configparser.ConfigParser) -> None:
        if config.has_section('provides'):
            raise WrapException('Unexpected "[provides]" section, did you mean "[provide]"?')
        if config.has_section('provide'):
            for k, v in config['provide'].items():
                if k == 'dependency_names':
                    # A comma separated list of dependency names that does not
                    # need a variable name; must be lowercase for consistency with
                    # dep=variable assignment
                    names_dict = {n.strip().lower(): None for n in v.split(',')}
                    self.provided_deps.update(names_dict)
                    continue
                if k == 'program_names':
                    # A comma separated list of program names
                    names_list = [n.strip() for n in v.split(',')]
                    self.provided_programs += names_list
                    continue
                if not v:
                    m = (f'Empty dependency variable name for {k!r} in {self.name}.wrap. '
                         'If the subproject uses meson.override_dependency() '
                         'it can be added in the "dependency_names" special key.')
                    raise WrapException(m)
                self.provided_deps[k] = v

    def get(self, key: str) -> str:
        try:
            return self.values[key]
        except KeyError:
            raise WrapException(f'Missing key {key!r} in {self.name}.wrap')

    @staticmethod
    def get_hashfile(subproject_directory: str) -> str:
        return os.path.join(subproject_directory, '.meson-subproject-wrap-hash.txt')

    def update_hash_cache(self, subproject_directory: str) -> None:
        if self.wrapfile_hash:
            with open(self.get_hashfile(subproject_directory), 'w', encoding='utf-8') as file:
                file.write(self.wrapfile_hash + '\n')

def get_directory(subdir_root: str, packagename: str) -> str:
    fname = os.path.join(subdir_root, packagename + '.wrap')
    if os.path.isfile(fname):
        wrap = PackageDefinition.from_wrap_file(fname)
        return wrap.directory
    return packagename

def verbose_git(cmd: T.List[str], workingdir: str, check: bool = False) -> bool:
    '''
    Wrapper to convert GitException to WrapException caught in interpreter.
    '''
    try:
        return mesonlib.verbose_git(cmd, workingdir, check=check)
    except mesonlib.GitException as e:
        raise WrapException(str(e))

@dataclass(eq=False)
class Resolver:
    source_dir: str
    subdir: str
    subproject: SubProject = SubProject('')
    wrap_mode: WrapMode = WrapMode.default
    wrap_frontend: bool = False
    allow_insecure: bool = False
    silent: bool = False

    def __post_init__(self) -> None:
        self.subdir_root = os.path.join(self.source_dir, self.subdir)
        self.cachedir = os.environ.get('MESON_PACKAGE_CACHE_DIR') or os.path.join(self.subdir_root, 'packagecache')
        self.wraps: T.Dict[str, PackageDefinition] = {}
        self.netrc: T.Optional[netrc] = None
        self.provided_deps: T.Dict[str, PackageDefinition] = {}
        self.provided_programs: T.Dict[str, PackageDefinition] = {}
        self.wrapdb: T.Dict[str, T.Any] = {}
        self.wrapdb_provided_deps: T.Dict[str, str] = {}
        self.wrapdb_provided_programs: T.Dict[str, str] = {}
        self.loaded_dirs: T.Set[str] = set()
        self.load_wraps()
        self.load_netrc()
        self.load_wrapdb()

    def load_netrc(self) -> None:
        try:
            self.netrc = netrc()
        except FileNotFoundError:
            return
        except Exception as e:
            mlog.warning(f'failed to process netrc file: {e}.', fatal=False)

    def load_wraps(self) -> None:
        # Load Cargo.lock at the root of source tree
        source_dir = os.path.dirname(self.subdir_root)
        if os.path.exists(os.path.join(source_dir, 'Cargo.lock')):
            from .. import cargo
            for wrap in cargo.load_wraps(source_dir, self.subdir_root):
                self.wraps[wrap.name] = wrap
        # Load subprojects/*.wrap
        if os.path.isdir(self.subdir_root):
            root, dirs, files = next(os.walk(self.subdir_root))
            for i in files:
                if not i.endswith('.wrap'):
                    continue
                fname = os.path.join(self.subdir_root, i)
                wrap = PackageDefinition.from_wrap_file(fname, self.subproject)
                self.wraps[wrap.name] = wrap
            # Add dummy package definition for directories not associated with a wrap file.
            ignore_dirs = {'packagecache', 'packagefiles'}
            for wrap in self.wraps.values():
                ignore_dirs |= {wrap.directory, wrap.name}
            for i in dirs:
                if i in ignore_dirs:
                    continue
                fname = os.path.join(self.subdir_root, i)
                wrap = PackageDefinition.from_directory(fname)
                self.wraps[wrap.name] = wrap
        # Add provided deps and programs into our lookup tables
        for wrap in self.wraps.values():
            self.add_wrap(wrap)
        self.loaded_dirs.add(self.subdir)

    def add_wrap(self, wrap: PackageDefinition) -> None:
        for k in wrap.provided_deps.keys():
            if k in self.provided_deps:
                prev_wrap = self.provided_deps[k]
                m = f'Multiple wrap files provide {k!r} dependency: {wrap.name} and {prev_wrap.name}'
                raise WrapException(m)
            self.provided_deps[k] = wrap
        for k in wrap.provided_programs:
            if k in self.provided_programs:
                prev_wrap = self.provided_programs[k]
                m = f'Multiple wrap files provide {k!r} program: {wrap.name} and {prev_wrap.name}'
                raise WrapException(m)
            self.provided_programs[k] = wrap

    def load_wrapdb(self) -> None:
        try:
            with Path(self.subdir_root, 'wrapdb.json').open('r', encoding='utf-8') as f:
                self.wrapdb = json.load(f)
        except FileNotFoundError:
            return
        for name, info in self.wrapdb.items():
            self.wrapdb_provided_deps.update({i: name for i in info.get('dependency_names', [])})
            self.wrapdb_provided_programs.update({i: name for i in info.get('program_names', [])})

    def get_from_wrapdb(self, subp_name: str) -> T.Optional[PackageDefinition]:
        info = self.wrapdb.get(subp_name)
        if not info:
            return None
        self.check_can_download()
        latest_version = info['versions'][0]
        version, revision = latest_version.rsplit('-', 1)
        url = open_wrapdburl(f'https://wrapdb.mesonbuild.com/v2/{subp_name}_{version}-{revision}/{subp_name}.wrap', allow_compression=True)
        fname = Path(self.subdir_root, f'{subp_name}.wrap')
        with fname.open('wb') as f:
            f.write(read_and_decompress(url))
        mlog.log(f'Installed {subp_name} version {version} revision {revision}')
        wrap = PackageDefinition.from_wrap_file(str(fname))
        self.wraps[wrap.name] = wrap
        self.add_wrap(wrap)
        return wrap

    def _merge_wraps(self, other_resolver: 'Resolver') -> None:
        for k, v in other_resolver.wraps.items():
            prev_wrap = self.wraps.get(v.directory)
            if prev_wrap and prev_wrap.type is None and v.type is not None:
                # This happens when a subproject has been previously downloaded
                # using a wrap from another subproject and the wrap-redirect got
                # deleted. In that case, the main project created a bare wrap
                # for the download directory, but now we have a proper wrap.
                # It also happens for wraps coming from Cargo.lock files, which
                # don't create wrap-redirect.
                del self.wraps[v.directory]
                del self.provided_deps[v.directory]
            if k not in self.wraps:
                self.wraps[k] = v
                self.add_wrap(v)

    def load_and_merge(self, subdir: str, subproject: SubProject) -> None:
        if self.wrap_mode != WrapMode.nopromote and subdir not in self.loaded_dirs:
            other_resolver = Resolver(self.source_dir, subdir, subproject, self.wrap_mode, self.wrap_frontend, self.allow_insecure, self.silent)
            self._merge_wraps(other_resolver)
            self.loaded_dirs.add(subdir)

    def find_dep_provider(self, packagename: str) -> T.Tuple[T.Optional[str], T.Optional[str]]:
        # Python's ini parser converts all key values to lowercase.
        # Thus the query name must also be in lower case.
        packagename = packagename.lower()
        wrap = self.provided_deps.get(packagename)
        if wrap:
            dep_var = wrap.provided_deps.get(packagename)
            return wrap.name, dep_var
        wrap_name = self.wrapdb_provided_deps.get(packagename)
        return wrap_name, None

    def get_varname(self, subp_name: str, depname: str) -> T.Optional[str]:
        wrap = self.wraps.get(subp_name)
        return wrap.provided_deps.get(depname) if wrap else None

    def find_program_provider(self, names: T.List[str]) -> T.Optional[str]:
        for name in names:
            wrap = self.provided_programs.get(name)
            if wrap:
                return wrap.name
            wrap_name = self.wrapdb_provided_programs.get(name)
            if wrap_name:
                return wrap_name
        return None

    def _resolve(self, packagename: str, force_method: T.Optional[Method] = None) -> T.Tuple[str, Method]:
        wrap = self.wraps.get(packagename)
        if wrap is None:
            wrap = self.get_from_wrapdb(packagename)
            if wrap is None:
                raise WrapNotFoundException(f'Neither a subproject directory nor a {packagename}.wrap file was found.')
        self.wrap = wrap
        self.directory = self.wrap.directory
        self.dirname = os.path.join(self.wrap.subprojects_dir, self.wrap.directory)
        if not os.path.exists(self.dirname):
            self.dirname = os.path.join(self.subdir_root, self.directory)
        rel_path = os.path.relpath(self.dirname, self.source_dir)

        if self.wrap.original_filename:
            # If the original wrap file is not in main project's subproject_dir,
            # write a wrap-redirect.
            basename = os.path.basename(self.wrap.original_filename)
            main_fname = os.path.join(self.subdir_root, basename)
            if self.wrap.original_filename != main_fname:
                rel = os.path.relpath(self.wrap.original_filename, self.source_dir)
                mlog.log('Using', mlog.bold(rel))
                # Write a dummy wrap file in main project that redirect to the
                # wrap we picked.
                with open(main_fname, 'w', encoding='utf-8') as f:
                    f.write(textwrap.dedent(f'''\
                        [wrap-redirect]
                        filename = {PurePath(os.path.relpath(self.wrap.original_filename, self.subdir_root)).as_posix()}
                        '''))

        # Map each supported method to a file that must exist at the root of source tree.
        methods_map: T.Dict[Method, str] = {
            'meson': 'meson.build',
            'cmake': 'CMakeLists.txt',
            'cargo': 'Cargo.toml',
        }

        # Check if this wrap forces a specific method, use meson otherwise.
        method = T.cast('T.Optional[Method]', self.wrap.values.get('method', force_method))
        if method and method not in methods_map:
            allowed_methods = ', '.join(methods_map.keys())
            raise WrapException(f'Wrap method {method!r} is not supported, must be one of: {allowed_methods}')
        if force_method and method != force_method:
            raise WrapException(f'Wrap method is {method!r} but we are trying to configure it with {force_method}')
        method = method or 'meson'

        def has_buildfile() -> bool:
            return os.path.exists(os.path.join(self.dirname, methods_map[method]))

        # The directory is there and has meson.build? Great, use it.
        if has_buildfile():
            self.validate()
            return rel_path, method

        # Check if the subproject is a git submodule
        self.resolve_git_submodule()

        if os.path.exists(self.dirname):
            if not os.path.isdir(self.dirname):
                raise WrapException('Path already exists but is not a directory')
        else:
            # Check first if we have the extracted directory in our cache. This can
            # happen for example when MESON_PACKAGE_CACHE_DIR=/usr/share/cargo/registry
            # on distros that ships Rust source code.
            # TODO: We don't currently clone git repositories into the cache
            # directory, but we should to avoid cloning multiple times the same
            # repository. In that case, we could do something smarter than
            # copy_tree() here.
            cached_directory = os.path.join(self.cachedir, self.directory)
            if os.path.isdir(cached_directory):
                self.copy_tree(cached_directory, self.dirname)
            elif self.wrap.type == 'file':
                self._get_file(packagename)
            else:
                self.check_can_download()
                if self.wrap.type == 'git':
                    self._get_git(packagename)
                elif self.wrap.type == "hg":
                    self._get_hg()
                elif self.wrap.type == "svn":
                    self._get_svn()
                else:
                    raise WrapException(f'Unknown wrap type {self.wrap.type!r}')
            try:
                self.apply_patch(packagename)
                self.apply_diff_files()
            except Exception:
                windows_proof_rmtree(self.dirname)
                raise

        if not has_buildfile():
            raise WrapException(f'Subproject exists but has no {methods_map[method]} file.')

        # At this point, the subproject has been successfully resolved for the
        # first time so save off the hash of the entire wrap file for future
        # reference.
        self.wrap.update_hash_cache(self.dirname)
        return rel_path, method

    def resolve(self, packagename: str, force_method: T.Optional[Method] = None) -> T.Tuple[str, Method]:
        try:
            with DirectoryLock(self.subdir_root, '.wraplock',
                               DirectoryLockAction.WAIT,
                               'Failed to lock subprojects directory', optional=True):
                return self._resolve(packagename, force_method)
        except FileNotFoundError:
            raise WrapNotFoundException('Attempted to resolve subproject without subprojects directory present.')

    def check_can_download(self) -> None:
        # Don't download subproject data based on wrap file if requested.
        # Git submodules are ok (see above)!
        if self.wrap_mode is WrapMode.nodownload:
            m = 'Automatic wrap-based subproject downloading is disabled'
            raise WrapException(m)

    def resolve_git_submodule(self) -> bool:
        # Is git installed? If not, we're probably not in a git repository and
        # definitely cannot try to conveniently set up a submodule.
        if not GIT:
            return False
        # Does the directory exist? Even uninitialised submodules checkout an
        # empty directory to work in
        if not os.path.isdir(self.dirname):
            return False
        # Are we in a git repository?
        ret, out = quiet_git(['rev-parse'], Path(self.dirname).parent)
        if not ret:
            return False
        # Is `dirname` a submodule?
        ret, out = quiet_git(['submodule', 'status', '.'], self.dirname)
        if not ret:
            return False
        # Submodule has not been added, add it
        if out.startswith('+'):
            mlog.warning('git submodule might be out of date')
            return True
        elif out.startswith('U'):
            raise WrapException('git submodule has merge conflicts')
        # Submodule exists, but is deinitialized or wasn't initialized
        elif out.startswith('-'):
            if verbose_git(['submodule', 'update', '--init', '.'], self.dirname):
                return True
            raise WrapException('git submodule failed to init')
        # Submodule looks fine, but maybe it wasn't populated properly. Do a checkout.
        elif out.startswith(' '):
            verbose_git(['submodule', 'update', '.'], self.dirname)
            verbose_git(['checkout', '.'], self.dirname)
            # Even if checkout failed, try building it anyway and let the user
            # handle any problems manually.
            return True
        elif out == '':
            # It is not a submodule, just a folder that exists in the main repository.
            return False
        raise WrapException(f'Unknown git submodule output: {out!r}')

    def _get_file(self, packagename: str) -> None:
        path = self._get_file_internal('source', packagename)
        extract_dir = self.subdir_root
        # Some upstreams ship packages that do not have a leading directory.
        # Create one for them.
        if 'lead_directory_missing' in self.wrap.values:
            os.mkdir(self.dirname)
            extract_dir = self.dirname
        try:
            shutil.unpack_archive(path, extract_dir)
        except OSError as e:
            raise WrapException(f'failed to unpack archive with error: {str(e)}') from e

    def _get_git(self, packagename: str) -> None:
        if not GIT:
            raise WrapException(f'Git program not found, cannot download {packagename}.wrap via git.')
        revno = self.wrap.get('revision')
        checkout_cmd = ['-c', 'advice.detachedHead=false', 'checkout', revno, '--']
        is_shallow = False
        depth_option: T.List[str] = []
        if self.wrap.values.get('depth', '') != '':
            is_shallow = True
            depth_option = ['--depth', self.wrap.values.get('depth')]
        # for some reason git only allows commit ids to be shallowly fetched by fetch not with clone
        if is_shallow and self.is_git_full_commit_id(revno):
            # git doesn't support directly cloning shallowly for commits,
            # so we follow https://stackoverflow.com/a/43136160
            verbose_git(['-c', 'init.defaultBranch=meson-dummy-branch', 'init', self.directory], self.subdir_root, check=True)
            verbose_git(['remote', 'add', 'origin', self.wrap.get('url')], self.dirname, check=True)
            revno = self.wrap.get('revision')
            verbose_git(['fetch', *depth_option, 'origin', revno], self.dirname, check=True)
            verbose_git(checkout_cmd, self.dirname, check=True)
        else:
            if not is_shallow:
                verbose_git(['clone', self.wrap.get('url'), self.directory], self.subdir_root, check=True)
                if revno.lower() != 'head':
                    if not verbose_git(checkout_cmd, self.dirname):
                        verbose_git(['fetch', self.wrap.get('url'), revno], self.dirname, check=True)
                        verbose_git(checkout_cmd, self.dirname, check=True)
            else:
                args = ['-c', 'advice.detachedHead=false', 'clone', *depth_option]
                if revno.lower() != 'head':
                    args += ['--branch', revno]
                args += [self.wrap.get('url'), self.directory]
                verbose_git(args, self.subdir_root, check=True)
        if self.wrap.values.get('clone-recursive', '').lower() == 'true':
            verbose_git(['submodule', 'update', '--init', '--checkout', '--recursive', *depth_option],
                        self.dirname, check=True)
        push_url = self.wrap.values.get('push-url')
        if push_url:
            verbose_git(['remote', 'set-url', '--push', 'origin', push_url], self.dirname, check=True)

    def validate(self) -> None:
        # This check is only for subprojects with wraps.
        if not self.wrap.wrapfile_hash:
            return

        # Retrieve original hash, if it exists.
        hashfile = self.wrap.get_hashfile(self.dirname)
        if os.path.isfile(hashfile):
            with open(hashfile, 'r', encoding='utf-8') as file:
                expected_hash = file.read().strip()
        else:
            # If stored hash doesn't exist then don't warn.
            return

        # Compare hashes and warn the user if they don't match.
        if expected_hash != self.wrap.wrapfile_hash:
            mlog.warning(f'Subproject {self.wrap.name}\'s revision may be out of date; its wrap file has changed since it was first configured')

    def is_git_full_commit_id(self, revno: str) -> bool:
        result = False
        if len(revno) in {40, 64}: # 40 for sha1, 64 for upcoming sha256
            result = all(ch in '0123456789AaBbCcDdEeFf' for ch in revno)
        return result

    def _get_hg(self) -> None:
        revno = self.wrap.get('revision')
        hg = shutil.which('hg')
        if not hg:
            raise WrapException('Mercurial program not found.')
        subprocess.check_call([hg, 'clone', self.wrap.get('url'),
                               self.directory], cwd=self.subdir_root)
        if revno.lower() != 'tip':
            subprocess.check_call([hg, 'checkout', revno],
                                  cwd=self.dirname)

    def _get_svn(self) -> None:
        revno = self.wrap.get('revision')
        svn = shutil.which('svn')
        if not svn:
            raise WrapException('SVN program not found.')
        subprocess.check_call([svn, 'checkout', '-r', revno, self.wrap.get('url'),
                               self.directory], cwd=self.subdir_root)

    def get_netrc_credentials(self, netloc: str) -> T.Optional[T.Tuple[str, str]]:
        if self.netrc is None or netloc not in self.netrc.hosts:
            return None

        login, account, password = self.netrc.authenticators(netloc)
        if account:
            login = account

        return login, password

    def get_data(self, urlstring: str) -> T.Tuple[str, str]:
        blocksize = 10 * 1024
        h = hashlib.sha256()
        tmpfile = tempfile.NamedTemporaryFile(mode='wb', dir=self.cachedir, delete=False)
        url = urllib.parse.urlparse(urlstring)
        if url.hostname and url.hostname.endswith(WHITELIST_SUBDOMAIN):
            resp = open_wrapdburl(urlstring, allow_insecure=self.allow_insecure, have_opt=self.wrap_frontend)
        elif WHITELIST_SUBDOMAIN in urlstring:
            raise WrapException(f'{urlstring} may be a WrapDB-impersonating URL')
        elif url.scheme == 'sftp':
            sftp = shutil.which('sftp')
            if sftp is None:
                raise WrapException('Scheme sftp is not available. Install sftp to enable it.')
            with tempfile.TemporaryDirectory() as workdir, \
                    tempfile.NamedTemporaryFile(mode='wb', dir=self.cachedir, delete=False) as tmpfile:
                args = []
                # Older versions of the sftp client cannot handle URLs, hence the splitting of url below
                if url.port:
                    args += ['-P', f'{url.port}']
                user = f'{url.username}@' if url.username else ''
                command = [sftp, '-o', 'KbdInteractiveAuthentication=no', *args, f'{user}{url.hostname}:{url.path[1:]}']
                subprocess.run(command, cwd=workdir, check=True)
                downloaded = os.path.join(workdir, os.path.basename(url.path))
                tmpfile.close()
                shutil.move(downloaded, tmpfile.name)
                return self.hash_file(tmpfile.name), tmpfile.name
        else:
            headers = {
                'User-Agent': f'mesonbuild/{coredata.version}',
                'Accept-Language': '*',
                'Accept-Encoding': '*',
            }
            creds = self.get_netrc_credentials(url.netloc)

            if creds is not None and '@' not in url.netloc:
                login, password = creds
                if url.scheme == 'https':
                    enc_creds = b64encode(f'{login}:{password}'.encode()).decode()
                    headers.update({'Authorization': f'Basic {enc_creds}'})
                elif url.scheme == 'ftp':
                    urlstring = urllib.parse.urlunparse(url._replace(netloc=f'{login}:{password}@{url.netloc}'))
                else:
                    mlog.warning('Meson is not going to use netrc credentials for protocols other than https/ftp',
                                 fatal=False)

            try:
                req = urllib.request.Request(urlstring, headers=headers)
                resp = urllib.request.urlopen(req, timeout=REQ_TIMEOUT)
            except OSError as e:
                mlog.log(str(e))
                raise WrapException(f'could not get {urlstring}; is the internet available?')
        with contextlib.closing(resp) as resp, tmpfile as tmpfile:
            try:
                dlsize = int(resp.info()['Content-Length'])
            except TypeError:
                dlsize = None
            if dlsize is None:
                print('Downloading file of unknown size.')
                while True:
                    block = resp.read(blocksize)
                    if block == b'':
                        break
                    h.update(block)
                    tmpfile.write(block)
                hashvalue = h.hexdigest()
                return hashvalue, tmpfile.name
            sys.stdout.flush()
            progress_bar = ProgressBar(bar_type='download', total=dlsize,
                                       desc='Downloading',
                                       disable=(self.silent or None))
            while True:
                block = resp.read(blocksize)
                if block == b'':
                    break
                h.update(block)
                tmpfile.write(block)
                progress_bar.update(len(block))
            progress_bar.close()
            hashvalue = h.hexdigest()
        return hashvalue, tmpfile.name

    def hash_file(self, path: str) -> str:
        h = hashlib.sha256()
        with open(path, 'rb') as f:
            h.update(f.read())
        return h.hexdigest()

    def check_hash(self, what: str, path: str, hash_required: bool = True) -> None:
        if what + '_hash' not in self.wrap.values and not hash_required:
            return
        expected = self.wrap.get(what + '_hash').lower()
        dhash = self.hash_file(path)
        if dhash != expected:
            raise WrapException(f'Incorrect hash for {what}:\n {expected} expected\n {dhash} actual.')

    def get_data_with_backoff(self, urlstring: str) -> T.Tuple[str, str]:
        delays = [1, 2, 4, 8, 16]
        for d in delays:
            try:
                return self.get_data(urlstring)
            except Exception as e:
                mlog.warning(f'failed to download with error: {e}. Trying after a delay...', fatal=False)
                time.sleep(d)
        return self.get_data(urlstring)

    def _download(self, what: str, ofname: str, packagename: str, fallback: bool = False) -> None:
        self.check_can_download()
        srcurl = self.wrap.get(what + ('_fallback_url' if fallback else '_url'))
        mlog.log('Downloading', mlog.bold(packagename), what, 'from', mlog.bold(srcurl))
        try:
            dhash, tmpfile = self.get_data_with_backoff(srcurl)
            expected = self.wrap.get(what + '_hash').lower()
            if dhash != expected:
                os.remove(tmpfile)
                raise WrapException(f'Incorrect hash for {what}:\n {expected} expected\n {dhash} actual.')
        except WrapException:
            if not fallback:
                if what + '_fallback_url' in self.wrap.values:
                    return self._download(what, ofname, packagename, fallback=True)
                mlog.log('A fallback URL could be specified using',
                         mlog.bold(what + '_fallback_url'), 'key in the wrap file')
            raise
        os.rename(tmpfile, ofname)

    def _get_file_internal(self, what: str, packagename: str) -> str:
        filename = self.wrap.get(what + '_filename')
        if what + '_url' in self.wrap.values:
            cache_path = os.path.join(self.cachedir, filename)

            if os.path.exists(cache_path):
                self.check_hash(what, cache_path)
                mlog.log('Using', mlog.bold(packagename), what, 'from cache.')
                return cache_path

            os.makedirs(self.cachedir, exist_ok=True)
            self._download(what, cache_path, packagename)
            return cache_path
        else:
            path = Path(self.wrap.filesdir) / filename

            if not path.exists():
                raise WrapException(f'File "{path}" does not exist')
            self.check_hash(what, path.as_posix(), hash_required=False)

            return path.as_posix()

    def apply_patch(self, packagename: str) -> None:
        if 'patch_filename' in self.wrap.values and 'patch_directory' in self.wrap.values:
            m = f'Wrap file {self.wrap.name!r} must not have both "patch_filename" and "patch_directory"'
            raise WrapException(m)
        if 'patch_filename' in self.wrap.values:
            path = self._get_file_internal('patch', packagename)
            try:
                shutil.unpack_archive(path, self.subdir_root)
            except Exception:
                with tempfile.TemporaryDirectory() as workdir:
                    shutil.unpack_archive(path, workdir)
                    self.copy_tree(workdir, self.subdir_root)
        elif 'patch_directory' in self.wrap.values:
            patch_dir = self.wrap.values['patch_directory']
            src_dir = os.path.join(self.wrap.filesdir, patch_dir)
            if not os.path.isdir(src_dir):
                raise WrapException(f'patch directory does not exist: {patch_dir}')
            self.copy_tree(src_dir, self.dirname)

    def apply_diff_files(self) -> None:
        for filename in self.wrap.diff_files:
            mlog.log(f'Applying diff file "{filename}"')
            path = Path(self.wrap.filesdir) / filename
            if not path.exists():
                raise WrapException(f'Diff file "{path}" does not exist')
            relpath = os.path.relpath(str(path), self.dirname)
            if PATCH:
                # Always pass a POSIX path to patch, because on Windows it's MSYS
                # Ignore whitespace when applying patches to workaround
                # line-ending differences
                cmd = [PATCH, '-l', '-f', '-p1', '-i', str(Path(relpath).as_posix())]
            elif GIT:
                # If the `patch` command is not available, fall back to `git
                # apply`. The `--work-tree` is necessary in case we're inside a
                # Git repository: by default, Git will try to apply the patch to
                # the repository root.
                cmd = [GIT, '--work-tree', '.', 'apply', '--ignore-whitespace', '-p1', relpath]
            else:
                raise WrapException('Missing "patch" or "git" commands to apply diff files')

            p, out, _ = Popen_safe(cmd, cwd=self.dirname, stderr=subprocess.STDOUT)
            if p.returncode != 0:
                mlog.log(out.strip())
                raise WrapException(f'Failed to apply diff file "{filename}"')

    def copy_tree(self, root_src_dir: str, root_dst_dir: str) -> None:
        """
        Copy directory tree. Overwrites also read only files.
        """
        for src_dir, _, files in os.walk(root_src_dir):
            dst_dir = src_dir.replace(root_src_dir, root_dst_dir, 1)
            if not os.path.exists(dst_dir):
                os.makedirs(dst_dir)
            for file_ in files:
                src_file = os.path.join(src_dir, file_)
                dst_file = os.path.join(dst_dir, file_)
                if os.path.exists(dst_file):
                    try:
                        os.remove(dst_file)
                    except PermissionError:
                        os.chmod(dst_file, stat.S_IWUSR)
                        os.remove(dst_file)
                shutil.copy2(src_file, dst_dir, follow_symlinks=False)
