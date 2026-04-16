# SPDX-License-Identifier: Apache-2.0
# Copyright 2013-2021 The Meson development team

from __future__ import annotations

from .base import ExternalDependency, DependencyException, DependencyTypeName
from .pkgconfig import PkgConfigDependency
from ..mesonlib import (Popen_safe, join_args, version_compare, version_compare_many)
from ..options import OptionKey
from ..programs import ExternalProgram
from .. import mlog
from enum import Enum
import re
import os
import json
import typing as T

if T.TYPE_CHECKING:
    from typing_extensions import TypedDict

    from ..environment import Environment
    from .base import DependencyObjectKWs

    # Definition of what `dub describe` returns (only the fields used by Meson)
    class DubDescription(TypedDict):
        platform: T.List[str]
        architecture: T.List[str]
        buildType: str
        packages: T.List[DubPackDesc]
        targets: T.List[DubTargetDesc]

    class DubPackDesc(TypedDict):
        name: str
        version: str
        active: bool
        configuration: str
        path: str
        targetType: str
        targetFileName: str

    class DubTargetDesc(TypedDict):
        rootPackage: str
        linkDependencies: T.List[str]
        buildSettings: DubBuildSettings
        cacheArtifactPath: str

    class DubBuildSettings(TypedDict):
        importPaths: T.List[str]
        stringImportPaths: T.List[str]
        versions: T.List[str]
        mainSourceFile: str
        sourceFiles: T.List[str]
        dflags: T.List[str]
        libs: T.List[str]
        lflags: T.List[str]

    class FindTargetEntry(TypedDict):
        search: str
        artifactPath: str

class DubDescriptionSource(Enum):
    Local = 'local'
    External = 'external'

class DubDependency(ExternalDependency):
    # dub program and version
    class_dubbin: T.Optional[T.Tuple[ExternalProgram, str]] = None
    class_dubbin_searched = False
    class_cache_dir = ''

    type_name = DependencyTypeName('dub')

    # Map Meson Compiler ID's to Dub Compiler ID's
    _ID_MAP: T.Mapping[str, str] = {
        'dmd': 'dmd',
        'gcc': 'gdc',
        'llvm': 'ldc',
    }

    def __init__(self, name: str, environment: 'Environment', kwargs: DependencyObjectKWs):
        kwargs['language'] = 'd'
        super().__init__(name, environment, kwargs)
        from ..compilers.d import DCompiler, d_feature_args

        _temp_comp = super().get_compiler()
        assert isinstance(_temp_comp, DCompiler)
        self.compiler = _temp_comp

        if kwargs.get('required') is not None:
            self.required = kwargs['required']

        if DubDependency.class_dubbin is None and not DubDependency.class_dubbin_searched:
            DubDependency.class_dubbin = self._check_dub()
            DubDependency.class_dubbin_searched = True
        if DubDependency.class_dubbin is None:
            if self.required:
                raise DependencyException('DUB not found.')
            return

        (self.dubbin, dubver) = DubDependency.class_dubbin  # pylint: disable=unpacking-non-sequence

        assert isinstance(self.dubbin, ExternalProgram)

        # Check Dub's compatibility with Meson
        self._search_in_cache = version_compare(dubver, '<=1.31.1')
        self._use_cache_describe = version_compare(dubver, '>=1.35.0')
        self._dub_has_build_deep = version_compare(dubver, '>=1.35.0')

        if not self._search_in_cache and not self._use_cache_describe:
            if self.required:
                raise DependencyException(
                    f'DUB version {dubver} is not compatible with Meson'
                    " (can't locate artifacts in DUB's cache). Upgrade to Dub >= 1.35.")
            else:
                mlog.warning(f'DUB dependency {name} not found because Dub {dubver} '
                             "is not compatible with Meson. (Can't locate artifacts in DUB's cache)."
                             ' Upgrade to Dub >= 1.35')
            return

        mlog.debug('Determining dependency {!r} with DUB executable '
                   '{!r}'.format(name, self.dubbin.get_path()))

        # we need to know the target architecture
        dub_arch = self.compiler.arch

        # we need to know the build type as well
        dub_buildtype = str(environment.coredata.optstore.get_value_for(OptionKey('buildtype')))
        # MESON types: choices=['plain', 'debug', 'debugoptimized', 'release', 'minsize', 'custom'])),
        # DUB types: debug (default), plain, release, release-debug, release-nobounds, unittest, profile, profile-gc,
        # docs, ddox, cov, unittest-cov, syntax and custom
        if dub_buildtype == 'debugoptimized':
            dub_buildtype = 'release-debug'
        elif dub_buildtype == 'minsize':
            dub_buildtype = 'release'

        result = self._get_dub_description(dub_arch, dub_buildtype)
        if result is None:
            return
        description, build_cmd, description_source = result
        dub_comp_id = self._ID_MAP[self.compiler.get_id()]

        self.compile_args = []
        self.link_args = self.raw_link_args = []

        show_buildtype_warning = False

        # collect all targets
        targets = {t['rootPackage']: t for t in description['targets']}

        def find_package_target(pkg: DubPackDesc) -> bool:
            nonlocal show_buildtype_warning
            # try to find a static library in a DUB folder corresponding to
            # version, configuration, compiler, arch and build-type
            # if can find, add to link_args.
            # link_args order is meaningful, so this function MUST be called in the right order
            pack_id = f'{pkg["name"]}@{pkg["version"]}'
            tgt_desc = targets[pkg['name']]
            (tgt_file, compatibilities) = self._find_target_in_cache(description, pkg, tgt_desc, dub_comp_id)
            if tgt_file is None:
                if not compatibilities:
                    mlog.error(mlog.bold(pack_id), 'not found')
                elif 'compiler' not in compatibilities:
                    mlog.error(mlog.bold(pack_id), 'found but not compiled with ', mlog.bold(dub_comp_id))
                elif dub_comp_id != 'gdc' and 'compiler_version' not in compatibilities:
                    mlog.error(mlog.bold(pack_id), 'found but not compiled with',
                               mlog.bold(f'{dub_comp_id}-{self.compiler.version}'))
                elif 'arch' not in compatibilities:
                    mlog.error(mlog.bold(pack_id), 'found but not compiled for', mlog.bold(dub_arch))
                elif 'platform' not in compatibilities:
                    mlog.error(mlog.bold(pack_id), 'found but not compiled for',
                               mlog.bold('.'.join(description['platform'])))
                elif 'configuration' not in compatibilities:
                    mlog.error(mlog.bold(pack_id), 'found but not compiled for the',
                               mlog.bold(pkg['configuration']), 'configuration')
                else:
                    mlog.error(mlog.bold(pack_id), 'not found')

                mlog.log('You may try the following command to install the necessary DUB libraries:')
                mlog.log(mlog.bold(build_cmd))

                return False

            if 'build_type' not in compatibilities:
                mlog.warning(mlog.bold(pack_id), 'found but not compiled as', mlog.bold(dub_buildtype))
                show_buildtype_warning = True

            self.link_args.append(tgt_file)
            return True

        # Main algorithm:
        # 1. Ensure that the target is a compatible library type (not dynamic)
        # 2. Find a compatible built library for the main dependency
        # 3. Do the same for each sub-dependency.
        #    link_args MUST be in the same order than the "linkDependencies" of the main target
        # 4. Add other build settings (imports, versions etc.)

        # 1
        packages: T.Dict[str, DubPackDesc] = {}
        found_it = False
        for pkg in description['packages']:
            packages[pkg['name']] = pkg

            if not pkg['active']:
                continue

            # check that the main dependency is indeed a library
            if pkg['name'] == name:
                if pkg['targetType'] not in ['library', 'sourceLibrary', 'staticLibrary']:
                    mlog.error(mlog.bold(name), "found but it isn't a static library, it is:",
                               pkg['targetType'])
                    return

                if self.version_reqs is not None:
                    ver = pkg['version']
                    if not version_compare_many(ver, self.version_reqs)[0]:
                        mlog.error(mlog.bold(f'{name}@{ver}'),
                                   'does not satisfy all version requirements of:',
                                   ' '.join(self.version_reqs))
                        return

                found_it = True
                self.version = pkg['version']
                self.pkg = pkg

        if not found_it:
            mlog.error(f'Could not find {name} in DUB description.')
            if description_source is DubDescriptionSource.Local:
                mlog.log('Make sure that the dependency is registered for your dub project by running:')
                mlog.log(mlog.bold(f'dub add {name}'))
            elif description_source is DubDescriptionSource.External:
                # `dub describe pkg` did not contain the pkg
                raise RuntimeError(f'`dub describe` succeeded but it does not contains {name}')
            return

        if name not in targets:
            if self.pkg['targetType'] == 'sourceLibrary':
                # source libraries have no associated targets,
                # but some build settings like import folders must be found from the package object.
                # Current algo only get these from "buildSettings" in the target object.
                # Let's save this for a future PR.
                # (See openssl DUB package for example of sourceLibrary)
                mlog.error('DUB targets of type', mlog.bold('sourceLibrary'), 'are not supported.')
            else:
                mlog.error('Could not find target description for', mlog.bold(self.name))
            return

        # Current impl only supports static libraries
        self.static = True

        # 2
        if not find_package_target(self.pkg):
            return

        # 3
        for link_dep in targets[name]['linkDependencies']:
            pkg = packages[link_dep]
            if not find_package_target(pkg):
                return

        if show_buildtype_warning:
            mlog.log('If it is not suitable, try the following command and reconfigure Meson with', mlog.bold('--clearcache'))
            mlog.log(mlog.bold(build_cmd))

        # 4
        bs = targets[name]['buildSettings']

        for flag in bs['dflags']:
            self.compile_args.append(flag)

        for path in bs['importPaths']:
            self.compile_args.append('-I' + path)

        for path in bs['stringImportPaths']:
            if 'import_dir' not in d_feature_args[self.compiler.id]:
                break
            flag = d_feature_args[self.compiler.id]['import_dir']
            self.compile_args.append(f'{flag}={path}')

        for ver in bs['versions']:
            if 'version' not in d_feature_args[self.compiler.id]:
                break
            flag = d_feature_args[self.compiler.id]['version']
            self.compile_args.append(f'{flag}={ver}')

        if bs['mainSourceFile']:
            self.compile_args.append(bs['mainSourceFile'])

        # pass static libraries
        # linkerFiles are added during step 3
        # for file in bs['linkerFiles']:
        #     self.link_args.append(file)

        for file in bs['sourceFiles']:
            # sourceFiles may contain static libraries
            if file.endswith('.lib') or file.endswith('.a'):
                self.link_args.append(file)

        for flag in bs['lflags']:
            self.link_args.append(flag)

        is_windows = self.env.machines.host.is_windows()
        if is_windows:
            winlibs = ['kernel32', 'user32', 'gdi32', 'winspool', 'shell32', 'ole32',
                       'oleaut32', 'uuid', 'comdlg32', 'advapi32', 'ws2_32']

        for lib in bs['libs']:
            if os.name != 'nt':
                # trying to add system libraries by pkg-config
                pkgdep = PkgConfigDependency(lib, environment, {'required': True, 'silent': True, 'native': self.for_machine})
                if pkgdep.is_found:
                    for arg in pkgdep.get_compile_args():
                        self.compile_args.append(arg)
                    for arg in pkgdep.get_link_args():
                        self.link_args.append(arg)
                    for arg in pkgdep.get_link_args(raw=True):
                        self.raw_link_args.append(arg)
                    continue

            if is_windows and lib in winlibs:
                self.link_args.append(lib + '.lib')
                continue

            # fallback
            self.link_args.append('-l'+lib)

        self.is_found = True

    # Get the dub description needed to resolve the dependency and a
    # build command that can be used to build the dependency in case it is
    # not present.
    def _get_dub_description(self, dub_arch: str, dub_buildtype: str) -> T.Optional[T.Tuple[DubDescription, str, DubDescriptionSource]]:
        def get_build_command() -> T.List[str]:
            if self._dub_has_build_deep:
                cmd = ['dub', 'build', '--deep']
            else:
                cmd = ['dub', 'run', '--yes', 'dub-build-deep', '--']

            return cmd + [
                '--arch=' + dub_arch,
                '--compiler=' + self.compiler.get_exelist()[-1],
                '--build=' + dub_buildtype,
            ]

        # Ask dub for the package
        describe_cmd = [
            'describe', '--arch=' + dub_arch,
            '--build=' + dub_buildtype, '--compiler=' + self.compiler.get_exelist()[-1]
        ]
        helper_build = join_args(get_build_command())
        source = DubDescriptionSource.Local
        ret, res, err = self._call_dubbin(describe_cmd)
        if ret == 0:
            return (json.loads(res), helper_build, source)

        pack_spec = self.name
        if self.version_reqs is not None:
            if len(self.version_reqs) > 1:
                mlog.error('Multiple version requirements are not supported for raw dub dependencies.')
                mlog.error("Please specify only an exact version like '1.2.3'")
                raise DependencyException('Multiple version requirements are not solvable for raw dub depencies')
            elif len(self.version_reqs) == 1:
                pack_spec += '@' + self.version_reqs[0]

        describe_cmd = [
            'describe', pack_spec, '--arch=' + dub_arch,
            '--build=' + dub_buildtype, '--compiler=' + self.compiler.get_exelist()[-1]
        ]
        helper_build = join_args(get_build_command() + [pack_spec])
        source = DubDescriptionSource.External
        ret, res, err = self._call_dubbin(describe_cmd)
        if ret == 0:
            return (json.loads(res), helper_build, source)

        mlog.debug('DUB describe failed: ' + err)
        if 'locally' in err:
            mlog.error(mlog.bold(pack_spec), 'is not present locally. You may try the following command:')
            mlog.log(mlog.bold(helper_build))
        return None

    # This function finds the target of the provided JSON package, built for the right
    # compiler, architecture, configuration...
    # It returns (target|None, {compatibilities})
    # If None is returned for target, compatibilities will list what other targets were found without full compatibility
    def _find_target_in_cache(self, desc: DubDescription, pkg_desc: DubPackDesc,
                              tgt_desc: DubTargetDesc, dub_comp_id: str
                              ) -> T.Tuple[T.Optional[str], T.Set[str]]:
        mlog.debug('Searching in DUB cache for compatible', pkg_desc['targetFileName'])

        # recent DUB versions include a direct path to a compatible cached artifact
        if self._use_cache_describe:
            tgt_file = tgt_desc['cacheArtifactPath']
            if os.path.exists(tgt_file):
                return (tgt_file, {'configuration', 'platform', 'arch', 'compiler', 'compiler_version', 'build_type'})
            else:
                return (None, set())

        assert self._search_in_cache

        # try to find a string like library-debug-linux.posix-x86_64-ldc_2081-EF934983A3319F8F8FF2F0E107A363BA

        # fields are:
        #  - configuration
        #  - build type
        #  - platform
        #  - architecture
        #  - compiler id (dmd, ldc, gdc)
        #  - compiler version or frontend id or frontend version?

        comp_versions = self._get_comp_versions_to_find(dub_comp_id)

        # build_type is not in check_list because different build types might be compatible.
        # We do show a WARNING that the build type is not the same.
        # It might be critical in release builds, and acceptable otherwise
        check_list = {'configuration', 'platform', 'arch', 'compiler', 'compiler_version'}
        compatibilities: T.Set[str] = set()

        for entry in self._cache_entries(pkg_desc):
            target = entry['artifactPath']
            if not os.path.exists(target):
                # unless Dub and Meson are racing, the target file should be present
                # when the directory is present
                mlog.debug("WARNING: Could not find a Dub target: " + target)
                continue

            # we build a new set for each entry, because if this target is returned
            # we want to return only the compatibilities associated to this target
            # otherwise we could miss the WARNING about build_type
            comps: T.Set[str] = set()

            search = entry['search']

            mlog.debug('searching compatibility in ' + search)
            mlog.debug('compiler_versions', comp_versions)

            if pkg_desc['configuration'] in search:
                comps.add('configuration')

            if desc['buildType'] in search:
                comps.add('build_type')

            if all(platform in search for platform in desc['platform']):
                comps.add('platform')

            if all(arch in search for arch in desc['architecture']):
                comps.add('arch')

            if dub_comp_id in search:
                comps.add('compiler')

            if not comp_versions or any(cv in search for cv in comp_versions):
                comps.add('compiler_version')

            if check_list.issubset(comps):
                mlog.debug('Found', target)
                return (target, comps)
            else:
                compatibilities = set.union(compatibilities, comps)

        return (None, compatibilities)

    def _cache_entries(self, pkg_desc: DubPackDesc) -> T.List[FindTargetEntry]:
        # the "old" cache is the `.dub` directory in every package of ~/.dub/packages
        dub_build_path = os.path.join(pkg_desc['path'], '.dub', 'build')

        if not os.path.exists(dub_build_path):
            mlog.warning('No such cache folder:', dub_build_path)
            return []

        mlog.debug('Checking in DUB cache folder', dub_build_path)

        return [
            {
                'search': dir_entry,
                'artifactPath': os.path.join(dub_build_path, dir_entry, pkg_desc['targetFileName'])
            }
            for dir_entry in os.listdir(dub_build_path)
        ]

    def _get_comp_versions_to_find(self, dub_comp_id: str) -> T.List[str]:
        # Get D frontend version implemented in the compiler, or the compiler version itself
        # gdc doesn't support this

        if dub_comp_id == 'gdc':
            return []

        comp_versions = [self.compiler.version]

        ret, res = self._call_compbin(['--version'])[0:2]
        if ret != 0:
            mlog.error('Failed to run', mlog.bold(' '.join(self.dubbin.get_command() + ['--version'])))
            return []
        d_ver_reg = re.search('v[0-9].[0-9][0-9][0-9].[0-9]', res)  # Ex.: v2.081.2

        if d_ver_reg is not None:
            frontend_version = d_ver_reg.group()
            frontend_id = frontend_version.rsplit('.', 1)[0].replace(
                'v', '').replace('.', '')  # Fix structure. Ex.: 2081
            comp_versions.extend([frontend_version, frontend_id])

        return comp_versions

    def _call_dubbin(self, args: T.List[str], env: T.Optional[T.Dict[str, str]] = None) -> T.Tuple[int, str, str]:
        assert isinstance(self.dubbin, ExternalProgram)
        p, out, err = Popen_safe(self.dubbin.get_command() + args, env=env, cwd=self.env.get_source_dir())
        return p.returncode, out.strip(), err.strip()

    def _call_compbin(self, args: T.List[str], env: T.Optional[T.Dict[str, str]] = None) -> T.Tuple[int, str, str]:
        p, out, err = Popen_safe(self.compiler.get_exelist() + args, env=env)
        return p.returncode, out.strip(), err.strip()

    def _check_dub(self) -> T.Optional[T.Tuple[ExternalProgram, str]]:

        def find() -> T.Optional[T.Tuple[ExternalProgram, str]]:
            dubbin = ExternalProgram('dub', silent=True)

            if not dubbin.found():
                return None

            try:
                p, out = Popen_safe(dubbin.get_command() + ['--version'])[0:2]
                if p.returncode != 0:
                    mlog.warning('Found dub {!r} but couldn\'t run it'
                                 ''.format(' '.join(dubbin.get_command())))
                    return None

            except (FileNotFoundError, PermissionError):
                return None

            vermatch = re.search(r'DUB version (\d+\.\d+\.\d+.*), ', out.strip())
            if vermatch:
                dubver = vermatch.group(1)
            else:
                mlog.warning(f"Found dub {' '.join(dubbin.get_command())} but couldn't parse version in {out.strip()}")
                return None

            return (dubbin, dubver)

        found = find()

        if found is None:
            mlog.log('Found DUB:', mlog.red('NO'))
        else:
            (dubbin, dubver) = found
            mlog.log('Found DUB:', mlog.bold(dubbin.get_path()),
                     '(version %s)' % dubver)

        return found
