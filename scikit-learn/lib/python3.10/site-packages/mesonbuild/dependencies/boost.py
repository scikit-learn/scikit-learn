# SPDX-License-Identifier: Apache-2.0
# Copyright 2013-2020 The Meson development team

from __future__ import annotations

import re
import dataclasses
import functools
import typing as T
from pathlib import Path

from .. import mlog
from .. import mesonlib
from ..options import OptionKey

from .base import DependencyException, SystemDependency
from .detect import packages
from .pkgconfig import PkgConfigDependency
from .misc import threads_factory

if T.TYPE_CHECKING:
    from ..envconfig import Properties
    from ..environment import Environment

# On windows 3 directory layouts are supported:
# * The default layout (versioned) installed:
#   - $BOOST_ROOT/include/boost-x_x/boost/*.hpp
#   - $BOOST_ROOT/lib/*.lib
# * The non-default layout (system) installed:
#   - $BOOST_ROOT/include/boost/*.hpp
#   - $BOOST_ROOT/lib/*.lib
# * The pre-built binaries from sf.net:
#   - $BOOST_ROOT/boost/*.hpp
#   - $BOOST_ROOT/lib<arch>-<compiler>/*.lib where arch=32/64 and compiler=msvc-14.1
#
# Note that we should also try to support:
# mingw-w64 / Windows : libboost_<module>-mt.a            (location = <prefix>/mingw64/lib/)
#                       libboost_<module>-mt.dll.a
#
# The `modules` argument accept library names. This is because every module that
# has libraries to link against also has multiple options regarding how to
# link. See for example:
# * http://www.boost.org/doc/libs/1_65_1/libs/test/doc/html/boost_test/usage_variants.html
# * http://www.boost.org/doc/libs/1_65_1/doc/html/stacktrace/configuration_and_build.html
# * http://www.boost.org/doc/libs/1_65_1/libs/math/doc/html/math_toolkit/main_tr1.html

# **On Unix**, official packaged versions of boost libraries follow the following schemes:
#
# Linux / Debian:   libboost_<module>.so -> libboost_<module>.so.1.66.0
# Linux / Red Hat:  libboost_<module>.so -> libboost_<module>.so.1.66.0
# Linux / OpenSuse: libboost_<module>.so -> libboost_<module>.so.1.66.0
# Win   / Cygwin:   libboost_<module>.dll.a                                 (location = /usr/lib)
#                   libboost_<module>.a
#                   cygboost_<module>_1_64.dll                              (location = /usr/bin)
# Win   / VS:       boost_<module>-vc<ver>-mt[-gd]-<arch>-1_67.dll          (location = C:/local/boost_1_67_0)
# Mac   / homebrew: libboost_<module>.dylib + libboost_<module>-mt.dylib    (location = /usr/local/lib)
# Mac   / macports: libboost_<module>.dylib + libboost_<module>-mt.dylib    (location = /opt/local/lib)
#
# It's not clear that any other abi tags (e.g. -gd) are used in official packages.
#
# On Linux systems, boost libs have multithreading support enabled, but without the -mt tag.
#
# Boost documentation recommends using complex abi tags like "-lboost_regex-gcc34-mt-d-1_36".
# (See http://www.boost.org/doc/libs/1_66_0/more/getting_started/unix-variants.html#library-naming)
# However, its not clear that any Unix distribution follows this scheme.
# Furthermore, the boost documentation for unix above uses examples from windows like
#   "libboost_regex-vc71-mt-d-x86-1_34.lib", so apparently the abi tags may be more aimed at windows.
#
# We follow the following strategy for finding modules:
# A) Detect potential boost root directories (uses also BOOST_ROOT env var)
# B) Foreach candidate
#   1. Look for the boost headers (boost/version.pp)
#   2. Find all boost libraries
#     2.1 Add all libraries in lib*
#     2.2 Filter out non boost libraries
#     2.3 Filter the remaining libraries based on the meson requirements (static/shared, etc.)
#     2.4 Ensure that all libraries have the same boost tag (and are thus compatible)
#   3. Select the libraries matching the requested modules

@dataclasses.dataclass(eq=False, order=False)
class UnknownFileException(Exception):
    path: Path

@functools.total_ordering
class BoostIncludeDir():
    def __init__(self, path: Path, version_int: int):
        self.path = path
        self.version_int = version_int
        major = int(self.version_int / 100000)
        minor = int((self.version_int / 100) % 1000)
        patch = int(self.version_int % 100)
        self.version = f'{major}.{minor}.{patch}'
        self.version_lib = f'{major}_{minor}'

    def __repr__(self) -> str:
        return f'<BoostIncludeDir: {self.version} -- {self.path}>'

    def __lt__(self, other: object) -> bool:
        if isinstance(other, BoostIncludeDir):
            return (self.version_int, self.path) < (other.version_int, other.path)
        return NotImplemented

@functools.total_ordering
class BoostLibraryFile():
    # Python libraries are special because of the included
    # minor version in the module name.
    boost_python_libs = ['boost_python', 'boost_numpy']
    reg_python_mod_split = re.compile(r'(boost_[a-zA-Z]+)([0-9]*)')

    reg_abi_tag = re.compile(r'^s?g?y?d?p?n?$')
    reg_ver_tag = re.compile(r'^[0-9_]+$')

    def __init__(self, path: Path):
        self.path = path
        self.name = self.path.name

        # Initialize default properties
        self.static = False
        self.toolset = ''
        self.arch = ''
        self.version_lib = ''
        self.mt = True

        self.runtime_static = False
        self.runtime_debug = False
        self.python_debug = False
        self.debug = False
        self.stlport = False
        self.deprecated_iostreams = False

        # Post process the library name
        name_parts = self.name.split('.')
        self.basename = name_parts[0]
        self.suffixes = name_parts[1:]
        self.vers_raw = [x for x in self.suffixes if x.isdigit()]
        self.suffixes = [x for x in self.suffixes if not x.isdigit()]
        self.nvsuffix = '.'.join(self.suffixes)  # Used for detecting the library type
        self.nametags = self.basename.split('-')
        self.mod_name = self.nametags[0]
        if self.mod_name.startswith('lib'):
            self.mod_name = self.mod_name[3:]

        # Set library version if possible
        if len(self.vers_raw) >= 2:
            self.version_lib = '{}_{}'.format(self.vers_raw[0], self.vers_raw[1])

        # Detecting library type
        if self.nvsuffix in {'so', 'dll', 'dll.a', 'dll.lib', 'dylib'}:
            self.static = False
        elif self.nvsuffix in {'a', 'lib'}:
            self.static = True
        else:
            raise UnknownFileException(self.path)

        # boost_.lib is the dll import library
        if self.basename.startswith('boost_') and self.nvsuffix == 'lib':
            self.static = False

        # Process tags
        tags = self.nametags[1:]
        # Filter out the python version tag and fix modname
        if self.is_python_lib():
            tags = self.fix_python_name(tags)
        if not tags:
            return

        # Without any tags mt is assumed, however, an absence of mt in the name
        # with tags present indicates that the lib was built without mt support
        self.mt = False
        for i in tags:
            if i == 'mt':
                self.mt = True
            elif len(i) == 3 and i[1:] in {'32', '64'}:
                self.arch = i
            elif BoostLibraryFile.reg_abi_tag.match(i):
                self.runtime_static = 's' in i
                self.runtime_debug = 'g' in i
                self.python_debug = 'y' in i
                self.debug = 'd' in i
                self.stlport = 'p' in i
                self.deprecated_iostreams = 'n' in i
            elif BoostLibraryFile.reg_ver_tag.match(i):
                self.version_lib = i
            else:
                self.toolset = i

    def __repr__(self) -> str:
        return f'<LIB: {self.abitag} {self.mod_name:<32} {self.path}>'

    def __lt__(self, other: object) -> bool:
        if isinstance(other, BoostLibraryFile):
            return (
                self.mod_name, self.static, self.version_lib, self.arch,
                not self.mt, not self.runtime_static,
                not self.debug, self.runtime_debug, self.python_debug,
                self.stlport, self.deprecated_iostreams,
                self.name,
            ) < (
                other.mod_name, other.static, other.version_lib, other.arch,
                not other.mt, not other.runtime_static,
                not other.debug, other.runtime_debug, other.python_debug,
                other.stlport, other.deprecated_iostreams,
                other.name,
            )
        return NotImplemented

    def __eq__(self, other: object) -> bool:
        if isinstance(other, BoostLibraryFile):
            return self.name == other.name
        return NotImplemented

    def __hash__(self) -> int:
        return hash(self.name)

    @property
    def abitag(self) -> str:
        abitag = ''
        abitag += 'S' if self.static else '-'
        abitag += 'M' if self.mt else '-'
        abitag += ' '
        abitag += 's' if self.runtime_static else '-'
        abitag += 'g' if self.runtime_debug else '-'
        abitag += 'y' if self.python_debug else '-'
        abitag += 'd' if self.debug else '-'
        abitag += 'p' if self.stlport else '-'
        abitag += 'n' if self.deprecated_iostreams else '-'
        abitag += ' ' + (self.arch or '???')
        abitag += ' ' + (self.toolset or '?')
        abitag += ' ' + (self.version_lib or 'x_xx')
        return abitag

    def is_boost(self) -> bool:
        return any(self.name.startswith(x) for x in ['libboost_', 'boost_'])

    def is_python_lib(self) -> bool:
        return any(self.mod_name.startswith(x) for x in BoostLibraryFile.boost_python_libs)

    def fix_python_name(self, tags: T.List[str]) -> T.List[str]:
        # Handle the boost_python naming madness.
        # See https://github.com/mesonbuild/meson/issues/4788 for some distro
        # specific naming variations.
        other_tags: T.List[str] = []

        # Split the current modname into the base name and the version
        m_cur = BoostLibraryFile.reg_python_mod_split.match(self.mod_name)
        cur_name = m_cur.group(1)
        cur_vers = m_cur.group(2)

        # Update the current version string if the new version string is longer
        def update_vers(new_vers: str) -> None:
            nonlocal cur_vers
            new_vers = new_vers.replace('_', '')
            new_vers = new_vers.replace('.', '')
            if not new_vers.isdigit():
                return
            if len(new_vers) > len(cur_vers):
                cur_vers = new_vers

        for i in tags:
            if i.startswith('py'):
                update_vers(i[2:])
            elif i.isdigit():
                update_vers(i)
            elif len(i) >= 3 and i[0].isdigit() and i[2].isdigit() and i[1] == '.':
                update_vers(i)
            else:
                other_tags += [i]

        self.mod_name = cur_name + cur_vers
        return other_tags

    def mod_name_matches(self, mod_name: str) -> bool:
        if self.mod_name == mod_name:
            return True
        if not self.is_python_lib():
            return False

        m_cur = BoostLibraryFile.reg_python_mod_split.match(self.mod_name)
        m_arg = BoostLibraryFile.reg_python_mod_split.match(mod_name)

        if not m_cur or not m_arg:
            return False

        if m_cur.group(1) != m_arg.group(1):
            return False

        cur_vers = m_cur.group(2)
        arg_vers = m_arg.group(2)

        # Always assume python 2 if nothing is specified
        if not arg_vers:
            arg_vers = '2'

        return cur_vers.startswith(arg_vers)

    def version_matches(self, version_lib: str) -> bool:
        # If no version tag is present, assume that it fits
        if not self.version_lib or not version_lib:
            return True
        return self.version_lib == version_lib

    def arch_matches(self, arch: str) -> bool:
        # If no version tag is present, assume that it fits
        if not self.arch or not arch:
            return True
        return self.arch == arch

    def vscrt_matches(self, vscrt: str) -> bool:
        # If no vscrt tag present, assume that it fits  ['/MD', '/MDd', '/MT', '/MTd']
        if not vscrt:
            return True
        if vscrt in {'/MD', '-MD'}:
            return not self.runtime_static and not self.runtime_debug
        elif vscrt in {'/MDd', '-MDd'}:
            return not self.runtime_static and self.runtime_debug
        elif vscrt in {'/MT', '-MT'}:
            return (self.runtime_static or not self.static) and not self.runtime_debug
        elif vscrt in {'/MTd', '-MTd'}:
            return (self.runtime_static or not self.static) and self.runtime_debug

        mlog.warning(f'Boost: unknown vscrt tag {vscrt}. This may cause the compilation to fail. Please consider reporting this as a bug.', once=True)
        return True

    def get_compiler_args(self) -> T.List[str]:
        args: T.List[str] = []
        if self.mod_name in boost_libraries:
            libdef = boost_libraries[self.mod_name]
            if self.static:
                args += libdef.static
            else:
                args += libdef.shared
            if self.mt:
                args += libdef.multi
            else:
                args += libdef.single
        return args

    def get_link_args(self) -> T.List[str]:
        return [self.path.as_posix()]

class BoostDependency(SystemDependency):
    def __init__(self, environment: Environment, kwargs: T.Dict[str, T.Any]) -> None:
        super().__init__('boost', environment, kwargs, language='cpp')
        buildtype = environment.coredata.optstore.get_value_for(OptionKey('buildtype'))
        assert isinstance(buildtype, str)
        self.debug = buildtype.startswith('debug')
        self.multithreading = kwargs.get('threading', 'multi') == 'multi'

        self.boost_root: T.Optional[Path] = None
        self.explicit_static = 'static' in kwargs

        # Extract and validate modules
        self.modules: T.List[str] = mesonlib.extract_as_list(kwargs, 'modules')
        for i in self.modules:
            if not isinstance(i, str):
                raise DependencyException('Boost module argument is not a string.')
            if i.startswith('boost_'):
                raise DependencyException('Boost modules must be passed without the boost_ prefix')

        self.modules_found: T.List[str] = []
        self.modules_missing: T.List[str] = []

        # Do we need threads?
        if 'thread' in self.modules:
            if not self._add_sub_dependency(threads_factory(environment, self.for_machine, {})):
                self.is_found = False
                return

        # Try figuring out the architecture tag
        self.arch = environment.machines[self.for_machine].cpu_family
        self.arch = boost_arch_map.get(self.arch, None)

        # First, look for paths specified in a machine file
        props = self.env.properties[self.for_machine]
        if any(x in self.env.properties[self.for_machine] for x in
               ['boost_includedir', 'boost_librarydir', 'boost_root']):
            self.detect_boost_machine_file(props)
            return

        # Finally, look for paths from .pc files and from searching the filesystem
        self.detect_roots()

    def check_and_set_roots(self, roots: T.List[Path], use_system: bool) -> None:
        roots = list(mesonlib.OrderedSet(roots))
        for j in roots:
            #   1. Look for the boost headers (boost/version.hpp)
            mlog.debug(f'Checking potential boost root {j.as_posix()}')
            inc_dirs = self.detect_inc_dirs(j)
            inc_dirs = sorted(inc_dirs, reverse=True)  # Prefer the newer versions

            # Early abort when boost is not found
            if not inc_dirs:
                continue

            lib_dirs = self.detect_lib_dirs(j, use_system)
            self.is_found = self.run_check(inc_dirs, lib_dirs)
            if self.is_found:
                self.boost_root = j
                break

    def detect_boost_machine_file(self, props: 'Properties') -> None:
        """Detect boost with values in the machine file or environment.

        The machine file values are defaulted to the environment values.
        """
        # XXX: if we had a TypedDict we wouldn't need this
        incdir = props.get('boost_includedir')
        assert incdir is None or isinstance(incdir, str)
        libdir = props.get('boost_librarydir')
        assert libdir is None or isinstance(libdir, str)

        if incdir and libdir:
            inc_dir = Path(incdir)
            lib_dir = Path(libdir)

            if not inc_dir.is_absolute() or not lib_dir.is_absolute():
                raise DependencyException('Paths given for boost_includedir and boost_librarydir in machine file must be absolute')

            mlog.debug('Trying to find boost with:')
            mlog.debug(f'  - boost_includedir = {inc_dir}')
            mlog.debug(f'  - boost_librarydir = {lib_dir}')

            return self.detect_split_root(inc_dir, lib_dir)

        elif incdir or libdir:
            raise DependencyException('Both boost_includedir *and* boost_librarydir have to be set in your machine file (one is not enough)')

        rootdir = props.get('boost_root')
        # It shouldn't be possible to get here without something in boost_root
        assert rootdir

        raw_paths = mesonlib.stringlistify(rootdir)
        paths = [Path(x) for x in raw_paths]
        if paths and any(not x.is_absolute() for x in paths):
            raise DependencyException('boost_root path given in machine file must be absolute')

        self.check_and_set_roots(paths, use_system=False)

    def run_check(self, inc_dirs: T.List[BoostIncludeDir], lib_dirs: T.List[Path]) -> bool:
        mlog.debug('  - potential library dirs: {}'.format([x.as_posix() for x in lib_dirs]))
        mlog.debug('  - potential include dirs: {}'.format([x.path.as_posix() for x in inc_dirs]))

        must_have_library = ['boost_python']

        #   2. Find all boost libraries
        libs: T.List[BoostLibraryFile] = []
        for i in lib_dirs:
            libs = self.detect_libraries(i)
            if libs:
                mlog.debug(f'  - found boost library dir: {i}')
                # mlog.debug('  - raw library list:')
                # for j in libs:
                #     mlog.debug('    - {}'.format(j))
                break
        libs = sorted(set(libs))

        any_libs_found = len(libs) > 0
        if not any_libs_found:
            return False

        modules = ['boost_' + x for x in self.modules]
        for inc in inc_dirs:
            mlog.debug(f'  - found boost {inc.version} include dir: {inc.path}')
            f_libs = self.filter_libraries(libs, inc.version_lib)

            mlog.debug('  - filtered library list:')
            for j in f_libs:
                mlog.debug(f'    - {j}')

            #   3. Select the libraries matching the requested modules
            not_found_as_libs: T.List[str] = []
            selected_modules: T.List[BoostLibraryFile] = []
            for mod in modules:
                found = False
                for l in f_libs:
                    if l.mod_name_matches(mod):
                        selected_modules += [l]
                        found = True
                        break
                if not found:
                    not_found_as_libs += [mod]

            # If a lib is not found, but an include directory exists,
            # assume it is a header only module.
            not_found: T.List[str] = []
            for boost_modulename in not_found_as_libs:
                assert boost_modulename.startswith('boost_')
                if boost_modulename in must_have_library:
                    not_found.append(boost_modulename)
                    continue
                include_subdir = boost_modulename.replace('boost_', 'boost/', 1)
                headerdir_found = False
                for inc_dir in inc_dirs:
                    if (inc_dir.path / include_subdir).is_dir():
                        headerdir_found = True
                        break
                if not headerdir_found:
                    not_found.append(boost_modulename)

            # log the result
            mlog.debug('  - found:')
            comp_args: T.List[str] = []
            link_args: T.List[str] = []
            for j in selected_modules:
                c_args = j.get_compiler_args()
                l_args = j.get_link_args()
                mlog.debug('    - {:<24} link={} comp={}'.format(j.mod_name, str(l_args), str(c_args)))
                comp_args += c_args
                link_args += l_args

            comp_args = list(mesonlib.OrderedSet(comp_args))
            link_args = list(mesonlib.OrderedSet(link_args))

            self.modules_found = [x.mod_name for x in selected_modules]
            self.modules_found = [x[6:] for x in self.modules_found]
            self.modules_found = sorted(set(self.modules_found))
            self.modules_missing = not_found
            self.modules_missing = [x[6:] for x in self.modules_missing]
            self.modules_missing = sorted(set(self.modules_missing))

            # if we found all modules we are done
            if not not_found:
                self.version = inc.version
                self.compile_args = ['-I' + inc.path.as_posix()]
                self.compile_args += comp_args
                self.compile_args += self._extra_compile_args()
                self.compile_args = list(mesonlib.OrderedSet(self.compile_args))
                self.link_args = link_args
                mlog.debug(f'  - final compile args: {self.compile_args}')
                mlog.debug(f'  - final link args:    {self.link_args}')
                return True

            # in case we missed something log it and try again
            mlog.debug('  - NOT found:')
            for mod in not_found:
                mlog.debug(f'    - {mod}')

        return False

    def detect_inc_dirs(self, root: Path) -> T.List[BoostIncludeDir]:
        candidates: T.List[Path] = []
        inc_root = root / 'include'

        candidates += [root / 'boost']
        candidates += [inc_root / 'boost']
        if inc_root.is_dir():
            for i in inc_root.iterdir():
                if not i.is_dir() or not i.name.startswith('boost-'):
                    continue
                candidates += [i / 'boost']
        candidates = [x for x in candidates if x.is_dir()]
        candidates = [x / 'version.hpp' for x in candidates]
        candidates = [x for x in candidates if x.exists()]
        return [self._include_dir_from_version_header(x) for x in candidates]

    def detect_lib_dirs(self, root: Path, use_system: bool) -> T.List[Path]:
        # First check the system include paths. Only consider those within the
        # given root path

        if use_system:
            system_dirs_t = self.clib_compiler.get_library_dirs(self.env)
            system_dirs = [Path(x) for x in system_dirs_t]
            system_dirs = [x.resolve() for x in system_dirs if x.exists()]
            system_dirs = [x for x in system_dirs if mesonlib.path_is_in_root(x, root)]
            system_dirs = list(mesonlib.OrderedSet(system_dirs))

            if system_dirs:
                return system_dirs

        # No system include paths were found --> fall back to manually looking
        # for library dirs in root
        dirs: T.List[Path] = []
        subdirs: T.List[Path] = []
        for i in root.iterdir():
            if i.is_dir() and i.name.startswith('lib'):
                dirs += [i]

        # Some distros put libraries not directly inside /usr/lib but in /usr/lib/x86_64-linux-gnu
        for i in dirs:
            for j in i.iterdir():
                if j.is_dir() and j.name.endswith('-linux-gnu'):
                    subdirs += [j]

        # Filter out paths that don't match the target arch to avoid finding
        # the wrong libraries. See https://github.com/mesonbuild/meson/issues/7110
        if not self.arch:
            return dirs + subdirs

        arch_list_32 = ['32', 'i386']
        arch_list_64 = ['64']

        raw_list = dirs + subdirs
        no_arch = [x for x in raw_list if not any(y in x.name for y in arch_list_32 + arch_list_64)]

        matching_arch: T.List[Path] = []
        if '32' in self.arch:
            matching_arch = [x for x in raw_list if any(y in x.name for y in arch_list_32)]
        elif '64' in self.arch:
            matching_arch = [x for x in raw_list if any(y in x.name for y in arch_list_64)]

        return sorted(matching_arch) + sorted(no_arch)

    def filter_libraries(self, libs: T.List[BoostLibraryFile], lib_vers: str) -> T.List[BoostLibraryFile]:
        # MSVC is very picky with the library tags
        vscrt = ''
        try:
            crt_val = self.env.coredata.optstore.get_value('b_vscrt')
            assert isinstance(crt_val, str)
            buildtype = self.env.coredata.optstore.get_value('buildtype')
            assert isinstance(buildtype, str)
            vscrt = self.clib_compiler.get_crt_compile_args(crt_val, buildtype)[0]
        except (KeyError, IndexError, AttributeError):
            pass

        # mlog.debug('    - static: {}'.format(self.static))
        # mlog.debug('    - not explicit static: {}'.format(not self.explicit_static))
        # mlog.debug('    - mt: {}'.format(self.multithreading))
        # mlog.debug('    - version: {}'.format(lib_vers))
        # mlog.debug('    - arch: {}'.format(self.arch))
        # mlog.debug('    - vscrt: {}'.format(vscrt))
        libs = [x for x in libs if x.static == self.static or not self.explicit_static]
        libs = [x for x in libs if x.mt == self.multithreading]
        if not self.env.machines[self.for_machine].is_openbsd():
            libs = [x for x in libs if x.version_matches(lib_vers)]
        libs = [x for x in libs if x.arch_matches(self.arch)]
        libs = [x for x in libs if x.vscrt_matches(vscrt)]
        libs = [x for x in libs if x.nvsuffix != 'dll']  # Only link to import libraries

        # Only filter by debug when we are building in release mode. Debug
        # libraries are automatically preferred through sorting otherwise.
        if not self.debug:
            libs = [x for x in libs if not x.debug]

        # Take the abitag from the first library and filter by it. This
        # ensures that we have a set of libraries that are always compatible.
        if not libs:
            return []
        abitag = libs[0].abitag
        libs = [x for x in libs if x.abitag == abitag]

        return libs

    def detect_libraries(self, libdir: Path) -> T.List[BoostLibraryFile]:
        libs: T.Set[BoostLibraryFile] = set()
        for i in libdir.iterdir():
            if not i.is_file():
                continue
            if not any(i.name.startswith(x) for x in ['libboost_', 'boost_']):
                continue
            # Windows binaries from SourceForge ship with PDB files alongside
            # DLLs (#8325).  Ignore them.
            if i.name.endswith('.pdb'):
                continue

            try:
                libs.add(BoostLibraryFile(i.resolve()))
            except UnknownFileException as e:
                mlog.warning('Boost: ignoring unknown file {} under lib directory'.format(e.path.name))

        return [x for x in libs if x.is_boost()]  # Filter out no boost libraries

    def detect_split_root(self, inc_dir: Path, lib_dir: Path) -> None:
        boost_inc_dir = None
        for j in [inc_dir / 'version.hpp', inc_dir / 'boost' / 'version.hpp']:
            if j.is_file():
                boost_inc_dir = self._include_dir_from_version_header(j)
                break
        if not boost_inc_dir:
            self.is_found = False
            return

        self.is_found = self.run_check([boost_inc_dir], [lib_dir])

    def detect_roots(self) -> None:
        roots: T.List[Path] = []

        # Try getting the BOOST_ROOT from a boost.pc if it exists. This primarily
        # allows BoostDependency to find boost from Conan. See #5438
        try:
            boost_pc = PkgConfigDependency('boost', self.env, {'required': False})
            if boost_pc.found():
                boost_lib_dir = boost_pc.get_variable(pkgconfig='libdir')
                boost_inc_dir = boost_pc.get_variable(pkgconfig='includedir')
                if boost_lib_dir and boost_inc_dir:
                    mlog.debug('Trying to find boost with:')
                    mlog.debug(f'  - boost_includedir = {Path(boost_inc_dir)}')
                    mlog.debug(f'  - boost_librarydir = {Path(boost_lib_dir)}')

                    self.detect_split_root(Path(boost_inc_dir), Path(boost_lib_dir))
                    return
                else:
                    boost_root = boost_pc.get_variable(pkgconfig='prefix')
                    if boost_root:
                        roots += [Path(boost_root)]
        except DependencyException:
            pass

        # Add roots from system paths
        inc_paths = [Path(x) for x in self.clib_compiler.get_default_include_dirs()]
        inc_paths = [x.parent for x in inc_paths if x.exists()]
        inc_paths = [x.resolve() for x in inc_paths]
        roots += inc_paths

        m = self.env.machines[self.for_machine]
        # Add system paths
        if m.is_windows():
            # Where boost built from source actually installs it
            c_root = Path('C:/Boost')
            if c_root.is_dir():
                roots += [c_root]

            # Where boost documentation says it should be
            prog_files = Path('C:/Program Files/boost')
            # Where boost prebuilt binaries are
            local_boost = Path('C:/local')

            candidates: T.List[Path] = []
            if prog_files.is_dir():
                candidates += [*prog_files.iterdir()]
            if local_boost.is_dir():
                candidates += [*local_boost.iterdir()]

            roots += [x for x in candidates if x.name.lower().startswith('boost') and x.is_dir()]
        else:
            tmp: T.List[Path] = []

            # Add some default system paths
            if m.is_darwin():
                tmp.extend([
                    Path('/opt/homebrew/'),        # for Apple Silicon MacOS
                    Path('/usr/local/opt/boost'),  # for Intel Silicon MacOS
                ])
            tmp += [Path('/opt/local')]
            tmp += [Path('/usr/local')]
            tmp += [Path('/usr')]

            # Cleanup paths
            tmp = [x for x in tmp if x.is_dir()]
            tmp = [x.resolve() for x in tmp]
            roots += tmp

        self.check_and_set_roots(roots, use_system=True)

    def log_details(self) -> str:
        res = ''
        if self.modules_found:
            res += 'found: ' + ', '.join(self.modules_found)
        if self.modules_missing:
            if res:
                res += ' | '
            res += 'missing: ' + ', '.join(self.modules_missing)
        return res

    def log_info(self) -> str:
        if self.boost_root:
            return self.boost_root.as_posix()
        return ''

    def _include_dir_from_version_header(self, hfile: Path) -> BoostIncludeDir:
        # Extract the version with a regex. Using clib_compiler.get_define would
        # also work, however, this is slower (since it the compiler has to be
        # invoked) and overkill since the layout of the header is always the same.
        assert hfile.exists()
        raw = hfile.read_text(encoding='utf-8')
        m = re.search(r'#define\s+BOOST_VERSION\s+([0-9]+)', raw)
        if not m:
            mlog.debug(f'Failed to extract version information from {hfile}')
            return BoostIncludeDir(hfile.parents[1], 0)
        return BoostIncludeDir(hfile.parents[1], int(m.group(1)))

    def _extra_compile_args(self) -> T.List[str]:
        # BOOST_ALL_DYN_LINK should not be required with the known defines below
        return ['-DBOOST_ALL_NO_LIB']  # Disable automatic linking

packages['boost'] = BoostDependency

# See https://www.boost.org/doc/libs/1_72_0/more/getting_started/unix-variants.html#library-naming
# See https://mesonbuild.com/Reference-tables.html#cpu-families
boost_arch_map = {
    'aarch64': 'a64',
    'arc': 'a32',
    'arm': 'a32',
    'ia64': 'i64',
    'mips': 'm32',
    'mips64': 'm64',
    'ppc': 'p32',
    'ppc64': 'p64',
    'sparc': 's32',
    'sparc64': 's64',
    'x86': 'x32',
    'x86_64': 'x64',
}


####      ---- BEGIN GENERATED ----      ####
#                                           #
# Generated with tools/boost_names.py:
#  - boost version:   1.73.0
#  - modules found:   159
#  - libraries found: 43
#

class BoostLibrary():
    def __init__(self, name: str, shared: T.List[str], static: T.List[str], single: T.List[str], multi: T.List[str]):
        self.name = name
        self.shared = shared
        self.static = static
        self.single = single
        self.multi = multi

class BoostModule():
    def __init__(self, name: str, key: str, desc: str, libs: T.List[str]):
        self.name = name
        self.key = key
        self.desc = desc
        self.libs = libs


# dict of all know libraries with additional compile options
boost_libraries = {
    'boost_atomic': BoostLibrary(
        name='boost_atomic',
        shared=['-DBOOST_ATOMIC_DYN_LINK=1'],
        static=['-DBOOST_ATOMIC_STATIC_LINK=1'],
        single=[],
        multi=[],
    ),
    'boost_chrono': BoostLibrary(
        name='boost_chrono',
        shared=['-DBOOST_CHRONO_DYN_LINK=1'],
        static=['-DBOOST_CHRONO_STATIC_LINK=1'],
        single=['-DBOOST_CHRONO_THREAD_DISABLED'],
        multi=[],
    ),
    'boost_container': BoostLibrary(
        name='boost_container',
        shared=['-DBOOST_CONTAINER_DYN_LINK=1'],
        static=['-DBOOST_CONTAINER_STATIC_LINK=1'],
        single=[],
        multi=[],
    ),
    'boost_context': BoostLibrary(
        name='boost_context',
        shared=['-DBOOST_CONTEXT_DYN_LINK=1'],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_contract': BoostLibrary(
        name='boost_contract',
        shared=['-DBOOST_CONTRACT_DYN_LINK'],
        static=['-DBOOST_CONTRACT_STATIC_LINK'],
        single=['-DBOOST_CONTRACT_DISABLE_THREADS'],
        multi=[],
    ),
    'boost_coroutine': BoostLibrary(
        name='boost_coroutine',
        shared=['-DBOOST_COROUTINES_DYN_LINK=1'],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_date_time': BoostLibrary(
        name='boost_date_time',
        shared=['-DBOOST_DATE_TIME_DYN_LINK=1'],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_exception': BoostLibrary(
        name='boost_exception',
        shared=[],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_fiber': BoostLibrary(
        name='boost_fiber',
        shared=['-DBOOST_FIBERS_DYN_LINK=1'],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_fiber_numa': BoostLibrary(
        name='boost_fiber_numa',
        shared=['-DBOOST_FIBERS_DYN_LINK=1'],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_filesystem': BoostLibrary(
        name='boost_filesystem',
        shared=['-DBOOST_FILESYSTEM_DYN_LINK=1'],
        static=['-DBOOST_FILESYSTEM_STATIC_LINK=1'],
        single=[],
        multi=[],
    ),
    'boost_graph': BoostLibrary(
        name='boost_graph',
        shared=[],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_iostreams': BoostLibrary(
        name='boost_iostreams',
        shared=['-DBOOST_IOSTREAMS_DYN_LINK=1'],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_locale': BoostLibrary(
        name='boost_locale',
        shared=[],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_log': BoostLibrary(
        name='boost_log',
        shared=['-DBOOST_LOG_DYN_LINK=1'],
        static=[],
        single=['-DBOOST_LOG_NO_THREADS'],
        multi=[],
    ),
    'boost_log_setup': BoostLibrary(
        name='boost_log_setup',
        shared=['-DBOOST_LOG_SETUP_DYN_LINK=1'],
        static=[],
        single=['-DBOOST_LOG_NO_THREADS'],
        multi=[],
    ),
    'boost_math_c99': BoostLibrary(
        name='boost_math_c99',
        shared=[],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_math_c99f': BoostLibrary(
        name='boost_math_c99f',
        shared=[],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_math_c99l': BoostLibrary(
        name='boost_math_c99l',
        shared=[],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_math_tr1': BoostLibrary(
        name='boost_math_tr1',
        shared=[],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_math_tr1f': BoostLibrary(
        name='boost_math_tr1f',
        shared=[],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_math_tr1l': BoostLibrary(
        name='boost_math_tr1l',
        shared=[],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_mpi': BoostLibrary(
        name='boost_mpi',
        shared=[],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_nowide': BoostLibrary(
        name='boost_nowide',
        shared=['-DBOOST_NOWIDE_DYN_LINK=1'],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_prg_exec_monitor': BoostLibrary(
        name='boost_prg_exec_monitor',
        shared=['-DBOOST_TEST_DYN_LINK=1'],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_program_options': BoostLibrary(
        name='boost_program_options',
        shared=[],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_random': BoostLibrary(
        name='boost_random',
        shared=['-DBOOST_RANDOM_DYN_LINK'],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_regex': BoostLibrary(
        name='boost_regex',
        shared=[],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_serialization': BoostLibrary(
        name='boost_serialization',
        shared=[],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_stacktrace_addr2line': BoostLibrary(
        name='boost_stacktrace_addr2line',
        shared=[],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_stacktrace_backtrace': BoostLibrary(
        name='boost_stacktrace_backtrace',
        shared=[],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_stacktrace_basic': BoostLibrary(
        name='boost_stacktrace_basic',
        shared=[],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_stacktrace_noop': BoostLibrary(
        name='boost_stacktrace_noop',
        shared=[],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_stacktrace_windbg': BoostLibrary(
        name='boost_stacktrace_windbg',
        shared=[],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_stacktrace_windbg_cached': BoostLibrary(
        name='boost_stacktrace_windbg_cached',
        shared=[],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_system': BoostLibrary(
        name='boost_system',
        shared=['-DBOOST_SYSTEM_DYN_LINK=1'],
        static=['-DBOOST_SYSTEM_STATIC_LINK=1'],
        single=[],
        multi=[],
    ),
    'boost_test_exec_monitor': BoostLibrary(
        name='boost_test_exec_monitor',
        shared=['-DBOOST_TEST_DYN_LINK=1'],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_thread': BoostLibrary(
        name='boost_thread',
        shared=['-DBOOST_THREAD_BUILD_DLL=1', '-DBOOST_THREAD_USE_DLL=1'],
        static=['-DBOOST_THREAD_BUILD_LIB=1', '-DBOOST_THREAD_USE_LIB=1'],
        single=[],
        multi=[],
    ),
    'boost_timer': BoostLibrary(
        name='boost_timer',
        shared=['-DBOOST_TIMER_DYN_LINK=1'],
        static=['-DBOOST_TIMER_STATIC_LINK=1'],
        single=[],
        multi=[],
    ),
    'boost_type_erasure': BoostLibrary(
        name='boost_type_erasure',
        shared=['-DBOOST_TYPE_ERASURE_DYN_LINK'],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_unit_test_framework': BoostLibrary(
        name='boost_unit_test_framework',
        shared=['-DBOOST_TEST_DYN_LINK=1'],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_wave': BoostLibrary(
        name='boost_wave',
        shared=[],
        static=[],
        single=[],
        multi=[],
    ),
    'boost_wserialization': BoostLibrary(
        name='boost_wserialization',
        shared=[],
        static=[],
        single=[],
        multi=[],
    ),
}

#                                           #
####       ---- END GENERATED ----       ####
