# SPDX-License-Identifier: Apache-2.0
# Copyright 2013-2025 The Meson development team
# Copyright © 2023-2025 Intel Corporation

from __future__ import annotations

import copy

from . import mlog, options
import argparse
import pickle, os, uuid
import sys
from functools import lru_cache
from itertools import chain
from collections import OrderedDict
import textwrap

from .mesonlib import (
    MesonException, MachineChoice, PerMachine,
    PerMachineDefaultable,
    default_prefix,
    pickle_load
)

from .options import OptionKey

from .machinefile import CmdLineFileParser

import ast
import enum
import shlex
import typing as T

if T.TYPE_CHECKING:
    from typing_extensions import Protocol

    from . import dependencies
    from .compilers.compilers import Compiler, CompileResult, RunResult, CompileCheckMode
    from .dependencies.detect import TV_DepID
    from .mesonlib import FileOrString
    from .cmake.traceparser import CMakeCacheEntry
    from .interpreterbase import SubProject
    from .options import ElementaryOptionValues, MutableKeyedOptionDictType
    from .build import BuildTarget

    class SharedCMDOptions(Protocol):

        """Representation of command line options from Meson setup, configure,
        and dist.

        :param cmd_line_options: command line options parsed into an OptionKey:
            str mapping
        """

        cmd_line_options: T.Dict[OptionKey, T.Optional[str]]
        cross_file: T.List[str]
        native_file: T.List[str]

    OptionDictType = T.Dict[str, options.AnyOptionType]
    CompilerCheckCacheKey = T.Tuple[T.Tuple[str, ...], str, FileOrString, T.Tuple[str, ...], CompileCheckMode]
    # code, args
    RunCheckCacheKey = T.Tuple[str, T.Tuple[str, ...]]

    # typeshed
    StrOrBytesPath = T.Union[str, bytes, os.PathLike[str], os.PathLike[bytes]]

# Check major_versions_differ() if changing versioning scheme.
#
# Pip requires that RCs are named like this: '0.1.0.rc1'
# But the corresponding Git tag needs to be '0.1.0rc1'
version = '1.9.1'

# The next stable version when we are in dev. This is used to allow projects to
# require meson version >=1.2.0 when using 1.1.99. FeatureNew won't warn when
# using a feature introduced in 1.2.0 when using Meson 1.1.99.
stable_version = version
if stable_version.endswith('.99'):
    stable_version_array = stable_version.split('.')
    stable_version_array[-1] = '0'
    stable_version_array[-2] = str(int(stable_version_array[-2]) + 1)
    stable_version = '.'.join(stable_version_array)


def get_genvs_default_buildtype_list() -> list[str]:
    # just debug, debugoptimized, and release for now
    # but this should probably be configurable through some extra option, alongside --genvslite.
    return options.buildtypelist[1:-2]


class MesonVersionMismatchException(MesonException):
    '''Build directory generated with Meson version is incompatible with current version'''
    def __init__(self, old_version: str, current_version: str, extra_msg: str = '') -> None:
        super().__init__(f'Build directory has been generated with Meson version {old_version}, '
                         f'which is incompatible with the current version {current_version}.'
                         + extra_msg)
        self.old_version = old_version
        self.current_version = current_version


class DependencyCacheType(enum.Enum):

    OTHER = 0
    PKG_CONFIG = 1
    CMAKE = 2

    @classmethod
    def from_type(cls, dep: 'dependencies.Dependency') -> 'DependencyCacheType':
        # As more types gain search overrides they'll need to be added here
        if dep.type_name == 'pkgconfig':
            return cls.PKG_CONFIG
        if dep.type_name == 'cmake':
            return cls.CMAKE
        return cls.OTHER


class DependencySubCache:

    def __init__(self, type_: DependencyCacheType):
        self.types = [type_]
        self.__cache: T.Dict[T.Tuple[str, ...], 'dependencies.Dependency'] = {}

    def __getitem__(self, key: T.Tuple[str, ...]) -> 'dependencies.Dependency':
        return self.__cache[key]

    def __setitem__(self, key: T.Tuple[str, ...], value: 'dependencies.Dependency') -> None:
        self.__cache[key] = value

    def __contains__(self, key: T.Tuple[str, ...]) -> bool:
        return key in self.__cache

    def values(self) -> T.Iterable['dependencies.Dependency']:
        return self.__cache.values()


class DependencyCache:

    """Class that stores a cache of dependencies.

    This class is meant to encapsulate the fact that we need multiple keys to
    successfully lookup by providing a simple get/put interface.
    """

    def __init__(self, builtins: options.OptionStore, for_machine: MachineChoice):
        self.__cache: T.MutableMapping[TV_DepID, DependencySubCache] = OrderedDict()
        self.__builtins = builtins
        self.__pkg_conf_key = options.OptionKey('pkg_config_path', machine=for_machine)
        self.__cmake_key = options.OptionKey('cmake_prefix_path', machine=for_machine)

    def __calculate_subkey(self, type_: DependencyCacheType) -> T.Tuple[str, ...]:
        data: T.Dict[DependencyCacheType, T.List[str]] = {
            DependencyCacheType.PKG_CONFIG: T.cast('T.List[str]', self.__builtins.get_value_for(self.__pkg_conf_key)),
            DependencyCacheType.CMAKE: T.cast('T.List[str]', self.__builtins.get_value_for(self.__cmake_key)),
            DependencyCacheType.OTHER: [],
        }
        assert type_ in data, 'Someone forgot to update subkey calculations for a new type'
        return tuple(data[type_])

    def __iter__(self) -> T.Iterator['TV_DepID']:
        return self.keys()

    def put(self, key: 'TV_DepID', dep: 'dependencies.Dependency') -> None:
        t = DependencyCacheType.from_type(dep)
        if key not in self.__cache:
            self.__cache[key] = DependencySubCache(t)
        subkey = self.__calculate_subkey(t)
        self.__cache[key][subkey] = dep

    def get(self, key: 'TV_DepID') -> T.Optional['dependencies.Dependency']:
        """Get a value from the cache.

        If there is no cache entry then None will be returned.
        """
        try:
            val = self.__cache[key]
        except KeyError:
            return None

        for t in val.types:
            subkey = self.__calculate_subkey(t)
            try:
                return val[subkey]
            except KeyError:
                pass
        return None

    def values(self) -> T.Iterator['dependencies.Dependency']:
        for c in self.__cache.values():
            yield from c.values()

    def keys(self) -> T.Iterator['TV_DepID']:
        return iter(self.__cache.keys())

    def items(self) -> T.Iterator[T.Tuple['TV_DepID', T.List['dependencies.Dependency']]]:
        for k, v in self.__cache.items():
            vs: T.List[dependencies.Dependency] = []
            for t in v.types:
                subkey = self.__calculate_subkey(t)
                if subkey in v:
                    vs.append(v[subkey])
            yield k, vs

    def clear(self) -> None:
        self.__cache.clear()


class CMakeStateCache:
    """Class that stores internal CMake compiler states.

    This cache is used to reduce the startup overhead of CMake by caching
    all internal CMake compiler variables.
    """

    def __init__(self) -> None:
        self.__cache: T.Dict[str, T.Dict[str, T.List[str]]] = {}
        self.cmake_cache: T.Dict[str, 'CMakeCacheEntry'] = {}

    def __iter__(self) -> T.Iterator[T.Tuple[str, T.Dict[str, T.List[str]]]]:
        return iter(self.__cache.items())

    def items(self) -> T.Iterator[T.Tuple[str, T.Dict[str, T.List[str]]]]:
        return iter(self.__cache.items())

    def update(self, language: str, variables: T.Dict[str, T.List[str]]):
        if language not in self.__cache:
            self.__cache[language] = {}
        self.__cache[language].update(variables)

    @property
    def languages(self) -> T.Set[str]:
        return set(self.__cache.keys())


# Can't bind this near the class method it seems, sadly.
_V = T.TypeVar('_V')

# This class contains all data that must persist over multiple
# invocations of Meson. It is roughly the same thing as
# cmakecache.

class CoreData:

    def __init__(self, cmd_options: SharedCMDOptions, scratch_dir: str, meson_command: T.List[str]):
        self.lang_guids = {
            'default': '8BC9CEB8-8B4A-11D0-8D11-00A0C91BC942',
            'c': '8BC9CEB8-8B4A-11D0-8D11-00A0C91BC942',
            'cpp': '8BC9CEB8-8B4A-11D0-8D11-00A0C91BC942',
            'masm': '8BC9CEB8-8B4A-11D0-8D11-00A0C91BC942',
            'test': '3AC096D0-A1C2-E12C-1390-A8335801FDAB',
            'directory': '2150E333-8FDC-42A3-9474-1A3956D46DE8',
        }
        self.test_guid = str(uuid.uuid4()).upper()
        self.regen_guid = str(uuid.uuid4()).upper()
        self.install_guid = str(uuid.uuid4()).upper()
        self.meson_command = meson_command
        self.target_guids = {}
        self.version = version
        self.cross_files = self.__load_config_files(cmd_options, scratch_dir, 'cross')
        self.compilers: PerMachine[T.Dict[str, Compiler]] = PerMachine(OrderedDict(), OrderedDict())
        self.optstore = options.OptionStore(self.is_cross_build())

        # Stores the (name, hash) of the options file, The name will be either
        # "meson_options.txt" or "meson.options".
        # This is used by mconf to reload the option file if it's changed.
        self.options_files: T.Dict[SubProject, T.Optional[T.Tuple[str, str]]] = {}

        # Set of subprojects that have already been initialized once, this is
        # required to be stored and reloaded with the coredata, as we don't
        # want to overwrite options for such subprojects.
        self.initialized_subprojects: T.Set[str] = set()

        # For host == build configurations these caches should be the same.
        self.deps: PerMachine[DependencyCache] = PerMachineDefaultable.default(
            self.is_cross_build(),
            DependencyCache(self.optstore, MachineChoice.BUILD),
            DependencyCache(self.optstore, MachineChoice.HOST))

        self.compiler_check_cache: T.Dict['CompilerCheckCacheKey', 'CompileResult'] = OrderedDict()
        self.run_check_cache: T.Dict['RunCheckCacheKey', 'RunResult'] = OrderedDict()

        # CMake cache
        self.cmake_cache: PerMachine[CMakeStateCache] = PerMachine(CMakeStateCache(), CMakeStateCache())

        # Only to print a warning if it changes between Meson invocations.
        self.config_files = self.__load_config_files(cmd_options, scratch_dir, 'native')
        self.builtin_options_libdir_cross_fixup()
        self.init_builtins()

    @staticmethod
    def __load_config_files(cmd_options: SharedCMDOptions, scratch_dir: str, ftype: str) -> T.List[str]:
        # Need to try and make the passed filenames absolute because when the
        # files are parsed later we'll have chdir()d.
        if ftype == 'cross':
            filenames = cmd_options.cross_file
        else:
            filenames = cmd_options.native_file

        if not filenames:
            return []

        found_invalid: T.List[str] = []
        missing: T.List[str] = []
        real: T.List[str] = []
        for i, f in enumerate(filenames):
            f = os.path.expanduser(os.path.expandvars(f))
            if os.path.exists(f):
                if os.path.isfile(f):
                    real.append(os.path.abspath(f))
                    continue
                elif os.path.isdir(f):
                    found_invalid.append(os.path.abspath(f))
                else:
                    # in this case we've been passed some kind of pipe, copy
                    # the contents of that file into the meson private (scratch)
                    # directory so that it can be re-read when wiping/reconfiguring
                    fcopy = os.path.join(scratch_dir, f'{uuid.uuid4()}.{ftype}.ini')
                    with open(f, encoding='utf-8') as rf:
                        with open(fcopy, 'w', encoding='utf-8') as wf:
                            wf.write(rf.read())
                    real.append(fcopy)

                    # Also replace the command line argument, as the pipe
                    # probably won't exist on reconfigure
                    filenames[i] = fcopy
                    continue
            if sys.platform != 'win32':
                paths = [
                    os.environ.get('XDG_DATA_HOME', os.path.expanduser('~/.local/share')),
                ] + os.environ.get('XDG_DATA_DIRS', '/usr/local/share:/usr/share').split(':')
                for path in paths:
                    path_to_try = os.path.join(path, 'meson', ftype, f)
                    if os.path.isfile(path_to_try):
                        real.append(path_to_try)
                        break
                else:
                    missing.append(f)
            else:
                missing.append(f)

        if missing:
            if found_invalid:
                mlog.log('Found invalid candidates for', ftype, 'file:', *found_invalid)
            mlog.log('Could not find any valid candidate for', ftype, 'files:', *missing)
            raise MesonException(f'Cannot find specified {ftype} file: {f}')
        return real

    def builtin_options_libdir_cross_fixup(self) -> None:
        # By default set libdir to "lib" when cross compiling since
        # getting the "system default" is always wrong on multiarch
        # platforms as it gets a value like lib/x86_64-linux-gnu.
        if self.cross_files:
            options.BUILTIN_OPTIONS[OptionKey('libdir')].default = 'lib'

    def init_builtins(self) -> None:
        # Create builtin options with default values
        for key, opt in options.BUILTIN_OPTIONS.items():
            self.add_builtin_option(self.optstore, key, opt)
        for for_machine in iter(MachineChoice):
            for key, opt in options.BUILTIN_OPTIONS_PER_MACHINE.items():
                self.add_builtin_option(self.optstore, key.evolve(machine=for_machine), opt)

    @staticmethod
    def add_builtin_option(optstore: options.OptionStore, key: OptionKey,
                           opt: options.AnyOptionType) -> None:
        # Create a copy of the object, as we're going to mutate it
        opt = copy.copy(opt)
        if key.subproject:
            if opt.yielding:
                # This option is global and not per-subproject
                return
        else:
            new_value = options.argparse_prefixed_default(
                opt, key, default_prefix())
            opt.set_value(new_value)

        modulename = key.get_module_prefix()
        if modulename:
            optstore.add_module_option(modulename, key, opt)
        else:
            optstore.add_system_option(key, opt)

    def init_backend_options(self, backend_name: str) -> None:
        if backend_name == 'ninja':
            self.optstore.add_system_option('backend_max_links', options.UserIntegerOption(
                'backend_max_links',
                'Maximum number of linker processes to run or 0 for no '
                'limit',
                0,
                min_value=0))
        elif backend_name.startswith('vs'):
            self.optstore.add_system_option('backend_startup_project', options.UserStringOption(
                'backend_startup_project',
                'Default project to execute in Visual Studio',
                ''))

    def get_option_for_target(self, target: 'BuildTarget', key: T.Union[str, OptionKey]) -> ElementaryOptionValues:
        if isinstance(key, str):
            assert ':' not in key
            newkey = OptionKey(key, target.subproject)
        else:
            newkey = key
        if newkey.subproject != target.subproject:
            # FIXME: this should be an error. The caller needs to ensure that
            # key and target have the same subproject for consistency.
            # Now just do this to get things going.
            newkey = newkey.evolve(subproject=target.subproject)
        (option_object, value) = self.optstore.get_value_object_and_value_for(newkey)
        override = target.get_override(newkey.name)
        if override is not None:
            return option_object.validate_value(override)
        return value

    def set_from_configure_command(self, options: SharedCMDOptions) -> bool:
        return self.optstore.set_from_configure_command(options.cmd_line_options)

    def set_option(self, key: OptionKey, value, first_invocation: bool = False) -> bool:
        dirty = False
        try:
            changed = self.optstore.set_option(key, value, first_invocation)
        except KeyError:
            raise MesonException(f'Tried to set unknown builtin option {str(key)}')
        dirty |= changed

        if key.name == 'buildtype':
            dirty |= self._set_others_from_buildtype(value)

        return dirty

    def clear_cache(self) -> None:
        self.deps.host.clear()
        self.deps.build.clear()
        self.compiler_check_cache.clear()
        self.run_check_cache.clear()

    def get_nondefault_buildtype_args(self) -> T.List[T.Union[T.Tuple[str, str, str], T.Tuple[str, bool, bool]]]:
        result: T.List[T.Union[T.Tuple[str, str, str], T.Tuple[str, bool, bool]]] = []
        value = self.optstore.get_value_for('buildtype')
        if value == 'plain':
            opt = 'plain'
            debug = False
        elif value == 'debug':
            opt = '0'
            debug = True
        elif value == 'debugoptimized':
            opt = '2'
            debug = True
        elif value == 'release':
            opt = '3'
            debug = False
        elif value == 'minsize':
            opt = 's'
            debug = True
        else:
            assert value == 'custom'
            return []
        actual_opt = self.optstore.get_value_for('optimization')
        actual_debug = self.optstore.get_value_for('debug')
        if actual_opt != opt:
            result.append(('optimization', actual_opt, opt))
        if actual_debug != debug:
            result.append(('debug', actual_debug, debug))
        return result

    def _set_others_from_buildtype(self, value: str) -> bool:
        dirty = False

        if value == 'plain':
            opt = 'plain'
            debug = False
        elif value == 'debug':
            opt = '0'
            debug = True
        elif value == 'debugoptimized':
            opt = '2'
            debug = True
        elif value == 'release':
            opt = '3'
            debug = False
        elif value == 'minsize':
            opt = 's'
            debug = True
        else:
            assert value == 'custom'
            return False

        dirty |= self.optstore.set_option(OptionKey('optimization'), opt)
        dirty |= self.optstore.set_option(OptionKey('debug'), debug)

        return dirty

    def get_external_args(self, for_machine: MachineChoice, lang: str) -> T.List[str]:
        # mypy cannot analyze type of OptionKey
        key = OptionKey(f'{lang}_args', machine=for_machine)
        return T.cast('T.List[str]', self.optstore.get_value(key))

    @lru_cache(maxsize=None)
    def get_external_link_args(self, for_machine: MachineChoice, lang: str) -> T.List[str]:
        # mypy cannot analyze type of OptionKey
        linkkey = OptionKey(f'{lang}_link_args', machine=for_machine)
        return T.cast('T.List[str]', self.optstore.get_value_for(linkkey))

    def is_cross_build(self, when_building_for: MachineChoice = MachineChoice.HOST) -> bool:
        if when_building_for == MachineChoice.BUILD:
            return False
        return len(self.cross_files) > 0

    def copy_build_options_from_regular_ones(self, shut_up_pylint: bool = True) -> bool:
        # FIXME, needs cross compilation support.
        if shut_up_pylint:
            return False
        dirty = False
        assert not self.is_cross_build()
        for k in options.BUILTIN_OPTIONS_PER_MACHINE:
            o = self.optstore.get_value_object_for(k.name)
            dirty |= self.optstore.set_option(k, o.value, True)
        for bk, bv in self.optstore.items():
            if bk.machine is MachineChoice.BUILD:
                hk = bk.as_host()
                try:
                    hv = self.optstore.get_value_object(hk)
                    dirty |= bv.set_value(hv.value)
                except KeyError:
                    continue

        return dirty

    def set_options(self, opts_to_set: T.Dict[OptionKey, T.Any], subproject: str = '', first_invocation: bool = False) -> bool:
        dirty = False
        if not self.is_cross_build():
            opts_to_set = {k: v for k, v in opts_to_set.items() if k.machine is not MachineChoice.BUILD}
        # Set prefix first because it's needed to sanitize other options
        pfk = OptionKey('prefix')
        if pfk in opts_to_set:
            prefix = self.optstore.sanitize_prefix(opts_to_set[pfk])
            for key in options.BUILTIN_DIR_NOPREFIX_OPTIONS:
                if key not in opts_to_set:
                    val = options.BUILTIN_OPTIONS[key].prefixed_default(key, prefix)
                    dirty |= self.optstore.set_option(key, val)

        unknown_options: T.List[OptionKey] = []
        for k, v in opts_to_set.items():
            if k == pfk:
                continue
            elif k.evolve(subproject=None) in self.optstore:
                dirty |= self.set_option(k, v, first_invocation)
            elif k.machine != MachineChoice.BUILD and not self.optstore.is_compiler_option(k):
                unknown_options.append(k)
        if unknown_options:
            if subproject:
                # The subproject may have top-level options that should be used
                # when it is not a subproject. Ignore those for now. With option
                # refactor they will get per-subproject values.
                really_unknown = []
                for uo in unknown_options:
                    topkey = uo.as_root()
                    if topkey not in self.optstore:
                        really_unknown.append(uo)
                unknown_options = really_unknown
            if unknown_options:
                unknown_options_str = ', '.join(sorted(str(s) for s in unknown_options))
                sub = f'In subproject {subproject}: ' if subproject else ''
                raise MesonException(f'{sub}Unknown options: "{unknown_options_str}"')

        if not self.is_cross_build():
            dirty |= self.copy_build_options_from_regular_ones()

        return dirty

    def add_compiler_options(self, c_options: MutableKeyedOptionDictType, lang: str, for_machine: MachineChoice,
                             subproject: str) -> None:
        for k, o in c_options.items():
            assert k.subproject is None and k.machine is for_machine
            if subproject:
                k = k.evolve(subproject=subproject)
            if lang == 'objc' and k.name == 'c_std':
                # For objective C, always fall back to c_std.
                self.optstore.add_compiler_option('c', k, o)
            elif lang == 'objcpp' and k.name == 'cpp_std':
                self.optstore.add_compiler_option('cpp', k, o)
            else:
                self.optstore.add_compiler_option(lang, k, o)

    def process_compiler_options(self, lang: str, comp: Compiler, subproject: str) -> None:
        self.add_compiler_options(comp.get_options(), lang, comp.for_machine, subproject)

        for key in [OptionKey(f'{lang}_args'), OptionKey(f'{lang}_link_args')]:
            if self.is_cross_build():
                key = key.evolve(machine=comp.for_machine)
            # the global option is already there, but any augment is still
            # sitting in pending_options has to be taken into account
            assert key in self.optstore
            if subproject:
                skey = key.evolve(subproject=subproject)
                self.optstore.add_compiler_option(lang, skey, self.optstore.get_value_object(key))

        for key in comp.base_options:
            if subproject:
                skey = key.evolve(subproject=subproject)
            else:
                skey = key
            if skey not in self.optstore:
                self.optstore.add_system_option(skey, copy.deepcopy(options.COMPILER_BASE_OPTIONS[key]))

        self.emit_base_options_warnings()

    def emit_base_options_warnings(self) -> None:
        bcodekey = OptionKey('b_bitcode')
        if bcodekey in self.optstore and self.optstore.get_value(bcodekey):
            msg = textwrap.dedent('''Base option 'b_bitcode' is enabled, which is incompatible with many linker options.
                                     Incompatible options such as \'b_asneeded\' have been disabled.'
                                     Please see https://mesonbuild.com/Builtin-options.html#Notes_about_Apple_Bitcode_support for more details.''')
            mlog.warning(msg, once=True, fatal=False)

def get_cmd_line_file(build_dir: str) -> str:
    return os.path.join(build_dir, 'meson-private', 'cmd_line.txt')

def read_cmd_line_file(build_dir: str, options: SharedCMDOptions) -> None:
    filename = get_cmd_line_file(build_dir)
    if not os.path.isfile(filename):
        return

    config = CmdLineFileParser()
    config.read(filename)

    # Do a copy because config is not really a dict. options.cmd_line_options
    # overrides values from the file.
    d = {OptionKey.from_string(k): v for k, v in config['options'].items()}
    d.update(options.cmd_line_options)
    options.cmd_line_options = d

    properties = config['properties']
    if not options.cross_file:
        options.cross_file = ast.literal_eval(properties.get('cross_file', '[]'))
    if not options.native_file:
        # This will be a string in the form: "['first', 'second', ...]", use
        # literal_eval to get it into the list of strings.
        options.native_file = ast.literal_eval(properties.get('native_file', '[]'))

def write_cmd_line_file(build_dir: str, options: SharedCMDOptions) -> None:
    filename = get_cmd_line_file(build_dir)
    config = CmdLineFileParser()

    properties: OrderedDict[str, str] = OrderedDict()
    if options.cross_file:
        properties['cross_file'] = options.cross_file
    if options.native_file:
        properties['native_file'] = options.native_file

    config['options'] = {str(k): str(v) for k, v in options.cmd_line_options.items()}
    config['properties'] = properties
    with open(filename, 'w', encoding='utf-8') as f:
        config.write(f)

def update_cmd_line_file(build_dir: str, options: SharedCMDOptions) -> None:
    filename = get_cmd_line_file(build_dir)
    config = CmdLineFileParser()
    config.read(filename)
    for k, v in options.cmd_line_options.items():
        keystr = str(k)
        if v is not None:
            config['options'][keystr] = str(v)
        elif keystr in config['options']:
            del config['options'][keystr]

    with open(filename, 'w', encoding='utf-8') as f:
        config.write(f)

def format_cmd_line_options(options: SharedCMDOptions) -> str:
    cmdline = ['-D{}={}'.format(str(k), v) for k, v in options.cmd_line_options.items()]
    if options.cross_file:
        cmdline += [f'--cross-file={f}' for f in options.cross_file]
    if options.native_file:
        cmdline += [f'--native-file={f}' for f in options.native_file]
    return ' '.join([shlex.quote(x) for x in cmdline])

def major_versions_differ(v1: str, v2: str) -> bool:
    v1_major, v1_minor = v1.rsplit('.', 1)
    v2_major, v2_minor = v2.rsplit('.', 1)
    # Major version differ, or one is development version but not the other.
    return v1_major != v2_major or ('99' in {v1_minor, v2_minor} and v1_minor != v2_minor)

def load(build_dir: str, suggest_reconfigure: bool = True) -> CoreData:
    filename = os.path.join(build_dir, 'meson-private', 'coredata.dat')
    return pickle_load(filename, 'Coredata', CoreData, suggest_reconfigure)


def save(obj: CoreData, build_dir: str) -> str:
    filename = os.path.join(build_dir, 'meson-private', 'coredata.dat')
    prev_filename = filename + '.prev'
    tempfilename = filename + '~'
    if major_versions_differ(obj.version, version):
        raise MesonException('Fatal version mismatch corruption.')
    if os.path.exists(filename):
        import shutil
        shutil.copyfile(filename, prev_filename)
    with open(tempfilename, 'wb') as f:
        pickle.dump(obj, f)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tempfilename, filename)
    return filename


class KeyNoneAction(argparse.Action):
    """
    Custom argparse Action that stores values in a dictionary as keys with value None.
    """

    def __init__(self, option_strings, dest, nargs=None, **kwargs: object) -> None:
        assert nargs is None or nargs == 1
        super().__init__(option_strings, dest, nargs=1, **kwargs)

    def __call__(self, parser: argparse.ArgumentParser, namespace: argparse.Namespace,
                 arg: T.List[str], option_string: str = None) -> None:
        current_dict = getattr(namespace, self.dest)
        if current_dict is None:
            current_dict = {}
            setattr(namespace, self.dest, current_dict)

        key = OptionKey.from_string(arg[0])
        current_dict[key] = None


class KeyValueAction(argparse.Action):
    """
    Custom argparse Action that parses KEY=VAL arguments and stores them in a dictionary.
    """

    def __init__(self, option_strings, dest, nargs=None, **kwargs: object) -> None:
        assert nargs is None or nargs == 1
        super().__init__(option_strings, dest, nargs=1, **kwargs)

    def __call__(self, parser: argparse.ArgumentParser, namespace: argparse.Namespace,
                 arg: T.List[str], option_string: str = None) -> None:
        current_dict = getattr(namespace, self.dest)
        if current_dict is None:
            current_dict = {}
            setattr(namespace, self.dest, current_dict)

        try:
            keystr, value = arg[0].split('=', 1)
            key = OptionKey.from_string(keystr)
            current_dict[key] = value
        except ValueError:
            parser.error(f'The argument for option {option_string!r} must be in OPTION=VALUE format.')


def register_builtin_arguments(parser: argparse.ArgumentParser) -> None:
    for n, b in options.BUILTIN_OPTIONS.items():
        options.option_to_argparse(b, n, parser, '')
    for n, b in options.BUILTIN_OPTIONS_PER_MACHINE.items():
        options.option_to_argparse(b, n, parser, ' (just for host machine)')
        options.option_to_argparse(b, n.as_build(), parser, ' (just for build machine)')
    parser.add_argument('-D', action=KeyValueAction, dest='cmd_line_options', default={}, metavar="option=value",
                        help='Set the value of an option, can be used several times to set multiple options.')

def parse_cmd_line_options(args: SharedCMDOptions) -> None:
    # Merge builtin options set with --option into the dict.
    for key in chain(
            options.BUILTIN_OPTIONS.keys(),
            (k.as_build() for k in options.BUILTIN_OPTIONS_PER_MACHINE.keys()),
            options.BUILTIN_OPTIONS_PER_MACHINE.keys(),
    ):
        name = str(key)
        value = getattr(args, name, None)
        if value is not None:
            if key in args.cmd_line_options:
                cmdline_name = options.argparse_name_to_arg(name)
                raise MesonException(
                    f'Got argument {name} as both -D{name} and {cmdline_name}. Pick one.')
            args.cmd_line_options[key] = value
            delattr(args, name)


FORBIDDEN_TARGET_NAMES = frozenset({
    'clean',
    'clean-ctlist',
    'clean-gcno',
    'clean-gcda',
    'coverage',
    'coverage-text',
    'coverage-xml',
    'coverage-html',
    'phony',
    'PHONY',
    'all',
    'test',
    'benchmark',
    'install',
    'uninstall',
    'build.ninja',
    'scan-build',
    'reconfigure',
    'dist',
    'distcheck',
})
