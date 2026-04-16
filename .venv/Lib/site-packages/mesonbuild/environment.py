# SPDX-License-Identifier: Apache-2.0
# Copyright 2012-2020 The Meson development team
# Copyright © 2023-2025 Intel Corporation

from __future__ import annotations

import itertools
import os, re
import typing as T
import collections

from . import cmdline
from . import coredata
from . import mesonlib
from . import machinefile
from . import options

from .mesonlib import (
    MesonException, MachineChoice, Popen_safe, PerMachine,
    PerMachineDefaultable, PerThreeMachineDefaultable, split_args,
    MesonBugException
)
from .options import OptionKey
from . import mlog
from .programs import ExternalProgram

from .envconfig import (
    BinaryTable, MachineInfo, Properties, CMakeVariables,
    detect_machine_info, machine_info_can_run
)
from . import compilers

from mesonbuild import envconfig

if T.TYPE_CHECKING:
    from .compilers.compilers import Compiler, CompilerDict, Language
    from .options import OptionDict, ElementaryOptionValues
    from .wrap.wrap import Resolver


NON_LANG_ENV_OPTIONS = [
    ('PKG_CONFIG_PATH', 'pkg_config_path'),
    ('CMAKE_PREFIX_PATH', 'cmake_prefix_path'),
    ('LDFLAGS', 'ldflags'),
    ('CPPFLAGS', 'cppflags'),
]

build_filename = 'meson.build'


def _as_str(val: object) -> str:
    assert isinstance(val, str), 'for mypy'
    return val


def _get_env_var(for_machine: MachineChoice, is_cross: bool, var_name: str) -> T.Optional[str]:
    """
    Returns the exact env var and the value.
    """
    candidates = PerMachine(
        # The prefixed build version takes priority, but if we are native
        # compiling we fall back on the unprefixed host version. This
        # allows native builds to never need to worry about the 'BUILD_*'
        # ones.
        ([var_name + '_FOR_BUILD'] if is_cross else [var_name]),
        # Always just the unprefixed host versions
        [var_name]
    )[for_machine]
    for var in candidates:
        value = os.environ.get(var)
        if value is not None:
            break
    else:
        formatted = ', '.join([f'{var!r}' for var in candidates])
        mlog.debug(f'None of {formatted} are defined in the environment, not changing global flags.')
        return None
    mlog.debug(f'Using {var!r} from environment with value: {value!r}')
    return value


class Environment:
    private_dir = 'meson-private'
    log_dir = 'meson-logs'
    info_dir = 'meson-info'

    def __init__(self, source_dir: str, build_dir: T.Optional[str], cmd_options: cmdline.SharedCMDOptions) -> None:
        self.source_dir = source_dir
        # Do not try to create build directories when build_dir is none.
        # This reduced mode is used by the --buildoptions introspector
        if build_dir is not None:
            self.build_dir = build_dir
            self.scratch_dir = os.path.join(build_dir, Environment.private_dir)
            self.log_dir = os.path.join(build_dir, Environment.log_dir)
            self.info_dir = os.path.join(build_dir, Environment.info_dir)
            os.makedirs(self.scratch_dir, exist_ok=True)
            os.makedirs(self.log_dir, exist_ok=True)
            os.makedirs(self.info_dir, exist_ok=True)
            try:
                self.coredata: coredata.CoreData = coredata.load(self.get_build_dir(), suggest_reconfigure=False)
                self.first_invocation = False
            except FileNotFoundError:
                self.create_new_coredata(cmd_options)
            except coredata.MesonVersionMismatchException as e:
                # This is routine, but tell the user the update happened
                mlog.log('Regenerating configuration from scratch:', str(e))
                cmdline.read_cmd_line_file(self.build_dir, cmd_options)
                self.create_new_coredata(cmd_options)
            except MesonException as e:
                # If we stored previous command line options, we can recover from
                # a broken/outdated coredata.
                if os.path.isfile(cmdline.get_cmd_line_file(self.build_dir)):
                    mlog.warning('Regenerating configuration from scratch.', fatal=False)
                    mlog.log('Reason:', mlog.red(str(e)))
                    cmdline.read_cmd_line_file(self.build_dir, cmd_options)
                    self.create_new_coredata(cmd_options)
                else:
                    raise MesonException(f'{str(e)} Try regenerating using "meson setup --wipe".')
        else:
            # Just create a fresh coredata in this case
            self.build_dir = ''
            self.scratch_dir = ''
            self.create_new_coredata(cmd_options)

        ## locally bind some unfrozen configuration

        # Stores machine infos, the only *three* machine one because we have a
        # target machine info on for the user (Meson never cares about the
        # target machine.)
        machines: PerThreeMachineDefaultable[MachineInfo] = PerThreeMachineDefaultable()

        # Similar to coredata.compilers, but lower level in that there is no
        # meta data, only names/paths.
        binaries: PerMachineDefaultable[BinaryTable] = PerMachineDefaultable()

        # Misc other properties about each machine.
        properties: PerMachineDefaultable[Properties] = PerMachineDefaultable()

        # CMake toolchain variables
        cmakevars: PerMachineDefaultable[CMakeVariables] = PerMachineDefaultable()

        ## Setup build machine defaults

        # Will be fully initialized later using compilers later.
        machines.build = detect_machine_info()

        # Just uses hard-coded defaults and environment variables. Might be
        # overwritten by a native file.
        binaries.build = BinaryTable()
        properties.build = Properties()

        # Options with the key parsed into an OptionKey type.
        #
        # Note that order matters because of 'buildtype', if it is after
        # 'optimization' and 'debug' keys, it override them.
        self.options: OptionDict = collections.OrderedDict()

        # Environment variables with the name converted into an OptionKey type.
        # These have subtly different behavior compared to machine files, so do
        # not store them in self.options.  See _set_default_options_from_env.
        self.env_opts: OptionDict = {}

        self.machinestore = machinefile.MachineFileStore(self.coredata.config_files, self.coredata.cross_files, self.source_dir)

        ## Read in native file(s) to override build machine configuration

        if self.coredata.config_files is not None:
            config = machinefile.parse_machine_files(self.coredata.config_files, self.source_dir)
            binaries.build = BinaryTable(config.get('binaries', {}))
            properties.build = Properties(config.get('properties', {}))
            cmakevars.build = CMakeVariables(config.get('cmake', {}))
            self._load_machine_file_options(
                config, properties.build,
                MachineChoice.BUILD if self.coredata.cross_files else MachineChoice.HOST)

        ## Read in cross file(s) to override host machine configuration

        if self.coredata.cross_files:
            config = machinefile.parse_machine_files(self.coredata.cross_files, self.source_dir)
            properties.host = Properties(config.get('properties', {}))
            binaries.host = BinaryTable(config.get('binaries', {}))
            cmakevars.host = CMakeVariables(config.get('cmake', {}))
            if 'host_machine' in config:
                machines.host = MachineInfo.from_literal(config['host_machine'])
            if 'target_machine' in config:
                machines.target = MachineInfo.from_literal(config['target_machine'])
            # Keep only per machine options from the native file. The cross
            # file takes precedence over all other options.
            for key, value in list(self.options.items()):
                if self.coredata.optstore.is_per_machine_option(key):
                    self.options[key.as_build()] = value
            self._load_machine_file_options(config, properties.host, MachineChoice.HOST)

        ## "freeze" now initialized configuration, and "save" to the class.

        self.machines = machines.default_missing()
        self.binaries = binaries.default_missing()
        self.properties = properties.default_missing()
        self.cmakevars = cmakevars.default_missing()

        # Set host machine info for machine-aware handling of directory options
        self.coredata.optstore.set_host_machine(self.machines.host)

        # Take default value from env if not set in cross/native files or command line.
        self._set_default_options_from_env()
        self._set_default_binaries_from_env()
        self._set_default_properties_from_env()

        # Warn if the user is using two different ways of setting build-type
        # options that override each other
        bt = OptionKey('buildtype')
        db = OptionKey('debug')
        op = OptionKey('optimization')
        if bt in self.options and (db in self.options or op in self.options):
            mlog.warning('Recommend using either -Dbuildtype or -Doptimization + -Ddebug. '
                         'Using both is redundant since they override each other. '
                         'See: https://mesonbuild.com/Builtin-options.html#build-type-options',
                         fatal=False)

        # Filter out build machine options that are not valid per-project.
        # We allow this in the file because it makes the machine files more
        # useful (ie, the same file can be used for host == build configuration
        # a host != build configuration)
        self.options = {k: v for k, v in self.options.items()
                        if k.machine is MachineChoice.HOST or self.coredata.optstore.is_per_machine_option(k)}

        exe_wrapper = self.lookup_binary_entry(MachineChoice.HOST, 'exe_wrapper')
        if exe_wrapper is not None:
            self.exe_wrapper = ExternalProgram.from_bin_list(self, MachineChoice.HOST, 'exe_wrapper')
        else:
            self.exe_wrapper = None

        self.default_cmake = ['cmake']
        self.default_pkgconfig = ['pkg-config']
        self.wrap_resolver: T.Optional['Resolver'] = None

    def mfilestr2key(self, machine_file_string: str, section: T.Optional[str], section_subproject: T.Optional[str], machine: MachineChoice) -> OptionKey:
        key = OptionKey.from_string(machine_file_string)
        if key.subproject:
            suggestion = section if section == 'project options' else 'built-in options'
            raise MesonException(f'Do not set subproject options in [{section}] section, use [subproject:{suggestion}] instead.')
        if section_subproject:
            key = key.evolve(subproject=section_subproject)
        if machine == MachineChoice.BUILD:
            if key.machine == MachineChoice.BUILD:
                mlog.deprecation('Setting build machine options in the native file does not need the "build." prefix', once=True)
            return key.evolve(machine=machine)
        return key

    def _load_machine_file_options(self, config: T.Mapping[str, T.Mapping[str, ElementaryOptionValues]],
                                   properties: Properties, machine: MachineChoice) -> None:
        """Read the contents of a Machine file and put it in the options store."""

        # Look for any options in the deprecated paths section, warn about
        # those, then assign them. They will be overwritten by the ones in the
        # "built-in options" section if they're in both sections.
        paths = config.get('paths')
        if paths:
            mlog.deprecation('The [paths] section is deprecated, use the [built-in options] section instead.')
            for strk, v in paths.items():
                k = self.mfilestr2key(strk, 'paths', None, machine)
                self.options[k] = v

        # Next look for compiler options in the "properties" section, this is
        # also deprecated, and these will also be overwritten by the "built-in
        # options" section. We need to remove these from this section, as well.
        deprecated_properties: T.Set[str] = set()
        for lang in compilers.all_languages:
            deprecated_properties.add(lang + '_args')
            deprecated_properties.add(lang + '_link_args')
        for strk, v in properties.properties.copy().items():
            if strk in deprecated_properties:
                mlog.deprecation(f'{strk} in the [properties] section of the machine file is deprecated, use the [built-in options] section.')
                k = self.mfilestr2key(strk, 'properties', None, machine)
                self.options[k] = v
                del properties.properties[strk]

        for section, values in config.items():
            if ':' in section:
                section_subproject, section = section.split(':', 1)
            else:
                section_subproject = ''
            if section == 'built-in options':
                for strk, v in values.items():
                    key = self.mfilestr2key(strk, section, section_subproject, machine)
                    # If we're in the cross file, and there is a `build.foo` warn about that. Later we'll remove it.
                    if machine is MachineChoice.HOST and key.machine is not machine:
                        mlog.deprecation('Setting build machine options in cross files, please use a native file instead, this will be removed in meson 2.0', once=True)
                    self.options[key] = v
            elif section == 'project options' and machine is MachineChoice.HOST:
                # Project options are only for the host machine, we don't want
                # to read these from the native file
                for strk, v in values.items():
                    # Project options are always for the host machine
                    key = self.mfilestr2key(strk, section, section_subproject, machine)
                    self.options[key] = v
            elif ':' in section:
                correct_subproject, correct_section = section.split(':')[-2:]
                raise MesonException(
                    'Subproject options should always be set as '
                    '`[subproject:section]`, even if the options are from a '
                    'nested subproject. '
                    f'Replace `[{section_subproject}:{section}]` with `[{correct_subproject}:{correct_section}]`')

    def _set_default_options_from_env(self) -> None:
        opts: T.List[T.Tuple[str, str]] = (
            [(v, f'{k}_args') for k, v in compilers.compilers.CFLAGS_MAPPING.items()] +
            NON_LANG_ENV_OPTIONS
        )

        env_opts: T.DefaultDict[OptionKey, T.List[str]] = collections.defaultdict(list)

        for (evar, keyname), for_machine in itertools.product(opts, MachineChoice):
            p_env = _get_env_var(for_machine, self.is_cross_build(), evar)
            if p_env is not None:
                # these may contain duplicates, which must be removed, else
                # a duplicates-in-array-option warning arises.
                if keyname == 'cmake_prefix_path':
                    if self.machines[for_machine].is_windows():
                        # Cannot split on ':' on Windows because its in the drive letter
                        _p_env = p_env.split(os.pathsep)
                    else:
                        # https://github.com/mesonbuild/meson/issues/7294
                        _p_env = re.split(r':|;', p_env)
                    p_list = list(mesonlib.OrderedSet(_p_env))
                elif keyname == 'pkg_config_path':
                    p_list = list(mesonlib.OrderedSet(p_env.split(os.pathsep)))
                else:
                    p_list = split_args(p_env)
                p_list = [e for e in p_list if e]  # filter out any empty elements

                # Take env vars only on first invocation, if the env changes when
                # reconfiguring it gets ignored.
                # FIXME: We should remember if we took the value from env to warn
                # if it changes on future invocations.
                if self.first_invocation:
                    if keyname == 'ldflags':
                        for lang in compilers.compilers.LANGUAGES_USING_LDFLAGS:
                            key = OptionKey(name=f'{lang}_link_args', machine=for_machine)
                            env_opts[key].extend(p_list)
                    elif keyname == 'cppflags':
                        for lang in compilers.compilers.LANGUAGES_USING_CPPFLAGS:
                            key = OptionKey(f'{lang}_args', machine=for_machine)
                            env_opts[key].extend(p_list)
                    else:
                        key = OptionKey.from_string(keyname).evolve(machine=for_machine)
                        env_opts[key].extend(p_list)

        # If this is an environment variable, we have to
        # store it separately until the compiler is
        # instantiated, as we don't know whether the
        # compiler will want to use these arguments at link
        # time and compile time (instead of just at compile
        # time) until we're instantiating that `Compiler`
        # object. This is required so that passing
        # `-Dc_args=` on the command line and `$CFLAGS`
        # have subtly different behavior. `$CFLAGS` will be
        # added to the linker command line if the compiler
        # acts as a linker driver, `-Dc_args` will not.
        for (_, keyname), for_machine in itertools.product(NON_LANG_ENV_OPTIONS, MachineChoice):
            key = OptionKey.from_string(keyname).evolve(machine=for_machine)
            # Only store options that are not already in self.options,
            # otherwise we'd override the machine files
            if key in env_opts and key not in self.options:
                self.options[key] = env_opts[key]
                del env_opts[key]

        self.env_opts.update(env_opts)

    def _set_default_binaries_from_env(self) -> None:
        """Set default binaries from the environment.

        For example, pkg-config can be set via PKG_CONFIG, or in the machine
        file. We want to set the default to the env variable.
        """
        opts = itertools.chain(envconfig.DEPRECATED_ENV_PROG_MAP.items(),
                               envconfig.ENV_VAR_PROG_MAP.items())

        for (name, evars), for_machine in itertools.product(opts, MachineChoice):
            for evar in evars:
                p_env = _get_env_var(for_machine, self.is_cross_build(), evar)
                if p_env is not None:
                    if os.path.exists(p_env):
                        self.binaries[for_machine].binaries.setdefault(name, [p_env])
                    else:
                        self.binaries[for_machine].binaries.setdefault(name, mesonlib.split_args(p_env))
                    break

    def _set_default_properties_from_env(self) -> None:
        """Properties which can also be set from the environment."""
        # name, evar, split
        opts: T.List[T.Tuple[str, T.List[str], bool]] = [
            ('boost_includedir', ['BOOST_INCLUDEDIR'], False),
            ('boost_librarydir', ['BOOST_LIBRARYDIR'], False),
            ('boost_root', ['BOOST_ROOT', 'BOOSTROOT'], True),
            ('java_home', ['JAVA_HOME'], False),
        ]

        for (name, evars, split), for_machine in itertools.product(opts, MachineChoice):
            for evar in evars:
                p_env = _get_env_var(for_machine, self.is_cross_build(), evar)
                if p_env is not None:
                    if split:
                        self.properties[for_machine].properties.setdefault(name, p_env.split(os.pathsep))
                    else:
                        self.properties[for_machine].properties.setdefault(name, p_env)
                    break

    def create_new_coredata(self, options: cmdline.SharedCMDOptions) -> None:
        # WARNING: Don't use any values from coredata in __init__. It gets
        # re-initialized with project options by the interpreter during
        # build file parsing.
        # meson_command is used by the regenchecker script, which runs meson
        meson_command = mesonlib.get_meson_command()
        if meson_command is None:
            meson_command = []
        else:
            meson_command = meson_command.copy()
        self.coredata = coredata.CoreData(options, self.scratch_dir, meson_command)
        self.first_invocation = True

    def init_backend_options(self, backend_name: str) -> None:
        # Only init backend options on first invocation otherwise it would
        # override values previously set from command line.
        if not self.first_invocation:
            return

        self.coredata.init_backend_options(backend_name)
        for k, v in self.options.items():
            if self.coredata.optstore.is_backend_option(k):
                self.coredata.optstore.set_option(k, v)

    def is_cross_build(self, when_building_for: MachineChoice = MachineChoice.HOST) -> bool:
        return self.coredata.is_cross_build(when_building_for)

    def dump_coredata(self) -> str:
        return coredata.save(self.coredata, self.get_build_dir())

    def get_log_dir(self) -> str:
        return self.log_dir

    def get_coredata(self) -> coredata.CoreData:
        return self.coredata

    @staticmethod
    def get_build_command(unbuffered: bool = False) -> T.List[str]:
        cmd = mesonlib.get_meson_command()
        if cmd is None:
            raise MesonBugException('No command?')
        cmd = cmd.copy()
        if unbuffered and 'python' in os.path.basename(cmd[0]):
            cmd.insert(1, '-u')
        return cmd

    def lookup_binary_entry(self, for_machine: MachineChoice, name: str) -> T.Optional[T.List[str]]:
        return self.binaries[for_machine].lookup_entry(name)

    def get_scratch_dir(self) -> str:
        return self.scratch_dir

    def get_source_dir(self) -> str:
        return self.source_dir

    def get_build_dir(self) -> str:
        return self.build_dir

    def get_import_lib_dir(self) -> str:
        "Install dir for the import library (library used for linking)"
        return self.get_libdir()

    def get_shared_module_dir(self) -> str:
        "Install dir for shared modules that are loaded at runtime"
        return self.get_libdir()

    def get_shared_lib_dir(self) -> str:
        "Install dir for the shared library"
        m = self.machines.host
        # Windows has no RPATH or similar, so DLLs must be next to EXEs.
        if m.is_windows() or m.is_cygwin():
            return self.get_bindir()
        return self.get_libdir()

    def get_jar_dir(self) -> str:
        """Install dir for JAR files"""
        return f"{self.get_datadir()}/java"

    def get_static_lib_dir(self) -> str:
        "Install dir for the static library"
        return self.get_libdir()

    def get_prefix(self) -> str:
        return _as_str(self.coredata.optstore.get_value_for(OptionKey('prefix')))

    def get_libdir(self) -> str:
        return _as_str(self.coredata.optstore.get_value_for(OptionKey('libdir')))

    def get_libexecdir(self) -> str:
        return _as_str(self.coredata.optstore.get_value_for(OptionKey('libexecdir')))

    def get_bindir(self) -> str:
        return _as_str(self.coredata.optstore.get_value_for(OptionKey('bindir')))

    def get_includedir(self) -> str:
        return _as_str(self.coredata.optstore.get_value_for(OptionKey('includedir')))

    def get_mandir(self) -> str:
        return _as_str(self.coredata.optstore.get_value_for(OptionKey('mandir')))

    def get_datadir(self) -> str:
        return _as_str(self.coredata.optstore.get_value_for(OptionKey('datadir')))

    def get_compiler_system_lib_dirs(self, for_machine: MachineChoice) -> T.List[str]:
        for comp in self.coredata.compilers[for_machine].values():
            if comp.id == 'clang':
                index = 1
                break
            elif comp.id == 'gcc':
                index = 2
                break
        else:
            # This option is only supported by gcc and clang. If we don't get a
            # GCC or Clang compiler return and empty list.
            return []

        p, out, _ = Popen_safe(comp.get_exelist() + ['-print-search-dirs'])
        if p.returncode != 0:
            raise mesonlib.MesonException('Could not calculate system search dirs')
        split = out.split('\n')[index].lstrip('libraries: =').split(':')
        return [os.path.normpath(p) for p in split]

    def get_compiler_system_include_dirs(self, for_machine: MachineChoice) -> T.List[str]:
        for comp in self.coredata.compilers[for_machine].values():
            if comp.id == 'clang':
                break
            elif comp.id == 'gcc':
                break
        else:
            # This option is only supported by gcc and clang. If we don't get a
            # GCC or Clang compiler return and empty list.
            return []
        return comp.get_default_include_dirs()

    def need_exe_wrapper(self, for_machine: MachineChoice = MachineChoice.HOST) -> bool:
        value = self.properties[for_machine].get('needs_exe_wrapper', None)
        if value is not None:
            assert isinstance(value, bool), 'for mypy'
            return value
        if not self.is_cross_build():
            return False
        return not machine_info_can_run(self.machines[for_machine])

    def get_exe_wrapper(self) -> T.Optional[ExternalProgram]:
        if not self.need_exe_wrapper():
            return None
        return self.exe_wrapper

    def has_exe_wrapper(self) -> bool:
        return self.exe_wrapper is not None and self.exe_wrapper.found()

    def get_env_for_paths(self, library_paths: T.Set[str], extra_paths: T.Set[str]) -> mesonlib.EnvironmentVariables:
        env = mesonlib.EnvironmentVariables()
        need_wine = not self.machines.build.is_windows() and self.machines.host.is_windows()
        if need_wine:
            # Executable paths should be in both PATH and WINEPATH.
            # - Having them in PATH makes bash completion find it,
            #   and make running "foo.exe" find it when wine-binfmt is installed.
            # - Having them in WINEPATH makes "wine foo.exe" find it.
            library_paths.update(extra_paths)
        if library_paths:
            if need_wine:
                env.prepend('WINEPATH', list(library_paths), separator=';')
            elif self.machines.host.is_windows() or self.machines.host.is_cygwin():
                extra_paths.update(library_paths)
            elif self.machines.host.is_darwin():
                env.prepend('DYLD_LIBRARY_PATH', list(library_paths))
            else:
                env.prepend('LD_LIBRARY_PATH', list(library_paths))
        if extra_paths:
            env.prepend('PATH', list(extra_paths))
        return env

    def add_lang_args(self, lang: Language, comp: T.Type['Compiler'],
                      for_machine: MachineChoice) -> None:
        """Add global language arguments that are needed before compiler/linker detection."""
        description = f'Extra arguments passed to the {lang}'
        argkey = OptionKey(f'{lang}_args', machine=for_machine)
        largkey = OptionKey(f'{lang}_link_args', machine=for_machine)

        comp_args_from_envvar = False
        comp_options = self.coredata.optstore.get_pending_value(argkey)
        if comp_options is None:
            comp_args_from_envvar = True
            comp_options = self.env_opts.get(argkey, [])

        link_options = self.coredata.optstore.get_pending_value(largkey)
        if link_options is None:
            link_options = self.env_opts.get(largkey, [])

        assert isinstance(comp_options, (str, list)), 'for mypy'
        assert isinstance(link_options, (str, list)), 'for mypy'

        cargs = options.UserStringArrayOption(
            argkey.name,
            description + ' compiler',
            comp_options, split_args=True, allow_dups=True)

        largs = options.UserStringArrayOption(
            largkey.name,
            description + ' linker',
            link_options, split_args=True, allow_dups=True)

        self.coredata.optstore.add_compiler_option(lang, argkey, cargs)
        self.coredata.optstore.add_compiler_option(lang, largkey, largs)

        if comp.USED_FOR_SEPARATE_LINKING_STEP and comp_args_from_envvar:
            # If the compiler acts as a linker driver, and we're using the
            # environment variable flags for both the compiler and linker
            # arguments, then put the compiler flags in the linker flags as well.
            # This is how autotools works, and the env vars feature is for
            # autotools compatibility.
            largs.extend_value(comp_options)

    def update_build_machine(self, compilers: T.Optional[CompilerDict] = None) -> None:
        """Redetect the build machine and update the machine definitions

        :compilers: An optional dictionary of compilers to use instead of the coredata dict.
        """
        compilers = compilers or self.coredata.compilers.build

        machines = self.machines.miss_defaulting()
        machines.build = detect_machine_info(compilers)
        self.machines = machines.default_missing()
