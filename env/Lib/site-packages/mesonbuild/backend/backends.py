# SPDX-License-Identifier: Apache-2.0
# Copyright 2012-2016 The Meson development team

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass, InitVar
from functools import lru_cache
from itertools import chain
from pathlib import Path
import copy
import enum
import json
import os
import pickle
import re
import shlex
import shutil
import typing as T
import hashlib

from .. import build
from .. import dependencies
from .. import programs
from .. import mesonlib
from .. import mlog
from ..compilers import LANGUAGES_USING_LDFLAGS, detect
from ..mesonlib import (
    File, MachineChoice, MesonException, OrderedSet,
    ExecutableSerialisation, classify_unity_sources,
)
from ..options import OptionKey


if T.TYPE_CHECKING:
    from .._typing import ImmutableListProtocol
    from ..arglist import CompilerArgs
    from ..compilers import Compiler
    from ..environment import Environment
    from ..interpreter import Interpreter, Test
    from ..linkers.linkers import StaticLinker
    from ..mesonlib import FileMode, FileOrString

    from typing_extensions import TypedDict

    _ALL_SOURCES_TYPE = T.List[T.Union[File, build.CustomTarget, build.CustomTargetIndex, build.GeneratedList]]

    class TargetIntrospectionData(TypedDict):

        language: str
        compiler: T.List[str]
        parameters: T.List[str]
        sources: T.List[str]
        generated_sources: T.List[str]


# Languages that can mix with C or C++ but don't support unity builds yet
# because the syntax we use for unity builds is specific to C/++/ObjC/++.
# Assembly files cannot be unitified and neither can LLVM IR files
LANGS_CANT_UNITY = ('d', 'fortran', 'vala')

@dataclass(eq=False)
class RegenInfo:
    source_dir: str
    build_dir: str
    depfiles: T.List[str]

class TestProtocol(enum.Enum):

    EXITCODE = 0
    TAP = 1
    GTEST = 2
    RUST = 3

    @classmethod
    def from_str(cls, string: str) -> 'TestProtocol':
        if string == 'exitcode':
            return cls.EXITCODE
        elif string == 'tap':
            return cls.TAP
        elif string == 'gtest':
            return cls.GTEST
        elif string == 'rust':
            return cls.RUST
        raise MesonException(f'unknown test format {string}')

    def __str__(self) -> str:
        cls = type(self)
        if self is cls.EXITCODE:
            return 'exitcode'
        elif self is cls.GTEST:
            return 'gtest'
        elif self is cls.RUST:
            return 'rust'
        return 'tap'


@dataclass(eq=False)
class CleanTrees:
    '''
    Directories outputted by custom targets that have to be manually cleaned
    because on Linux `ninja clean` only deletes empty directories.
    '''
    build_dir: str
    trees: T.List[str]

@dataclass(eq=False)
class InstallData:
    source_dir: str
    build_dir: str
    prefix: str
    libdir: str
    strip_bin: T.List[str]
    # TODO: in python 3.8 or with typing_Extensions this could be:
    # `T.Union[T.Literal['preserve'], int]`, which would be more accurate.
    install_umask: T.Union[str, int]
    mesonintrospect: T.List[str]
    version: str

    def __post_init__(self) -> None:
        self.targets: T.List[TargetInstallData] = []
        self.headers: T.List[InstallDataBase] = []
        self.man: T.List[InstallDataBase] = []
        self.emptydir: T.List[InstallEmptyDir] = []
        self.data: T.List[InstallDataBase] = []
        self.symlinks: T.List[InstallSymlinkData] = []
        self.install_scripts: T.List[ExecutableSerialisation] = []
        self.install_subdirs: T.List[SubdirInstallData] = []

@dataclass(eq=False)
class TargetInstallData:
    fname: str
    outdir: str
    outdir_name: InitVar[T.Optional[str]]
    strip: bool
    install_name_mappings: T.Mapping[str, str]
    rpath_dirs_to_remove: T.Set[bytes]
    install_rpath: str
    # TODO: install_mode should just always be a FileMode object
    install_mode: T.Optional['FileMode']
    subproject: str
    optional: bool = False
    tag: T.Optional[str] = None
    can_strip: bool = False

    def __post_init__(self, outdir_name: T.Optional[str]) -> None:
        if outdir_name is None:
            outdir_name = os.path.join('{prefix}', self.outdir)
        self.out_name = os.path.join(outdir_name, os.path.basename(self.fname))

@dataclass(eq=False)
class InstallEmptyDir:
    path: str
    install_mode: 'FileMode'
    subproject: str
    tag: T.Optional[str] = None

@dataclass(eq=False)
class InstallDataBase:
    path: str
    install_path: str
    install_path_name: str
    install_mode: 'FileMode'
    subproject: str
    tag: T.Optional[str] = None
    data_type: T.Optional[str] = None
    follow_symlinks: T.Optional[bool] = None

@dataclass(eq=False)
class InstallSymlinkData:
    target: str
    name: str
    install_path: str
    subproject: str
    tag: T.Optional[str] = None
    allow_missing: bool = False

# cannot use dataclass here because "exclude" is out of order
class SubdirInstallData(InstallDataBase):
    def __init__(self, path: str, install_path: str, install_path_name: str,
                 install_mode: 'FileMode', exclude: T.Tuple[T.Set[str], T.Set[str]],
                 subproject: str, tag: T.Optional[str] = None, data_type: T.Optional[str] = None,
                 follow_symlinks: T.Optional[bool] = None):
        super().__init__(path, install_path, install_path_name, install_mode, subproject, tag, data_type, follow_symlinks)
        self.exclude = exclude


@dataclass(eq=False)
class TestSerialisation:
    name: str
    project_name: str
    suite: T.List[str]
    fname: T.List[str]
    is_cross_built: bool
    exe_wrapper: T.Optional[programs.ExternalProgram]
    needs_exe_wrapper: bool
    is_parallel: bool
    cmd_args: T.List[str]
    env: mesonlib.EnvironmentVariables
    should_fail: bool
    timeout: T.Optional[int]
    workdir: T.Optional[str]
    extra_paths: T.List[str]
    protocol: TestProtocol
    priority: int
    cmd_is_built: bool
    cmd_is_exe: bool
    depends: T.List[str]
    version: str
    verbose: bool

    def __post_init__(self) -> None:
        if self.exe_wrapper is not None:
            assert isinstance(self.exe_wrapper, programs.ExternalProgram)


def get_backend_from_name(backend: str, build: T.Optional[build.Build] = None, interpreter: T.Optional['Interpreter'] = None) -> T.Optional['Backend']:
    if backend == 'ninja':
        from . import ninjabackend
        return ninjabackend.NinjaBackend(build, interpreter)
    elif backend == 'vs':
        from . import vs2010backend
        return vs2010backend.autodetect_vs_version(build, interpreter)
    elif backend == 'vs2010':
        from . import vs2010backend
        return vs2010backend.Vs2010Backend(build, interpreter)
    elif backend == 'vs2012':
        from . import vs2012backend
        return vs2012backend.Vs2012Backend(build, interpreter)
    elif backend == 'vs2013':
        from . import vs2013backend
        return vs2013backend.Vs2013Backend(build, interpreter)
    elif backend == 'vs2015':
        from . import vs2015backend
        return vs2015backend.Vs2015Backend(build, interpreter)
    elif backend == 'vs2017':
        from . import vs2017backend
        return vs2017backend.Vs2017Backend(build, interpreter)
    elif backend == 'vs2019':
        from . import vs2019backend
        return vs2019backend.Vs2019Backend(build, interpreter)
    elif backend == 'vs2022':
        from . import vs2022backend
        return vs2022backend.Vs2022Backend(build, interpreter)
    elif backend == 'xcode':
        from . import xcodebackend
        return xcodebackend.XCodeBackend(build, interpreter)
    elif backend == 'none':
        from . import nonebackend
        return nonebackend.NoneBackend(build, interpreter)
    return None


def get_genvslite_backend(genvsname: str, build: T.Optional[build.Build] = None, interpreter: T.Optional['Interpreter'] = None) -> T.Optional['Backend']:
    if genvsname == 'vs2022':
        from . import vs2022backend
        return vs2022backend.Vs2022Backend(build, interpreter, gen_lite = True)
    return None

# This class contains the basic functionality that is needed by all backends.
# Feel free to move stuff in and out of it as you see fit.
class Backend:

    environment: T.Optional['Environment']
    name = '<UNKNOWN>'

    def __init__(self, build: T.Optional[build.Build], interpreter: T.Optional['Interpreter']):
        # Make it possible to construct a dummy backend
        # This is used for introspection without a build directory
        if build is None:
            self.environment = None
            return
        self.build = build
        self.interpreter = interpreter
        self.environment = build.environment
        self.processed_targets: T.Set[str] = set()
        self.build_dir = self.environment.get_build_dir()
        self.source_dir = self.environment.get_source_dir()
        self.build_to_src = mesonlib.relpath(self.environment.get_source_dir(),
                                             self.environment.get_build_dir())
        self.src_to_build = mesonlib.relpath(self.environment.get_build_dir(),
                                             self.environment.get_source_dir())

    # If requested via 'capture = True', returns captured compile args per
    # target (e.g. captured_args[target]) that can be used later, for example,
    # to populate things like intellisense fields in generated visual studio
    # projects (as is the case when using '--genvslite').
    #
    # 'vslite_ctx' is only provided when
    # we expect this backend setup/generation to make use of previously captured
    # compile args (as is the case when using '--genvslite').
    def generate(self, capture: bool = False, vslite_ctx: T.Optional[T.Dict] = None) -> T.Optional[T.Dict]:
        raise RuntimeError(f'generate is not implemented in {type(self).__name__}')

    def get_target_filename(self, t: T.Union[build.Target, build.CustomTargetIndex], *, warn_multi_output: bool = True) -> str:
        if isinstance(t, build.CustomTarget):
            if warn_multi_output and len(t.get_outputs()) != 1:
                mlog.warning(f'custom_target {t.name!r} has more than one output! '
                             f'Using the first one. Consider using `{t.name}[0]`.')
            filename = t.get_outputs()[0]
        elif isinstance(t, build.CustomTargetIndex):
            filename = t.get_outputs()[0]
        else:
            assert isinstance(t, build.BuildTarget), t
            filename = t.get_filename()
        return os.path.join(self.get_target_dir(t), filename)

    def get_target_filename_abs(self, target: T.Union[build.Target, build.CustomTargetIndex]) -> str:
        return os.path.join(self.environment.get_build_dir(), self.get_target_filename(target))

    def get_target_debug_filename(self, target: build.BuildTarget) -> T.Optional[str]:
        assert isinstance(target, build.BuildTarget), target
        if target.get_debug_filename():
            debug_filename = target.get_debug_filename()
            return os.path.join(self.get_target_dir(target), debug_filename)
        else:
            return None

    def get_target_debug_filename_abs(self, target: build.BuildTarget) -> T.Optional[str]:
        assert isinstance(target, build.BuildTarget), target
        if not target.get_debug_filename():
            return None
        return os.path.join(self.environment.get_build_dir(), self.get_target_debug_filename(target))

    def get_source_dir_include_args(self, target: build.BuildTarget, compiler: 'Compiler', *, absolute_path: bool = False) -> T.List[str]:
        curdir = target.get_subdir()
        if absolute_path:
            lead = self.source_dir
        else:
            lead = self.build_to_src
        tmppath = os.path.normpath(os.path.join(lead, curdir))
        return compiler.get_include_args(tmppath, False)

    def get_build_dir_include_args(self, target: build.BuildTarget, compiler: 'Compiler', *, absolute_path: bool = False) -> T.List[str]:
        if absolute_path:
            curdir = os.path.join(self.build_dir, target.get_subdir())
        else:
            curdir = target.get_subdir()
            if curdir == '':
                curdir = '.'
        return compiler.get_include_args(curdir, False)

    def get_target_filename_for_linking(self, target: T.Union[build.Target, build.CustomTargetIndex]) -> T.Optional[str]:
        # On some platforms (msvc for instance), the file that is used for
        # dynamic linking is not the same as the dynamic library itself. This
        # file is called an import library, and we want to link against that.
        # On all other platforms, we link to the library directly.
        if isinstance(target, build.SharedLibrary):
            link_lib = target.get_import_filename() or target.get_filename()
            # In AIX, if we archive .so, the blibpath must link to archived shared library otherwise to the .so file.
            if mesonlib.is_aix() and target.aix_so_archive:
                link_lib = re.sub('[.][a]([.]?([0-9]+))*([.]?([a-z]+))*', '.a', link_lib.replace('.so', '.a'))
            return Path(self.get_target_dir(target), link_lib).as_posix()
        elif isinstance(target, build.StaticLibrary):
            return Path(self.get_target_dir(target), target.get_filename()).as_posix()
        elif isinstance(target, (build.CustomTarget, build.CustomTargetIndex)):
            if not target.is_linkable_target():
                raise MesonException(f'Tried to link against custom target "{target.name}", which is not linkable.')
            return Path(self.get_target_dir(target), target.get_filename()).as_posix()
        elif isinstance(target, build.Executable):
            if target.import_filename:
                return Path(self.get_target_dir(target), target.get_import_filename()).as_posix()
            else:
                return None
        raise AssertionError(f'BUG: Tried to link to {target!r} which is not linkable')

    @lru_cache(maxsize=None)
    def get_target_dir(self, target: T.Union[build.Target, build.CustomTargetIndex]) -> str:
        if isinstance(target, build.RunTarget):
            # this produces no output, only a dummy top-level name
            dirname = ''
        elif self.environment.coredata.get_option(OptionKey('layout')) == 'mirror':
            dirname = target.get_subdir()
        else:
            dirname = 'meson-out'
        return dirname

    def get_target_dir_relative_to(self, t: build.Target, o: build.Target) -> str:
        '''Get a target dir relative to another target's directory'''
        target_dir = os.path.join(self.environment.get_build_dir(), self.get_target_dir(t))
        othert_dir = os.path.join(self.environment.get_build_dir(), self.get_target_dir(o))
        return os.path.relpath(target_dir, othert_dir)

    def get_target_source_dir(self, target: build.Target) -> str:
        # if target dir is empty, avoid extraneous trailing / from os.path.join()
        target_dir = self.get_target_dir(target)
        if target_dir:
            return os.path.join(self.build_to_src, target_dir)
        return self.build_to_src

    def get_target_private_dir(self, target: T.Union[build.BuildTarget, build.CustomTarget, build.CustomTargetIndex]) -> str:
        return os.path.join(self.get_target_filename(target, warn_multi_output=False) + '.p')

    def get_target_private_dir_abs(self, target: T.Union[build.BuildTarget, build.CustomTarget, build.CustomTargetIndex]) -> str:
        return os.path.join(self.environment.get_build_dir(), self.get_target_private_dir(target))

    @lru_cache(maxsize=None)
    def get_target_generated_dir(
            self, target: T.Union[build.BuildTarget, build.CustomTarget, build.CustomTargetIndex],
            gensrc: T.Union[build.CustomTarget, build.CustomTargetIndex, build.GeneratedList],
            src: str) -> str:
        """
        Takes a BuildTarget, a generator source (CustomTarget or GeneratedList),
        and a generated source filename.
        Returns the full path of the generated source relative to the build root
        """
        # CustomTarget generators output to the build dir of the CustomTarget
        if isinstance(gensrc, (build.CustomTarget, build.CustomTargetIndex)):
            return os.path.join(self.get_target_dir(gensrc), src)
        # GeneratedList generators output to the private build directory of the
        # target that the GeneratedList is used in
        return os.path.join(self.get_target_private_dir(target), src)

    def get_unity_source_file(self, target: T.Union[build.BuildTarget, build.CustomTarget, build.CustomTargetIndex],
                              suffix: str, number: int) -> mesonlib.File:
        # There is a potential conflict here, but it is unlikely that
        # anyone both enables unity builds and has a file called foo-unity.cpp.
        osrc = f'{target.name}-unity{number}.{suffix}'
        return mesonlib.File.from_built_file(self.get_target_private_dir(target), osrc)

    def generate_unity_files(self, target: build.BuildTarget, unity_src: str) -> T.List[mesonlib.File]:
        abs_files: T.List[str] = []
        result: T.List[mesonlib.File] = []
        compsrcs = classify_unity_sources(target.compilers.values(), unity_src)
        unity_size = target.get_option(OptionKey('unity_size'))
        assert isinstance(unity_size, int), 'for mypy'

        def init_language_file(suffix: str, unity_file_number: int) -> T.TextIO:
            unity_src = self.get_unity_source_file(target, suffix, unity_file_number)
            outfileabs = unity_src.absolute_path(self.environment.get_source_dir(),
                                                 self.environment.get_build_dir())
            outfileabs_tmp = outfileabs + '.tmp'
            abs_files.append(outfileabs)
            outfileabs_tmp_dir = os.path.dirname(outfileabs_tmp)
            if not os.path.exists(outfileabs_tmp_dir):
                os.makedirs(outfileabs_tmp_dir)
            result.append(unity_src)
            return open(outfileabs_tmp, 'w', encoding='utf-8')

        # For each language, generate unity source files and return the list
        for comp, srcs in compsrcs.items():
            files_in_current = unity_size + 1
            unity_file_number = 0
            # TODO: this could be simplified with an algorithm that pre-sorts
            # the sources into the size of chunks we want
            ofile = None
            for src in srcs:
                if files_in_current >= unity_size:
                    if ofile:
                        ofile.close()
                    ofile = init_language_file(comp.get_default_suffix(), unity_file_number)
                    unity_file_number += 1
                    files_in_current = 0
                ofile.write(f'#include<{src}>\n')
                files_in_current += 1
            if ofile:
                ofile.close()

        for x in abs_files:
            mesonlib.replace_if_different(x, x + '.tmp')
        return result

    @staticmethod
    def relpath(todir: str, fromdir: str) -> str:
        return os.path.relpath(os.path.join('dummyprefixdir', todir),
                               os.path.join('dummyprefixdir', fromdir))

    def flatten_object_list(self, target: build.BuildTarget, proj_dir_to_build_root: str = ''
                            ) -> T.Tuple[T.List[str], T.List[build.BuildTargetTypes]]:
        obj_list, deps = self._flatten_object_list(target, target.get_objects(), proj_dir_to_build_root)
        return list(dict.fromkeys(obj_list)), deps

    def determine_ext_objs(self, objects: build.ExtractedObjects, proj_dir_to_build_root: str = '') -> T.List[str]:
        obj_list, _ = self._flatten_object_list(objects.target, [objects], proj_dir_to_build_root)
        return list(dict.fromkeys(obj_list))

    def _flatten_object_list(self, target: build.BuildTarget,
                             objects: T.Sequence[T.Union[str, 'File', build.ExtractedObjects]],
                             proj_dir_to_build_root: str) -> T.Tuple[T.List[str], T.List[build.BuildTargetTypes]]:
        obj_list: T.List[str] = []
        deps: T.List[build.BuildTargetTypes] = []
        for obj in objects:
            if isinstance(obj, str):
                o = os.path.join(proj_dir_to_build_root,
                                 self.build_to_src, target.get_subdir(), obj)
                obj_list.append(o)
            elif isinstance(obj, mesonlib.File):
                if obj.is_built:
                    o = os.path.join(proj_dir_to_build_root,
                                     obj.rel_to_builddir(self.build_to_src))
                    obj_list.append(o)
                else:
                    o = os.path.join(proj_dir_to_build_root,
                                     self.build_to_src)
                    obj_list.append(obj.rel_to_builddir(o))
            elif isinstance(obj, build.ExtractedObjects):
                if obj.recursive:
                    objs, d = self._flatten_object_list(obj.target, obj.objlist, proj_dir_to_build_root)
                    obj_list.extend(objs)
                    deps.extend(d)
                obj_list.extend(self._determine_ext_objs(obj, proj_dir_to_build_root))
                deps.append(obj.target)
            else:
                raise MesonException('Unknown data type in object list.')
        return obj_list, deps

    @staticmethod
    def is_swift_target(target: build.BuildTarget) -> bool:
        for s in target.sources:
            if s.endswith('swift'):
                return True
        return False

    def determine_swift_dep_dirs(self, target: build.BuildTarget) -> T.List[str]:
        result: T.List[str] = []
        for l in target.link_targets:
            result.append(self.get_target_private_dir_abs(l))
        return result

    def get_executable_serialisation(
            self, cmd: T.Sequence[T.Union[programs.ExternalProgram, build.BuildTarget, build.CustomTarget, File, str]],
            workdir: T.Optional[str] = None,
            extra_bdeps: T.Optional[T.List[build.BuildTarget]] = None,
            capture: T.Optional[str] = None,
            feed: T.Optional[str] = None,
            env: T.Optional[mesonlib.EnvironmentVariables] = None,
            tag: T.Optional[str] = None,
            verbose: bool = False,
            installdir_map: T.Optional[T.Dict[str, str]] = None) -> 'ExecutableSerialisation':

        # XXX: cmd_args either need to be lowered to strings, or need to be checked for non-string arguments, right?
        exe, *raw_cmd_args = cmd
        if isinstance(exe, programs.ExternalProgram):
            exe_cmd = exe.get_command()
            exe_for_machine = exe.for_machine
        elif isinstance(exe, build.BuildTarget):
            exe_cmd = [self.get_target_filename_abs(exe)]
            exe_for_machine = exe.for_machine
        elif isinstance(exe, build.CustomTarget):
            # The output of a custom target can either be directly runnable
            # or not, that is, a script, a native binary or a cross compiled
            # binary when exe wrapper is available and when it is not.
            # This implementation is not exhaustive but it works in the
            # common cases.
            exe_cmd = [self.get_target_filename_abs(exe)]
            exe_for_machine = MachineChoice.BUILD
        elif isinstance(exe, mesonlib.File):
            exe_cmd = [exe.rel_to_builddir(self.environment.source_dir)]
            exe_for_machine = MachineChoice.BUILD
        else:
            exe_cmd = [exe]
            exe_for_machine = MachineChoice.BUILD

        cmd_args: T.List[str] = []
        for c in raw_cmd_args:
            if isinstance(c, programs.ExternalProgram):
                p = c.get_path()
                assert isinstance(p, str)
                cmd_args.append(p)
            elif isinstance(c, (build.BuildTarget, build.CustomTarget)):
                cmd_args.append(self.get_target_filename_abs(c))
            elif isinstance(c, mesonlib.File):
                cmd_args.append(c.rel_to_builddir(self.environment.source_dir))
            else:
                cmd_args.append(c)

        machine = self.environment.machines[exe_for_machine]
        if machine.is_windows() or machine.is_cygwin():
            extra_paths = self.determine_windows_extra_paths(exe, extra_bdeps or [])
        else:
            extra_paths = []

        is_cross_built = not self.environment.machines.matches_build_machine(exe_for_machine)
        if is_cross_built and self.environment.need_exe_wrapper():
            if not self.environment.has_exe_wrapper():
                msg = 'An exe_wrapper is needed but was not found. Please define one ' \
                      'in cross file and check the command and/or add it to PATH.'
                raise MesonException(msg)
            exe_wrapper = self.environment.get_exe_wrapper()
        else:
            if exe_cmd[0].endswith('.jar'):
                exe_cmd = ['java', '-jar'] + exe_cmd
            elif exe_cmd[0].endswith('.exe') and not (mesonlib.is_windows() or mesonlib.is_cygwin() or mesonlib.is_wsl()):
                exe_cmd = ['mono'] + exe_cmd
            exe_wrapper = None

        workdir = workdir or self.environment.get_build_dir()
        return ExecutableSerialisation(exe_cmd + cmd_args, env,
                                       exe_wrapper, workdir,
                                       extra_paths, capture, feed, tag, verbose, installdir_map)

    def as_meson_exe_cmdline(self, exe: T.Union[str, mesonlib.File, build.BuildTarget, build.CustomTarget, programs.ExternalProgram],
                             cmd_args: T.Sequence[T.Union[str, mesonlib.File, build.BuildTarget, build.CustomTarget, programs.ExternalProgram]],
                             workdir: T.Optional[str] = None,
                             extra_bdeps: T.Optional[T.List[build.BuildTarget]] = None,
                             capture: T.Optional[str] = None,
                             feed: T.Optional[str] = None,
                             force_serialize: bool = False,
                             env: T.Optional[mesonlib.EnvironmentVariables] = None,
                             verbose: bool = False) -> T.Tuple[T.List[str], str]:
        '''
        Serialize an executable for running with a generator or a custom target
        '''
        cmd: T.List[T.Union[str, mesonlib.File, build.BuildTarget, build.CustomTarget, programs.ExternalProgram]] = []
        cmd.append(exe)
        cmd.extend(cmd_args)
        es = self.get_executable_serialisation(cmd, workdir, extra_bdeps, capture, feed, env, verbose=verbose)
        reasons: T.List[str] = []
        if es.extra_paths:
            reasons.append('to set PATH')

        if es.exe_wrapper:
            reasons.append('to use exe_wrapper')

        if workdir:
            reasons.append('to set workdir')

        if any('\n' in c for c in es.cmd_args):
            reasons.append('because command contains newlines')

        if env and env.varnames:
            reasons.append('to set env')

        # force_serialize passed to this function means that the VS backend has
        # decided it absolutely cannot use real commands. This is "always",
        # because it's not clear what will work (other than compilers) and so
        # we don't bother to handle a variety of common cases that probably do
        # work.
        #
        # It's also overridden for a few conditions that can't be handled
        # inside a command line

        can_use_env = env.can_use_env and not force_serialize
        force_serialize = force_serialize or bool(reasons)

        if capture:
            reasons.append('to capture output')
        if feed:
            reasons.append('to feed input')

        if can_use_env and reasons == ['to set env'] and shutil.which('env'):
            envlist = []
            for k, v in env.get_env({}).items():
                envlist.append(f'{k}={v}')
            return ['env'] + envlist + es.cmd_args, ', '.join(reasons)

        if not force_serialize:
            if not capture and not feed:
                return es.cmd_args, ''
            args: T.List[str] = []
            if capture:
                args += ['--capture', capture]
            if feed:
                args += ['--feed', feed]

            return (
                self.environment.get_build_command() + ['--internal', 'exe'] + args + ['--'] + es.cmd_args,
                ', '.join(reasons)
            )

        if isinstance(exe, (programs.ExternalProgram,
                            build.BuildTarget, build.CustomTarget)):
            basename = os.path.basename(exe.name)
        elif isinstance(exe, mesonlib.File):
            basename = os.path.basename(exe.fname)
        else:
            basename = os.path.basename(exe)

        # Can't just use exe.name here; it will likely be run more than once
        # Take a digest of the cmd args, env, workdir, capture, and feed. This
        # avoids collisions and also makes the name deterministic over
        # regenerations which avoids a rebuild by Ninja because the cmdline
        # stays the same.
        hasher = hashlib.sha1()
        if es.env:
            es.env.hash(hasher)
        hasher.update(bytes(str(es.cmd_args), encoding='utf-8'))
        hasher.update(bytes(str(es.workdir), encoding='utf-8'))
        hasher.update(bytes(str(capture), encoding='utf-8'))
        hasher.update(bytes(str(feed), encoding='utf-8'))
        digest = hasher.hexdigest()
        scratch_file = f'meson_exe_{basename}_{digest}.dat'
        exe_data = os.path.join(self.environment.get_scratch_dir(), scratch_file)
        with open(exe_data, 'wb') as f:
            pickle.dump(es, f)
        return (self.environment.get_build_command() + ['--internal', 'exe', '--unpickle', exe_data],
                ', '.join(reasons))

    def serialize_tests(self) -> T.Tuple[str, str]:
        test_data = os.path.join(self.environment.get_scratch_dir(), 'meson_test_setup.dat')
        with open(test_data, 'wb') as datafile:
            self.write_test_file(datafile)
        benchmark_data = os.path.join(self.environment.get_scratch_dir(), 'meson_benchmark_setup.dat')
        with open(benchmark_data, 'wb') as datafile:
            self.write_benchmark_file(datafile)
        return test_data, benchmark_data

    def determine_linker_and_stdlib_args(self, target: build.BuildTarget) -> T.Tuple[T.Union['Compiler', 'StaticLinker'], T.List[str]]:
        '''
        If we're building a static library, there is only one static linker.
        Otherwise, we query the target for the dynamic linker.
        '''
        if isinstance(target, build.StaticLibrary):
            return self.build.static_linker[target.for_machine], []
        l, stdlib_args = target.get_clink_dynamic_linker_and_stdlibs()
        return l, stdlib_args

    @staticmethod
    def _libdir_is_system(libdir: str, compilers: T.Mapping[str, 'Compiler'], env: 'Environment') -> bool:
        libdir = os.path.normpath(libdir)
        for cc in compilers.values():
            if libdir in cc.get_library_dirs(env):
                return True
        return False

    def get_external_rpath_dirs(self, target: build.BuildTarget) -> T.Set[str]:
        args: T.List[str] = []
        for lang in LANGUAGES_USING_LDFLAGS:
            try:
                e = self.environment.coredata.get_external_link_args(target.for_machine, lang)
                if isinstance(e, str):
                    args.append(e)
                else:
                    args.extend(e)
            except Exception:
                pass
        return self.get_rpath_dirs_from_link_args(args)

    @staticmethod
    def get_rpath_dirs_from_link_args(args: T.List[str]) -> T.Set[str]:
        dirs: T.Set[str] = set()
        # Match rpath formats:
        # -Wl,-rpath=
        # -Wl,-rpath,
        rpath_regex = re.compile(r'-Wl,-rpath[=,]([^,]+)')
        # Match solaris style compat runpath formats:
        # -Wl,-R
        # -Wl,-R,
        runpath_regex = re.compile(r'-Wl,-R[,]?([^,]+)')
        # Match symbols formats:
        # -Wl,--just-symbols=
        # -Wl,--just-symbols,
        symbols_regex = re.compile(r'-Wl,--just-symbols[=,]([^,]+)')
        for arg in args:
            rpath_match = rpath_regex.match(arg)
            if rpath_match:
                for dir in rpath_match.group(1).split(':'):
                    dirs.add(dir)
            runpath_match = runpath_regex.match(arg)
            if runpath_match:
                for dir in runpath_match.group(1).split(':'):
                    # The symbols arg is an rpath if the path is a directory
                    if Path(dir).is_dir():
                        dirs.add(dir)
            symbols_match = symbols_regex.match(arg)
            if symbols_match:
                for dir in symbols_match.group(1).split(':'):
                    # Prevent usage of --just-symbols to specify rpath
                    if Path(dir).is_dir():
                        raise MesonException(f'Invalid arg for --just-symbols, {dir} is a directory.')
        return dirs

    @lru_cache(maxsize=None)
    def rpaths_for_non_system_absolute_shared_libraries(self, target: build.BuildTarget, exclude_system: bool = True) -> 'ImmutableListProtocol[str]':
        paths: OrderedSet[str] = OrderedSet()
        srcdir = self.environment.get_source_dir()

        for dep in target.external_deps:
            if dep.type_name not in {'library', 'pkgconfig', 'cmake'}:
                continue
            for libpath in dep.link_args:
                # For all link args that are absolute paths to a library file, add RPATH args
                if not os.path.isabs(libpath):
                    continue
                libdir = os.path.dirname(libpath)
                if exclude_system and self._libdir_is_system(libdir, target.compilers, self.environment):
                    # No point in adding system paths.
                    continue
                # Don't remove rpaths specified in LDFLAGS.
                if libdir in self.get_external_rpath_dirs(target):
                    continue
                # Windows doesn't support rpaths, but we use this function to
                # emulate rpaths by setting PATH
                # .dll is there for mingw gcc
                # .so's may be extended with version information, e.g. libxyz.so.1.2.3
                if not (
                    os.path.splitext(libpath)[1] in {'.dll', '.lib', '.so', '.dylib'}
                    or re.match(r'.+\.so(\.|$)', os.path.basename(libpath))
                ):
                    continue

                try:
                    commonpath = os.path.commonpath((libdir, srcdir))
                except ValueError: # when paths are on different drives on Windows
                    commonpath = ''

                if commonpath == srcdir:
                    rel_to_src = libdir[len(srcdir) + 1:]
                    assert not os.path.isabs(rel_to_src), f'rel_to_src: {rel_to_src} is absolute'
                    paths.add(os.path.join(self.build_to_src, rel_to_src))
                else:
                    paths.add(libdir)
            # Don't remove rpaths specified by the dependency
            paths.difference_update(self.get_rpath_dirs_from_link_args(dep.link_args))
        for i in chain(target.link_targets, target.link_whole_targets):
            if isinstance(i, build.BuildTarget):
                paths.update(self.rpaths_for_non_system_absolute_shared_libraries(i, exclude_system))
        return list(paths)

    # This may take other types
    def determine_rpath_dirs(self, target: T.Union[build.BuildTarget, build.CustomTarget, build.CustomTargetIndex]
                             ) -> T.Tuple[str, ...]:
        result: OrderedSet[str]
        if self.environment.coredata.get_option(OptionKey('layout')) == 'mirror':
            # Need a copy here
            result = OrderedSet(target.get_link_dep_subdirs())
        else:
            result = OrderedSet()
            result.add('meson-out')
        if isinstance(target, build.BuildTarget):
            result.update(self.rpaths_for_non_system_absolute_shared_libraries(target))
            target.rpath_dirs_to_remove.update([d.encode('utf-8') for d in result])
        return tuple(result)

    @staticmethod
    def canonicalize_filename(fname: str) -> str:
        parts = Path(fname).parts
        hashed = ''
        if len(parts) > 5:
            temp = '/'.join(parts[-5:])
            # is it shorter to hash the beginning of the path?
            if len(fname) > len(temp) + 41:
                hashed = hashlib.sha1(fname.encode('utf-8')).hexdigest() + '_'
                fname = temp
        for ch in ('/', '\\', ':'):
            fname = fname.replace(ch, '_')
        return hashed + fname

    def object_filename_from_source(self, target: build.BuildTarget, source: 'FileOrString', targetdir: T.Optional[str] = None) -> str:
        assert isinstance(source, mesonlib.File)
        if isinstance(target, build.CompileTarget):
            return target.sources_map[source]
        build_dir = self.environment.get_build_dir()
        rel_src = source.rel_to_builddir(self.build_to_src)

        # foo.vala files compile down to foo.c and then foo.c.o, not foo.vala.o
        if rel_src.endswith(('.vala', '.gs')):
            # See description in generate_vala_compile for this logic.
            if source.is_built:
                if os.path.isabs(rel_src):
                    rel_src = rel_src[len(build_dir) + 1:]
                rel_src = os.path.relpath(rel_src, self.get_target_private_dir(target))
            else:
                rel_src = os.path.basename(rel_src)
            # A meson- prefixed directory is reserved; hopefully no-one creates a file name with such a weird prefix.
            gen_source = 'meson-generated_' + rel_src[:-5] + '.c'
        elif source.is_built:
            if os.path.isabs(rel_src):
                rel_src = rel_src[len(build_dir) + 1:]
            # A meson- prefixed directory is reserved; hopefully no-one creates a file name with such a weird prefix.
            gen_source = 'meson-generated_' + os.path.relpath(rel_src, self.get_target_private_dir(target))
        else:
            if os.path.isabs(rel_src):
                # Use the absolute path directly to avoid file name conflicts
                gen_source = rel_src
            else:
                gen_source = os.path.relpath(os.path.join(build_dir, rel_src),
                                             os.path.join(self.environment.get_source_dir(), target.get_subdir()))
        machine = self.environment.machines[target.for_machine]
        ret = self.canonicalize_filename(gen_source) + '.' + machine.get_object_suffix()
        if targetdir is not None:
            return os.path.join(targetdir, ret)
        return ret

    def _determine_ext_objs(self, extobj: 'build.ExtractedObjects', proj_dir_to_build_root: str) -> T.List[str]:
        result: T.List[str] = []

        targetdir = self.get_target_private_dir(extobj.target)

        # Merge sources and generated sources
        raw_sources = list(extobj.srclist)
        for gensrc in extobj.genlist:
            for r in gensrc.get_outputs():
                path = self.get_target_generated_dir(extobj.target, gensrc, r)
                dirpart, fnamepart = os.path.split(path)
                raw_sources.append(File(True, dirpart, fnamepart))

        # Filter out headers and all non-source files
        sources: T.List['FileOrString'] = []
        for s in raw_sources:
            if self.environment.is_source(s):
                sources.append(s)
            elif self.environment.is_object(s):
                result.append(s.relative_name())

        # MSVC generate an object file for PCH
        if extobj.pch and self.target_uses_pch(extobj.target):
            for lang, pch in extobj.target.pch.items():
                compiler = extobj.target.compilers[lang]
                if compiler.get_argument_syntax() == 'msvc':
                    objname = self.get_msvc_pch_objname(lang, pch)
                    result.append(os.path.join(proj_dir_to_build_root, targetdir, objname))

        # extobj could contain only objects and no sources
        if not sources:
            return result

        # With unity builds, sources don't map directly to objects,
        # we only support extracting all the objects in this mode,
        # so just return all object files.
        if extobj.target.is_unity:
            compsrcs = classify_unity_sources(extobj.target.compilers.values(), sources)
            sources = []
            unity_size = extobj.target.get_option(OptionKey('unity_size'))
            assert isinstance(unity_size, int), 'for mypy'

            for comp, srcs in compsrcs.items():
                if comp.language in LANGS_CANT_UNITY:
                    sources += srcs
                    continue
                for i in range((len(srcs) + unity_size - 1) // unity_size):
                    _src = self.get_unity_source_file(extobj.target,
                                                      comp.get_default_suffix(), i)
                    sources.append(_src)

        for osrc in sources:
            objname = self.object_filename_from_source(extobj.target, osrc, targetdir)
            objpath = os.path.join(proj_dir_to_build_root, objname)
            result.append(objpath)

        return result

    def get_pch_include_args(self, compiler: 'Compiler', target: build.BuildTarget) -> T.List[str]:
        args: T.List[str] = []
        pchpath = self.get_target_private_dir(target)
        includeargs = compiler.get_include_args(pchpath, False)
        p = target.get_pch(compiler.get_language())
        if p:
            args += compiler.get_pch_use_args(pchpath, p[0])
        return includeargs + args

    def get_msvc_pch_objname(self, lang: str, pch: T.List[str]) -> str:
        if len(pch) == 1:
            # Same name as in create_msvc_pch_implementation() below.
            return f'meson_pch-{lang}.obj'
        return os.path.splitext(pch[1])[0] + '.obj'

    def create_msvc_pch_implementation(self, target: build.BuildTarget, lang: str, pch_header: str) -> str:
        # We have to include the language in the file name, otherwise
        # pch.c and pch.cpp will both end up as pch.obj in VS backends.
        impl_name = f'meson_pch-{lang}.{lang}'
        pch_rel_to_build = os.path.join(self.get_target_private_dir(target), impl_name)
        # Make sure to prepend the build dir, since the working directory is
        # not defined. Otherwise, we might create the file in the wrong path.
        pch_file = os.path.join(self.build_dir, pch_rel_to_build)
        os.makedirs(os.path.dirname(pch_file), exist_ok=True)

        content = f'#include "{os.path.basename(pch_header)}"'
        pch_file_tmp = pch_file + '.tmp'
        with open(pch_file_tmp, 'w', encoding='utf-8') as f:
            f.write(content)
        mesonlib.replace_if_different(pch_file, pch_file_tmp)
        return pch_rel_to_build

    def target_uses_pch(self, target: build.BuildTarget) -> bool:
        try:
            return T.cast('bool', target.get_option(OptionKey('b_pch')))
        except (KeyError, AttributeError):
            return False

    @staticmethod
    def escape_extra_args(args: T.List[str]) -> T.List[str]:
        # all backslashes in defines are doubly-escaped
        extra_args: T.List[str] = []
        for arg in args:
            if arg.startswith(('-D', '/D')):
                arg = arg.replace('\\', '\\\\')
            extra_args.append(arg)

        return extra_args

    def get_no_stdlib_args(self, target: 'build.BuildTarget', compiler: 'Compiler') -> T.List[str]:
        if compiler.language in self.build.stdlibs[target.for_machine]:
            return compiler.get_no_stdinc_args()
        return []

    def generate_basic_compiler_args(self, target: build.BuildTarget, compiler: 'Compiler') -> 'CompilerArgs':
        # Create an empty commands list, and start adding arguments from
        # various sources in the order in which they must override each other
        # starting from hard-coded defaults followed by build options and so on.
        commands = compiler.compiler_args()

        copt_proxy = target.get_options()
        # First, the trivial ones that are impossible to override.
        #
        # Add -nostdinc/-nostdinc++ if needed; can't be overridden
        commands += self.get_no_stdlib_args(target, compiler)
        # Add things like /NOLOGO or -pipe; usually can't be overridden
        commands += compiler.get_always_args()
        # warning_level is a string, but mypy can't determine that
        commands += compiler.get_warn_args(T.cast('str', target.get_option(OptionKey('warning_level'))))
        # Add -Werror if werror=true is set in the build options set on the
        # command-line or default_options inside project(). This only sets the
        # action to be done for warnings if/when they are emitted, so it's ok
        # to set it after or get_warn_args().
        if target.get_option(OptionKey('werror')):
            commands += compiler.get_werror_args()
        # Add compile args for c_* or cpp_* build options set on the
        # command-line or default_options inside project().
        commands += compiler.get_option_compile_args(copt_proxy)

        optimization = target.get_option(OptionKey('optimization'))
        assert isinstance(optimization, str), 'for mypy'
        commands += compiler.get_optimization_args(optimization)

        debug = target.get_option(OptionKey('debug'))
        assert isinstance(debug, bool), 'for mypy'
        commands += compiler.get_debug_args(debug)

        # Add compile args added using add_project_arguments()
        commands += self.build.get_project_args(compiler, target.subproject, target.for_machine)
        # Add compile args added using add_global_arguments()
        # These override per-project arguments
        commands += self.build.get_global_args(compiler, target.for_machine)
        # Compile args added from the env: CFLAGS/CXXFLAGS, etc, or the cross
        # file. We want these to override all the defaults, but not the
        # per-target compile args.
        commands += self.environment.coredata.get_external_args(target.for_machine, compiler.get_language())
        # Using both /Z7 or /ZI and /Zi at the same times produces a compiler warning.
        # We do not add /Z7 or /ZI by default. If it is being used it is because the user has explicitly enabled it.
        # /Zi needs to be removed in that case to avoid cl's warning to that effect (D9025 : overriding '/Zi' with '/ZI')
        if ('/Zi' in commands) and (('/ZI' in commands) or ('/Z7' in commands)):
            commands.remove('/Zi')
        # Always set -fPIC for shared libraries
        if isinstance(target, build.SharedLibrary):
            commands += compiler.get_pic_args()
        # Set -fPIC for static libraries by default unless explicitly disabled
        if isinstance(target, build.StaticLibrary) and target.pic:
            commands += compiler.get_pic_args()
        elif isinstance(target, (build.StaticLibrary, build.Executable)) and target.pie:
            commands += compiler.get_pie_args()
        # Add compile args needed to find external dependencies. Link args are
        # added while generating the link command.
        # NOTE: We must preserve the order in which external deps are
        # specified, so we reverse the list before iterating over it.
        for dep in reversed(target.get_external_deps()):
            if not dep.found():
                continue

            if compiler.language == 'vala':
                if dep.type_name == 'pkgconfig':
                    assert isinstance(dep, dependencies.ExternalDependency)
                    if dep.name == 'glib-2.0' and dep.version_reqs is not None:
                        for req in dep.version_reqs:
                            if req.startswith(('>=', '==')):
                                commands += ['--target-glib', req[2:]]
                                break
                    commands += ['--pkg', dep.name]
                elif isinstance(dep, dependencies.ExternalLibrary):
                    commands += dep.get_link_args('vala')
            else:
                commands += compiler.get_dependency_compile_args(dep)
            # Qt needs -fPIC for executables
            # XXX: We should move to -fPIC for all executables
            if isinstance(target, build.Executable):
                commands += dep.get_exe_args(compiler)
            # For 'automagic' deps: Boost and GTest. Also dependency('threads').
            # pkg-config puts the thread flags itself via `Cflags:`
        # Fortran requires extra include directives.
        if compiler.language == 'fortran':
            for lt in chain(target.link_targets, target.link_whole_targets):
                priv_dir = self.get_target_private_dir(lt)
                commands += compiler.get_include_args(priv_dir, False)
        return commands

    def build_target_link_arguments(self, compiler: 'Compiler', deps: T.List[build.Target]) -> T.List[str]:
        args: T.List[str] = []
        for d in deps:
            if not d.is_linkable_target():
                raise RuntimeError(f'Tried to link with a non-library target "{d.get_basename()}".')
            arg = self.get_target_filename_for_linking(d)
            if not arg:
                continue
            if compiler.get_language() == 'd':
                arg = '-Wl,' + arg
            elif compiler.get_linker_id() == 'nvlink' and arg.endswith('.a'):
                # We need to pass static archives without -Xlinker= to nvcc,
                # since they may contain relocatable device code. When passing
                # the static archive to nvcc with -Xlinker=, we bypass the
                # frontend which means we lose the opportunity to perform device
                # linking. We only need to do this for static archives, since
                # nvcc doesn't support device linking with dynamic libraries:
                # https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#libraries
                pass
            else:
                arg = compiler.get_linker_lib_prefix() + arg
            args.append(arg)
        return args

    def get_mingw_extra_paths(self, target: build.BuildTarget) -> T.List[str]:
        paths: OrderedSet[str] = OrderedSet()
        # The cross bindir
        root = self.environment.properties[target.for_machine].get_root()
        if root:
            paths.add(os.path.join(root, 'bin'))
        # The toolchain bindir
        sys_root = self.environment.properties[target.for_machine].get_sys_root()
        if sys_root:
            paths.add(os.path.join(sys_root, 'bin'))
        # Get program and library dirs from all target compilers
        if isinstance(target, build.BuildTarget):
            for cc in target.compilers.values():
                paths.update(cc.get_program_dirs(self.environment))
                paths.update(cc.get_library_dirs(self.environment))
        return list(paths)

    @staticmethod
    @lru_cache(maxsize=None)
    def search_dll_path(link_arg: str) -> T.Optional[str]:
        if link_arg.startswith(('-l', '-L')):
            link_arg = link_arg[2:]

        p = Path(link_arg)
        if not p.is_absolute():
            return None

        try:
            p = p.resolve(strict=True)
        except FileNotFoundError:
            return None

        for f in p.parent.glob('*.dll'):
            # path contains dlls
            return str(p.parent)

        if p.is_file():
            p = p.parent
        # Heuristic: replace *last* occurrence of '/lib'
        binpath = Path('/bin'.join(p.as_posix().rsplit('/lib', maxsplit=1)))
        for _ in binpath.glob('*.dll'):
            return str(binpath)

        return None

    @classmethod
    @lru_cache(maxsize=None)
    def extract_dll_paths(cls, target: build.BuildTarget) -> T.Set[str]:
        """Find paths to all DLLs needed for a given target, since
        we link against import libs, and we don't know the actual
        path of the DLLs.

        1. If there are DLLs in the same directory than the .lib dir, use it
        2. If there is a sibbling directory named 'bin' with DLLs in it, use it
        """
        results = set()
        for dep in target.external_deps:

            if dep.type_name == 'pkgconfig':
                # If by chance pkg-config knows the bin dir...
                bindir = dep.get_variable(pkgconfig='bindir', default_value='')
                if bindir:
                    results.add(bindir)
                    continue

            results.update(filter(None, map(cls.search_dll_path, dep.link_args)))  # pylint: disable=bad-builtin

        for i in chain(target.link_targets, target.link_whole_targets):
            if isinstance(i, build.BuildTarget):
                results.update(cls.extract_dll_paths(i))

        return results

    def determine_windows_extra_paths(
            self, target: T.Union[build.BuildTarget, build.CustomTarget, build.CustomTargetIndex, programs.ExternalProgram, mesonlib.File, str],
            extra_bdeps: T.Sequence[T.Union[build.BuildTarget, build.CustomTarget, build.CustomTargetIndex]]) -> T.List[str]:
        """On Windows there is no such thing as an rpath.

        We must determine all locations of DLLs that this exe
        links to and return them so they can be used in unit
        tests.
        """
        result: T.Set[str] = set()
        prospectives: T.Set[build.BuildTargetTypes] = set()
        if isinstance(target, build.BuildTarget):
            prospectives.update(target.get_transitive_link_deps())
            # External deps
            result.update(self.extract_dll_paths(target))

        for bdep in extra_bdeps:
            prospectives.add(bdep)
            if isinstance(bdep, build.BuildTarget):
                prospectives.update(bdep.get_transitive_link_deps())
        # Internal deps
        for ld in prospectives:
            dirseg = os.path.join(self.environment.get_build_dir(), self.get_target_dir(ld))
            result.add(dirseg)
        if (isinstance(target, build.BuildTarget) and
                not self.environment.machines.matches_build_machine(target.for_machine)):
            result.update(self.get_mingw_extra_paths(target))
        return list(result)

    def write_benchmark_file(self, datafile: T.BinaryIO) -> None:
        self.write_test_serialisation(self.build.get_benchmarks(), datafile)

    def write_test_file(self, datafile: T.BinaryIO) -> None:
        self.write_test_serialisation(self.build.get_tests(), datafile)

    def create_test_serialisation(self, tests: T.List['Test']) -> T.List[TestSerialisation]:
        arr: T.List[TestSerialisation] = []
        for t in sorted(tests, key=lambda tst: -1 * tst.priority):
            exe = t.get_exe()
            if isinstance(exe, programs.ExternalProgram):
                cmd = exe.get_command()
            else:
                cmd = [os.path.join(self.environment.get_build_dir(), self.get_target_filename(exe))]
            if isinstance(exe, (build.BuildTarget, programs.ExternalProgram)):
                test_for_machine = exe.for_machine
            else:
                # E.g. an external verifier or simulator program run on a generated executable.
                # Can always be run without a wrapper.
                test_for_machine = MachineChoice.BUILD

            # we allow passing compiled executables to tests, which may be cross built.
            # We need to consider these as well when considering whether the target is cross or not.
            for a in t.cmd_args:
                if isinstance(a, build.BuildTarget):
                    if a.for_machine is MachineChoice.HOST:
                        test_for_machine = MachineChoice.HOST
                        break

            is_cross = self.environment.is_cross_build(test_for_machine)
            exe_wrapper = self.environment.get_exe_wrapper()
            machine = self.environment.machines[exe.for_machine]
            if machine.is_windows() or machine.is_cygwin():
                extra_bdeps: T.List[T.Union[build.BuildTarget, build.CustomTarget, build.CustomTargetIndex]] = []
                if isinstance(exe, build.CustomTarget):
                    extra_bdeps = list(exe.get_transitive_build_target_deps())
                extra_paths = self.determine_windows_extra_paths(exe, extra_bdeps)
                for a in t.cmd_args:
                    if isinstance(a, build.BuildTarget):
                        for p in self.determine_windows_extra_paths(a, []):
                            if p not in extra_paths:
                                extra_paths.append(p)
            else:
                extra_paths = []

            cmd_args: T.List[str] = []
            depends: T.Set[build.Target] = set(t.depends)
            if isinstance(exe, build.Target):
                depends.add(exe)
            for a in t.cmd_args:
                if isinstance(a, build.Target):
                    depends.add(a)
                elif isinstance(a, build.CustomTargetIndex):
                    depends.add(a.target)

                if isinstance(a, mesonlib.File):
                    a = os.path.join(self.environment.get_build_dir(), a.rel_to_builddir(self.build_to_src))
                    cmd_args.append(a)
                elif isinstance(a, str):
                    cmd_args.append(a)
                elif isinstance(a, (build.Target, build.CustomTargetIndex)):
                    cmd_args.extend(self.construct_target_rel_paths(a, t.workdir))
                elif isinstance(a, programs.ExternalProgram):
                    cmd_args.extend(a.get_command())
                else:
                    raise MesonException('Bad object in test command.')

            t_env = copy.deepcopy(t.env)
            if not machine.is_windows() and not machine.is_cygwin() and not machine.is_darwin():
                ld_lib_path_libs: T.Set[build.SharedLibrary] = set()
                for d in depends:
                    if isinstance(d, build.BuildTarget):
                        for l in d.get_all_link_deps():
                            if isinstance(l, build.SharedLibrary):
                                ld_lib_path_libs.add(l)

                env_build_dir = self.environment.get_build_dir()
                ld_lib_path: T.Set[str] = set(os.path.join(env_build_dir, l.get_subdir()) for l in ld_lib_path_libs)

                if ld_lib_path:
                    t_env.prepend('LD_LIBRARY_PATH', list(ld_lib_path), ':')

            ts = TestSerialisation(t.get_name(), t.project_name, t.suite, cmd, is_cross,
                                   exe_wrapper, self.environment.need_exe_wrapper(),
                                   t.is_parallel, cmd_args, t_env,
                                   t.should_fail, t.timeout, t.workdir,
                                   extra_paths, t.protocol, t.priority,
                                   isinstance(exe, (build.Target, build.CustomTargetIndex)),
                                   isinstance(exe, build.Executable),
                                   [x.get_id() for x in depends],
                                   self.environment.coredata.version,
                                   t.verbose)
            arr.append(ts)
        return arr

    def write_test_serialisation(self, tests: T.List['Test'], datafile: T.BinaryIO) -> None:
        pickle.dump(self.create_test_serialisation(tests), datafile)

    def construct_target_rel_paths(self, t: T.Union[build.Target, build.CustomTargetIndex], workdir: T.Optional[str]) -> T.List[str]:
        target_dir = self.get_target_dir(t)
        # ensure that test executables can be run when passed as arguments
        if isinstance(t, build.Executable) and workdir is None:
            target_dir = target_dir or '.'

        if isinstance(t, build.BuildTarget):
            outputs = [t.get_filename()]
        else:
            assert isinstance(t, (build.CustomTarget, build.CustomTargetIndex))
            outputs = t.get_outputs()

        outputs = [os.path.join(target_dir, x) for x in outputs]
        if workdir is not None:
            assert os.path.isabs(workdir)
            outputs = [os.path.join(self.environment.get_build_dir(), x) for x in outputs]
            outputs = [os.path.relpath(x, workdir) for x in outputs]
        return outputs

    def generate_depmf_install(self, d: InstallData) -> None:
        depmf_path = self.build.dep_manifest_name
        if depmf_path is None:
            option_dir = self.environment.coredata.get_option(OptionKey('licensedir'))
            assert isinstance(option_dir, str), 'for mypy'
            if option_dir:
                depmf_path = os.path.join(option_dir, 'depmf.json')
            else:
                return
        ifilename = os.path.join(self.environment.get_build_dir(), 'depmf.json')
        ofilename = os.path.join(self.environment.get_prefix(), depmf_path)
        odirname = os.path.join(self.environment.get_prefix(), os.path.dirname(depmf_path))
        out_name = os.path.join('{prefix}', depmf_path)
        out_dir = os.path.join('{prefix}', os.path.dirname(depmf_path))
        mfobj = {'type': 'dependency manifest', 'version': '1.0',
                 'projects': {k: v.to_json() for k, v in self.build.dep_manifest.items()}}
        with open(ifilename, 'w', encoding='utf-8') as f:
            f.write(json.dumps(mfobj))
        # Copy file from, to, and with mode unchanged
        d.data.append(InstallDataBase(ifilename, ofilename, out_name, None, '',
                                      tag='devel', data_type='depmf'))
        for m in self.build.dep_manifest.values():
            for ifilename, name in m.license_files:
                ofilename = os.path.join(odirname, name.relative_name())
                out_name = os.path.join(out_dir, name.relative_name())
                d.data.append(InstallDataBase(ifilename, ofilename, out_name, None,
                                              m.subproject, tag='devel', data_type='depmf'))

    def get_regen_filelist(self) -> T.List[str]:
        '''List of all files whose alteration means that the build
        definition needs to be regenerated.'''
        deps = OrderedSet([str(Path(self.build_to_src) / df)
                           for df in self.interpreter.get_build_def_files()])
        if self.environment.is_cross_build():
            deps.update(self.environment.coredata.cross_files)
        deps.update(self.environment.coredata.config_files)
        deps.add('meson-private/coredata.dat')
        self.check_clock_skew(deps)
        return list(deps)

    def generate_regen_info(self) -> None:
        deps = self.get_regen_filelist()
        regeninfo = RegenInfo(self.environment.get_source_dir(),
                              self.environment.get_build_dir(),
                              deps)
        filename = os.path.join(self.environment.get_scratch_dir(),
                                'regeninfo.dump')
        with open(filename, 'wb') as f:
            pickle.dump(regeninfo, f)

    def check_clock_skew(self, file_list: T.Iterable[str]) -> None:
        # If a file that leads to reconfiguration has a time
        # stamp in the future, it will trigger an eternal reconfigure
        # loop.
        import time
        now = time.time()
        for f in file_list:
            absf = os.path.join(self.environment.get_build_dir(), f)
            ftime = os.path.getmtime(absf)
            delta = ftime - now
            # On Windows disk time stamps sometimes point
            # to the future by a minuscule amount, less than
            # 0.001 seconds. I don't know why.
            if delta > 0.001:
                raise MesonException(f'Clock skew detected. File {absf} has a time stamp {delta:.4f}s in the future.')

    def build_target_to_cmd_array(self, bt: T.Union[build.BuildTarget, programs.ExternalProgram]) -> T.List[str]:
        if isinstance(bt, build.BuildTarget):
            arr = [os.path.join(self.environment.get_build_dir(), self.get_target_filename(bt))]
        else:
            arr = bt.get_command()
        return arr

    def replace_extra_args(self, args: T.List[str], genlist: 'build.GeneratedList') -> T.List[str]:
        final_args: T.List[str] = []
        for a in args:
            if a == '@EXTRA_ARGS@':
                final_args += genlist.get_extra_args()
            else:
                final_args.append(a)
        return final_args

    def replace_outputs(self, args: T.List[str], private_dir: str, output_list: T.List[str]) -> T.List[str]:
        newargs: T.List[str] = []
        regex = re.compile(r'@OUTPUT(\d+)@')
        for arg in args:
            m = regex.search(arg)
            while m is not None:
                index = int(m.group(1))
                src = f'@OUTPUT{index}@'
                arg = arg.replace(src, os.path.join(private_dir, output_list[index]))
                m = regex.search(arg)
            newargs.append(arg)
        return newargs

    def get_build_by_default_targets(self) -> 'T.OrderedDict[str, T.Union[build.BuildTarget, build.CustomTarget]]':
        result: 'T.OrderedDict[str, T.Union[build.BuildTarget, build.CustomTarget]]' = OrderedDict()
        # Get all build and custom targets that must be built by default
        for name, b in self.build.get_targets().items():
            if b.build_by_default:
                result[name] = b
        return result

    def get_testlike_targets(self, benchmark: bool = False) -> T.OrderedDict[str, T.Union[build.BuildTarget, build.CustomTarget]]:
        result: T.OrderedDict[str, T.Union[build.BuildTarget, build.CustomTarget]] = OrderedDict()
        targets = self.build.get_benchmarks() if benchmark else self.build.get_tests()
        for t in targets:
            exe = t.exe
            if isinstance(exe, (build.CustomTarget, build.BuildTarget)):
                result[exe.get_id()] = exe
            for arg in t.cmd_args:
                if not isinstance(arg, (build.CustomTarget, build.BuildTarget)):
                    continue
                result[arg.get_id()] = arg
            for dep in t.depends:
                assert isinstance(dep, (build.CustomTarget, build.BuildTarget, build.CustomTargetIndex))
                result[dep.get_id()] = dep
        return result

    @lru_cache(maxsize=None)
    def get_custom_target_provided_by_generated_source(self, generated_source: build.CustomTarget) -> 'ImmutableListProtocol[str]':
        libs: T.List[str] = []
        for f in generated_source.get_outputs():
            if self.environment.is_library(f):
                libs.append(os.path.join(self.get_target_dir(generated_source), f))
        return libs

    @lru_cache(maxsize=None)
    def get_custom_target_provided_libraries(self, target: T.Union[build.BuildTarget, build.CustomTarget]) -> 'ImmutableListProtocol[str]':
        libs: T.List[str] = []
        for t in target.get_generated_sources():
            if not isinstance(t, build.CustomTarget):
                continue
            libs.extend(self.get_custom_target_provided_by_generated_source(t))
        return libs

    def get_custom_target_sources(self, target: build.CustomTarget) -> T.List[str]:
        '''
        Custom target sources can be of various object types; strings, File,
        BuildTarget, even other CustomTargets.
        Returns the path to them relative to the build root directory.
        '''
        srcs: T.List[str] = []
        for i in target.get_sources():
            if isinstance(i, str):
                fname = [os.path.join(self.build_to_src, target.subdir, i)]
            elif isinstance(i, build.BuildTarget):
                fname = [self.get_target_filename(i)]
            elif isinstance(i, (build.CustomTarget, build.CustomTargetIndex)):
                fname = [os.path.join(self.get_custom_target_output_dir(i), p) for p in i.get_outputs()]
            elif isinstance(i, build.GeneratedList):
                fname = [os.path.join(self.get_target_private_dir(target), p) for p in i.get_outputs()]
            elif isinstance(i, build.ExtractedObjects):
                fname = self.determine_ext_objs(i)
            elif isinstance(i, programs.ExternalProgram):
                assert i.found(), "This shouldn't be possible"
                assert i.path is not None, 'for mypy'
                fname = [i.path]
            else:
                fname = [i.rel_to_builddir(self.build_to_src)]
            if target.absolute_paths:
                fname = [os.path.join(self.environment.get_build_dir(), f) for f in fname]
            srcs += fname
        return srcs

    def get_target_depend_files(self, target: T.Union[build.CustomTarget, build.BuildTarget], absolute_paths: bool = False) -> T.List[str]:
        deps: T.List[str] = []
        for i in target.depend_files:
            if isinstance(i, mesonlib.File):
                if absolute_paths:
                    deps.append(i.absolute_path(self.environment.get_source_dir(),
                                                self.environment.get_build_dir()))
                else:
                    deps.append(i.rel_to_builddir(self.build_to_src))
            else:
                if absolute_paths:
                    deps.append(os.path.join(self.environment.get_source_dir(), target.subdir, i))
                else:
                    deps.append(os.path.join(self.build_to_src, target.subdir, i))
        return deps

    def get_custom_target_output_dir(self, target: T.Union[build.Target, build.CustomTargetIndex]) -> str:
        # The XCode backend is special. A target foo/bar does
        # not go to ${BUILDDIR}/foo/bar but instead to
        # ${BUILDDIR}/${BUILDTYPE}/foo/bar.
        # Currently we set the include dir to be the former,
        # and not the latter. Thus we need this extra customisation
        # point. If in the future we make include dirs et al match
        # ${BUILDDIR}/${BUILDTYPE} instead, this becomes unnecessary.
        return self.get_target_dir(target)

    @lru_cache(maxsize=None)
    def get_normpath_target(self, source: str) -> str:
        return os.path.normpath(source)

    def get_custom_target_dirs(self, target: build.CustomTarget, compiler: 'Compiler', *,
                               absolute_path: bool = False) -> T.List[str]:
        custom_target_include_dirs: T.List[str] = []
        for i in target.get_generated_sources():
            # Generator output goes into the target private dir which is
            # already in the include paths list. Only custom targets have their
            # own target build dir.
            if not isinstance(i, (build.CustomTarget, build.CustomTargetIndex)):
                continue
            idir = self.get_normpath_target(self.get_custom_target_output_dir(i))
            if not idir:
                idir = '.'
            if absolute_path:
                idir = os.path.join(self.environment.get_build_dir(), idir)
            if idir not in custom_target_include_dirs:
                custom_target_include_dirs.append(idir)
        return custom_target_include_dirs

    def get_custom_target_dir_include_args(
            self, target: build.CustomTarget, compiler: 'Compiler', *,
            absolute_path: bool = False) -> T.List[str]:
        incs: T.List[str] = []
        for i in self.get_custom_target_dirs(target, compiler, absolute_path=absolute_path):
            incs += compiler.get_include_args(i, False)
        return incs

    def eval_custom_target_command(
            self, target: build.CustomTarget, absolute_outputs: bool = False) -> \
            T.Tuple[T.List[str], T.List[str], T.List[str]]:
        # We want the outputs to be absolute only when using the VS backend
        # XXX: Maybe allow the vs backend to use relative paths too?
        source_root = self.build_to_src
        build_root = '.'
        outdir = self.get_custom_target_output_dir(target)
        if absolute_outputs:
            source_root = self.environment.get_source_dir()
            build_root = self.environment.get_build_dir()
            outdir = os.path.join(self.environment.get_build_dir(), outdir)
        outputs = [os.path.join(outdir, i) for i in target.get_outputs()]
        inputs = self.get_custom_target_sources(target)
        # Evaluate the command list
        cmd: T.List[str] = []
        for i in target.command:
            if isinstance(i, build.BuildTarget):
                cmd += self.build_target_to_cmd_array(i)
                continue
            elif isinstance(i, build.CustomTarget):
                # GIR scanner will attempt to execute this binary but
                # it assumes that it is in path, so always give it a full path.
                tmp = i.get_outputs()[0]
                i = os.path.join(self.get_custom_target_output_dir(i), tmp)
            elif isinstance(i, mesonlib.File):
                i = i.rel_to_builddir(self.build_to_src)
                if target.absolute_paths or absolute_outputs:
                    i = os.path.join(self.environment.get_build_dir(), i)
            # FIXME: str types are blindly added ignoring 'target.absolute_paths'
            # because we can't know if they refer to a file or just a string
            elif isinstance(i, str):
                if '@SOURCE_ROOT@' in i:
                    i = i.replace('@SOURCE_ROOT@', source_root)
                if '@BUILD_ROOT@' in i:
                    i = i.replace('@BUILD_ROOT@', build_root)
                if '@CURRENT_SOURCE_DIR@' in i:
                    i = i.replace('@CURRENT_SOURCE_DIR@', os.path.join(source_root, target.subdir))
                if '@DEPFILE@' in i:
                    if target.depfile is None:
                        msg = f'Custom target {target.name!r} has @DEPFILE@ but no depfile ' \
                              'keyword argument.'
                        raise MesonException(msg)
                    dfilename = os.path.join(outdir, target.depfile)
                    i = i.replace('@DEPFILE@', dfilename)
                if '@PRIVATE_DIR@' in i:
                    pdir = self.get_target_private_dir_abs(target)
                    os.makedirs(pdir, exist_ok=True)
                    if not target.absolute_paths:
                        pdir = self.get_target_private_dir(target)
                    i = i.replace('@PRIVATE_DIR@', pdir)
            else:
                raise RuntimeError(f'Argument {i} is of unknown type {type(i)}')
            cmd.append(i)
        # Substitute the rest of the template strings
        values = mesonlib.get_filenames_templates_dict(inputs, outputs)
        cmd = mesonlib.substitute_values(cmd, values)
        # This should not be necessary but removing it breaks
        # building GStreamer on Windows. The underlying issue
        # is problems with quoting backslashes on Windows
        # which is the seventh circle of hell. The downside is
        # that this breaks custom targets whose command lines
        # have backslashes. If you try to fix this be sure to
        # check that it does not break GST.
        #
        # The bug causes file paths such as c:\foo to get escaped
        # into c:\\foo.
        #
        # Unfortunately we have not been able to come up with an
        # isolated test case for this so unless you manage to come up
        # with one, the only way is to test the building with Gst's
        # setup. Note this in your MR or ping us and we will get it
        # fixed.
        #
        # https://github.com/mesonbuild/meson/pull/737
        cmd = [i.replace('\\', '/') for i in cmd]
        return inputs, outputs, cmd

    def get_introspect_command(self) -> str:
        return ' '.join(shlex.quote(x) for x in self.environment.get_build_command() + ['introspect'])

    def get_run_target_env(self, target: build.RunTarget) -> mesonlib.EnvironmentVariables:
        env = target.env if target.env else mesonlib.EnvironmentVariables()
        if target.default_env:
            env.set('MESON_SOURCE_ROOT', [self.environment.get_source_dir()])
            env.set('MESON_BUILD_ROOT', [self.environment.get_build_dir()])
            env.set('MESON_SUBDIR', [target.subdir])
            env.set('MESONINTROSPECT', [self.get_introspect_command()])
        return env

    def run_postconf_scripts(self) -> None:
        from ..scripts.meson_exe import run_exe
        env = {'MESON_SOURCE_ROOT': self.environment.get_source_dir(),
               'MESON_BUILD_ROOT': self.environment.get_build_dir(),
               'MESONINTROSPECT': self.get_introspect_command(),
               }

        for s in self.build.postconf_scripts:
            name = ' '.join(s.cmd_args)
            mlog.log(f'Running postconf script {name!r}')
            rc = run_exe(s, env)
            if rc != 0:
                raise MesonException(f'Postconf script \'{name}\' failed with exit code {rc}.')

    def create_install_data(self) -> InstallData:
        strip_bin = self.environment.lookup_binary_entry(MachineChoice.HOST, 'strip')
        if strip_bin is None:
            if self.environment.is_cross_build():
                mlog.warning('Cross file does not specify strip binary, result will not be stripped.')
            else:
                # TODO go through all candidates, like others
                strip_bin = [detect.defaults['strip'][0]]

        umask = self.environment.coredata.get_option(OptionKey('install_umask'))
        assert isinstance(umask, (str, int)), 'for mypy'

        d = InstallData(self.environment.get_source_dir(),
                        self.environment.get_build_dir(),
                        self.environment.get_prefix(),
                        self.environment.get_libdir(),
                        strip_bin,
                        umask,
                        self.environment.get_build_command() + ['introspect'],
                        self.environment.coredata.version)
        self.generate_depmf_install(d)
        self.generate_target_install(d)
        self.generate_header_install(d)
        self.generate_man_install(d)
        self.generate_emptydir_install(d)
        self.generate_data_install(d)
        self.generate_symlink_install(d)
        self.generate_custom_install_script(d)
        self.generate_subdir_install(d)
        return d

    def create_install_data_files(self) -> None:
        install_data_file = os.path.join(self.environment.get_scratch_dir(), 'install.dat')
        with open(install_data_file, 'wb') as ofile:
            pickle.dump(self.create_install_data(), ofile)

    def guess_install_tag(self, fname: str, outdir: T.Optional[str] = None) -> T.Optional[str]:
        prefix = self.environment.get_prefix()
        bindir = Path(prefix, self.environment.get_bindir())
        libdir = Path(prefix, self.environment.get_libdir())
        incdir = Path(prefix, self.environment.get_includedir())
        _ldir = self.environment.coredata.get_option(OptionKey('localedir'))
        assert isinstance(_ldir, str), 'for mypy'
        localedir = Path(prefix, _ldir)
        dest_path = Path(prefix, outdir, Path(fname).name) if outdir else Path(prefix, fname)
        if bindir in dest_path.parents:
            return 'runtime'
        elif libdir in dest_path.parents:
            if dest_path.suffix in {'.a', '.pc'}:
                return 'devel'
            elif dest_path.suffix in {'.so', '.dll'}:
                return 'runtime'
        elif incdir in dest_path.parents:
            return 'devel'
        elif localedir in dest_path.parents:
            return 'i18n'
        elif 'installed-tests' in dest_path.parts:
            return 'tests'
        elif 'systemtap' in dest_path.parts:
            return 'systemtap'
        mlog.debug('Failed to guess install tag for', dest_path)
        return None

    def generate_target_install(self, d: InstallData) -> None:
        for t in self.build.get_targets().values():
            if not t.should_install():
                continue
            outdirs, install_dir_names, custom_install_dir = t.get_install_dir()
            # Sanity-check the outputs and install_dirs
            num_outdirs, num_out = len(outdirs), len(t.get_outputs())
            if num_outdirs not in {1, num_out}:
                m = 'Target {!r} has {} outputs: {!r}, but only {} "install_dir"s were found.\n' \
                    "Pass 'false' for outputs that should not be installed and 'true' for\n" \
                    'using the default installation directory for an output.'
                raise MesonException(m.format(t.name, num_out, t.get_outputs(), num_outdirs))
            assert len(t.install_tag) == num_out
            install_mode = t.get_custom_install_mode()
            # because mypy gets confused type narrowing in lists
            first_outdir = outdirs[0]
            first_outdir_name = install_dir_names[0]

            # Install the target output(s)
            if isinstance(t, build.BuildTarget):
                # In general, stripping static archives is tricky and full of pitfalls.
                # Wholesale stripping of static archives with a command such as
                #
                #   strip libfoo.a
                #
                # is broken, as GNU's strip will remove *every* symbol in a static
                # archive. One solution to this nonintuitive behaviour would be
                # to only strip local/debug symbols. Unfortunately, strip arguments
                # are not specified by POSIX and therefore not portable. GNU's `-g`
                # option (i.e. remove debug symbols) is equivalent to Apple's `-S`.
                #
                # TODO: Create GNUStrip/AppleStrip/etc. hierarchy for more
                #       fine-grained stripping of static archives.
                can_strip = not isinstance(t, build.StaticLibrary)
                should_strip = can_strip and t.get_option(OptionKey('strip'))
                assert isinstance(should_strip, bool), 'for mypy'
                # Install primary build output (library/executable/jar, etc)
                # Done separately because of strip/aliases/rpath
                if first_outdir is not False:
                    tag = t.install_tag[0] or ('devel' if isinstance(t, build.StaticLibrary) else 'runtime')
                    mappings = t.get_link_deps_mapping(d.prefix)
                    i = TargetInstallData(self.get_target_filename(t), first_outdir,
                                          first_outdir_name,
                                          should_strip, mappings, t.rpath_dirs_to_remove,
                                          t.install_rpath, install_mode, t.subproject,
                                          tag=tag, can_strip=can_strip)
                    d.targets.append(i)

                    for alias, to, tag in t.get_aliases():
                        alias = os.path.join(first_outdir, alias)
                        s = InstallSymlinkData(to, alias, first_outdir, t.subproject, tag, allow_missing=True)
                        d.symlinks.append(s)

                    if isinstance(t, (build.SharedLibrary, build.SharedModule, build.Executable)):
                        # On toolchains/platforms that use an import library for
                        # linking (separate from the shared library with all the
                        # code), we need to install that too (dll.a/.lib).
                        if t.get_import_filename():
                            if custom_install_dir:
                                # If the DLL is installed into a custom directory,
                                # install the import library into the same place so
                                # it doesn't go into a surprising place
                                implib_install_dir = first_outdir
                            else:
                                implib_install_dir = self.environment.get_import_lib_dir()
                            # Install the import library; may not exist for shared modules
                            i = TargetInstallData(self.get_target_filename_for_linking(t),
                                                  implib_install_dir, first_outdir_name,
                                                  False, {}, set(), '', install_mode,
                                                  t.subproject, optional=isinstance(t, build.SharedModule),
                                                  tag='devel')
                            d.targets.append(i)

                        if not should_strip and t.get_debug_filename():
                            debug_file = os.path.join(self.get_target_dir(t), t.get_debug_filename())
                            i = TargetInstallData(debug_file, first_outdir,
                                                  first_outdir_name,
                                                  False, {}, set(), '',
                                                  install_mode, t.subproject,
                                                  optional=True, tag='devel')
                            d.targets.append(i)
                # Install secondary outputs. Only used for Vala right now.
                if num_outdirs > 1:
                    for output, outdir, outdir_name, tag in zip(t.get_outputs()[1:], outdirs[1:], install_dir_names[1:], t.install_tag[1:]):
                        # User requested that we not install this output
                        if outdir is False:
                            continue
                        f = os.path.join(self.get_target_dir(t), output)
                        i = TargetInstallData(f, outdir, outdir_name, False, {}, set(), None,
                                              install_mode, t.subproject,
                                              tag=tag)
                        d.targets.append(i)
            elif isinstance(t, build.CustomTarget):
                # If only one install_dir is specified, assume that all
                # outputs will be installed into it. This is for
                # backwards-compatibility and because it makes sense to
                # avoid repetition since this is a common use-case.
                #
                # To selectively install only some outputs, pass `false` as
                # the install_dir for the corresponding output by index
                #
                # XXX: this wouldn't be needed if we just always matches outdirs
                # to the length of outputs…
                if num_outdirs == 1 and num_out > 1:
                    if first_outdir is not False:
                        for output, tag in zip(t.get_outputs(), t.install_tag):
                            tag = tag or self.guess_install_tag(output, first_outdir)
                            f = os.path.join(self.get_target_dir(t), output)
                            i = TargetInstallData(f, first_outdir, first_outdir_name,
                                                  False, {}, set(), None, install_mode,
                                                  t.subproject, optional=not t.build_by_default,
                                                  tag=tag)
                            d.targets.append(i)
                else:
                    for output, outdir, outdir_name, tag in zip(t.get_outputs(), outdirs, install_dir_names, t.install_tag):
                        # User requested that we not install this output
                        if outdir is False:
                            continue
                        tag = tag or self.guess_install_tag(output, outdir)
                        f = os.path.join(self.get_target_dir(t), output)
                        i = TargetInstallData(f, outdir, outdir_name,
                                              False, {}, set(), None, install_mode,
                                              t.subproject, optional=not t.build_by_default,
                                              tag=tag)
                        d.targets.append(i)

    def generate_custom_install_script(self, d: InstallData) -> None:
        d.install_scripts = self.build.install_scripts
        for i in d.install_scripts:
            if not i.tag:
                mlog.debug('Failed to guess install tag for install script:', ' '.join(i.cmd_args))

    def generate_header_install(self, d: InstallData) -> None:
        incroot = self.environment.get_includedir()
        headers = self.build.get_headers()

        srcdir = self.environment.get_source_dir()
        builddir = self.environment.get_build_dir()
        for h in headers:
            outdir = outdir_name = h.get_custom_install_dir()
            if outdir is None:
                subdir = h.get_install_subdir()
                if subdir is None:
                    outdir = incroot
                    outdir_name = '{includedir}'
                else:
                    outdir = os.path.join(incroot, subdir)
                    outdir_name = os.path.join('{includedir}', subdir)

            for f in h.get_sources():
                abspath = f.absolute_path(srcdir, builddir)
                i = InstallDataBase(abspath, outdir, outdir_name, h.get_custom_install_mode(), h.subproject, tag='devel', follow_symlinks=h.follow_symlinks)
                d.headers.append(i)

    def generate_man_install(self, d: InstallData) -> None:
        manroot = self.environment.get_mandir()
        man = self.build.get_man()
        for m in man:
            for f in m.get_sources():
                num = f.split('.')[-1]
                subdir = m.get_custom_install_dir()
                if subdir is None:
                    if m.locale:
                        subdir = os.path.join('{mandir}', m.locale, 'man' + num)
                    else:
                        subdir = os.path.join('{mandir}', 'man' + num)
                fname = f.fname
                if m.locale: # strip locale from file name
                    fname = fname.replace(f'.{m.locale}', '')
                srcabs = f.absolute_path(self.environment.get_source_dir(), self.environment.get_build_dir())
                dstname = os.path.join(subdir, os.path.basename(fname))
                dstabs = dstname.replace('{mandir}', manroot)
                i = InstallDataBase(srcabs, dstabs, dstname, m.get_custom_install_mode(), m.subproject, tag='man')
                d.man.append(i)

    def generate_emptydir_install(self, d: InstallData) -> None:
        emptydir: T.List[build.EmptyDir] = self.build.get_emptydir()
        for e in emptydir:
            tag = e.install_tag or self.guess_install_tag(e.path)
            i = InstallEmptyDir(e.path, e.install_mode, e.subproject, tag)
            d.emptydir.append(i)

    def generate_data_install(self, d: InstallData) -> None:
        data = self.build.get_data()
        srcdir = self.environment.get_source_dir()
        builddir = self.environment.get_build_dir()
        for de in data:
            assert isinstance(de, build.Data)
            subdir = de.install_dir
            subdir_name = de.install_dir_name
            for src_file, dst_name in zip(de.sources, de.rename):
                assert isinstance(src_file, mesonlib.File)
                dst_abs = os.path.join(subdir, dst_name)
                dstdir_name = os.path.join(subdir_name, dst_name)
                tag = de.install_tag or self.guess_install_tag(dst_abs)
                i = InstallDataBase(src_file.absolute_path(srcdir, builddir), dst_abs, dstdir_name,
                                    de.install_mode, de.subproject, tag=tag, data_type=de.data_type,
                                    follow_symlinks=de.follow_symlinks)
                d.data.append(i)

    def generate_symlink_install(self, d: InstallData) -> None:
        links: T.List[build.SymlinkData] = self.build.get_symlinks()
        for l in links:
            assert isinstance(l, build.SymlinkData)
            install_dir = l.install_dir
            name_abs = os.path.join(install_dir, l.name)
            tag = l.install_tag or self.guess_install_tag(name_abs)
            s = InstallSymlinkData(l.target, name_abs, install_dir, l.subproject, tag)
            d.symlinks.append(s)

    def generate_subdir_install(self, d: InstallData) -> None:
        for sd in self.build.get_install_subdirs():
            if sd.from_source_dir:
                from_dir = self.environment.get_source_dir()
            else:
                from_dir = self.environment.get_build_dir()
            src_dir = os.path.join(from_dir,
                                   sd.source_subdir,
                                   sd.installable_subdir).rstrip('/')
            dst_dir = os.path.join(self.environment.get_prefix(),
                                   sd.install_dir)
            dst_name = os.path.join('{prefix}', sd.install_dir)
            if sd.install_dir != sd.install_dir_name:
                dst_name = sd.install_dir_name
            if not sd.strip_directory:
                dst_dir = os.path.join(dst_dir, os.path.basename(src_dir))
                dst_name = os.path.join(dst_name, os.path.basename(src_dir))
            tag = sd.install_tag or self.guess_install_tag(os.path.join(sd.install_dir, 'dummy'))
            i = SubdirInstallData(src_dir, dst_dir, dst_name, sd.install_mode, sd.exclude, sd.subproject, tag,
                                  follow_symlinks=sd.follow_symlinks)
            d.install_subdirs.append(i)

    def get_introspection_data(self, target_id: str, target: build.Target) -> T.List['TargetIntrospectionData']:
        '''
        Returns a list of source dicts with the following format for a given target:
        [
            {
                "language": "<LANG>",
                "compiler": ["result", "of", "comp.get_exelist()"],
                "parameters": ["list", "of", "compiler", "parameters],
                "sources": ["list", "of", "all", "<LANG>", "source", "files"],
                "generated_sources": ["list", "of", "generated", "source", "files"]
            }
        ]

        This is a limited fallback / reference implementation. The backend should override this method.
        '''
        if isinstance(target, (build.CustomTarget, build.BuildTarget)):
            source_list_raw = target.sources
            source_list = []
            for j in source_list_raw:
                if isinstance(j, mesonlib.File):
                    source_list += [j.absolute_path(self.source_dir, self.build_dir)]
                elif isinstance(j, str):
                    source_list += [os.path.join(self.source_dir, j)]
                elif isinstance(j, (build.CustomTarget, build.BuildTarget)):
                    source_list += [os.path.join(self.build_dir, j.get_subdir(), o) for o in j.get_outputs()]
            source_list = [os.path.normpath(s) for s in source_list]

            compiler: T.List[str] = []
            if isinstance(target, build.CustomTarget):
                tmp_compiler = target.command
                for j in tmp_compiler:
                    if isinstance(j, mesonlib.File):
                        compiler += [j.absolute_path(self.source_dir, self.build_dir)]
                    elif isinstance(j, str):
                        compiler += [j]
                    elif isinstance(j, (build.BuildTarget, build.CustomTarget)):
                        compiler += j.get_outputs()
                    else:
                        raise RuntimeError(f'Type "{type(j).__name__}" is not supported in get_introspection_data. This is a bug')

            return [{
                'language': 'unknown',
                'compiler': compiler,
                'parameters': [],
                'sources': source_list,
                'generated_sources': []
            }]

        return []

    def get_devenv(self) -> mesonlib.EnvironmentVariables:
        env = mesonlib.EnvironmentVariables()
        extra_paths = set()
        library_paths = set()
        build_machine = self.environment.machines[MachineChoice.BUILD]
        host_machine = self.environment.machines[MachineChoice.HOST]
        need_wine = not build_machine.is_windows() and host_machine.is_windows()
        for t in self.build.get_targets().values():
            in_default_dir = t.should_install() and not t.get_install_dir()[2]
            if t.for_machine != MachineChoice.HOST or not in_default_dir:
                continue
            tdir = os.path.join(self.environment.get_build_dir(), self.get_target_dir(t))
            if isinstance(t, build.Executable):
                # Add binaries that are going to be installed in bindir into PATH
                # so they get used by default instead of searching on system when
                # in developer environment.
                extra_paths.add(tdir)
                if host_machine.is_windows() or host_machine.is_cygwin():
                    # On windows we cannot rely on rpath to run executables from build
                    # directory. We have to add in PATH the location of every DLL needed.
                    library_paths.update(self.determine_windows_extra_paths(t, []))
            elif isinstance(t, build.SharedLibrary):
                # Add libraries that are going to be installed in libdir into
                # LD_LIBRARY_PATH. This allows running system applications using
                # that library.
                library_paths.add(tdir)
        if need_wine:
            # Executable paths should be in both PATH and WINEPATH.
            # - Having them in PATH makes bash completion find it,
            #   and make running "foo.exe" find it when wine-binfmt is installed.
            # - Having them in WINEPATH makes "wine foo.exe" find it.
            library_paths.update(extra_paths)
        if library_paths:
            if need_wine:
                env.prepend('WINEPATH', list(library_paths), separator=';')
            elif host_machine.is_windows() or host_machine.is_cygwin():
                extra_paths.update(library_paths)
            elif host_machine.is_darwin():
                env.prepend('DYLD_LIBRARY_PATH', list(library_paths))
            else:
                env.prepend('LD_LIBRARY_PATH', list(library_paths))
        if extra_paths:
            env.prepend('PATH', list(extra_paths))
        return env

    def compiler_to_generator_args(self, target: build.BuildTarget,
                                   compiler: 'Compiler', output: str = '@OUTPUT@',
                                   depfile: T.Union[str, None] = '@DEPFILE@',
                                   extras: T.Union[T.List[str], None] = None,
                                   input: str = '@INPUT@') -> CompilerArgs:
        '''
        The VS and Xcode backends need the full set of arguments for making a
        custom build rule. This is a convenience method to convert a Compiler
        to its arguments, for later concatenation.
        '''
        # FIXME: There are many other args missing
        commands = self.generate_basic_compiler_args(target, compiler)
        if depfile:
            commands += compiler.get_dependency_gen_args(output, depfile)
        commands += compiler.get_output_args(output)
        commands += self.get_source_dir_include_args(target, compiler)
        commands += self.get_build_dir_include_args(target, compiler)
        commands += compiler.get_compile_only_args()
        # Add per-target compile args, f.ex, `c_args : ['-DFOO']`. We set these
        # near the end since these are supposed to override everything else.
        commands += self.escape_extra_args(target.get_extra_args(compiler.get_language()))
        # Do not escape this one, it is interpreted by the build system
        # (Xcode considers these as variables to expand at build time)
        if extras is not None:
            commands += extras
        commands += [input]
        return commands

    def compiler_to_generator(self, target: build.BuildTarget,
                              compiler: 'Compiler',
                              sources: _ALL_SOURCES_TYPE,
                              output_templ: str,
                              depends: T.Optional[T.List[T.Union[build.BuildTarget, build.CustomTarget, build.CustomTargetIndex]]] = None,
                              ) -> build.GeneratedList:
        '''
        Some backends don't support custom compilers. This is a convenience
        method to convert a Compiler to a Generator.
        '''
        exelist = compiler.get_exelist()
        exe = programs.ExternalProgram(exelist[0])
        args = exelist[1:]
        commands = self.compiler_to_generator_args(target, compiler)
        generator = build.Generator(exe, args + commands.to_native(),
                                    [output_templ], depfile='@PLAINNAME@.d',
                                    depends=depends)
        return generator.process_files(sources, self.interpreter)

    def compile_target_to_generator(self, target: build.CompileTarget) -> build.GeneratedList:
        all_sources = T.cast('_ALL_SOURCES_TYPE', target.sources) + T.cast('_ALL_SOURCES_TYPE', target.generated)
        return self.compiler_to_generator(target, target.compiler, all_sources,
                                          target.output_templ, target.depends)
