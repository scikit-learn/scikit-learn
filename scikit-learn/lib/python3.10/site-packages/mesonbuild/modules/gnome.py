# SPDX-License-Identifier: Apache-2.0
# Copyright 2015-2016 The Meson development team
# Copyright © 2023-2024 Intel Corporation

'''This module provides helper functions for Gnome/GLib related
functionality such as gobject-introspection, gresources and gtk-doc'''
from __future__ import annotations

import copy
import itertools
import functools
import os
import subprocess
import textwrap
import typing as T

from . import (
    ExtensionModule, GirTarget, GResourceHeaderTarget, GResourceTarget, ModuleInfo,
    ModuleReturnValue, TypelibTarget, VapiTarget,
)
from .. import build
from .. import interpreter
from .. import mesonlib
from .. import mlog
from ..build import CustomTarget, CustomTargetIndex, Executable, GeneratedList, InvalidArguments
from ..dependencies import Dependency, InternalDependency
from ..dependencies.pkgconfig import PkgConfigDependency, PkgConfigInterface
from ..interpreter.type_checking import DEPENDS_KW, DEPEND_FILES_KW, ENV_KW, INSTALL_DIR_KW, INSTALL_KW, NoneType, DEPENDENCY_SOURCES_KW, in_set_validator
from ..interpreterbase import noPosargs, noKwargs, FeatureNew, FeatureDeprecated
from ..interpreterbase import typed_kwargs, KwargInfo, ContainerTypeInfo
from ..interpreterbase.decorators import typed_pos_args
from ..mesonlib import (
    MachineChoice, MesonException, OrderedSet, Popen_safe, join_args, quote_arg
)
from ..options import OptionKey
from ..programs import OverrideProgram
from ..scripts.gettext import read_linguas

if T.TYPE_CHECKING:
    from typing_extensions import Literal, TypeAlias, TypedDict

    from . import ModuleState
    from ..build import BuildTarget
    from ..compilers import Compiler
    from ..interpreter import Interpreter
    from ..interpreterbase import TYPE_var, TYPE_kwargs
    from ..mesonlib import FileOrString
    from ..programs import ExternalProgram

    class PostInstall(TypedDict):
        glib_compile_schemas: bool
        gio_querymodules: T.List[str]
        gtk_update_icon_cache: bool
        update_desktop_database: bool
        update_mime_database: bool

    class CompileSchemas(TypedDict):

        build_by_default: bool
        depend_files: T.List[FileOrString]

    class Yelp(TypedDict):

        languages: T.List[str]
        media: T.List[str]
        sources: T.List[str]
        symlink_media: bool

    class CompileResources(TypedDict):

        build_by_default: bool
        c_name: T.Optional[str]
        dependencies: T.List[T.Union[mesonlib.File, CustomTarget, CustomTargetIndex]]
        export: bool
        extra_args: T.List[str]
        gresource_bundle: bool
        install: bool
        install_dir: T.Optional[str]
        install_header: bool
        source_dir: T.List[str]

    class GenerateGir(TypedDict):

        build_by_default: bool
        dependencies: T.List[Dependency]
        doc_format: T.Optional[str]
        export_packages: T.List[str]
        extra_args: T.List[str]
        fatal_warnings: bool
        header: T.List[str]
        identifier_prefix: T.List[str]
        include_directories: T.List[T.Union[build.IncludeDirs, str]]
        includes: T.List[T.Union[str, GirTarget]]
        install: bool
        install_dir_gir: T.Optional[str]
        install_dir_typelib: T.Optional[str]
        link_with: T.List[T.Union[build.SharedLibrary, build.StaticLibrary]]
        namespace: str
        nsversion: str
        sources: T.List[T.Union[FileOrString, build.GeneratedTypes]]
        symbol_prefix: T.List[str]

    class GtkDoc(TypedDict):

        src_dir: T.List[T.Union[str, build.IncludeDirs]]
        main_sgml: str
        main_xml: str
        module_version: str
        namespace: str
        mode: Literal['xml', 'smgl', 'auto', 'none']
        html_args: T.List[str]
        scan_args: T.List[str]
        scanobjs_args: T.List[str]
        fixxref_args: T.List[str]
        mkdb_args: T.List[str]
        content_files: T.List[T.Union[build.GeneratedTypes, FileOrString]]
        ignore_headers: T.List[str]
        install_dir: T.List[str]
        check: bool
        install: bool
        gobject_typesfile: T.List[FileOrString]
        html_assets: T.List[FileOrString]
        expand_content_files: T.List[FileOrString]
        c_args: T.List[str]
        include_directories: T.List[T.Union[str, build.IncludeDirs]]
        dependencies: T.List[T.Union[Dependency, build.SharedLibrary, build.StaticLibrary]]

    class GdbusCodegen(TypedDict):

        sources: T.List[FileOrString]
        extra_args: T.List[str]
        interface_prefix: T.Optional[str]
        namespace: T.Optional[str]
        object_manager: bool
        build_by_default: bool
        annotations: T.List[T.List[str]]
        install_header: bool
        install_dir: T.Optional[str]
        docbook: T.Optional[str]
        rst: T.Optional[str]
        markdown: T.Optional[str]
        autocleanup: Literal['all', 'none', 'objects', 'default']

    class GenMarshal(TypedDict):

        build_always: T.Optional[str]
        build_always_stale: T.Optional[bool]
        build_by_default: T.Optional[bool]
        depend_files: T.List[mesonlib.File]
        extra_args: T.List[str]
        install_dir: T.Optional[str]
        install_header: bool
        internal: bool
        nostdinc: bool
        prefix: T.Optional[str]
        skip_source: bool
        sources: T.List[FileOrString]
        stdinc: bool
        valist_marshallers: bool

    class GenerateVapi(TypedDict):

        sources: T.List[T.Union[str, GirTarget]]
        install_dir: T.Optional[str]
        install: bool
        vapi_dirs: T.List[str]
        metadata_dirs: T.List[str]
        gir_dirs: T.List[str]
        packages: T.List[T.Union[str, InternalDependency]]

    class _MkEnumsCommon(TypedDict):

        install_header: bool
        install_dir: T.Optional[str]
        identifier_prefix: T.Optional[str]
        symbol_prefix: T.Optional[str]

    class MkEnumsSimple(_MkEnumsCommon):

        sources: T.List[FileOrString]
        header_prefix: str
        decorator: str
        function_prefix: str
        body_prefix: str

    class MkEnums(_MkEnumsCommon):

        sources: T.List[T.Union[FileOrString, build.GeneratedTypes]]
        c_template: T.Optional[FileOrString]
        h_template: T.Optional[FileOrString]
        comments: T.Optional[str]
        eprod: T.Optional[str]
        fhead: T.Optional[str]
        fprod: T.Optional[str]
        ftail: T.Optional[str]
        vhead: T.Optional[str]
        vprod: T.Optional[str]
        vtail: T.Optional[str]
        depends: T.List[T.Union[BuildTarget, CustomTarget, CustomTargetIndex]]

    ToolType: TypeAlias = T.Union[Executable, ExternalProgram, OverrideProgram]


# Differs from the CustomTarget version in that it straight defaults to True
_BUILD_BY_DEFAULT: KwargInfo[bool] = KwargInfo(
    'build_by_default', bool, default=True,
)

_EXTRA_ARGS_KW: KwargInfo[T.List[str]] = KwargInfo(
    'extra_args',
    ContainerTypeInfo(list, str),
    default=[],
    listify=True,
)

_MK_ENUMS_COMMON_KWS: T.List[KwargInfo] = [
    INSTALL_KW.evolve(name='install_header'),
    INSTALL_DIR_KW,
    KwargInfo('identifier_prefix', (str, NoneType)),
    KwargInfo('symbol_prefix', (str, NoneType)),
]

def annotations_validator(annotations: T.List[T.Union[str, T.List[str]]]) -> T.Optional[str]:
    """Validate gdbus-codegen annotations argument"""

    badlist = 'must be made up of 3 strings for ELEMENT, KEY, and VALUE'

    if not annotations:
        return None
    elif all(isinstance(annot, str) for annot in annotations):
        if len(annotations) == 3:
            return None
        else:
            return badlist
    elif not all(isinstance(annot, list) for annot in annotations):
        for c, annot in enumerate(annotations):
            if not isinstance(annot, list):
                return f'element {c+1} must be a list'
    else:
        for c, annot in enumerate(annotations):
            if len(annot) != 3 or not all(isinstance(i, str) for i in annot):
                return f'element {c+1} {badlist}'
    return None

# gresource compilation is broken due to the way
# the resource compiler and Ninja clash about it
#
# https://github.com/ninja-build/ninja/issues/1184
# https://bugzilla.gnome.org/show_bug.cgi?id=774368
gresource_dep_needed_version = '>= 2.51.1'

class GnomeModule(ExtensionModule):

    INFO = ModuleInfo('gnome')

    def __init__(self, interpreter: 'Interpreter') -> None:
        super().__init__(interpreter)
        self.gir_dep: T.Optional[Dependency] = None
        self.giscanner: T.Optional[ToolType] = None
        self.gicompiler: T.Optional[ToolType] = None
        self.install_glib_compile_schemas = False
        self.install_gio_querymodules: T.List[str] = []
        self.install_gtk_update_icon_cache = False
        self.install_update_desktop_database = False
        self.install_update_mime_database = False
        self.devenv: T.Optional[mesonlib.EnvironmentVariables] = None
        self.native_glib_version: T.Optional[str] = None
        self.methods.update({
            'post_install': self.post_install,
            'compile_resources': self.compile_resources,
            'generate_gir': self.generate_gir,
            'compile_schemas': self.compile_schemas,
            'yelp': self.yelp,
            'gtkdoc': self.gtkdoc,
            'gtkdoc_html_dir': self.gtkdoc_html_dir,
            'gdbus_codegen': self.gdbus_codegen,
            'mkenums': self.mkenums,
            'mkenums_simple': self.mkenums_simple,
            'genmarshal': self.genmarshal,
            'generate_vapi': self.generate_vapi,
        })

    def _get_native_glib_version(self, state: 'ModuleState') -> str:
        if self.native_glib_version is None:
            glib_dep = PkgConfigDependency('glib-2.0', state.environment,
                                           {'native': True, 'required': False})
            if glib_dep.found():
                self.native_glib_version = glib_dep.get_version()
            else:
                mlog.warning('Could not detect glib version, assuming 2.54. '
                             'You may get build errors if your glib is older.')
                self.native_glib_version = '2.54'
        return self.native_glib_version

    @mesonlib.run_once
    def __print_gresources_warning(self, state: 'ModuleState') -> None:
        if not mesonlib.version_compare(self._get_native_glib_version(state),
                                        gresource_dep_needed_version):
            mlog.warning('GLib compiled dependencies do not work reliably with \n'
                         'the current version of GLib. See the following upstream issue:',
                         mlog.bold('https://bugzilla.gnome.org/show_bug.cgi?id=774368'),
                         once=True, fatal=False)

    @staticmethod
    def _print_gdbus_warning() -> None:
        mlog.warning('Code generated with gdbus_codegen() requires the root directory be added to\n'
                     '  include_directories of targets with GLib < 2.51.3:',
                     mlog.bold('https://github.com/mesonbuild/meson/issues/1387'),
                     once=True, fatal=False)

    @staticmethod
    def _find_tool(state: 'ModuleState', tool: str) -> 'ToolType':
        tool_map = {
            'gio-querymodules': 'gio-2.0',
            'glib-compile-schemas': 'gio-2.0',
            'glib-compile-resources': 'gio-2.0',
            'gdbus-codegen': 'gio-2.0',
            'glib-genmarshal': 'glib-2.0',
            'glib-mkenums': 'glib-2.0',
            'g-ir-scanner': 'gobject-introspection-1.0',
            'g-ir-compiler': 'gobject-introspection-1.0',
        }
        depname = tool_map[tool]
        varname = tool.replace('-', '_')
        return state.find_tool(tool, depname, varname)

    @typed_kwargs(
        'gnome.post_install',
        KwargInfo('glib_compile_schemas', bool, default=False),
        KwargInfo('gio_querymodules', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo('gtk_update_icon_cache', bool, default=False),
        KwargInfo('update_desktop_database', bool, default=False, since='0.59.0'),
        KwargInfo('update_mime_database', bool, default=False, since='0.64.0'),
    )
    @noPosargs
    @FeatureNew('gnome.post_install', '0.57.0')
    def post_install(self, state: 'ModuleState', args: T.List['TYPE_var'], kwargs: 'PostInstall') -> ModuleReturnValue:
        rv: T.List['mesonlib.ExecutableSerialisation'] = []
        datadir_abs = os.path.join(state.environment.get_prefix(), state.environment.get_datadir())
        if kwargs['glib_compile_schemas'] and not self.install_glib_compile_schemas:
            self.install_glib_compile_schemas = True
            prog = self._find_tool(state, 'glib-compile-schemas')
            schemasdir = os.path.join(datadir_abs, 'glib-2.0', 'schemas')
            script = state.backend.get_executable_serialisation([prog, schemasdir])
            script.skip_if_destdir = True
            rv.append(script)
        for d in kwargs['gio_querymodules']:
            if d not in self.install_gio_querymodules:
                self.install_gio_querymodules.append(d)
                prog = self._find_tool(state, 'gio-querymodules')
                moduledir = os.path.join(state.environment.get_prefix(), d)
                script = state.backend.get_executable_serialisation([prog, moduledir])
                script.skip_if_destdir = True
                rv.append(script)
        if kwargs['gtk_update_icon_cache'] and not self.install_gtk_update_icon_cache:
            self.install_gtk_update_icon_cache = True
            prog = state.find_program('gtk4-update-icon-cache', required=False)
            found = isinstance(prog, Executable) or prog.found()
            if not found:
                prog = state.find_program('gtk-update-icon-cache')
            icondir = os.path.join(datadir_abs, 'icons', 'hicolor')
            script = state.backend.get_executable_serialisation([prog, '-q', '-t', '-f', icondir])
            script.skip_if_destdir = True
            rv.append(script)
        if kwargs['update_desktop_database'] and not self.install_update_desktop_database:
            self.install_update_desktop_database = True
            prog = state.find_program('update-desktop-database')
            appdir = os.path.join(datadir_abs, 'applications')
            script = state.backend.get_executable_serialisation([prog, '-q', appdir])
            script.skip_if_destdir = True
            rv.append(script)
        if kwargs['update_mime_database'] and not self.install_update_mime_database:
            self.install_update_mime_database = True
            prog = state.find_program('update-mime-database')
            appdir = os.path.join(datadir_abs, 'mime')
            script = state.backend.get_executable_serialisation([prog, appdir])
            script.skip_if_destdir = True
            rv.append(script)
        return ModuleReturnValue(None, rv)

    @typed_pos_args('gnome.compile_resources', str, (str, mesonlib.File, CustomTarget, CustomTargetIndex, GeneratedList))
    @typed_kwargs(
        'gnome.compile_resources',
        _BUILD_BY_DEFAULT,
        _EXTRA_ARGS_KW,
        INSTALL_KW,
        INSTALL_KW.evolve(name='install_header', since='0.37.0'),
        INSTALL_DIR_KW,
        KwargInfo('c_name', (str, NoneType)),
        KwargInfo('dependencies', ContainerTypeInfo(list, (mesonlib.File, CustomTarget, CustomTargetIndex)), default=[], listify=True),
        KwargInfo('export', bool, default=False, since='0.37.0'),
        KwargInfo('gresource_bundle', bool, default=False, since='0.37.0'),
        KwargInfo('source_dir', ContainerTypeInfo(list, str), default=[], listify=True),
    )
    def compile_resources(self, state: 'ModuleState', args: T.Tuple[str, 'FileOrString'],
                          kwargs: 'CompileResources') -> 'ModuleReturnValue':
        self.__print_gresources_warning(state)
        glib_version = self._get_native_glib_version(state)

        glib_compile_resources = self._find_tool(state, 'glib-compile-resources')
        cmd: T.List[T.Union['ToolType', str]] = [glib_compile_resources, '@INPUT@']

        source_dirs = kwargs['source_dir']
        dependencies = kwargs['dependencies']

        target_name, input_file = args

        # Validate dependencies
        subdirs: T.List[str] = []
        depends: T.List[T.Union[CustomTarget, CustomTargetIndex]] = []
        for dep in dependencies:
            if isinstance(dep, mesonlib.File):
                subdirs.append(dep.subdir)
            else:
                depends.append(dep)
                subdirs.append(dep.get_subdir())
                if not mesonlib.version_compare(glib_version, gresource_dep_needed_version):
                    m = 'The "dependencies" argument of gnome.compile_resources() cannot\n' \
                        'be used with the current version of glib-compile-resources due to\n' \
                        '<https://bugzilla.gnome.org/show_bug.cgi?id=774368>'
                    raise MesonException(m)

        if not mesonlib.version_compare(glib_version, gresource_dep_needed_version):
            # Resource xml files generated at build-time cannot be used with
            # gnome.compile_resources() because we need to scan the xml for
            # dependencies. Use configure_file() instead to generate it at
            # configure-time
            if isinstance(input_file, mesonlib.File):
                # glib-compile-resources will be run inside the source dir,
                # so we need either 'src_to_build' or the absolute path.
                # Absolute path is the easiest choice.
                if input_file.is_built:
                    ifile = os.path.join(state.environment.get_build_dir(), input_file.subdir, input_file.fname)
                else:
                    ifile = os.path.join(input_file.subdir, input_file.fname)

            elif isinstance(input_file, (CustomTarget, CustomTargetIndex, GeneratedList)):
                raise MesonException('Resource xml files generated at build-time cannot be used with '
                                     'gnome.compile_resources() in the current version of glib-compile-resources '
                                     'because we need to scan the xml for dependencies due to '
                                     '<https://bugzilla.gnome.org/show_bug.cgi?id=774368>\nUse '
                                     'configure_file() instead to generate it at configure-time.')
            else:
                ifile = os.path.join(state.subdir, input_file)

            depend_files, depends, subdirs = self._get_gresource_dependencies(
                state, ifile, source_dirs, dependencies)
        else:
            depend_files = []

        # Make source dirs relative to build dir now
        source_dirs = [os.path.join(state.build_to_src, state.subdir, d) for d in source_dirs]
        # Ensure build directories of generated deps are included
        source_dirs += subdirs
        # Always include current directory, but after paths set by user
        source_dirs.append(os.path.join(state.build_to_src, state.subdir))

        # Clean up duplicate directories
        source_dirs = list(OrderedSet(os.path.normpath(dir) for dir in source_dirs))

        for source_dir in source_dirs:
            cmd += ['--sourcedir', source_dir]

        if kwargs['c_name']:
            cmd += ['--c-name', kwargs['c_name']]
        if not kwargs['export']:
            cmd += ['--internal']

        cmd += ['--generate', '--target', '@OUTPUT@']
        cmd += kwargs['extra_args']

        gresource = kwargs['gresource_bundle']
        if gresource:
            output = f'{target_name}.gresource'
            name = f'{target_name}_gresource'
        else:
            if 'c' in state.environment.coredata.compilers.host:
                output = f'{target_name}.c'
                name = f'{target_name}_c'
            elif 'cpp' in state.environment.coredata.compilers.host:
                output = f'{target_name}.cpp'
                name = f'{target_name}_cpp'
            else:
                raise MesonException('Compiling GResources into code is only supported in C and C++ projects')

        if kwargs['install']:
            if not gresource:
                raise MesonException('The install kwarg only applies to gresource bundles, see install_header')
            elif not kwargs['install_dir']:
                raise MesonException('gnome.compile_resources: "install_dir" keyword argument must be set when "install" is true.')

        install_header = kwargs['install_header']
        if install_header and gresource:
            raise MesonException('The install_header kwarg does not apply to gresource bundles')
        if install_header and not kwargs['export']:
            raise MesonException('GResource header is installed yet export is not enabled')

        depfile: T.Optional[str] = None
        target_cmd: T.List[T.Union['ToolType', str]]
        if not mesonlib.version_compare(glib_version, gresource_dep_needed_version):
            # This will eventually go out of sync if dependencies are added
            target_cmd = cmd
        else:
            depfile = f'{output}.d'
            target_cmd = copy.copy(cmd) + ['--dependency-file', '@DEPFILE@']
        target_c = GResourceTarget(
            name,
            state.subdir,
            state.subproject,
            state.environment,
            target_cmd,
            [input_file],
            [output],
            build_by_default=kwargs['build_by_default'],
            depfile=depfile,
            depend_files=depend_files,
            extra_depends=depends,
            install=kwargs['install'],
            install_dir=[kwargs['install_dir']] if kwargs['install_dir'] else [],
            install_tag=['runtime'],
        )
        target_c.source_dirs = source_dirs

        if gresource: # Only one target for .gresource files
            return ModuleReturnValue(target_c, [target_c])

        install_dir = kwargs['install_dir'] or state.environment.coredata.optstore.get_value_for(OptionKey('includedir'))
        assert isinstance(install_dir, str), 'for mypy'
        target_h = GResourceHeaderTarget(
            f'{target_name}_h',
            state.subdir,
            state.subproject,
            state.environment,
            cmd,
            [input_file],
            [f'{target_name}.h'],
            build_by_default=kwargs['build_by_default'],
            extra_depends=depends,
            install=install_header,
            install_dir=[install_dir],
            install_tag=['devel'],
        )
        rv = [target_c, target_h]
        return ModuleReturnValue(rv, rv)

    @staticmethod
    def _get_gresource_dependencies(
            state: 'ModuleState', input_file: str, source_dirs: T.List[str],
            dependencies: T.Sequence[T.Union[mesonlib.File, CustomTarget, CustomTargetIndex]]
            ) -> T.Tuple[T.List[mesonlib.FileOrString], T.List[T.Union[CustomTarget, CustomTargetIndex]], T.List[str]]:

        cmd = ['glib-compile-resources',
               input_file,
               '--generate-dependencies']

        # Prefer generated files over source files
        cmd += ['--sourcedir', state.subdir] # Current build dir
        for source_dir in source_dirs:
            cmd += ['--sourcedir', os.path.join(state.subdir, source_dir)]

        try:
            pc, stdout, stderr = Popen_safe(cmd, cwd=state.environment.get_source_dir())
        except (FileNotFoundError, PermissionError):
            raise MesonException('Could not execute glib-compile-resources.')
        if pc.returncode != 0:
            m = f'glib-compile-resources failed to get dependencies for {cmd[1]}:\n{stderr}'
            mlog.warning(m)
            raise subprocess.CalledProcessError(pc.returncode, cmd)

        raw_dep_files: T.List[str] = stdout.split('\n')[:-1]

        depends: T.List[T.Union[CustomTarget, CustomTargetIndex]] = []
        subdirs: T.List[str] = []
        dep_files: T.List[mesonlib.FileOrString] = []
        for resfile in raw_dep_files.copy():
            resbasename = os.path.basename(resfile)
            for dep in dependencies:
                if isinstance(dep, mesonlib.File):
                    if dep.fname != resbasename:
                        continue
                    raw_dep_files.remove(resfile)
                    dep_files.append(dep)
                    subdirs.append(dep.subdir)
                    break
                elif isinstance(dep, (CustomTarget, CustomTargetIndex)):
                    fname = None
                    outputs = {(o, os.path.basename(o)) for o in dep.get_outputs()}
                    for o, baseo in outputs:
                        if baseo == resbasename:
                            fname = o
                            break
                    if fname is not None:
                        raw_dep_files.remove(resfile)
                        depends.append(dep)
                        subdirs.append(dep.get_subdir())
                        break
            else:
                # In generate-dependencies mode, glib-compile-resources doesn't raise
                # an error for missing resources but instead prints whatever filename
                # was listed in the input file.  That's good because it means we can
                # handle resource files that get generated as part of the build, as
                # follows.
                #
                # If there are multiple generated resource files with the same basename
                # then this code will get confused.
                try:
                    f = mesonlib.File.from_source_file(state.environment.get_source_dir(),
                                                       ".", resfile)
                except MesonException:
                    raise MesonException(
                        f'Resource "{resfile}" listed in "{input_file}" was not found. '
                        'If this is a generated file, pass the target that generates '
                        'it to gnome.compile_resources() using the "dependencies" '
                        'keyword argument.')
                raw_dep_files.remove(resfile)
                dep_files.append(f)
        dep_files.extend(raw_dep_files)
        return dep_files, depends, subdirs

    def _get_link_args(self, state: 'ModuleState',
                       lib: T.Union[build.SharedLibrary, build.StaticLibrary],
                       depends: T.Sequence[T.Union[build.BuildTarget, 'build.GeneratedTypes', 'FileOrString', build.StructuredSources]],
                       include_rpath: bool = False,
                       use_gir_args: bool = False
                       ) -> T.Tuple[T.List[str], T.List[T.Union[build.BuildTarget, 'build.GeneratedTypes', 'FileOrString', build.StructuredSources]]]:
        link_command: T.List[str] = []
        new_depends = list(depends)
        # Construct link args
        if isinstance(lib, build.SharedLibrary):
            libdir = os.path.join(state.environment.get_build_dir(), state.backend.get_target_dir(lib))
            link_command.append('-L' + libdir)
            if include_rpath:
                link_command.append('-Wl,-rpath,' + libdir)
            new_depends.append(lib)
            # Needed for the following binutils bug:
            # https://github.com/mesonbuild/meson/issues/1911
            # However, g-ir-scanner does not understand -Wl,-rpath
            # so we need to use -L instead
            for d in lib.determine_rpath_dirs():
                d = os.path.join(state.environment.get_build_dir(), d)
                link_command.append('-L' + d)
                if include_rpath:
                    link_command.append('-Wl,-rpath,' + d)
        if use_gir_args and self._gir_has_option('--extra-library'):
            link_command.append('--extra-library=' + lib.name)
        else:
            link_command.append('-l' + lib.name)
        return link_command, new_depends

    def _get_dependencies_flags_raw(
            self, deps: T.Sequence[T.Union['Dependency', build.BuildTarget, CustomTarget, CustomTargetIndex]],
            state: 'ModuleState',
            depends: T.Sequence[T.Union[build.BuildTarget, 'build.GeneratedTypes', 'FileOrString', build.StructuredSources]],
            include_rpath: bool,
            use_gir_args: bool,
            ) -> T.Tuple[OrderedSet[str], OrderedSet[T.Union[str, T.Tuple[str, str]]], OrderedSet[T.Union[str, T.Tuple[str, str]]], OrderedSet[str],
                         T.List[T.Union[build.BuildTarget, 'build.GeneratedTypes', 'FileOrString', build.StructuredSources]]]:
        cflags: OrderedSet[str] = OrderedSet()
        # External linker flags that can't be de-duped reliably because they
        # require two args in order, such as -framework AVFoundation will be stored as a tuple.
        internal_ldflags: OrderedSet[T.Union[str, T.Tuple[str, str]]] = OrderedSet()
        external_ldflags: OrderedSet[T.Union[str, T.Tuple[str, str]]] = OrderedSet()
        gi_includes: OrderedSet[str] = OrderedSet()
        deps = mesonlib.listify(deps)
        depends = list(depends)

        for dep in deps:
            if isinstance(dep, Dependency):
                girdir = dep.get_variable(pkgconfig='girdir', internal='girdir', default_value='')
                if girdir:
                    assert isinstance(girdir, str), 'for mypy'
                    gi_includes.update([girdir])
            if isinstance(dep, InternalDependency):
                cflags.update(dep.get_compile_args())
                cflags.update(state.get_include_args(dep.include_directories))
                for lib in dep.libraries:
                    if isinstance(lib, build.SharedLibrary):
                        _ld, depends = self._get_link_args(state, lib, depends, include_rpath)
                        internal_ldflags.update(_ld)
                        libdepflags = self._get_dependencies_flags_raw(lib.get_external_deps(), state, depends, include_rpath,
                                                                       use_gir_args)
                        cflags.update(libdepflags[0])
                        internal_ldflags.update(libdepflags[1])
                        external_ldflags.update(libdepflags[2])
                        gi_includes.update(libdepflags[3])
                        depends = libdepflags[4]
                extdepflags = self._get_dependencies_flags_raw(dep.ext_deps, state, depends, include_rpath,
                                                               use_gir_args)
                cflags.update(extdepflags[0])
                internal_ldflags.update(extdepflags[1])
                external_ldflags.update(extdepflags[2])
                gi_includes.update(extdepflags[3])
                depends = extdepflags[4]
                for source in dep.sources:
                    if isinstance(source, GirTarget):
                        gi_includes.update([os.path.join(state.environment.get_build_dir(),
                                            source.get_subdir())])
            # This should be any dependency other than an internal one.
            elif isinstance(dep, Dependency):
                cflags.update(dep.get_compile_args())
                ldflags = iter(dep.get_link_args(raw=True))
                for flag in ldflags:
                    if (os.path.isabs(flag) and
                            # For PkgConfigDependency only:
                            getattr(dep, 'is_libtool', False)):
                        lib_dir = os.path.dirname(flag)
                        external_ldflags.update([f'-L{lib_dir}'])
                        if include_rpath:
                            external_ldflags.update([f'-Wl,-rpath,{lib_dir}'])
                        libname = os.path.basename(flag)
                        if libname.startswith("lib"):
                            libname = libname[3:]
                        libname = libname.split(".so")[0]
                        flag = f"-l{libname}"
                    # FIXME: Hack to avoid passing some compiler options in
                    if flag.startswith("-W") and not flag.startswith('-Wl,-rpath,'):
                        continue
                    # If it's a framework arg, slurp the framework name too
                    # to preserve the order of arguments
                    if flag == '-framework':
                        external_ldflags.update([(flag, next(ldflags))])
                    else:
                        external_ldflags.update([flag])
            elif isinstance(dep, (build.StaticLibrary, build.SharedLibrary)):
                cflags.update(state.get_include_args(dep.get_include_dirs()))
                depends.append(dep)
            else:
                mlog.log(f'dependency {dep!r} not handled to build gir files')
                continue

        if use_gir_args and self._gir_has_option('--extra-library'):
            def fix_ldflags(ldflags: T.Iterable[T.Union[str, T.Tuple[str, str]]]) -> OrderedSet[T.Union[str, T.Tuple[str, str]]]:
                fixed_ldflags: OrderedSet[T.Union[str, T.Tuple[str, str]]] = OrderedSet()
                for ldflag in ldflags:
                    if isinstance(ldflag, str) and ldflag.startswith("-l"):
                        ldflag = ldflag.replace('-l', '--extra-library=', 1)
                    fixed_ldflags.add(ldflag)
                return fixed_ldflags
            internal_ldflags = fix_ldflags(internal_ldflags)
            external_ldflags = fix_ldflags(external_ldflags)
        return cflags, internal_ldflags, external_ldflags, gi_includes, depends

    def _get_dependencies_flags(
            self, deps: T.Sequence[T.Union['Dependency', build.BuildTarget, CustomTarget, CustomTargetIndex]],
            state: 'ModuleState',
            depends: T.Sequence[T.Union[build.BuildTarget, 'build.GeneratedTypes', 'FileOrString', build.StructuredSources]],
            include_rpath: bool = False,
            use_gir_args: bool = False,
            ) -> T.Tuple[OrderedSet[str], T.List[str], T.List[str], OrderedSet[str],
                         T.List[T.Union[build.BuildTarget, 'build.GeneratedTypes', 'FileOrString', build.StructuredSources]]]:

        cflags, internal_ldflags_raw, external_ldflags_raw, gi_includes, depends = self._get_dependencies_flags_raw(deps, state, depends, include_rpath, use_gir_args)
        internal_ldflags: T.List[str] = []
        external_ldflags: T.List[str] = []

        # Extract non-deduplicable argument groups out of the tuples.
        for ldflag in internal_ldflags_raw:
            if isinstance(ldflag, str):
                internal_ldflags.append(ldflag)
            else:
                internal_ldflags.extend(ldflag)
        for ldflag in external_ldflags_raw:
            if isinstance(ldflag, str):
                external_ldflags.append(ldflag)
            else:
                external_ldflags.extend(ldflag)

        return cflags, internal_ldflags, external_ldflags, gi_includes, depends

    def _unwrap_gir_target(self, girtarget: T.Union[Executable, build.StaticLibrary, build.SharedLibrary], state: 'ModuleState'
                           ) -> T.Union[Executable, build.StaticLibrary, build.SharedLibrary]:
        if not isinstance(girtarget, (Executable, build.SharedLibrary,
                                      build.StaticLibrary)):
            raise MesonException(f'Gir target must be an executable or library but is "{girtarget}" of type {type(girtarget).__name__}')

        STATIC_BUILD_REQUIRED_VERSION = ">=1.58.1"
        if isinstance(girtarget, (build.StaticLibrary)) and \
           not self._giscanner_version_compare(state, STATIC_BUILD_REQUIRED_VERSION):
            raise MesonException('Static libraries can only be introspected with GObject-Introspection ' + STATIC_BUILD_REQUIRED_VERSION)

        return girtarget

    def _devenv_prepend(self, varname: str, value: str) -> None:
        if self.devenv is None:
            self.devenv = mesonlib.EnvironmentVariables()
        self.devenv.prepend(varname, [value])

    def postconf_hook(self, b: build.Build) -> None:
        if self.devenv is not None:
            b.devenv.append(self.devenv)

    def _get_gir_dep(self, state: 'ModuleState') -> T.Tuple[Dependency, T.Union[Executable, 'ExternalProgram', 'OverrideProgram'],
                                                            T.Union[Executable, 'ExternalProgram', 'OverrideProgram']]:
        if not self.gir_dep:
            self.gir_dep = state.dependency('gobject-introspection-1.0')
            self.giscanner = self._find_tool(state, 'g-ir-scanner')
            self.gicompiler = self._find_tool(state, 'g-ir-compiler')
        return self.gir_dep, self.giscanner, self.gicompiler

    def _giscanner_version_compare(self, state: 'ModuleState', cmp: str) -> bool:
        # Support for --version was introduced in g-i 1.58, but Ubuntu
        # Bionic shipped 1.56.1. As all our version checks are greater
        # than 1.58, we can just return False if get_version fails.
        try:
            giscanner, _, _ = self._get_gir_dep(state)
            return mesonlib.version_compare(giscanner.get_version(), cmp)
        except MesonException:
            return False

    @functools.lru_cache(maxsize=None)
    def _gir_has_option(self, option: str) -> bool:
        exe = self.giscanner
        if isinstance(exe, (Executable, OverrideProgram)):
            # Handle overridden g-ir-scanner
            assert option in {'--extra-library', '--sources-top-dirs'}
            return True
        p, o, _ = Popen_safe(exe.get_command() + ['--help'], stderr=subprocess.STDOUT)
        return p.returncode == 0 and option in o

    # May mutate depends and gir_inc_dirs
    @staticmethod
    def _scan_include(state: 'ModuleState', includes: T.List[T.Union[str, GirTarget]]
                      ) -> T.Tuple[T.List[str], T.List[str], T.List[GirTarget]]:
        ret: T.List[str] = []
        gir_inc_dirs: T.List[str] = []
        depends: T.List[GirTarget] = []

        for inc in includes:
            if isinstance(inc, str):
                ret += [f'--include={inc}']
            elif isinstance(inc, GirTarget):
                gir_inc_dirs .append(os.path.join(state.environment.get_build_dir(), inc.get_subdir()))
                ret.append(f"--include-uninstalled={os.path.join(inc.get_subdir(), inc.get_basename())}")
                depends.append(inc)

        return ret, gir_inc_dirs, depends

    @staticmethod
    def _scan_langs(state: 'ModuleState', langs: T.Iterable[str]) -> T.List[str]:
        ret: T.List[str] = []

        for lang in langs:
            link_args = state.environment.coredata.get_external_link_args(MachineChoice.HOST, lang)
            for link_arg in link_args:
                if link_arg.startswith('-L'):
                    ret.append(link_arg)

        return ret

    @staticmethod
    def _scan_gir_targets(state: 'ModuleState', girtargets: T.Sequence[build.BuildTarget]) -> T.List[T.Union[str, Executable]]:
        ret: T.List[T.Union[str, Executable]] = []

        for girtarget in girtargets:
            if isinstance(girtarget, Executable):
                ret += ['--program', girtarget]
            else:
                # Because of https://gitlab.gnome.org/GNOME/gobject-introspection/merge_requests/72
                # we can't use the full path until this is merged.
                libpath = os.path.join(girtarget.get_subdir(), girtarget.get_filename())
                # Must use absolute paths here because g-ir-scanner will not
                # add them to the runtime path list if they're relative. This
                # means we cannot use @BUILD_ROOT@
                build_root = state.environment.get_build_dir()
                if isinstance(girtarget, build.SharedLibrary):
                    # need to put our output directory first as we need to use the
                    # generated libraries instead of any possibly installed system/prefix
                    # ones.
                    ret += ["-L{}/{}".format(build_root, os.path.dirname(libpath))]
                    libname = girtarget.get_basename()
                else:
                    libname = os.path.join(f"{build_root}/{libpath}")
                ret += ['--library', libname]
                # Needed for the following binutils bug:
                # https://github.com/mesonbuild/meson/issues/1911
                # However, g-ir-scanner does not understand -Wl,-rpath
                # so we need to use -L instead
                for d in girtarget.determine_rpath_dirs():
                    d = os.path.join(state.environment.get_build_dir(), d)
                    ret.append('-L' + d)

        return ret

    @staticmethod
    def _get_girtargets_langs_compilers(girtargets: T.Sequence[build.BuildTarget]) -> T.List[T.Tuple[str, 'Compiler']]:
        ret: T.List[T.Tuple[str, 'Compiler']] = []
        for girtarget in girtargets:
            for lang, compiler in girtarget.compilers.items():
                # XXX: Can you use g-i with any other language?
                if lang in {'c', 'cpp', 'objc', 'objcpp', 'd'}:
                    ret.append((lang, compiler))
                    break

        return ret

    @staticmethod
    def _get_gir_targets_deps(girtargets: T.Sequence[build.BuildTarget]
                              ) -> T.List[T.Union[build.BuildTarget, CustomTarget, CustomTargetIndex, Dependency]]:
        ret: T.List[T.Union[build.BuildTarget, CustomTarget, CustomTargetIndex, Dependency]] = []
        for girtarget in girtargets:
            ret += girtarget.get_all_link_deps()
            ret += girtarget.get_external_deps()
        return ret

    @staticmethod
    def _get_gir_targets_inc_dirs(girtargets: T.Sequence[build.BuildTarget]) -> OrderedSet[build.IncludeDirs]:
        ret: OrderedSet = OrderedSet()
        for girtarget in girtargets:
            ret.update(girtarget.get_include_dirs())
        return ret

    @staticmethod
    def _get_langs_compilers_flags(state: 'ModuleState', langs_compilers: T.List[T.Tuple[str, 'Compiler']]
                                   ) -> T.Tuple[T.List[str], T.List[str], T.List[str]]:
        cflags: T.List[str] = []
        internal_ldflags: T.List[str] = []
        external_ldflags: T.List[str] = []

        for lang, compiler in langs_compilers:
            if state.global_args.get(lang):
                cflags += state.global_args[lang]
            if state.project_args.get(lang):
                cflags += state.project_args[lang]
            if OptionKey('b_sanitize') in compiler.base_options:
                sanitize = state.environment.coredata.optstore.get_value('b_sanitize')
                assert isinstance(sanitize, list)
                cflags += compiler.sanitizer_compile_args(sanitize)
                # These must be first in ldflags
                if 'address' in sanitize:
                    internal_ldflags += ['-lasan']
                if 'thread' in sanitize:
                    internal_ldflags += ['-ltsan']
                if 'undefined' in sanitize:
                    internal_ldflags += ['-lubsan']
                # FIXME: Linking directly to lib*san is not recommended but g-ir-scanner
                # does not understand -f LDFLAGS. https://bugzilla.gnome.org/show_bug.cgi?id=783892
                # ldflags += compiler.sanitizer_link_args(sanitize)

        return cflags, internal_ldflags, external_ldflags

    @staticmethod
    def _make_gir_filelist(state: 'ModuleState', srcdir: str, ns: str,
                           nsversion: str, girtargets: T.Sequence[build.BuildTarget],
                           libsources: T.Sequence[T.Union[
                               str, mesonlib.File, GeneratedList,
                               CustomTarget, CustomTargetIndex]]
                           ) -> str:
        gir_filelist_dir = state.backend.get_target_private_dir_abs(girtargets[0])
        if not os.path.isdir(gir_filelist_dir):
            os.mkdir(gir_filelist_dir)
        gir_filelist_filename = os.path.join(gir_filelist_dir, f'{ns}_{nsversion}_gir_filelist')

        with open(gir_filelist_filename, 'w', encoding='utf-8') as gir_filelist:
            for s in libsources:
                if isinstance(s, (CustomTarget, CustomTargetIndex)):
                    for custom_output in s.get_outputs():
                        gir_filelist.write(os.path.join(state.environment.get_build_dir(),
                                                        state.backend.get_target_dir(s),
                                                        custom_output) + '\n')
                elif isinstance(s, mesonlib.File):
                    gir_filelist.write(s.rel_to_builddir(state.build_to_src) + '\n')
                elif isinstance(s, GeneratedList):
                    for gen_src in s.get_outputs():
                        gir_filelist.write(os.path.join(srcdir, gen_src) + '\n')
                else:
                    gir_filelist.write(os.path.join(srcdir, s) + '\n')

        return gir_filelist_filename

    def _make_gir_target(
            self,
            state: 'ModuleState',
            girfile: str,
            scan_command: T.Sequence[T.Union['FileOrString', Executable, ExternalProgram, OverrideProgram]],
            generated_files: T.Sequence[T.Union[str, mesonlib.File, CustomTarget, CustomTargetIndex, GeneratedList]],
            depends: T.Sequence[T.Union['FileOrString', build.BuildTarget, 'build.GeneratedTypes', build.StructuredSources]],
            env_flags: T.Sequence[str],
            kwargs: T.Dict[str, T.Any]) -> GirTarget:
        install = kwargs['install_gir']
        if install is None:
            install = kwargs['install']

        install_dir = kwargs['install_dir_gir']
        if install_dir is None:
            install_dir = os.path.join(state.environment.get_datadir(), 'gir-1.0')
        elif install_dir is False:
            install = False

        # g-ir-scanner uses pkg-config to find libraries such as glib. They could
        # be built as subproject in which case we need to trick it to use
        # -uninstalled.pc files Meson generated. It also must respect pkgconfig
        # settings user could have set in machine file, like PKG_CONFIG_LIBDIR,
        # SYSROOT, etc.
        run_env = PkgConfigInterface.get_env(state.environment, MachineChoice.HOST, uninstalled=True)
        # g-ir-scanner uses Python's distutils to find the compiler, which uses 'CC'
        cc_exelist = state.environment.coredata.compilers.host['c'].get_exelist()
        run_env.set('CC', [quote_arg(x) for x in cc_exelist], ' ')
        run_env.set('CFLAGS', [quote_arg(x) for x in env_flags], ' ')
        run_env.merge(kwargs['env'])

        # response file supported?
        rspable = self._giscanner_version_compare(state, '>= 1.85.0')

        return GirTarget(
            girfile,
            state.subdir,
            state.subproject,
            state.environment,
            scan_command,
            generated_files,
            [girfile],
            build_by_default=kwargs['build_by_default'],
            extra_depends=depends,
            install=install,
            install_dir=[install_dir],
            install_tag=['devel'],
            env=run_env,
            rspable=rspable,
        )

    @staticmethod
    def _make_typelib_target(state: 'ModuleState', typelib_output: str,
                             typelib_cmd: T.Sequence[T.Union[str, Executable, ExternalProgram, CustomTarget]],
                             generated_files: T.Sequence[T.Union[str, mesonlib.File, CustomTarget, CustomTargetIndex, GeneratedList]],
                             kwargs: T.Dict[str, T.Any]) -> TypelibTarget:
        install = kwargs['install_typelib']
        if install is None:
            install = kwargs['install']

        install_dir = kwargs['install_dir_typelib']
        if install_dir is None:
            install_dir = os.path.join(state.environment.get_libdir(), 'girepository-1.0')
        elif install_dir is False:
            install = False

        return TypelibTarget(
            typelib_output,
            state.subdir,
            state.subproject,
            state.environment,
            typelib_cmd,
            generated_files,
            [typelib_output],
            install=install,
            install_dir=[install_dir],
            install_tag=['typelib'],
            build_by_default=kwargs['build_by_default'],
            env=kwargs['env'],
        )

    @staticmethod
    def _gather_typelib_includes_and_update_depends(
            state: 'ModuleState',
            deps: T.Sequence[T.Union[Dependency, build.BuildTarget, CustomTarget, CustomTargetIndex]],
            depends: T.Sequence[T.Union[build.BuildTarget, 'build.GeneratedTypes', 'FileOrString', build.StructuredSources]]
            ) -> T.Tuple[T.List[str], T.List[T.Union[build.BuildTarget, 'build.GeneratedTypes', 'FileOrString', build.StructuredSources]]]:
        # Need to recursively add deps on GirTarget sources from our
        # dependencies and also find the include directories needed for the
        # typelib generation custom target below.
        typelib_includes: T.List[str] = []
        new_depends = list(depends)
        for dep in deps:
            # Add a dependency on each GirTarget listed in dependencies and add
            # the directory where it will be generated to the typelib includes
            if isinstance(dep, InternalDependency):
                for source in dep.sources:
                    if isinstance(source, GirTarget) and source not in depends:
                        new_depends.append(source)
                        subdir = os.path.join(state.environment.get_build_dir(),
                                              source.get_subdir())
                        if subdir not in typelib_includes:
                            typelib_includes.append(subdir)
            # Do the same, but for dependencies of dependencies. These are
            # stored in the list of generated sources for each link dep (from
            # girtarget.get_all_link_deps() above).
            # FIXME: Store this in the original form from declare_dependency()
            # so it can be used here directly.
            elif isinstance(dep, build.SharedLibrary):
                for g_source in dep.generated:
                    if isinstance(g_source, GirTarget):
                        subdir = os.path.join(state.environment.get_build_dir(),
                                              g_source.get_subdir())
                        if subdir not in typelib_includes:
                            typelib_includes.append(subdir)
            if isinstance(dep, Dependency):
                girdir = dep.get_variable(pkgconfig='girdir', internal='girdir', default_value='')
                assert isinstance(girdir, str), 'for mypy'
                if girdir and girdir not in typelib_includes:
                    typelib_includes.append(girdir)
        return typelib_includes, new_depends

    @staticmethod
    def _get_external_args_for_langs(state: 'ModuleState', langs: T.List[str]) -> T.List[str]:
        ret: T.List[str] = []
        for lang in langs:
            ret += mesonlib.listify(state.environment.coredata.get_external_args(MachineChoice.HOST, lang))
        return ret

    @staticmethod
    def _get_scanner_cflags(cflags: T.Iterable[str]) -> T.Iterable[str]:
        'g-ir-scanner only accepts -I/-D/-U; must ignore all other flags'
        for f in cflags:
            # _FORTIFY_SOURCE depends on / works together with -O, on the other hand this
            # just invokes the preprocessor anyway
            if f.startswith(('-D', '-U', '-I')) and not f.startswith('-D_FORTIFY_SOURCE'):
                yield f

    @staticmethod
    def _get_scanner_ldflags(ldflags: T.Iterable[str]) -> tuple[list[str], list[str]]:
        'g-ir-scanner only accepts -L/-l; must ignore -F and other linker flags'
        return (
            [f for f in ldflags if f.startswith(('-L', '-l', '--extra-library'))],
            [f for f in ldflags if f.startswith(('-Wl,-rpath,'))],
        )

    @typed_pos_args('gnome.generate_gir', varargs=(Executable, build.SharedLibrary, build.StaticLibrary), min_varargs=1)
    @typed_kwargs(
        'gnome.generate_gir',
        INSTALL_KW,
        _BUILD_BY_DEFAULT.evolve(since='0.40.0'),
        _EXTRA_ARGS_KW,
        ENV_KW.evolve(since='1.2.0'),
        KwargInfo('dependencies', ContainerTypeInfo(list, Dependency), default=[], listify=True),
        KwargInfo('doc_format', (str, NoneType), since='1.8.0'),
        KwargInfo('export_packages', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo('fatal_warnings', bool, default=False, since='0.55.0'),
        KwargInfo('header', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo('identifier_prefix', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo('include_directories', ContainerTypeInfo(list, (str, build.IncludeDirs)), default=[], listify=True),
        KwargInfo('includes', ContainerTypeInfo(list, (str, GirTarget)), default=[], listify=True),
        KwargInfo('install_gir', (bool, NoneType), since='0.61.0'),
        KwargInfo('install_dir_gir', (str, bool, NoneType),
                  deprecated_values={False: ('0.61.0', 'Use install_gir to disable installation')},
                  validator=lambda x: 'as boolean can only be false' if x is True else None),
        KwargInfo('install_typelib', (bool, NoneType), since='0.61.0'),
        KwargInfo('install_dir_typelib', (str, bool, NoneType),
                  deprecated_values={False: ('0.61.0', 'Use install_typelib to disable installation')},
                  validator=lambda x: 'as boolean can only be false' if x is True else None),
        KwargInfo('link_with', ContainerTypeInfo(list, (build.SharedLibrary, build.StaticLibrary)), default=[], listify=True),
        KwargInfo('namespace', str, required=True),
        KwargInfo('nsversion', str, required=True),
        KwargInfo('sources', ContainerTypeInfo(list, (str, mesonlib.File, GeneratedList, CustomTarget, CustomTargetIndex)), default=[], listify=True),
        KwargInfo('symbol_prefix', ContainerTypeInfo(list, str), default=[], listify=True),
    )
    def generate_gir(self, state: 'ModuleState', args: T.Tuple[T.List[T.Union[Executable, build.SharedLibrary, build.StaticLibrary]]],
                     kwargs: 'GenerateGir') -> ModuleReturnValue:
        # Ensure we have a C compiler even in C++ projects.
        state.add_language('c', MachineChoice.HOST)

        girtargets = [self._unwrap_gir_target(arg, state) for arg in args[0]]
        if len(girtargets) > 1 and any(isinstance(el, Executable) for el in girtargets):
            raise MesonException('generate_gir only accepts a single argument when one of the arguments is an executable')

        gir_dep, giscanner, gicompiler = self._get_gir_dep(state)

        ns = kwargs['namespace']
        nsversion = kwargs['nsversion']
        libsources = kwargs['sources']

        girfile = f'{ns}-{nsversion}.gir'
        srcdir = os.path.join(state.environment.get_source_dir(), state.subdir)
        builddir = os.path.join(state.environment.get_build_dir(), state.subdir)

        depends: T.List[T.Union['FileOrString', 'build.GeneratedTypes', build.BuildTarget, build.StructuredSources]] = []
        depends.extend(gir_dep.sources)
        depends.extend(girtargets)

        langs_compilers = self._get_girtargets_langs_compilers(girtargets)
        cflags, internal_ldflags, external_ldflags = self._get_langs_compilers_flags(state, langs_compilers)
        deps = self._get_gir_targets_deps(girtargets)
        deps += kwargs['dependencies']
        deps += [gir_dep]
        typelib_includes, depends = self._gather_typelib_includes_and_update_depends(state, deps, depends)
        # ldflags will be misinterpreted by gir scanner (showing
        # spurious dependencies) but building GStreamer fails if they
        # are not used here.
        dep_cflags, dep_internal_ldflags, dep_external_ldflags, gi_includes, depends = \
            self._get_dependencies_flags(deps, state, depends, use_gir_args=True)
        scan_cflags = []
        scan_cflags += list(self._get_scanner_cflags(cflags))
        scan_cflags += list(self._get_scanner_cflags(dep_cflags))
        scan_cflags += list(self._get_scanner_cflags(self._get_external_args_for_langs(state, [lc[0] for lc in langs_compilers])))
        scan_internal_ldflags = []
        scan_external_ldflags = []
        scan_env_ldflags = state.environment.coredata.get_external_link_args(MachineChoice.HOST, 'c')
        for cli_flags, env_flags in (self._get_scanner_ldflags(internal_ldflags), self._get_scanner_ldflags(dep_internal_ldflags)):
            scan_internal_ldflags += cli_flags
            scan_env_ldflags += env_flags
        for cli_flags, env_flags in (self._get_scanner_ldflags(external_ldflags), self._get_scanner_ldflags(dep_external_ldflags)):
            scan_external_ldflags += cli_flags
            scan_env_ldflags += env_flags
        girtargets_inc_dirs = self._get_gir_targets_inc_dirs(girtargets)
        inc_dirs = kwargs['include_directories']

        gir_inc_dirs: T.List[str] = []

        scan_command: T.List[T.Union[str, Executable, 'ExternalProgram', 'OverrideProgram']] = [giscanner]
        scan_command += ['--quiet']
        scan_command += ['--no-libtool']
        scan_command += ['--namespace=' + ns, '--nsversion=' + nsversion]
        scan_command += ['--warn-all']
        scan_command += ['--output', '@OUTPUT@']
        scan_command += [f'--c-include={h}' for h in kwargs['header']]
        scan_command += kwargs['extra_args']
        scan_command += ['-I' + srcdir, '-I' + builddir]
        scan_command += state.get_include_args(girtargets_inc_dirs)
        scan_command += ['--filelist=' + self._make_gir_filelist(state, srcdir, ns, nsversion, girtargets, libsources)]
        for l in kwargs['link_with']:
            _cflags, depends = self._get_link_args(state, l, depends, use_gir_args=True)
            scan_command.extend(_cflags)
        _cmd, _ginc, _deps = self._scan_include(state, kwargs['includes'])
        scan_command.extend(_cmd)
        gir_inc_dirs.extend(_ginc)
        depends.extend(_deps)

        scan_command += [f'--symbol-prefix={p}' for p in kwargs['symbol_prefix']]
        scan_command += [f'--identifier-prefix={p}' for p in kwargs['identifier_prefix']]
        scan_command += [f'--pkg-export={p}' for p in kwargs['export_packages']]
        scan_command += ['--cflags-begin']
        scan_command += scan_cflags
        scan_command += ['--cflags-end']
        scan_command += state.get_include_args(inc_dirs)
        scan_command += state.get_include_args(itertools.chain(gi_includes, gir_inc_dirs, inc_dirs), prefix='--add-include-path=')
        scan_command += list(scan_internal_ldflags)
        scan_command += self._scan_gir_targets(state, girtargets)
        scan_command += self._scan_langs(state, [lc[0] for lc in langs_compilers])
        scan_command += list(scan_external_ldflags)

        if self._gir_has_option('--sources-top-dirs'):
            scan_command += ['--sources-top-dirs', os.path.join(state.environment.get_source_dir(), state.root_subdir)]
            scan_command += ['--sources-top-dirs', os.path.join(state.environment.get_build_dir(), state.root_subdir)]

        if kwargs['doc_format'] is not None and self._gir_has_option('--doc-format'):
            scan_command += ['--doc-format', kwargs['doc_format']]

        if '--warn-error' in scan_command:
            FeatureDeprecated.single_use('gnome.generate_gir argument --warn-error', '0.55.0',
                                         state.subproject, 'Use "fatal_warnings" keyword argument', state.current_node)
        if kwargs['fatal_warnings']:
            scan_command.append('--warn-error')

        generated_files = [f for f in libsources if isinstance(f, (GeneratedList, CustomTarget, CustomTargetIndex))]

        scan_target = self._make_gir_target(
            state, girfile, scan_command, generated_files, depends, scan_env_ldflags,
            # We have to cast here because mypy can't figure this out
            T.cast('T.Dict[str, T.Any]', kwargs))

        typelib_output = f'{ns}-{nsversion}.typelib'
        typelib_cmd = [gicompiler, scan_target, '--output', '@OUTPUT@']
        typelib_cmd += state.get_include_args(gir_inc_dirs, prefix='--includedir=')

        for incdir in typelib_includes:
            typelib_cmd += ["--includedir=" + incdir]

        typelib_target = self._make_typelib_target(state, typelib_output, typelib_cmd, generated_files, T.cast('T.Dict[str, T.Any]', kwargs))

        self._devenv_prepend('GI_TYPELIB_PATH', os.path.join(state.environment.get_build_dir(), state.subdir))

        rv = [scan_target, typelib_target]

        return ModuleReturnValue(rv, rv)

    @noPosargs
    @typed_kwargs('gnome.compile_schemas', _BUILD_BY_DEFAULT.evolve(since='0.40.0'), DEPEND_FILES_KW)
    def compile_schemas(self, state: 'ModuleState', args: T.List['TYPE_var'], kwargs: 'CompileSchemas') -> ModuleReturnValue:
        srcdir = os.path.join(state.build_to_src, state.subdir)
        outdir = state.subdir

        cmd: T.List[T.Union['ToolType', str]] = [self._find_tool(state, 'glib-compile-schemas'), '--targetdir', outdir, srcdir]
        if state.subdir == '':
            targetname = 'gsettings-compile'
        else:
            targetname = 'gsettings-compile-' + state.subdir.replace('/', '_')
        target_g = CustomTarget(
            targetname,
            state.subdir,
            state.subproject,
            state.environment,
            cmd,
            [],
            ['gschemas.compiled'],
            build_by_default=kwargs['build_by_default'],
            depend_files=kwargs['depend_files'],
            description='Compiling gschemas {}',
        )
        self._devenv_prepend('GSETTINGS_SCHEMA_DIR', os.path.join(state.environment.get_build_dir(), state.subdir))
        return ModuleReturnValue(target_g, [target_g])

    @typed_pos_args('gnome.yelp', str, varargs=str)
    @typed_kwargs(
        'gnome.yelp',
        KwargInfo(
            'languages', ContainerTypeInfo(list, str),
            listify=True, default=[],
            deprecated='0.43.0',
            deprecated_message='Use a LINGUAS file in the source directory instead',
        ),
        KwargInfo('media', ContainerTypeInfo(list, str), listify=True, default=[]),
        KwargInfo('sources', ContainerTypeInfo(list, str), listify=True, default=[]),
        KwargInfo('symlink_media', bool, default=True),
    )
    def yelp(self, state: 'ModuleState', args: T.Tuple[str, T.List[str]], kwargs: 'Yelp') -> ModuleReturnValue:
        project_id = args[0]
        sources = kwargs['sources']
        if args[1]:
            FeatureDeprecated.single_use('gnome.yelp more than one positional argument', '0.60.0',
                                         state.subproject, 'use the "sources" keyword argument instead.', state.current_node)
        if not sources:
            sources = args[1]
            if not sources:
                raise MesonException('Yelp requires a list of sources')
        elif args[1]:
            mlog.warning('"gnome.yelp" ignores positional sources arguments when the "sources" keyword argument is set')
        sources_files = [mesonlib.File.from_source_file(state.environment.source_dir,
                                                        os.path.join(state.subdir, 'C'),
                                                        s) for s in sources]

        langs = kwargs['languages']
        if not langs:
            langs = read_linguas(os.path.join(state.environment.source_dir, state.subdir))

        media = kwargs['media']
        symlinks = kwargs['symlink_media']
        targets: T.List[T.Union['build.Target', build.Data, build.SymlinkData]] = []
        potargets: T.List[build.RunTarget] = []

        itstool = state.find_program('itstool')
        msgmerge = state.find_program('msgmerge')
        msgfmt = state.find_program('msgfmt')

        install_dir = os.path.join(state.environment.get_datadir(), 'help')
        c_install_dir = os.path.join(install_dir, 'C', project_id)
        c_data = build.Data(sources_files, c_install_dir, c_install_dir,
                            mesonlib.FileMode(), state.subproject, install_tag='doc')
        targets.append(c_data)

        media_files: T.List[mesonlib.File] = []
        for m in media:
            f = mesonlib.File.from_source_file(state.environment.source_dir,
                                               os.path.join(state.subdir, 'C'), m)
            media_files.append(f)
            m_install_dir = os.path.join(c_install_dir, os.path.dirname(m))
            m_data = build.Data([f], m_install_dir, m_install_dir,
                                mesonlib.FileMode(), state.subproject, install_tag='doc')
            targets.append(m_data)

        pot_file = os.path.join('@SOURCE_ROOT@', state.subdir, 'C', project_id + '.pot')
        pot_sources = [os.path.join('@SOURCE_ROOT@', state.subdir, 'C', s) for s in sources]
        pot_args: T.List[T.Union[ExternalProgram, Executable, OverrideProgram, str]] = [itstool, '-o', pot_file]
        pot_args.extend(pot_sources)
        pottarget = build.RunTarget(f'help-{project_id}-pot', pot_args, [],
                                    os.path.join(state.subdir, 'C'), state.subproject,
                                    state.environment)
        targets.append(pottarget)

        for l in langs:
            l_subdir = os.path.join(state.subdir, l)
            l_install_dir = os.path.join(install_dir, l, project_id)

            for i, m in enumerate(media):
                m_dir = os.path.dirname(m)
                m_install_dir = os.path.join(l_install_dir, m_dir)
                try:
                    m_file: T.Optional[mesonlib.File] = mesonlib.File.from_source_file(state.environment.source_dir, l_subdir, m)
                except MesonException:
                    m_file = None

                l_data: T.Union[build.Data, build.SymlinkData]
                if symlinks and not m_file:
                    link_target = os.path.join(os.path.relpath(c_install_dir, start=m_install_dir), m)
                    l_data = build.SymlinkData(link_target, os.path.basename(m),
                                               m_install_dir, state.subproject, install_tag='doc')
                else:
                    if not m_file:
                        m_file = media_files[i]
                    l_data = build.Data([m_file], m_install_dir, m_install_dir,
                                        mesonlib.FileMode(), state.subproject, install_tag='doc')
                targets.append(l_data)

            po_file = l + '.po'
            po_args: T.List[T.Union[ExternalProgram, Executable, OverrideProgram, str]] = [
                msgmerge, '-q', '-o',
                os.path.join('@SOURCE_ROOT@', l_subdir, po_file),
                os.path.join('@SOURCE_ROOT@', l_subdir, po_file), pot_file]
            potarget = build.RunTarget(f'help-{project_id}-{l}-update-po',
                                       po_args, [pottarget], l_subdir, state.subproject,
                                       state.environment)
            targets.append(potarget)
            potargets.append(potarget)

            gmo_file = project_id + '-' + l + '.gmo'
            gmotarget = CustomTarget(
                f'help-{project_id}-{l}-gmo',
                l_subdir,
                state.subproject,
                state.environment,
                [msgfmt, '@INPUT@', '-o', '@OUTPUT@'],
                [po_file],
                [gmo_file],
                install_tag=['doc'],
                description='Generating yelp doc {}',
            )
            targets.append(gmotarget)

            mergetarget = CustomTarget(
                f'help-{project_id}-{l}',
                l_subdir,
                state.subproject,
                state.environment,
                [itstool, '-m', os.path.join(l_subdir, gmo_file), '--lang', l, '-o', '@OUTDIR@', '@INPUT@'],
                sources_files,
                sources,
                extra_depends=[gmotarget],
                install=True,
                install_dir=[l_install_dir],
                install_tag=['doc'],
                description='Generating yelp doc {}',
            )
            targets.append(mergetarget)

        allpotarget = build.AliasTarget(f'help-{project_id}-update-po', potargets,
                                        state.subdir, state.subproject, state.environment)
        targets.append(allpotarget)

        return ModuleReturnValue(None, targets)

    @typed_pos_args('gnome.gtkdoc', str)
    @typed_kwargs(
        'gnome.gtkdoc',
        KwargInfo('c_args', ContainerTypeInfo(list, str), since='0.48.0', default=[], listify=True),
        KwargInfo('check', bool, default=False, since='0.52.0'),
        KwargInfo('content_files', ContainerTypeInfo(list, (str, mesonlib.File, GeneratedList, CustomTarget, CustomTargetIndex)), default=[], listify=True),
        KwargInfo(
            'dependencies',
            ContainerTypeInfo(list, (Dependency, build.SharedLibrary, build.StaticLibrary)),
            listify=True, default=[]),
        KwargInfo('expand_content_files', ContainerTypeInfo(list, (str, mesonlib.File)), default=[], listify=True),
        KwargInfo('fixxref_args', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo('gobject_typesfile', ContainerTypeInfo(list, (str, mesonlib.File)), default=[], listify=True),
        KwargInfo('html_args', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo('html_assets', ContainerTypeInfo(list, (str, mesonlib.File)), default=[], listify=True),
        KwargInfo('ignore_headers', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo(
            'include_directories',
            ContainerTypeInfo(list, (str, build.IncludeDirs)),
            listify=True, default=[]),
        KwargInfo('install', bool, default=True),
        KwargInfo('install_dir', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo('main_sgml', (str, NoneType)),
        KwargInfo('main_xml', (str, NoneType)),
        KwargInfo('mkdb_args', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo(
            'mode', str, default='auto', since='0.37.0',
            validator=in_set_validator({'xml', 'sgml', 'none', 'auto'})),
        KwargInfo('module_version', str, default='', since='0.48.0'),
        KwargInfo('namespace', str, default='', since='0.37.0'),
        KwargInfo('scan_args', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo('scanobjs_args', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo('src_dir', ContainerTypeInfo(list, (str, build.IncludeDirs)), listify=True, required=True),
    )
    def gtkdoc(self, state: 'ModuleState', args: T.Tuple[str], kwargs: 'GtkDoc') -> ModuleReturnValue:
        modulename = args[0]
        main_file = kwargs['main_sgml']
        main_xml = kwargs['main_xml']
        if main_xml is not None:
            if main_file is not None:
                raise InvalidArguments('gnome.gtkdoc: main_xml and main_sgml are exclusive arguments')
            main_file = main_xml
        moduleversion = kwargs['module_version']
        targetname = modulename + ('-' + moduleversion if moduleversion else '') + '-doc'
        command = state.environment.get_build_command()

        namespace = kwargs['namespace']

        # Ensure we have a C compiler even in C++ projects.
        state.add_language('c', MachineChoice.HOST)

        def abs_filenames(files: T.Iterable['FileOrString']) -> T.Iterator[str]:
            for f in files:
                if isinstance(f, mesonlib.File):
                    yield f.absolute_path(state.environment.get_source_dir(), state.environment.get_build_dir())
                else:
                    yield os.path.join(state.environment.get_source_dir(), state.subdir, f)

        src_dirs = kwargs['src_dir']
        header_dirs: T.List[str] = []
        for src_dir in src_dirs:
            if isinstance(src_dir, build.IncludeDirs):
                header_dirs.extend(src_dir.to_string_list(state.environment.get_source_dir(),
                                                          state.environment.get_build_dir()))
            else:
                header_dirs.append(src_dir)

        t_args: T.List[str] = [
            '--internal', 'gtkdoc',
            '--sourcedir=' + state.environment.get_source_dir(),
            '--builddir=' + state.environment.get_build_dir(),
            '--subdir=' + state.subdir,
            '--headerdirs=' + '@@'.join(header_dirs),
            '--mainfile=' + main_file,
            '--modulename=' + modulename,
            '--moduleversion=' + moduleversion,
            '--mode=' + kwargs['mode']]
        for tool in ['scan', 'scangobj', 'mkdb', 'mkhtml', 'fixxref']:
            program_name = 'gtkdoc-' + tool
            program = state.find_program(program_name)
            path = program.get_path()
            assert path is not None, "This shouldn't be possible since program should be found"
            t_args.append(f'--{program_name}={path}')
        if namespace:
            t_args.append('--namespace=' + namespace)
        exe_wrapper = state.environment.get_exe_wrapper()
        if exe_wrapper:
            t_args.append('--run=' + ' '.join(exe_wrapper.get_command()))
        t_args.append(f'--htmlargs={"@@".join(kwargs["html_args"])}')
        t_args.append(f'--scanargs={"@@".join(kwargs["scan_args"])}')
        t_args.append(f'--scanobjsargs={"@@".join(kwargs["scanobjs_args"])}')
        t_args.append(f'--gobjects-types-file={"@@".join(abs_filenames(kwargs["gobject_typesfile"]))}')
        t_args.append(f'--fixxrefargs={"@@".join(kwargs["fixxref_args"])}')
        t_args.append(f'--mkdbargs={"@@".join(kwargs["mkdb_args"])}')
        t_args.append(f'--html-assets={"@@".join(abs_filenames(kwargs["html_assets"]))}')

        depends: T.List['build.GeneratedTypes'] = []
        content_files = []
        for s in kwargs['content_files']:
            if isinstance(s, (CustomTarget, CustomTargetIndex)):
                depends.append(s)
                for o in s.get_outputs():
                    content_files.append(os.path.join(state.environment.get_build_dir(),
                                                      state.backend.get_target_dir(s),
                                                      o))
            elif isinstance(s, mesonlib.File):
                content_files.append(s.absolute_path(state.environment.get_source_dir(),
                                                     state.environment.get_build_dir()))
            elif isinstance(s, GeneratedList):
                depends.append(s)
                for gen_src in s.get_outputs():
                    content_files.append(os.path.join(state.environment.get_source_dir(),
                                                      state.subdir,
                                                      gen_src))
            else:
                content_files.append(os.path.join(state.environment.get_source_dir(),
                                                  state.subdir,
                                                  s))
        t_args += ['--content-files=' + '@@'.join(content_files)]

        t_args.append(f'--expand-content-files={"@@".join(abs_filenames(kwargs["expand_content_files"]))}')
        t_args.append(f'--ignore-headers={"@@".join(kwargs["ignore_headers"])}')
        t_args.append(f'--installdir={"@@".join(kwargs["install_dir"])}')
        build_args, new_depends = self._get_build_args(kwargs['c_args'], kwargs['include_directories'],
                                                       kwargs['dependencies'], state, depends)
        t_args.extend(build_args)
        new_depends.extend(depends)
        custom_target = CustomTarget(
            targetname,
            state.subdir,
            state.subproject,
            state.environment,
            command + t_args,
            [],
            [f'{modulename}-decl.txt'],
            build_always_stale=True,
            extra_depends=new_depends,
            description='Generating gtkdoc {}',
        )
        alias_target = build.AliasTarget(targetname, [custom_target], state.subdir, state.subproject, state.environment)
        if kwargs['check']:
            check_cmd = state.find_program('gtkdoc-check')
            check_env = ['DOC_MODULE=' + modulename,
                         'DOC_MAIN_SGML_FILE=' + main_file]
            check_args = (targetname + '-check', check_cmd)
            check_workdir = os.path.join(state.environment.get_build_dir(), state.subdir)
            state.test(check_args, env=check_env, workdir=check_workdir, depends=[custom_target])
        res: T.List[T.Union[build.Target, mesonlib.ExecutableSerialisation]] = [custom_target, alias_target]
        if kwargs['install']:
            res.append(state.backend.get_executable_serialisation(command + t_args, tag='doc'))
        return ModuleReturnValue(custom_target, res)

    def _get_build_args(self, c_args: T.List[str], inc_dirs: T.List[T.Union[str, build.IncludeDirs]],
                        deps: T.List[T.Union[Dependency, build.SharedLibrary, build.StaticLibrary]],
                        state: 'ModuleState',
                        depends: T.Sequence[T.Union[build.BuildTarget, 'build.GeneratedTypes']]) -> T.Tuple[
                                T.List[str], T.List[T.Union[build.BuildTarget, 'build.GeneratedTypes', 'FileOrString', build.StructuredSources]]]:
        args: T.List[str] = []
        cflags = c_args.copy()
        deps_cflags, internal_ldflags, external_ldflags, _gi_includes, new_depends = \
            self._get_dependencies_flags(deps, state, depends, include_rpath=True)

        cflags.extend(deps_cflags)
        cflags.extend(state.get_include_args(inc_dirs))
        ldflags: T.List[str] = []
        ldflags.extend(internal_ldflags)
        ldflags.extend(external_ldflags)

        cflags.extend(state.environment.coredata.get_external_args(MachineChoice.HOST, 'c'))
        ldflags.extend(state.environment.coredata.get_external_link_args(MachineChoice.HOST, 'c'))
        compiler = state.environment.coredata.compilers[MachineChoice.HOST]['c']

        compiler_flags = self._get_langs_compilers_flags(state, [('c', compiler)])
        cflags.extend(compiler_flags[0])
        ldflags.extend(compiler_flags[1])
        ldflags.extend(compiler_flags[2])
        if compiler:
            args += ['--cc=%s' % join_args(compiler.get_exelist())]
            args += ['--ld=%s' % join_args(compiler.get_linker_exelist())]
        if cflags:
            args += ['--cflags=%s' % join_args(cflags)]
        if ldflags:
            args += ['--ldflags=%s' % join_args(ldflags)]

        return args, new_depends

    @noKwargs
    @typed_pos_args('gnome.gtkdoc_html_dir', str)
    def gtkdoc_html_dir(self, state: 'ModuleState', args: T.Tuple[str], kwargs: 'TYPE_kwargs') -> str:
        return os.path.join('share/gtk-doc/html', args[0])

    @typed_pos_args('gnome.gdbus_codegen', str, optargs=[(str, mesonlib.File, CustomTarget, CustomTargetIndex, GeneratedList)])
    @typed_kwargs(
        'gnome.gdbus_codegen',
        _BUILD_BY_DEFAULT.evolve(since='0.40.0'),
        DEPENDENCY_SOURCES_KW.evolve(since='0.46.0'),
        KwargInfo('extra_args', ContainerTypeInfo(list, str), since='0.47.0', default=[], listify=True),
        KwargInfo('interface_prefix', (str, NoneType)),
        KwargInfo('namespace', (str, NoneType)),
        KwargInfo('object_manager', bool, default=False),
        KwargInfo(
            'annotations', ContainerTypeInfo(list, (list, str)),
            default=[],
            validator=annotations_validator,
            convertor=lambda x: [x] if x and isinstance(x[0], str) else x,
        ),
        KwargInfo('install_header', bool, default=False, since='0.46.0'),
        KwargInfo('docbook', (str, NoneType)),
        KwargInfo('rst', (str, NoneType), since='1.9.0'),
        KwargInfo('markdown', (str, NoneType), since='1.9.0'),
        KwargInfo(
            'autocleanup', str, default='default', since='0.47.0',
            validator=in_set_validator({'all', 'none', 'objects'})),
        INSTALL_DIR_KW.evolve(since='0.46.0')
    )
    def gdbus_codegen(self, state: 'ModuleState', args: T.Tuple[str, T.Optional[T.Union['FileOrString', build.GeneratedTypes]]],
                      kwargs: 'GdbusCodegen') -> ModuleReturnValue:
        namebase = args[0]
        xml_files: T.List[T.Union['FileOrString', build.GeneratedTypes]] = [args[1]] if args[1] else []
        cmd: T.List[T.Union['ToolType', str]] = [self._find_tool(state, 'gdbus-codegen')]
        cmd.extend(kwargs['extra_args'])

        # Autocleanup supported?
        glib_version = self._get_native_glib_version(state)
        if not mesonlib.version_compare(glib_version, '>= 2.49.1'):
            # Warn if requested, silently disable if not
            if kwargs['autocleanup'] != 'default':
                mlog.warning(f'Glib version ({glib_version}) is too old to support the \'autocleanup\' '
                             'kwarg, need 2.49.1 or newer')
        else:
            # Handle legacy glib versions that don't have autocleanup
            ac = kwargs['autocleanup']
            if ac == 'default':
                ac = 'all'
            cmd.extend(['--c-generate-autocleanup', ac])

        if kwargs['interface_prefix'] is not None:
            cmd.extend(['--interface-prefix', kwargs['interface_prefix']])
        if kwargs['namespace'] is not None:
            cmd.extend(['--c-namespace', kwargs['namespace']])
        if kwargs['object_manager']:
            cmd.extend(['--c-generate-object-manager'])
        xml_files.extend(kwargs['sources'])
        build_by_default = kwargs['build_by_default']

        # Annotations are a bit ugly in that they are a list of lists of strings...
        for annot in kwargs['annotations']:
            cmd.append('--annotate')
            cmd.extend(annot)

        targets = []
        install_header = kwargs['install_header']
        install_dir = kwargs['install_dir'] or state.environment.coredata.optstore.get_value_for(OptionKey('includedir'))
        assert isinstance(install_dir, str), 'for mypy'

        output = namebase + '.c'
        # Added in https://gitlab.gnome.org/GNOME/glib/commit/e4d68c7b3e8b01ab1a4231bf6da21d045cb5a816 (2.55.2)
        # Fixed in https://gitlab.gnome.org/GNOME/glib/commit/cd1f82d8fc741a2203582c12cc21b4dacf7e1872 (2.56.2)
        if mesonlib.version_compare(glib_version, '>= 2.56.2'):
            c_cmd = cmd + ['--body', '--output', '@OUTPUT@', '@INPUT@']
        else:
            if kwargs['docbook'] is not None:
                docbook = kwargs['docbook']

                cmd += ['--generate-docbook', docbook]

            if kwargs['rst'] is not None:
                if not mesonlib.version_compare(glib_version, '>= 2.71.1'):
                    mlog.error(f'Glib version ({glib_version}) is too old to '
                               'support the \'rst\' kwarg, need 2.71.1 or '
                               'newer')

                rst = kwargs['rst']

                cmd += ['--generate-rst', rst]

            if kwargs['markdown'] is not None:
                if not mesonlib.version_compare(glib_version, '>= 2.75.2'):
                    mlog.error(f'Glib version ({glib_version}) is too old to '
                               'support the \'markdown\' kwarg, need 2.75.2 '
                               'or newer')

                markdown = kwargs['markdown']

                cmd += ['--generate-md', markdown]

            # https://git.gnome.org/browse/glib/commit/?id=ee09bb704fe9ccb24d92dd86696a0e6bb8f0dc1a
            if mesonlib.version_compare(glib_version, '>= 2.51.3'):
                cmd += ['--output-directory', '@OUTDIR@', '--generate-c-code', namebase, '@INPUT@']
            else:
                self._print_gdbus_warning()
                cmd += ['--generate-c-code', '@OUTDIR@/' + namebase, '@INPUT@']
            c_cmd = cmd

        cfile_custom_target = CustomTarget(
            output,
            state.subdir,
            state.subproject,
            state.environment,
            c_cmd,
            xml_files,
            [output],
            build_by_default=build_by_default,
            description='Generating gdbus source {}',
        )
        targets.append(cfile_custom_target)

        output = namebase + '.h'
        if mesonlib.version_compare(glib_version, '>= 2.56.2'):
            hfile_cmd = cmd + ['--header', '--output', '@OUTPUT@', '@INPUT@']
            depends = []
        else:
            hfile_cmd = cmd
            depends = [cfile_custom_target]

        hfile_custom_target = CustomTarget(
            output,
            state.subdir,
            state.subproject,
            state.environment,
            hfile_cmd,
            xml_files,
            [output],
            build_by_default=build_by_default,
            extra_depends=depends,
            install=install_header,
            install_dir=[install_dir],
            install_tag=['devel'],
            description='Generating gdbus header {}',
        )
        targets.append(hfile_custom_target)

        if kwargs['docbook'] is not None:
            docbook = kwargs['docbook']
            # The docbook output is always ${docbook}-${name_of_xml_file}
            output = namebase + '-docbook'
            outputs = []
            for f in xml_files:
                outputs.append('{}-{}'.format(docbook, os.path.basename(str(f))))

            if mesonlib.version_compare(glib_version, '>= 2.56.2'):
                docbook_cmd = cmd + ['--output-directory', '@OUTDIR@', '--generate-docbook', docbook, '@INPUT@']
                depends = []
            else:
                docbook_cmd = cmd
                depends = [cfile_custom_target]

            docbook_custom_target = CustomTarget(
                output,
                state.subdir,
                state.subproject,
                state.environment,
                docbook_cmd,
                xml_files,
                outputs,
                build_by_default=build_by_default,
                extra_depends=depends,
                description='Generating gdbus docbook {}',
            )
            targets.append(docbook_custom_target)

        if kwargs['rst'] is not None:
            rst = kwargs['rst']
            # The rst output is always ${rst}-${name_of_xml_file}
            output = namebase + '-rst'
            outputs = []
            for f in xml_files:
                outputs.append('{}-{}'.format(rst, os.path.basename(str(f))))

            rst_custom_target = CustomTarget(
                output,
                state.subdir,
                state.subproject,
                state.environment,
                cmd + ['--output-directory', '@OUTDIR@', '--generate-rst', rst, '@INPUT@'],
                xml_files,
                outputs,
                build_by_default=build_by_default,
                description='Generating gdbus reStructuredText {}',
            )
            targets.append(rst_custom_target)

        if kwargs['markdown'] is not None:
            markdown = kwargs['markdown']
            # The markdown output is always ${markdown}-${name_of_xml_file}
            output = namebase + '-markdown'
            outputs = []
            for f in xml_files:
                outputs.append('{}-{}'.format(markdown, os.path.basename(str(f))))

            markdown_custom_target = CustomTarget(
                output,
                state.subdir,
                state.subproject,
                state.environment,
                cmd + ['--output-directory', '@OUTDIR@', '--generate-md', markdown, '@INPUT@'],
                xml_files,
                outputs,
                build_by_default=build_by_default,
                description='Generating gdbus markdown {}',
            )
            targets.append(markdown_custom_target)

        return ModuleReturnValue(targets, targets)

    @typed_pos_args('gnome.mkenums', str)
    @typed_kwargs(
        'gnome.mkenums',
        *_MK_ENUMS_COMMON_KWS,
        DEPENDS_KW,
        KwargInfo(
            'sources',
            ContainerTypeInfo(list, (str, mesonlib.File, CustomTarget, CustomTargetIndex,
                                     GeneratedList)),
            listify=True,
            required=True,
        ),
        KwargInfo('c_template', (str, mesonlib.File, NoneType)),
        KwargInfo('h_template', (str, mesonlib.File, NoneType)),
        KwargInfo('comments', (str, NoneType)),
        KwargInfo('eprod', (str, NoneType)),
        KwargInfo('fhead', (str, NoneType)),
        KwargInfo('fprod', (str, NoneType)),
        KwargInfo('ftail', (str, NoneType)),
        KwargInfo('vhead', (str, NoneType)),
        KwargInfo('vprod', (str, NoneType)),
        KwargInfo('vtail', (str, NoneType)),
    )
    def mkenums(self, state: 'ModuleState', args: T.Tuple[str], kwargs: 'MkEnums') -> ModuleReturnValue:
        basename = args[0]

        c_template = kwargs['c_template']
        if isinstance(c_template, mesonlib.File):
            c_template = c_template.absolute_path(state.environment.source_dir, state.environment.build_dir)
        h_template = kwargs['h_template']
        if isinstance(h_template, mesonlib.File):
            h_template = h_template.absolute_path(state.environment.source_dir, state.environment.build_dir)

        cmd: T.List[str] = []
        known_kwargs = ['comments', 'eprod', 'fhead', 'fprod', 'ftail',
                        'identifier_prefix', 'symbol_prefix',
                        'vhead', 'vprod', 'vtail']
        for arg in known_kwargs:
            # Mypy can't figure out that this TypedDict index is correct, without repeating T.Literal for the entire list
            if kwargs[arg]:                                         # type: ignore
                cmd += ['--' + arg.replace('_', '-'), kwargs[arg]]  # type: ignore

        targets: T.List[CustomTarget] = []

        h_target: T.Optional[CustomTarget] = None
        if h_template is not None:
            h_output = os.path.basename(os.path.splitext(h_template)[0])
            # We always set template as the first element in the source array
            # so --template consumes it.
            h_cmd = cmd + ['--template', '@INPUT@']
            h_sources: T.List[T.Union[FileOrString, 'build.GeneratedTypes']] = [h_template]
            h_sources.extend(kwargs['sources'])
            h_target = self._make_mkenum_impl(
                state, h_sources, h_output, h_cmd, install=kwargs['install_header'],
                install_dir=kwargs['install_dir'])
            targets.append(h_target)

        if c_template is not None:
            c_output = os.path.basename(os.path.splitext(c_template)[0])
            # We always set template as the first element in the source array
            # so --template consumes it.
            c_cmd = cmd + ['--template', '@INPUT@']
            c_sources: T.List[T.Union[FileOrString, 'build.GeneratedTypes']] = [c_template]
            c_sources.extend(kwargs['sources'])

            depends = kwargs['depends'].copy()
            if h_target is not None:
                depends.append(h_target)
            c_target = self._make_mkenum_impl(
                state, c_sources, c_output, c_cmd, depends=depends)
            targets.insert(0, c_target)

        if c_template is None and h_template is None:
            generic_cmd = cmd + ['@INPUT@']
            target = self._make_mkenum_impl(
                state, kwargs['sources'], basename, generic_cmd,
                install=kwargs['install_header'],
                install_dir=kwargs['install_dir'])
            return ModuleReturnValue(target, [target])
        else:
            return ModuleReturnValue(targets, targets)

    @FeatureNew('gnome.mkenums_simple', '0.42.0')
    @typed_pos_args('gnome.mkenums_simple', str)
    @typed_kwargs(
        'gnome.mkenums_simple',
        *_MK_ENUMS_COMMON_KWS,
        KwargInfo(
            'sources',
            ContainerTypeInfo(list, (str, mesonlib.File)),
            listify=True,
            required=True,
        ),
        KwargInfo('header_prefix', str, default=''),
        KwargInfo('function_prefix', str, default=''),
        KwargInfo('body_prefix', str, default=''),
        KwargInfo('decorator', str, default=''),
    )
    def mkenums_simple(self, state: 'ModuleState', args: T.Tuple[str], kwargs: 'MkEnumsSimple') -> ModuleReturnValue:
        hdr_filename = f'{args[0]}.h'
        body_filename = f'{args[0]}.c'

        header_prefix = kwargs['header_prefix']
        decl_decorator = kwargs['decorator']
        func_prefix = kwargs['function_prefix']
        body_prefix = kwargs['body_prefix']

        cmd: T.List[str] = []
        if kwargs['identifier_prefix']:
            cmd.extend(['--identifier-prefix', kwargs['identifier_prefix']])
        if kwargs['symbol_prefix']:
            cmd.extend(['--symbol-prefix', kwargs['symbol_prefix']])

        c_cmd = cmd.copy()
        # Maybe we should write our own template files into the build dir
        # instead, but that seems like much more work, nice as it would be.
        fhead = ''
        if body_prefix != '':
            fhead += '%s\n' % body_prefix
        fhead += '#include "%s"\n' % hdr_filename
        for hdr in self.interpreter.source_strings_to_files(kwargs['sources']):
            hdr_path = os.path.relpath(hdr.relative_name(), state.subdir)
            fhead += f'#include "{hdr_path}"\n'
        fhead += textwrap.dedent(
            '''
            #define C_ENUM(v) ((gint) v)
            #define C_FLAGS(v) ((guint) v)
            ''')
        c_cmd.extend(['--fhead', fhead])

        c_cmd.append('--fprod')
        c_cmd.append(textwrap.dedent(
            '''
            /* enumerations from "@basename@" */
            '''))

        c_cmd.append('--vhead')
        c_cmd.append(textwrap.dedent(
            f'''
            GType
            {func_prefix}@enum_name@_get_type (void)
            {{
                static gsize gtype_id = 0;
                static const G@Type@Value values[] = {{'''))

        c_cmd.extend(['--vprod', '        { C_@TYPE@ (@VALUENAME@), "@VALUENAME@", "@valuenick@" },'])

        c_cmd.append('--vtail')
        c_cmd.append(textwrap.dedent(
            '''\
                    { 0, NULL, NULL }
                };
                if (g_once_init_enter (&gtype_id)) {
                    GType new_type = g_@type@_register_static (g_intern_static_string ("@EnumName@"), values);
                    g_once_init_leave (&gtype_id, new_type);
                }
                return (GType) gtype_id;
            }'''))
        c_cmd.append('@INPUT@')

        c_file = self._make_mkenum_impl(state, kwargs['sources'], body_filename, c_cmd)

        # .h file generation
        h_cmd = cmd.copy()

        if header_prefix and not header_prefix.endswith('\n'):
            header_prefix += '\n'  # Extra trailing newline for style

        h_cmd.append('--fhead')
        h_cmd.append(textwrap.dedent(
            f'''\
            #pragma once

            #include <glib-object.h>
            {header_prefix}
            G_BEGIN_DECLS
            '''))

        h_cmd.append('--fprod')
        h_cmd.append(textwrap.dedent(
            '''
            /* enumerations from "@basename@" */
            '''))

        extra_newline = ''
        if decl_decorator:
            extra_newline = '\n'  # Extra leading newline for style

        h_cmd.append('--vhead')
        h_cmd.append(extra_newline + textwrap.dedent(
            f'''\
            {decl_decorator}
            GType {func_prefix}@enum_name@_get_type (void);
            #define @ENUMPREFIX@_TYPE_@ENUMSHORT@ ({func_prefix}@enum_name@_get_type())'''))

        h_cmd.append('--ftail')
        h_cmd.append(textwrap.dedent(
            '''
            G_END_DECLS'''))
        h_cmd.append('@INPUT@')

        h_file = self._make_mkenum_impl(
            state, kwargs['sources'], hdr_filename, h_cmd,
            install=kwargs['install_header'],
            install_dir=kwargs['install_dir'])

        return ModuleReturnValue([c_file, h_file], [c_file, h_file])

    def _make_mkenum_impl(
            self,
            state: 'ModuleState',
            sources: T.Sequence[T.Union[str, mesonlib.File, CustomTarget, CustomTargetIndex, GeneratedList]],
            output: str,
            cmd: T.List[str],
            *,
            install: bool = False,
            install_dir: T.Optional[T.Sequence[T.Union[str, bool]]] = None,
            depends: T.Optional[T.Sequence[T.Union[CustomTarget, CustomTargetIndex, BuildTarget]]] = None
            ) -> build.CustomTarget:
        real_cmd: T.List[T.Union[str, 'ToolType']] = [self._find_tool(state, 'glib-mkenums')]
        real_cmd.extend(cmd)
        _install_dir = install_dir or state.environment.coredata.optstore.get_value_for(OptionKey('includedir'))
        assert isinstance(_install_dir, str), 'for mypy'

        return CustomTarget(
            output,
            state.subdir,
            state.subproject,
            state.environment,
            real_cmd,
            sources,
            [output],
            capture=True,
            install=install,
            install_dir=[_install_dir],
            install_tag=['devel'],
            extra_depends=depends,
            # https://github.com/mesonbuild/meson/issues/973
            absolute_paths=True,
            rspable=mesonlib.is_windows() or mesonlib.is_cygwin(),
            description='Generating GObject enum file {}',
        )

    @typed_pos_args('gnome.genmarshal', str)
    @typed_kwargs(
        'gnome.genmarshal',
        DEPEND_FILES_KW.evolve(since='0.61.0'),
        DEPENDS_KW.evolve(since='0.61.0'),
        INSTALL_KW.evolve(name='install_header'),
        INSTALL_DIR_KW,
        KwargInfo('extra_args', ContainerTypeInfo(list, str), listify=True, default=[]),
        KwargInfo('internal', bool, default=False),
        KwargInfo('nostdinc', bool, default=False),
        KwargInfo('prefix', (str, NoneType)),
        KwargInfo('skip_source', bool, default=False),
        KwargInfo('sources', ContainerTypeInfo(list, (str, mesonlib.File), allow_empty=False), listify=True, required=True),
        KwargInfo('stdinc', bool, default=False),
        KwargInfo('valist_marshallers', bool, default=False),
    )
    def genmarshal(self, state: 'ModuleState', args: T.Tuple[str], kwargs: 'GenMarshal') -> ModuleReturnValue:
        output = args[0]
        sources = kwargs['sources']

        new_genmarshal = mesonlib.version_compare(self._get_native_glib_version(state), '>= 2.53.3')

        cmd: T.List[T.Union['ToolType', str]] = [self._find_tool(state, 'glib-genmarshal'), '--quiet']
        if kwargs['prefix']:
            cmd.extend(['--prefix', kwargs['prefix']])
        if kwargs['extra_args']:
            if new_genmarshal:
                cmd.extend(kwargs['extra_args'])
            else:
                mlog.warning('The current version of GLib does not support extra arguments \n'
                             'for glib-genmarshal. You need at least GLib 2.53.3. See ',
                             mlog.bold('https://github.com/mesonbuild/meson/pull/2049'),
                             once=True, fatal=False)
        for k in ['internal', 'nostdinc', 'skip_source', 'stdinc', 'valist_marshallers']:
            # Mypy can't figure out that this TypedDict index is correct, without repeating T.Literal for the entire list
            if kwargs[k]:                                            # type: ignore
                cmd.append(f'--{k.replace("_", "-")}')

        install_header = kwargs['install_header']
        capture = False

        # https://github.com/GNOME/glib/commit/0fbc98097fac4d3e647684f344e508abae109fdf
        if mesonlib.version_compare(self._get_native_glib_version(state), '>= 2.51.0'):
            cmd += ['--output', '@OUTPUT@']
        else:
            capture = True

        header_file = output + '.h'
        h_cmd = cmd + ['--header', '@INPUT@']
        if new_genmarshal:
            h_cmd += ['--pragma-once']
        header = CustomTarget(
            output + '_h',
            state.subdir,
            state.subproject,
            state.environment,
            h_cmd,
            sources,
            [header_file],
            install=install_header,
            install_dir=[kwargs['install_dir']] if kwargs['install_dir'] else [],
            install_tag=['devel'],
            capture=capture,
            depend_files=kwargs['depend_files'],
            description='Generating glib marshaller header {}',
        )

        c_cmd = cmd + ['--body', '@INPUT@']
        extra_deps: T.List[CustomTarget] = []
        if mesonlib.version_compare(self._get_native_glib_version(state), '>= 2.53.4'):
            # Silence any warnings about missing prototypes
            c_cmd += ['--include-header', header_file]
            extra_deps.append(header)
        body = CustomTarget(
            output + '_c',
            state.subdir,
            state.subproject,
            state.environment,
            c_cmd,
            sources,
            [f'{output}.c'],
            capture=capture,
            depend_files=kwargs['depend_files'],
            extra_depends=extra_deps,
            description='Generating glib marshaller source {}',
        )

        rv = [body, header]
        return ModuleReturnValue(rv, rv)

    def _extract_vapi_packages(self, state: 'ModuleState', packages: T.List[T.Union[InternalDependency, str]],
                               ) -> T.Tuple[T.List[str], T.List[VapiTarget], T.List[str], T.List[str], T.List[str]]:
        '''
        Packages are special because we need to:
        - Get a list of packages for the .deps file
        - Get a list of depends for any VapiTargets
        - Get package name from VapiTargets
        - Add include dirs for any VapiTargets
        '''
        if not packages:
            return [], [], [], [], []
        vapi_depends: T.List[VapiTarget] = []
        vapi_packages: T.List[str] = []
        vapi_includes: T.List[str] = []
        vapi_args: T.List[str] = []
        remaining_args = []
        for arg in packages:
            if isinstance(arg, InternalDependency):
                targets = [t for t in arg.sources if isinstance(t, VapiTarget)]
                for target in targets:
                    srcdir = os.path.join(state.environment.get_source_dir(),
                                          target.get_subdir())
                    outdir = os.path.join(state.environment.get_build_dir(),
                                          target.get_subdir())
                    outfile = target.get_outputs()[0][:-5] # Strip .vapi
                    vapi_args.append('--vapidir=' + outdir)
                    vapi_args.append('--girdir=' + outdir)
                    vapi_args.append('--pkg=' + outfile)
                    vapi_depends.append(target)
                    vapi_packages.append(outfile)
                    vapi_includes.append(srcdir)
            else:
                assert isinstance(arg, str), 'for mypy'
                vapi_args.append(f'--pkg={arg}')
                vapi_packages.append(arg)
                remaining_args.append(arg)

        # TODO: this is supposed to take IncludeDirs, but it never worked
        return vapi_args, vapi_depends, vapi_packages, vapi_includes, remaining_args

    def _generate_deps(self, state: 'ModuleState', library: str, packages: T.List[str], install_dir: str) -> build.Data:
        outdir = state.environment.scratch_dir
        fname = os.path.join(outdir, library + '.deps')
        with open(fname, 'w', encoding='utf-8') as ofile:
            for package in packages:
                ofile.write(package + '\n')
        return build.Data([mesonlib.File(True, outdir, fname)], install_dir, install_dir, mesonlib.FileMode(), state.subproject)

    def _get_vapi_link_with(self, target: CustomTarget) -> T.List[build.LibTypes]:
        link_with: T.List[build.LibTypes] = []
        for dep in target.get_target_dependencies():
            if isinstance(dep, build.SharedLibrary):
                link_with.append(dep)
            elif isinstance(dep, GirTarget):
                link_with += self._get_vapi_link_with(dep)
        return link_with

    @typed_pos_args('gnome.generate_vapi', str)
    @typed_kwargs(
        'gnome.generate_vapi',
        INSTALL_KW,
        INSTALL_DIR_KW,
        KwargInfo(
            'sources',
            ContainerTypeInfo(list, (str, GirTarget), allow_empty=False),
            listify=True,
            required=True,
        ),
        KwargInfo('vapi_dirs', ContainerTypeInfo(list, str), listify=True, default=[]),
        KwargInfo('metadata_dirs', ContainerTypeInfo(list, str), listify=True, default=[]),
        KwargInfo('gir_dirs', ContainerTypeInfo(list, str), listify=True, default=[]),
        KwargInfo('packages', ContainerTypeInfo(list, (str, InternalDependency)), listify=True, default=[]),
    )
    def generate_vapi(self, state: 'ModuleState', args: T.Tuple[str], kwargs: 'GenerateVapi') -> ModuleReturnValue:
        created_values: T.List[T.Union[Dependency, build.Data]] = []
        library = args[0]
        build_dir = os.path.join(state.environment.get_build_dir(), state.subdir)
        source_dir = os.path.join(state.environment.get_source_dir(), state.subdir)
        pkg_cmd, vapi_depends, vapi_packages, vapi_includes, packages = self._extract_vapi_packages(state, kwargs['packages'])
        cmd: T.List[T.Union[ExternalProgram, Executable, OverrideProgram, str]]
        cmd = [state.find_program('vapigen'), '--quiet', f'--library={library}', f'--directory={build_dir}']
        cmd.extend([f'--vapidir={d}' for d in kwargs['vapi_dirs']])
        cmd.extend([f'--metadatadir={d}' for d in kwargs['metadata_dirs']])
        cmd.extend([f'--girdir={d}' for d in kwargs['gir_dirs']])
        cmd += pkg_cmd
        cmd += ['--metadatadir=' + source_dir]

        inputs = kwargs['sources']

        link_with: T.List[build.LibTypes] = []
        for i in inputs:
            if isinstance(i, str):
                cmd.append(os.path.join(source_dir, i))
            elif isinstance(i, GirTarget):
                link_with += self._get_vapi_link_with(i)
                subdir = os.path.join(state.environment.get_build_dir(),
                                      i.get_subdir())
                gir_file = os.path.join(subdir, i.get_outputs()[0])
                cmd.append(gir_file)

        vapi_output = library + '.vapi'
        datadir = state.environment.coredata.optstore.get_value_for(OptionKey('datadir'))
        assert isinstance(datadir, str), 'for mypy'
        install_dir = kwargs['install_dir'] or os.path.join(datadir, 'vala', 'vapi')

        if kwargs['install']:
            # We shouldn't need this locally but we install it
            deps_target = self._generate_deps(state, library, vapi_packages, install_dir)
            created_values.append(deps_target)
        vapi_target = VapiTarget(
            vapi_output,
            state.subdir,
            state.subproject,
            state.environment,
            command=cmd,
            sources=inputs,
            outputs=[vapi_output],
            extra_depends=vapi_depends,
            install=kwargs['install'],
            install_dir=[install_dir],
            install_tag=['devel'],
        )

        # So to try our best to get this to just work we need:
        # - link with the correct library
        # - include the vapi and dependent vapi files in sources
        # - add relevant directories to include dirs
        incs = [build.IncludeDirs(state.subdir, ['.'] + vapi_includes, False)]
        sources = [vapi_target] + vapi_depends
        rv = InternalDependency(None, incs, [], [], link_with, [], sources, [], [], {}, [], [], [])
        created_values.append(rv)
        return ModuleReturnValue(rv, created_values)

def initialize(interp: 'Interpreter') -> GnomeModule:
    mod = GnomeModule(interp)
    mod.interpreter.append_holder_map(GResourceTarget, interpreter.CustomTargetHolder)
    mod.interpreter.append_holder_map(GResourceHeaderTarget, interpreter.CustomTargetHolder)
    mod.interpreter.append_holder_map(GirTarget, interpreter.CustomTargetHolder)
    mod.interpreter.append_holder_map(TypelibTarget, interpreter.CustomTargetHolder)
    mod.interpreter.append_holder_map(VapiTarget, interpreter.CustomTargetHolder)
    return mod
