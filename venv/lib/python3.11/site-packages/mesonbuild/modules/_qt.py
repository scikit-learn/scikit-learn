# SPDX-License-Identifier: Apache-2.0
# Copyright 2015 The Meson development team
# Copyright © 2021-2023 Intel Corporation

from __future__ import annotations

import os
import shutil
import typing as T
import xml.etree.ElementTree as ET
import re

from . import ModuleReturnValue, ExtensionModule
from .. import build
from .. import options
from .. import mlog
from ..dependencies import find_external_dependency, Dependency, ExternalLibrary, InternalDependency
from ..mesonlib import MesonException, File, FileMode, version_compare, Popen_safe
from ..interpreter import extract_required_kwarg
from ..interpreter.type_checking import INSTALL_DIR_KW, INSTALL_KW, NoneType
from ..interpreterbase import ContainerTypeInfo, FeatureDeprecated, KwargInfo, noPosargs, FeatureNew, typed_kwargs, typed_pos_args
from ..programs import NonExistingExternalProgram

if T.TYPE_CHECKING:
    from . import ModuleState
    from ..dependencies.qt import QtPkgConfigDependency, QmakeQtDependency
    from ..interpreter import Interpreter
    from ..interpreter import kwargs
    from ..mesonlib import FileOrString
    from ..programs import ExternalProgram
    from typing_extensions import Literal

    QtDependencyType = T.Union[QtPkgConfigDependency, QmakeQtDependency]

    from typing_extensions import TypedDict

    class ResourceCompilerKwArgs(TypedDict):

        """Keyword arguments for the Resource Compiler method."""

        name: T.Optional[str]
        sources: T.Sequence[T.Union[FileOrString, build.CustomTarget, build.CustomTargetIndex, build.GeneratedList]]
        extra_args: T.List[str]
        method: str

    class UICompilerKwArgs(TypedDict):

        """Keyword arguments for the Ui Compiler method."""

        sources: T.Sequence[T.Union[FileOrString, build.CustomTarget, build.CustomTargetIndex, build.GeneratedList]]
        extra_args: T.List[str]
        method: str
        preserve_paths: bool

    class MocCompilerKwArgs(TypedDict):

        """Keyword arguments for the Moc Compiler method."""

        sources: T.Sequence[T.Union[FileOrString, build.CustomTarget, build.CustomTargetIndex, build.GeneratedList]]
        headers: T.Sequence[T.Union[FileOrString, build.CustomTarget, build.CustomTargetIndex, build.GeneratedList]]
        extra_args: T.List[str]
        method: str
        include_directories: T.List[T.Union[str, build.IncludeDirs]]
        dependencies: T.List[T.Union[Dependency, ExternalLibrary]]
        preserve_paths: bool
        output_json: bool

    class PreprocessKwArgs(TypedDict):

        sources: T.List[FileOrString]
        moc_sources: T.List[T.Union[FileOrString, build.CustomTarget]]
        moc_headers: T.List[T.Union[FileOrString, build.CustomTarget]]
        qresources: T.List[FileOrString]
        ui_files: T.List[T.Union[FileOrString, build.CustomTarget]]
        moc_extra_arguments: T.List[str]
        rcc_extra_arguments: T.List[str]
        uic_extra_arguments: T.List[str]
        moc_output_json: bool
        include_directories: T.List[T.Union[str, build.IncludeDirs]]
        dependencies: T.List[T.Union[Dependency, ExternalLibrary]]
        method: str
        preserve_paths: bool

    class HasToolKwArgs(kwargs.ExtractRequired):

        method: str
        tools: T.List[Literal['moc', 'uic', 'rcc', 'lrelease', 'qmlcachegen', 'qmltyperegistrar']]

    class CompileTranslationsKwArgs(TypedDict):

        build_by_default: bool
        install: bool
        install_dir: T.Optional[str]
        method: str
        qresource: T.Optional[str]
        rcc_extra_arguments: T.List[str]
        ts_files: T.List[T.Union[str, File, build.CustomTarget, build.CustomTargetIndex, build.GeneratedList]]

    class GenQrcKwArgs(TypedDict):

        sources: T.Sequence[File]
        aliases: T.Sequence[str]
        prefix: str
        output: str

    class GenQmldirKwArgs(TypedDict):

        module_name: str
        module_version: str
        module_prefix: str
        qml_sources: T.Sequence[T.Union[FileOrString, build.GeneratedTypes]]
        qml_singletons: T.Sequence[T.Union[FileOrString, build.GeneratedTypes]]
        qml_internals: T.Sequence[T.Union[FileOrString, build.GeneratedTypes]]
        designer_supported: bool
        imports: T.List[str]
        optional_imports: T.List[str]
        default_imports: T.List[str]
        depends_imports: T.List[str]
        typeinfo: str
        output: str

    class GenQmlCachegenKwArgs(TypedDict):

        target_name: str
        qml_sources: T.Sequence[T.Union[FileOrString, build.GeneratedTypes]]
        qml_qrc: T.Union[FileOrString, build.GeneratedTypes]
        extra_args: T.List[str]
        module_prefix: str
        method: str

    class GenQmlTypeRegistrarKwArgs(TypedDict):

        target_name: str
        import_name: str
        major_version: str
        minor_version: str
        namespace: str
        typeinfo: str
        generate_qmltype: bool
        collected_json: T.Optional[T.Union[FileOrString, build.CustomTarget]]
        extra_args: T.List[str]
        method: str
        install: bool
        install_dir: T.Optional[str]

    class MocJsonCollectKwArgs(TypedDict):

        target_name: str
        moc_json: T.Sequence[build.GeneratedList]
        method: str

    class QmlModuleKwArgs(TypedDict):

        version: str
        qml_sources: T.List[T.Union[FileOrString, build.GeneratedTypes]]
        qml_singletons: T.List[T.Union[FileOrString, build.GeneratedTypes]]
        qml_internals: T.List[T.Union[FileOrString, build.GeneratedTypes]]
        resources_prefix: str
        moc_headers: T.List[T.Union[FileOrString, build.GeneratedTypes]]
        include_directories: T.List[T.Union[str, build.IncludeDirs]]
        imports: T.List[str]
        optional_imports: T.List[str]
        default_imports: T.List[str]
        depends_imports: T.List[str]
        designer_supported: bool
        namespace: str
        typeinfo: str
        moc_extra_arguments: T.List[str]
        rcc_extra_arguments: T.List[str]
        qmlcachegen_extra_arguments: T.List[str]
        qmltyperegistrar_extra_arguments: T.List[str]
        generate_qmldir: bool
        generate_qmltype: bool
        cachegen: bool
        dependencies: T.List[T.Union[Dependency, ExternalLibrary]]
        method: str
        preserve_paths: bool
        install_dir: str
        install: bool

def _list_in_set_validator(choices: T.Set[str]) -> T.Callable[[T.List[str]], T.Optional[str]]:
    """Check that the choice given was one of the given set."""
    def inner(checklist: T.List[str]) -> T.Optional[str]:
        invalid = set(checklist).difference(choices)
        if invalid:
            return f"invalid selections {', '.join(sorted(invalid))}, valid elements are {', '.join(sorted(choices))}."
        return None

    return inner

#While Qt recomment module name to a be dot separated alphanum, it can technically be
#any well-formed ECMAScript Identifier Name.
#As best effort here we just check for illegal characters
#see https://doc.qt.io/qt-6/qtqml-modules-identifiedmodules.html
_MODULE_NAME_PUNCT = r'- {}<>()[\].:;~%?&,+^=|!\/*"\''
_MODULE_NAME_RE = f'[^{_MODULE_NAME_PUNCT}0-9][^{_MODULE_NAME_PUNCT}]*(\\.[^{_MODULE_NAME_PUNCT}0-9][^{_MODULE_NAME_PUNCT}]*)*'

class QtBaseModule(ExtensionModule):
    _tools_detected = False
    _rcc_supports_depfiles = False
    _moc_supports_depfiles = False
    _set_of_qt_tools = {'moc', 'uic', 'rcc', 'lrelease', 'qmlcachegen', 'qmltyperegistrar'}
    _moc_supports_json = False
    _support_qml_module = False

    def __init__(self, interpreter: Interpreter, qt_version: int = 5):
        ExtensionModule.__init__(self, interpreter)
        self.qt_version = qt_version
        # It is important that this list does not change order as the order of
        # the returned ExternalPrograms will change as well
        self.tools: T.Dict[str, T.Union[ExternalProgram, build.Executable]] = {
            tool: NonExistingExternalProgram(tool) for tool in self._set_of_qt_tools
        }
        self.methods.update({
            'has_tools': self.has_tools,
            'preprocess': self.preprocess,
            'compile_translations': self.compile_translations,
            'compile_resources': self.compile_resources,
            'compile_ui': self.compile_ui,
            'compile_moc': self.compile_moc,
            'qml_module': self.qml_module,
        })

    def compilers_detect(self, state: ModuleState, qt_dep: QtDependencyType) -> None:
        """Detect Qt (4 or 5) moc, uic, rcc in the specified bindir or in PATH"""
        wanted = f'== {qt_dep.version}'

        def gen_bins() -> T.Generator[T.Tuple[str, str], None, None]:
            for b in self.tools:
                if qt_dep.bindir:
                    yield os.path.join(qt_dep.bindir, b), b
                if qt_dep.libexecdir:
                    yield os.path.join(qt_dep.libexecdir, b), b
                # prefer the (official) <tool><version> or (unofficial) <tool>-qt<version>
                # of the tool to the plain one, as we
                # don't know what the unsuffixed one points to without calling it.
                yield f'{b}{qt_dep.qtver}', b
                yield f'{b}-qt{qt_dep.qtver}', b
                yield b, b

        for b, name in gen_bins():
            if self.tools[name].found():
                continue

            if name == 'lrelease':
                arg = ['-version']
            elif version_compare(qt_dep.version, '>= 5'):
                arg = ['--version']
            else:
                arg = ['-v']

            # Ensure that the version of qt and each tool are the same
            def get_version(p: T.Union[ExternalProgram, build.Executable]) -> str:
                _, out, err = Popen_safe(p.get_command() + arg)
                if name == 'lrelease' or not qt_dep.version.startswith('4'):
                    care = out
                else:
                    care = err
                return care.rsplit(' ', maxsplit=1)[-1].replace(')', '').strip()

            p = state.find_program(b, required=False,
                                   version_func=get_version,
                                   wanted=wanted)
            if p.found():
                self.tools[name] = p

    def _detect_tools(self, state: ModuleState, method: str, required: bool = True) -> None:
        if self._tools_detected:
            return
        self._tools_detected = True
        mlog.log(f'Detecting Qt{self.qt_version} tools')
        kwargs = {'required': required, 'modules': 'Core', 'method': method}
        # Just pick one to make mypy happy
        qt = T.cast('QtPkgConfigDependency', find_external_dependency(f'qt{self.qt_version}', state.environment, kwargs))
        if qt.found():
            # Get all tools and then make sure that they are the right version
            self.compilers_detect(state, qt)
            if version_compare(qt.version, '>=6.2.0'):
                #5.1x supports qmlcachegen and other tools to some extend, but arguments/build process marginally differs
                self._support_qml_module = True
            if version_compare(qt.version, '>=5.15.0'):
                self._moc_supports_depfiles = True
                self._moc_supports_json = True
            else:
                mlog.warning('moc dependencies will not work properly until you move to Qt >= 5.15', fatal=False)
            if version_compare(qt.version, '>=5.14.0'):
                self._rcc_supports_depfiles = True
            else:
                mlog.warning('rcc dependencies will not work properly until you move to Qt >= 5.14:',
                             mlog.bold('https://bugreports.qt.io/browse/QTBUG-45460'), fatal=False)
        else:
            suffix = f'-qt{self.qt_version}'
            self.tools['moc'] = NonExistingExternalProgram(name='moc' + suffix)
            self.tools['uic'] = NonExistingExternalProgram(name='uic' + suffix)
            self.tools['rcc'] = NonExistingExternalProgram(name='rcc' + suffix)
            self.tools['lrelease'] = NonExistingExternalProgram(name='lrelease' + suffix)

    @staticmethod
    def _qrc_nodes(state: ModuleState, rcc_file: FileOrString) -> T.Tuple[str, T.List[str]]:
        abspath: str
        if isinstance(rcc_file, str):
            abspath = os.path.join(state.environment.source_dir, state.subdir, rcc_file)
        else:
            abspath = rcc_file.absolute_path(state.environment.source_dir, state.environment.build_dir)
        rcc_dirname = os.path.dirname(abspath)

        # FIXME: what error are we actually trying to check here? (probably parse errors?)
        try:
            tree = ET.parse(abspath)
            root = tree.getroot()
            result: T.List[str] = []
            for child in root[0]:
                if child.tag != 'file':
                    mlog.warning("malformed rcc file: ", os.path.join(state.subdir, str(rcc_file)))
                    break
                elif child.text is None:
                    raise MesonException(f'<file> element without a path in {os.path.join(state.subdir, str(rcc_file))}')
                else:
                    result.append(child.text)

            return rcc_dirname, result
        except MesonException:
            raise
        except Exception:
            raise MesonException(f'Unable to parse resource file {abspath}')

    def _parse_qrc_deps(self, state: ModuleState,
                        rcc_file_: T.Union[FileOrString, build.CustomTarget, build.CustomTargetIndex, build.GeneratedList]) -> T.List[File]:
        result: T.List[File] = []
        inputs: T.Sequence['FileOrString'] = []
        if isinstance(rcc_file_, (str, File)):
            inputs = [rcc_file_]
        else:
            inputs = rcc_file_.get_outputs()

        for rcc_file in inputs:
            rcc_dirname, nodes = self._qrc_nodes(state, rcc_file)
            for resource_path in nodes:
                # We need to guess if the pointed resource is:
                #   a) in build directory -> implies a generated file
                #   b) in source directory
                #   c) somewhere else external dependency file to bundle
                #
                # Also from qrc documentation: relative path are always from qrc file
                # So relative path must always be computed from qrc file !
                if os.path.isabs(resource_path):
                    # a)
                    if resource_path.startswith(os.path.abspath(state.environment.build_dir)):
                        resource_relpath = os.path.relpath(resource_path, state.environment.build_dir)
                        result.append(File(is_built=True, subdir='', fname=resource_relpath))
                    # either b) or c)
                    else:
                        result.append(File(is_built=False, subdir=state.subdir, fname=resource_path))
                else:
                    path_from_rcc = os.path.normpath(os.path.join(rcc_dirname, resource_path))
                    # a)
                    if path_from_rcc.startswith(state.environment.build_dir):
                        result.append(File(is_built=True, subdir=state.subdir, fname=resource_path))
                    # b)
                    else:
                        result.append(File(is_built=False, subdir=state.subdir, fname=path_from_rcc))
        return result

    @FeatureNew('qt.has_tools', '0.54.0')
    @noPosargs
    @typed_kwargs(
        'qt.has_tools',
        KwargInfo('required', (bool, options.UserFeatureOption), default=False),
        KwargInfo('method', str, default='auto'),
        KwargInfo('tools', ContainerTypeInfo(list, str), listify=True,
                  default=['moc', 'uic', 'rcc', 'lrelease'],
                  validator=_list_in_set_validator(_set_of_qt_tools),
                  since='1.6.0'),
    )
    def has_tools(self, state: ModuleState, args: T.Tuple, kwargs: HasToolKwArgs) -> bool:
        method = kwargs.get('method', 'auto')
        # We have to cast here because TypedDicts are invariant, even though
        # ExtractRequiredKwArgs is a subset of HasToolKwArgs, type checkers
        # will insist this is wrong
        disabled, required, feature = extract_required_kwarg(kwargs, state.subproject, default=False)
        if disabled:
            mlog.log('qt.has_tools skipped: feature', mlog.bold(feature), 'disabled')
            return False
        self._detect_tools(state, method, required=False)
        for tool in kwargs['tools']:
            assert tool in self._set_of_qt_tools, f'tools must be in {self._set_of_qt_tools}'
            if not self.tools[tool].found():
                if required:
                    raise MesonException('Qt tools not found')
                return False
        return True

    @FeatureNew('qt.compile_resources', '0.59.0')
    @noPosargs
    @typed_kwargs(
        'qt.compile_resources',
        KwargInfo('name', (str, NoneType)),
        KwargInfo(
            'sources',
            ContainerTypeInfo(list, (File, str, build.CustomTarget, build.CustomTargetIndex, build.GeneratedList), allow_empty=False),
            listify=True,
            required=True,
        ),
        KwargInfo('extra_args', ContainerTypeInfo(list, str), listify=True, default=[]),
        KwargInfo('method', str, default='auto')
    )
    def compile_resources(self, state: 'ModuleState', args: T.Tuple, kwargs: 'ResourceCompilerKwArgs') -> ModuleReturnValue:
        """Compile Qt resources files.

        Uses CustomTargets to generate .cpp files from .qrc files.
        """
        if any(isinstance(s, (build.CustomTarget, build.CustomTargetIndex, build.GeneratedList)) for s in kwargs['sources']):
            FeatureNew.single_use('qt.compile_resources: custom_target or generator for "sources" keyword argument',
                                  '0.60.0', state.subproject, location=state.current_node)
        out = self._compile_resources_impl(state, kwargs)
        return ModuleReturnValue(out, [out])

    def _compile_resources_impl(self, state: 'ModuleState', kwargs: 'ResourceCompilerKwArgs') -> T.List[build.CustomTarget]:
        # Avoid the FeatureNew when dispatching from preprocess
        self._detect_tools(state, kwargs['method'])
        if not self.tools['rcc'].found():
            err_msg = ("{0} sources specified and couldn't find {1}, "
                       "please check your qt{2} installation")
            raise MesonException(err_msg.format('RCC', f'rcc-qt{self.qt_version}', self.qt_version))

        # List of generated CustomTargets
        targets: T.List[build.CustomTarget] = []

        # depfile arguments
        DEPFILE_ARGS: T.List[str] = ['--depfile', '@DEPFILE@'] if self._rcc_supports_depfiles else []

        name = kwargs['name']
        sources: T.List['FileOrString'] = []
        for s in kwargs['sources']:
            if isinstance(s, (str, File)):
                sources.append(s)
            else:
                sources.extend(s.get_outputs())
        extra_args = kwargs['extra_args']

        # If a name was set generate a single .cpp file from all of the qrc
        # files, otherwise generate one .cpp file per qrc file.
        if name:
            qrc_deps: T.List[File] = []
            for s in sources:
                qrc_deps.extend(self._parse_qrc_deps(state, s))

            res_target = build.CustomTarget(
                name,
                state.subdir,
                state.subproject,
                state.environment,
                self.tools['rcc'].get_command() + ['-name', name, '-o', '@OUTPUT@'] + extra_args + ['@INPUT@'] + DEPFILE_ARGS,
                sources,
                [f'{name}.cpp'],
                depend_files=qrc_deps,
                depfile=f'{name}.d',
                description='Compiling Qt resources {}',
            )
            targets.append(res_target)
        else:
            for rcc_file in sources:
                qrc_deps = self._parse_qrc_deps(state, rcc_file)
                if isinstance(rcc_file, str):
                    basename = os.path.basename(rcc_file)
                else:
                    basename = os.path.basename(rcc_file.fname)
                name = f'qt{self.qt_version}-{basename.replace(".", "_")}'
                res_target = build.CustomTarget(
                    name,
                    state.subdir,
                    state.subproject,
                    state.environment,
                    self.tools['rcc'].get_command() + ['-name', '@BASENAME@', '-o', '@OUTPUT@'] + extra_args + ['@INPUT@'] + DEPFILE_ARGS,
                    [rcc_file],
                    [f'{name}.cpp'],
                    depend_files=qrc_deps,
                    depfile=f'{name}.d',
                    description='Compiling Qt resources {}',
                )
                targets.append(res_target)

        return targets

    @FeatureNew('qt.compile_ui', '0.59.0')
    @noPosargs
    @typed_kwargs(
        'qt.compile_ui',
        KwargInfo(
            'sources',
            ContainerTypeInfo(list, (File, str, build.CustomTarget, build.CustomTargetIndex, build.GeneratedList), allow_empty=False),
            listify=True,
            required=True,
        ),
        KwargInfo('extra_args', ContainerTypeInfo(list, str), listify=True, default=[]),
        KwargInfo('method', str, default='auto'),
        KwargInfo('preserve_paths', bool, default=False, since='1.4.0'),
    )
    def compile_ui(self, state: ModuleState, args: T.Tuple, kwargs: UICompilerKwArgs) -> ModuleReturnValue:
        """Compile UI resources into cpp headers."""
        if any(isinstance(s, (build.CustomTarget, build.CustomTargetIndex, build.GeneratedList)) for s in kwargs['sources']):
            FeatureNew.single_use('qt.compile_ui: custom_target or generator for "sources" keyword argument',
                                  '0.60.0', state.subproject, location=state.current_node)
        out = self._compile_ui_impl(state, kwargs)
        return ModuleReturnValue(out, [out])

    def _compile_ui_impl(self, state: ModuleState, kwargs: UICompilerKwArgs) -> build.GeneratedList:
        # Avoid the FeatureNew when dispatching from preprocess
        self._detect_tools(state, kwargs['method'])
        if not self.tools['uic'].found():
            err_msg = ("{0} sources specified and couldn't find {1}, "
                       "please check your qt{2} installation")
            raise MesonException(err_msg.format('UIC', f'uic-qt{self.qt_version}', self.qt_version))

        preserve_path_from = os.path.join(state.source_root, state.subdir) if kwargs['preserve_paths'] else None
        # TODO: This generator isn't added to the generator list in the Interpreter
        gen = build.Generator(
            self.tools['uic'],
            kwargs['extra_args'] + ['-o', '@OUTPUT@', '@INPUT@'],
            ['ui_@BASENAME@.h'],
            name=f'Qt{self.qt_version} ui')
        return gen.process_files(kwargs['sources'], state, preserve_path_from)

    @FeatureNew('qt.compile_moc', '0.59.0')
    @noPosargs
    @typed_kwargs(
        'qt.compile_moc',
        KwargInfo(
            'sources',
            ContainerTypeInfo(list, (File, str, build.CustomTarget, build.CustomTargetIndex, build.GeneratedList)),
            listify=True,
            default=[],
        ),
        KwargInfo(
            'headers',
            ContainerTypeInfo(list, (File, str, build.CustomTarget, build.CustomTargetIndex, build.GeneratedList)),
            listify=True,
            default=[]
        ),
        KwargInfo('extra_args', ContainerTypeInfo(list, str), listify=True, default=[]),
        KwargInfo('method', str, default='auto'),
        KwargInfo('include_directories', ContainerTypeInfo(list, (build.IncludeDirs, str)), listify=True, default=[]),
        KwargInfo('dependencies', ContainerTypeInfo(list, (Dependency, ExternalLibrary)), listify=True, default=[]),
        KwargInfo('preserve_paths', bool, default=False, since='1.4.0'),
        KwargInfo('output_json', bool, default=False, since='1.7.0'),
    )
    def compile_moc(self, state: ModuleState, args: T.Tuple, kwargs: MocCompilerKwArgs) -> ModuleReturnValue:
        if any(isinstance(s, (build.CustomTarget, build.CustomTargetIndex, build.GeneratedList)) for s in kwargs['headers']):
            FeatureNew.single_use('qt.compile_moc: custom_target or generator for "headers" keyword argument',
                                  '0.60.0', state.subproject, location=state.current_node)
        if any(isinstance(s, (build.CustomTarget, build.CustomTargetIndex, build.GeneratedList)) for s in kwargs['sources']):
            FeatureNew.single_use('qt.compile_moc: custom_target or generator for "sources" keyword argument',
                                  '0.60.0', state.subproject, location=state.current_node)
        out = self._compile_moc_impl(state, kwargs)
        return ModuleReturnValue(out, [out])

    def _compile_moc_impl(self, state: ModuleState, kwargs: MocCompilerKwArgs) -> T.List[build.GeneratedList]:
        # Avoid the FeatureNew when dispatching from preprocess
        self._detect_tools(state, kwargs['method'])
        if not self.tools['moc'].found():
            err_msg = ("{0} sources specified and couldn't find {1}, "
                       "please check your qt{2} installation")
            raise MesonException(err_msg.format('MOC', f'uic-qt{self.qt_version}', self.qt_version))

        if not (kwargs['headers'] or kwargs['sources']):
            raise build.InvalidArguments('At least one of the "headers" or "sources" keyword arguments must be provided and not empty')

        inc = state.get_include_args(include_dirs=kwargs['include_directories'])
        compile_args: T.List[str] = []
        for dep in kwargs['dependencies']:
            compile_args.extend(a for a in dep.get_all_compile_args() if a.startswith(('-I', '-D')))
            if isinstance(dep, InternalDependency):
                for incl in dep.include_directories:
                    compile_args.extend(f'-I{i}' for i in incl.to_string_list(self.interpreter.source_root, self.interpreter.environment.build_dir))

        output: T.List[build.GeneratedList] = []

        do_output_json: bool = kwargs['output_json']
        if do_output_json and not self._moc_supports_json:
            raise MesonException(f'moc-qt{self.qt_version} doesn\'t support "output_json" option')

        # depfile arguments (defaults to <output-name>.d)
        DEPFILE_ARGS: T.List[str] = ['--output-dep-file'] if self._moc_supports_depfiles else []
        JSON_ARGS: T.List[str] = ['--output-json'] if do_output_json else []

        arguments = kwargs['extra_args'] + DEPFILE_ARGS + JSON_ARGS + inc + compile_args + ['@INPUT@', '-o', '@OUTPUT0@']
        preserve_path_from = os.path.join(state.source_root, state.subdir) if kwargs['preserve_paths'] else None
        if kwargs['headers']:
            header_gen_output: T.List[str] = ['moc_@BASENAME@.cpp']
            if do_output_json:
                header_gen_output.append('moc_@BASENAME@.cpp.json')
            moc_gen = build.Generator(
                self.tools['moc'], arguments, header_gen_output,
                depfile='moc_@BASENAME@.cpp.d',
                name=f'Qt{self.qt_version} moc header')
            output.append(moc_gen.process_files(kwargs['headers'], state, preserve_path_from))
        if kwargs['sources']:
            source_gen_output: T.List[str] = ['@BASENAME@.moc']
            if do_output_json:
                source_gen_output.append('@BASENAME@.moc.json')
            moc_gen = build.Generator(
                self.tools['moc'], arguments, source_gen_output,
                depfile='@BASENAME@.moc.d',
                name=f'Qt{self.qt_version} moc source')
            output.append(moc_gen.process_files(kwargs['sources'], state, preserve_path_from))

        return output

    # We can't use typed_pos_args here, the signature is ambiguous
    @typed_kwargs(
        'qt.preprocess',
        KwargInfo('sources', ContainerTypeInfo(list, (File, str)), listify=True, default=[], deprecated='0.59.0'),
        KwargInfo('qresources', ContainerTypeInfo(list, (File, str)), listify=True, default=[]),
        KwargInfo('ui_files', ContainerTypeInfo(list, (File, str, build.CustomTarget)), listify=True, default=[]),
        KwargInfo('moc_sources', ContainerTypeInfo(list, (File, str, build.CustomTarget)), listify=True, default=[]),
        KwargInfo('moc_headers', ContainerTypeInfo(list, (File, str, build.CustomTarget)), listify=True, default=[]),
        KwargInfo('moc_extra_arguments', ContainerTypeInfo(list, str), listify=True, default=[], since='0.44.0'),
        KwargInfo('rcc_extra_arguments', ContainerTypeInfo(list, str), listify=True, default=[], since='0.49.0'),
        KwargInfo('uic_extra_arguments', ContainerTypeInfo(list, str), listify=True, default=[], since='0.49.0'),
        KwargInfo('method', str, default='auto'),
        KwargInfo('include_directories', ContainerTypeInfo(list, (build.IncludeDirs, str)), listify=True, default=[]),
        KwargInfo('dependencies', ContainerTypeInfo(list, (Dependency, ExternalLibrary)), listify=True, default=[]),
        KwargInfo('preserve_paths', bool, default=False, since='1.4.0'),
        KwargInfo('moc_output_json', bool, default=False, since='1.7.0'),
    )
    def preprocess(self, state: ModuleState, args: T.List[T.Union[str, File]], kwargs: PreprocessKwArgs) -> ModuleReturnValue:
        _sources = args[1:]
        if _sources:
            FeatureDeprecated.single_use('qt.preprocess positional sources', '0.59', state.subproject, location=state.current_node)
        # List is invariant, os we have to cast...
        sources = T.cast('T.List[T.Union[str, File, build.GeneratedList, build.CustomTarget]]',
                         _sources + kwargs['sources'])
        for s in sources:
            if not isinstance(s, (str, File)):
                raise build.InvalidArguments('Variadic arguments to qt.preprocess must be Strings or Files')
        method = kwargs['method']

        if kwargs['qresources']:
            # custom output name set? -> one output file, multiple otherwise
            rcc_kwargs: ResourceCompilerKwArgs = {'name': '', 'sources': kwargs['qresources'], 'extra_args': kwargs['rcc_extra_arguments'], 'method': method}
            if args:
                name = args[0]
                if not isinstance(name, str):
                    raise build.InvalidArguments('First argument to qt.preprocess must be a string')
                rcc_kwargs['name'] = name
            sources.extend(self._compile_resources_impl(state, rcc_kwargs))

        if kwargs['ui_files']:
            ui_kwargs: UICompilerKwArgs = {
                'sources': kwargs['ui_files'],
                'extra_args': kwargs['uic_extra_arguments'],
                'method': method,
                'preserve_paths': kwargs['preserve_paths'],
            }
            sources.append(self._compile_ui_impl(state, ui_kwargs))

        if kwargs['moc_headers'] or kwargs['moc_sources']:
            moc_kwargs: MocCompilerKwArgs = {
                'extra_args': kwargs['moc_extra_arguments'],
                'sources': kwargs['moc_sources'],
                'headers': kwargs['moc_headers'],
                'include_directories': kwargs['include_directories'],
                'dependencies': kwargs['dependencies'],
                'method': method,
                'preserve_paths': kwargs['preserve_paths'],
                'output_json': kwargs['moc_output_json']
            }
            sources.extend(self._compile_moc_impl(state, moc_kwargs))

        return ModuleReturnValue(sources, [sources])

    @FeatureNew('qt.compile_translations', '0.44.0')
    @noPosargs
    @typed_kwargs(
        'qt.compile_translations',
        KwargInfo('build_by_default', bool, default=False),
        INSTALL_KW,
        INSTALL_DIR_KW,
        KwargInfo('method', str, default='auto'),
        KwargInfo('qresource', (str, NoneType), since='0.56.0'),
        KwargInfo('rcc_extra_arguments', ContainerTypeInfo(list, str), listify=True, default=[], since='0.56.0'),
        KwargInfo('ts_files', ContainerTypeInfo(list, (str, File, build.CustomTarget, build.CustomTargetIndex, build.GeneratedList)), listify=True, default=[]),
    )
    def compile_translations(self, state: ModuleState, args: T.Tuple, kwargs: CompileTranslationsKwArgs) -> ModuleReturnValue:
        ts_files = kwargs['ts_files']
        if any(isinstance(s, (build.CustomTarget, build.CustomTargetIndex, build.GeneratedList)) for s in ts_files):
            FeatureNew.single_use('qt.compile_translations: custom_target or generator for "ts_files" keyword argument',
                                  '0.60.0', state.subproject, location=state.current_node)
        if kwargs['install'] and not kwargs['install_dir']:
            raise MesonException('qt.compile_translations: "install_dir" keyword argument must be set when "install" is true.')
        qresource = kwargs['qresource']
        if qresource:
            if ts_files:
                raise MesonException('qt.compile_translations: Cannot specify both ts_files and qresource')
            if os.path.dirname(qresource) != '':
                raise MesonException('qt.compile_translations: qresource file name must not contain a subdirectory.')
            qresource_file = File.from_built_file(state.subdir, qresource)
            infile_abs = os.path.join(state.environment.source_dir, qresource_file.relative_name())
            outfile_abs = os.path.join(state.environment.build_dir, qresource_file.relative_name())
            os.makedirs(os.path.dirname(outfile_abs), exist_ok=True)
            shutil.copy2(infile_abs, outfile_abs)
            self.interpreter.add_build_def_file(infile_abs)

            _, nodes = self._qrc_nodes(state, qresource_file)
            for c in nodes:
                if c.endswith('.qm'):
                    ts_files.append(c.rstrip('.qm') + '.ts')
                else:
                    raise MesonException(f'qt.compile_translations: qresource can only contain qm files, found {c}')
            results = self.preprocess(state, [], {'qresources': qresource_file, 'rcc_extra_arguments': kwargs['rcc_extra_arguments']})
        self._detect_tools(state, kwargs['method'])
        translations: T.List[build.CustomTarget] = []
        for ts in ts_files:
            if not self.tools['lrelease'].found():
                raise MesonException('qt.compile_translations: ' +
                                     self.tools['lrelease'].name + ' not found')
            if qresource:
                # In this case we know that ts_files is always a List[str], as
                # it's generated above and no ts_files are passed in. However,
                # mypy can't figure that out so we use assert to assure it that
                # what we're doing is safe
                assert isinstance(ts, str), 'for mypy'
                outdir = os.path.dirname(os.path.normpath(os.path.join(state.subdir, ts)))
                ts = os.path.basename(ts)
            else:
                outdir = state.subdir
            cmd: T.List[T.Union[ExternalProgram, build.Executable, str]] = [self.tools['lrelease'], '@INPUT@', '-qm', '@OUTPUT@']
            lrelease_target = build.CustomTarget(
                f'qt{self.qt_version}-compile-{ts}',
                outdir,
                state.subproject,
                state.environment,
                cmd,
                [ts],
                ['@BASENAME@.qm'],
                install=kwargs['install'],
                install_dir=[kwargs['install_dir']],
                install_tag=['i18n'],
                build_by_default=kwargs['build_by_default'],
                description='Compiling Qt translations {}',
            )
            translations.append(lrelease_target)
        if qresource:
            return ModuleReturnValue(results.return_value[0], [results.new_objects, translations])
        else:
            return ModuleReturnValue(translations, [translations])

    def _source_to_files(self, state: ModuleState, sources: T.Sequence[T.Union[FileOrString, build.GeneratedTypes]]) -> T.List[File]:

        content_files = []
        for s in sources:
            if isinstance(s, (build.CustomTarget, build.CustomTargetIndex)):
                for o in s.get_outputs():
                    content_files.append(File.from_built_file(state.backend.get_target_dir(s), o))
            elif isinstance(s, File):
                content_files.append(s)
            elif isinstance(s, build.GeneratedList):
                for gen_src in s.get_outputs():
                    content_files.append(File.from_built_file(state.subdir, gen_src))
            else:
                content_files.append(File.from_source_file(
                    state.environment.get_source_dir(),
                    state.subdir,
                    s
                ))
        return content_files

    def _gen_qrc(self, state: ModuleState, kwargs: GenQrcKwArgs) -> File:

        fileout = File.from_built_file(state.subdir, kwargs['output'])
        fileout_abs = os.path.join(state.environment.build_dir, fileout.relative_name())
        if not os.path.isdir(state.environment.build_dir):
            os.mkdir(state.environment.build_dir)

        rcc = ET.Element('RCC')
        qresource = ET.SubElement(rcc, 'qresource', prefix='/' + kwargs['prefix'])
        assert (len(kwargs['sources']) == len(kwargs['aliases']))
        for source, alias in zip(kwargs['sources'], kwargs['aliases']):
            filenode = ET.SubElement(qresource, 'file', alias=alias)
            filenode.text = source.absolute_path(
                state.environment.get_source_dir(),
                state.environment.get_build_dir()
            )

        tree = ET.ElementTree(rcc)
        tree.write(fileout_abs)
        return fileout

    def _gen_qmldir(self, state: ModuleState, kwargs: GenQmldirKwArgs) -> File:
        module_name: str = kwargs['module_name']
        module_version: str = kwargs['module_version']
        module_prefix: str = kwargs['module_prefix']
        designer_supported: bool = kwargs['designer_supported']
        typeinfo_file: str = kwargs['typeinfo']

        #Foo.Bar/1.0 foo.bar/auto foo.bar
        import_re = re.compile(r'^(' + _MODULE_NAME_RE + r')(/((\d+(\.\d+)?)|auto))?$')

        fileout = File.from_built_file(state.subdir, kwargs['output'])
        fileout_abs = os.path.join(state.environment.build_dir, fileout.relative_name())
        if not os.path.isdir(state.environment.build_dir):
            os.mkdir(state.environment.build_dir)

        with open(fileout_abs, 'w', encoding='utf-8') as fd:

            def __gen_import(import_type: str, importlist: T.Sequence[str]) -> None:
                for import_string in importlist:
                    match = import_re.match(import_string)
                    if not match:
                        raise MesonException(f'invalid syntax for qml import {import_string}')
                    module: str = match.group(1)
                    version: str = match.group(4) or ''
                    fd.write(f'{import_type} {module} {version}\n')

            def __gen_declaration(qualifier: str, version: str, importlist: T.Sequence[T.Union[FileOrString, build.GeneratedTypes]]) -> None:
                importpathlist = self._source_to_files(state, importlist)
                for s in importpathlist:
                    basename: str = os.path.basename(s.fname)
                    classname: str = basename.rsplit('.', maxsplit=1)[0]

                    if not basename.endswith(('.qml', '.js', '.mjs')):
                        raise MesonException(f'unexpected file type declared in qml sources {s}')

                    if not classname or '.' in classname or classname[0].islower():
                        raise MesonException(f'{basename} is not a valid QML file name')
                    if version:
                        fd.write(f'{qualifier}{classname} {version} {basename}\n')
                    else:
                        fd.write(f'{qualifier}{classname} {basename}\n')

            fd.write(f'module {module_name}\n')
            fd.write(f'prefer :/{module_prefix}/\n')

            __gen_import('import', kwargs['imports'])
            __gen_import('optional import', kwargs['optional_imports'])
            __gen_import('default import', kwargs['default_imports'])
            __gen_import('depends', kwargs['depends_imports'])
            __gen_declaration('', module_version, kwargs['qml_sources'])
            __gen_declaration('singleton ', module_version, kwargs['qml_singletons'])
            __gen_declaration('internal ', '', kwargs['qml_internals'])

            if typeinfo_file:
                fd.write(f'typeinfo {typeinfo_file}\n')

            if designer_supported:
                fd.write('designersupported\n')
        return fileout

    def _moc_json_collect(self, state: ModuleState, kwargs: MocJsonCollectKwArgs) -> build.CustomTarget:
        self._detect_tools(state, kwargs['method'])
        if not self.tools['moc'].found():
            raise MesonException('qt.qml_module: ' +
                                 self.tools['moc'].name + ' not found')

        target_name: str = kwargs['target_name']
        moc_json: T.Sequence[build.GeneratedList] = kwargs['moc_json']

        #there may be a better way :-/
        input_args: T.List[str] = []
        input_counter = 0
        for g in moc_json:
            for fname in g.get_outputs():
                if fname.endswith('.json'):
                    input_args.append(f'@INPUT{input_counter}@')
                input_counter += 1

        return build.CustomTarget(
            f'moc_collect_json_{target_name}',
            state.subdir,
            state.subproject,
            state.environment,
            self.tools['moc'].get_command() + ['--collect-json', '-o', '@OUTPUT@'] + input_args,
            moc_json,
            [f'{target_name}_json_collect.json'],
            description=f'Collecting json type information for {target_name}',
        )

    def _gen_qml_cachegen(self, state: ModuleState, kwargs: GenQmlCachegenKwArgs) -> T.List[T.Union[build.CustomTarget, build.GeneratedList]]:
        self._detect_tools(state, kwargs['method'])
        if not self.tools['qmlcachegen'].found():
            raise MesonException('qt.qml_module: ' +
                                 self.tools['qmlcachegen'].name + ' not found')

        target_name: str = kwargs['target_name']

        command_args = ['-o', '@OUTPUT@'] + kwargs['extra_args']
        for qrc in self._source_to_files(state, [kwargs['qml_qrc']]):
            command_args.extend(['--resource', qrc.absolute_path(
                state.environment.get_source_dir(),
                state.environment.get_build_dir()
            )])

        command_args.append('@INPUT@')

        cache_gen = build.Generator(
            self.tools['qmlcachegen'],
            command_args,
            [f'{target_name}_@BASENAME@.cpp'],
            name=f'Qml cache generation for {target_name}')

        output: T.List[T.Union[build.CustomTarget, build.GeneratedList]] = []
        output.append(cache_gen.process_files(kwargs['qml_sources'], state))

        cachegen_inputs: T.List[str] = []
        qml_sources_paths = self._source_to_files(state, kwargs['qml_sources'])
        for s in qml_sources_paths:
            source_basename = os.path.basename(s.fname)
            ressource_path = os.path.join('/', kwargs['module_prefix'], source_basename)
            cachegen_inputs.append(ressource_path)

        cacheloader_target = build.CustomTarget(
            f'cacheloader_{target_name}',
            state.subdir,
            state.subproject,
            state.environment,
            self.tools['qmlcachegen'].get_command() + ['-o', '@OUTPUT@'] + ['--resource-name', f'qmlcache_{target_name}'] + kwargs['extra_args'] + ['--resource=@INPUT@'] + cachegen_inputs,
            [kwargs['qml_qrc']],
            #output name format matters here
            [f'{target_name}_qmlcache_loader.cpp'],
            description=f'Qml cache loader for {target_name}',
        )
        output.append(cacheloader_target)
        return output

    def _qml_type_registrar(self, state: ModuleState, kwargs: GenQmlTypeRegistrarKwArgs) -> build.CustomTarget:
        self._detect_tools(state, kwargs['method'])
        if not self.tools['qmltyperegistrar'].found():
            raise MesonException('qt.qml_module: ' +
                                 self.tools['qmltyperegistrar'].name + ' not found')

        import_name: str = kwargs['import_name']
        major_version: str = kwargs['major_version']
        minor_version: str = kwargs['minor_version']
        namespace: str = kwargs['namespace']
        typeinfo: str = kwargs['typeinfo']
        target_name: str = kwargs['target_name']
        collected_json: T.Optional[T.Union[FileOrString, build.CustomTarget]] = kwargs['collected_json']

        inputs: T.Sequence[T.Union[FileOrString, build.CustomTarget]] = [collected_json] if collected_json else []
        outputs: T.List[str] = [f'{target_name}_qmltyperegistrations.cpp']
        install_dir: T.List[T.Union[str, Literal[False]]] = [False]
        install_tag: T.List[T.Union[str, None]] = [None]

        cmd = self.tools['qmltyperegistrar'].get_command() + [
            '--import-name', import_name,
            '--major-version', major_version,
            '--minor-version', minor_version,
            '-o', '@OUTPUT0@',
        ]

        cmd.extend(kwargs['extra_args'])

        if namespace:
            cmd.extend(['--namespace', namespace])

        if kwargs['generate_qmltype']:
            cmd.extend(['--generate-qmltypes', '@OUTPUT1@'])
            if typeinfo == '':
                outputs.append(f'{target_name}.qmltypes')
            else:
                outputs.append(f'{typeinfo}')
            install_dir.append(kwargs['install_dir'])
            install_tag.append('devel')

        if collected_json:
            cmd.append('@INPUT@')

        return build.CustomTarget(
            f'typeregistrar_{target_name}',
            state.subdir,
            state.subproject,
            state.environment,
            cmd,
            inputs,
            outputs,
            install=kwargs['install'],
            install_dir=install_dir,
            install_tag=install_tag,
            description=f'Qml type registration for {target_name}',
        )

    @FeatureNew('qt.qml_module', '1.7')
    @typed_pos_args('qt.qml_module', str)
    @typed_kwargs(
        'qt.qml_module',
        KwargInfo('version', str, default='254.254'),
        #qml sources
        KwargInfo('qml_sources', ContainerTypeInfo(list, (File, str, build.CustomTarget)), listify=True, default=[]),
        KwargInfo('qml_singletons', ContainerTypeInfo(list, (File, str, build.CustomTarget)), listify=True, default=[]),
        KwargInfo('qml_internals', ContainerTypeInfo(list, (File, str, build.CustomTarget)), listify=True, default=[]),
        KwargInfo('resources_prefix', str, default='qt/qml'),
        #qmldir generation
        KwargInfo('imports', ContainerTypeInfo(list, (str)), default=[]),
        KwargInfo('optional_imports', ContainerTypeInfo(list, (str)), default=[]),
        KwargInfo('default_imports', ContainerTypeInfo(list, (str)), default=[]),
        #match DEPENDENCIES argument from CMake, but dependencies keyword is already taken
        KwargInfo('depends_imports', ContainerTypeInfo(list, (str)), default=[]),
        KwargInfo('designer_supported', bool, default=False),
        #for type registration, same arguments as moc
        #moc_sources is voluntary ommited as typeregistrar needs to import a header
        KwargInfo('moc_headers', ContainerTypeInfo(list, (File, str, build.CustomTarget)), listify=True, default=[]),
        KwargInfo('include_directories', ContainerTypeInfo(list, (build.IncludeDirs, str)), listify=True, default=[]),
        KwargInfo('namespace', str, default=''),
        KwargInfo('typeinfo', str, default=''),

        KwargInfo('moc_extra_arguments', ContainerTypeInfo(list, str), listify=True, default=[]),
        KwargInfo('rcc_extra_arguments', ContainerTypeInfo(list, str), listify=True, default=[]),
        KwargInfo('qmlcachegen_extra_arguments', ContainerTypeInfo(list, str), listify=True, default=[]),
        KwargInfo('qmltyperegistrar_extra_arguments', ContainerTypeInfo(list, str), listify=True, default=[]),

        KwargInfo('generate_qmldir', bool, default=True),
        KwargInfo('generate_qmltype', bool, default=True),
        KwargInfo('cachegen', bool, default=True),

        KwargInfo('dependencies', ContainerTypeInfo(list, (Dependency, ExternalLibrary)), listify=True, default=[]),
        INSTALL_DIR_KW,
        INSTALL_KW,
        KwargInfo('method', str, default='auto'),
        KwargInfo('preserve_paths', bool, default=False),
    )
    def qml_module(self, state: ModuleState, args: T.Tuple[str], kwargs: QmlModuleKwArgs) -> ModuleReturnValue:

        self._detect_tools(state, kwargs['method'])
        if not self._support_qml_module:
            raise MesonException('qt.qml_module is not suppported for this version of Qt')

        #Major.Minor(.Patch)
        version_re = re.compile(r'^(\d+)\.(\d+)(\.(\d+))?$')
        module_name_re = re.compile(_MODULE_NAME_RE)

        output: T.List[T.Union[build.CustomTarget, build.GeneratedList]] = []

        module_name: str = args[0]
        if not module_name_re.fullmatch(module_name):
            raise MesonException(f'qml module URI should be in the form Foo.Bar.xxx, got {module_name}')

        module_version: str = kwargs['version']
        module_version_match = version_re.match(module_version)
        if not module_version_match:
            raise MesonException(f'qml module version should be in the form Major.Minor, got {module_version}')
        module_version_major: str = module_version_match.group(1)
        module_version_minor: str = module_version_match.group(2)
        #qt ignores .patch version
        module_version_short = f'{module_version_major}.{module_version_minor}'

        module_prefix_list: T.List[str] = module_name.split('.')
        module_prefix: str = os.path.join(*module_prefix_list)
        module_prefix_full: str = os.path.join(*(kwargs['resources_prefix'].split('/') + module_prefix_list))

        #same format as the one derived from qmltyperegistrar
        target_name = re.sub(r'[^A-Za-z0-9]', '_', module_name)

        qrc_resouces: T.List[T.Union[FileOrString, build.GeneratedTypes]] = []
        all_qml: T.Sequence[T.Union[FileOrString, build.GeneratedTypes]] = kwargs['qml_sources'] + kwargs['qml_singletons'] + kwargs['qml_internals']
        all_qml_files: T.List[File] = self._source_to_files(state, all_qml)
        all_qml_basename: T.List[str] = [os.path.basename(p.fname) for p in all_qml_files]

        install_dir: str = kwargs['install_dir'] or 'qml'
        module_install_dir: str = os.path.join(install_dir, module_prefix)

        if len(all_qml) != 0:
            qml_qrc_kwargs: GenQrcKwArgs = {
                'output': f'{target_name}_qml.qrc',
                'sources': all_qml_files,
                'aliases': all_qml_basename,
                'prefix': module_prefix_full,
            }
            qml_qrc = self._gen_qrc(state, qml_qrc_kwargs)

            if not kwargs['cachegen']:
                qrc_resouces.append(qml_qrc)
            else:
                cachegen_kwargs: GenQmlCachegenKwArgs = {
                    'target_name': target_name,
                    'qml_qrc': qml_qrc,
                    'qml_sources': all_qml,
                    'module_prefix': module_prefix_full,
                    'extra_args': kwargs['qmlcachegen_extra_arguments'],
                    'method': kwargs['method'],
                }
                output.extend(self._gen_qml_cachegen(state, cachegen_kwargs))

        #copy QML files for Qt tools
        if kwargs['install']:
            self.interpreter.install_data_impl(all_qml_files, module_install_dir,
                                               FileMode(), all_qml_basename, 'devel')

        collected_json: T.Optional[T.Union[FileOrString, build.CustomTarget]] = None
        if kwargs['moc_headers']:
            compile_moc_kwargs: MocCompilerKwArgs = {
                'sources': [],
                'headers': kwargs['moc_headers'],
                'extra_args': kwargs['moc_extra_arguments'],
                'method': kwargs['method'],
                'include_directories': kwargs['include_directories'],
                'dependencies': kwargs['dependencies'],
                'preserve_paths': kwargs['preserve_paths'],
                'output_json': True,
            }
            moc_output = self._compile_moc_impl(state, compile_moc_kwargs)
            output.extend(moc_output)

            moc_collect_json_kwargs: MocJsonCollectKwArgs = {
                'target_name': target_name,
                'moc_json': moc_output,
                'method': kwargs['method'],
            }
            collected_json = self._moc_json_collect(state, moc_collect_json_kwargs)
            output.append(collected_json)

        typeinfo_file: str = ''
        #cmake NO_GENERATE_QMLTYPE disable the whole type registration, not just the .qmltype generation
        if kwargs['generate_qmltype']:
            qmltyperegistrar_kwargs: GenQmlTypeRegistrarKwArgs = {
                'target_name': target_name,
                'import_name': module_name,
                'major_version': module_version_major,
                'minor_version': module_version_minor,
                'collected_json': collected_json,
                'namespace': kwargs['namespace'],
                'generate_qmltype': True,
                'extra_args': kwargs['qmltyperegistrar_extra_arguments'],
                'typeinfo': kwargs['typeinfo'],
                'method': kwargs['method'],
                'install': kwargs['install'],
                'install_dir': module_install_dir,
            }
            type_registrar_output = self._qml_type_registrar(state, qmltyperegistrar_kwargs)
            output.append(type_registrar_output)
            if len(type_registrar_output.get_outputs()) == 2:
                typeinfo_file = type_registrar_output.get_outputs()[1]

        if kwargs['generate_qmldir']:
            qmldir_kwargs: GenQmldirKwArgs = {
                'output': f'{target_name}_qmldir',
                'module_name': module_name,
                'module_version': module_version_short,
                'qml_sources': kwargs['qml_sources'],
                'qml_singletons': kwargs['qml_singletons'],
                'qml_internals': kwargs['qml_internals'],
                'imports': kwargs['imports'],
                'optional_imports': kwargs['optional_imports'],
                'default_imports': kwargs['default_imports'],
                'depends_imports': kwargs['depends_imports'],
                'designer_supported': kwargs['designer_supported'],
                'typeinfo': typeinfo_file,
                'module_prefix': module_prefix_full,
            }
            qmldir_file: File = self._gen_qmldir(state, qmldir_kwargs)

            qmldir_qrc_kwargs: GenQrcKwArgs = {
                'output': f'{target_name}_qmldir.qrc',
                'sources': self._source_to_files(state, [qmldir_file]),
                'aliases': ['qmldir'],
                'prefix': module_prefix_full,
            }
            qrc_resouces.append(self._gen_qrc(state, qmldir_qrc_kwargs))

            if kwargs['install']:
                self.interpreter.install_data_impl([qmldir_file], module_install_dir,
                                                   FileMode(), ['qmldir'], 'devel')

        if qrc_resouces:
            compile_resource_kwargs: ResourceCompilerKwArgs = {
                'name': target_name,
                'sources': qrc_resouces,
                'extra_args': kwargs['rcc_extra_arguments'],
                'method': kwargs['method'],
            }
            output.extend(self._compile_resources_impl(state, compile_resource_kwargs))

        return ModuleReturnValue(output, [output])
