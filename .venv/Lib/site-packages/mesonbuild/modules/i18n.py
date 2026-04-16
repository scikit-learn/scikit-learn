# SPDX-License-Identifier: Apache-2.0
# Copyright 2016 The Meson development team

from __future__ import annotations

from os import path
from pathlib import Path
import shlex
import typing as T

from . import ExtensionModule, ModuleReturnValue, ModuleInfo
from .. import build
from .. import mesonlib
from ..options import OptionKey
from .. import mlog
from ..interpreter.primitives import OptionString
from ..interpreter.type_checking import CT_BUILD_BY_DEFAULT, CT_INPUT_KW, INSTALL_TAG_KW, OUTPUT_KW, INSTALL_DIR_KW, INSTALL_KW, NoneType, in_set_validator
from ..interpreterbase import FeatureNew
from ..interpreterbase.exceptions import InvalidArguments
from ..interpreterbase.decorators import ContainerTypeInfo, KwargInfo, noPosargs, typed_kwargs, typed_pos_args
from ..programs import ExternalProgram
from ..scripts.gettext import read_linguas

if T.TYPE_CHECKING:
    from typing_extensions import Literal, TypedDict

    from . import ModuleState
    from ..build import Target
    from ..interpreter import Interpreter
    from ..interpreterbase import TYPE_var
    from ..programs import Program

    class MergeFile(TypedDict):

        input: T.List[T.Union[
            str, build.BuildTarget, build.CustomTarget, build.CustomTargetIndex,
            build.ExtractedObjects, build.GeneratedList, Program,
            mesonlib.File]]
        output: str
        build_by_default: bool
        install: bool
        install_dir: T.Optional[str]
        install_tag: T.Optional[str]
        args: T.List[str]
        data_dirs: T.List[str]
        po_dir: str
        type: Literal['xml', 'desktop']

    class Gettext(TypedDict):

        args: T.List[str]
        data_dirs: T.List[str]
        install: bool
        install_dir: T.Optional[str]
        languages: T.List[str]
        preset: T.Optional[str]

    class ItsJoinFile(TypedDict):

        input: T.List[T.Union[
            str, build.BuildTarget, build.GeneratedTypes,
            build.ExtractedObjects, Program, mesonlib.File]]
        output: str
        build_by_default: bool
        install: bool
        install_dir: T.Optional[str]
        install_tag: T.Optional[str]
        its_files: T.List[str]
        mo_targets: T.List[build.BuildTargetTypes]

    class XgettextProgramT(TypedDict):

        args: T.List[str]
        recursive: bool
        install: bool
        install_dir: T.Optional[str]
        install_tag: T.Optional[str]

    SourcesType = T.Union[str, mesonlib.File, build.BuildTargetTypes, build.BothLibraries]


_ARGS: KwargInfo[T.List[str]] = KwargInfo(
    'args',
    ContainerTypeInfo(list, str),
    default=[],
    listify=True,
)

_DATA_DIRS: KwargInfo[T.List[str]] = KwargInfo(
    'data_dirs',
    ContainerTypeInfo(list, str),
    default=[],
    listify=True
)

PRESET_ARGS = {
    'glib': [
        '--from-code=UTF-8',
        '--add-comments',

        # https://developer.gnome.org/glib/stable/glib-I18N.html
        '--keyword=_',
        '--keyword=N_',
        '--keyword=C_:1c,2',
        '--keyword=NC_:1c,2',
        '--keyword=g_dcgettext:2',
        '--keyword=g_dngettext:2,3',
        '--keyword=g_dpgettext2:2c,3',

        '--flag=N_:1:pass-c-format',
        '--flag=C_:2:pass-c-format',
        '--flag=NC_:2:pass-c-format',
        '--flag=g_dngettext:2:pass-c-format',
        '--flag=g_strdup_printf:1:c-format',
        '--flag=g_string_printf:2:c-format',
        '--flag=g_string_append_printf:2:c-format',
        '--flag=g_error_new:3:c-format',
        '--flag=g_set_error:4:c-format',
        '--flag=g_markup_printf_escaped:1:c-format',
        '--flag=g_log:3:c-format',
        '--flag=g_print:1:c-format',
        '--flag=g_printerr:1:c-format',
        '--flag=g_printf:1:c-format',
        '--flag=g_fprintf:2:c-format',
        '--flag=g_sprintf:2:c-format',
        '--flag=g_snprintf:3:c-format',
    ]
}


class XgettextProgram:

    pot_files: T.Dict[str, build.CustomTarget] = {}

    def __init__(self, xgettext: ExternalProgram, interpreter: Interpreter):
        self.xgettext = xgettext
        self.interpreter = interpreter

    def extract(self,
                name: str,
                sources: T.List[SourcesType],
                args: T.List[str],
                recursive: bool,
                install: bool,
                install_dir: T.Optional[str],
                install_tag: T.Optional[str]) -> build.CustomTarget:

        if not name.endswith('.pot'):
            name += '.pot'

        source_files = self._get_source_files(sources)

        command = self.xgettext.command + args
        command.append(f'--directory={self.interpreter.environment.get_source_dir()}')
        command.append(f'--directory={self.interpreter.environment.get_build_dir()}')
        command.append('--output=@OUTPUT@')

        depends = list(self._get_depends(sources)) if recursive else []
        rsp_file = self._get_rsp_file(name, source_files, depends, command)
        inputs: T.List[T.Union[mesonlib.File, build.CustomTarget]]
        if rsp_file:
            inputs = [rsp_file]
            depend_files = list(source_files)
            command.append('--files-from=@INPUT@')
        else:
            inputs = list(source_files) + depends
            depends = None
            depend_files = None
            command.append('@INPUT@')

        ct = build.CustomTarget(
            '',
            self.interpreter.subdir,
            self.interpreter.subproject,
            self.interpreter.environment,
            command,
            inputs,
            [name],
            depend_files = depend_files,
            extra_depends = depends,
            install = install,
            install_dir = [install_dir] if install_dir else None,
            install_tag = [install_tag] if install_tag else None,
            description = 'Extracting translations to {}',
        )

        for source_id in self._get_source_id(sources):
            self.pot_files[source_id] = ct
        self.pot_files[ct.get_id()] = ct

        self.interpreter.add_target(ct.name, ct)
        return ct

    def _get_source_files(self, sources: T.Iterable[SourcesType]) -> T.Set[mesonlib.File]:
        source_files = set()
        for source in sources:
            if isinstance(source, mesonlib.File):
                source_files.add(source)
            elif isinstance(source, str):
                mesonlib.check_direntry_issues(source)
                source_files.add(mesonlib.File.from_source_file(self.interpreter.source_root, self.interpreter.subdir, source))
            elif isinstance(source, build.BuildTarget):
                source_files.update(source.get_sources())
            elif isinstance(source, build.BothLibraries):
                source_files.update(source.get('shared').get_sources())
            elif isinstance(source, (build.CustomTarget, build.CustomTargetIndex)):
                source_files.update(mesonlib.File.from_built_file(source.get_subdir(), f) for f in source.get_outputs())
        return source_files

    def _get_depends(self, sources: T.Iterable[SourcesType]) -> T.Set[build.CustomTarget]:
        depends = set()
        for source in sources:
            if isinstance(source, build.BuildTarget):
                for source_id in self._get_source_id(source.get_dependencies()):
                    if source_id in self.pot_files:
                        depends.add(self.pot_files[source_id])
            elif isinstance(source, build.CustomTarget):
                # Dependency on another extracted pot file
                source_id = source.get_id()
                if source_id in self.pot_files:
                    depends.add(self.pot_files[source_id])
        return depends

    def _get_rsp_file(self,
                      name: str,
                      source_files: T.Iterable[mesonlib.File],
                      depends: T.Iterable[build.CustomTarget],
                      arguments: T.List[str]) -> T.Optional[mesonlib.File]:
        source_list = '\n'.join(source.relative_name() for source in source_files)
        for dep in depends:
            source_list += '\n' + path.join(dep.subdir, dep.get_filename())

        estimated_cmdline_length = len(source_list) + sum(len(arg) + 1 for arg in arguments) + 1
        if estimated_cmdline_length < mesonlib.get_rsp_threshold():
            return None

        rsp_file = Path(self.interpreter.environment.build_dir, self.interpreter.subdir, name+'.rsp')
        rsp_file.write_text(source_list, encoding='utf-8')

        return mesonlib.File.from_built_file(self.interpreter.subdir, rsp_file.name)

    @staticmethod
    def _get_source_id(sources: T.Iterable[SourcesType]) -> T.Iterable[str]:
        for source in sources:
            if isinstance(source, build.Target):
                yield source.get_id()
            elif isinstance(source, build.BothLibraries):
                yield source.get('static').get_id()
                yield source.get('shared').get_id()


class I18nModule(ExtensionModule):

    INFO = ModuleInfo('i18n')

    def __init__(self, interpreter: 'Interpreter'):
        super().__init__(interpreter)
        self.methods.update({
            'merge_file': self.merge_file,
            'gettext': self.gettext,
            'itstool_join': self.itstool_join,
            'xgettext': self.xgettext,
        })
        self.tools: T.Dict[str, T.Optional[Program]] = {
            'itstool': None,
            'msgfmt': None,
            'msginit': None,
            'msgmerge': None,
            'xgettext': None,
        }

    @staticmethod
    def _get_data_dirs(state: 'ModuleState', dirs: T.Iterable[str]) -> T.List[str]:
        """Returns source directories of relative paths"""
        src_dir = path.join(state.environment.get_source_dir(), state.subdir)
        return [path.join(src_dir, d) for d in dirs]

    @FeatureNew('i18n.merge_file', '0.37.0')
    @noPosargs
    @typed_kwargs(
        'i18n.merge_file',
        CT_BUILD_BY_DEFAULT,
        CT_INPUT_KW,
        KwargInfo('install_dir', (str, NoneType)),
        INSTALL_TAG_KW,
        OUTPUT_KW,
        INSTALL_KW,
        _ARGS.evolve(since='0.51.0'),
        _DATA_DIRS.evolve(since='0.41.0'),
        KwargInfo('po_dir', str, required=True),
        KwargInfo('type', str, default='xml', validator=in_set_validator({'xml', 'desktop'})),
    )
    def merge_file(self, state: 'ModuleState', args: T.List['TYPE_var'], kwargs: 'MergeFile') -> ModuleReturnValue:
        if kwargs['install'] and not kwargs['install_dir']:
            raise InvalidArguments('i18n.merge_file: "install_dir" keyword argument must be set when "install" is true.')

        if self.tools['msgfmt'] is None or not self.tools['msgfmt'].found():
            self.tools['msgfmt'] = state.find_program('msgfmt', for_machine=mesonlib.MachineChoice.BUILD)
        if isinstance(self.tools['msgfmt'], ExternalProgram):
            try:
                have_version = self.tools['msgfmt'].get_version()
            except mesonlib.MesonException as e:
                raise mesonlib.MesonException('i18n.merge_file requires GNU msgfmt') from e
            want_version = '>=0.19' if kwargs['type'] == 'desktop' else '>=0.19.7'
            if not mesonlib.version_compare(have_version, want_version):
                msg = f'i18n.merge_file requires GNU msgfmt {want_version} to produce files of type: ' + kwargs['type'] + f' (got: {have_version})'
                raise mesonlib.MesonException(msg)
        podir = path.join(state.build_to_src, state.subdir, kwargs['po_dir'])

        ddirs = self._get_data_dirs(state, kwargs['data_dirs'])
        datadirs = '--datadirs=' + ':'.join(ddirs) if ddirs else None

        command: T.List[T.Union[str, build.BuildTargetTypes, Program, mesonlib.File]] = []
        command.extend(state.environment.get_build_command())
        command.extend([
            '--internal', 'msgfmthelper',
            '--msgfmt=' + self.tools['msgfmt'].get_path(),
        ])
        if datadirs:
            command.append(datadirs)
        command.extend(['@INPUT@', '@OUTPUT@', kwargs['type'], podir])
        if kwargs['args']:
            command.append('--')
            command.extend(kwargs['args'])

        build_by_default = kwargs['build_by_default']
        if build_by_default is None:
            build_by_default = kwargs['install']

        install_tag = [kwargs['install_tag']] if kwargs['install_tag'] is not None else None

        ct = build.CustomTarget(
            '',
            state.subdir,
            state.subproject,
            state.environment,
            command,
            kwargs['input'],
            [kwargs['output']],
            build_by_default=build_by_default,
            install=kwargs['install'],
            install_dir=[kwargs['install_dir']] if kwargs['install_dir'] is not None else None,
            install_tag=install_tag,
            description='Merging translations for {}',
        )

        return ModuleReturnValue(ct, [ct])

    @typed_pos_args('i18n.gettext', str)
    @typed_kwargs(
        'i18n.gettext',
        _ARGS,
        _DATA_DIRS.evolve(since='0.36.0'),
        INSTALL_KW.evolve(default=True),
        INSTALL_DIR_KW.evolve(since='0.50.0'),
        KwargInfo('languages', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo(
            'preset',
            (str, NoneType),
            validator=in_set_validator(set(PRESET_ARGS)),
            since='0.37.0',
        ),
    )
    def gettext(self, state: 'ModuleState', args: T.Tuple[str], kwargs: 'Gettext') -> ModuleReturnValue:
        for tool, strict in [('msgfmt', True), ('msginit', False), ('msgmerge', False), ('xgettext', False)]:
            if self.tools[tool] is None:
                self.tools[tool] = state.find_program(tool, required=False, for_machine=mesonlib.MachineChoice.BUILD)
            # still not found?
            if not self.tools[tool].found():
                if strict:
                    mlog.warning('Gettext not found, all translation (po) targets will be ignored.',
                                 once=True, location=state.current_node)
                    return ModuleReturnValue(None, [])
                else:
                    mlog.warning(f'{tool!r} not found, maintainer targets will not work',
                                 once=True, fatal=False, location=state.current_node)
        packagename = args[0]
        pkg_arg = f'--pkgname={packagename}'

        languages = kwargs['languages']
        lang_arg = '--langs=' + '@@'.join(languages) if languages else None

        _datadirs = ':'.join(self._get_data_dirs(state, kwargs['data_dirs']))
        datadirs = f'--datadirs={_datadirs}' if _datadirs else None

        extra_args = kwargs['args']
        targets: T.List['Target'] = []
        gmotargets: T.List['build.CustomTarget'] = []

        preset = kwargs['preset']
        if preset:
            preset_args = PRESET_ARGS[preset]
            extra_args = list(mesonlib.OrderedSet(preset_args + extra_args))

        extra_arg = '--extra-args=' + '@@'.join(extra_args) if extra_args else None

        source_root = path.join(state.source_root, state.root_subdir)
        subdir = path.relpath(state.subdir, start=state.root_subdir) if state.subdir else None

        potargs = state.environment.get_build_command() + ['--internal', 'gettext', 'pot', pkg_arg]
        potargs.append(f'--source-root={source_root}')
        if subdir:
            potargs.append(f'--subdir={subdir}')
        if datadirs:
            potargs.append(datadirs)
        if extra_arg:
            potargs.append(extra_arg)
        if self.tools['xgettext'].found():
            potargs.append('--xgettext=' + self.tools['xgettext'].get_path())
        pottarget = build.RunTarget(packagename + '-pot', potargs, [], state.subdir, state.subproject,
                                    state.environment, default_env=False)
        targets.append(pottarget)

        install = kwargs['install']
        install_dir = kwargs['install_dir'] or state.environment.coredata.optstore.get_value_for(OptionKey('localedir'))
        assert isinstance(install_dir, str), 'for mypy'
        if not languages:
            languages = read_linguas(path.join(state.environment.source_dir, state.subdir))
        for l in languages:
            po_file = mesonlib.File.from_source_file(state.environment.source_dir, state.subdir, l+'.po')
            mo_install_dir = path.join(install_dir, l, 'LC_MESSAGES')
            if isinstance(install_dir, OptionString):
                name = path.join(install_dir.optname, l, 'LC_MESSAGES')
                mo_install_dir = OptionString(mo_install_dir, name)
            gmotarget = build.CustomTarget(
                f'{packagename}-{l}.mo',
                path.join(state.subdir, l, 'LC_MESSAGES'),
                state.subproject,
                state.environment,
                [self.tools['msgfmt'], '-o', '@OUTPUT@', '@INPUT@'],
                [po_file],
                [f'{packagename}.mo'],
                install=install,
                # We have multiple files all installed as packagename+'.mo' in different install subdirs.
                # What we really wanted to do, probably, is have a rename: kwarg, but that's not available
                # to custom_targets. Crude hack: set the build target's subdir manually.
                # Bonus: the build tree has something usable as an uninstalled bindtextdomain() target dir.
                install_dir=[mo_install_dir],
                install_tag=['i18n'],
                description='Building translation {}',
            )
            targets.append(gmotarget)
            gmotargets.append(gmotarget)

        allgmotarget = build.AliasTarget(packagename + '-gmo', gmotargets, state.subdir, state.subproject,
                                         state.environment)
        targets.append(allgmotarget)

        updatepoargs = state.environment.get_build_command() + ['--internal', 'gettext', 'update_po', pkg_arg]
        updatepoargs.append(f'--source-root={source_root}')
        if subdir:
            updatepoargs.append(f'--subdir={subdir}')
        if lang_arg:
            updatepoargs.append(lang_arg)
        if datadirs:
            updatepoargs.append(datadirs)
        if extra_arg:
            updatepoargs.append(extra_arg)
        for tool in ['msginit', 'msgmerge']:
            if self.tools[tool].found():
                updatepoargs.append(f'--{tool}=' + self.tools[tool].get_path())
        updatepotarget = build.RunTarget(packagename + '-update-po', updatepoargs, [], state.subdir, state.subproject,
                                         state.environment, default_env=False)
        targets.append(updatepotarget)

        return ModuleReturnValue([gmotargets, pottarget, updatepotarget], targets)

    @FeatureNew('i18n.itstool_join', '0.62.0')
    @noPosargs
    @typed_kwargs(
        'i18n.itstool_join',
        CT_BUILD_BY_DEFAULT,
        CT_INPUT_KW,
        KwargInfo('install_dir', (str, NoneType)),
        INSTALL_TAG_KW,
        OUTPUT_KW,
        INSTALL_KW,
        _ARGS.evolve(),
        KwargInfo('its_files', ContainerTypeInfo(list, str)),
        KwargInfo('mo_targets', ContainerTypeInfo(list, build.CustomTarget), required=True),
    )
    def itstool_join(self, state: 'ModuleState', args: T.List['TYPE_var'], kwargs: 'ItsJoinFile') -> ModuleReturnValue:
        if kwargs['install'] and not kwargs['install_dir']:
            raise InvalidArguments('i18n.itstool_join: "install_dir" keyword argument must be set when "install" is true.')

        if self.tools['itstool'] is None:
            self.tools['itstool'] = state.find_program('itstool', for_machine=mesonlib.MachineChoice.BUILD)
        mo_targets = kwargs['mo_targets']
        its_files = kwargs.get('its_files', [])

        mo_fnames = []
        for target in mo_targets:
            mo_fnames.append(path.join(target.get_builddir(), target.get_outputs()[0]))

        command: T.List[T.Union[str, build.BuildTargetTypes, Program, mesonlib.File]] = []
        command.extend(state.environment.get_build_command())

        itstool_cmd = self.tools['itstool'].get_command()
        # TODO: python 3.8 can use shlex.join()
        command.extend([
            '--internal', 'itstool', 'join',
            '-i', '@INPUT@',
            '-o', '@OUTPUT@',
            '--itstool=' + ' '.join(shlex.quote(c) for c in itstool_cmd),
        ])
        if its_files:
            for fname in its_files:
                if not path.isabs(fname):
                    fname = path.join(state.environment.source_dir, state.subdir, fname)
                command.extend(['--its', fname])
        command.extend(mo_fnames)

        build_by_default = kwargs['build_by_default']
        if build_by_default is None:
            build_by_default = kwargs['install']

        install_tag = [kwargs['install_tag']] if kwargs['install_tag'] is not None else None

        ct = build.CustomTarget(
            '',
            state.subdir,
            state.subproject,
            state.environment,
            command,
            kwargs['input'],
            [kwargs['output']],
            build_by_default=build_by_default,
            extra_depends=mo_targets,
            install=kwargs['install'],
            install_dir=[kwargs['install_dir']] if kwargs['install_dir'] is not None else None,
            install_tag=install_tag,
            description='Merging translations for {}',
        )

        return ModuleReturnValue(ct, [ct])

    @FeatureNew('i18n.xgettext', '1.8.0')
    @typed_pos_args('i18n.xgettext', str, varargs=(str, mesonlib.File, build.BuildTarget, build.BothLibraries, build.CustomTarget, build.CustomTargetIndex), min_varargs=1)
    @typed_kwargs(
        'i18n.xgettext',
        _ARGS,
        KwargInfo('recursive', bool, default=False),
        INSTALL_KW,
        INSTALL_DIR_KW,
        INSTALL_TAG_KW,
    )
    def xgettext(self, state: ModuleState, args: T.Tuple[str, T.List[SourcesType]], kwargs: XgettextProgramT) -> build.CustomTarget:
        if any(isinstance(a, build.CustomTarget) for a in args[1]):
            FeatureNew.single_use('i18n.xgettext with custom_target is broken until 1.10', '1.10.0', self.interpreter.subproject, location=self.interpreter.current_node)
        if any(isinstance(a, build.CustomTargetIndex) for a in args[1]):
            FeatureNew.single_use('i18n.xgettext with custom_target index', '1.10.0', self.interpreter.subproject, location=self.interpreter.current_node)

        toolname = 'xgettext'
        if self.tools[toolname] is None or not self.tools[toolname].found():
            self.tools[toolname] = state.find_program(toolname, required=True, for_machine=mesonlib.MachineChoice.BUILD)

        if kwargs['install'] and not kwargs['install_dir']:
            raise InvalidArguments('i18n.xgettext: "install_dir" keyword argument must be set when "install" is true.')

        xgettext_program = XgettextProgram(T.cast('ExternalProgram', self.tools[toolname]), self.interpreter)
        return xgettext_program.extract(*args, **kwargs)


def initialize(interp: 'Interpreter') -> I18nModule:
    return I18nModule(interp)
