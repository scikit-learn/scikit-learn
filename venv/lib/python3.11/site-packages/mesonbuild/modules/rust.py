# SPDX-License-Identifier: Apache-2.0
# Copyright © 2020-2024 Intel Corporation

from __future__ import annotations
import itertools
import os
import typing as T

from mesonbuild.interpreterbase.decorators import FeatureNew

from . import ExtensionModule, ModuleReturnValue, ModuleInfo
from .. import mesonlib, mlog
from ..build import (BothLibraries, BuildTarget, CustomTargetIndex, Executable, ExtractedObjects, GeneratedList,
                     CustomTarget, InvalidArguments, Jar, StructuredSources, SharedLibrary)
from ..compilers.compilers import are_asserts_disabled, lang_suffixes
from ..interpreter.type_checking import (
    DEPENDENCIES_KW, LINK_WITH_KW, SHARED_LIB_KWS, TEST_KWS, OUTPUT_KW,
    INCLUDE_DIRECTORIES, SOURCES_VARARGS, NoneType, in_set_validator
)
from ..interpreterbase import ContainerTypeInfo, InterpreterException, KwargInfo, typed_kwargs, typed_pos_args, noPosargs, permittedKwargs
from ..mesonlib import File
from ..programs import ExternalProgram

if T.TYPE_CHECKING:
    from . import ModuleState
    from ..build import IncludeDirs, LibTypes
    from ..dependencies import Dependency, ExternalLibrary
    from ..interpreter import Interpreter
    from ..interpreter import kwargs as _kwargs
    from ..interpreter.interpreter import SourceInputs, SourceOutputs
    from ..programs import OverrideProgram
    from ..interpreter.type_checking import SourcesVarargsType

    from typing_extensions import TypedDict, Literal

    class FuncTest(_kwargs.BaseTest):

        dependencies: T.List[T.Union[Dependency, ExternalLibrary]]
        is_parallel: bool
        link_with: T.List[LibTypes]
        rust_args: T.List[str]

    class FuncBindgen(TypedDict):

        args: T.List[str]
        c_args: T.List[str]
        include_directories: T.List[IncludeDirs]
        input: T.List[SourceInputs]
        output: str
        output_inline_wrapper: str
        dependencies: T.List[T.Union[Dependency, ExternalLibrary]]
        language: T.Optional[Literal['c', 'cpp']]
        bindgen_version: T.List[str]


class RustModule(ExtensionModule):

    """A module that holds helper functions for rust."""

    INFO = ModuleInfo('rust', '0.57.0', stabilized='1.0.0')

    def __init__(self, interpreter: Interpreter) -> None:
        super().__init__(interpreter)
        self._bindgen_bin: T.Optional[T.Union[ExternalProgram, Executable, OverrideProgram]] = None
        if 'rust' in interpreter.compilers.host:
            self._bindgen_rust_target: T.Optional[str] = interpreter.compilers.host['rust'].version
        else:
            self._bindgen_rust_target = None
        self.methods.update({
            'test': self.test,
            'bindgen': self.bindgen,
            'proc_macro': self.proc_macro,
        })

    @typed_pos_args('rust.test', str, BuildTarget)
    @typed_kwargs(
        'rust.test',
        *TEST_KWS,
        DEPENDENCIES_KW,
        LINK_WITH_KW.evolve(since='1.2.0'),
        KwargInfo(
            'rust_args',
            ContainerTypeInfo(list, str),
            listify=True,
            default=[],
            since='1.2.0',
        ),
        KwargInfo('is_parallel', bool, default=False),
    )
    def test(self, state: ModuleState, args: T.Tuple[str, BuildTarget], kwargs: FuncTest) -> ModuleReturnValue:
        """Generate a rust test target from a given rust target.

        Rust puts its unitests inside its main source files, unlike most
        languages that put them in external files. This means that normally
        you have to define two separate targets with basically the same
        arguments to get tests:

        ```meson
        rust_lib_sources = [...]
        rust_lib = static_library(
            'rust_lib',
            rust_lib_sources,
        )

        rust_lib_test = executable(
            'rust_lib_test',
            rust_lib_sources,
            rust_args : ['--test'],
        )

        test(
            'rust_lib_test',
            rust_lib_test,
            protocol : 'rust',
        )
        ```

        This is all fine, but not very DRY. This method makes it much easier
        to define rust tests:

        ```meson
        rust = import('unstable-rust')

        rust_lib = static_library(
            'rust_lib',
            [sources],
        )

        rust.test('rust_lib_test', rust_lib)
        ```
        """
        if any(isinstance(t, Jar) for t in kwargs.get('link_with', [])):
            raise InvalidArguments('Rust tests cannot link with Jar targets')

        name = args[0]
        base_target: BuildTarget = args[1]
        if not base_target.uses_rust():
            raise InterpreterException('Second positional argument to rustmod.test() must be a rust based target')
        extra_args = kwargs['args']

        # Delete any arguments we don't want passed
        if '--test' in extra_args:
            mlog.warning('Do not add --test to rustmod.test arguments')
            extra_args.remove('--test')
        if '--format' in extra_args:
            mlog.warning('Do not add --format to rustmod.test arguments')
            i = extra_args.index('--format')
            # Also delete the argument to --format
            del extra_args[i + 1]
            del extra_args[i]
        for i, a in enumerate(extra_args):
            if isinstance(a, str) and a.startswith('--format='):
                del extra_args[i]
                break

        # We need to cast here, as currently these don't have protocol in them, but test itself does.
        tkwargs = T.cast('_kwargs.FuncTest', kwargs.copy())

        tkwargs['args'] = extra_args + ['--test', '--format', 'pretty']
        tkwargs['protocol'] = 'rust'

        new_target_kwargs = base_target.original_kwargs.copy()
        # Don't mutate the shallow copied list, instead replace it with a new
        # one
        new_target_kwargs['install'] = False
        new_target_kwargs['dependencies'] = new_target_kwargs.get('dependencies', []) + kwargs['dependencies']
        new_target_kwargs['link_with'] = new_target_kwargs.get('link_with', []) + kwargs['link_with']
        del new_target_kwargs['rust_crate_type']
        for kw in ['pic', 'prelink', 'rust_abi', 'version', 'soversion', 'darwin_versions']:
            if kw in new_target_kwargs:
                del new_target_kwargs[kw]

        lang_args = base_target.extra_args.copy()
        lang_args['rust'] = base_target.extra_args['rust'] + kwargs['rust_args'] + ['--test']
        new_target_kwargs['language_args'] = lang_args

        sources = T.cast('T.List[SourceOutputs]', base_target.sources.copy())
        sources.extend(base_target.generated)

        new_target = Executable(
            name, base_target.subdir, state.subproject, base_target.for_machine,
            sources, base_target.structured_sources,
            base_target.objects, base_target.environment, base_target.compilers,
            new_target_kwargs
        )

        test = self.interpreter.make_test(
            self.interpreter.current_node, (name, new_target), tkwargs)

        return ModuleReturnValue(None, [new_target, test])

    @noPosargs
    @typed_kwargs(
        'rust.bindgen',
        KwargInfo('c_args', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo('args', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo(
            'input',
            ContainerTypeInfo(list, (File, GeneratedList, BuildTarget, BothLibraries, ExtractedObjects, CustomTargetIndex, CustomTarget, str), allow_empty=False),
            default=[],
            listify=True,
            required=True,
        ),
        KwargInfo('language', (str, NoneType), since='1.4.0', validator=in_set_validator({'c', 'cpp'})),
        KwargInfo('bindgen_version', ContainerTypeInfo(list, str), default=[], listify=True, since='1.4.0'),
        INCLUDE_DIRECTORIES.evolve(since_values={ContainerTypeInfo(list, str): '1.0.0'}),
        OUTPUT_KW,
        KwargInfo(
            'output_inline_wrapper',
            str,
            default='',
            since='1.4.0',
        ),
        DEPENDENCIES_KW.evolve(since='1.0.0'),
    )
    def bindgen(self, state: ModuleState, args: T.List, kwargs: FuncBindgen) -> ModuleReturnValue:
        """Wrapper around bindgen to simplify its use.

        The main thing this simplifies is the use of `include_directory`
        objects, instead of having to pass a plethora of `-I` arguments.
        """
        header, *_deps = self.interpreter.source_strings_to_files(kwargs['input'])

        # Split File and Target dependencies to add pass to CustomTarget
        depends: T.List[SourceOutputs] = []
        depend_files: T.List[File] = []
        for d in _deps:
            if isinstance(d, File):
                depend_files.append(d)
            else:
                depends.append(d)

        # Copy to avoid subsequent calls mutating the original
        # TODO: if we want this to be per-machine we'll need a native kwarg
        clang_args = state.environment.properties.host.get_bindgen_clang_args().copy()

        for i in state.process_include_dirs(kwargs['include_directories']):
            # bindgen always uses clang, so it's safe to hardcode -I here
            clang_args.extend([f'-I{x}' for x in i.to_string_list(
                state.environment.get_source_dir(), state.environment.get_build_dir())])
        if are_asserts_disabled(state.environment.coredata.optstore):
            clang_args.append('-DNDEBUG')

        for de in kwargs['dependencies']:
            for i in de.get_include_dirs():
                clang_args.extend([f'-I{x}' for x in i.to_string_list(
                    state.environment.get_source_dir(), state.environment.get_build_dir())])
            clang_args.extend(de.get_all_compile_args())
            for s in de.get_sources():
                if isinstance(s, File):
                    depend_files.append(s)
                elif isinstance(s, CustomTarget):
                    depends.append(s)

        if self._bindgen_bin is None:
            self._bindgen_bin = state.find_program('bindgen', wanted=kwargs['bindgen_version'])
            if self._bindgen_rust_target is not None:
                # ExternalCommand.command's type is bonkers
                _, _, err = mesonlib.Popen_safe(
                    T.cast('T.List[str]', self._bindgen_bin.get_command()) +
                    ['--rust-target', self._bindgen_rust_target])
                # Sometimes this is "invalid Rust target" and sometimes "invalid
                # rust target"
                if 'Got an invalid' in err:
                    self._bindgen_rust_target = None

        name: str
        if isinstance(header, File):
            name = header.fname
        elif isinstance(header, (BuildTarget, BothLibraries, ExtractedObjects, StructuredSources)):
            raise InterpreterException('bindgen source file must be a C header, not an object or build target')
        else:
            name = header.get_outputs()[0]

        # bindgen assumes that C++ headers will be called .hpp. We want to
        # ensure that anything Meson considers a C++ header is treated as one.
        language = kwargs['language']
        if language is None:
            ext = os.path.splitext(name)[1][1:]
            if ext in lang_suffixes['cpp']:
                language = 'cpp'
            elif ext == 'h':
                language = 'c'
            else:
                raise InterpreterException(f'Unknown file type extension for: {name}')

        # We only want include directories and defines, other things may not be valid
        cargs = state.get_option(f'{language}_args', state.subproject)
        assert isinstance(cargs, list), 'for mypy'
        for a in itertools.chain(state.global_args.get(language, []), state.project_args.get(language, []), cargs):
            if a.startswith(('-I', '/I', '-D', '/D', '-U', '/U')):
                clang_args.append(a)

        if language == 'cpp':
            clang_args.extend(['-x', 'c++'])

        # Add the C++ standard to the clang arguments. Attempt to translate VS
        # extension versions into the nearest standard version
        std = state.get_option(f'{language}_std')
        assert isinstance(std, str), 'for mypy'
        if std.startswith('vc++'):
            if std.endswith('latest'):
                mlog.warning('Attempting to translate vc++latest into a clang compatible version.',
                             'Currently this is hardcoded for c++20', once=True, fatal=False)
                std = 'c++20'
            else:
                mlog.debug('The current C++ standard is a Visual Studio extension version.',
                           'bindgen will use a the nearest C++ standard instead')
                std = std[1:]

        if std != 'none':
            clang_args.append(f'-std={std}')

        inline_wrapper_args: T.List[str] = []
        outputs = [kwargs['output']]
        if kwargs['output_inline_wrapper']:
            # Todo drop this isinstance once Executable supports version_compare
            if isinstance(self._bindgen_bin, ExternalProgram):
                if mesonlib.version_compare(self._bindgen_bin.get_version(), '< 0.65'):
                    raise InterpreterException('\'output_inline_wrapper\' parameter of rust.bindgen requires bindgen-0.65 or newer')

            outputs.append(kwargs['output_inline_wrapper'])
            inline_wrapper_args = [
                '--experimental', '--wrap-static-fns',
                '--wrap-static-fns-path', os.path.join(state.environment.build_dir, '@OUTPUT1@')
            ]

        cmd = self._bindgen_bin.get_command() + \
            [
                '@INPUT@', '--output',
                os.path.join(state.environment.build_dir, '@OUTPUT0@')
            ] + \
            kwargs['args'] + inline_wrapper_args
        if self._bindgen_rust_target and '--rust-target' not in cmd:
            cmd.extend(['--rust-target', self._bindgen_rust_target])
        cmd.append('--')
        cmd.extend(kwargs['c_args'])
        cmd.extend(clang_args)
        cmd.extend(['-MD', '-MQ', '@INPUT@', '-MF', '@DEPFILE@'])

        target = CustomTarget(
            f'rustmod-bindgen-{name}'.replace('/', '_'),
            state.subdir,
            state.subproject,
            state.environment,
            cmd,
            [header],
            outputs,
            depfile='@PLAINNAME@.d',
            extra_depends=depends,
            depend_files=depend_files,
            backend=state.backend,
            description='Generating bindings for Rust {}',
        )

        return ModuleReturnValue(target, [target])

    # Allow a limited set of kwargs, but still use the full set of typed_kwargs()
    # because it could be setting required default values.
    @FeatureNew('rust.proc_macro', '1.3.0')
    @permittedKwargs({'rust_args', 'rust_dependency_map', 'sources', 'dependencies', 'extra_files',
                      'link_args', 'link_depends', 'link_with', 'override_options'})
    @typed_pos_args('rust.proc_macro', str, varargs=SOURCES_VARARGS)
    @typed_kwargs('rust.proc_macro', *SHARED_LIB_KWS, allow_unknown=True)
    def proc_macro(self, state: ModuleState, args: T.Tuple[str, SourcesVarargsType], kwargs: _kwargs.SharedLibrary) -> SharedLibrary:
        kwargs['native'] = True  # type: ignore
        kwargs['rust_crate_type'] = 'proc-macro'  # type: ignore
        kwargs['rust_args'] = kwargs['rust_args'] + ['--extern', 'proc_macro']
        target = state._interpreter.build_target(state.current_node, args, kwargs, SharedLibrary)
        return target


def initialize(interp: Interpreter) -> RustModule:
    return RustModule(interp)
