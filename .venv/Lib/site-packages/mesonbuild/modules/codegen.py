# SPDX-License-Identifier: Apache-2.0
# Copyright © 2024-2025 Intel Corporation

from __future__ import annotations
import dataclasses
import os
import typing as T

from . import ExtensionModule, ModuleInfo
from ..build import CustomTarget, CustomTargetIndex, GeneratedList
from ..compilers.compilers import lang_suffixes
from ..interpreter.interpreterobjects import extract_required_kwarg
from ..interpreter.type_checking import NoneType, REQUIRED_KW, DISABLER_KW, NATIVE_KW
from ..interpreterbase import (
    ContainerTypeInfo, ObjectHolder, KwargInfo, typed_pos_args, typed_kwargs,
    noPosargs, noKwargs, disablerIfNotFound, InterpreterObject
)
from ..mesonlib import File, MesonException, Popen_safe, version_compare
from ..programs import Program, ExternalProgram, NonExistingExternalProgram
from ..utils.core import HoldableObject
from .. import mlog

if T.TYPE_CHECKING:
    from typing_extensions import Literal, TypedDict

    from . import ModuleState
    from .._typing import ImmutableListProtocol
    from ..interpreter import Interpreter
    from ..interpreter.kwargs import ExtractRequired
    from ..interpreterbase import TYPE_var, TYPE_kwargs
    from ..mesonlib import MachineChoice
    from ..programs import CommandList

    LexImpls = Literal['lex', 'flex', 'reflex', 'win_flex']
    YaccImpls = Literal['yacc', 'byacc', 'bison', 'win_bison']

    class LexGenerateKwargs(TypedDict):

        args: T.List[str]
        source: T.Optional[str]
        header: T.Optional[str]
        table: T.Optional[str]
        plainname: bool

    class FindLexKwargs(ExtractRequired):

        lex_version: T.List[str]
        flex_version: T.List[str]
        reflex_version: T.List[str]
        win_flex_version: T.List[str]
        implementations: T.List[LexImpls]
        native: MachineChoice

    class YaccGenerateKWargs(TypedDict):

        args: T.List[str]
        source: T.Optional[str]
        header: T.Optional[str]
        locations: T.Optional[str]
        plainname: bool

    class FindYaccKwargs(ExtractRequired):

        yacc_version: T.List[str]
        byacc_version: T.List[str]
        bison_version: T.List[str]
        win_bison_version: T.List[str]
        implementations: T.List[YaccImpls]
        native: MachineChoice


def is_subset_validator(choices: T.Set[str]) -> T.Callable[[T.List[str]], T.Optional[str]]:

    def inner(check: T.List[str]) -> T.Optional[str]:
        if not set(check).issubset(choices):
            invalid = ', '.join(sorted(set(check).difference(choices)))
            valid = ', '.join(sorted(choices))
            return f"valid members are '{valid}', not '{invalid}'"
        return None

    return inner


@dataclasses.dataclass
class _CodeGenerator(HoldableObject):

    name: str
    program: Program
    arguments: ImmutableListProtocol[str] = dataclasses.field(default_factory=list)

    def command(self) -> CommandList:
        return T.cast('CommandList', [self.program]) + T.cast('CommandList', self.arguments)

    def found(self) -> bool:
        return self.program.found()


@dataclasses.dataclass
class LexGenerator(_CodeGenerator):
    pass


class LexHolder(ObjectHolder[LexGenerator]):

    @noPosargs
    @noKwargs
    @InterpreterObject.method('implementation')
    def implementation_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return self.held_object.name

    @noPosargs
    @noKwargs
    @InterpreterObject.method('found')
    def found_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> bool:
        return self.held_object.found()

    @typed_pos_args('codegen.lex.generate', (str, File, GeneratedList, CustomTarget, CustomTargetIndex))
    @typed_kwargs(
        'codegen.lex.generate',
        KwargInfo('args', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo('source', (str, NoneType)),
        KwargInfo('header', (str, NoneType)),
        KwargInfo('table', (str, NoneType)),
        KwargInfo('plainname', bool, default=False),
    )
    @InterpreterObject.method('generate')
    def generate_method(self, args: T.Tuple[T.Union[str, File, GeneratedList, CustomTarget, CustomTargetIndex]], kwargs: LexGenerateKwargs) -> CustomTarget:
        if not self.held_object.found():
            raise MesonException('Attempted to call generate without finding a lex implementation')

        input = self.interpreter.source_strings_to_files([args[0]])[0]
        if isinstance(input, File):
            is_cpp = input.endswith(".ll")
            name = os.path.splitext(input.fname)[0]
        else:
            gen_input = input.get_outputs()
            if len(gen_input) != 1:
                raise MesonException('codegen.lex.generate: generated type inputs must have exactly one output, index into them to select the correct input')
            is_cpp = gen_input[0].endswith('.ll')
            name = os.path.splitext(gen_input[0])[0]
        name = os.path.basename(name)

        # If an explicit source was given, use that to determine whether the
        # user expects this to be a C or C++ source.
        if kwargs['source'] is not None:
            ext = kwargs['source'].rsplit('.', 1)[1]
            is_cpp = ext in lang_suffixes['cpp']

        for_machine = self.held_object.program.for_machine

        # Flex uses FlexLexer.h for C++ code
        if is_cpp and self.held_object.name in {'flex', 'win_flex'}:
            try:
                comp = self.interpreter.environment.coredata.compilers[for_machine]['cpp']
            except KeyError:
                raise MesonException(f"Could not find a C++ compiler for {for_machine} to search for FlexLexer.h")
            found, _ = comp.has_header('FlexLexer.h', '')
            if not found:
                raise MesonException('Could not find FlexLexer.h, which is required for Flex with C++')

        if kwargs['source'] is None:
            outputs = ['@{}@.{}'.format(
                'PLAINNAME' if kwargs['plainname'] else 'BASENAME',
                'cpp' if is_cpp else 'c')]
        else:
            outputs = [kwargs['source']]

        command = self.held_object.command()
        if kwargs['header'] is not None:
            outputs.append(kwargs['header'])
            command.append(f'--header-file=@OUTPUT{len(outputs) - 1}@')
        if kwargs['table'] is not None:
            outputs.append(kwargs['table'])
            command.append(f'--tables-file=@OUTPUT{len(outputs) - 1}@')
        command.extend(kwargs['args'])
        # Flex, at least, seems to require that input be the last argument given
        command.append('@INPUT@')

        target = CustomTarget(
            f'codegen-lex-{name}-{for_machine.get_lower_case_name()}',
            self.interpreter.subdir,
            self.interpreter.subproject,
            self.interpreter.environment,
            command,
            [input],
            outputs,
            backend=self.interpreter.backend,
            description='Generating lexer {{}} with {}'.format(self.held_object.name),
        )
        self.interpreter.add_target(target.name, target)

        return target


@dataclasses.dataclass
class YaccGenerator(_CodeGenerator):
    pass


class YaccHolder(ObjectHolder[YaccGenerator]):

    @noPosargs
    @noKwargs
    @InterpreterObject.method('implementation')
    def implementation_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return self.held_object.name

    @noPosargs
    @noKwargs
    @InterpreterObject.method('found')
    def found_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> bool:
        return self.held_object.found()

    @typed_pos_args('codegen.yacc.generate', (str, File, GeneratedList, CustomTarget, CustomTargetIndex))
    @typed_kwargs(
        'codegen.yacc.generate',
        KwargInfo('args', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo('source', (str, NoneType)),
        KwargInfo('header', (str, NoneType)),
        KwargInfo('locations', (str, NoneType)),
        KwargInfo('plainname', bool, default=False),
    )
    @InterpreterObject.method('generate')
    def generate_method(self, args: T.Tuple[T.Union[str, File, CustomTarget, CustomTargetIndex, GeneratedList]], kwargs: YaccGenerateKWargs) -> CustomTarget:
        if not self.held_object.found():
            raise MesonException('Attempted to call generate without finding a yacc implementation')

        input = self.interpreter.source_strings_to_files([args[0]])[0]
        if isinstance(input, File):
            is_cpp = input.endswith(".yy")
            name = os.path.splitext(input.fname)[0]
        else:
            gen_input = input.get_outputs()
            if len(gen_input) != 1:
                raise MesonException('codegen.lex.generate: generated type inputs must have exactly one output, index into them to select the correct input')
            is_cpp = gen_input[0].endswith('.yy')
            name = os.path.splitext(gen_input[0])[0]
        name = os.path.basename(name)

        command = self.held_object.command()
        command.extend(kwargs['args'])

        source_ext = 'cpp' if is_cpp else 'c'
        header_ext = 'hpp' if is_cpp else 'h'

        base = '@PLAINNAME@' if kwargs['plainname'] else '@BASENAME@'
        outputs: T.List[str] = []
        outputs.append(f'{base}.{source_ext}' if kwargs['source'] is None else kwargs['source'])
        outputs.append(f'{base}.{header_ext}' if kwargs['header'] is None else kwargs['header'])
        if kwargs['locations'] is not None:
            outputs.append(kwargs['locations'])

        for_machine = self.held_object.program.for_machine
        target = CustomTarget(
            f'codegen-yacc-{name}-{for_machine.get_lower_case_name()}',
            self.interpreter.subdir,
            self.interpreter.subproject,
            self.interpreter.environment,
            command,
            [input],
            outputs,
            backend=self.interpreter.backend,
            description='Generating parser {{}} with {}'.format(self.held_object.name),
        )
        self.interpreter.add_target(target.name, target)
        return target


class CodeGenModule(ExtensionModule):

    """Module with helpers for codegen wrappers."""

    INFO = ModuleInfo('codegen', '1.10.0', unstable=True)

    def __init__(self, interpreter: Interpreter) -> None:
        super().__init__(interpreter)
        self.methods.update({
            'lex': self.lex_method,
            'yacc': self.yacc_method,
        })

    @noPosargs
    @typed_kwargs(
        'codegen.lex',
        KwargInfo('lex_version', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo('flex_version', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo('reflex_version', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo('win_flex_version', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo(
            'implementations',
            ContainerTypeInfo(list, str),
            default=[],
            listify=True,
            validator=is_subset_validator({'lex', 'flex', 'reflex', 'win_flex'})
        ),
        REQUIRED_KW,
        DISABLER_KW,
        NATIVE_KW
    )
    @disablerIfNotFound
    def lex_method(self, state: ModuleState, args: T.Tuple, kwargs: FindLexKwargs) -> LexGenerator:
        disabled, required, feature = extract_required_kwarg(kwargs, state.subproject)
        if disabled:
            mlog.log('generator lex skipped: feature', mlog.bold(feature), 'disabled')
            return LexGenerator('lex', NonExistingExternalProgram('lex'))

        names: T.List[LexImpls] = []
        if kwargs['implementations']:
            names = kwargs['implementations']
        else:
            assert state.environment.machines[kwargs['native']] is not None, 'for mypy'
            if state.environment.machines[kwargs['native']].system == 'windows':
                names.append('win_flex')
            names.extend(['flex', 'reflex', 'lex'])

        versions: T.Mapping[str, T.List[str]] = {
            'lex': kwargs['lex_version'],
            'flex': kwargs['flex_version'],
            'reflex': kwargs['reflex_version'],
            'win_flex': kwargs['win_flex_version']
        }

        for name in names:
            bin = state.find_program(
                name, wanted=versions[name], for_machine=kwargs['native'], required=False)
            if bin.found():
                # If you're building reflex as a subproject, we consider that you
                # know what you're doing.
                if name == 'reflex' and isinstance(bin, ExternalProgram):
                    # there are potentially 3 programs called "reflex":
                    # 1. https://invisible-island.net/reflex/, an alternate fork
                    #    of the original flex, this is supported
                    # 2. https://www.genivia.com/doc/reflex/html/, an
                    #    alternative implementation for generating C++ scanners.
                    #    Not supported
                    # 3. https://github.com/cespare/reflex, which is not a lex
                    #    implementation at all, but a file watcher
                    _, out, err = Popen_safe(bin.get_command() + ['--version'])
                    if 'unknown flag: --version' in err:
                        mlog.debug('Skipping cespare/reflex, which is not a lexer and is not supported')
                        continue
                    if 'Written by Robert van Engelen' in out:
                        mlog.debug('Skipping RE/flex, which is not compatible with POSIX lex.')
                        continue
                break
        else:
            if required:
                raise MesonException.from_node(
                    'Could not find a lex implementation. Tried: ', ", ".join(names),
                    node=state.current_node)
            return LexGenerator(name, bin)

        lex_args: T.List[str] = []
        # This option allows compiling with MSVC
        # https://github.com/lexxmark/winflexbison/blob/master/UNISTD_ERROR.readme
        if bin.name == 'win_flex' and state.environment.machines[kwargs['native']].is_windows():
            lex_args.append('--wincompat')
        lex_args.extend(['-o', '@OUTPUT0@'])
        return LexGenerator(name, bin, T.cast('ImmutableListProtocol[str]', lex_args))

    @noPosargs
    @typed_kwargs(
        'codegen.yacc',
        KwargInfo('yacc_version', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo('byacc_version', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo('bison_version', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo('win_bison_version', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo(
            'implementations',
            ContainerTypeInfo(list, str),
            default=[],
            listify=True,
            validator=is_subset_validator({'yacc', 'byacc', 'bison', 'win_bison'})
        ),
        REQUIRED_KW,
        DISABLER_KW,
        NATIVE_KW,
    )
    @disablerIfNotFound
    def yacc_method(self, state: ModuleState, args: T.Tuple, kwargs: FindYaccKwargs) -> YaccGenerator:
        disabled, required, feature = extract_required_kwarg(kwargs, state.subproject)
        if disabled:
            mlog.log('generator yacc skipped: feature', mlog.bold(feature), 'disabled')
            return YaccGenerator('yacc', NonExistingExternalProgram('yacc'))
        names: T.List[YaccImpls]
        if kwargs['implementations']:
            names = kwargs['implementations']
        else:
            assert state.environment.machines[kwargs['native']] is not None, 'for mypy'
            if state.environment.machines[kwargs['native']].system == 'windows':
                names = ['win_bison', 'bison', 'yacc']
            else:
                names = ['bison', 'byacc', 'yacc']

        versions: T.Mapping[YaccImpls, T.List[str]] = {
            'yacc': kwargs['yacc_version'],
            'byacc': kwargs['byacc_version'],
            'bison': kwargs['bison_version'],
            'win_bison': kwargs['win_bison_version'],
        }

        for name in names:
            bin = state.find_program(
                name, wanted=versions[name], for_machine=kwargs['native'], required=False)
            if bin.found():
                break
        else:
            if required:
                raise MesonException.from_node(
                    'Could not find a yacc implementation. Tried: ', ", ".join(names),
                    node=state.current_node)
            return YaccGenerator(name, bin)

        yacc_args: T.List[str] = ['@INPUT@', '-o', '@OUTPUT0@']

        impl = T.cast('YaccImpls', bin.name)
        if impl == 'yacc' and isinstance(bin, ExternalProgram):
            _, out, _ = Popen_safe(bin.get_command() + ['--version'])
            if 'GNU Bison' in out:
                impl = 'bison'
            elif out.startswith('yacc - 2'):
                impl = 'byacc'

        if impl in {'bison', 'win_bison'}:
            yacc_args.append('--defines=@OUTPUT1@')
            if isinstance(bin, ExternalProgram) and version_compare(bin.get_version(), '>= 3.4'):
                yacc_args.append('--color=always')
        elif impl == 'byacc':
            yacc_args.extend(['-H', '@OUTPUT1@'])
        else:
            mlog.warning('This yacc does not appear to be bison or byacc, the '
                         'POSIX specification does not require that header '
                         'output location be configurable, and may not work.',
                         fatal=False)
            yacc_args.append('-H')
        return YaccGenerator(name, bin, T.cast('ImmutableListProtocol[str]', yacc_args))


def initialize(interpreter: Interpreter) -> CodeGenModule:
    interpreter.append_holder_map(LexGenerator, LexHolder)
    interpreter.append_holder_map(YaccGenerator, YaccHolder)
    return CodeGenModule(interpreter)
