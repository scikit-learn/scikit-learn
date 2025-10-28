# SPDX-License-Identifier: Apache-2.0
# Copyright 2021 The Meson development team
from __future__ import annotations

import re
import os

import typing as T

from ...mesonlib import version_compare, version_compare_many
from ...interpreterbase import (
    InterpreterObject,
    MesonOperator,
    ObjectHolder,
    FeatureNew,
    typed_operator,
    noArgsFlattening,
    noKwargs,
    noPosargs,
    typed_pos_args,
    InvalidArguments,
    FeatureBroken,
    stringifyUserArguments,
)


if T.TYPE_CHECKING:
    from ...interpreterbase import TYPE_var, TYPE_kwargs

class StringHolder(ObjectHolder[str]):
    TRIVIAL_OPERATORS = {
        # Arithmetic
        MesonOperator.PLUS: (str, lambda obj, x: obj.held_object + x),

        # Comparison
        MesonOperator.EQUALS: (str, lambda obj, x: obj.held_object == x),
        MesonOperator.NOT_EQUALS: (str, lambda obj, x: obj.held_object != x),
        MesonOperator.GREATER: (str, lambda obj, x: obj.held_object > x),
        MesonOperator.LESS: (str, lambda obj, x: obj.held_object < x),
        MesonOperator.GREATER_EQUALS: (str, lambda obj, x: obj.held_object >= x),
        MesonOperator.LESS_EQUALS: (str, lambda obj, x: obj.held_object <= x),
    }

    def display_name(self) -> str:
        return 'str'

    @noKwargs
    @typed_pos_args('str.contains', str)
    @InterpreterObject.method('contains')
    def contains_method(self, args: T.Tuple[str], kwargs: TYPE_kwargs) -> bool:
        return self.held_object.find(args[0]) >= 0

    @noKwargs
    @typed_pos_args('str.startswith', str)
    @InterpreterObject.method('startswith')
    def startswith_method(self, args: T.Tuple[str], kwargs: TYPE_kwargs) -> bool:
        return self.held_object.startswith(args[0])

    @noKwargs
    @typed_pos_args('str.endswith', str)
    @InterpreterObject.method('endswith')
    def endswith_method(self, args: T.Tuple[str], kwargs: TYPE_kwargs) -> bool:
        return self.held_object.endswith(args[0])

    @noArgsFlattening
    @noKwargs
    @typed_pos_args('str.format', varargs=object)
    @InterpreterObject.method('format')
    def format_method(self, args: T.Tuple[T.List[TYPE_var]], kwargs: TYPE_kwargs) -> str:
        arg_strings: T.List[str] = []
        for arg in args[0]:
            try:
                arg_strings.append(stringifyUserArguments(arg, self.subproject))
            except InvalidArguments as e:
                FeatureBroken.single_use(f'str.format: {str(e)}', '1.3.0', self.subproject, location=self.current_node)
                arg_strings.append(str(arg))

        def arg_replace(match: T.Match[str]) -> str:
            idx = int(match.group(1))
            if idx >= len(arg_strings):
                raise InvalidArguments(f'Format placeholder @{idx}@ out of range.')
            return arg_strings[idx]

        return re.sub(r'@(\d+)@', arg_replace, self.held_object)

    @noKwargs
    @noPosargs
    @FeatureNew('str.splitlines', '1.2.0')
    @InterpreterObject.method('splitlines')
    def splitlines_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> T.List[str]:
        return self.held_object.splitlines()

    @noKwargs
    @typed_pos_args('str.join', varargs=str)
    @InterpreterObject.method('join')
    def join_method(self, args: T.Tuple[T.List[str]], kwargs: TYPE_kwargs) -> str:
        return self.held_object.join(args[0])

    @noKwargs
    @FeatureNew('str.replace', '0.58.0')
    @typed_pos_args('str.replace', str, str)
    @InterpreterObject.method('replace')
    def replace_method(self, args: T.Tuple[str, str], kwargs: TYPE_kwargs) -> str:
        return self.held_object.replace(args[0], args[1])

    @noKwargs
    @typed_pos_args('str.split', optargs=[str])
    @InterpreterObject.method('split')
    def split_method(self, args: T.Tuple[T.Optional[str]], kwargs: TYPE_kwargs) -> T.List[str]:
        return self.held_object.split(args[0])

    @noKwargs
    @typed_pos_args('str.strip', optargs=[str])
    @InterpreterObject.method('strip')
    def strip_method(self, args: T.Tuple[T.Optional[str]], kwargs: TYPE_kwargs) -> str:
        if args[0]:
            FeatureNew.single_use('str.strip with a positional argument', '0.43.0', self.subproject, location=self.current_node)
        return self.held_object.strip(args[0])

    @noKwargs
    @FeatureNew('str.substring', '0.56.0')
    @typed_pos_args('str.substring', optargs=[int, int])
    @InterpreterObject.method('substring')
    def substring_method(self, args: T.Tuple[T.Optional[int], T.Optional[int]], kwargs: TYPE_kwargs) -> str:
        start = args[0] if args[0] is not None else 0
        end = args[1] if args[1] is not None else len(self.held_object)
        return self.held_object[start:end]

    @noKwargs
    @noPosargs
    @InterpreterObject.method('to_int')
    def to_int_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> int:
        try:
            return int(self.held_object)
        except ValueError:
            raise InvalidArguments(f'String {self.held_object!r} cannot be converted to int')

    @noKwargs
    @noPosargs
    @InterpreterObject.method('to_lower')
    def to_lower_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return self.held_object.lower()

    @noKwargs
    @noPosargs
    @InterpreterObject.method('to_upper')
    def to_upper_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return self.held_object.upper()

    @noKwargs
    @noPosargs
    @InterpreterObject.method('underscorify')
    def underscorify_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return re.sub(r'[^a-zA-Z0-9]', '_', self.held_object)

    @noKwargs
    @InterpreterObject.method('version_compare')
    @typed_pos_args('str.version_compare', varargs=str, min_varargs=1)
    def version_compare_method(self, args: T.Tuple[T.List[str]], kwargs: TYPE_kwargs) -> bool:
        if len(args[0]) > 1:
            FeatureNew.single_use('version_compare() with multiple arguments', '1.8.0', self.subproject, location=self.current_node)
        return version_compare_many(self.held_object, args[0])[0]

    @staticmethod
    def _op_div(this: str, other: str) -> str:
        return os.path.join(this, other).replace('\\', '/')

    @FeatureNew('/ with string arguments', '0.49.0')
    @typed_operator(MesonOperator.DIV, str)
    @InterpreterObject.operator(MesonOperator.DIV)
    def op_div(self, other: str) -> str:
        return self._op_div(self.held_object, other)

    @typed_operator(MesonOperator.INDEX, int)
    @InterpreterObject.operator(MesonOperator.INDEX)
    def op_index(self, other: int) -> str:
        try:
            return self.held_object[other]
        except IndexError:
            raise InvalidArguments(f'Index {other} out of bounds of string of size {len(self.held_object)}.')

    @FeatureNew('"in" string operator', '1.0.0')
    @typed_operator(MesonOperator.IN, str)
    @InterpreterObject.operator(MesonOperator.IN)
    def op_in(self, other: str) -> bool:
        return other in self.held_object

    @FeatureNew('"not in" string operator', '1.0.0')
    @typed_operator(MesonOperator.NOT_IN, str)
    @InterpreterObject.operator(MesonOperator.NOT_IN)
    def op_notin(self, other: str) -> bool:
        return other not in self.held_object


class MesonVersionString(str):
    pass

class MesonVersionStringHolder(StringHolder):
    @noKwargs
    @typed_pos_args('str.version_compare', str)
    @InterpreterObject.method('version_compare')
    def version_compare_method(self, args: T.Tuple[str], kwargs: TYPE_kwargs) -> bool:
        self.interpreter.tmp_meson_version = args[0]
        return version_compare(self.held_object, args[0])

# These special subclasses of string exist to cover the case where a dependency
# exports a string variable interchangeable with a system dependency. This
# matters because a dependency can only have string-type get_variable() return
# values. If at any time dependencies start supporting additional variable
# types, this class could be deprecated.
class DependencyVariableString(str):
    pass

class DependencyVariableStringHolder(StringHolder):
    @InterpreterObject.operator(MesonOperator.DIV)
    def op_div(self, other: str) -> T.Union[str, DependencyVariableString]:
        ret = super().op_div(other)
        if '..' in other:
            return ret
        return DependencyVariableString(ret)


class OptionString(str):
    optname: str

    def __new__(cls, value: str, name: str) -> 'OptionString':
        obj = str.__new__(cls, value)
        obj.optname = name
        return obj

    def __getnewargs__(self) -> T.Tuple[str, str]: # type: ignore # because the entire point of this is to diverge
        return (str(self), self.optname)


class OptionStringHolder(StringHolder):
    held_object: OptionString

    @InterpreterObject.operator(MesonOperator.DIV)
    def op_div(self, other: str) -> T.Union[str, OptionString]:
        ret = super().op_div(other)
        name = self._op_div(self.held_object.optname, other)
        return OptionString(ret, name)
