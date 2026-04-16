# Copyright 2021 The Meson development team
# SPDX-License-Identifier: Apache-2.0
from __future__ import annotations

from ...interpreterbase import (
    InterpreterObject,
    MesonOperator,
    ObjectHolder,
    typed_pos_args,
    noKwargs,
    noPosargs,

    InvalidArguments
)

import typing as T

if T.TYPE_CHECKING:
    from ...interpreterbase import TYPE_var, TYPE_kwargs

class BooleanHolder(ObjectHolder[bool]):
    TRIVIAL_OPERATORS = {
        MesonOperator.BOOL: (None, lambda obj, x: obj.held_object),
        MesonOperator.NOT: (None, lambda obj, x: not obj.held_object),
        MesonOperator.EQUALS: (bool, lambda obj, x: obj.held_object == x),
        MesonOperator.NOT_EQUALS: (bool, lambda obj, x: obj.held_object != x),
    }

    def display_name(self) -> str:
        return 'bool'

    @noKwargs
    @noPosargs
    @InterpreterObject.method('to_int')
    def to_int_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> int:
        return 1 if self.held_object else 0

    @noKwargs
    @typed_pos_args('bool.to_string', optargs=[str, str])
    @InterpreterObject.method('to_string')
    def to_string_method(self, args: T.Tuple[T.Optional[str], T.Optional[str]], kwargs: TYPE_kwargs) -> str:
        true_str = args[0] or 'true'
        false_str = args[1] or 'false'
        if any(x is not None for x in args) and not all(x is not None for x in args):
            raise InvalidArguments('bool.to_string() must have either no arguments or exactly two string arguments that signify what values to return for true and false.')
        return true_str if self.held_object else false_str
