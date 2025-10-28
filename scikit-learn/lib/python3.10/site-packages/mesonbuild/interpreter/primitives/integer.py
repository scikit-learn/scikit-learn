# SPDX-License-Identifier: Apache-2.0
# Copyright 2021 The Meson development team
from __future__ import annotations

from ...interpreterbase import (
    InterpreterObject, MesonOperator, ObjectHolder,
    FeatureBroken, InvalidArguments, KwargInfo,
    noKwargs, noPosargs, typed_operator, typed_kwargs
)

import typing as T

if T.TYPE_CHECKING:
    from ...interpreterbase import TYPE_var, TYPE_kwargs

class IntegerHolder(ObjectHolder[int]):
    # Operators that only require type checks
    TRIVIAL_OPERATORS = {
        # Arithmetic
        MesonOperator.UMINUS: (None, lambda obj, x: -obj.held_object),
        MesonOperator.PLUS: (int, lambda obj, x: obj.held_object + x),
        MesonOperator.MINUS: (int, lambda obj, x: obj.held_object - x),
        MesonOperator.TIMES: (int, lambda obj, x: obj.held_object * x),

        # Comparison
        MesonOperator.EQUALS: (int, lambda obj, x: obj.held_object == x),
        MesonOperator.NOT_EQUALS: (int, lambda obj, x: obj.held_object != x),
        MesonOperator.GREATER: (int, lambda obj, x: obj.held_object > x),
        MesonOperator.LESS: (int, lambda obj, x: obj.held_object < x),
        MesonOperator.GREATER_EQUALS: (int, lambda obj, x: obj.held_object >= x),
        MesonOperator.LESS_EQUALS: (int, lambda obj, x: obj.held_object <= x),
    }

    def display_name(self) -> str:
        return 'int'

    def operator_call(self, operator: MesonOperator, other: TYPE_var) -> TYPE_var:
        if isinstance(other, bool):
            FeatureBroken.single_use('int operations with non-int', '1.2.0', self.subproject,
                                     'It is not commutative and only worked because of leaky Python abstractions.',
                                     location=self.current_node)
        return super().operator_call(operator, other)

    @noKwargs
    @noPosargs
    @InterpreterObject.method('is_even')
    def is_even_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> bool:
        return self.held_object % 2 == 0

    @noKwargs
    @noPosargs
    @InterpreterObject.method('is_odd')
    def is_odd_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> bool:
        return self.held_object % 2 != 0

    @typed_kwargs(
        'to_string',
        KwargInfo('fill', int, default=0, since='1.3.0')
    )
    @noPosargs
    @InterpreterObject.method('to_string')
    def to_string_method(self, args: T.List[TYPE_var], kwargs: T.Dict[str, T.Any]) -> str:
        return str(self.held_object).zfill(kwargs['fill'])

    @typed_operator(MesonOperator.DIV, int)
    @InterpreterObject.operator(MesonOperator.DIV)
    def op_div(self, other: int) -> int:
        if other == 0:
            raise InvalidArguments('Tried to divide by 0')
        return self.held_object // other

    @typed_operator(MesonOperator.MOD, int)
    @InterpreterObject.operator(MesonOperator.MOD)
    def op_mod(self, other: int) -> int:
        if other == 0:
            raise InvalidArguments('Tried to divide by 0')
        return self.held_object % other
