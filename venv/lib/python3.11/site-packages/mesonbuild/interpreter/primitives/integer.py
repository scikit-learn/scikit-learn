# SPDX-License-Identifier: Apache-2.0
# Copyright 2021 The Meson development team
from __future__ import annotations

from ...interpreterbase import (
    FeatureBroken, InvalidArguments, MesonOperator, ObjectHolder, KwargInfo,
    noKwargs, noPosargs, typed_operator, typed_kwargs
)

import typing as T

if T.TYPE_CHECKING:
    # Object holders need the actual interpreter
    from ...interpreter import Interpreter
    from ...interpreterbase import TYPE_var, TYPE_kwargs

class IntegerHolder(ObjectHolder[int]):
    def __init__(self, obj: int, interpreter: 'Interpreter') -> None:
        super().__init__(obj, interpreter)
        self.methods.update({
            'is_even': self.is_even_method,
            'is_odd': self.is_odd_method,
            'to_string': self.to_string_method,
        })

        self.trivial_operators.update({
            # Arithmetic
            MesonOperator.UMINUS: (None, lambda x: -self.held_object),
            MesonOperator.PLUS: (int, lambda x: self.held_object + x),
            MesonOperator.MINUS: (int, lambda x: self.held_object - x),
            MesonOperator.TIMES: (int, lambda x: self.held_object * x),

            # Comparison
            MesonOperator.EQUALS: (int, lambda x: self.held_object == x),
            MesonOperator.NOT_EQUALS: (int, lambda x: self.held_object != x),
            MesonOperator.GREATER: (int, lambda x: self.held_object > x),
            MesonOperator.LESS: (int, lambda x: self.held_object < x),
            MesonOperator.GREATER_EQUALS: (int, lambda x: self.held_object >= x),
            MesonOperator.LESS_EQUALS: (int, lambda x: self.held_object <= x),
        })

        # Use actual methods for functions that require additional checks
        self.operators.update({
            MesonOperator.DIV: self.op_div,
            MesonOperator.MOD: self.op_mod,
        })

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
    def is_even_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> bool:
        return self.held_object % 2 == 0

    @noKwargs
    @noPosargs
    def is_odd_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> bool:
        return self.held_object % 2 != 0

    @typed_kwargs(
        'to_string',
        KwargInfo('fill', int, default=0, since='1.3.0')
    )
    @noPosargs
    def to_string_method(self, args: T.List[TYPE_var], kwargs: T.Dict[str, T.Any]) -> str:
        return str(self.held_object).zfill(kwargs['fill'])

    @typed_operator(MesonOperator.DIV, int)
    def op_div(self, other: int) -> int:
        if other == 0:
            raise InvalidArguments('Tried to divide by 0')
        return self.held_object // other

    @typed_operator(MesonOperator.MOD, int)
    def op_mod(self, other: int) -> int:
        if other == 0:
            raise InvalidArguments('Tried to divide by 0')
        return self.held_object % other
