# SPDX-License-Identifier: Apache-2.0
# Copyright 2021 The Meson development team
from __future__ import annotations

import typing as T

from ...interpreterbase import (
    ObjectHolder,
    IterableObject,
    MesonOperator,
    typed_operator,
    noKwargs,
    noPosargs,
    noArgsFlattening,
    typed_pos_args,

    TYPE_var,

    InvalidArguments,
)

if T.TYPE_CHECKING:
    # Object holders need the actual interpreter
    from ...interpreter import Interpreter
    from ...interpreterbase import TYPE_kwargs

class DictHolder(ObjectHolder[T.Dict[str, TYPE_var]], IterableObject):
    def __init__(self, obj: T.Dict[str, TYPE_var], interpreter: 'Interpreter') -> None:
        super().__init__(obj, interpreter)
        self.methods.update({
            'has_key': self.has_key_method,
            'keys': self.keys_method,
            'get': self.get_method,
        })

        self.trivial_operators.update({
            # Arithmetic
            MesonOperator.PLUS: (dict, lambda x: {**self.held_object, **x}),

            # Comparison
            MesonOperator.EQUALS: (dict, lambda x: self.held_object == x),
            MesonOperator.NOT_EQUALS: (dict, lambda x: self.held_object != x),
            MesonOperator.IN: (str, lambda x: x in self.held_object),
            MesonOperator.NOT_IN: (str, lambda x: x not in self.held_object),
        })

        # Use actual methods for functions that require additional checks
        self.operators.update({
            MesonOperator.INDEX: self.op_index,
        })

    def display_name(self) -> str:
        return 'dict'

    def iter_tuple_size(self) -> int:
        return 2

    def iter_self(self) -> T.Iterator[T.Tuple[str, TYPE_var]]:
        return iter(self.held_object.items())

    def size(self) -> int:
        return len(self.held_object)

    @noKwargs
    @typed_pos_args('dict.has_key', str)
    def has_key_method(self, args: T.Tuple[str], kwargs: TYPE_kwargs) -> bool:
        return args[0] in self.held_object

    @noKwargs
    @noPosargs
    def keys_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> T.List[str]:
        return sorted(self.held_object)

    @noArgsFlattening
    @noKwargs
    @typed_pos_args('dict.get', str, optargs=[object])
    def get_method(self, args: T.Tuple[str, T.Optional[TYPE_var]], kwargs: TYPE_kwargs) -> TYPE_var:
        if args[0] in self.held_object:
            return self.held_object[args[0]]
        if args[1] is not None:
            return args[1]
        raise InvalidArguments(f'Key {args[0]!r} is not in the dictionary.')

    @typed_operator(MesonOperator.INDEX, str)
    def op_index(self, other: str) -> TYPE_var:
        if other not in self.held_object:
            raise InvalidArguments(f'Key {other} is not in the dictionary.')
        return self.held_object[other]
