# SPDX-License-Identifier: Apache-2.0
# Copyright 2021 The Meson development team
from __future__ import annotations

import typing as T

from ...interpreterbase import (
    InterpreterObject,
    IterableObject,
    MesonOperator,
    ObjectHolder,
    FeatureNew,
    typed_operator,
    noKwargs,
    noPosargs,
    noArgsFlattening,
    typed_pos_args,

    TYPE_var,

    InvalidArguments,
)

if T.TYPE_CHECKING:
    from ...interpreterbase import TYPE_kwargs

class DictHolder(ObjectHolder[T.Dict[str, TYPE_var]], IterableObject):
    # Operators that only require type checks
    TRIVIAL_OPERATORS = {
        # Arithmetic
        MesonOperator.PLUS: (dict, lambda obj, x: {**obj.held_object, **x}),

        # Comparison
        MesonOperator.EQUALS: (dict, lambda obj, x: obj.held_object == x),
        MesonOperator.NOT_EQUALS: (dict, lambda obj, x: obj.held_object != x),
        MesonOperator.IN: (str, lambda obj, x: x in obj.held_object),
        MesonOperator.NOT_IN: (str, lambda obj, x: x not in obj.held_object),
    }

    def display_name(self) -> str:
        return 'dict'

    def iter_tuple_size(self) -> int:
        return 2

    def iter_self(self) -> T.Iterator[T.Tuple[str, TYPE_var]]:
        return iter(self.held_object.items())

    def size(self) -> int:
        return len(self.held_object)

    def _keys_getter(self) -> T.List[str]:
        return sorted(self.held_object)

    @noKwargs
    @typed_pos_args('dict.has_key', str)
    @InterpreterObject.method('has_key')
    def has_key_method(self, args: T.Tuple[str], kwargs: TYPE_kwargs) -> bool:
        return args[0] in self.held_object

    @noKwargs
    @noPosargs
    @InterpreterObject.method('keys')
    def keys_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> T.List[str]:
        return self._keys_getter()

    @noKwargs
    @noPosargs
    @InterpreterObject.method('values')
    @FeatureNew('dict.values', '1.10.0')
    def values_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> T.List[TYPE_var]:
        return [self.held_object[k] for k in self._keys_getter()]

    @noArgsFlattening
    @noKwargs
    @typed_pos_args('dict.get', str, optargs=[object])
    @InterpreterObject.method('get')
    def get_method(self, args: T.Tuple[str, T.Optional[TYPE_var]], kwargs: TYPE_kwargs) -> TYPE_var:
        if args[0] in self.held_object:
            return self.held_object[args[0]]
        if args[1] is not None:
            return args[1]
        raise InvalidArguments(f'Key {args[0]!r} is not in the dictionary.')

    @typed_operator(MesonOperator.INDEX, str)
    @InterpreterObject.operator(MesonOperator.INDEX)
    def op_index(self, other: str) -> TYPE_var:
        if other not in self.held_object:
            raise InvalidArguments(f'Key {other} is not in the dictionary.')
        return self.held_object[other]
