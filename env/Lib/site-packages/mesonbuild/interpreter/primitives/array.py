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
    FeatureNew,

    TYPE_var,

    InvalidArguments,
)
from ...mparser import PlusAssignmentNode

if T.TYPE_CHECKING:
    # Object holders need the actual interpreter
    from ...interpreter import Interpreter
    from ...interpreterbase import TYPE_kwargs

class ArrayHolder(ObjectHolder[T.List[TYPE_var]], IterableObject):
    def __init__(self, obj: T.List[TYPE_var], interpreter: 'Interpreter') -> None:
        super().__init__(obj, interpreter)
        self.methods.update({
            'contains': self.contains_method,
            'length': self.length_method,
            'get': self.get_method,
        })

        self.trivial_operators.update({
            MesonOperator.EQUALS: (list, lambda x: self.held_object == x),
            MesonOperator.NOT_EQUALS: (list, lambda x: self.held_object != x),
            MesonOperator.IN: (object, lambda x: x in self.held_object),
            MesonOperator.NOT_IN: (object, lambda x: x not in self.held_object),
        })

        # Use actual methods for functions that require additional checks
        self.operators.update({
            MesonOperator.PLUS: self.op_plus,
            MesonOperator.INDEX: self.op_index,
        })

    def display_name(self) -> str:
        return 'array'

    def iter_tuple_size(self) -> None:
        return None

    def iter_self(self) -> T.Iterator[TYPE_var]:
        return iter(self.held_object)

    def size(self) -> int:
        return len(self.held_object)

    @noArgsFlattening
    @noKwargs
    @typed_pos_args('array.contains', object)
    def contains_method(self, args: T.Tuple[object], kwargs: TYPE_kwargs) -> bool:
        def check_contains(el: T.List[TYPE_var]) -> bool:
            for element in el:
                if isinstance(element, list):
                    found = check_contains(element)
                    if found:
                        return True
                if element == args[0]:
                    return True
            return False
        return check_contains(self.held_object)

    @noKwargs
    @noPosargs
    def length_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> int:
        return len(self.held_object)

    @noArgsFlattening
    @noKwargs
    @typed_pos_args('array.get', int, optargs=[object])
    def get_method(self, args: T.Tuple[int, T.Optional[TYPE_var]], kwargs: TYPE_kwargs) -> TYPE_var:
        index = args[0]
        if index < -len(self.held_object) or index >= len(self.held_object):
            if args[1] is None:
                raise InvalidArguments(f'Array index {index} is out of bounds for array of size {len(self.held_object)}.')
            return args[1]
        return self.held_object[index]

    @typed_operator(MesonOperator.PLUS, object)
    def op_plus(self, other: TYPE_var) -> T.List[TYPE_var]:
        if not isinstance(other, list):
            if not isinstance(self.current_node, PlusAssignmentNode):
                FeatureNew.single_use('list.<plus>', '0.60.0', self.subproject, 'The right hand operand was not a list.',
                                      location=self.current_node)
            other = [other]
        return self.held_object + other

    @typed_operator(MesonOperator.INDEX, int)
    def op_index(self, other: int) -> TYPE_var:
        try:
            return self.held_object[other]
        except IndexError:
            raise InvalidArguments(f'Index {other} out of bounds of array of size {len(self.held_object)}.')
