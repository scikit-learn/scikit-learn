# SPDX-License-Identifier: Apache-2.0
# Copyright 2021 The Meson development team
from __future__ import annotations

import typing as T

from ...interpreterbase import (
    InterpreterObject,
    IterableObject,
    KwargInfo,
    MesonOperator,
    ObjectHolder,
    typed_operator,
    noKwargs,
    noPosargs,
    noArgsFlattening,
    typed_kwargs,
    typed_pos_args,
    FeatureNew,

    TYPE_var,

    InvalidArguments,
)
from ...mparser import PlusAssignmentNode

if T.TYPE_CHECKING:
    from ...interpreterbase import TYPE_kwargs

class ArrayHolder(ObjectHolder[T.List[TYPE_var]], IterableObject):
    # Operators that only require type checks
    TRIVIAL_OPERATORS = {
        MesonOperator.EQUALS: (list, lambda obj, x: obj.held_object == x),
        MesonOperator.NOT_EQUALS: (list, lambda obj, x: obj.held_object != x),
        MesonOperator.IN: (object, lambda obj, x: x in obj.held_object),
        MesonOperator.NOT_IN: (object, lambda obj, x: x not in obj.held_object),
    }

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
    @InterpreterObject.method('contains')
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
    @InterpreterObject.method('length')
    def length_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> int:
        return len(self.held_object)

    @noArgsFlattening
    @noKwargs
    @typed_pos_args('array.get', int, optargs=[object])
    @InterpreterObject.method('get')
    def get_method(self, args: T.Tuple[int, T.Optional[TYPE_var]], kwargs: TYPE_kwargs) -> TYPE_var:
        index = args[0]
        if index < -len(self.held_object) or index >= len(self.held_object):
            if args[1] is None:
                raise InvalidArguments(f'Array index {index} is out of bounds for array of size {len(self.held_object)}.')
            return args[1]
        return self.held_object[index]

    @FeatureNew('array.slice', '1.10.0')
    @typed_kwargs('array.slice', KwargInfo('step', int, default=1))
    @typed_pos_args('array.slice', optargs=[int, int])
    @InterpreterObject.method('slice')
    def slice_method(self, args: T.Tuple[T.Optional[int], T.Optional[int]], kwargs: T.Dict[str, int]) -> TYPE_var:
        start, stop = args
        if start is not None and stop is None:
            raise InvalidArguments('Providing only one positional slice argument is ambiguous.')
        if kwargs['step'] == 0:
            raise InvalidArguments('Slice step cannot be zero.')
        return self.held_object[start:stop:kwargs['step']]

    @typed_operator(MesonOperator.PLUS, object)
    @InterpreterObject.operator(MesonOperator.PLUS)
    def op_plus(self, other: TYPE_var) -> T.List[TYPE_var]:
        if not isinstance(other, list):
            if not isinstance(self.current_node, PlusAssignmentNode):
                FeatureNew.single_use('list.<plus>', '0.60.0', self.subproject, 'The right hand operand was not a list.',
                                      location=self.current_node)
            other = [other]
        return self.held_object + other

    @typed_operator(MesonOperator.INDEX, int)
    @InterpreterObject.operator(MesonOperator.INDEX)
    def op_index(self, other: int) -> TYPE_var:
        try:
            return self.held_object[other]
        except IndexError:
            raise InvalidArguments(f'Index {other} out of bounds of array of size {len(self.held_object)}.')

    @noPosargs
    @noKwargs
    @FeatureNew('array.flatten', '1.9.0')
    @InterpreterObject.method('flatten')
    def flatten_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> TYPE_var:
        def flatten(obj: TYPE_var) -> T.Iterable[TYPE_var]:
            if isinstance(obj, list):
                for o in obj:
                    yield from flatten(o)
            else:
                yield obj

        return list(flatten(self.held_object))
