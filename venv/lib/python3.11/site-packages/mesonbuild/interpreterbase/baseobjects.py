# SPDX-License-Identifier: Apache-2.0
# Copyright 2013-2021 The Meson development team

from __future__ import annotations

from .. import mparser
from .exceptions import InvalidCode, InvalidArguments
from .helpers import flatten, resolve_second_level_holders
from .operator import MesonOperator
from ..mesonlib import HoldableObject, MesonBugException
import textwrap

import typing as T
from abc import ABCMeta
from contextlib import AbstractContextManager

if T.TYPE_CHECKING:
    from typing_extensions import Protocol, TypeAlias

    # Object holders need the actual interpreter
    from ..interpreter import Interpreter

    __T = T.TypeVar('__T', bound='TYPE_var', contravariant=True)

    class OperatorCall(Protocol[__T]):
        def __call__(self, other: __T) -> 'TYPE_var': ...


TV_func = T.TypeVar('TV_func', bound=T.Callable[..., T.Any])

TYPE_elementary: TypeAlias = T.Union[str, int, bool, T.Sequence['TYPE_elementary'], T.Dict[str, 'TYPE_elementary']]
TYPE_var: TypeAlias = T.Union[TYPE_elementary, HoldableObject, 'MesonInterpreterObject', T.Sequence['TYPE_var'], T.Dict[str, 'TYPE_var']]
TYPE_nvar = T.Union[TYPE_var, mparser.BaseNode]
TYPE_kwargs = T.Dict[str, TYPE_var]
TYPE_nkwargs = T.Dict[str, TYPE_nvar]
TYPE_key_resolver = T.Callable[[mparser.BaseNode], str]

SubProject = T.NewType('SubProject', str)

class InterpreterObject:
    def __init__(self, *, subproject: T.Optional['SubProject'] = None) -> None:
        self.methods: T.Dict[
            str,
            T.Callable[[T.List[TYPE_var], TYPE_kwargs], TYPE_var]
        ] = {}
        self.operators: T.Dict[MesonOperator, 'OperatorCall'] = {}
        self.trivial_operators: T.Dict[
            MesonOperator,
            T.Tuple[
                T.Union[T.Type, T.Tuple[T.Type, ...]],
                'OperatorCall'
            ]
        ] = {}
        # Current node set during a method call. This can be used as location
        # when printing a warning message during a method call.
        self.current_node:  mparser.BaseNode = None
        self.subproject = subproject or SubProject('')

        # Some default operators supported by all objects
        self.operators.update({
            MesonOperator.EQUALS: self.op_equals,
            MesonOperator.NOT_EQUALS: self.op_not_equals,
        })

    # The type of the object that can be printed to the user
    def display_name(self) -> str:
        return type(self).__name__

    def method_call(
                self,
                method_name: str,
                args: T.List[TYPE_var],
                kwargs: TYPE_kwargs
            ) -> TYPE_var:
        if method_name in self.methods:
            method = self.methods[method_name]
            if not getattr(method, 'no-args-flattening', False):
                args = flatten(args)
            if not getattr(method, 'no-second-level-holder-flattening', False):
                args, kwargs = resolve_second_level_holders(args, kwargs)
            return method(args, kwargs)
        raise InvalidCode(f'Unknown method "{method_name}" in object {self} of type {type(self).__name__}.')

    def operator_call(self, operator: MesonOperator, other: TYPE_var) -> TYPE_var:
        if operator in self.trivial_operators:
            op = self.trivial_operators[operator]
            if op[0] is None and other is not None:
                raise MesonBugException(f'The unary operator `{operator.value}` of {self.display_name()} was passed the object {other} of type {type(other).__name__}')
            if op[0] is not None and not isinstance(other, op[0]):
                raise InvalidArguments(f'The `{operator.value}` operator of {self.display_name()} does not accept objects of type {type(other).__name__} ({other})')
            return op[1](other)
        if operator in self.operators:
            return self.operators[operator](other)
        raise InvalidCode(f'Object {self} of type {self.display_name()} does not support the `{operator.value}` operator.')

    # Default comparison operator support
    def _throw_comp_exception(self, other: TYPE_var, opt_type: str) -> T.NoReturn:
        raise InvalidArguments(textwrap.dedent(
            f'''
                Trying to compare values of different types ({self.display_name()}, {type(other).__name__}) using {opt_type}.
                This was deprecated and undefined behavior previously and is as of 0.60.0 a hard error.
            '''
        ))

    def op_equals(self, other: TYPE_var) -> bool:
        # We use `type(...) == type(...)` here to enforce an *exact* match for comparison. We
        # don't want comparisons to be possible where `isinstance(derived_obj, type(base_obj))`
        # would pass because this comparison must never be true: `derived_obj == base_obj`
        if type(self) is not type(other):
            self._throw_comp_exception(other, '==')
        return self == other

    def op_not_equals(self, other: TYPE_var) -> bool:
        if type(self) is not type(other):
            self._throw_comp_exception(other, '!=')
        return self != other

class MesonInterpreterObject(InterpreterObject):
    ''' All non-elementary objects and non-object-holders should be derived from this '''

class MutableInterpreterObject:
    ''' Dummy class to mark the object type as mutable '''

HoldableTypes = (HoldableObject, int, bool, str, list, dict)
TYPE_HoldableTypes = T.Union[TYPE_var, HoldableObject]
InterpreterObjectTypeVar = T.TypeVar('InterpreterObjectTypeVar', bound=TYPE_HoldableTypes)

class ObjectHolder(InterpreterObject, T.Generic[InterpreterObjectTypeVar]):
    def __init__(self, obj: InterpreterObjectTypeVar, interpreter: 'Interpreter') -> None:
        super().__init__(subproject=interpreter.subproject)
        # This causes some type checkers to assume that obj is a base
        # HoldableObject, not the specialized type, so only do this assert in
        # non-type checking situations
        if not T.TYPE_CHECKING:
            assert isinstance(obj, HoldableTypes), f'This is a bug: Trying to hold object of type `{type(obj).__name__}` that is not in `{HoldableTypes}`'
        self.held_object = obj
        self.interpreter = interpreter
        self.env = self.interpreter.environment

    # Hide the object holder abstraction from the user
    def display_name(self) -> str:
        return type(self.held_object).__name__

    # Override default comparison operators for the held object
    def op_equals(self, other: TYPE_var) -> bool:
        # See the comment from InterpreterObject why we are using `type()` here.
        if type(self.held_object) is not type(other):
            self._throw_comp_exception(other, '==')
        return self.held_object == other

    def op_not_equals(self, other: TYPE_var) -> bool:
        if type(self.held_object) is not type(other):
            self._throw_comp_exception(other, '!=')
        return self.held_object != other

    def __repr__(self) -> str:
        return f'<[{type(self).__name__}] holds [{type(self.held_object).__name__}]: {self.held_object!r}>'

class IterableObject(metaclass=ABCMeta):
    '''Base class for all objects that can be iterated over in a foreach loop'''

    def iter_tuple_size(self) -> T.Optional[int]:
        '''Return the size of the tuple for each iteration. Returns None if only a single value is returned.'''
        raise MesonBugException(f'iter_tuple_size not implemented for {self.__class__.__name__}')

    def iter_self(self) -> T.Iterator[T.Union[TYPE_var, T.Tuple[TYPE_var, ...]]]:
        raise MesonBugException(f'iter not implemented for {self.__class__.__name__}')

    def size(self) -> int:
        raise MesonBugException(f'size not implemented for {self.__class__.__name__}')

class ContextManagerObject(MesonInterpreterObject, AbstractContextManager):
    def __init__(self, subproject: 'SubProject') -> None:
        super().__init__(subproject=subproject)
