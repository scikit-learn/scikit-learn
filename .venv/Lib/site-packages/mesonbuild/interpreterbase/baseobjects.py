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
    from typing_extensions import TypeAlias

    # Object holders need the actual interpreter
    from ..interpreter import Interpreter


TV_func = T.TypeVar('TV_func', bound=T.Callable[..., T.Any])

TYPE_elementary: TypeAlias = T.Union[str, int, bool, T.Sequence['TYPE_elementary'], T.Dict[str, 'TYPE_elementary']]
TYPE_var: TypeAlias = T.Union[TYPE_elementary, HoldableObject, 'MesonInterpreterObject', T.Sequence['TYPE_var'], T.Dict[str, 'TYPE_var']]
TYPE_nvar = T.Union[TYPE_var, mparser.BaseNode]
TYPE_kwargs = T.Dict[str, TYPE_var]
TYPE_nkwargs = T.Dict[str, TYPE_nvar]
TYPE_key_resolver = T.Callable[[mparser.BaseNode], str]
TYPE_op_arg = T.TypeVar('TYPE_op_arg', bound='TYPE_var', contravariant=True)
TYPE_op_func = T.Callable[[TYPE_op_arg, TYPE_op_arg], TYPE_var]
TYPE_method_func = T.Callable[['InterpreterObject', T.List[TYPE_var], TYPE_kwargs], TYPE_var]


SubProject = T.NewType('SubProject', str)

class InterpreterObject:
    TRIVIAL_OPERATORS: T.Dict[
        MesonOperator,
        T.Tuple[
            T.Union[T.Type, T.Tuple[T.Type, ...]],
            TYPE_op_func
        ]
    ] = {}

    OPERATORS: T.Dict[MesonOperator, TYPE_op_func] = {}

    METHODS: T.Dict[
        str,
        TYPE_method_func,
    ] = {}

    def __init_subclass__(cls: T.Type[InterpreterObject], **kwargs: T.Any) -> None:
        super().__init_subclass__(**kwargs)
        saved_trivial_operators = cls.TRIVIAL_OPERATORS

        cls.METHODS = {}
        cls.OPERATORS = {}
        cls.TRIVIAL_OPERATORS = {}

        # Compute inherited operators and methods according to the Python resolution
        # order.  Reverse the result of mro() because update() will overwrite entries
        # that are set by the superclass with those that are set by the subclass.
        for superclass in reversed(cls.mro()[1:]):
            if superclass is InterpreterObject:
                # InterpreterObject cannot use @InterpreterObject.operator because
                # __init_subclass__ does not operate on InterpreterObject itself
                cls.OPERATORS.update({
                    MesonOperator.EQUALS: InterpreterObject.op_equals,
                    MesonOperator.NOT_EQUALS: InterpreterObject.op_not_equals
                })

            elif issubclass(superclass, InterpreterObject):
                cls.METHODS.update(superclass.METHODS)
                cls.OPERATORS.update(superclass.OPERATORS)
                cls.TRIVIAL_OPERATORS.update(superclass.TRIVIAL_OPERATORS)

        for name, method in cls.__dict__.items():
            if hasattr(method, 'meson_method'):
                cls.METHODS[method.meson_method] = method
            if hasattr(method, 'meson_operator'):
                cls.OPERATORS[method.meson_operator] = method
        cls.TRIVIAL_OPERATORS.update(saved_trivial_operators)

    @staticmethod
    def method(name: str) -> T.Callable[[TV_func], TV_func]:
        '''Decorator that tags a Python method as the implementation of a method
           for the Meson interpreter'''
        def decorator(f: TV_func) -> TV_func:
            f.meson_method = name    # type: ignore[attr-defined]
            return f
        return decorator

    @staticmethod
    def operator(op: MesonOperator) -> T.Callable[[TV_func], TV_func]:
        '''Decorator that tags a method as the implementation of an operator
           for the Meson interpreter'''
        def decorator(f: TV_func) -> TV_func:
            f.meson_operator = op    # type: ignore[attr-defined]
            return f
        return decorator

    def __init__(self, *, subproject: T.Optional['SubProject'] = None) -> None:
        # Current node set during a method call. This can be used as location
        # when printing a warning message during a method call.
        self.current_node:  mparser.BaseNode = None
        self.subproject = subproject or SubProject('')

    # The type of the object that can be printed to the user
    def display_name(self) -> str:
        return type(self).__name__

    def method_call(
                self,
                method_name: str,
                args: T.List[TYPE_var],
                kwargs: TYPE_kwargs
            ) -> TYPE_var:
        if method_name in self.METHODS:
            method = self.METHODS[method_name]
            if not getattr(method, 'no-args-flattening', False):
                args = flatten(args)
            if not getattr(method, 'no-second-level-holder-flattening', False):
                args, kwargs = resolve_second_level_holders(args, kwargs)
            return method(self, args, kwargs)
        raise InvalidCode(f'Unknown method "{method_name}" in object {self} of type {type(self).__name__}.')

    def operator_call(self, operator: MesonOperator, other: TYPE_var) -> TYPE_var:
        if operator in self.TRIVIAL_OPERATORS:
            op = self.TRIVIAL_OPERATORS[operator]
            if op[0] is None and other is not None:
                raise MesonBugException(f'The unary operator `{operator.value}` of {self.display_name()} was passed the object {other} of type {type(other).__name__}')
            if op[0] is not None and not isinstance(other, op[0]):
                raise InvalidArguments(f'The `{operator.value}` operator of {self.display_name()} does not accept objects of type {type(other).__name__} ({other})')
            return op[1](self, other)
        if operator in self.OPERATORS:
            return self.OPERATORS[operator](self, other)

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

class UnknownValue(MesonInterpreterObject):
    '''This class is only used for the rewriter/static introspection tool and
    indicates that a value cannot be determined statically, either because of
    limitations in our code or because the value differs from machine to
    machine.'''

class UndefinedVariable(MesonInterpreterObject):
    '''This class is only used for the rewriter/static introspection tool and
    represents the `value` a meson-variable has if it was never written to.'''

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
    @InterpreterObject.operator(MesonOperator.EQUALS)
    def op_equals(self, other: TYPE_var) -> bool:
        # See the comment from InterpreterObject why we are using `type()` here.
        if type(self.held_object) is not type(other):
            self._throw_comp_exception(other, '==')
        return self.held_object == other

    @InterpreterObject.operator(MesonOperator.NOT_EQUALS)
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
