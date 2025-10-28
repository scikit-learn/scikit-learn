"""Helpers for inspecting Python modules."""

from __future__ import annotations

import ast
import builtins
import contextlib
import enum
import inspect
import re
import sys
import types
import typing
from collections.abc import Mapping
from functools import cached_property, partial, partialmethod, singledispatchmethod
from importlib import import_module
from inspect import Parameter, Signature
from io import StringIO
from types import ClassMethodDescriptorType, MethodDescriptorType, WrapperDescriptorType
from typing import TYPE_CHECKING, Any, ForwardRef

from sphinx.pycode.ast import unparse as ast_unparse
from sphinx.util import logging
from sphinx.util.typing import stringify_annotation

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence
    from inspect import _ParameterKind
    from types import MethodType, ModuleType
    from typing import Final, Protocol, TypeAlias

    from typing_extensions import TypeIs

    class _SupportsGet(Protocol):
        def __get__(self, __instance: Any, __owner: type | None = ...) -> Any: ...  # NoQA: E704

    class _SupportsSet(Protocol):
        # instance and value are contravariants but we do not need that precision
        def __set__(self, __instance: Any, __value: Any) -> None: ...  # NoQA: E704

    class _SupportsDelete(Protocol):
        # instance is contravariant but we do not need that precision
        def __delete__(self, __instance: Any) -> None: ...  # NoQA: E704

    _RoutineType: TypeAlias = (
        types.FunctionType
        | types.LambdaType
        | types.MethodType
        | types.BuiltinFunctionType
        | types.BuiltinMethodType
        | types.WrapperDescriptorType
        | types.MethodDescriptorType
        | types.ClassMethodDescriptorType
    )
    _SignatureType: TypeAlias = Callable[..., Any] | staticmethod | classmethod

logger = logging.getLogger(__name__)

memory_address_re = re.compile(r' at 0x[0-9a-f]{8,16}(?=>)', re.IGNORECASE)

# re-export as is
isasyncgenfunction = inspect.isasyncgenfunction
ismethod = inspect.ismethod
ismethoddescriptor = inspect.ismethoddescriptor
isclass = inspect.isclass
ismodule = inspect.ismodule


def unwrap(obj: Any) -> Any:
    """Get an original object from wrapped object (wrapped functions).

    Mocked objects are returned as is.
    """
    if hasattr(obj, '__sphinx_mock__'):
        # Skip unwrapping mock object to avoid RecursionError
        return obj

    try:
        return inspect.unwrap(obj)
    except ValueError:
        # might be a mock object
        return obj


def unwrap_all(obj: Any, *, stop: Callable[[Any], bool] | None = None) -> Any:
    """Get an original object from wrapped object.

    Unlike :func:`unwrap`, this unwraps partial functions, wrapped functions,
    class methods and static methods.

    When specified, *stop* is a predicate indicating whether an object should
    be unwrapped or not.
    """
    if callable(stop):
        while not stop(obj):
            if ispartial(obj):
                obj = obj.func
            elif inspect.isroutine(obj) and hasattr(obj, '__wrapped__'):
                obj = obj.__wrapped__
            elif isclassmethod(obj) or isstaticmethod(obj):
                obj = obj.__func__
            else:
                return obj
        return obj  # in case the while loop never starts

    while True:
        if ispartial(obj):
            obj = obj.func
        elif inspect.isroutine(obj) and hasattr(obj, '__wrapped__'):
            obj = obj.__wrapped__
        elif isclassmethod(obj) or isstaticmethod(obj):
            obj = obj.__func__
        else:
            return obj


def getall(obj: Any) -> Sequence[str] | None:
    """Get the ``__all__`` attribute of an object as a sequence.

    This returns ``None`` if the given ``obj.__all__`` does not exist and
    raises :exc:`ValueError` if ``obj.__all__`` is not a list or tuple of
    strings.
    """
    __all__ = safe_getattr(obj, '__all__', None)
    if __all__ is None:
        return None
    if isinstance(__all__, list | tuple) and all(isinstance(e, str) for e in __all__):
        return __all__
    raise ValueError(__all__)


def getannotations(obj: Any) -> Mapping[str, Any]:
    """Safely get the ``__annotations__`` attribute of an object."""
    __annotations__ = safe_getattr(obj, '__annotations__', None)
    if isinstance(__annotations__, Mapping):
        return __annotations__
    return {}


def getglobals(obj: Any) -> Mapping[str, Any]:
    """Safely get :attr:`obj.__globals__ <function.__globals__>`."""
    __globals__ = safe_getattr(obj, '__globals__', None)
    if isinstance(__globals__, Mapping):
        return __globals__
    return {}


def getmro(obj: Any) -> tuple[type, ...]:
    """Safely get :attr:`obj.__mro__ <class.__mro__>`."""
    __mro__ = safe_getattr(obj, '__mro__', None)
    if isinstance(__mro__, tuple):
        return __mro__
    return ()


def getorigbases(obj: Any) -> tuple[Any, ...] | None:
    """Safely get ``obj.__orig_bases__``.

    This returns ``None`` if the object is not a class or if ``__orig_bases__``
    is not well-defined (e.g., a non-tuple object or an empty sequence).
    """
    if not isclass(obj):
        return None

    # Get __orig_bases__ from obj.__dict__ to avoid accessing the parent's __orig_bases__.
    # refs: https://github.com/sphinx-doc/sphinx/issues/9607
    __dict__ = safe_getattr(obj, '__dict__', {})
    __orig_bases__ = __dict__.get('__orig_bases__')
    if isinstance(__orig_bases__, tuple) and len(__orig_bases__) > 0:
        return __orig_bases__
    return None


def getslots(obj: Any) -> dict[str, Any] | dict[str, None] | None:
    """Safely get :term:`obj.__slots__ <__slots__>` as a dictionary if any.

    - This returns ``None`` if ``obj.__slots__`` does not exist.
    - This raises a :exc:`TypeError` if *obj* is not a class.
    - This raises a :exc:`ValueError` if ``obj.__slots__`` is invalid.
    """
    if not isclass(obj):
        raise TypeError

    __slots__ = safe_getattr(obj, '__slots__', None)
    if __slots__ is None:
        return None
    elif isinstance(__slots__, dict):
        return __slots__
    elif isinstance(__slots__, str):
        return {__slots__: None}
    elif isinstance(__slots__, list | tuple):
        return dict.fromkeys(__slots__)
    else:
        raise ValueError


def isenumclass(x: Any) -> TypeIs[type[enum.Enum]]:
    """Check if the object is an :class:`enumeration class <enum.Enum>`."""
    return isclass(x) and issubclass(x, enum.Enum)


def isenumattribute(x: Any) -> TypeIs[enum.Enum]:
    """Check if the object is an enumeration attribute."""
    return isinstance(x, enum.Enum)


def unpartial(obj: Any) -> Any:
    """Get an original object from a partial-like object.

    If *obj* is not a partial object, it is returned as is.

    .. seealso:: :func:`ispartial`
    """
    while ispartial(obj):
        obj = obj.func
    return obj


def ispartial(obj: Any) -> TypeIs[partial | partialmethod]:
    """Check if the object is a partial function or method."""
    return isinstance(obj, partial | partialmethod)


def isclassmethod(
    obj: Any,
    cls: Any = None,
    name: str | None = None,
) -> TypeIs[classmethod]:
    """Check if the object is a :class:`classmethod`."""
    if isinstance(obj, classmethod):
        return True
    if ismethod(obj) and obj.__self__ is not None and isclass(obj.__self__):
        return True
    if cls and name:
        # trace __mro__ if the method is defined in parent class
        sentinel = object()
        for basecls in getmro(cls):
            meth = basecls.__dict__.get(name, sentinel)
            if meth is not sentinel:
                return isclassmethod(meth)
    return False


def isstaticmethod(
    obj: Any,
    cls: Any = None,
    name: str | None = None,
) -> TypeIs[staticmethod]:
    """Check if the object is a :class:`staticmethod`."""
    if isinstance(obj, staticmethod):
        return True
    if cls and name:
        # trace __mro__ if the method is defined in parent class
        sentinel = object()
        for basecls in getattr(cls, '__mro__', [cls]):
            meth = basecls.__dict__.get(name, sentinel)
            if meth is not sentinel:
                return isinstance(meth, staticmethod)
    return False


def isdescriptor(x: Any) -> TypeIs[_SupportsGet | _SupportsSet | _SupportsDelete]:
    """Check if the object is a :external+python:term:`descriptor`."""
    return any(
        callable(safe_getattr(x, item, None))
        for item in ('__get__', '__set__', '__delete__')
    )


def isabstractmethod(obj: Any) -> bool:
    """Check if the object is an :func:`abstractmethod`."""
    return safe_getattr(obj, '__isabstractmethod__', False) is True


def isboundmethod(method: MethodType) -> bool:
    """Check if the method is a bound method."""
    return safe_getattr(method, '__self__', None) is not None


def is_cython_function_or_method(obj: Any) -> bool:
    """Check if the object is a function or method in cython."""
    try:
        return obj.__class__.__name__ == 'cython_function_or_method'
    except AttributeError:
        return False


_DESCRIPTOR_LIKE: Final[tuple[type, ...]] = (
    ClassMethodDescriptorType,
    MethodDescriptorType,
    WrapperDescriptorType,
)


def isattributedescriptor(obj: Any) -> bool:
    """Check if the object is an attribute-like descriptor."""
    if inspect.isdatadescriptor(obj):
        # data descriptor is kind of attribute
        return True
    if isdescriptor(obj):
        # non data descriptor
        unwrapped = unwrap(obj)
        if isfunction(unwrapped) or isbuiltin(unwrapped) or ismethod(unwrapped):
            # attribute must not be either function, builtin and method
            return False
        if is_cython_function_or_method(unwrapped):
            # attribute must not be either function and method (for cython)
            return False
        if isclass(unwrapped):
            # attribute must not be a class
            return False
        if isinstance(unwrapped, _DESCRIPTOR_LIKE):
            # attribute must not be a method descriptor
            return False
        # attribute must not be an instancemethod (C-API)
        return type(unwrapped).__name__ != 'instancemethod'
    return False


def is_singledispatch_function(obj: Any) -> bool:
    """Check if the object is a :func:`~functools.singledispatch` function."""
    return (
        inspect.isfunction(obj)
        and hasattr(obj, 'dispatch')
        and hasattr(obj, 'register')
        and obj.dispatch.__module__ == 'functools'
    )


def is_singledispatch_method(obj: Any) -> TypeIs[singledispatchmethod]:
    """Check if the object is a :class:`~functools.singledispatchmethod`."""
    return isinstance(obj, singledispatchmethod)


def isfunction(obj: Any) -> TypeIs[types.FunctionType]:
    """Check if the object is a user-defined function.

    Partial objects are unwrapped before checking them.

    .. seealso:: :external+python:func:`inspect.isfunction`
    """
    return inspect.isfunction(unpartial(obj))


def isbuiltin(obj: Any) -> TypeIs[types.BuiltinFunctionType]:
    """Check if the object is a built-in function or method.

    Partial objects are unwrapped before checking them.

    .. seealso:: :external+python:func:`inspect.isbuiltin`
    """
    return inspect.isbuiltin(unpartial(obj))


def isroutine(obj: Any) -> TypeIs[_RoutineType]:
    """Check if the object is a kind of function or method.

    Partial objects are unwrapped before checking them.

    .. seealso:: :external+python:func:`inspect.isroutine`
    """
    return inspect.isroutine(unpartial(obj))


def iscoroutinefunction(obj: Any) -> TypeIs[Callable[..., types.CoroutineType]]:
    """Check if the object is a :external+python:term:`coroutine` function."""
    obj = unwrap_all(obj, stop=_is_wrapped_coroutine)
    return inspect.iscoroutinefunction(obj)


def _is_wrapped_coroutine(obj: Any) -> bool:
    """Check if the object is wrapped coroutine-function."""
    if isstaticmethod(obj) or isclassmethod(obj) or ispartial(obj):
        # staticmethod, classmethod and partial method are not a wrapped coroutine-function
        # Note: Since 3.10, staticmethod and classmethod becomes a kind of wrappers
        return False
    return hasattr(obj, '__wrapped__')


def isproperty(obj: Any) -> TypeIs[property | cached_property]:
    """Check if the object is property (possibly cached)."""
    return isinstance(obj, property | cached_property)


def isgenericalias(obj: Any) -> TypeIs[types.GenericAlias]:
    """Check if the object is a generic alias."""
    return isinstance(obj, types.GenericAlias | typing._BaseGenericAlias)  # type: ignore[attr-defined]


def safe_getattr(obj: Any, name: str, *defargs: Any) -> Any:
    """A getattr() that turns all exceptions into AttributeErrors."""
    if len(defargs) > 1:
        msg = f'safe_getattr expected at most 3 arguments, got {len(defargs)}'
        raise TypeError(msg)

    try:
        return getattr(obj, name, *defargs)
    except Exception as exc:
        # sometimes accessing a property raises an exception (e.g.
        # NotImplementedError), so let's try to read the attribute directly
        try:
            # In case the object does weird things with attribute access
            # such that accessing `obj.__dict__` may raise an exception
            return obj.__dict__[name]
        except Exception:
            pass

        # this is a catch-all for all the weird things that some modules do
        # with attribute access
        if defargs:
            return defargs[0]

        raise AttributeError(name) from exc


def object_description(obj: Any, *, _seen: frozenset[int] = frozenset()) -> str:
    """A repr() implementation that returns text safe to use in reST context.

    Maintains a set of 'seen' object IDs to detect and avoid infinite recursion.
    """
    seen = _seen
    if isinstance(obj, dict):
        if id(obj) in seen:
            return 'dict(...)'
        seen |= {id(obj)}
        try:
            sorted_keys = sorted(obj)
        except TypeError:
            # Cannot sort dict keys, fall back to using descriptions as a sort key
            sorted_keys = sorted(obj, key=lambda k: object_description(k, _seen=seen))

        items = (
            (
                object_description(key, _seen=seen),
                object_description(obj[key], _seen=seen),
            )
            for key in sorted_keys
        )
        return '{%s}' % ', '.join(f'{key}: {value}' for (key, value) in items)
    elif isinstance(obj, set):
        if id(obj) in seen:
            return 'set(...)'
        seen |= {id(obj)}
        try:
            sorted_values = sorted(obj)
        except TypeError:
            # Cannot sort set values, fall back to using descriptions as a sort key
            sorted_values = sorted(obj, key=lambda x: object_description(x, _seen=seen))
        return '{%s}' % ', '.join(
            object_description(x, _seen=seen) for x in sorted_values
        )
    elif isinstance(obj, frozenset):
        if id(obj) in seen:
            return 'frozenset(...)'
        seen |= {id(obj)}
        try:
            sorted_values = sorted(obj)
        except TypeError:
            # Cannot sort frozenset values, fall back to using descriptions as a sort key
            sorted_values = sorted(obj, key=lambda x: object_description(x, _seen=seen))
        return 'frozenset({%s})' % ', '.join(
            object_description(x, _seen=seen) for x in sorted_values
        )
    elif isinstance(obj, enum.Enum):
        if obj.__repr__.__func__ is not enum.Enum.__repr__:  # type: ignore[attr-defined]
            return repr(obj)
        return f'{obj.__class__.__name__}.{obj.name}'
    elif isinstance(obj, tuple):
        if id(obj) in seen:
            return 'tuple(...)'
        seen |= frozenset([id(obj)])
        return '({}{})'.format(
            ', '.join(object_description(x, _seen=seen) for x in obj),
            ',' * (len(obj) == 1),
        )
    elif isinstance(obj, list):
        if id(obj) in seen:
            return 'list(...)'
        seen |= {id(obj)}
        return '[%s]' % ', '.join(object_description(x, _seen=seen) for x in obj)

    try:
        s = repr(obj)
    except Exception as exc:
        raise ValueError from exc
    # Strip non-deterministic memory addresses such as
    # ``<__main__.A at 0x7f68cb685710>``
    s = memory_address_re.sub('', s)
    return s.replace('\n', ' ')


def is_builtin_class_method(obj: Any, attr_name: str) -> bool:
    """Check whether *attr_name* is implemented on a builtin class.

        >>> is_builtin_class_method(int, '__init__')
        True


    This function is needed since CPython implements ``int.__init__`` via
    descriptors, but PyPy implementation is written in pure Python code.
    """
    mro = getmro(obj)

    try:
        cls = next(c for c in mro if attr_name in safe_getattr(c, '__dict__', {}))
    except StopIteration:
        return False

    try:
        name = safe_getattr(cls, '__name__')
    except AttributeError:
        return False

    return getattr(builtins, name, None) is cls


class DefaultValue:
    """A simple wrapper for default value of the parameters of overload functions."""

    def __init__(self, value: str) -> None:
        self.value = value

    def __eq__(self, other: object) -> bool:
        return self.value == other

    def __repr__(self) -> str:
        return self.value


class TypeAliasForwardRef:
    """Pseudo typing class for :confval:`autodoc_type_aliases`.

    This avoids the error on evaluating the type inside :func:`typing.get_type_hints()`.
    """

    def __init__(self, name: str) -> None:
        self.name = name

    def __call__(self) -> None:
        # Dummy method to imitate special typing classes
        pass

    def __eq__(self, other: Any) -> bool:
        return self.name == other

    def __hash__(self) -> int:
        return hash(self.name)

    def __repr__(self) -> str:
        return self.name


class TypeAliasModule:
    """Pseudo module class for :confval:`autodoc_type_aliases`."""

    def __init__(self, modname: str, mapping: Mapping[str, str]) -> None:
        self.__modname = modname
        self.__mapping = mapping

        self.__module: ModuleType | None = None

    def __getattr__(self, name: str) -> Any:
        fullname = '.'.join(filter(None, [self.__modname, name]))
        if fullname in self.__mapping:
            # exactly matched
            return TypeAliasForwardRef(self.__mapping[fullname])
        else:
            prefix = fullname + '.'
            nested = {k: v for k, v in self.__mapping.items() if k.startswith(prefix)}
            if nested:
                # sub modules or classes found
                return TypeAliasModule(fullname, nested)
            else:
                # no sub modules or classes found.
                try:
                    # return the real submodule if exists
                    return import_module(fullname)
                except ImportError:
                    # return the real class
                    if self.__module is None:
                        self.__module = import_module(self.__modname)

                    return getattr(self.__module, name)


class TypeAliasNamespace(dict[str, Any]):
    """Pseudo namespace class for :confval:`autodoc_type_aliases`.

    Useful for looking up nested objects via ``namespace.foo.bar.Class``.
    """

    def __init__(self, mapping: Mapping[str, str]) -> None:
        super().__init__()
        self.__mapping = mapping

    def __getitem__(self, key: str) -> Any:
        if key in self.__mapping:
            # exactly matched
            return TypeAliasForwardRef(self.__mapping[key])
        else:
            prefix = key + '.'
            nested = {k: v for k, v in self.__mapping.items() if k.startswith(prefix)}
            if nested:
                # sub modules or classes found
                return TypeAliasModule(key, nested)
            else:
                raise KeyError


def _should_unwrap(subject: _SignatureType) -> bool:
    """Check the function should be unwrapped on getting signature."""
    __globals__ = getglobals(subject)
    # contextmanger should be unwrapped
    return (
        __globals__.get('__name__') == 'contextlib'
        and __globals__.get('__file__') == contextlib.__file__
    )


def signature(
    subject: _SignatureType,
    bound_method: bool = False,
    type_aliases: Mapping[str, str] | None = None,
) -> Signature:
    """Return a Signature object for the given *subject*.

    :param bound_method: Specify *subject* is a bound method or not
    """
    if type_aliases is None:
        type_aliases = {}

    try:
        if _should_unwrap(subject):
            signature = inspect.signature(subject)  # type: ignore[arg-type]
        else:
            signature = inspect.signature(subject, follow_wrapped=True)  # type: ignore[arg-type]
    except ValueError:
        # follow built-in wrappers up (ex. functools.lru_cache)
        signature = inspect.signature(subject)  # type: ignore[arg-type]
    parameters = list(signature.parameters.values())
    return_annotation = signature.return_annotation

    try:
        # Resolve annotations using ``get_type_hints()`` and type_aliases.
        localns = TypeAliasNamespace(type_aliases)
        annotations = typing.get_type_hints(subject, None, localns, include_extras=True)
        for i, param in enumerate(parameters):
            if param.name in annotations:
                annotation = annotations[param.name]
                if isinstance(annotation, TypeAliasForwardRef):
                    annotation = annotation.name
                parameters[i] = param.replace(annotation=annotation)
        if 'return' in annotations:
            if isinstance(annotations['return'], TypeAliasForwardRef):
                return_annotation = annotations['return'].name
            else:
                return_annotation = annotations['return']
    except Exception:
        # ``get_type_hints()`` does not support some kind of objects like partial,
        # ForwardRef and so on.
        pass

    if bound_method:
        if inspect.ismethod(subject):
            # ``inspect.signature()`` considers the subject is a bound method and removes
            # first argument from signature.  Therefore no skips are needed here.
            pass
        else:
            if len(parameters) > 0:
                parameters.pop(0)

    # To allow to create signature object correctly for pure python functions,
    # pass an internal parameter __validate_parameters__=False to Signature
    #
    # For example, this helps a function having a default value `inspect._empty`.
    # refs: https://github.com/sphinx-doc/sphinx/issues/7935
    return Signature(
        parameters, return_annotation=return_annotation, __validate_parameters__=False
    )


def evaluate_signature(
    sig: Signature,
    globalns: dict[str, Any] | None = None,
    localns: dict[str, Any] | None = None,
) -> Signature:
    """Evaluate unresolved type annotations in a signature object."""
    if globalns is None:
        globalns = {}
    if localns is None:
        localns = globalns

    parameters = list(sig.parameters.values())
    for i, param in enumerate(parameters):
        if param.annotation:
            annotation = _evaluate(param.annotation, globalns, localns)
            parameters[i] = param.replace(annotation=annotation)

    return_annotation = sig.return_annotation
    if return_annotation:
        return_annotation = _evaluate(return_annotation, globalns, localns)

    return sig.replace(parameters=parameters, return_annotation=return_annotation)


def _evaluate_forwardref(
    ref: ForwardRef,
    globalns: dict[str, Any] | None,
    localns: dict[str, Any] | None,
) -> Any:
    """Evaluate a forward reference."""
    if sys.version_info >= (3, 12, 4):
        # ``type_params`` were added in 3.13 and the signature of _evaluate()
        # is not backward-compatible (it was backported to 3.12.4, so anything
        # before 3.12.4 still has the old signature).
        #
        # See: https://github.com/python/cpython/pull/118104.
        return ref._evaluate(globalns, localns, {}, recursive_guard=frozenset())  # type: ignore[arg-type, misc]
    return ref._evaluate(globalns, localns, frozenset())


def _evaluate(
    annotation: Any,
    globalns: dict[str, Any],
    localns: dict[str, Any],
) -> Any:
    """Evaluate unresolved type annotation."""
    try:
        if isinstance(annotation, str):
            ref = ForwardRef(annotation, True)
            annotation = _evaluate_forwardref(ref, globalns, localns)

            if isinstance(annotation, ForwardRef):
                annotation = _evaluate_forwardref(ref, globalns, localns)
            elif isinstance(annotation, str):
                # might be a ForwardRef'ed annotation in overloaded functions
                ref = ForwardRef(annotation, True)
                annotation = _evaluate_forwardref(ref, globalns, localns)
    except (NameError, TypeError):
        # failed to evaluate type. skipped.
        pass

    return annotation


def stringify_signature(
    sig: Signature,
    show_annotation: bool = True,
    show_return_annotation: bool = True,
    unqualified_typehints: bool = False,
) -> str:
    """Stringify a :class:`~inspect.Signature` object.

    :param show_annotation: If enabled, show annotations on the signature
    :param show_return_annotation: If enabled, show annotation of the return value
    :param unqualified_typehints: If enabled, show annotations as unqualified
                                  (ex. io.StringIO -> StringIO)
    """
    if unqualified_typehints:
        mode = 'smart'
    else:
        mode = 'fully-qualified'

    EMPTY = Parameter.empty

    args = []
    last_kind = None
    for param in sig.parameters.values():
        if (
            param.kind != Parameter.POSITIONAL_ONLY
            and last_kind == Parameter.POSITIONAL_ONLY
        ):
            # PEP-570: Separator for Positional Only Parameter: /
            args.append('/')
        if param.kind == Parameter.KEYWORD_ONLY and last_kind in (
            Parameter.POSITIONAL_OR_KEYWORD,
            Parameter.POSITIONAL_ONLY,
            None,
        ):
            # PEP-3102: Separator for Keyword Only Parameter: *
            args.append('*')

        arg = StringIO()
        if param.kind is Parameter.VAR_POSITIONAL:
            arg.write('*' + param.name)
        elif param.kind is Parameter.VAR_KEYWORD:
            arg.write('**' + param.name)
        else:
            arg.write(param.name)

        if show_annotation and param.annotation is not EMPTY:
            arg.write(': ')
            arg.write(stringify_annotation(param.annotation, mode))  # type: ignore[arg-type]
        if param.default is not EMPTY:
            if show_annotation and param.annotation is not EMPTY:
                arg.write(' = ')
            else:
                arg.write('=')
            arg.write(object_description(param.default))

        args.append(arg.getvalue())
        last_kind = param.kind

    if last_kind is Parameter.POSITIONAL_ONLY:
        # PEP-570: Separator for Positional Only Parameter: /
        args.append('/')

    concatenated_args = ', '.join(args)
    if (
        sig.return_annotation is EMPTY
        or not show_annotation
        or not show_return_annotation
    ):
        return f'({concatenated_args})'
    else:
        retann = stringify_annotation(sig.return_annotation, mode)  # type: ignore[arg-type]
        return f'({concatenated_args}) -> {retann}'


def signature_from_str(signature: str) -> Signature:
    """Create a :class:`~inspect.Signature` object from a string."""
    code = 'def func' + signature + ': pass'
    module = ast.parse(code)
    function = typing.cast(ast.FunctionDef, module.body[0])

    return signature_from_ast(function, code)


def signature_from_ast(node: ast.FunctionDef, code: str = '') -> Signature:
    """Create a :class:`~inspect.Signature` object from an AST node."""
    EMPTY = Parameter.empty

    args: ast.arguments = node.args
    defaults: tuple[ast.expr | None, ...] = tuple(args.defaults)
    pos_only_offset = len(args.posonlyargs)
    defaults_offset = pos_only_offset + len(args.args) - len(defaults)
    # The sequence ``D = args.defaults`` contains non-None AST expressions,
    # so we can use ``None`` as a sentinel value for that to indicate that
    # there is no default value for a specific parameter.
    #
    # Let *p* be the number of positional-only and positional-or-keyword
    # arguments. Note that ``0 <= len(D) <= p`` and ``D[0]`` is the default
    # value corresponding to a positional-only *or* a positional-or-keyword
    # argument. Since a non-default argument cannot follow a default argument,
    # the sequence *D* can be completed on the left by adding None sentinels
    # so that ``len(D) == p`` and ``D[i]`` is the *i*-th default argument.
    defaults = (None,) * defaults_offset + defaults

    # construct the parameter list
    params: list[Parameter] = []

    # positional-only arguments (introduced in Python 3.8)
    for arg, defexpr in zip(args.posonlyargs, defaults, strict=False):
        params.append(_define(Parameter.POSITIONAL_ONLY, arg, code, defexpr=defexpr))

    # normal arguments
    for arg, defexpr in zip(args.args, defaults[pos_only_offset:], strict=False):
        params.append(
            _define(Parameter.POSITIONAL_OR_KEYWORD, arg, code, defexpr=defexpr)
        )

    # variadic positional argument (no possible default expression)
    if args.vararg:
        params.append(
            _define(Parameter.VAR_POSITIONAL, args.vararg, code, defexpr=None)
        )

    # keyword-only arguments
    for arg, defexpr in zip(args.kwonlyargs, args.kw_defaults, strict=False):
        params.append(_define(Parameter.KEYWORD_ONLY, arg, code, defexpr=defexpr))

    # variadic keyword argument (no possible default expression)
    if args.kwarg:
        params.append(_define(Parameter.VAR_KEYWORD, args.kwarg, code, defexpr=None))

    return_annotation = ast_unparse(node.returns, code) or EMPTY
    return Signature(params, return_annotation=return_annotation)


def _define(
    kind: _ParameterKind,
    arg: ast.arg,
    code: str,
    *,
    defexpr: ast.expr | None,
) -> Parameter:
    EMPTY = Parameter.empty

    default = EMPTY if defexpr is None else DefaultValue(ast_unparse(defexpr, code))
    annotation = ast_unparse(arg.annotation, code) or EMPTY
    return Parameter(arg.arg, kind, default=default, annotation=annotation)


def getdoc(
    obj: Any,
    attrgetter: Callable = safe_getattr,
    allow_inherited: bool = False,
    cls: Any = None,
    name: str | None = None,
) -> str | None:
    """Get the docstring for the object.

    This tries to obtain the docstring for some kind of objects additionally:

    * partial functions
    * inherited docstring
    * inherited decorated methods
    """
    if cls and name and isclassmethod(obj, cls, name):
        for basecls in getmro(cls):
            meth = basecls.__dict__.get(name)
            if meth and hasattr(meth, '__func__'):
                doc: str | None = getdoc(meth.__func__)
                if doc is not None or not allow_inherited:
                    return doc

    doc = _getdoc_internal(obj)
    if ispartial(obj) and doc == obj.__class__.__doc__:
        return getdoc(obj.func)
    elif doc is None and allow_inherited:
        if cls and name:
            # Check a docstring of the attribute or method from super classes.
            for basecls in getmro(cls):
                meth = safe_getattr(basecls, name, None)
                if meth is not None:
                    doc = _getdoc_internal(meth)
                    if doc is not None:
                        break

            if doc is None:
                # retry using `inspect.getdoc()`
                for basecls in getmro(cls):
                    meth = safe_getattr(basecls, name, None)
                    if meth is not None:
                        doc = inspect.getdoc(meth)
                        if doc is not None:
                            break

        if doc is None:
            doc = inspect.getdoc(obj)

    return doc


def _getdoc_internal(
    obj: Any, attrgetter: Callable[[Any, str, Any], Any] = safe_getattr
) -> str | None:
    doc = attrgetter(obj, '__doc__', None)
    if isinstance(doc, str):
        return doc
    return None
