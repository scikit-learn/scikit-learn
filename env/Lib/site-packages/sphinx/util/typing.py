"""The composite types for Sphinx."""

from __future__ import annotations

import dataclasses
import sys
import types
import typing
from collections.abc import Callable, Sequence
from contextvars import Context, ContextVar, Token
from struct import Struct
from typing import (
    TYPE_CHECKING,
    Annotated,
    Any,
    ForwardRef,
    NewType,
    TypedDict,
    TypeVar,
    Union,
)

from docutils import nodes
from docutils.parsers.rst.states import Inliner

from sphinx.util import logging

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Final, Literal, Protocol, TypeAlias

    from typing_extensions import TypeIs

    from sphinx.application import Sphinx

    _RestifyMode: TypeAlias = Literal[
        'fully-qualified-except-typing',
        'smart',
    ]
    _StringifyMode: TypeAlias = Literal[
        'fully-qualified-except-typing',
        'fully-qualified',
        'smart',
    ]

logger = logging.getLogger(__name__)


# classes that have an incorrect .__module__ attribute
_INVALID_BUILTIN_CLASSES: Final[Mapping[object, str]] = {
    Context: 'contextvars.Context',  # Context.__module__ == '_contextvars'
    ContextVar: 'contextvars.ContextVar',  # ContextVar.__module__ == '_contextvars'
    Token: 'contextvars.Token',  # Token.__module__ == '_contextvars'
    Struct: 'struct.Struct',  # Struct.__module__ == '_struct'
    # types in 'types' with <type>.__module__ == 'builtins':
    types.AsyncGeneratorType: 'types.AsyncGeneratorType',
    types.BuiltinFunctionType: 'types.BuiltinFunctionType',
    types.BuiltinMethodType: 'types.BuiltinMethodType',
    types.CellType: 'types.CellType',
    types.ClassMethodDescriptorType: 'types.ClassMethodDescriptorType',
    types.CodeType: 'types.CodeType',
    types.CoroutineType: 'types.CoroutineType',
    types.FrameType: 'types.FrameType',
    types.FunctionType: 'types.FunctionType',
    types.GeneratorType: 'types.GeneratorType',
    types.GetSetDescriptorType: 'types.GetSetDescriptorType',
    types.LambdaType: 'types.LambdaType',
    types.MappingProxyType: 'types.MappingProxyType',
    types.MemberDescriptorType: 'types.MemberDescriptorType',
    types.MethodDescriptorType: 'types.MethodDescriptorType',
    types.MethodType: 'types.MethodType',
    types.MethodWrapperType: 'types.MethodWrapperType',
    types.ModuleType: 'types.ModuleType',
    types.TracebackType: 'types.TracebackType',
    types.WrapperDescriptorType: 'types.WrapperDescriptorType',
}


def is_invalid_builtin_class(obj: Any) -> bool:
    """Check *obj* is an invalid built-in class."""
    try:
        return obj in _INVALID_BUILTIN_CLASSES
    except TypeError:  # unhashable type
        return False


# Text like nodes which are initialized with text and rawsource
TextlikeNode: TypeAlias = nodes.Text | nodes.TextElement

# path matcher
PathMatcher: TypeAlias = Callable[[str], bool]

# common role functions
if TYPE_CHECKING:

    class RoleFunction(Protocol):
        def __call__(  # NoQA: E704
            self,
            name: str,
            rawtext: str,
            text: str,
            lineno: int,
            inliner: Inliner,
            /,
            options: dict[str, Any] | None = None,
            content: Sequence[str] = (),
        ) -> tuple[list[nodes.Node], list[nodes.system_message]]: ...

else:
    RoleFunction: TypeAlias = Callable[
        [str, str, str, int, Inliner, dict[str, Any], Sequence[str]],
        tuple[list[nodes.Node], list[nodes.system_message]],
    ]

# A option spec for directive
OptionSpec: TypeAlias = dict[str, Callable[[str], Any]]

# title getter functions for enumerable nodes (see sphinx.domains.std)
TitleGetter: TypeAlias = Callable[[nodes.Node], str]

# Readable file stream for inventory loading
if TYPE_CHECKING:
    from types import TracebackType

    from typing_extensions import Self

    _T_co = TypeVar('_T_co', str, bytes, covariant=True)

    class _ReadableStream(Protocol[_T_co]):
        def read(self, size: int = ...) -> _T_co: ...  # NoQA: E704

        def __enter__(self) -> Self: ...  # NoQA: E704

        def __exit__(  # NoQA: E704
            self,
            exc_type: type[BaseException] | None,
            exc_val: BaseException | None,
            exc_tb: TracebackType | None,
        ) -> None: ...


# inventory data on memory
InventoryItem: TypeAlias = tuple[
    str,  # project name
    str,  # project version
    str,  # URL
    str,  # display name
]
Inventory: TypeAlias = dict[str, dict[str, InventoryItem]]


class ExtensionMetadata(TypedDict, total=False):
    """The metadata returned by an extension's ``setup()`` function.

    See :ref:`ext-metadata`.
    """

    version: str
    """The extension version (default: ``'unknown version'``)."""
    env_version: int
    """An integer that identifies the version of env data added by the extension."""
    parallel_read_safe: bool
    """Indicate whether parallel reading of source files is supported
    by the extension.
    """
    parallel_write_safe: bool
    """Indicate whether parallel writing of output files is supported
    by the extension (default: ``True``).
    """


if TYPE_CHECKING:
    _ExtensionSetupFunc: TypeAlias = Callable[[Sphinx], ExtensionMetadata]


def get_type_hints(
    obj: Any,
    globalns: dict[str, Any] | None = None,
    localns: dict[str, Any] | None = None,
    include_extras: bool = False,
) -> dict[str, Any]:
    """Return a dictionary containing type hints for a function, method, module or class
    object.

    This is a simple wrapper of `typing.get_type_hints()` that does not raise an error on
    runtime.
    """
    from sphinx.util.inspect import safe_getattr  # lazy loading

    try:
        return typing.get_type_hints(
            obj, globalns, localns, include_extras=include_extras
        )
    except NameError:
        # Failed to evaluate ForwardRef (maybe TYPE_CHECKING)
        return safe_getattr(obj, '__annotations__', {})
    except AttributeError:
        # Failed to evaluate ForwardRef (maybe not runtime checkable)
        return safe_getattr(obj, '__annotations__', {})
    except TypeError:
        # Invalid object is given. But try to get __annotations__ as a fallback.
        return safe_getattr(obj, '__annotations__', {})
    except KeyError:
        # a broken class found (refs: https://github.com/sphinx-doc/sphinx/issues/8084)
        return {}


def is_system_TypeVar(typ: Any) -> bool:
    """Check *typ* is system defined TypeVar."""
    modname = getattr(typ, '__module__', '')
    return modname == 'typing' and isinstance(typ, TypeVar)


def _is_annotated_form(obj: Any) -> TypeIs[Annotated[Any, ...]]:
    """Check if *obj* is an annotated type."""
    return (
        typing.get_origin(obj) is Annotated
        or str(obj).startswith('typing.Annotated')
    )  # fmt: skip


def _is_unpack_form(obj: Any) -> bool:
    """Check if the object is :class:`typing.Unpack` or equivalent."""
    if sys.version_info >= (3, 11):
        from typing import Unpack

        # typing_extensions.Unpack != typing.Unpack for 3.11, but we assume
        # that typing_extensions.Unpack should not be used in that case
        return typing.get_origin(obj) is Unpack

    # Python 3.10 requires typing_extensions.Unpack
    origin = typing.get_origin(obj)
    return (
        getattr(origin, '__module__', None) == 'typing_extensions'
        and origin.__name__ == 'Unpack'
    )


def restify(cls: Any, mode: _RestifyMode = 'fully-qualified-except-typing') -> str:
    """Convert a type-like object to a reST reference.

    :param mode: Specify a method how annotations will be stringified.

                 'fully-qualified-except-typing'
                     Show the module name and qualified name of the annotation except
                     the "typing" module.
                 'smart'
                     Show the name of the annotation.
    """
    from sphinx.ext.autodoc.mock import ismock, ismockmodule  # lazy loading
    from sphinx.util.inspect import isgenericalias, object_description  # lazy loading

    valid_modes = {'fully-qualified-except-typing', 'smart'}
    if mode not in valid_modes:
        valid = ', '.join(map(repr, sorted(valid_modes)))
        msg = f'mode must be one of {valid}; got {mode!r}'
        raise ValueError(msg)

    # things that are not types
    if cls is None or cls == types.NoneType:
        return ':py:obj:`None`'
    if cls is Ellipsis:
        return '...'
    if isinstance(cls, str):
        return cls

    cls_module_is_typing = getattr(cls, '__module__', '') == 'typing'

    # If the mode is 'smart', we always use '~'.
    # If the mode is 'fully-qualified-except-typing',
    # we use '~' only for the objects in the ``typing`` module.
    module_prefix = '~' if mode == 'smart' or cls_module_is_typing else ''

    try:
        if ismockmodule(cls):
            return f':py:class:`{module_prefix}{cls.__name__}`'
        elif ismock(cls):
            return f':py:class:`{module_prefix}{cls.__module__}.{cls.__name__}`'
        elif is_invalid_builtin_class(cls):
            # The above predicate never raises TypeError but should not be
            # evaluated before determining whether *cls* is a mocked object
            # or not; instead of two try-except blocks, we keep it here.
            return f':py:class:`{module_prefix}{_INVALID_BUILTIN_CLASSES[cls]}`'
        elif _is_annotated_form(cls):
            args = restify(cls.__args__[0], mode)
            meta_args = []
            for m in cls.__metadata__:
                if isinstance(m, type):
                    meta_args.append(restify(m, mode))
                elif dataclasses.is_dataclass(m):
                    # use restify for the repr of field values rather than repr
                    d_fields = ', '.join([
                        rf'{f.name}=\ {restify(getattr(m, f.name), mode)}'
                        for f in dataclasses.fields(m)
                        if f.repr
                    ])
                    meta_args.append(rf'{restify(type(m), mode)}\ ({d_fields})')
                else:
                    meta_args.append(repr(m))
            meta = ', '.join(meta_args)
            if sys.version_info[:2] <= (3, 11):
                # Hardcoded to fix errors on Python 3.11 and earlier.
                return rf':py:class:`~typing.Annotated`\ [{args}, {meta}]'
            return (
                f':py:class:`{module_prefix}{cls.__module__}.{cls.__name__}`'
                rf'\ [{args}, {meta}]'
            )
        elif isinstance(cls, NewType):
            return f':py:class:`{module_prefix}{cls.__module__}.{cls.__name__}`'  # type: ignore[attr-defined]
        elif isinstance(cls, types.UnionType):
            # Union types (PEP 585) retain their definition order when they
            # are printed natively and ``None``-like types are kept as is.
            return ' | '.join(restify(a, mode) for a in cls.__args__)
        elif cls.__module__ in ('__builtin__', 'builtins'):
            if hasattr(cls, '__args__'):
                if not cls.__args__:  # Empty tuple, list, ...
                    return rf':py:class:`{cls.__name__}`\ [{cls.__args__!r}]'

                concatenated_args = ', '.join(
                    restify(arg, mode) for arg in cls.__args__
                )
                return rf':py:class:`{cls.__name__}`\ [{concatenated_args}]'
            return f':py:class:`{cls.__name__}`'
        elif isgenericalias(cls) and cls_module_is_typing and cls.__origin__ is Union:
            # *cls* is defined in ``typing``, and thus ``__args__`` must exist
            return ' | '.join(restify(a, mode) for a in cls.__args__)
        elif isgenericalias(cls):
            if isinstance(cls.__origin__, typing._SpecialForm):
                # ClassVar; Concatenate; Final; Literal; Unpack; TypeGuard; TypeIs
                # Required/NotRequired
                text = restify(cls.__origin__, mode)
            elif cls.__name__:
                text = f':py:class:`{module_prefix}{cls.__module__}.{cls.__name__}`'
            else:
                text = restify(cls.__origin__, mode)

            __args__ = getattr(cls, '__args__', ())
            if not __args__:
                return text
            if all(map(is_system_TypeVar, __args__)):
                # Don't print the arguments; they're all system defined type variables.
                return text

            # Callable has special formatting
            if (
                (cls_module_is_typing and cls.__name__ == 'Callable')
                or (cls.__module__ == 'collections.abc' and cls.__name__ == 'Callable')
            ):  # fmt: skip
                args = ', '.join(restify(a, mode) for a in __args__[:-1])
                returns = restify(__args__[-1], mode)
                return rf'{text}\ [[{args}], {returns}]'

            if cls_module_is_typing and cls.__origin__.__name__ == 'Literal':
                args = ', '.join(
                    _format_literal_arg_restify(a, mode=mode) for a in cls.__args__
                )
                return rf'{text}\ [{args}]'

            # generic representation of the parameters
            args = ', '.join(restify(a, mode) for a in __args__)
            return rf'{text}\ [{args}]'
        elif isinstance(cls, typing._SpecialForm):
            return f':py:obj:`~{cls.__module__}.{cls.__name__}`'  # type: ignore[attr-defined]
        elif sys.version_info[:2] >= (3, 11) and cls is typing.Any:
            # handle bpo-46998
            return f':py:obj:`~{cls.__module__}.{cls.__name__}`'
        elif hasattr(cls, '__qualname__'):
            return f':py:class:`{module_prefix}{cls.__module__}.{cls.__qualname__}`'
        elif isinstance(cls, ForwardRef):
            return f':py:class:`{cls.__forward_arg__}`'
        else:
            # not a class (ex. TypeVar) but should have a __name__
            return f':py:obj:`{module_prefix}{cls.__module__}.{cls.__name__}`'
    except (AttributeError, TypeError) as exc:
        logger.debug('restify on %r in mode %r failed: %r', cls, mode, exc)
        return object_description(cls)


def _format_literal_arg_restify(arg: Any, /, *, mode: str) -> str:
    from sphinx.util.inspect import isenumattribute  # lazy loading

    if isenumattribute(arg):
        enum_cls = arg.__class__
        if mode == 'smart' or enum_cls.__module__ == 'typing':
            # MyEnum.member
            return (
                f':py:attr:`~{enum_cls.__module__}.{enum_cls.__qualname__}.{arg.name}`'
            )
        # module.MyEnum.member
        return f':py:attr:`{enum_cls.__module__}.{enum_cls.__qualname__}.{arg.name}`'
    return repr(arg)


def stringify_annotation(
    annotation: Any,
    /,
    mode: _StringifyMode = 'fully-qualified-except-typing',
) -> str:
    """Stringify type annotation object.

    :param annotation: The annotation to stringified.
    :param mode: Specify a method how annotations will be stringified.

                 'fully-qualified-except-typing'
                     Show the module name and qualified name of the annotation except
                     the "typing" module.
                 'smart'
                     Show the name of the annotation.
                 'fully-qualified'
                     Show the module name and qualified name of the annotation.
    """
    from sphinx.ext.autodoc.mock import ismock, ismockmodule  # lazy loading

    valid_modes = {'fully-qualified-except-typing', 'fully-qualified', 'smart'}
    if mode not in valid_modes:
        valid = ', '.join(map(repr, sorted(valid_modes)))
        msg = f'mode must be one of {valid}; got {mode!r}'
        raise ValueError(msg)

    # things that are not types
    if annotation is None or annotation == types.NoneType:
        return 'None'
    if annotation is Ellipsis:
        return '...'
    if isinstance(annotation, str):
        if annotation.startswith("'") and annotation.endswith("'"):
            # Might be a double Forward-ref'ed type.  Go unquoting.
            return annotation[1:-1]
        return annotation
    if not annotation:
        return repr(annotation)

    module_prefix = '~' if mode == 'smart' else ''

    # The values below must be strings if the objects are well-formed.
    annotation_qualname: str = getattr(annotation, '__qualname__', '')
    annotation_module: str = getattr(annotation, '__module__', '')
    annotation_name: str = getattr(annotation, '__name__', '')
    annotation_module_is_typing = annotation_module == 'typing'

    # Extract the annotation's base type by considering formattable cases
    if isinstance(annotation, TypeVar) and not _is_unpack_form(annotation):
        # typing_extensions.Unpack is incorrectly determined as a TypeVar
        if annotation_module_is_typing and mode in {
            'fully-qualified-except-typing',
            'smart',
        }:
            return annotation_name
        return module_prefix + f'{annotation_module}.{annotation_name}'
    elif isinstance(annotation, NewType):
        return module_prefix + f'{annotation_module}.{annotation_name}'
    elif ismockmodule(annotation):
        return module_prefix + annotation_name
    elif ismock(annotation):
        return module_prefix + f'{annotation_module}.{annotation_name}'
    elif is_invalid_builtin_class(annotation):
        return module_prefix + _INVALID_BUILTIN_CLASSES[annotation]
    elif _is_annotated_form(annotation):  # for py310+
        pass
    elif annotation_module == 'builtins' and annotation_qualname:
        args = getattr(annotation, '__args__', None)
        if args is None:
            return annotation_qualname

        # PEP 585 generic
        if not args:  # Empty tuple, list, ...
            return repr(annotation)

        concatenated_args = ', '.join(stringify_annotation(arg, mode) for arg in args)
        return f'{annotation_qualname}[{concatenated_args}]'
    else:
        # add other special cases that can be directly formatted
        pass

    module_prefix = f'{annotation_module}.'
    annotation_forward_arg: str | None = getattr(annotation, '__forward_arg__', None)
    if annotation_qualname or (
        annotation_module_is_typing and not annotation_forward_arg
    ):
        if mode == 'smart':
            module_prefix = f'~{module_prefix}'
        if annotation_module_is_typing and mode == 'fully-qualified-except-typing':
            module_prefix = ''
    elif _is_unpack_form(annotation) and annotation_module == 'typing_extensions':
        module_prefix = '~' if mode == 'smart' else ''
    else:
        module_prefix = ''

    if annotation_module_is_typing:
        if annotation_forward_arg:
            # handle ForwardRefs
            qualname = annotation_forward_arg
        else:
            if annotation_name:
                qualname = annotation_name
            elif annotation_qualname:
                qualname = annotation_qualname
            else:
                # in this case, we know that the annotation is a member
                # of ``typing`` and all of them define ``__origin__``
                qualname = stringify_annotation(
                    annotation.__origin__,
                    'fully-qualified-except-typing',
                ).replace('typing.', '')  # ex. Union
    elif annotation_qualname:
        qualname = annotation_qualname
    elif hasattr(annotation, '__origin__'):
        # instantiated generic provided by a user
        qualname = stringify_annotation(annotation.__origin__, mode)
    elif isinstance(annotation, types.UnionType):
        qualname = 'types.UnionType'
    else:
        # we weren't able to extract the base type, appending arguments would
        # only make them appear twice
        return repr(annotation)

    # Process the generic arguments (if any).
    # They must be a list or a tuple, otherwise they are considered 'broken'.
    annotation_args = getattr(annotation, '__args__', ())
    if annotation_args and isinstance(annotation_args, list | tuple):
        if (
            qualname in {'Union', 'types.UnionType'}
            and all(getattr(a, '__origin__', ...) is typing.Literal for a in annotation_args)
        ):  # fmt: skip
            # special case to flatten a Union of Literals into a literal
            flattened_args = typing.Literal[annotation_args].__args__  # type: ignore[attr-defined]
            args = ', '.join(
                _format_literal_arg_stringify(a, mode=mode) for a in flattened_args
            )
            return f'{module_prefix}Literal[{args}]'
        if qualname in {'Optional', 'Union', 'types.UnionType'}:
            return ' | '.join(stringify_annotation(a, mode) for a in annotation_args)
        elif qualname == 'Callable':
            args = ', '.join(
                stringify_annotation(a, mode) for a in annotation_args[:-1]
            )
            returns = stringify_annotation(annotation_args[-1], mode)
            return f'{module_prefix}Callable[[{args}], {returns}]'
        elif qualname == 'Literal':
            args = ', '.join(
                _format_literal_arg_stringify(a, mode=mode) for a in annotation_args
            )
            return f'{module_prefix}Literal[{args}]'
        elif _is_annotated_form(annotation):  # for py310+
            args = stringify_annotation(annotation_args[0], mode)
            meta_args = []
            for m in annotation.__metadata__:
                if isinstance(m, type):
                    meta_args.append(stringify_annotation(m, mode))
                elif dataclasses.is_dataclass(m):
                    # use stringify_annotation for the repr of field values rather than repr
                    d_fields = ', '.join([
                        f'{f.name}={stringify_annotation(getattr(m, f.name), mode)}'
                        for f in dataclasses.fields(m)
                        if f.repr
                    ])
                    meta_args.append(
                        f'{stringify_annotation(type(m), mode)}({d_fields})'
                    )
                else:
                    meta_args.append(repr(m))
            meta = ', '.join(meta_args)
            if sys.version_info[:2] <= (3, 11):
                if mode == 'fully-qualified-except-typing':
                    return f'Annotated[{args}, {meta}]'
                module_prefix = module_prefix.replace('builtins', 'typing')
                return f'{module_prefix}Annotated[{args}, {meta}]'
            return f'{module_prefix}Annotated[{args}, {meta}]'
        elif all(is_system_TypeVar(a) for a in annotation_args):
            # Suppress arguments if all system defined TypeVars (ex. Dict[KT, VT])
            return module_prefix + qualname
        else:
            args = ', '.join(stringify_annotation(a, mode) for a in annotation_args)
            return f'{module_prefix}{qualname}[{args}]'

    return module_prefix + qualname


def _format_literal_arg_stringify(arg: Any, /, *, mode: str) -> str:
    from sphinx.util.inspect import isenumattribute  # lazy loading

    if isenumattribute(arg):
        enum_cls = arg.__class__
        if mode == 'smart' or enum_cls.__module__ == 'typing':
            # MyEnum.member
            return f'{enum_cls.__qualname__}.{arg.name}'
        # module.MyEnum.member
        return f'{enum_cls.__module__}.{enum_cls.__qualname__}.{arg.name}'
    return repr(arg)


# deprecated name -> (object to return, canonical path or empty string, removal version)
_DEPRECATED_OBJECTS: dict[str, tuple[Any, str, tuple[int, int]]] = {
}  # fmt: skip


def __getattr__(name: str) -> Any:
    if name not in _DEPRECATED_OBJECTS:
        msg = f'module {__name__!r} has no attribute {name!r}'
        raise AttributeError(msg)

    from sphinx.deprecation import _deprecation_warning

    deprecated_object, canonical_name, remove = _DEPRECATED_OBJECTS[name]
    _deprecation_warning(__name__, name, canonical_name, remove=remove)
    return deprecated_object
