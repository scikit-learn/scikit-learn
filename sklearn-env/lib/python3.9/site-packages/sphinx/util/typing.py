"""
    sphinx.util.typing
    ~~~~~~~~~~~~~~~~~~

    The composite types for Sphinx.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import sys
import typing
from struct import Struct
from types import TracebackType
from typing import Any, Callable, Dict, Generator, List, Optional, Tuple, Type, TypeVar, Union

from docutils import nodes
from docutils.parsers.rst.states import Inliner

from sphinx.deprecation import RemovedInSphinx60Warning, deprecated_alias

if sys.version_info > (3, 7):
    from typing import ForwardRef
else:
    from typing import _ForwardRef  # type: ignore

    class ForwardRef:
        """A pseudo ForwardRef class for py36."""
        def __init__(self, arg: Any, is_argument: bool = True) -> None:
            self.arg = arg

        def _evaluate(self, globalns: Dict, localns: Dict) -> Any:
            ref = _ForwardRef(self.arg)
            return ref._eval_type(globalns, localns)

try:
    from types import UnionType  # type: ignore  # python 3.10 or above
except ImportError:
    UnionType = None

if False:
    # For type annotation
    from typing import Type  # NOQA # for python3.5.1


# builtin classes that have incorrect __module__
INVALID_BUILTIN_CLASSES = {
    Struct: 'struct.Struct',  # Before Python 3.9
    TracebackType: 'types.TracebackType',
}


# Text like nodes which are initialized with text and rawsource
TextlikeNode = Union[nodes.Text, nodes.TextElement]

# type of None
NoneType = type(None)

# path matcher
PathMatcher = Callable[[str], bool]

# common role functions
RoleFunction = Callable[[str, str, str, int, Inliner, Dict[str, Any], List[str]],
                        Tuple[List[nodes.Node], List[nodes.system_message]]]

# A option spec for directive
OptionSpec = Dict[str, Callable[[str], Any]]

# title getter functions for enumerable nodes (see sphinx.domains.std)
TitleGetter = Callable[[nodes.Node], str]

# inventory data on memory
InventoryItem = Tuple[str, str, str, str]
Inventory = Dict[str, Dict[str, InventoryItem]]


def get_type_hints(obj: Any, globalns: Dict = None, localns: Dict = None) -> Dict[str, Any]:
    """Return a dictionary containing type hints for a function, method, module or class object.

    This is a simple wrapper of `typing.get_type_hints()` that does not raise an error on
    runtime.
    """
    from sphinx.util.inspect import safe_getattr  # lazy loading

    try:
        return typing.get_type_hints(obj, globalns, localns)
    except NameError:
        # Failed to evaluate ForwardRef (maybe TYPE_CHECKING)
        return safe_getattr(obj, '__annotations__', {})
    except AttributeError:
        # Failed to evaluate ForwardRef (maybe not runtime checkable)
        return safe_getattr(obj, '__annotations__', {})
    except TypeError:
        # Invalid object is given. But try to get __annotations__ as a fallback for
        # the code using type union operator (PEP 604) in python 3.9 or below.
        return safe_getattr(obj, '__annotations__', {})
    except KeyError:
        # a broken class found (refs: https://github.com/sphinx-doc/sphinx/issues/8084)
        return {}


def is_system_TypeVar(typ: Any) -> bool:
    """Check *typ* is system defined TypeVar."""
    modname = getattr(typ, '__module__', '')
    return modname == 'typing' and isinstance(typ, TypeVar)


def restify(cls: Optional[Type], mode: str = 'fully-qualified-except-typing') -> str:
    """Convert python class to a reST reference.

    :param mode: Specify a method how annotations will be stringified.

                 'fully-qualified-except-typing'
                     Show the module name and qualified name of the annotation except
                     the "typing" module.
                 'smart'
                     Show the name of the annotation.
    """
    from sphinx.util import inspect  # lazy loading

    if mode == 'smart':
        modprefix = '~'
    else:
        modprefix = ''

    try:
        if cls is None or cls is NoneType:
            return ':py:obj:`None`'
        elif cls is Ellipsis:
            return '...'
        elif isinstance(cls, str):
            return cls
        elif cls in INVALID_BUILTIN_CLASSES:
            return ':py:class:`%s%s`' % (modprefix, INVALID_BUILTIN_CLASSES[cls])
        elif inspect.isNewType(cls):
            if sys.version_info > (3, 10):
                # newtypes have correct module info since Python 3.10+
                return ':py:class:`%s%s.%s`' % (modprefix, cls.__module__, cls.__name__)
            else:
                return ':py:class:`%s`' % cls.__name__
        elif UnionType and isinstance(cls, UnionType):
            if len(cls.__args__) > 1 and None in cls.__args__:
                args = ' | '.join(restify(a, mode) for a in cls.__args__ if a)
                return 'Optional[%s]' % args
            else:
                return ' | '.join(restify(a, mode) for a in cls.__args__)
        elif cls.__module__ in ('__builtin__', 'builtins'):
            if hasattr(cls, '__args__'):
                return ':py:class:`%s`\\ [%s]' % (
                    cls.__name__,
                    ', '.join(restify(arg, mode) for arg in cls.__args__),
                )
            else:
                return ':py:class:`%s`' % cls.__name__
        else:
            if sys.version_info >= (3, 7):  # py37+
                return _restify_py37(cls, mode)
            else:
                return _restify_py36(cls, mode)
    except (AttributeError, TypeError):
        return inspect.object_description(cls)


def _restify_py37(cls: Optional[Type], mode: str = 'fully-qualified-except-typing') -> str:
    """Convert python class to a reST reference."""
    from sphinx.util import inspect  # lazy loading

    if mode == 'smart':
        modprefix = '~'
    else:
        modprefix = ''

    if (inspect.isgenericalias(cls) and
            cls.__module__ == 'typing' and cls.__origin__ is Union):
        # Union
        if len(cls.__args__) > 1 and cls.__args__[-1] is NoneType:
            if len(cls.__args__) > 2:
                args = ', '.join(restify(a, mode) for a in cls.__args__[:-1])
                return ':py:obj:`~typing.Optional`\\ [:obj:`~typing.Union`\\ [%s]]' % args
            else:
                return ':py:obj:`~typing.Optional`\\ [%s]' % restify(cls.__args__[0], mode)
        else:
            args = ', '.join(restify(a, mode) for a in cls.__args__)
            return ':py:obj:`~typing.Union`\\ [%s]' % args
    elif inspect.isgenericalias(cls):
        if isinstance(cls.__origin__, typing._SpecialForm):
            text = restify(cls.__origin__, mode)  # type: ignore
        elif getattr(cls, '_name', None):
            if cls.__module__ == 'typing':
                text = ':py:class:`~%s.%s`' % (cls.__module__, cls._name)
            else:
                text = ':py:class:`%s%s.%s`' % (modprefix, cls.__module__, cls._name)
        else:
            text = restify(cls.__origin__, mode)

        origin = getattr(cls, '__origin__', None)
        if not hasattr(cls, '__args__'):
            pass
        elif all(is_system_TypeVar(a) for a in cls.__args__):
            # Suppress arguments if all system defined TypeVars (ex. Dict[KT, VT])
            pass
        elif cls.__module__ == 'typing' and cls._name == 'Callable':
            args = ', '.join(restify(a, mode) for a in cls.__args__[:-1])
            text += r"\ [[%s], %s]" % (args, restify(cls.__args__[-1], mode))
        elif cls.__module__ == 'typing' and getattr(origin, '_name', None) == 'Literal':
            text += r"\ [%s]" % ', '.join(repr(a) for a in cls.__args__)
        elif cls.__args__:
            text += r"\ [%s]" % ", ".join(restify(a, mode) for a in cls.__args__)

        return text
    elif isinstance(cls, typing._SpecialForm):
        return ':py:obj:`~%s.%s`' % (cls.__module__, cls._name)
    elif hasattr(cls, '__qualname__'):
        if cls.__module__ == 'typing':
            return ':py:class:`~%s.%s`' % (cls.__module__, cls.__qualname__)
        else:
            return ':py:class:`%s%s.%s`' % (modprefix, cls.__module__, cls.__qualname__)
    elif isinstance(cls, ForwardRef):
        return ':py:class:`%s`' % cls.__forward_arg__
    else:
        # not a class (ex. TypeVar)
        if cls.__module__ == 'typing':
            return ':py:obj:`~%s.%s`' % (cls.__module__, cls.__name__)
        else:
            return ':py:obj:`%s%s.%s`' % (modprefix, cls.__module__, cls.__name__)


def _restify_py36(cls: Optional[Type], mode: str = 'fully-qualified-except-typing') -> str:
    if mode == 'smart':
        modprefix = '~'
    else:
        modprefix = ''

    module = getattr(cls, '__module__', None)
    if module == 'typing':
        if getattr(cls, '_name', None):
            qualname = cls._name
        elif getattr(cls, '__qualname__', None):
            qualname = cls.__qualname__
        elif getattr(cls, '__forward_arg__', None):
            qualname = cls.__forward_arg__
        elif getattr(cls, '__origin__', None):
            qualname = stringify(cls.__origin__)  # ex. Union
        else:
            qualname = repr(cls).replace('typing.', '')
    elif hasattr(cls, '__qualname__'):
        qualname = '%s%s.%s' % (modprefix, module, cls.__qualname__)
    else:
        qualname = repr(cls)

    if (isinstance(cls, typing.TupleMeta) and  # type: ignore
            not hasattr(cls, '__tuple_params__')):
        if module == 'typing':
            reftext = ':py:class:`~typing.%s`' % qualname
        else:
            reftext = ':py:class:`%s%s`' % (modprefix, qualname)

        params = cls.__args__
        if params:
            param_str = ', '.join(restify(p, mode) for p in params)
            return reftext + '\\ [%s]' % param_str
        else:
            return reftext
    elif isinstance(cls, typing.GenericMeta):
        if module == 'typing':
            reftext = ':py:class:`~typing.%s`' % qualname
        else:
            reftext = ':py:class:`%s%s`' % (modprefix, qualname)

        if cls.__args__ is None or len(cls.__args__) <= 2:
            params = cls.__args__
        elif cls.__origin__ == Generator:
            params = cls.__args__
        else:  # typing.Callable
            args = ', '.join(restify(arg, mode) for arg in cls.__args__[:-1])
            result = restify(cls.__args__[-1], mode)
            return reftext + '\\ [[%s], %s]' % (args, result)

        if params:
            param_str = ', '.join(restify(p, mode) for p in params)
            return reftext + '\\ [%s]' % (param_str)
        else:
            return reftext
    elif (hasattr(cls, '__origin__') and
          cls.__origin__ is typing.Union):
        params = cls.__args__
        if params is not None:
            if len(params) > 1 and params[-1] is NoneType:
                if len(params) > 2:
                    param_str = ", ".join(restify(p, mode) for p in params[:-1])
                    return (':py:obj:`~typing.Optional`\\ '
                            '[:py:obj:`~typing.Union`\\ [%s]]' % param_str)
                else:
                    return ':py:obj:`~typing.Optional`\\ [%s]' % restify(params[0], mode)
            else:
                param_str = ', '.join(restify(p, mode) for p in params)
                return ':py:obj:`~typing.Union`\\ [%s]' % param_str
        else:
            return ':py:obj:`Union`'
    elif hasattr(cls, '__qualname__'):
        if cls.__module__ == 'typing':
            return ':py:class:`~%s.%s`' % (cls.__module__, cls.__qualname__)
        else:
            return ':py:class:`%s%s.%s`' % (modprefix, cls.__module__, cls.__qualname__)
    elif hasattr(cls, '_name'):
        # SpecialForm
        if cls.__module__ == 'typing':
            return ':py:obj:`~%s.%s`' % (cls.__module__, cls._name)
        else:
            return ':py:obj:`%s%s.%s`' % (modprefix, cls.__module__, cls._name)
    elif hasattr(cls, '__name__'):
        # not a class (ex. TypeVar)
        if cls.__module__ == 'typing':
            return ':py:obj:`~%s.%s`' % (cls.__module__, cls.__name__)
        else:
            return ':py:obj:`%s%s.%s`' % (modprefix, cls.__module__, cls.__name__)
    else:
        # others (ex. Any)
        if cls.__module__ == 'typing':
            return ':py:obj:`~%s.%s`' % (cls.__module__, qualname)
        else:
            return ':py:obj:`%s%s.%s`' % (modprefix, cls.__module__, qualname)


def stringify(annotation: Any, mode: str = 'fully-qualified-except-typing') -> str:
    """Stringify type annotation object.

    :param mode: Specify a method how annotations will be stringified.

                 'fully-qualified-except-typing'
                     Show the module name and qualified name of the annotation except
                     the "typing" module.
                 'smart'
                     Show the name of the annotation.
                 'fully-qualified'
                     Show the module name and qualified name of the annotation.
    """
    from sphinx.util import inspect  # lazy loading

    if mode == 'smart':
        modprefix = '~'
    else:
        modprefix = ''

    if isinstance(annotation, str):
        if annotation.startswith("'") and annotation.endswith("'"):
            # might be a double Forward-ref'ed type.  Go unquoting.
            return annotation[1:-1]
        else:
            return annotation
    elif isinstance(annotation, TypeVar):
        if (annotation.__module__ == 'typing' and
                mode in ('fully-qualified-except-typing', 'smart')):
            return annotation.__name__
        else:
            return modprefix + '.'.join([annotation.__module__, annotation.__name__])
    elif inspect.isNewType(annotation):
        if sys.version_info > (3, 10):
            # newtypes have correct module info since Python 3.10+
            return modprefix + '%s.%s' % (annotation.__module__, annotation.__name__)
        else:
            return annotation.__name__
    elif not annotation:
        return repr(annotation)
    elif annotation is NoneType:
        return 'None'
    elif annotation in INVALID_BUILTIN_CLASSES:
        return modprefix + INVALID_BUILTIN_CLASSES[annotation]
    elif str(annotation).startswith('typing.Annotated'):  # for py310+
        pass
    elif (getattr(annotation, '__module__', None) == 'builtins' and
          getattr(annotation, '__qualname__', None)):
        if hasattr(annotation, '__args__'):  # PEP 585 generic
            return repr(annotation)
        else:
            return annotation.__qualname__
    elif annotation is Ellipsis:
        return '...'

    if sys.version_info >= (3, 7):  # py37+
        return _stringify_py37(annotation, mode)
    else:
        return _stringify_py36(annotation, mode)


def _stringify_py37(annotation: Any, mode: str = 'fully-qualified-except-typing') -> str:
    """stringify() for py37+."""
    module = getattr(annotation, '__module__', None)
    modprefix = ''
    if module == 'typing' and getattr(annotation, '__forward_arg__', None):
        qualname = annotation.__forward_arg__
    elif module == 'typing':
        if getattr(annotation, '_name', None):
            qualname = annotation._name
        elif getattr(annotation, '__qualname__', None):
            qualname = annotation.__qualname__
        else:
            qualname = stringify(annotation.__origin__).replace('typing.', '')  # ex. Union

        if mode == 'smart':
            modprefix = '~%s.' % module
        elif mode == 'fully-qualified':
            modprefix = '%s.' % module
    elif hasattr(annotation, '__qualname__'):
        if mode == 'smart':
            modprefix = '~%s.' % module
        else:
            modprefix = '%s.' % module
        qualname = annotation.__qualname__
    elif hasattr(annotation, '__origin__'):
        # instantiated generic provided by a user
        qualname = stringify(annotation.__origin__, mode)
    elif UnionType and isinstance(annotation, UnionType):  # types.Union (for py3.10+)
        qualname = 'types.Union'
    else:
        # we weren't able to extract the base type, appending arguments would
        # only make them appear twice
        return repr(annotation)

    if getattr(annotation, '__args__', None):
        if not isinstance(annotation.__args__, (list, tuple)):
            # broken __args__ found
            pass
        elif qualname in ('Optional', 'Union'):
            if len(annotation.__args__) > 1 and annotation.__args__[-1] is NoneType:
                if len(annotation.__args__) > 2:
                    args = ', '.join(stringify(a, mode) for a in annotation.__args__[:-1])
                    return '%sOptional[%sUnion[%s]]' % (modprefix, modprefix, args)
                else:
                    return '%sOptional[%s]' % (modprefix,
                                               stringify(annotation.__args__[0], mode))
            else:
                args = ', '.join(stringify(a, mode) for a in annotation.__args__)
                return '%sUnion[%s]' % (modprefix, args)
        elif qualname == 'types.Union':
            if len(annotation.__args__) > 1 and None in annotation.__args__:
                args = ' | '.join(stringify(a) for a in annotation.__args__ if a)
                return '%sOptional[%s]' % (modprefix, args)
            else:
                return ' | '.join(stringify(a) for a in annotation.__args__)
        elif qualname == 'Callable':
            args = ', '.join(stringify(a, mode) for a in annotation.__args__[:-1])
            returns = stringify(annotation.__args__[-1], mode)
            return '%s%s[[%s], %s]' % (modprefix, qualname, args, returns)
        elif qualname == 'Literal':
            args = ', '.join(repr(a) for a in annotation.__args__)
            return '%s%s[%s]' % (modprefix, qualname, args)
        elif str(annotation).startswith('typing.Annotated'):  # for py39+
            return stringify(annotation.__args__[0], mode)
        elif all(is_system_TypeVar(a) for a in annotation.__args__):
            # Suppress arguments if all system defined TypeVars (ex. Dict[KT, VT])
            return modprefix + qualname
        else:
            args = ', '.join(stringify(a, mode) for a in annotation.__args__)
            return '%s%s[%s]' % (modprefix, qualname, args)

    return modprefix + qualname


def _stringify_py36(annotation: Any, mode: str = 'fully-qualified-except-typing') -> str:
    """stringify() for py36."""
    module = getattr(annotation, '__module__', None)
    modprefix = ''
    if module == 'typing' and getattr(annotation, '__forward_arg__', None):
        qualname = annotation.__forward_arg__
    elif module == 'typing':
        if getattr(annotation, '_name', None):
            qualname = annotation._name
        elif getattr(annotation, '__qualname__', None):
            qualname = annotation.__qualname__
        elif getattr(annotation, '__origin__', None):
            qualname = stringify(annotation.__origin__)  # ex. Union
        else:
            qualname = repr(annotation).replace('typing.', '')

        if mode == 'smart':
            modprefix = '~%s.' % module
        elif mode == 'fully-qualified':
            modprefix = '%s.' % module
    elif hasattr(annotation, '__qualname__'):
        if mode == 'smart':
            modprefix = '~%s.' % module
        else:
            modprefix = '%s.' % module
        qualname = annotation.__qualname__
    else:
        qualname = repr(annotation)

    if (isinstance(annotation, typing.TupleMeta) and  # type: ignore
            not hasattr(annotation, '__tuple_params__')):  # for Python 3.6
        params = annotation.__args__
        if params:
            param_str = ', '.join(stringify(p, mode) for p in params)
            return '%s%s[%s]' % (modprefix, qualname, param_str)
        else:
            return modprefix + qualname
    elif isinstance(annotation, typing.GenericMeta):
        params = None
        if annotation.__args__ is None or len(annotation.__args__) <= 2:  # type: ignore  # NOQA
            params = annotation.__args__  # type: ignore
        elif annotation.__origin__ == Generator:  # type: ignore
            params = annotation.__args__  # type: ignore
        else:  # typing.Callable
            args = ', '.join(stringify(arg, mode) for arg
                             in annotation.__args__[:-1])  # type: ignore
            result = stringify(annotation.__args__[-1])  # type: ignore
            return '%s%s[[%s], %s]' % (modprefix, qualname, args, result)
        if params is not None:
            param_str = ', '.join(stringify(p, mode) for p in params)
            return '%s%s[%s]' % (modprefix, qualname, param_str)
    elif (hasattr(annotation, '__origin__') and
          annotation.__origin__ is typing.Union):
        params = annotation.__args__
        if params is not None:
            if len(params) > 1 and params[-1] is NoneType:
                if len(params) > 2:
                    param_str = ", ".join(stringify(p, mode) for p in params[:-1])
                    return '%sOptional[%sUnion[%s]]' % (modprefix, modprefix, param_str)
                else:
                    return '%sOptional[%s]' % (modprefix, stringify(params[0], mode))
            else:
                param_str = ', '.join(stringify(p, mode) for p in params)
                return '%sUnion[%s]' % (modprefix, param_str)

    return modprefix + qualname


deprecated_alias('sphinx.util.typing',
                 {
                     'DirectiveOption': Callable[[str], Any],
                 },
                 RemovedInSphinx60Warning)
