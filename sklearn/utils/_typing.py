import inspect
import numbers
import typing
from typing import Union
from typing import TypeVar
from typing import Iterator

import numpy as np


try:
    import typing_extensions  # noqa

    TYPING_EXTENSION_INSTALLED = True
except ImportError:
    TYPING_EXTENSION_INSTALLED = False


if typing.TYPE_CHECKING or TYPING_EXTENSION_INSTALLED:
    from typing_extensions import Literal  # noqa
else:

    class _SimpleLiteral:
        def __getitem__(self, values):
            return typing.Any

    Literal = _SimpleLiteral()


if typing.TYPE_CHECKING or TYPING_EXTENSION_INSTALLED:
    from typing_extensions import Protocol  # noqa

    class CVSplitter(Protocol):
        def get_n_splits(self):
            """Get the number of splits."""

        def split(self, X, y=None, groups=None):
            """Split data"""

else:
    CVSplitter = TypeVar("CVSplitter")  # typing: ignore

CVType = Union[int, CVSplitter, Iterator, None]
ArrayLike = TypeVar("ArrayLike")
NDArray = TypeVar("NDArray")
EstimatorType = TypeVar("EstimatorType")
JoblibMemory = TypeVar("JoblibMemory")
RandomStateType = Union[int, np.random.RandomState, None]
MemoryType = Union[str, JoblibMemory, None]


def get_annotation_class_name(annotation) -> str:
    # Special cases
    if annotation is None or annotation is type(None):  # noqa
        return "None"

    if getattr(annotation, "__qualname__", None):
        return annotation.__qualname__
    elif getattr(
        annotation, "_name", None
    ):  # Required for generic aliases on Python 3.7+
        return annotation._name

    origin = getattr(annotation, "__origin__", None)
    if origin:
        if getattr(origin, "__qualname__", None):
            # Required for Protocol subclasses
            return origin.__qualname__
        elif getattr(origin, "_name", None):
            # Required for Union on Python 3.7+
            return origin._name
        else:
            return origin.__class__.__qualname__.lstrip(
                "_"
            )  # Required for Union on Python < 3.7

    annotation_cls = (annotation
                      if inspect.isclass(annotation) else annotation.__class__)
    return annotation_cls.__qualname__.lstrip("_")


def format_docstring_annotation(annotation):
    """Convert annotation to docstring."""

    # handle some annotations directly
    if annotation == np.random.RandomState:
        return "RandomState instance"
    elif annotation == CVSplitter:
        return "cross-validation generator"
    elif annotation == Iterator:
        return "iterable"

    class_name = get_annotation_class_name(annotation)

    if class_name == "None":
        return "None"
    elif class_name == "TypeVar":
        name = annotation.__name__
        if name == "EstimatorType":
            return "estimator instance"
        elif name == "ArrayLike":
            return "array-like"
        elif name == "NDArray":
            return "ndarray"
        elif name == "JoblibMemory":
            return "object with the joblib.Memory interface"
        else:
            raise ValueError(f"Unrecognized TypeVar: {annotation}")

    elif class_name == 'Union':
        values = [format_docstring_annotation(t)
                  for t in annotation.__args__]
        if len(values) == 2:
            return ' or '.join(values)
        first = ', '.join(values[:-1])
        return f'{first} or {values[-1]}'

    elif class_name == 'Literal':
        if hasattr(annotation, '__values__'):
            # For Python == 3.6 support
            args = annotation.__values__
        else:
            args = annotation.__args__
        items = [repr(t) for t in args]
        if len(items) == 1:
            return items[0]
        values = ', '.join(items)
        return f'{{{values}}}'

    elif class_name == 'List':
        values = ', '.join(format_docstring_annotation(t)
                           for t in annotation.__args__)
        return f'list of {values}'

    return class_name


def get_docstring_annotations(obj):
    """Get human readable docstring for types for a obj with annotations.

    This function requires `typing_extensions` to be installed to run.

    Parameters
    ----------
    obj: object
        Object to get annotations from

    Returns
    -------
    output: dict
        dictionary mapping from name to human-readable docstring.
    """
    if not hasattr(obj, '__annotations__'):
        return {}

    annotations = typing.get_type_hints(obj)
    # get defaults
    params = inspect.signature(obj).parameters
    defaults = {p: v.default for p, v in params.items()
                if v.default != inspect.Parameter.empty}

    output = {}
    for name, annotation in annotations.items():
        anno = format_docstring_annotation(annotation)
        if name in defaults:
            default = defaults[name]
            if (isinstance(default, numbers.Real) and
                    not isinstance(default, numbers.Integral)):
                # For floats the representation can vary, i.e:
                # default=np.inf or default=1e-4
                anno += ", default="
            else:
                anno += f", default={repr(default)}"
        output[name] = anno
    return output
