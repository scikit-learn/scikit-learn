import inspect

from typing import Union
from typing import Any
from typing import TypeVar
from typing_extensions import Literal  # noqa
from typing_extensions import Annotated  # noqa

import numpy as np


RandomState = Union[int, np.random.RandomState, None]
ArrayLike = TypeVar('ArrayLike')


class Shape:
    def __init__(self, *shapes):
        if any(not isinstance(s, tuple) for s in shapes):
            raise ValueError("All shapes must be tuple")
        self.shapes = shapes

    def _join_one(self, shape):
        if len(shape) == 1:
            return f"({shape[0]},)"
        inner = ', '.join(str(s) for s in shape)
        return f"({inner})"

    def __repr__(self):
        output = ' or '.join(self._join_one(shape) for shape in self.shapes)
        return f"of shape {output}"


def get_annotation_class_name(annotation) -> str:
    if annotation is None:
        return 'None'
    elif annotation is Any:
        return 'Any'
    elif hasattr(annotation, '__metadata__'):
        return 'Annotated'

    if getattr(annotation, '__qualname__', None):
        return annotation.__qualname__
    elif getattr(annotation, '_name', None):
        #  generic for >= 3.7
        return annotation._name

    origin = getattr(annotation, '__origin__', None)
    if origin:
        return get_annotation_class_name(annotation.__origin__)

    if inspect.isclass(annotation):
        annotation = annotation.__class__
    return annotation.__qualname__.lstrip('_')


def format_annotation(annotation):
    """Convert annotation to docstring"""
    class_name = get_annotation_class_name(annotation)

    if class_name == 'BaseEstimator':
        return 'estimator instance'
    elif class_name == 'ArrayLike':
        return 'array-like'
    elif class_name == 'NoneType':
        return 'None'
    elif class_name == 'RandomState':
        return 'RandomState instance'
    elif class_name == 'Annotated':
        inner_annotation = format_annotation(annotation.__origin__)
        args = ', '.join(repr(t) for t in annotation.__metadata__)
        return f'{inner_annotation} {args}'
    elif class_name == 'Union':
        values = [format_annotation(t) for t in annotation.__args__]
        if len(values) == 2:
            return ' or '.join(values)
        # greater than 2
        first = ', '.join(values[:-1])
        return f'{first} or {values[-1]}'
    elif class_name == 'Literal':
        values = ', '.join(repr(t) for t in annotation.__args__)
        return f'{{{values}}}'
    elif class_name in ('list', 'List'):
        values = ', '.join(format_annotation(t)
                           for t in annotation.__args__)
        return f'list of {values}'

    return class_name


def get_annotations(instance):
    if not hasattr(instance, '__annotations__'):
        raise ValueError(f"{instance} does not have annotations")

    annotations = instance.__annotations__
    # get defaults
    params = inspect.signature(instance).parameters
    defaults = {p: v.default for p, v in params.items()
                if v.default != inspect.Parameter.empty}

    output = {}
    for name, annotation in annotations.items():
        anno = format_annotation(annotation)
        if name in defaults:
            anno += f", default={repr(defaults[name])}"
        output[name] = anno
    return output
