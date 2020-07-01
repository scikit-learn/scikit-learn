import re
import numpy as np

from typing import Union
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal  # noqa

try:
    from typing import Annotated  # python 3.9
except ImportError:
    from typing_extensions import Annotated  # noqa


RandomState = Union[int, np.random.RandomState, None]


class Shape:
    def __init__(self, *shapes):
        self.shapes = []
        for shape in shapes:
            if not isinstance(shape, (tuple, list)):
                self.shapes.append((shape, ))
            else:
                self.shapes.append(shape)

    def _join_one(self, shape):
        if len(shape) == 1:
            return f"({shape[0]},)"
        inner = ', '.join(str(s) for s in shape)
        return f"({inner})"

    def __repr__(self):
        output = ' or '.join(self._join_one(shape) for shape in self.shapes)
        return f"of shape {output}"


def format_annotation(annotation):
    """Convert annotation to docstring"""
    if annotation is None or annotation is type(None):  # noqa
        return 'None'

    if annotation in [int, bool, float, str, dict, np.ndarray]:
        return annotation.__qualname__

    if hasattr(annotation, '__name__'):
        name = annotation.__name__
        if name == 'BaseEstimator':
            return 'estimator instance'
        elif name == 'RandomState':
            return 'int, RandomState instance, or None'
        elif name == 'ArrayLike':
            return 'array-like'

    if hasattr(annotation, '__origin__'):
        origin = annotation.__origin__
        if hasattr(annotation, '__metadata__'):  # Annotated
            metadata = ', '.join(str(t) for t in annotation.__metadata__)
            type_info = format_annotation(origin)
            return f'{type_info} {metadata}'

        if getattr(origin, '__qualname__', None):
            name = origin.__qualname__
        elif getattr(origin, '_name', None):
            # Required for Union on Python 3.7+
            name = origin._name
        else:
            # Required for Union on Python < 3.7
            name = origin.__class__.__qualname__.lstrip('_')

        if name == 'Union':
            values = [format_annotation(t) for t in annotation.__args__]
            if len(values) == 2:
                return ' or '.join(values)
            # greater than 2
            first = ', '.join(values[:-1])
            return f'{first}, or {values[-1]}'

        elif name == "Literal":
            values = ', '.join(format_annotation(t)
                               for t in annotation.__args__)
            return f'{{{values}}}'
        elif name == 'list':
            values = ', '.join(format_annotation(t)
                               for t in annotation.__args__)
            return f'list of {values}'

    return repr(annotation)


def add_types_to_docstring(docstring, annotations, defaults):

    indent_regex = r"^( +)Parameters\s*\n +[-=]{10}"
    indent_match = re.search(indent_regex, docstring, flags=re.MULTILINE)
    if not indent_match:
        return docstring
    n_indent = len(indent_match.group(1))

    indent = " " * n_indent
    param_regex = re.compile(f"{indent}(\\w+) :")
    lines = docstring.split('\n')

    for lineno, line in enumerate(lines):
        found_param = param_regex.match(line)
        if not found_param:
            continue
        name = found_param.group(1)

        if name not in annotations:
            continue

        annotation = annotations[name]
        type_str = format_annotation(annotation)
        new_line = f"{indent}{name} : {type_str}"
        if name in defaults:
            new_line += f", (default={defaults[name]})"
        lines[lineno] = new_line

    return "\n".join(lines)
