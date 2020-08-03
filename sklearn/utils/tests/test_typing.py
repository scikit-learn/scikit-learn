from typing import Dict
from typing import Any
from typing import List
from typing import Union
from typing import Callable

import pytest

from sklearn.base import BaseEstimator
from sklearn.utils._typing import Literal
from sklearn.utils._typing import _get_annotation_class_name
# from sklearn.utils._typing import _format_annotation


@pytest.mark.parametrize("annotation, expected_class", [
    (None, 'None'),
    (Any, 'Any'),
    (str, 'str'),
    (int, 'int'),
    (float, 'float'),
    (list, 'list'),
    (BaseEstimator, 'BaseEstimator'),
    (List[int], 'List'),
    (Union[int, float], 'Union'),
    (Dict, 'Dict'),
    (Literal['a', 'b'], 'Literal'),
    (Callable, 'Callable'),
    (Callable[[str], str], 'Callable'),
])
def test_get_annotation_class_name(annotation, expected_class):
    assert _get_annotation_class_name(annotation) == expected_class


# @pytest.mark.parametrize("annotation", "expected_str", [
#     (None, 'None')
# ])
# def test_format_annotation(annotation, expected_str):
#     pass
