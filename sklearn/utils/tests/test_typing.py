from typing import Dict
from typing import Any
from typing import List
from typing import Union
from typing import Callable
from typing import Optional
from typing_extensions import Literal

import pytest
import numpy as np

from sklearn.base import BaseEstimator
from sklearn.utils._typing import RandomState
from sklearn.utils._typing import _get_annotation_class_name
from sklearn.utils._typing import _format_docstring_annotation
from sklearn.utils._typing import get_docstring_annotations


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


@pytest.mark.parametrize("annotation, expected_str", [
    (None, "None"),
    (BaseEstimator, "estimator instance"),
    (np.random.RandomState, "RandomState instance"),
    (int, "int"),
    (float, 'float'),
    (list, "list"),
    (str, "str"),
    (List[int], "list of int"),
    (Optional[List[int]], "list of int or None"),
    (List[BaseEstimator], "list of estimator instance"),
    (Optional[BaseEstimator], "estimator instance or None"),
    (Union[int, float], "int or float"),
    (Literal['cat', 'dog'], "{'cat', 'dog'}"),
    (Literal['cat'], "'cat'"),
    (RandomState, "int, RandomState instance or None")
])
def test_format_docstring_annotation(annotation, expected_str):
    assert _format_docstring_annotation(annotation) == expected_str


class TestObject:
    def __init__(self,
                 estimator: BaseEstimator,
                 num: int = 10, union_num: Union[int, float] = 1.4,
                 float_num: float = 1e-4,
                 pet: Literal['cat', 'dog'] = 'dog',
                 random_state: RandomState = None):
        pass


def test_get_docstring_annotations():
    annotations = get_docstring_annotations(TestObject.__init__)
    assert annotations['estimator'] == "estimator instance"
    assert annotations['num'] == "int, default=10"
    assert annotations['float_num'] == "float, default="
    assert annotations['union_num'] == "int or float, default="
    assert annotations['pet'] == "{'cat', 'dog'}, default='dog'"
    assert annotations['random_state'] == ("int, RandomState instance or None"
                                           ", default=None")
