from typing import Dict
from typing import List
from typing import Union
from typing import Callable
from typing import Optional

import pytest
import numpy as np

from sklearn.utils._typing import RandomStateType
from sklearn.utils._typing import EstimatorType
from sklearn.utils._typing import Literal
from sklearn.utils._typing import get_annotation_class_name
from sklearn.utils._typing import format_docstring_annotation
from sklearn.utils._typing import get_docstring_annotations


@pytest.mark.parametrize("annotation, expected_class", [
    (None, 'None'),
    (str, 'str'),
    (int, 'int'),
    (float, 'float'),
    (list, 'list'),
    (EstimatorType, 'TypeVar'),
    (List[int], 'List'),
    (Union[int, float], 'Union'),
    (Dict, 'Dict'),
    (Callable, 'Callable'),
    (Callable[[str], str], 'Callable'),
])
def test_get_annotation_class_name(annotation, expected_class):
    """Test annotation names are returned correct."""
    assert get_annotation_class_name(annotation) == expected_class


@pytest.mark.parametrize("annotation, expected_str", [
    (None, "None"),
    (EstimatorType, "estimator instance"),
    (np.random.RandomState, "RandomState instance"),
    (int, "int"),
    (float, 'float'),
    (list, "list"),
    (str, "str"),
    (List[int], "list of int"),
    (Optional[List[int]], "list of int or None"),
    (List[EstimatorType], "list of estimator instance"),
    (Optional[EstimatorType], "estimator instance or None"),
    (Union[int, float], "int or float"),
    (RandomStateType, "int, RandomState instance or None")
])
def test_format_docstring_annotation(annotation, expected_str):
    """Check format for auto generation annotations."""
    assert format_docstring_annotation(annotation) == expected_str


class _TypingObject:
    def __init__(self,
                 estimator: EstimatorType,
                 num: int = 10, union_num: Union[int, float] = 1.4,
                 float_num: float = 1e-4,
                 pet: Literal['dog'] = 'dog',
                 weather: Literal['sunny', 'cloudy'] = 'sunny',
                 random_state: RandomStateType = None):
        pass


def test_get_docstring_annotations():
    """Check docstring for annotations."""
    pytest.importorskip("typing_extensions")
    annotations = get_docstring_annotations(_TypingObject.__init__)

    assert annotations['estimator'] == "estimator instance"
    assert annotations['num'] == "int, default=10"
    assert annotations['float_num'] == "float, default="
    assert annotations['union_num'] == "int or float, default="
    assert annotations['pet'] == "'dog', default='dog'"
    assert annotations['weather'] == "{'sunny', 'cloudy'}, default='sunny'"
    assert annotations['random_state'] == ("int, RandomState instance or None"
                                           ", default=None")
