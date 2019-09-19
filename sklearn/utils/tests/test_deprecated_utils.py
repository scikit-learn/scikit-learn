import pytest
import numpy as np

from sklearn.dummy import DummyClassifier
from sklearn.utils.estimator_checks import choose_check_classifiers_labels
from sklearn.utils.estimator_checks import NotAnArray
from sklearn.utils.estimator_checks import enforce_estimator_tags_y
from sklearn.utils.estimator_checks import is_public_parameter
from sklearn.utils.estimator_checks import pairwise_estimator_convert_X
from sklearn.utils.estimator_checks import set_checking_parameters


# This file tests the utils that are deprecated


def test_choose_check_classifiers_labels_deprecated():
    with pytest.warns(DeprecationWarning, match="removed in version 0.24"):
        choose_check_classifiers_labels(None, None, None)


def test_enforce_estimator_tags_y():
    with pytest.warns(DeprecationWarning, match="removed in version 0.24"):
        enforce_estimator_tags_y(DummyClassifier(), np.array([0, 1]))


def test_notanarray():
    with pytest.warns(DeprecationWarning, match="removed in version 0.24"):
        NotAnArray([1, 2])


def test_is_public_parameter():
    with pytest.warns(DeprecationWarning, match="removed in version 0.24"):
        is_public_parameter('hello')


def test_pairwise_estimator_convert_X():
    with pytest.warns(DeprecationWarning, match="removed in version 0.24"):
        pairwise_estimator_convert_X([[1, 2]], DummyClassifier())


def test_set_checking_parameters():
    with pytest.warns(DeprecationWarning, match="removed in version 0.24"):
        set_checking_parameters(DummyClassifier())
