import pytest
import numpy as np

from sklearn.dummy import DummyClassifier
from sklearn.utils.estimator_checks import choose_check_classifiers_labels
from sklearn.utils.estimator_checks import enforce_estimator_tags_y


# This file tests the utils that are deprecated


def test_choose_check_classifiers_labels_deprecated():
    with pytest.warns(DeprecationWarning, match="removed in version 0.24"):
        choose_check_classifiers_labels(None, None, None)


def test_enforce_estimator_tags_y():
    with pytest.warns(DeprecationWarning, match="removed in version 0.24"):
        enforce_estimator_tags_y(DummyClassifier(), np.array([0, 1]))
