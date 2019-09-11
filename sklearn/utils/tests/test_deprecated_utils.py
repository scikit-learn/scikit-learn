import pytest

from sklearn.utils.estimator_checks import choose_check_classifiers_labels


# This file tests the utils that are deprecated


def test_choose_check_classifiers_labels_deprecated():
    with pytest.warns(DeprecationWarning, match="removed in version 0.24"):
        choose_check_classifiers_labels(None, None, None)
