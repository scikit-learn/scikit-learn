# Should be in test_estimator_checks. Putting this here temporarily because I'm
# too lazy to make it work without pytest.

import warnings

import pytest
from sklearn.experimental import enable_hist_gradient_boosting  # noqa
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.exceptions import SkipTestWarning
from sklearn.utils.estimator_checks import check_estimator


def test_strict_mode():

    est = HistGradientBoostingClassifier(max_iter=10, min_samples_leaf=1)
    with pytest.warns(SkipTestWarning):
        # some_strict_check() was skipped
        check_estimator(est, strict_mode=False)

    with warnings.catch_warnings(record=True) as w:
        # some_strict_check() was ran
        check_estimator(est, strict_mode=True)
    assert not w
