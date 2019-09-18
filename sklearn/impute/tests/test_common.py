import pytest

import numpy as np

from sklearn.utils.testing import assert_allclose
from sklearn.utils.testing import assert_array_equal

from sklearn.experimental import enable_iterative_imputer  # noqa

from sklearn.impute import IterativeImputer
from sklearn.impute import KNNImputer
from sklearn.impute import SimpleImputer


IMPUTERS = [IterativeImputer(), KNNImputer(), SimpleImputer()]


@pytest.mark.filterwarnings("ignore::sklearn.exceptions.ConvergenceWarning")
@pytest.mark.parametrize("imputer", IMPUTERS)
def test_imputation_missing_value_in_test_array(imputer):
    # [Non Regression Test for issue #13968] Missing value in test set should
    # not throw an error and return a finite dataset
    train = [[1], [2]]
    test = [[3], [np.nan]]
    imputer.set_params(add_indicator=True)
    imputer.fit(train).transform(test)


@pytest.mark.filterwarnings("ignore::sklearn.exceptions.ConvergenceWarning")
@pytest.mark.parametrize("marker", [np.nan, -1, 0])
@pytest.mark.parametrize("imputer", IMPUTERS)
def test_imputers_add_indicator(marker, imputer):
    X = np.array([
        [marker, 1,      5,      marker, 1],
        [2,      marker, 1,      marker, 2],
        [6,      3,      marker, marker, 3],
        [1,      2,      9,      marker, 4]
    ])
    X_true_indicator = np.array([
        [1., 0., 0., 1.],
        [0., 1., 0., 1.],
        [0., 0., 1., 1.],
        [0., 0., 0., 1.]
    ])
    imputer.set_params(missing_values=marker, add_indicator=True)

    X_trans = imputer.fit(X).transform(X)
    # The test is for testing the indicator,
    # that's why we're looking at the last 4 columns only.
    assert_allclose(X_trans[:, -4:], X_true_indicator)
    assert_array_equal(imputer.indicator_.features_, np.array([0, 1, 2, 3]))
