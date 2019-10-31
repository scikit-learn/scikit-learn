import pytest

import numpy as np
from scipy import sparse

from sklearn.utils._testing import assert_allclose
from sklearn.utils._testing import assert_allclose_dense_sparse
from sklearn.utils._testing import assert_array_equal

from sklearn.experimental import enable_iterative_imputer  # noqa

from sklearn.impute import IterativeImputer
from sklearn.impute import KNNImputer
from sklearn.impute import SimpleImputer


IMPUTERS = [IterativeImputer(), KNNImputer(), SimpleImputer()]
SPARSE_IMPUTERS = [SimpleImputer()]


# ConvergenceWarning will be raised by the IterativeImputer
@pytest.mark.filterwarnings("ignore::sklearn.exceptions.ConvergenceWarning")
@pytest.mark.parametrize("imputer", IMPUTERS)
def test_imputation_missing_value_in_test_array(imputer):
    # [Non Regression Test for issue #13968] Missing value in test set should
    # not throw an error and return a finite dataset
    train = [[1], [2]]
    test = [[3], [np.nan]]
    imputer.set_params(add_indicator=True)
    imputer.fit(train).transform(test)


# ConvergenceWarning will be raised by the IterativeImputer
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

    X_trans = imputer.fit_transform(X)
    assert_allclose(X_trans[:, -4:], X_true_indicator)
    assert_array_equal(imputer.indicator_.features_, np.array([0, 1, 2, 3]))

    imputer.set_params(add_indicator=False)
    X_trans_no_indicator = imputer.fit_transform(X)
    assert_allclose(X_trans[:, :-4], X_trans_no_indicator)


# ConvergenceWarning will be raised by the IterativeImputer
@pytest.mark.filterwarnings("ignore::sklearn.exceptions.ConvergenceWarning")
@pytest.mark.parametrize("marker", [np.nan, -1])
@pytest.mark.parametrize("imputer", SPARSE_IMPUTERS)
def test_imputers_add_indicator_sparse(imputer, marker):
    X = sparse.csr_matrix([
        [marker, 1,      5,      marker, 1],
        [2,      marker, 1,      marker, 2],
        [6,      3,      marker, marker, 3],
        [1,      2,      9,      marker, 4]
    ])
    X_true_indicator = sparse.csr_matrix([
        [1., 0., 0., 1.],
        [0., 1., 0., 1.],
        [0., 0., 1., 1.],
        [0., 0., 0., 1.]
    ])
    imputer.set_params(missing_values=marker, add_indicator=True)

    X_trans = imputer.fit_transform(X)
    assert_allclose_dense_sparse(X_trans[:, -4:], X_true_indicator)
    assert_array_equal(imputer.indicator_.features_, np.array([0, 1, 2, 3]))

    imputer.set_params(add_indicator=False)
    X_trans_no_indicator = imputer.fit_transform(X)
    assert_allclose_dense_sparse(X_trans[:, :-4], X_trans_no_indicator)
