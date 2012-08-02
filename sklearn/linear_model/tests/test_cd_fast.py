import warnings
from sys import version_info
import numpy as np
from sklearn.linear_model import cd_fast
from sklearn.linear_model.cd_fast import enet_coordinate_descent
from sklearn.linear_model.coordinate_descent import ElasticNet
from sklearn.linear_model.base import center_data
from numpy.testing import assert_array_almost_equal, assert_almost_equal, assert_equal
from sklearn.utils.testing import assert_greater
from sklearn.datasets.samples_generator import make_regression


def test_iteration_set():
    """
    Test iter_set option for enet_coordinate_descent.

    This is validated against enet_coordinate_descent, run only on the reduced
    data, corresponding to the iter_set.
    """
    X, y = make_regression(n_samples=40, n_features=20, n_informative=5,
                    random_state=0)
    n_samples, n_features = X.shape
    rho = 0.80
    alpha = 10
    l1_reg = alpha * rho * n_samples
    l2_reg = alpha * (1.0 - rho) * n_samples
    X = np.asfortranarray(X)
    w = np.zeros(n_features)

    iter_set = np.arange(n_features / 2, (n_features / 2) + 5, dtype=np.int32)
    X_red = X[:, iter_set]
    X_red = np.asfortranarray(X_red)
    w_red = w[iter_set]

    result_red, _, _ = cd_fast.enet_coordinate_descent(w_red, l1_reg, l2_reg,
                    X_red, y, max_iter=1000, tol=1e-9, positive=False)

    result_iter_set, _, _ = cd_fast.enet_coordinate_descent(w, l1_reg, l2_reg,
                        X, y, max_iter=1000, tol=1e-9, positive=False, iter_set=iter_set)

    assert_array_almost_equal(result_red, result_iter_set[iter_set], 7)
    assert_array_almost_equal(w[~iter_set], result_iter_set[~iter_set], 7)


if __name__ == '__main__':
    import nose
    nose.runmodule()
