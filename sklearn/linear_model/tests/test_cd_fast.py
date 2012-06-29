import warnings
from sys import version_info
import numpy as np
from sklearn.linear_model import cd_fast
from sklearn.linear_model.cd_fast import enet_coordinate_descent
from sklearn.linear_model.cd_fast import enet_coordinate_descent_old
from sklearn.linear_model.coordinate_descent import ElasticNet
from sklearn.linear_model.base import center_data
from numpy.testing import assert_array_almost_equal, assert_almost_equal, \
                          assert_equal
from sklearn.utils.testing import assert_greater
from sklearn.datasets.samples_generator import make_regression


def check_warnings():
    if version_info < (2, 6):
        raise SkipTest("Testing for warnings is not supported in versions \
        older than Python 2.6")

# ATTENTION does not pass with w = 0 as start value
def test_line():

    X = np.array([[-1], [0.], [1.]])
    y = np.array([-1.0, 0.0, 1.0]) # just a straight line
    n_samples, n_features = X.shape
    rho = 0.3
    alpha = 0.5
    alpha = alpha * rho * n_samples
    beta = alpha * (1.0 - rho) * n_samples
    w = np.array([0.2])
    old_result = cd_fast.enet_coordinate_descent(w, alpha, beta, X, y,
                                max_iter=100, tol=1e-4, positive=False)[0]

    assert_array_almost_equal(old_result,
                np.array([0.52631579]))
# cd_fast.enet_coordinate_descent(w, alpha, beta,
# X, y, max_iter=100, tol=1e-4, positive=False)[0])


def test_2d():
    X = np.array([[-1, 0.0], [0., 1.0], [1., -1.]])
    y = np.array([-1.0, 0.0, 1.0]) # just a straight line

    rho = 0.3
    alpha = 0.5

    n_samples, n_features = X.shape
    l2_reg = alpha * rho * n_samples
    l1_reg = alpha * (1.0 - rho) * n_samples
    w = np.zeros(n_features)
    X = np.asfortranarray(X)
    old_result, old_gap, old_tol = enet_coordinate_descent_old(w, l2_reg, l1_reg,
                                X, y, max_iter=100, tol=1e-7, positive=False)
    w = np.zeros(n_features)
    #print result_org
    my_result, gab, tol = enet_coordinate_descent(w, l2_reg, l1_reg,
                                X, y, max_iter=100, tol=1e-7, positive=False)
    assert_array_almost_equal(my_result, old_result, 9)
    # assert_array_almost_equal(my_result,
    # np.array([0.52323384, -0.00908868]),7)


def test_active_set():
    # test set with zeros in solution coef
    X, y = make_regression(n_samples=40, n_features=20, n_informative=5,
                    random_state=0)
    n_samples, n_features = X.shape
    rho = 0.80
    alpha = 10
    alpha = alpha * rho * n_samples
    beta = alpha * (1.0 - rho) * n_samples
    w = np.zeros(n_features)
    X = np.asfortranarray(X)
    old_result, old_gap, old_tol = cd_fast.enet_coordinate_descent_old(w, alpha, beta,
                                X, y, max_iter=1000, tol=1e-9, positive=False)
    w = np.zeros(n_features)
    result, gap, tol = cd_fast.enet_coordinate_descent(w, alpha, beta,
                                X, y, max_iter=1000, tol=1e-9, positive=False)

    assert_array_almost_equal(result, old_result, 7)
    assert_array_almost_equal(gap, old_gap, 7)

#    np.array([38.18037338, 18.4702112, 9.86198851, -1.46801215, 16.52490931
#     , 14.26861543, 18.15508878, 36.40871624, 0., 12.35964046
#     ,6.98213445, 30.17242224,7.07032768,4.42177579, -1.73831861
#     , -7.26278943,0.34912212, 48.84641316,8.05922053, 10.301779])


def test_not_enough_memory():
    X, y = make_regression(n_samples=40, n_features=20, n_informative=5,
                    random_state=0)
    n_samples, n_features = X.shape
    rho = 0.80
    alpha = 10
    alpha = alpha * rho * n_samples
    beta = alpha * (1.0 - rho) * n_samples
    w = np.zeros(n_features)
    X = np.asfortranarray(X)
    old_result, old_gap, old_tol = cd_fast.enet_coordinate_descent_old(w, alpha, beta,
                                X, y, max_iter=1000, tol=1e-9, positive=False)
    w = np.zeros(n_features)

    with warnings.catch_warnings():
        # Here we have a small number of iterations, and thus the
        # ElasticNet might not converge. This is to speed up tests
        warnings.simplefilter("ignore", UserWarning)

        result, gap, tol = cd_fast.enet_coordinate_descent(w, alpha, beta,
                        X, y, max_iter=1000, tol=1e-9, positive=False, memory_limit=10)

    assert_array_almost_equal(result, old_result, 7)
    assert_array_almost_equal(gap, old_gap, 7)


def test_memory_limit_warning():
    check_warnings()  # Skip if unsupported Python version
    with warnings.catch_warnings(record=True) as warning:
        warnings.simplefilter('always')

            # test set with zeros in solution coef
        X, y = make_regression(n_samples=40, n_features=20, n_informative=5,
                        random_state=0)
        n_samples, n_features = X.shape
        rho = 0.80
        alpha = 10
        alpha = alpha * rho * n_samples
        beta = alpha * (1.0 - rho) * n_samples
        w = np.zeros(n_features)
        X = np.asfortranarray(X)

        w = np.zeros(n_features)
        cd_fast.enet_coordinate_descent(w, alpha, beta,
                                    X, y, max_iter=1000, tol=1e-9, positive=False, 
                                    memory_limit=10)
        assert_greater(len(warning), 0)  # warnings should be raised


if __name__ == '__main__':
    import nose
    nose.runmodule()
