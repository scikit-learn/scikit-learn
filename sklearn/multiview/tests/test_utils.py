import numpy as np
from numpy.testing import assert_array_almost_equal
import sklearn.multiview.utils as utils


def test_hbeta():
    data = np.arange(25, dtype=float).reshape((5, 5))
    H, P = utils.Hbeta(data, 2)

    real_H = 100.145867478
    real_P = np.array([[8.64664717e-01, 1.17019644e-01, 1.58368867e-02,
                        2.14328955e-03, 2.90062698e-04],
                       [3.92557174e-05, 5.31268363e-06, 7.18993544e-07,
                        9.73051950e-08, 1.31688261e-08],
                       [1.78220681e-09, 2.41195464e-10, 3.26422564e-11,
                        4.41764902e-12, 5.97863781e-13],
                       [8.09120641e-14, 1.09502571e-14, 1.48195615e-15,
                        2.00560955e-16, 2.71429737e-17],
                       [3.67340203e-18, 4.97140904e-19, 6.72807051e-20,
                        9.10545327e-21, 1.23228910e-21]])

    assert_array_almost_equal(H, real_H, decimal=4)
    assert_array_almost_equal(P, real_P, decimal=4)


def test_x2p():
    data = np.arange(25, dtype=float).reshape((5, 5))
    P, beta = utils.x2p(data)

    real_P = np.array([[0., 0.25, 0.25, 0.25, 0.25],
                       [0.25, 0., 0.25, 0.25, 0.25],
                       [0.25, 0.25, 0., 0.25, 0.25],
                       [0.25, 0.25, 0.25, 0., 0.25],
                       [0.25, 0.25, 0.25, 0.25, 0.]])
    real_beta = np.array([8.88178420e-16, 8.88178420e-16, 8.88178420e-16,
                          8.88178420e-16, 8.88178420e-16])

    assert_array_almost_equal(P, real_P, decimal=4)
    assert_array_almost_equal(beta, real_beta, decimal=10)


def test_whiten():
    data = np.array([[1, 2, 3, 4], [4, 3, 2, 1], [2, 4, 1, 3], [1, 3, 2, 4]])
    whitened = utils.whiten(data)

    real_whitened = np.array([[9.63475981e-01, 1.11961253e+00, 1.49011612e-08,
                               0.00000000e+00],
                              [-1.55893688e+00, 6.91958598e-01, 0.00000000e+00,
                               0.00000000e+00],
                              [-1.84007539e-01, -1.46559183e+00,
                               -1.49011612e-08, 0.00000000e+00],
                              [7.79468442e-01, -3.45979299e-01, 0.00000000e+00,
                               0.00000000e+00]])

    assert_array_almost_equal(whitened, real_whitened, decimal=0)
