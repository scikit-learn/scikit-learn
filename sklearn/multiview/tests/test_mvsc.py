import numpy as np
from numpy.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_raises
import sklearn.multiview.mvsc as mvsc


def test_laplacian_ng():
    data = np.arange(25, dtype=float).reshape((5, 5))
    laplacian = mvsc.laplacian_ng(data)

    real_laplacian = np.array([[0., 0.0534, 0.0816, 0.1028, 0.1206],
                               [0.2672, 0.1714, 0.1527, 0.1466, 0.1450],
                               [0.4082, 0.2400, 0.2, 0.1820, 0.1723],
                               [0.5144, 0.2933, 0.2380, 0.2117, 0.1964],
                               [0.6030, 0.3384, 0.2708, 0.2378, 0.2181]])

    assert_array_almost_equal(laplacian, real_laplacian, decimal=4)


def test_suggested_sigma():
    data = np.arange(25, dtype=float).reshape((5, 5))
    s_sigma = mvsc.suggested_sigma(data)

    real_s_sigma = 7.0

    assert_array_almost_equal(s_sigma, real_s_sigma, decimal=4)


def test_gaussian_similarity():
    data = np.arange(25, dtype=float).reshape((5, 5))
    similarity = mvsc.distance_gaussian_similarity(data, 2)
    print(similarity)

    real_similarity = np.array([[1., 0.8824, 0.6065, 0.3246, 0.1353],
                                [0.0439, 0.0110, 0.0021, 0.0003, 0.],
                                [0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0.]])

    assert_array_almost_equal(similarity, real_similarity, decimal=4)


def test_mvsc_error():
    # Data and is_distane do not have the same length.
    one = np.arange(25, dtype=float).reshape((5, 5))
    two = np.arange(25, 50, dtype=float).reshape((5, 5))
    data = [one, two]
    is_distance = [False, False, False]

    mvsc_est = mvsc.MVSC(k=2)
    assert_raises(ValueError, mvsc_est.fit, data, is_distance)

    # Sample data matrices do not have the same number of rows
    one = np.arange(25, dtype=float).reshape((5, 5))
    two = np.arange(25, 49, dtype=float).reshape((4, 6))
    data = [one, two]
    is_distance = [False, False]

    mvsc_est = mvsc.MVSC(k=2)
    assert_raises(ValueError, mvsc_est.fit, data, is_distance)

    # k value cannot be negative
    one = np.arange(25, dtype=float).reshape((5, 5))
    two = np.arange(25, 50, dtype=float).reshape((5, 5))
    data = [one, two]
    is_distance = [False, False]

    mvsc_est = mvsc.MVSC(k=-2)
    assert_raises(ValueError, mvsc_est.fit, data, is_distance)


def test_mvsc():
    one = np.arange(25, dtype=float).reshape((5, 5))
    two = np.arange(25, 50, dtype=float).reshape((5, 5))
    data = [one, two]
    is_distance = [False, False]

    mvsc_est = mvsc.MVSC(k=2)
    mvsc_est.fit(data, is_distance)
