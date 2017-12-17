import numpy as np
from numpy.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_raises
import sklearn.multiview.mvmds as mvmds


def test_preprocess_mds():
    data = np.arange(25, dtype=float).reshape((5, 5))
    preprocessed_data = mvmds.preprocess_mvmds(data)

    sim = np.array([[40., 20., 0., -20., -40.],
                    [20., 10., 0., -10., -20.],
                    [0., 0., 0., 0., 0.],
                    [-20., -10., 0., 10., 20.],
                    [-40., -20., 0., 20., 40.]])

    assert_array_almost_equal(preprocessed_data, sim, decimal=3)


def test_mvmds_error():
    # Data and is_distane do not have the same length.
    one = np.arange(25, dtype=float).reshape((5, 5))
    two = np.arange(25, 50, dtype=float).reshape((5, 5))
    data = [one, two]
    is_distance = [False, False, False]

    mvmds_est = mvmds.MVMDS(k=2)
    assert_raises(ValueError, mvmds_est.fit, data, is_distance)

    # Sample data matrices do not have the same number of rows
    one = np.arange(25, dtype=float).reshape((5, 5))
    two = np.arange(25, 49, dtype=float).reshape((4, 6))
    data = [one, two]
    is_distance = [False, False]

    mvmds_est = mvmds.MVMDS(k=2)
    assert_raises(ValueError, mvmds_est.fit, data, is_distance)

    # k value cannot be negative
    one = np.arange(25, dtype=float).reshape((5, 5))
    two = np.arange(25, 50, dtype=float).reshape((5, 5))
    data = [one, two]
    is_distance = [False, False]

    mvmds_est = mvmds.MVMDS(k=-2)
    assert_raises(ValueError, mvmds_est.fit, data, is_distance)


def test_mvmds():
    one = np.arange(25, dtype=float).reshape((5, 5))
    two = np.arange(25, 50, dtype=float).reshape((5, 5))
    data = [one, two]
    is_distance = [False, False]

    mvmds_est = mvmds.MVMDS(k=2)
    mvmds_est.fit(data, is_distance)
