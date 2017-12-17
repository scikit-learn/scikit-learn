import numpy as np
from sklearn.utils.testing import assert_raises
import sklearn.multiview.mvtsne as mvtsne


def test_mvtsne_error():
    # Data and is_distane do not have the same length.
    one = np.arange(25, dtype=float).reshape((5, 5))
    two = np.arange(25, 50, dtype=float).reshape((5, 5))
    data = [one, two]
    is_distance = [False, False, False]

    mvtsne_est = mvtsne.MvtSNE(k=2)
    assert_raises(ValueError, mvtsne_est.fit, data, is_distance)

    # Sample data matrices do not have the same number of rows
    one = np.arange(25, dtype=float).reshape((5, 5))
    two = np.arange(25, 49, dtype=float).reshape((4, 6))
    data = [one, two]
    is_distance = [False, False]

    mvtsne_est = mvtsne.MvtSNE(k=2)
    assert_raises(ValueError, mvtsne_est.fit, data, is_distance)

    # k value cannot be negative
    one = np.arange(25, dtype=float).reshape((5, 5))
    two = np.arange(25, 50, dtype=float).reshape((5, 5))
    data = [one, two]
    is_distance = [False, False]

    mvtsne_est = mvtsne.MvtSNE(k=-2)
    assert_raises(ValueError, mvtsne_est.fit, data, is_distance)


def test_mvtsne():
    one = np.arange(25, dtype=float).reshape((5, 5))
    two = np.arange(25, 50, dtype=float).reshape((5, 5))
    data = [one, two]
    is_distance = [False, False]

    mvtsne_est = mvtsne.MvtSNE(k=2)
    mvtsne_est.fit(data, is_distance)
