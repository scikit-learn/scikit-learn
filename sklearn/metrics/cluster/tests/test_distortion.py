import numpy as np

from sklearn.utils.testing import assert_almost_equal

from sklearn.metrics.cluster.distortion import distortion


def test_distortion():
    X = np.array([[0, 0], [2, 2],
                  [5, 5], [6, 6]])
    labels = [0, 0, 1, 1]
    assert_almost_equal(distortion(X, labels), 2.5)
