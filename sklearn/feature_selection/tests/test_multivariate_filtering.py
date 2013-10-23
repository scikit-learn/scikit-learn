from sklearn.utils.testing import (assert_array_equal, assert_raises)

import numpy as np

from sklearn.feature_selection import MinRedundancyMaxRelevance

X = np.array([[1, 3, 1],
              [3, 3, 3],
              [1, 3, 1],
              [1, 3, 3],
              [1, 3, 1]])

y = np.array([3, 1, 3, 1, 3])

def test_mMRM():
    """
    Test MinRedundancyMaxRelevance with default setting.
    """

    m = MinRedundancyMaxRelevance().fit(X, y)

    assert_array_equal([2, 0], m.mask)

    assert_array_equal(0.6730116670092563, m.score[0])

    assert_raises(ValueError, MinRedundancyMaxRelevance, rule='none')
