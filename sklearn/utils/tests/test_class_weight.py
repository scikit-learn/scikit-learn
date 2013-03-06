import numpy as np

from sklearn.utils.class_weight import compute_class_weight

from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_true


def test_non_index_class_labels():
    """Test class_weight="auto" with unique(y) not the range [0, K)."""
    y = np.asarray([2, 2, 2, 3, 3, 4])
    cw = compute_class_weight("auto", np.asarray([2, 3, 4]), y)
    assert_almost_equal(cw.sum(), np.unique(y).shape)
    assert_true(cw[0] < cw[1] < cw[2])
