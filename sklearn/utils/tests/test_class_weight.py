import numpy as np

from sklearn.utils.class_weight import compute_class_weight
from sklearn.utils.fixes import unique

from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_true


def test_compute_class_weight():
    """Test (and demo) compute_class_weight."""
    classes, y = unique(np.asarray([2, 2, 2, 3, 3, 4]), return_inverse=True)
    cw = compute_class_weight("auto", classes, y)
    assert_almost_equal(cw.sum(), classes.shape)
    assert_true(cw[0] < cw[1] < cw[2])
