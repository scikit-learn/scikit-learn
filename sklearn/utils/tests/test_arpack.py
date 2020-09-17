import numpy as np

from sklearn.utils._arpack import _init_arpack_v0
from numpy.testing import assert_allclose


def test_init_arpack_v0():
    v0s = [_init_arpack_v0(1000, i) for i in range(100)]
    v0 = np.concatenate(v0s)
    assert np.all(v0 <= 1)
    assert np.all(v0 >= -1)

    assert_allclose(np.mean(v0), 0, atol=1e-2)
    assert_allclose(np.std(v0), 1 / np.sqrt(3), atol=1e-3)
