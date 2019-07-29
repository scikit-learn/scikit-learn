import pytest
import numpy as np

from sklearn.utils.testing import assert_array_equal
from sklearn.utils import NamedArray


def test_basics():
    x = NamedArray(np.random.rand(5, 3), feature_names=['a', 'b', 'c'])
    assert_array_equal(x.feature_names, ['a', 'b', 'c'])
    assert not isinstance(x + 1, NamedArray)
    assert not isinstance(x + x, NamedArray)
    assert not isinstance(x + np.ones(shape=(5, 3)), NamedArray)
    assert x[0, :].shape == (3,)
    assert x[:, 0].shape == (5,)
    assert x[0:2].shape == (2, 3)


def test_validation():
    with pytest.raises(ValueError, match="column names provided"):
        NamedArray(np.ones(shape=(3, 3)), feature_names=[1])

    # allow None as feature_names
    NamedArray(np.ones(shape=(3, 3)))


def test_getattr():
    x = NamedArray(np.random.rand(5, 3), feature_names=['a', 'b', 'c'])
    # these would fail if __getattr__ doesn't work
    x.ndim
    x.shape


def test_pandas():
    _ = pytest.importorskip("pandas")
    x = NamedArray(np.random.rand(5, 3), feature_names=['a', 'b', 'c'])
    assert all(x.todataframe().columns == ['a', 'b', 'c'])


def test_1d():
    x = NamedArray(np.array([1, 2]), feature_names=['a'])
    assert x.feature_names == ['a']


def test_repr():
    repr_ = ("NamedArray([[1 2]\n"
             "            [3 4]],\n"
             "           feature_names=['a' 'b'])")
    x = NamedArray([[1, 2], [3, 4]], feature_names=['a', 'b'])
    assert repr(x) == repr_
