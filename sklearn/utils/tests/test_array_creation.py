import numpy as np
import pytest
from unittest.mock import MagicMock

from sklearn.utils.array_creation import empty_like
from sklearn.utils.array_creation import np_version
import sklearn.utils.array_creation


@pytest.mark.skipif(np_version < (1, 17),
                    reason="NEP18 not supported before 1.17")
def test_empty_like_nep18():
    class ArrayLike:
        __array_function__ = MagicMock(return_value=42)

    # if NEP18 is supported, empty_like should be forwarded to us
    array_like = ArrayLike()
    value = empty_like(array_like, dtype=np.float32, shape=(4, 2))
    assert value == 42


def test_empty_like():
    # Normaly arrays should just work with all versions of numpy
    X = np.arange(8)
    Y = empty_like(X.reshape((4, 2)))
    assert isinstance(Y, np.ndarray)
    assert Y.shape == (4, 2)


def test_empty_like_no_nep18():
    class NotAnArray:
        def __array__(self):
            return np.arange(8, dtype=np.float64).reshape((4, 2))
    try:
        # we trick this module into thinking it is working with an older
        # version to also test/cover this branch with newer versions of numpy
        real_np_version = sklearn.utils.array_creation.np_version
        sklearn.utils.array_creation.np_version = (1, 16)

        no_array = NotAnArray()
        empty_like(no_array, dtype=np.float32, shape=(4, 2))
        # for numpy < 1.17, we should give an error msg, if we provide shape
        # with a non-numpy array, and no dtype
        with pytest.raises(NotImplementedError):
            empty_like(no_array, shape=(4, 2))

        # we can pass a non-ndarray object, but without shape
        no_array = NotAnArray()
        an_array = empty_like(no_array, dtype=np.float32)
        assert an_array.shape == (4, 2)
        assert an_array.dtype == np.float32

        # but with a ndarray, we can pass with shape
        second_array = empty_like(an_array, dtype=np.float64, shape=(3, 5))
        assert second_array.shape == (3, 5)
        assert second_array.dtype == np.float64

        # and the dtype is optional for ndarrays
        second_array_same_type = empty_like(an_array, shape=(3, 5))
        assert second_array_same_type.shape == (3, 5)
        assert second_array_same_type.dtype == np.float32

        c_like_array = empty_like(an_array.T, shape=(3, 5))
        assert c_like_array.flags['C_CONTIGUOUS']

        fortran_like_array = empty_like(an_array.T, order='A', shape=(3, 5))
        assert fortran_like_array.flags['F_CONTIGUOUS']

        # unlike numpy, we don't implement order=K
        with pytest.raises(NotImplementedError):
            empty_like(an_array, order='K', shape=(4, 2))

    finally:
        sklearn.utils.array_creation.np_version = real_np_version
