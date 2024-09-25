import numpy as np
import pytest

from scipy.conftest import array_api_compatible
from scipy._lib._array_api import (
    _GLOBAL_CONFIG, array_namespace, _asarray, copy, xp_assert_equal, is_numpy
)
import scipy._lib.array_api_compat.numpy as np_compat

skip_xp_backends = pytest.mark.skip_xp_backends


@pytest.mark.skipif(not _GLOBAL_CONFIG["SCIPY_ARRAY_API"],
        reason="Array API test; set environment variable SCIPY_ARRAY_API=1 to run it")
class TestArrayAPI:

    def test_array_namespace(self):
        x, y = np.array([0, 1, 2]), np.array([0, 1, 2])
        xp = array_namespace(x, y)
        assert 'array_api_compat.numpy' in xp.__name__

        _GLOBAL_CONFIG["SCIPY_ARRAY_API"] = False
        xp = array_namespace(x, y)
        assert 'array_api_compat.numpy' in xp.__name__
        _GLOBAL_CONFIG["SCIPY_ARRAY_API"] = True

    @array_api_compatible
    def test_asarray(self, xp):
        x, y = _asarray([0, 1, 2], xp=xp), _asarray(np.arange(3), xp=xp)
        ref = xp.asarray([0, 1, 2])
        xp_assert_equal(x, ref)
        xp_assert_equal(y, ref)

    @pytest.mark.filterwarnings("ignore: the matrix subclass")
    def test_raises(self):
        msg = "of type `numpy.ma.MaskedArray` are not supported"
        with pytest.raises(TypeError, match=msg):
            array_namespace(np.ma.array(1), np.array(1))

        msg = "of type `numpy.matrix` are not supported"
        with pytest.raises(TypeError, match=msg):
            array_namespace(np.array(1), np.matrix(1))

        msg = "only boolean and numerical dtypes are supported"
        with pytest.raises(TypeError, match=msg):
            array_namespace([object()])
        with pytest.raises(TypeError, match=msg):
            array_namespace('abc')

    def test_array_likes(self):
        # should be no exceptions
        array_namespace([0, 1, 2])
        array_namespace(1, 2, 3)
        array_namespace(1)

    @skip_xp_backends('jax.numpy',
                      reasons=["JAX arrays do not support item assignment"])
    @pytest.mark.usefixtures("skip_xp_backends")
    @array_api_compatible
    def test_copy(self, xp):
        for _xp in [xp, None]:
            x = xp.asarray([1, 2, 3])
            y = copy(x, xp=_xp)
            # with numpy we'd want to use np.shared_memory, but that's not specified
            # in the array-api
            x[0] = 10
            x[1] = 11
            x[2] = 12

            assert x[0] != y[0]
            assert x[1] != y[1]
            assert x[2] != y[2]
            assert id(x) != id(y)

    @array_api_compatible
    @pytest.mark.parametrize('dtype', ['int32', 'int64', 'float32', 'float64'])
    @pytest.mark.parametrize('shape', [(), (3,)])
    def test_strict_checks(self, xp, dtype, shape):
        # Check that `_strict_check` behaves as expected
        dtype = getattr(xp, dtype)
        x = xp.broadcast_to(xp.asarray(1, dtype=dtype), shape)
        x = x if shape else x[()]
        y = np_compat.asarray(1)[()]

        options = dict(check_namespace=True, check_dtype=False, check_shape=False)
        if xp == np:
            xp_assert_equal(x, y, **options)
        else:
            with pytest.raises(AssertionError, match="Namespaces do not match."):
                xp_assert_equal(x, y, **options)

        options = dict(check_namespace=False, check_dtype=True, check_shape=False)
        if y.dtype.name in str(x.dtype):
            xp_assert_equal(x, y, **options)
        else:
            with pytest.raises(AssertionError, match="dtypes do not match."):
                xp_assert_equal(x, y, **options)

        options = dict(check_namespace=False, check_dtype=False, check_shape=True)
        if x.shape == y.shape:
            xp_assert_equal(x, y, **options)
        else:
            with pytest.raises(AssertionError, match="Shapes do not match."):
                xp_assert_equal(x, y, **options)

    @array_api_compatible
    def test_check_scalar(self, xp):
        if not is_numpy(xp):
            pytest.skip("Scalars only exist in NumPy")

        if is_numpy(xp):
            with pytest.raises(AssertionError, match="Types do not match."):
                xp_assert_equal(xp.asarray(0.), xp.float64(0))
            xp_assert_equal(xp.float64(0), xp.asarray(0.))
