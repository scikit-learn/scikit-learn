import pytest
from hypothesis import given, strategies, reproduce_failure  # noqa: F401
import hypothesis.extra.numpy as npst

from scipy.special._support_alternative_backends import (get_array_special_func,
                                                         array_special_func_map)
from scipy.conftest import array_api_compatible
from scipy import special
from scipy._lib._array_api import xp_assert_close, is_jax
from scipy._lib.array_api_compat import numpy as np

try:
    import array_api_strict
    HAVE_ARRAY_API_STRICT = True
except ImportError:
    HAVE_ARRAY_API_STRICT = False


@pytest.mark.skipif(not HAVE_ARRAY_API_STRICT,
                    reason="`array_api_strict` not installed")
def test_dispatch_to_unrecognize_library():
    xp = array_api_strict
    f = get_array_special_func('ndtr', xp=xp, n_array_args=1)
    x = [1, 2, 3]
    res = f(xp.asarray(x))
    ref = xp.asarray(special.ndtr(np.asarray(x)))
    xp_assert_close(res, ref, xp=xp)


@pytest.mark.parametrize('dtype', ['float32', 'float64', 'int64'])
@pytest.mark.skipif(not HAVE_ARRAY_API_STRICT,
                    reason="`array_api_strict` not installed")
def test_rel_entr_generic(dtype):
    xp = array_api_strict
    f = get_array_special_func('rel_entr', xp=xp, n_array_args=2)
    dtype_np = getattr(np, dtype)
    dtype_xp = getattr(xp, dtype)
    x, y = [-1, 0, 0, 1], [1, 0, 2, 3]

    x_xp, y_xp = xp.asarray(x, dtype=dtype_xp), xp.asarray(y, dtype=dtype_xp)
    res = f(x_xp, y_xp)

    x_np, y_np = np.asarray(x, dtype=dtype_np), np.asarray(y, dtype=dtype_np)
    ref = special.rel_entr(x_np[-1], y_np[-1])
    ref = np.asarray([np.inf, 0, 0, ref], dtype=ref.dtype)

    xp_assert_close(res, xp.asarray(ref), xp=xp)


@pytest.mark.fail_slow(2)
@array_api_compatible
@given(data=strategies.data())
@pytest.mark.parametrize('f_name_n_args', array_special_func_map.items())
def test_support_alternative_backends(xp, data, f_name_n_args):
    f_name, n_args = f_name_n_args

    if is_jax(xp):
        if f_name in ['gammainc', 'gammaincc']:
            pytest.skip("google/jax#20507")
        if f_name == 'rel_entr':
            pytest.skip("google/jax#21265")

    f = getattr(special, f_name)

    mbs = npst.mutually_broadcastable_shapes(num_shapes=n_args)
    shapes, final_shape = data.draw(mbs)

    dtype = data.draw(strategies.sampled_from(['float32', 'float64']))
    dtype_np = getattr(np, dtype)
    dtype_xp = getattr(xp, dtype)

    elements = dict(min_value=dtype_np(-10), max_value=dtype_np(10),
                    allow_subnormal=False)
    args_np = [np.asarray(data.draw(npst.arrays(dtype_np, shape, elements=elements)))
               for shape in shapes]

    # `torch.asarray(np.asarray(1.))` produces
    # TypeError: can't convert np.ndarray of type numpy.object_.
    # So we extract the scalar from 0d arrays.
    args_xp = [xp.asarray(arg[()], dtype=dtype_xp) for arg in args_np]

    ref = np.asarray(f(*args_np))
    res = f(*args_xp)

    eps = np.finfo(dtype).eps
    # PyTorch seems to struggle with precision near the poles of `gammaln`,
    # so the tolerance needs to be quite loose (eps**0.2) - see gh-19935.
    # To compensate, we also check that the root-mean-square error is
    # less than eps**0.5.
    ref = xp.asarray(ref, dtype=dtype_xp)
    xp_assert_close(res, ref, rtol=eps**0.2, atol=eps*20,
                    check_namespace=True, check_shape=True, check_dtype=True,)
    xp_assert_close(xp.sqrt(xp.mean(res**2)), xp.sqrt(xp.mean(ref**2)),
                    rtol=eps**0.5, atol=eps*20,
                    check_namespace=False, check_shape=False, check_dtype=False,)
