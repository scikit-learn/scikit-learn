import pytest

from scipy.special._support_alternative_backends import (get_array_special_func,
                                                         array_special_func_map)
from scipy.conftest import array_api_compatible
from scipy import special
from scipy._lib._array_api_no_0d import xp_assert_close
from scipy._lib._array_api import is_jax, is_torch, SCIPY_DEVICE
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


@pytest.mark.fail_slow(5)
@array_api_compatible
# @pytest.mark.skip_xp_backends('numpy', reason='skip while debugging')
# @pytest.mark.usefixtures("skip_xp_backends")
# `reversed` is for developer convenience: test new function first = less waiting
@pytest.mark.parametrize('f_name_n_args', reversed(array_special_func_map.items()))
@pytest.mark.parametrize('dtype', ['float32', 'float64'])
@pytest.mark.parametrize('shapes', [[(0,)]*4, [tuple()]*4, [(10,)]*4,
                                    [(10,), (11, 1), (12, 1, 1), (13, 1, 1, 1)]])
def test_support_alternative_backends(xp, f_name_n_args, dtype, shapes):
    f_name, n_args = f_name_n_args

    if (SCIPY_DEVICE != 'cpu'
        and is_torch(xp)
        and f_name in {'stdtr', 'betaincc', 'betainc'}
    ):
        pytest.skip(f"`{f_name}` does not have an array-agnostic implementation "
                    f"and cannot delegate to PyTorch.")

    shapes = shapes[:n_args]
    f = getattr(special, f_name)

    dtype_np = getattr(np, dtype)
    dtype_xp = getattr(xp, dtype)

    # # To test the robustness of the alternative backend's implementation,
    # # use Hypothesis to generate arguments
    # from hypothesis import given, strategies, reproduce_failure, assume
    # import hypothesis.extra.numpy as npst
    # @given(data=strategies.data())
    # mbs = npst.mutually_broadcastable_shapes(num_shapes=n_args)
    # shapes, final_shape = data.draw(mbs)
    # elements = dict(allow_subnormal=False)  # consider min_value, max_value
    # args_np = [np.asarray(data.draw(npst.arrays(dtype_np, shape, elements=elements)),
    #                       dtype=dtype_np)
    #            for shape in shapes]

    # For CI, be a little more forgiving; just generate normally distributed arguments
    rng = np.random.default_rng(984254252920492019)
    args_np = [rng.standard_normal(size=shape, dtype=dtype_np) for shape in shapes]

    if (is_jax(xp) and f_name == 'gammaincc'  # google/jax#20699
            or f_name == 'chdtrc'):  # gh-20972
        args_np[0] = np.abs(args_np[0])
        args_np[1] = np.abs(args_np[1])

    args_xp = [xp.asarray(arg[()], dtype=dtype_xp) for arg in args_np]

    res = f(*args_xp)
    ref = xp.asarray(f(*args_np), dtype=dtype_xp)

    eps = np.finfo(dtype_np).eps
    xp_assert_close(res, ref, atol=10*eps)


@array_api_compatible
def test_chdtr_gh21311(xp):
    # the edge case behavior of generic chdtr was not right; see gh-21311
    # be sure to test at least these cases
    # should add `np.nan` into the mix when gh-21317 is resolved
    x = np.asarray([-np.inf, -1., 0., 1., np.inf])
    v = x.reshape(-1, 1)
    ref = special.chdtr(v, x)
    res = special.chdtr(xp.asarray(v), xp.asarray(x))
    xp_assert_close(res, xp.asarray(ref))
