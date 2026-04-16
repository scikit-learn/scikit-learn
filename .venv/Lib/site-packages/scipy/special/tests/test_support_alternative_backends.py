from functools import partial
import pickle

import pytest
from hypothesis import given, strategies
import hypothesis.extra.numpy as npst
from packaging import version

from scipy import special
from scipy.special._support_alternative_backends import _special_funcs
from scipy._lib._array_api_no_0d import xp_assert_close
from scipy._lib._array_api import (is_cupy, is_dask, is_jax, is_torch,
                                   make_xp_pytest_param, make_xp_test_case,
                                   get_native_namespace_name)
from scipy._lib.array_api_compat import numpy as np
import scipy._lib.array_api_extra as xpx


# Run all tests in this module in the Array API CI, including those without
# the xp fixture
pytestmark = pytest.mark.array_api_backends

lazy_xp_modules = [special]


def _skip_or_tweak_alternative_backends(xp, nfo, dtypes, int_only):
    """Skip tests for specific intersections of scipy.special functions 
    vs. backends vs. dtypes vs. devices.
    Also suggest bespoke tweaks.

    Returns
    -------
    positive_only : list[bool]
        Whether you should exclusively test positive inputs.
    dtypes : list[str]
        dtype strings 'float64', 'int32', 'int64', etc. with integer types
        mapped to the type of the NumPy default int.
    """
    f_name = nfo.name
    if isinstance(nfo.positive_only, dict):
        positive_only = nfo.positive_only.get(get_native_namespace_name(xp), False)
    else:
        positive_only = nfo.positive_only
    if isinstance(positive_only, bool):
        positive_only = [positive_only]*nfo.n_args

    dtypes = [np.intp.__name__ if dtype == "intp" else dtype for dtype in dtypes]

    if f_name in {'betaincinv'} and is_cupy(xp):
        pytest.xfail("CuPy uses different convention for out of domain input.")

    if (
            f_name == "sinc" and "float32" in dtypes
            and version.parse(np.__version__) < version.parse("2")
    ):
        pytest.xfail("https://github.com/numpy/numpy/issues/11204")

    if not any('int' in dtype for dtype in dtypes):
        return positive_only, dtypes

    # Integer-specific issues from this point onwards

    if f_name in {'gamma', 'gammasgn'} and is_cupy(xp):
        # CuPy has not yet updated gamma pole behavior to match
        # https://github.com/scipy/scipy/pull/21827.
        positive_only = [True]

    if f_name in {'poch'} and is_jax(xp):
        # Jax uses a different convention at gamma poles.
        positive_only = [True, True]

    if f_name == 'multigammaln':
        pytest.skip("multigammaln raises for out of domain inputs.")

    if ((is_torch(xp) and f_name in {'gammainc', 'gammaincc'})
        or (is_cupy(xp) and f_name in {'stdtr', 'i0e', 'i1e'})
        or (is_jax(xp) and f_name in {'stdtr', 'ndtr', 'ndtri', 'log_ndtr', 'hyp1f1',
                                      'hyp2f1', 'spence', 'kl_div'})
    ):
        pytest.skip(f"`{f_name}` does not support integer types")

    # int/float mismatched args support is sketchy
    if (any('float' in dtype for dtype in dtypes)
        and ((is_torch(xp) and f_name in ('rel_entr', 'xlogy', 'polygamma',
                                          'zeta', 'xlog1py'))
             or (is_jax(xp) and f_name in ('gammainc', 'gammaincc', 'expn',
                                           'rel_entr', 'xlogy', 'betaln',
                                           'polygamma', 'zeta', 'poch',
                                           'xlog1py')))
    ):
        pytest.xfail("dtypes do not match")

    if (is_torch(xp) and xpx.default_dtype(xp) == xp.float32):
        # On PyTorch with float32 default dtype, all ints are promoted
        # to float32, but when falling back to NumPy/SciPy int64 is promoted
        # instead to float64. Integer only parameters essentially do not
        # participate in determination of the result type in PyTorch with
        # float32 default dtype, but will impact the output dtype as if
        # they were float64 when falling back to NumPy/SciPy.
        if not nfo.torch_native:
            pytest.xfail("dtypes do not match")

    return positive_only, dtypes


@pytest.mark.filterwarnings("ignore:invalid value encountered:RuntimeWarning:dask")
@pytest.mark.filterwarnings("ignore:overflow encountered:RuntimeWarning")
@pytest.mark.parametrize('shapes', [[(0,)]*4, [tuple()]*4, [(10,)]*4,
                                    [(10,), (11, 1), (12, 1, 1), (13, 1, 1, 1)]])
@pytest.mark.parametrize('base_dtype', ['float32', 'float64', 'intp'])
@pytest.mark.parametrize(
    'func,nfo', [make_xp_pytest_param(i.wrapper, i) for i in _special_funcs])
def test_support_alternative_backends(xp, func, nfo, base_dtype, shapes):
    int_only = nfo.int_only
    if int_only is None:
        int_only = (False, ) * nfo.n_args
        dtypes = (base_dtype, ) * nfo.n_args
    else:
        dtypes = tuple(
            'intp' if needs_int else base_dtype for needs_int in int_only
        )

    positive_only, dtypes = _skip_or_tweak_alternative_backends(
        xp, nfo, dtypes, int_only
    )

    dtypes_np = [getattr(np, dtype) for dtype in dtypes]
    dtypes_xp = [getattr(xp, dtype) for dtype in dtypes]

    shapes = shapes[:nfo.n_args]
    rng = np.random.default_rng(984254252920492019)
    args_np = []

    # Handle cases where there's an argument which only takes scalar values.
    python_int_only = nfo.python_int_only
    if isinstance(python_int_only, dict):
        python_int_only = python_int_only.get(get_native_namespace_name(xp))
    scalar_or_0d_only = nfo.scalar_or_0d_only
    if isinstance(scalar_or_0d_only, dict):
        scalar_or_0d_only = scalar_or_0d_only.get(get_native_namespace_name(xp))

    test_large_ints = nfo.test_large_ints
    if isinstance(nfo.test_large_ints, dict):
        test_large_ints = test_large_ints.get(get_native_namespace_name(xp), False)

    if python_int_only is None:
        python_int_only = [False] * nfo.n_args
    if scalar_or_0d_only is None:
        scalar_or_0d_only = [False] * nfo.n_args

    no_shape = [
        cond1 or cond2 for cond1, cond2 in zip(python_int_only, scalar_or_0d_only)
    ]

    shapes = [shape if not cond else None for shape, cond in zip(shapes, no_shape)]

    for dtype, dtype_np, shape, needs_python_int in zip(
            dtypes, dtypes_np, shapes, python_int_only
    ):
        if 'int' in dtype and test_large_ints:
            iinfo = np.iinfo(dtype_np)
            rand = partial(rng.integers, iinfo.min, iinfo.max + 1)
        elif 'int' in dtype:
            rand = partial(rng.integers, -20, 21)
        else:
            rand = rng.standard_normal
        val = rand(size=shape, dtype=dtype_np)
        if needs_python_int:
            # The logic above for determining shapes guarantees that
            # shape will be None in the above line when a Python int is required,
            # so this can safely be converted to an int.
            val = int(val)
        args_np.append(val)

    args_np = [
        np.abs(arg) if cond else arg for arg, cond in zip(args_np, positive_only)
    ]

    args_xp = [
        xp.asarray(arg, dtype=dtype_xp) if not needs_python_int
        else arg
        for arg, dtype_xp, needs_python_int
        in zip(args_np, dtypes_xp, python_int_only)
    ]

    args_np = [
        np.asarray(arg, dtype=dtype_np) if not needs_python_int
        else arg
        for arg, dtype_np, needs_python_int
        in zip(args_np, dtypes_np, python_int_only)
    ]

    if is_dask(xp):
        # We're using map_blocks to dispatch the function to Dask.
        # This is the correct thing to do IF all tested functions are elementwise;
        # otherwise the output would change depending on chunking.
        # Try to trigger bugs related to having multiple chunks.
        args_xp = [arg.rechunk(5) for arg in args_xp]

    res = nfo.wrapper(*args_xp)  # Also wrapped by lazy_xp_function
    ref = nfo.func(*args_np)  # Unwrapped ufunc
    if (
            is_torch(xp)
            and xpx.default_dtype(xp) == xp.float32
            and "float64" not in dtypes
    ):
        # int64 promotes like float32 on torch with default dtype = float32
        # cast reference if needed
        ref = np.float32(ref)
    # When dtype_np is integer, the output dtype can be float
    atol = 0 if ref.dtype.kind in 'iu' else 10 * np.finfo(ref.dtype).eps
    rtol = None
    if is_torch(xp) and func.__name__ == 'j1':
        # If we end up needing more function/backend specific tolerance
        # adjustments, this should be factored out properly.
        atol = 1e-7
        rtol = 1e-5
    xp_assert_close(
        res, xp.asarray(ref), rtol=rtol, atol=atol, check_0d=nfo.produces_0d
    )


@pytest.mark.parametrize(
    'func, nfo',
    [make_xp_pytest_param(i.wrapper, i) for i in _special_funcs if i.n_args >= 2])
@pytest.mark.filterwarnings("ignore:invalid value encountered:RuntimeWarning:dask")
def test_support_alternative_backends_mismatched_dtypes(xp, func, nfo):
    """Test mix-n-match of int and float arguments"""
    if func.__name__ in {'expn', 'polygamma', 'multigammaln', 'bdtr', 'bdtrc', 'bdtri',
                         'nbdtr', 'nbdtrc', 'nbdtri', 'pdtri'}:
        pytest.skip(f"dtypes for {func.__name__} make it a bad fit for this test.")
    dtypes = ['intp', 'float32', 'float64', 'float64'][:nfo.n_args]

    positive_only, dtypes = _skip_or_tweak_alternative_backends(
        xp, nfo, dtypes, (False,)*nfo.n_args
    )
    dtypes_np = [getattr(np, dtype) for dtype in dtypes]
    dtypes_xp = [getattr(xp, dtype) for dtype in dtypes]

    rng = np.random.default_rng(984254252920492019)
    iinfo = np.iinfo(np.intp)
    if nfo.test_large_ints:
        randint = partial(rng.integers, iinfo.min, iinfo.max + 1)
    else:
        randint = partial(rng.integers, -20, 21)
    args_np = [
        randint(size=1, dtype=np.int64),
        rng.standard_normal(size=1, dtype=np.float32),
        rng.standard_normal(size=1, dtype=np.float64),
        rng.standard_normal(size=1, dtype=np.float64),
    ][:nfo.n_args]
    args_np = [
        np.abs(arg) if cond else arg for arg, cond in zip(args_np, positive_only)
    ]

    args_xp = [xp.asarray(arg, dtype=dtype_xp)
               for arg, dtype_xp in zip(args_np, dtypes_xp)]
    args_np = [np.asarray(arg, dtype=dtype_np)
               for arg, dtype_np in zip(args_np, dtypes_np)]

    res = nfo.wrapper(*args_xp)  # Also wrapped by lazy_xp_function
    ref = nfo.func(*args_np)  # Unwrapped ufunc
    if (
            is_torch(xp)
            and xpx.default_dtype(xp) == xp.float32
            and "float64" not in dtypes
    ):
        # int64 promotes like float32 on torch with default dtype = float32
        # cast reference if needed
        ref = np.float32(ref)

    atol = 10 * np.finfo(ref.dtype).eps
    xp_assert_close(res, xp.asarray(ref), atol=atol)


@pytest.mark.xslow
@given(data=strategies.data())
@pytest.mark.fail_slow(5)
@pytest.mark.parametrize(
    'func,nfo', [make_xp_pytest_param(nfo.wrapper, nfo) for nfo in _special_funcs])
@pytest.mark.filterwarnings("ignore:invalid value encountered:RuntimeWarning:dask")
@pytest.mark.filterwarnings("ignore:divide by zero encountered:RuntimeWarning:dask")
@pytest.mark.filterwarnings("ignore:overflow encountered:RuntimeWarning:dask")
@pytest.mark.filterwarnings(
    "ignore:overflow encountered:RuntimeWarning:array_api_strict"
)
def test_support_alternative_backends_hypothesis(xp, func, nfo, data):
    if func.__name__ in {'expn', 'polygamma', 'multigammaln', 'bdtr', 'bdtrc', 'bdtri',
                         'nbdtr', 'nbdtrc', 'nbdtri', 'pdtri'}:
        pytest.skip(f"dtypes for {func.__name__} make it a bad fit for this test.")
    dtype = data.draw(strategies.sampled_from(['float32', 'float64', 'intp']))
    positive_only, dtypes = _skip_or_tweak_alternative_backends(
        xp, nfo, [dtype], (False,)*nfo.n_args
    )
    dtype_np = getattr(np, dtypes[0])
    dtype_xp = getattr(xp, dtypes[0])

    elements = {'allow_subnormal': False}
    # Most failures are due to NaN or infinity; uncomment to suppress them
    # elements['allow_infinity'] = False
    # elements['allow_nan'] = False
    if any(positive_only):
        elements['min_value'] = 0

    shapes, _ = data.draw(
        npst.mutually_broadcastable_shapes(num_shapes=nfo.n_args))
    args_np = [data.draw(npst.arrays(dtype_np, shape, elements=elements))
               for shape in shapes]

    args_xp = [xp.asarray(arg, dtype=dtype_xp) for arg in args_np]
    args_np = [np.asarray(arg, dtype=dtype_np) for arg in args_np]

    res = nfo.wrapper(*args_xp)  # Also wrapped by lazy_xp_function
    ref = nfo.func(*args_np)  # Unwrapped ufunc
    if (
            is_torch(xp)
            and xpx.default_dtype(xp) == xp.float32
            and dtype != "float64"
    ):
        # int64 promotes like float32 on torch with default dtype = float32
        # cast reference if needed
        ref = np.float32(ref)

    # When dtype_np is integer, the output dtype can be float
    atol = 0 if ref.dtype.kind in 'iu' else 10 * np.finfo(ref.dtype).eps
    xp_assert_close(res, xp.asarray(ref), atol=atol)


@pytest.mark.filterwarnings("ignore:numpy.core is deprecated:DeprecationWarning")
@pytest.mark.parametrize("func", [nfo.wrapper for nfo in _special_funcs])
def test_pickle(func):
    roundtrip = pickle.loads(pickle.dumps(func))
    assert roundtrip is func


@pytest.mark.parametrize("func", [nfo.wrapper for nfo in _special_funcs])
def test_repr(func):
    assert func.__name__ in repr(func)
    assert "locals" not in repr(func)


@pytest.mark.skipif(
    version.parse(np.__version__) < version.parse("2.2"),
    reason="Can't update ufunc __doc__ when SciPy is compiled vs. NumPy < 2.2")
@pytest.mark.parametrize('func', [nfo.wrapper for nfo in _special_funcs])
def test_doc(func):
    """xp_capabilities updates the docstring in place. 
    Make sure it does so exactly once, including when SCIPY_ARRAY_API is not set.
    """
    match = "has experimental support for Python Array API"
    assert func.__doc__.count(match) == 1


@pytest.mark.parametrize(
    'func,n_args,int_only,is_ufunc',
    [(nfo.wrapper, nfo.n_args, nfo.int_only, nfo.is_ufunc)
     for nfo in _special_funcs]
)
def test_ufunc_kwargs(func, n_args, int_only, is_ufunc):
    """Test that numpy-specific out= and dtype= keyword arguments
    of ufuncs still work when SCIPY_ARRAY_API is set.
    """
    if not is_ufunc:
        pytest.skip(f"{func.__name__} is not a ufunc.")
    if int_only is None:
        int_only = (False, ) * n_args
    # out=
    args = [
        np.asarray([.1, .2]) if not needs_int
        else np.asarray([1, 2])
        for needs_int in int_only
    ]
    out = np.empty(2)
    y = func(*args, out=out)
    xp_assert_close(y, out)

    # out= with out.dtype != args.dtype
    out = np.empty(2, dtype=np.float32)
    y = func(*args, out=out)
    xp_assert_close(y, out)

    if func.__name__ in {"bdtr", "bdtrc", "bdtri"}:
        # The below function evaluation will trigger a deprecation warning
        # with dtype=np.float32. This will go away if the trigger is actually
        # pulled on the deprecation.
        return

    # dtype=
    y = func(*args, dtype=np.float32)
    assert y.dtype == np.float32


@make_xp_test_case(special.chdtr)
def test_chdtr_gh21311(xp):
    # the edge case behavior of generic chdtr was not right; see gh-21311
    # be sure to test at least these cases
    # should add `np.nan` into the mix when gh-21317 is resolved
    x = np.asarray([-np.inf, -1., 0., 1., np.inf])
    v = x.reshape(-1, 1)
    ref = special.chdtr(v, x)
    res = special.chdtr(xp.asarray(v), xp.asarray(x))
    xp_assert_close(res, xp.asarray(ref))
