from multiprocessing import Pool
from multiprocessing.pool import Pool as PWL
import re
import math
from fractions import Fraction

import numpy as np
from numpy.testing import assert_equal, assert_
import pytest
from pytest import raises as assert_raises
import hypothesis.extra.numpy as npst
from hypothesis import given, strategies, reproduce_failure  # noqa: F401
from scipy.conftest import array_api_compatible, skip_xp_invalid_arg

from scipy._lib._array_api import (xp_assert_equal, xp_assert_close, is_numpy,
                                   xp_copy, is_array_api_strict)
from scipy._lib._util import (_aligned_zeros, check_random_state, MapWrapper,
                              getfullargspec_no_self, FullArgSpec,
                              rng_integers, _validate_int, _rename_parameter,
                              _contains_nan, _rng_html_rewrite, _lazywhere)
from scipy import cluster, interpolate, linalg, optimize, sparse, spatial, stats

skip_xp_backends = pytest.mark.skip_xp_backends


@pytest.mark.slow
def test__aligned_zeros():
    niter = 10

    def check(shape, dtype, order, align):
        err_msg = repr((shape, dtype, order, align))
        x = _aligned_zeros(shape, dtype, order, align=align)
        if align is None:
            align = np.dtype(dtype).alignment
        assert_equal(x.__array_interface__['data'][0] % align, 0)
        if hasattr(shape, '__len__'):
            assert_equal(x.shape, shape, err_msg)
        else:
            assert_equal(x.shape, (shape,), err_msg)
        assert_equal(x.dtype, dtype)
        if order == "C":
            assert_(x.flags.c_contiguous, err_msg)
        elif order == "F":
            if x.size > 0:
                # Size-0 arrays get invalid flags on NumPy 1.5
                assert_(x.flags.f_contiguous, err_msg)
        elif order is None:
            assert_(x.flags.c_contiguous, err_msg)
        else:
            raise ValueError()

    # try various alignments
    for align in [1, 2, 3, 4, 8, 16, 32, 64, None]:
        for n in [0, 1, 3, 11]:
            for order in ["C", "F", None]:
                for dtype in [np.uint8, np.float64]:
                    for shape in [n, (1, 2, 3, n)]:
                        for j in range(niter):
                            check(shape, dtype, order, align)


def test_check_random_state():
    # If seed is None, return the RandomState singleton used by np.random.
    # If seed is an int, return a new RandomState instance seeded with seed.
    # If seed is already a RandomState instance, return it.
    # Otherwise raise ValueError.
    rsi = check_random_state(1)
    assert_equal(type(rsi), np.random.RandomState)
    rsi = check_random_state(rsi)
    assert_equal(type(rsi), np.random.RandomState)
    rsi = check_random_state(None)
    assert_equal(type(rsi), np.random.RandomState)
    assert_raises(ValueError, check_random_state, 'a')
    rg = np.random.Generator(np.random.PCG64())
    rsi = check_random_state(rg)
    assert_equal(type(rsi), np.random.Generator)


def test_getfullargspec_no_self():
    p = MapWrapper(1)
    argspec = getfullargspec_no_self(p.__init__)
    assert_equal(argspec, FullArgSpec(['pool'], None, None, (1,), [],
                                      None, {}))
    argspec = getfullargspec_no_self(p.__call__)
    assert_equal(argspec, FullArgSpec(['func', 'iterable'], None, None, None,
                                      [], None, {}))

    class _rv_generic:
        def _rvs(self, a, b=2, c=3, *args, size=None, **kwargs):
            return None

    rv_obj = _rv_generic()
    argspec = getfullargspec_no_self(rv_obj._rvs)
    assert_equal(argspec, FullArgSpec(['a', 'b', 'c'], 'args', 'kwargs',
                                      (2, 3), ['size'], {'size': None}, {}))


def test_mapwrapper_serial():
    in_arg = np.arange(10.)
    out_arg = np.sin(in_arg)

    p = MapWrapper(1)
    assert_(p._mapfunc is map)
    assert_(p.pool is None)
    assert_(p._own_pool is False)
    out = list(p(np.sin, in_arg))
    assert_equal(out, out_arg)

    with assert_raises(RuntimeError):
        p = MapWrapper(0)


def test_pool():
    with Pool(2) as p:
        p.map(math.sin, [1, 2, 3, 4])


def test_mapwrapper_parallel():
    in_arg = np.arange(10.)
    out_arg = np.sin(in_arg)

    with MapWrapper(2) as p:
        out = p(np.sin, in_arg)
        assert_equal(list(out), out_arg)

        assert_(p._own_pool is True)
        assert_(isinstance(p.pool, PWL))
        assert_(p._mapfunc is not None)

    # the context manager should've closed the internal pool
    # check that it has by asking it to calculate again.
    with assert_raises(Exception) as excinfo:
        p(np.sin, in_arg)

    assert_(excinfo.type is ValueError)

    # can also set a PoolWrapper up with a map-like callable instance
    with Pool(2) as p:
        q = MapWrapper(p.map)

        assert_(q._own_pool is False)
        q.close()

        # closing the PoolWrapper shouldn't close the internal pool
        # because it didn't create it
        out = p.map(np.sin, in_arg)
        assert_equal(list(out), out_arg)


def test_rng_integers():
    rng = np.random.RandomState()

    # test that numbers are inclusive of high point
    arr = rng_integers(rng, low=2, high=5, size=100, endpoint=True)
    assert np.max(arr) == 5
    assert np.min(arr) == 2
    assert arr.shape == (100, )

    # test that numbers are inclusive of high point
    arr = rng_integers(rng, low=5, size=100, endpoint=True)
    assert np.max(arr) == 5
    assert np.min(arr) == 0
    assert arr.shape == (100, )

    # test that numbers are exclusive of high point
    arr = rng_integers(rng, low=2, high=5, size=100, endpoint=False)
    assert np.max(arr) == 4
    assert np.min(arr) == 2
    assert arr.shape == (100, )

    # test that numbers are exclusive of high point
    arr = rng_integers(rng, low=5, size=100, endpoint=False)
    assert np.max(arr) == 4
    assert np.min(arr) == 0
    assert arr.shape == (100, )

    # now try with np.random.Generator
    try:
        rng = np.random.default_rng()
    except AttributeError:
        return

    # test that numbers are inclusive of high point
    arr = rng_integers(rng, low=2, high=5, size=100, endpoint=True)
    assert np.max(arr) == 5
    assert np.min(arr) == 2
    assert arr.shape == (100, )

    # test that numbers are inclusive of high point
    arr = rng_integers(rng, low=5, size=100, endpoint=True)
    assert np.max(arr) == 5
    assert np.min(arr) == 0
    assert arr.shape == (100, )

    # test that numbers are exclusive of high point
    arr = rng_integers(rng, low=2, high=5, size=100, endpoint=False)
    assert np.max(arr) == 4
    assert np.min(arr) == 2
    assert arr.shape == (100, )

    # test that numbers are exclusive of high point
    arr = rng_integers(rng, low=5, size=100, endpoint=False)
    assert np.max(arr) == 4
    assert np.min(arr) == 0
    assert arr.shape == (100, )


class TestValidateInt:

    @pytest.mark.parametrize('n', [4, np.uint8(4), np.int16(4), np.array(4)])
    def test_validate_int(self, n):
        n = _validate_int(n, 'n')
        assert n == 4

    @pytest.mark.parametrize('n', [4.0, np.array([4]), Fraction(4, 1)])
    def test_validate_int_bad(self, n):
        with pytest.raises(TypeError, match='n must be an integer'):
            _validate_int(n, 'n')

    def test_validate_int_below_min(self):
        with pytest.raises(ValueError, match='n must be an integer not '
                                             'less than 0'):
            _validate_int(-1, 'n', 0)


class TestRenameParameter:
    # check that wrapper `_rename_parameter` for backward-compatible
    # keyword renaming works correctly

    # Example method/function that still accepts keyword `old`
    @_rename_parameter("old", "new")
    def old_keyword_still_accepted(self, new):
        return new

    # Example method/function for which keyword `old` is deprecated
    @_rename_parameter("old", "new", dep_version="1.9.0")
    def old_keyword_deprecated(self, new):
        return new

    def test_old_keyword_still_accepted(self):
        # positional argument and both keyword work identically
        res1 = self.old_keyword_still_accepted(10)
        res2 = self.old_keyword_still_accepted(new=10)
        res3 = self.old_keyword_still_accepted(old=10)
        assert res1 == res2 == res3 == 10

        # unexpected keyword raises an error
        message = re.escape("old_keyword_still_accepted() got an unexpected")
        with pytest.raises(TypeError, match=message):
            self.old_keyword_still_accepted(unexpected=10)

        # multiple values for the same parameter raises an error
        message = re.escape("old_keyword_still_accepted() got multiple")
        with pytest.raises(TypeError, match=message):
            self.old_keyword_still_accepted(10, new=10)
        with pytest.raises(TypeError, match=message):
            self.old_keyword_still_accepted(10, old=10)
        with pytest.raises(TypeError, match=message):
            self.old_keyword_still_accepted(new=10, old=10)

    @pytest.fixture
    def kwarg_lock(self):
        from threading import Lock
        return Lock()

    def test_old_keyword_deprecated(self, kwarg_lock):
        # positional argument and both keyword work identically,
        # but use of old keyword results in DeprecationWarning
        dep_msg = "Use of keyword argument `old` is deprecated"
        res1 = self.old_keyword_deprecated(10)
        res2 = self.old_keyword_deprecated(new=10)
        # pytest warning filter is not thread-safe, enforce serialization
        with kwarg_lock:
            with pytest.warns(DeprecationWarning, match=dep_msg):
                    res3 = self.old_keyword_deprecated(old=10)
        assert res1 == res2 == res3 == 10

        # unexpected keyword raises an error
        message = re.escape("old_keyword_deprecated() got an unexpected")
        with pytest.raises(TypeError, match=message):
            self.old_keyword_deprecated(unexpected=10)

        # multiple values for the same parameter raises an error and,
        # if old keyword is used, results in DeprecationWarning
        message = re.escape("old_keyword_deprecated() got multiple")
        with pytest.raises(TypeError, match=message):
            self.old_keyword_deprecated(10, new=10)
        with kwarg_lock:
            with pytest.raises(TypeError, match=message), \
                    pytest.warns(DeprecationWarning, match=dep_msg):
                    # breakpoint()
                    self.old_keyword_deprecated(10, old=10)
        with kwarg_lock:
            with pytest.raises(TypeError, match=message), \
                    pytest.warns(DeprecationWarning, match=dep_msg):
                    self.old_keyword_deprecated(new=10, old=10)


class TestContainsNaNTest:

    def test_policy(self):
        data = np.array([1, 2, 3, np.nan])

        contains_nan, nan_policy = _contains_nan(data, nan_policy="propagate")
        assert contains_nan
        assert nan_policy == "propagate"

        contains_nan, nan_policy = _contains_nan(data, nan_policy="omit")
        assert contains_nan
        assert nan_policy == "omit"

        msg = "The input contains nan values"
        with pytest.raises(ValueError, match=msg):
            _contains_nan(data, nan_policy="raise")

        msg = "nan_policy must be one of"
        with pytest.raises(ValueError, match=msg):
            _contains_nan(data, nan_policy="nan")

    def test_contains_nan(self):
        data1 = np.array([1, 2, 3])
        assert not _contains_nan(data1)[0]

        data2 = np.array([1, 2, 3, np.nan])
        assert _contains_nan(data2)[0]

        data3 = np.array([np.nan, 2, 3, np.nan])
        assert _contains_nan(data3)[0]

        data4 = np.array([[1, 2], [3, 4]])
        assert not _contains_nan(data4)[0]

        data5 = np.array([[1, 2], [3, np.nan]])
        assert _contains_nan(data5)[0]

    @skip_xp_invalid_arg
    def test_contains_nan_with_strings(self):
        data1 = np.array([1, 2, "3", np.nan])  # converted to string "nan"
        assert not _contains_nan(data1)[0]

        data2 = np.array([1, 2, "3", np.nan], dtype='object')
        assert _contains_nan(data2)[0]

        data3 = np.array([["1", 2], [3, np.nan]])  # converted to string "nan"
        assert not _contains_nan(data3)[0]

        data4 = np.array([["1", 2], [3, np.nan]], dtype='object')
        assert _contains_nan(data4)[0]

    @skip_xp_backends('jax.numpy',
                      reason="JAX arrays do not support item assignment")
    @pytest.mark.usefixtures("skip_xp_backends")
    @array_api_compatible
    @pytest.mark.parametrize("nan_policy", ['propagate', 'omit', 'raise'])
    def test_array_api(self, xp, nan_policy):
        rng = np.random.default_rng(932347235892482)
        x0 = rng.random(size=(2, 3, 4))
        x = xp.asarray(x0)
        x_nan = xp_copy(x, xp=xp)
        x_nan[1, 2, 1] = np.nan

        contains_nan, nan_policy_out = _contains_nan(x, nan_policy=nan_policy)
        assert not contains_nan
        assert nan_policy_out == nan_policy

        if nan_policy == 'raise':
            message = 'The input contains...'
            with pytest.raises(ValueError, match=message):
                _contains_nan(x_nan, nan_policy=nan_policy)
        elif nan_policy == 'omit' and not is_numpy(xp):
            message = "`nan_policy='omit' is incompatible..."
            with pytest.raises(ValueError, match=message):
                _contains_nan(x_nan, nan_policy=nan_policy)
        elif nan_policy == 'propagate':
            contains_nan, nan_policy_out = _contains_nan(
                x_nan, nan_policy=nan_policy)
            assert contains_nan
            assert nan_policy_out == nan_policy


def test__rng_html_rewrite():
    def mock_str():
        lines = [
            'np.random.default_rng(8989843)',
            'np.random.default_rng(seed)',
            'np.random.default_rng(0x9a71b21474694f919882289dc1559ca)',
            ' bob ',
        ]
        return lines

    res = _rng_html_rewrite(mock_str)()
    ref = [
        'np.random.default_rng()',
        'np.random.default_rng(seed)',
        'np.random.default_rng()',
        ' bob ',
    ]

    assert res == ref


class TestTransitionToRNG:
    def kmeans(self, **kwargs):
        rng = np.random.default_rng(3458934594269824562)
        return cluster.vq.kmeans2(rng.random(size=(20, 3)), 3, **kwargs)

    def kmeans2(self, **kwargs):
        rng = np.random.default_rng(3458934594269824562)
        return cluster.vq.kmeans2(rng.random(size=(20, 3)), 3, **kwargs)

    def barycentric(self, **kwargs):
        rng = np.random.default_rng(3458934594269824562)
        x1, x2, y1 = rng.random((3, 10))
        f = interpolate.BarycentricInterpolator(x1, y1, **kwargs)
        return f(x2)

    def clarkson_woodruff_transform(self, **kwargs):
        rng = np.random.default_rng(3458934594269824562)
        return linalg.clarkson_woodruff_transform(rng.random((10, 10)), 3, **kwargs)

    def basinhopping(self, **kwargs):
        rng = np.random.default_rng(3458934594269824562)
        return optimize.basinhopping(optimize.rosen, rng.random(3), **kwargs).x

    def opt(self, fun, **kwargs):
        rng = np.random.default_rng(3458934594269824562)
        bounds = optimize.Bounds(-rng.random(3) * 10, rng.random(3) * 10)
        return fun(optimize.rosen, bounds, **kwargs).x

    def differential_evolution(self, **kwargs):
        return self.opt(optimize.differential_evolution, **kwargs)

    def dual_annealing(self, **kwargs):
        return self.opt(optimize.dual_annealing, **kwargs)

    def check_grad(self, **kwargs):
        rng = np.random.default_rng(3458934594269824562)
        x = rng.random(3)
        return optimize.check_grad(optimize.rosen, optimize.rosen_der, x,
                                   direction='random', **kwargs)

    def random_array(self, **kwargs):
        return sparse.random_array((10, 10), density=1.0, **kwargs).toarray()

    def random(self, **kwargs):
        return sparse.random(10, 10, density=1.0, **kwargs).toarray()

    def rand(self, **kwargs):
        return sparse.rand(10, 10, density=1.0, **kwargs).toarray()

    def svds(self, **kwargs):
        rng = np.random.default_rng(3458934594269824562)
        A = rng.random((10, 10))
        return sparse.linalg.svds(A, **kwargs)

    def random_rotation(self, **kwargs):
        return spatial.transform.Rotation.random(3, **kwargs).as_matrix()

    def goodness_of_fit(self, **kwargs):
        rng = np.random.default_rng(3458934594269824562)
        data = rng.random(100)
        return stats.goodness_of_fit(stats.laplace, data, **kwargs).pvalue

    def permutation_test(self, **kwargs):
        rng = np.random.default_rng(3458934594269824562)
        data = tuple(rng.random((2, 100)))
        def statistic(x, y, axis): return np.mean(x, axis=axis) - np.mean(y, axis=axis)
        return stats.permutation_test(data, statistic, **kwargs).pvalue

    def bootstrap(self, **kwargs):
        rng = np.random.default_rng(3458934594269824562)
        data = (rng.random(100),)
        return stats.bootstrap(data, np.mean, **kwargs).confidence_interval

    def dunnett(self, **kwargs):
        rng = np.random.default_rng(3458934594269824562)
        x, y, control = rng.random((3, 100))
        return stats.dunnett(x, y, control=control, **kwargs).pvalue

    def sobol_indices(self, **kwargs):
        def f_ishigami(x): return (np.sin(x[0]) + 7 * np.sin(x[1]) ** 2
                                   + 0.1 * (x[2] ** 4) * np.sin(x[0]))
        dists = [stats.uniform(loc=-np.pi, scale=2 * np.pi),
                 stats.uniform(loc=-np.pi, scale=2 * np.pi),
                 stats.uniform(loc=-np.pi, scale=2 * np.pi)]
        res = stats.sobol_indices(func=f_ishigami, n=1024, dists=dists, **kwargs)
        return res.first_order

    def qmc_engine(self, engine, **kwargs):
        qrng = engine(d=1, **kwargs)
        return qrng.random(4)

    def halton(self, **kwargs):
        return self.qmc_engine(stats.qmc.Halton, **kwargs)

    def sobol(self, **kwargs):
        return self.qmc_engine(stats.qmc.Sobol, **kwargs)

    def latin_hypercube(self, **kwargs):
        return self.qmc_engine(stats.qmc.LatinHypercube, **kwargs)

    def poisson_disk(self, **kwargs):
        return self.qmc_engine(stats.qmc.PoissonDisk, **kwargs)

    def multivariate_normal_qmc(self, **kwargs):
        X = stats.qmc.MultivariateNormalQMC([0], **kwargs)
        return X.random(4)

    def multinomial_qmc(self, **kwargs):
        X = stats.qmc.MultinomialQMC([0.5, 0.5], 4, **kwargs)
        return X.random(4)

    def permutation_method(self, **kwargs):
        rng = np.random.default_rng(3458934594269824562)
        data = tuple(rng.random((2, 100)))
        method = stats.PermutationMethod(**kwargs)
        return stats.pearsonr(*data, method=method).pvalue

    def bootstrap_method(self, **kwargs):
        rng = np.random.default_rng(3458934594269824562)
        data = tuple(rng.random((2, 100)))
        res = stats.pearsonr(*data)
        method = stats.BootstrapMethod(**kwargs)
        return res.confidence_interval(method=method)

    @pytest.mark.fail_slow(10)
    @pytest.mark.slow
    @pytest.mark.parametrize("method, arg_name", [
        (kmeans, "seed"),
        (kmeans2, "seed"),
        (barycentric, "random_state"),
        (clarkson_woodruff_transform, "seed"),
        (basinhopping, "seed"),
        (differential_evolution, "seed"),
        (dual_annealing, "seed"),
        (check_grad, "seed"),
        (random_array, 'random_state'),
        (random, 'random_state'),
        (rand, 'random_state'),
        (svds, "random_state"),
        (random_rotation, "random_state"),
        (goodness_of_fit, "random_state"),
        (permutation_test, "random_state"),
        (bootstrap, "random_state"),
        (permutation_method, "random_state"),
        (bootstrap_method, "random_state"),
        (dunnett, "random_state"),
        (sobol_indices, "random_state"),
        (halton, "seed"),
        (sobol, "seed"),
        (latin_hypercube, "seed"),
        (poisson_disk, "seed"),
        (multivariate_normal_qmc, "seed"),
        (multinomial_qmc, "seed"),
    ])
    def test_rng_deterministic(self, method, arg_name):
        np.random.seed(None)
        seed = 2949672964

        rng = np.random.default_rng(seed)
        message = "got multiple values for argument now known as `rng`"
        with pytest.raises(TypeError, match=message):
            method(self, **{'rng': rng, arg_name: seed})

        rng = np.random.default_rng(seed)
        res1 = method(self, rng=rng)
        res2 = method(self, rng=seed)
        assert_equal(res2, res1)

        if method.__name__ in {"dunnett", "sobol_indices"}:
            # the two kwargs have essentially the same behavior for these functions
            res3 = method(self, **{arg_name: seed})
            assert_equal(res3, res1)
            return

        rng = np.random.RandomState(seed)
        res1 = method(self, **{arg_name: rng})
        res2 = method(self, **{arg_name: seed})

        if method.__name__ in {"halton", "sobol", "latin_hypercube", "poisson_disk",
                               "multivariate_normal_qmc", "multinomial_qmc"}:
            # For these, passing `random_state=RandomState(seed)` is not the same as
            # passing integer `seed`.
            res1b = method(self, **{arg_name: np.random.RandomState(seed)})
            assert_equal(res1b, res1)
            res2b = method(self, **{arg_name: seed})
            assert_equal(res2b, res2)
            return

        np.random.seed(seed)
        res3 = method(self, **{arg_name: None})
        assert_equal(res2, res1)
        assert_equal(res3, res1)


class TestLazywhere:
    n_arrays = strategies.integers(min_value=1, max_value=3)
    rng_seed = strategies.integers(min_value=1000000000, max_value=9999999999)
    dtype = strategies.sampled_from((np.float32, np.float64))
    p = strategies.floats(min_value=0, max_value=1)
    data = strategies.data()

    @pytest.mark.fail_slow(10)
    @pytest.mark.filterwarnings('ignore::RuntimeWarning')  # overflows, etc.
    @skip_xp_backends('jax.numpy',
                      reason="JAX arrays do not support item assignment")
    @pytest.mark.usefixtures("skip_xp_backends")
    @array_api_compatible
    @given(n_arrays=n_arrays, rng_seed=rng_seed, dtype=dtype, p=p, data=data)
    @pytest.mark.thread_unsafe
    def test_basic(self, n_arrays, rng_seed, dtype, p, data, xp):
        mbs = npst.mutually_broadcastable_shapes(num_shapes=n_arrays+1,
                                                 min_side=0)
        input_shapes, result_shape = data.draw(mbs)
        cond_shape, *shapes = input_shapes
        elements = {'allow_subnormal': False}  # cupy/cupy#8382
        fillvalue = xp.asarray(data.draw(npst.arrays(dtype=dtype, shape=tuple(),
                                                     elements=elements)))
        float_fillvalue = float(fillvalue)
        arrays = [xp.asarray(data.draw(npst.arrays(dtype=dtype, shape=shape)))
                  for shape in shapes]

        def f(*args):
            return sum(arg for arg in args)

        def f2(*args):
            return sum(arg for arg in args) / 2

        rng = np.random.default_rng(rng_seed)
        cond = xp.asarray(rng.random(size=cond_shape) > p)

        res1 = _lazywhere(cond, arrays, f, fillvalue)
        res2 = _lazywhere(cond, arrays, f, f2=f2)
        if not is_array_api_strict(xp):
            res3 = _lazywhere(cond, arrays, f, float_fillvalue)

        # Ensure arrays are at least 1d to follow sane type promotion rules.
        # This can be removed when minimum supported NumPy is 2.0
        if xp == np:
            cond, fillvalue, *arrays = np.atleast_1d(cond, fillvalue, *arrays)

        ref1 = xp.where(cond, f(*arrays), fillvalue)
        ref2 = xp.where(cond, f(*arrays), f2(*arrays))
        if not is_array_api_strict(xp):
            # Array API standard doesn't currently define behavior when fillvalue is a
            # Python scalar. When it does, test can be run with array_api_strict, too.
            ref3 = xp.where(cond, f(*arrays), float_fillvalue)

        if xp == np:  # because we ensured arrays are at least 1d
            ref1 = ref1.reshape(result_shape)
            ref2 = ref2.reshape(result_shape)
            ref3 = ref3.reshape(result_shape)

        xp_assert_close(res1, ref1, rtol=2e-16)
        xp_assert_equal(res2, ref2)
        if not is_array_api_strict(xp):
            xp_assert_equal(res3, ref3)
