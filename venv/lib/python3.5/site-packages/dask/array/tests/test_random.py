import pytest
pytest.importorskip('numpy')

import numpy as np

import dask.array as da
from dask.array.core import Array
from dask.array.random import random, exponential, normal
from dask.array.utils import assert_eq
from dask.multiprocessing import get as mpget
from dask.multiprocessing import _dumps, _loads


def test_RandomState():
    state = da.random.RandomState(5)
    x = state.normal(10, 1, size=10, chunks=5)
    assert_eq(x, x)

    state = da.random.RandomState(5)
    y = state.normal(10, 1, size=10, chunks=5)
    assert_eq(x, y)


def test_concurrency():
    state = da.random.RandomState(5)
    x = state.normal(10, 1, size=10, chunks=2)

    state = da.random.RandomState(5)
    y = state.normal(10, 1, size=10, chunks=2)
    assert (x.compute(get=mpget) == y.compute(get=mpget)).all()


def test_doc_randomstate():
    assert 'mean' in da.random.RandomState(5).normal.__doc__


def test_serializability():
    state = da.random.RandomState(5)
    x = state.normal(10, 1, size=10, chunks=5)

    y = _loads(_dumps(x))

    assert_eq(x, y)


def test_determinisim_through_dask_values():
    samples_1 = da.random.RandomState(42).normal(size=1000, chunks=10)
    samples_2 = da.random.RandomState(42).normal(size=1000, chunks=10)

    assert set(samples_1.dask) == set(samples_2.dask)
    assert_eq(samples_1, samples_2)


def test_randomstate_consistent_names():
    state1 = da.random.RandomState(42)
    state2 = da.random.RandomState(42)
    assert (sorted(state1.normal(size=(100, 100), chunks=(10, 10)).dask) ==
            sorted(state2.normal(size=(100, 100), chunks=(10, 10)).dask))
    assert (sorted(state1.normal(size=100, loc=4.5, scale=5.0, chunks=10).dask) ==
            sorted(state2.normal(size=100, loc=4.5, scale=5.0, chunks=10).dask))


def test_random():
    a = random((10, 10), chunks=(5, 5))
    assert isinstance(a, Array)
    assert isinstance(a.name, str) and a.name
    assert a.shape == (10, 10)
    assert a.chunks == ((5, 5), (5, 5))

    x = set(np.array(a).flat)

    assert len(x) > 90


def test_parametrized_random_function():
    a = exponential(1000, (10, 10), chunks=(5, 5))
    assert isinstance(a, Array)
    assert isinstance(a.name, str) and a.name
    assert a.shape == (10, 10)
    assert a.chunks == ((5, 5), (5, 5))

    x = np.array(a)
    assert 10 < x.mean() < 100000

    y = set(x.flat)
    assert len(y) > 90


def test_kwargs():
    a = normal(loc=10.0, scale=0.1, size=(10, 10), chunks=(5, 5))
    assert isinstance(a, Array)
    x = np.array(a)
    assert 8 < x.mean() < 12


def test_unique_names():
    a = random((10, 10), chunks=(5, 5))
    b = random((10, 10), chunks=(5, 5))

    assert a.name != b.name


def test_docs():
    assert 'exponential' in exponential.__doc__
    assert 'exponential' in exponential.__name__


def test_can_make_really_big_random_array():
    normal(10, 1, (1000000, 1000000), chunks=(100000, 100000))


def test_random_seed():
    da.random.seed(123)
    x = da.random.normal(size=10, chunks=5)
    y = da.random.normal(size=10, chunks=5)

    da.random.seed(123)
    a = da.random.normal(size=10, chunks=5)
    b = da.random.normal(size=10, chunks=5)

    assert_eq(x, a)
    assert_eq(y, b)


def test_consistent_across_sizes():
    x1 = da.random.RandomState(123).random(20, chunks=20)
    x2 = da.random.RandomState(123).random(100, chunks=20)[:20]
    x3 = da.random.RandomState(123).random(200, chunks=20)[:20]
    assert_eq(x1, x2)
    assert_eq(x1, x3)


def test_random_all():
    da.random.beta(1, 2, size=5, chunks=3).compute()
    da.random.binomial(10, 0.5, size=5, chunks=3).compute()
    da.random.chisquare(1, size=5, chunks=3).compute()
    da.random.exponential(1, size=5, chunks=3).compute()
    da.random.f(1, 2, size=5, chunks=3).compute()
    da.random.gamma(5, 1, size=5, chunks=3).compute()
    da.random.geometric(1, size=5, chunks=3).compute()
    da.random.gumbel(1, size=5, chunks=3).compute()
    da.random.hypergeometric(1, 2, 3, size=5, chunks=3).compute()
    da.random.laplace(size=5, chunks=3).compute()
    da.random.logistic(size=5, chunks=3).compute()
    da.random.lognormal(size=5, chunks=3).compute()
    da.random.logseries(0.5, size=5, chunks=3).compute()
    da.random.multinomial(20, [1 / 6.] * 6, size=5, chunks=3).compute()
    da.random.negative_binomial(5, 0.5, size=5, chunks=3).compute()
    da.random.noncentral_chisquare(2, 2, size=5, chunks=3).compute()

    da.random.noncentral_f(2, 2, 3, size=5, chunks=3).compute()
    da.random.normal(2, 2, size=5, chunks=3).compute()
    da.random.pareto(1, size=5, chunks=3).compute()
    da.random.poisson(size=5, chunks=3).compute()

    da.random.power(1, size=5, chunks=3).compute()
    da.random.rayleigh(size=5, chunks=3).compute()
    da.random.random_sample(size=5, chunks=3).compute()

    da.random.triangular(1, 2, 3, size=5, chunks=3).compute()
    da.random.uniform(size=5, chunks=3).compute()
    da.random.vonmises(2, 3, size=5, chunks=3).compute()
    da.random.wald(1, 2, size=5, chunks=3).compute()

    da.random.weibull(2, size=5, chunks=3).compute()
    da.random.zipf(2, size=5, chunks=3).compute()

    da.random.standard_cauchy(size=5, chunks=3).compute()
    da.random.standard_exponential(size=5, chunks=3).compute()
    da.random.standard_gamma(2, size=5, chunks=3).compute()
    da.random.standard_normal(size=5, chunks=3).compute()
    da.random.standard_t(2, size=5, chunks=3).compute()


@pytest.mark.skipif(not hasattr(np,'broadcast_to'),
                    reason='requires numpy 1.10 method "broadcast_to"')
def test_array_broadcasting():
    arr = np.arange(6).reshape((2, 3))
    daones = da.ones((2, 3, 4), chunks=3)
    assert da.random.poisson(arr, chunks=3).compute().shape == (2, 3)

    for x in (arr, daones):
        y = da.random.normal(x, 2, chunks=3)
        assert y.shape == x.shape
        assert y.compute().shape == x.shape

    y = da.random.normal(daones, 2, chunks=3)
    assert set(daones.dask).issubset(set(y.dask))

    assert da.random.normal(np.ones((1, 4)),
                            da.ones((2, 3, 4), chunks=(2, 3, 4)),
                            chunks=(2, 3, 4)).compute().shape == (2, 3, 4)
    assert da.random.normal(scale=np.ones((1, 4)),
                            loc=da.ones((2, 3, 4), chunks=(2, 3, 4)),
                            size=(2, 2, 3, 4),
                            chunks=(2, 2, 3, 4)).compute().shape == (2, 2, 3, 4)

    with pytest.raises(ValueError):
        da.random.normal(arr, np.ones((3, 1)), size=(2, 3, 4), chunks=3)

    for o in (np.ones(100), da.ones(100, chunks=(50,)), 1):
        a = da.random.normal(1000 * o, 0.01, chunks=(50,))
        assert 800 < a.mean().compute() < 1200

    # ensure that mis-matched chunks align well
    x = np.arange(10)**3
    y = da.from_array(x, chunks=(1,))
    z = da.random.normal(y, 0.01, chunks=(10,))

    assert 0.8 < z.mean().compute() / x.mean() < 1.2


def test_multinomial():
    for size, chunks in [(5, 3), ((5, 4), (2, 3))]:
        x = da.random.multinomial(20, [1 / 6.] * 6, size=size, chunks=chunks)
        y = np.random.multinomial(20, [1 / 6.] * 6, size=size)

        assert x.shape == y.shape == x.compute().shape


def test_choice():
    np_dtype = np.random.choice(1, size=()).dtype
    size = (10, 3)
    chunks = 4
    x = da.random.choice(3, size=size, chunks=chunks)
    assert x.dtype == np_dtype
    assert x.shape == size
    res = x.compute()
    assert res.dtype == np_dtype
    assert res.shape == size

    np_a = np.array([1, 3, 5, 7, 9], dtype='f8')
    da_a = da.from_array(np_a, chunks=2)

    for a in [np_a, da_a]:
        x = da.random.choice(a, size=size, chunks=chunks)
        res = x.compute()
        assert x.dtype == np_a.dtype
        assert res.dtype == np_a.dtype
        assert set(np.unique(res)).issubset(np_a)

    np_p = np.array([0, 0.2, 0.2, 0.3, 0.3])
    da_p = da.from_array(np_p, chunks=2)

    for a, p in [(da_a, np_p), (np_a, da_p)]:
        x = da.random.choice(a, size=size, chunks=chunks, p=p)
        res = x.compute()
        assert x.dtype == np_a.dtype
        assert res.dtype == np_a.dtype
        assert set(np.unique(res)).issubset(np_a[1:])

    np_dtype = np.random.choice(1, size=(), p=np.array([1])).dtype
    x = da.random.choice(5, size=size, chunks=chunks, p=np_p)
    res = x.compute()
    assert x.dtype == np_dtype
    assert res.dtype == np_dtype

    errs = [(-1, None),             # negative a
            (np_a[:, None], None),  # a must be 1D
            (np_a, np_p[:, None]),  # p must be 1D
            (np_a, np_p[:-2]),      # a and p must match
            (3, np_p),              # a and p must match
            (4, [0.2, 0.2, 0.3])]   # p must sum to 1

    for (a, p) in errs:
        with pytest.raises(ValueError):
            da.random.choice(a, size=size, chunks=chunks, p=p)
