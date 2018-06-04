# -*- coding: utf-8 -*-

import os
import pytest
from operator import add, mul
import subprocess
import sys
from toolz import merge

import dask
from dask import delayed
from dask.base import (compute, tokenize, normalize_token, normalize_function,
                       visualize, persist, function_cache, is_dask_collection,
                       DaskMethodsMixin, optimize, unpack_collections)
from dask.delayed import Delayed
from dask.utils import tmpdir, tmpfile, ignoring
from dask.utils_test import inc, dec
from dask.compatibility import long, unicode


def import_or_none(path):
    with ignoring():
        return pytest.importorskip(path)
    return None


tz = pytest.importorskip('toolz')
da = import_or_none('dask.array')
db = import_or_none('dask.bag')
dd = import_or_none('dask.dataframe')
np = import_or_none('numpy')
sp = import_or_none('scipy.sparse')
pd = import_or_none('pandas')


def f1(a, b, c=1):
    pass


def f2(a, b=1, c=2):
    pass


def f3(a):
    pass


def test_normalize_function():

    assert normalize_function(f2)

    assert normalize_function(lambda a: a)

    assert (normalize_function(tz.partial(f2, b=2)) ==
            normalize_function(tz.partial(f2, b=2)))

    assert (normalize_function(tz.partial(f2, b=2)) !=
            normalize_function(tz.partial(f2, b=3)))

    assert (normalize_function(tz.partial(f1, b=2)) !=
            normalize_function(tz.partial(f2, b=2)))

    assert (normalize_function(tz.compose(f2, f3)) ==
            normalize_function(tz.compose(f2, f3)))

    assert (normalize_function(tz.compose(f2, f3)) !=
            normalize_function(tz.compose(f2, f1)))

    assert normalize_function(tz.curry(f2)) == normalize_function(tz.curry(f2))
    assert normalize_function(tz.curry(f2)) != normalize_function(tz.curry(f1))
    assert (normalize_function(tz.curry(f2, b=1)) ==
            normalize_function(tz.curry(f2, b=1)))
    assert (normalize_function(tz.curry(f2, b=1)) !=
            normalize_function(tz.curry(f2, b=2)))


def test_tokenize():
    a = (1, 2, 3)
    assert isinstance(tokenize(a), (str, bytes))


@pytest.mark.skipif('not np')
def test_tokenize_numpy_array_consistent_on_values():
    assert (tokenize(np.random.RandomState(1234).random_sample(1000)) ==
            tokenize(np.random.RandomState(1234).random_sample(1000)))


@pytest.mark.skipif('not np')
def test_tokenize_numpy_array_supports_uneven_sizes():
    tokenize(np.random.random(7).astype(dtype='i2'))


@pytest.mark.skipif('not np')
def test_tokenize_discontiguous_numpy_array():
    tokenize(np.random.random(8)[::2])


@pytest.mark.skipif('not np')
def test_tokenize_numpy_datetime():
    tokenize(np.array(['2000-01-01T12:00:00'], dtype='M8[ns]'))


@pytest.mark.skipif('not np')
def test_tokenize_numpy_scalar():
    assert tokenize(np.array(1.0, dtype='f8')) == tokenize(np.array(1.0, dtype='f8'))
    assert (tokenize(np.array([(1, 2)], dtype=[('a', 'i4'), ('b', 'i8')])[0]) ==
            tokenize(np.array([(1, 2)], dtype=[('a', 'i4'), ('b', 'i8')])[0]))


@pytest.mark.skipif('not np')
def test_tokenize_numpy_array_on_object_dtype():
    assert (tokenize(np.array(['a', 'aa', 'aaa'], dtype=object)) ==
            tokenize(np.array(['a', 'aa', 'aaa'], dtype=object)))
    assert (tokenize(np.array(['a', None, 'aaa'], dtype=object)) ==
            tokenize(np.array(['a', None, 'aaa'], dtype=object)))
    assert (tokenize(np.array([(1, 'a'), (1, None), (1, 'aaa')], dtype=object)) ==
            tokenize(np.array([(1, 'a'), (1, None), (1, 'aaa')], dtype=object)))
    if sys.version_info[0] == 2:
        assert (tokenize(np.array([unicode("Rebeca Alón", encoding="utf-8")], dtype=object)) ==
                tokenize(np.array([unicode("Rebeca Alón", encoding="utf-8")], dtype=object)))


@pytest.mark.skipif('not np')
def test_tokenize_numpy_memmap():
    with tmpfile('.npy') as fn:
        x = np.arange(5)
        np.save(fn, x)
        y = tokenize(np.load(fn, mmap_mode='r'))

    with tmpfile('.npy') as fn:
        x = np.arange(5)
        np.save(fn, x)
        z = tokenize(np.load(fn, mmap_mode='r'))

    assert y != z

    with tmpfile('.npy') as fn:
        x = np.random.normal(size=(10, 10))
        np.save(fn, x)
        mm = np.load(fn, mmap_mode='r')
        mm2 = np.load(fn, mmap_mode='r')
        a = tokenize(mm[0, :])
        b = tokenize(mm[1, :])
        c = tokenize(mm[0:3, :])
        d = tokenize(mm[:, 0])
        assert len(set([a, b, c, d])) == 4
        assert tokenize(mm) == tokenize(mm2)
        assert tokenize(mm[1, :]) == tokenize(mm2[1, :])


@pytest.mark.skipif('not np')
def test_tokenize_numpy_memmap_no_filename():
    # GH 1562:
    with tmpfile('.npy') as fn1, tmpfile('.npy') as fn2:
        x = np.arange(5)
        np.save(fn1, x)
        np.save(fn2, x)

        a = np.load(fn1, mmap_mode='r')
        b = a + a
        assert tokenize(b) == tokenize(b)


@pytest.mark.skipif('not np')
def test_tokenize_numpy_ufunc_consistent():
    assert tokenize(np.sin) == '02106e2c67daf452fb480d264e0dac21'
    assert tokenize(np.cos) == 'c99e52e912e4379882a9a4b387957a0b'

    # Make a ufunc that isn't in the numpy namespace. Similar to
    # any found in other packages.
    inc = np.frompyfunc(lambda x: x + 1, 1, 1)
    assert tokenize(inc) == tokenize(inc)


def test_tokenize_partial_func_args_kwargs_consistent():
    f = tz.partial(f3, f2, c=f1)
    res = normalize_token(f)
    sol = (b'cdask.tests.test_base\nf3\np0\n.',
           (b'cdask.tests.test_base\nf2\np0\n.',),
           (('c', b'cdask.tests.test_base\nf1\np0\n.'),))
    assert res == sol


def test_normalize_base():
    for i in [1, long(1), 1.1, '1', slice(1, 2, 3)]:
        assert normalize_token(i) is i


@pytest.mark.skipif('not pd')
def test_tokenize_pandas():
    a = pd.DataFrame({'x': [1, 2, 3], 'y': ['4', 'asd', None]}, index=[1, 2, 3])
    b = pd.DataFrame({'x': [1, 2, 3], 'y': ['4', 'asd', None]}, index=[1, 2, 3])

    assert tokenize(a) == tokenize(b)
    b.index.name = 'foo'
    assert tokenize(a) != tokenize(b)

    a = pd.DataFrame({'x': [1, 2, 3], 'y': ['a', 'b', 'a']})
    b = pd.DataFrame({'x': [1, 2, 3], 'y': ['a', 'b', 'a']})
    a['z'] = a.y.astype('category')
    assert tokenize(a) != tokenize(b)
    b['z'] = a.y.astype('category')
    assert tokenize(a) == tokenize(b)


def test_tokenize_kwargs():
    assert tokenize(5, x=1) == tokenize(5, x=1)
    assert tokenize(5) != tokenize(5, x=1)
    assert tokenize(5, x=1) != tokenize(5, x=2)
    assert tokenize(5, x=1) != tokenize(5, y=1)


def test_tokenize_same_repr():
    class Foo(object):

        def __init__(self, x):
            self.x = x

        def __repr__(self):
            return 'a foo'

    assert tokenize(Foo(1)) != tokenize(Foo(2))


def test_tokenize_method():
    class Foo(object):
        def __init__(self, x):
            self.x = x

        def __dask_tokenize__(self):
            return self.x

    a, b = Foo(1), Foo(2)
    assert tokenize(a) == tokenize(a)
    assert tokenize(a) != tokenize(b)

    # dispatch takes precedence
    before = tokenize(a)
    normalize_token.register(Foo, lambda self: self.x + 1)
    after = tokenize(a)
    assert before != after


@pytest.mark.skipif('not np')
def test_tokenize_sequences():
    assert tokenize([1]) != tokenize([2])
    assert tokenize([1]) != tokenize((1,))
    assert tokenize([1]) == tokenize([1])

    x = np.arange(2000)  # long enough to drop information in repr
    y = np.arange(2000)
    y[1000] = 0  # middle isn't printed in repr
    assert tokenize([x]) != tokenize([y])


def test_tokenize_dict():
    assert tokenize({'x': 1, 1: 'x'}) == tokenize({'x': 1, 1: 'x'})


def test_tokenize_set():
    assert tokenize({1, 2, 'x', (1, 'x')}) == tokenize({1, 2, 'x', (1, 'x')})


def test_tokenize_ordered_dict():
    with ignoring(ImportError):
        from collections import OrderedDict
        a = OrderedDict([('a', 1), ('b', 2)])
        b = OrderedDict([('a', 1), ('b', 2)])
        c = OrderedDict([('b', 2), ('a', 1)])

        assert tokenize(a) == tokenize(b)
        assert tokenize(a) != tokenize(c)


@pytest.mark.skipif('not np')
def test_tokenize_object_array_with_nans():
    a = np.array([u'foo', u'Jos\xe9', np.nan], dtype='O')
    assert tokenize(a) == tokenize(a)


@pytest.mark.parametrize('x', [1, True, 'a', b'a', 1.0, 1j, 1.0j,
                               [], (), {}, None, str, int])
def test_tokenize_base_types(x):
    assert tokenize(x) == tokenize(x), x


@pytest.mark.skipif('not np')
def test_tokenize_numpy_matrix():
    rng = np.random.RandomState(1234)
    a = np.asmatrix(rng.rand(100))
    b = a.copy()
    assert tokenize(a) == tokenize(b)

    b[:10] = 1
    assert tokenize(a) != tokenize(b)


@pytest.mark.skipif('not sp')
@pytest.mark.parametrize('cls_name',
                         ('dia', 'bsr', 'coo', 'csc', 'csr', 'dok', 'lil'))
def test_tokenize_dense_sparse_array(cls_name):
    rng = np.random.RandomState(1234)

    with pytest.warns(None):
        # ignore scipy.sparse.SparseEfficiencyWarning
        a = sp.rand(10, 10000, random_state=rng).asformat(cls_name)
    b = a.copy()

    assert tokenize(a) == tokenize(b)

    # modifying the data values
    if hasattr(b, 'data'):
        b.data[:10] = 1
    elif cls_name == 'dok':
        b[3, 3] = 1
    else:
        raise ValueError

    assert tokenize(a) != tokenize(b)

    # modifying the data indices
    with pytest.warns(None):
        b = a.copy().asformat('coo')
        b.row[:10] = np.arange(10)
        b = b.asformat(cls_name)
    assert tokenize(a) != tokenize(b)


def test_is_dask_collection():
    class DummyCollection(object):
        def __init__(self, dsk=None):
            self.dask = dsk

        def __dask_graph__(self):
            return self.dask

    x = delayed(1) + 2
    assert is_dask_collection(x)
    assert not is_dask_collection(2)
    assert is_dask_collection(DummyCollection({}))
    assert not is_dask_collection(DummyCollection())
    assert not is_dask_collection(DummyCollection)


def test_unpack_collections():
    a = delayed(1) + 5
    b = a + 1
    c = a + 2

    def build(a, b, c, iterator):
        return (a, b,               # Top-level collections
                {'a': a,            # dict
                 a: b,              # collections as keys
                 'b': [1, 2, [b]],  # list
                 'c': 10,           # other builtins pass through unchanged
                 'd': (c, 2),       # tuple
                 'e': {a, 2, 3}},   # set
                iterator)           # Iterator

    args = build(a, b, c, (i for i in [a, b, c]))

    collections, repack = unpack_collections(*args)
    assert len(collections) == 3

    # Replace collections with `'~a'` strings
    result = repack(['~a', '~b', '~c'])
    sol = build('~a', '~b', '~c', ['~a', '~b', '~c'])
    assert result == sol

    # traverse=False
    collections, repack = unpack_collections(*args, traverse=False)
    assert len(collections) == 2  # just a and b
    assert repack(collections) == args

    # No collections
    collections, repack = unpack_collections(1, 2, {'a': 3})
    assert not collections
    assert repack(collections) == (1, 2, {'a': 3})

    # Result that looks like a task
    def fail(*args):
        raise ValueError("Shouldn't have been called")

    collections, repack = unpack_collections(a, (fail, 1), [(fail, 2, 3)],
                                             traverse=False)
    repack(collections)  # Smoketest task literals
    repack([(fail, 1)])  # Smoketest results that look like tasks


class Tuple(DaskMethodsMixin):
    __slots__ = ('_dask', '_keys')
    __dask_scheduler__ = staticmethod(dask.threaded.get)

    def __init__(self, dsk, keys):
        self._dask = dsk
        self._keys = keys

    def __add__(self, other):
        if isinstance(other, Tuple):
            return Tuple(merge(self._dask, other._dask),
                         self._keys + other._keys)
        return NotImplemented

    def __dask_graph__(self):
        return self._dask

    def __dask_keys__(self):
        return self._keys

    def __dask_tokenize__(self):
        return self._keys

    def __dask_postcompute__(self):
        return tuple, ()

    def __dask_postpersist__(self):
        return Tuple, (self._keys,)


def test_custom_collection():
    dsk = {'a': 1, 'b': 2}
    dsk2 = {'c': (add, 'a', 'b'),
            'd': (add, 'c', 1)}
    dsk2.update(dsk)
    dsk3 = {'e': (add, 'a', 4),
            'f': (inc, 'e')}
    dsk3.update(dsk)

    x = Tuple(dsk, ['a', 'b'])
    y = Tuple(dsk2, ['c', 'd'])
    z = Tuple(dsk3, ['e', 'f'])

    # __slots__ defined on base mixin class propogates
    with pytest.raises(AttributeError):
        x.foo = 1

    # is_dask_collection
    assert is_dask_collection(x)

    # tokenize
    assert tokenize(x) == tokenize(x)
    assert tokenize(x) != tokenize(y)

    # compute
    assert x.compute() == (1, 2)
    assert dask.compute(x, [y, z]) == ((1, 2), [(3, 4), (5, 6)])
    t = x + y + z
    assert t.compute() == (1, 2, 3, 4, 5, 6)

    # persist
    t2 = t.persist()
    assert isinstance(t2, Tuple)
    assert t2._dask == dict(zip('abcdef', range(1, 7)))
    assert t2.compute() == (1, 2, 3, 4, 5, 6)
    x2, y2, z2 = dask.persist(x, y, z)
    t3 = x2 + y2 + z2
    assert t2._dask == t3._dask


@pytest.mark.skipif('not db')
def test_compute_no_opt():
    # Bag does `fuse` by default. Test that with `optimize_graph=False` that
    # doesn't get called. We check this by using a callback to track the keys
    # that are computed.
    from dask.callbacks import Callback
    b = db.from_sequence(range(100), npartitions=4)
    add1 = tz.partial(add, 1)
    mul2 = tz.partial(mul, 2)
    o = b.map(add1).map(mul2)
    # Check that with the kwarg, the optimization doesn't happen
    keys = []
    with Callback(pretask=lambda key, *args: keys.append(key)):
        o.compute(get=dask.get, optimize_graph=False)
    assert len([k for k in keys if '-mul-' in k[0]]) == 4
    assert len([k for k in keys if '-add-' in k[0]]) == 4
    # Check that without the kwarg, the optimization does happen
    keys = []
    with Callback(pretask=lambda key, *args: keys.append(key)):
        o.compute(get=dask.get)
    # Names of fused tasks have been merged, and the original key is an alias.
    # Otherwise, the lengths below would be 4 and 0.
    assert len([k for k in keys if '-mul-' in k[0]]) == 8
    assert len([k for k in keys if '-add-' in k[0]]) == 4
    assert len([k for k in keys if 'add-map-mul' in k[0]]) == 4  # See? Renamed


@pytest.mark.skipif('not da')
def test_compute_array():
    arr = np.arange(100).reshape((10, 10))
    darr = da.from_array(arr, chunks=(5, 5))
    darr1 = darr + 1
    darr2 = darr + 2
    out1, out2 = compute(darr1, darr2)
    assert np.allclose(out1, arr + 1)
    assert np.allclose(out2, arr + 2)


@pytest.mark.skipif('not da')
def test_persist_array():
    from dask.array.utils import assert_eq
    arr = np.arange(100).reshape((10, 10))
    x = da.from_array(arr, chunks=(5, 5))
    x = (x + 1) - x.mean(axis=0)
    y = x.persist()

    assert_eq(x, y)
    assert set(y.dask).issubset(x.dask)
    assert len(y.dask) == y.npartitions


@pytest.mark.skipif('not dd')
def test_compute_dataframe():
    df = pd.DataFrame({'a': [1, 2, 3, 4], 'b': [5, 5, 3, 3]})
    ddf = dd.from_pandas(df, npartitions=2)
    ddf1 = ddf.a + 1
    ddf2 = ddf.a + ddf.b
    out1, out2 = compute(ddf1, ddf2)
    pd.util.testing.assert_series_equal(out1, df.a + 1)
    pd.util.testing.assert_series_equal(out2, df.a + df.b)


@pytest.mark.skipif('not dd or not da')
def test_compute_array_dataframe():
    arr = np.arange(100).reshape((10, 10))
    darr = da.from_array(arr, chunks=(5, 5)) + 1
    df = pd.DataFrame({'a': [1, 2, 3, 4], 'b': [5, 5, 3, 3]})
    ddf = dd.from_pandas(df, npartitions=2).a + 2
    arr_out, df_out = compute(darr, ddf)
    assert np.allclose(arr_out, arr + 1)
    pd.util.testing.assert_series_equal(df_out, df.a + 2)


@pytest.mark.skipif('not da or not db')
def test_compute_array_bag():
    x = da.arange(5, chunks=2)
    b = db.from_sequence([1, 2, 3])

    pytest.raises(ValueError, lambda: compute(x, b))

    xx, bb = compute(x, b, get=dask.get)
    assert np.allclose(xx, np.arange(5))
    assert bb == [1, 2, 3]


@pytest.mark.skipif('not da')
def test_compute_with_literal():
    x = da.arange(5, chunks=2)
    y = 10

    xx, yy = compute(x, y)
    assert (xx == x.compute()).all()
    assert yy == y

    assert compute(5) == (5,)


def test_compute_nested():
    a = delayed(1) + 5
    b = a + 1
    c = a + 2
    assert (compute({'a': a, 'b': [1, 2, b]}, (c, 2)) ==
            ({'a': 6, 'b': [1, 2, 7]}, (8, 2)))

    res = compute([a, b], c, traverse=False)
    assert res[0][0] is a
    assert res[0][1] is b
    assert res[1] == 8


@pytest.mark.skipif('not da')
@pytest.mark.skipif(sys.flags.optimize,
                    reason="graphviz exception with Python -OO flag")
def test_visualize():
    pytest.importorskip('graphviz')
    with tmpdir() as d:
        x = da.arange(5, chunks=2)
        x.visualize(filename=os.path.join(d, 'mydask'))
        assert os.path.exists(os.path.join(d, 'mydask.png'))

        x.visualize(filename=os.path.join(d, 'mydask.pdf'))
        assert os.path.exists(os.path.join(d, 'mydask.pdf'))

        visualize(x, 1, 2, filename=os.path.join(d, 'mydask.png'))
        assert os.path.exists(os.path.join(d, 'mydask.png'))

        dsk = {'a': 1, 'b': (add, 'a', 2), 'c': (mul, 'a', 1)}
        visualize(x, dsk, filename=os.path.join(d, 'mydask.png'))
        assert os.path.exists(os.path.join(d, 'mydask.png'))

        x = Tuple(dsk, ['a', 'b', 'c'])
        visualize(x, filename=os.path.join(d, 'mydask.png'))
        assert os.path.exists(os.path.join(d, 'mydask.png'))


@pytest.mark.skipif('not da')
@pytest.mark.skipif(sys.flags.optimize,
                    reason="graphviz exception with Python -OO flag")
def test_visualize_order():
    pytest.importorskip('matplotlib')
    x = da.arange(5, chunks=2)
    with tmpfile(extension='dot') as fn:
        x.visualize(color='order', filename=fn, cmap='RdBu')
        with open(fn) as f:
            text = f.read()
        assert 'color="#' in text


def test_use_cloudpickle_to_tokenize_functions_in__main__():
    import sys
    from textwrap import dedent

    defn = dedent("""
    def inc():
        return x
    """)

    __main__ = sys.modules['__main__']
    exec(compile(defn, '<test>', 'exec'), __main__.__dict__)
    f = __main__.inc

    t = normalize_token(f)
    assert b'cloudpickle' in t


def inc_to_dec(dsk, keys):
    for key in dsk:
        if dsk[key][0] == inc:
            dsk[key] = (dec,) + dsk[key][1:]
    return dsk


def test_optimizations_keyword():
    x = dask.delayed(inc)(1)
    assert x.compute() == 2

    with dask.set_options(optimizations=[inc_to_dec]):
        assert x.compute() == 0

    assert x.compute() == 2


def test_optimize():
    x = dask.delayed(inc)(1)
    y = dask.delayed(inc)(x)
    z = x + y

    x2, y2, z2, constant = optimize(x, y, z, 1)
    assert constant == 1

    # Same graphs for each
    dsk = dict(x2.dask)
    assert dict(y2.dask) == dsk
    assert dict(z2.dask) == dsk

    # Computationally equivalent
    assert dask.compute(x2, y2, z2) == dask.compute(x, y, z)

    # Applying optimizations before compute and during compute gives
    # same results. Shows optimizations are occurring.
    sols = dask.compute(x, y, z, optimizations=[inc_to_dec])
    x3, y3, z3 = optimize(x, y, z, optimizations=[inc_to_dec])
    assert dask.compute(x3, y3, z3) == sols

    # Optimize respects global optimizations as well
    with dask.set_options(optimizations=[inc_to_dec]):
        x4, y4, z4 = optimize(x, y, z)
    for a, b in zip([x3, y3, z3], [x4, y4, z4]):
        assert dict(a.dask) == dict(b.dask)


def test_optimize_nested():
    a = dask.delayed(inc)(1)
    b = dask.delayed(inc)(a)
    c = a + b

    result = optimize({'a': a, 'b': [1, 2, b]}, (c, 2))

    a2 = result[0]['a']
    b2 = result[0]['b'][2]
    c2 = result[1][0]

    assert isinstance(a2, Delayed)
    assert isinstance(b2, Delayed)
    assert isinstance(c2, Delayed)
    assert dict(a2.dask) == dict(b2.dask) == dict(c2.dask)
    assert compute(*result) == ({'a': 2, 'b': [1, 2, 3]}, (5, 2))

    res = optimize([a, b], c, traverse=False)
    assert res[0][0] is a
    assert res[0][1] is b
    assert res[1].compute() == 5


# TODO: remove after deprecation cycle of `dask.optimize` module is completed
def test_optimize_has_deprecated_module_functions_as_attributes():
    import dask.optimize as deprecated_optimize
    # Function has method attributes
    assert dask.optimize.cull is deprecated_optimize.cull
    assert dask.optimize.inline is deprecated_optimize.inline
    with pytest.warns(UserWarning):
        dask.optimize.cull({}, [])


def test_default_imports():
    """
    Startup time: `import dask` should not import too many modules.
    """
    code = """if 1:
        import dask
        import sys

        print(sorted(sys.modules))
        """

    out = subprocess.check_output([sys.executable, '-c', code])
    modules = set(eval(out.decode()))
    assert 'dask' in modules
    blacklist = ['dask.array', 'dask.dataframe', 'numpy', 'pandas',
                 'partd', 's3fs', 'distributed']
    for mod in blacklist:
        assert mod not in modules


def test_persist_literals():
    assert persist(1, 2, 3) == (1, 2, 3)


def test_persist_nested():
    a = delayed(1) + 5
    b = a + 1
    c = a + 2
    result = persist({'a': a, 'b': [1, 2, b]}, (c, 2))
    assert isinstance(result[0]['a'], Delayed)
    assert isinstance(result[0]['b'][2], Delayed)
    assert isinstance(result[1][0], Delayed)
    assert compute(*result) == ({'a': 6, 'b': [1, 2, 7]}, (8, 2))

    res = persist([a, b], c, traverse=False)
    assert res[0][0] is a
    assert res[0][1] is b
    assert res[1].compute() == 8


def test_persist_delayed():
    x1 = delayed(1)
    x2 = delayed(inc)(x1)
    x3 = delayed(inc)(x2)

    xx, = persist(x3)
    assert isinstance(xx, Delayed)
    assert xx.key == x3.key
    assert len(xx.dask) == 1

    assert x3.compute() == xx.compute()


@pytest.mark.skipif('not da or not db')
def test_persist_array_bag():
    x = da.arange(5, chunks=2) + 1
    b = db.from_sequence([1, 2, 3]).map(inc)

    with pytest.raises(ValueError):
        persist(x, b)

    xx, bb = persist(x, b, get=dask.get)

    assert isinstance(xx, da.Array)
    assert isinstance(bb, db.Bag)

    assert xx.name == x.name
    assert bb.name == b.name
    assert len(xx.dask) == xx.npartitions < len(x.dask)
    assert len(bb.dask) == bb.npartitions < len(b.dask)

    assert np.allclose(x, xx)
    assert list(b) == list(bb)


def test_normalize_function_limited_size():
    for i in range(1000):
        normalize_function(lambda x: x)

    assert 50 < len(function_cache) < 600


def test_optimize_globals():
    da = pytest.importorskip('dask.array')
    db = pytest.importorskip('dask.bag')

    x = da.ones(10, chunks=(5,))

    def optimize_double(dsk, keys):
        return {k: (mul, 2, v) for k, v in dsk.items()}

    from dask.array.utils import assert_eq

    assert_eq(x + 1, np.ones(10) + 1)

    with dask.set_options(array_optimize=optimize_double):
        assert_eq(x + 1, (np.ones(10) * 2 + 1) * 2)

    assert_eq(x + 1, np.ones(10) + 1)

    b = db.range(10, npartitions=2)

    with dask.set_options(array_optimize=optimize_double):
        xx, bb = dask.compute(x + 1, b.map(inc), get=dask.get)
        assert_eq(xx, (np.ones(10) * 2 + 1) * 2)


def test_optimize_None():
    da = pytest.importorskip('dask.array')

    x = da.ones(10, chunks=(5,))
    y = x[:9][1:8][::2] + 1  # normally these slices would be fused

    def my_get(dsk, keys):
        assert dsk == dict(y.dask)  # but they aren't
        return dask.get(dsk, keys)

    with dask.set_options(array_optimize=None, get=my_get):
        y.compute()
