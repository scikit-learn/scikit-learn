from itertools import product
import warnings

import pytest
np = pytest.importorskip('numpy')

import dask
from dask.utils import funcname
from dask.array.utils import assert_eq
from dask.array.rechunk import intersect_chunks, rechunk, normalize_chunks
from dask.array.rechunk import cumdims_label, _breakpoints, _intersect_1d, _old_to_new
from dask.array.rechunk import plan_rechunk, divide_to_width, merge_to_number
import dask.array as da


def test_rechunk_internals_1():
    """ Test the cumdims_label and _breakpoints and
    _intersect_1d internal funcs to rechunk."""
    new = cumdims_label(((1,1,2),(1,5,1)),'n')
    old = cumdims_label(((4, ),(1,) * 5),'o')
    breaks = tuple(_breakpoints(o, n) for o, n in zip(old, new))
    answer = (('o', 0), ('n', 0), ('n', 1), ('n', 2), ('o', 4), ('n', 4))
    assert breaks[0] == answer
    answer2 = (('o', 0), ('n', 0), ('o', 1), ('n', 1), ('o', 2),
               ('o', 3), ('o', 4), ('o', 5), ('n', 6), ('n', 7))
    assert breaks[1] == answer2
    i1d = [_intersect_1d(b) for b in breaks]
    answer3 = [[(0, slice(0, 1))],
               [(0, slice(1, 2))],
               [(0, slice(2, 4))]]
    assert i1d[0] == answer3
    answer4 = [[(0, slice(0, 1))],
               [(1, slice(0, 1)),
                (2, slice(0, 1)),
                (3, slice(0, 1)),
                (4, slice(0, 1)),
                (5, slice(0, 1))],
               [(5, slice(1, 2))]]
    assert i1d[1] == answer4


def test_intersect_1():
    """ Convert 1 D chunks"""
    old = ((10, 10, 10, 10, 10), )
    new = ((25, 5, 20), )
    answer = [(((0, slice(0, 10)), ),
              ((1, slice(0, 10)), ),
              ((2, slice(0, 5)), )),
              (((2, slice(5, 10)), ), ),
              (((3, slice(0, 10)), ), ((4, slice(0, 10)), ))
              ]
    cross = list(intersect_chunks(old_chunks=old, new_chunks=new))
    assert answer == cross


def test_intersect_2():
    """ Convert 1 D chunks"""
    old = ((20, 20, 20, 20, 20), )
    new = ((58, 4, 20, 18),)
    answer = [(((0, slice(0, 20)), ),
              ((1, slice(0, 20)), ),
              ((2, slice(0, 18)), )),
              (((2, slice(18, 20)), ), ((3, slice(0, 2)), )),
              (((3, slice(2, 20)), ), ((4, slice(0, 2)), )),
              (((4, slice(2, 20)), ), )
              ]
    cross = list(intersect_chunks(old_chunks=old, new_chunks=new))
    assert answer == cross


def test_rechunk_1d():
    """Try rechunking a random 1d matrix"""
    a = np.random.uniform(0, 1, 300)
    x = da.from_array(a, chunks=((100, ) * 3, ))
    new = ((50, ) * 6,)
    x2 = rechunk(x, chunks=new)
    assert x2.chunks == new
    assert np.all(x2.compute() == a)


def test_rechunk_2d():
    """Try rechunking a random 2d matrix"""
    a = np.random.uniform(0, 1, 300).reshape((10, 30))
    x = da.from_array(a, chunks=((1, 2, 3, 4), (5, ) * 6))
    new = ((5, 5), (15, ) * 2)
    x2 = rechunk(x, chunks=new)
    assert x2.chunks == new
    assert np.all(x2.compute() == a)


def test_rechunk_4d():
    """Try rechunking a random 4d matrix"""
    old = ((5, 5), ) * 4
    a = np.random.uniform(0, 1, 10000).reshape((10, ) * 4)
    x = da.from_array(a, chunks=old)
    new = ((10, ), ) * 4
    x2 = rechunk(x, chunks=new)
    assert x2.chunks == new
    assert np.all(x2.compute() == a)


def test_rechunk_expand():
    a = np.random.uniform(0, 1, 100).reshape((10, 10))
    x = da.from_array(a, chunks=(5, 5))
    y = x.rechunk(chunks=((3, 3, 3, 1), (3, 3, 3, 1)))
    assert np.all(y.compute() == a)


def test_rechunk_expand2():
    (a, b) = (3, 2)
    orig = np.random.uniform(0, 1, a ** b).reshape((a,) * b)
    for off, off2 in product(range(1, a - 1), range(1, a - 1)):
        old = ((a - off, off), ) * b
        x = da.from_array(orig, chunks=old)
        new = ((a - off2, off2), ) * b
        assert np.all(x.rechunk(chunks=new).compute() == orig)
        if a - off - off2 > 0:
            new = ((off, a - off2 - off, off2), ) * b
            y = x.rechunk(chunks=new).compute()
            assert np.all(y == orig)


def test_rechunk_method():
    """ Test rechunking can be done as a method of dask array."""
    old = ((5, 2, 3), ) * 4
    new = ((3, 3, 3, 1), ) * 4
    a = np.random.uniform(0, 1, 10000).reshape((10, ) * 4)
    x = da.from_array(a, chunks=old)
    x2 = x.rechunk(chunks=new)
    assert x2.chunks == new
    assert np.all(x2.compute() == a)


def test_rechunk_blockshape():
    """ Test that blockshape can be used."""
    new_shape, new_chunks = (10, 10), (4, 3)
    new_blockdims = normalize_chunks(new_chunks, new_shape)
    old_chunks = ((4, 4, 2), (3, 3, 3, 1))
    a = np.random.uniform(0,1,100).reshape((10, 10))
    x = da.from_array(a, chunks=old_chunks)
    check1 = rechunk(x, chunks=new_chunks)
    assert check1.chunks == new_blockdims
    assert np.all(check1.compute() == a)


def test_dtype():
    x = da.ones(5, chunks=(2,))
    assert x.rechunk(chunks=(1,)).dtype == x.dtype


def test_rechunk_with_dict():
    x = da.ones((24, 24), chunks=(4, 8))
    y = x.rechunk(chunks={0: 12})
    assert y.chunks == ((12, 12), (8, 8, 8))

    x = da.ones((24, 24), chunks=(4, 8))
    y = x.rechunk(chunks={0: (12, 12)})
    assert y.chunks == ((12, 12), (8, 8, 8))

    x = da.ones((24, 24), chunks=(4, 8))
    y = x.rechunk(chunks={0: -1})
    assert y.chunks == ((24,), (8, 8, 8))


def test_rechunk_with_empty_input():
    x = da.ones((24, 24), chunks=(4, 8))
    assert x.rechunk(chunks={}).chunks == x.chunks
    pytest.raises(ValueError, lambda: x.rechunk(chunks=()))


def test_rechunk_with_null_dimensions():
    x = da.from_array(np.ones((24, 24)), chunks=(4, 8))
    assert (x.rechunk(chunks=(None, 4)).chunks ==
            da.ones((24, 24), chunks=(4, 4)).chunks)


def test_rechunk_with_integer():
    x = da.from_array(np.arange(5), chunks=4)
    y = x.rechunk(3)
    assert y.chunks == ((3, 2),)
    assert (x.compute() == y.compute()).all()


def test_rechunk_0d():
    a = np.array(42)
    x = da.from_array(a, chunks=())
    y = x.rechunk(())
    assert y.chunks == ()
    assert y.compute() == a


def test_rechunk_empty():
    x = da.ones((0, 10), chunks=(5, 5))
    y = x.rechunk((2, 2))
    assert y.chunks == ((0,), (2,) * 5)
    assert_eq(x, y)


def test_rechunk_same():
    x = da.ones((24, 24), chunks=(4, 8))
    y = x.rechunk(x.chunks)
    assert x is y


def test_rechunk_minus_one():
    x = da.ones((24, 24), chunks=(4, 8))
    y = x.rechunk((-1, 8))
    assert y.chunks == ((24,), (8, 8, 8))
    assert_eq(x, y)


def test_rechunk_intermediates():
    x = da.random.normal(10, 0.1, (10, 10), chunks=(10, 1))
    y = x.rechunk((1, 10))
    assert len(y.dask) > 30


def test_divide_to_width():
    chunks = divide_to_width((8, 9, 10), 10)
    assert chunks == (8, 9, 10)
    chunks = divide_to_width((8, 2, 9, 10, 11, 12), 4)
    # Note how 9 gives (3, 3, 3), not (4, 4, 1) or whatever
    assert chunks == (4, 4,
                      2,
                      3, 3, 3,
                      3, 3, 4,
                      3, 4, 4,
                      4, 4, 4,
                      )


def test_merge_to_number():
    chunks = merge_to_number((10,) * 4, 5)
    assert chunks == (10, 10, 10, 10)
    chunks = merge_to_number((10,) * 4, 4)
    assert chunks == (10, 10, 10, 10)
    chunks = merge_to_number((10,) * 4, 3)
    assert chunks == (20, 10, 10)
    chunks = merge_to_number((10,) * 4, 2)
    assert chunks == (20, 20)
    chunks = merge_to_number((10,) * 4, 1)
    assert chunks == (40,)

    chunks = merge_to_number((10,) * 10, 2)
    assert chunks == (50,) * 2
    chunks = merge_to_number((10,) * 10, 3)
    assert chunks == (40, 30, 30)

    chunks = merge_to_number((5, 1, 1, 15, 10), 4)
    assert chunks == (5, 2, 15, 10)
    chunks = merge_to_number((5, 1, 1, 15, 10), 3)
    assert chunks == (7, 15, 10)
    chunks = merge_to_number((5, 1, 1, 15, 10), 2)
    assert chunks == (22, 10)
    chunks = merge_to_number((5, 1, 1, 15, 10), 1)
    assert chunks == (32,)

    chunks = merge_to_number((1, 1, 1, 1, 3, 1, 1), 6)
    assert chunks == (2, 1, 1, 3, 1, 1)
    chunks = merge_to_number((1, 1, 1, 1, 3, 1, 1), 5)
    assert chunks == (2, 2, 3, 1, 1)
    chunks = merge_to_number((1, 1, 1, 1, 3, 1, 1), 4)
    assert chunks == (2, 2, 3, 2)
    chunks = merge_to_number((1, 1, 1, 1, 3, 1, 1), 3)
    assert chunks == (4, 3, 2)
    chunks = merge_to_number((1, 1, 1, 1, 3, 1, 1), 2)
    assert chunks == (4, 5)
    chunks = merge_to_number((1, 1, 1, 1, 3, 1, 1), 1)
    assert chunks == (9,)


def _plan(old_chunks, new_chunks, itemsize=1, block_size_limit=1e7):
    return plan_rechunk(old_chunks, new_chunks,
                        itemsize=itemsize,
                        block_size_limit=block_size_limit)


def _assert_steps(steps, expected):
    assert len(steps) == len(expected)
    assert steps == expected


def test_plan_rechunk():
    c = ((20,) * 2)              # coarse
    f = ((2,) * 20)              # fine
    nc = ((float('nan'),) * 2)   # nan-coarse
    nf = ((float('nan'),) * 20)  # nan-fine

    # Trivial cases
    steps = _plan((), ())
    _assert_steps(steps, [()])
    steps = _plan((c, ()), (f, ()))
    _assert_steps(steps, [(f, ())])

    # No intermediate required
    steps = _plan((c,), (f,))
    _assert_steps(steps, [(f,)])
    steps = _plan((f,), (c,))
    _assert_steps(steps, [(c,)])
    steps = _plan((c, c), (f, f))
    _assert_steps(steps, [(f, f)])
    steps = _plan((f, f), (c, c))
    _assert_steps(steps, [(c, c)])
    steps = _plan((f, c), (c, c))
    _assert_steps(steps, [(c, c)])

    # An intermediate is used to reduce graph size
    steps = _plan((f, c), (c, f))
    _assert_steps(steps, [(c, c), (c, f)])

    steps = _plan((c + c, c + f), (f + f, c + c))
    _assert_steps(steps, [(c + c, c + c), (f + f, c + c)])

    # Same, with unknown dim
    steps = _plan((nc + nf, c + c, c + f), (nc + nf, f + f, c + c))
    _assert_steps(steps, steps)

    # Just at the memory limit => an intermediate is used
    steps = _plan((f, c), (c, f), block_size_limit=400)
    _assert_steps(steps, [(c, c), (c, f)])

    # Hitting the memory limit => partial merge
    m = ((10,) * 4)   # mid

    steps = _plan((f, c), (c, f), block_size_limit=399)
    _assert_steps(steps, [(m, c), (c, f)])

    steps2 = _plan((f, c), (c, f), block_size_limit=3999, itemsize=10)
    _assert_steps(steps2, steps)

    # Larger problem size => more intermediates
    c = ((1000,) * 2)   # coarse
    f = ((2,) * 1000)   # fine

    steps = _plan((f, c), (c, f), block_size_limit=99999)
    assert len(steps) == 3
    assert steps[-1] == (c, f)
    for i in range(len(steps) - 1):
        prev = steps[i]
        succ = steps[i + 1]
        # Merging on the first dim, splitting on the second dim
        assert len(succ[0]) <= len(prev[0]) / 2.0
        assert len(succ[1]) >= len(prev[1]) * 2.0


def test_plan_rechunk_5d():
    # 5d problem
    c = ((10,) * 1)    # coarse
    f = ((1,) * 10)    # fine

    steps = _plan((c, c, c, c, c), (f, f, f, f, f))
    _assert_steps(steps, [(f, f, f, f, f)])
    steps = _plan((f, f, f, f, c), (c, c, c, f, f))
    _assert_steps(steps, [(c, c, c, f, c), (c, c, c, f, f)])
    # Only 1 dim can be merged at first
    steps = _plan((c, c, f, f, c), (c, c, c, f, f), block_size_limit=2e4)
    _assert_steps(steps, [(c, c, c, f, c), (c, c, c, f, f)])


def test_plan_rechunk_heterogenous():
    c = ((10,) * 1)    # coarse
    f = ((1,) * 10)    # fine
    cf = c + f
    cc = c + c
    ff = f + f
    fc = f + c

    # No intermediate required
    steps = _plan((cc, cf), (ff, ff))
    _assert_steps(steps, [(ff, ff)])
    steps = _plan((cf, fc), (ff, cf))
    _assert_steps(steps, [(ff, cf)])

    # An intermediate is used to reduce graph size
    steps = _plan((cc, cf), (ff, cc))
    _assert_steps(steps, [(cc, cc), (ff, cc)])

    steps = _plan((cc, cf, cc), (ff, cc, cf))
    _assert_steps(steps, [(cc, cc, cc), (ff, cc, cf)])

    # Imposing a memory limit => the first intermediate is constrained:
    #  * cc -> ff would increase the graph size: no
    #  * ff -> cf would increase the block size too much: no
    #  * cf -> cc fits the bill (graph size /= 10, block size neutral)
    #  * cf -> fc also fits the bill (graph size and block size neutral)
    steps = _plan((cc, ff, cf), (ff, cf, cc), block_size_limit=100)
    _assert_steps(steps, [(cc, ff, cc), (ff, cf, cc)])


def test_plan_rechunk_asymmetric():
    a = ((1,) * 1000, (80000000,))
    b = ((1000,), (80000,) * 1000)
    steps = plan_rechunk(a, b, itemsize=8)
    assert len(steps) > 1

    x = da.ones((1000, 80000000), chunks=(1, 80000000))
    y = x.rechunk((1000, x.shape[1] // 1000))
    assert len(y.dask) < 100000


def test_rechunk_warning():
    N = 20
    x = da.random.normal(size=(N, N, 100), chunks=(1, N, 100))
    with warnings.catch_warnings(record=True) as w:
        x = x.rechunk((N, 1, 100))

    assert not w


@pytest.mark.parametrize('shape,chunks', [[(4,), (2,)],
                                          [(4, 4), (2, 2)],
                                          [(4, 4), (4, 2)]])
def test_dont_concatenate_single_chunks(shape, chunks):
    x = da.ones(shape, chunks=shape)
    y = x.rechunk(chunks)
    dsk = dict(y.dask)
    assert not any(funcname(task[0]).startswith('concat')
                   for task in dsk.values()
                   if dask.istask(task))


def test_intersect_nan():
    old_chunks = ((float('nan'), float('nan')), (8,))
    new_chunks = ((float('nan'), float('nan')), (4, 4))

    result = list(intersect_chunks(old_chunks, new_chunks))
    expected = [
        (((0, slice(0, None, None)), (0, slice(0, 4, None))),),
        (((0, slice(0, None, None)), (0, slice(4, 8, None))),),
        (((1, slice(0, None, None)), (0, slice(0, 4, None))),),
        (((1, slice(0, None, None)), (0, slice(4, 8, None))),)
    ]
    assert result == expected


def test_intersect_nan_single():
    old_chunks = ((float('nan'),), (10,))
    new_chunks = ((float('nan'),), (5, 5))

    result = list(intersect_chunks(old_chunks, new_chunks))
    expected = [(((0, slice(0, None, None)), (0, slice(0, 5, None))),),
                (((0, slice(0, None, None)), (0, slice(5, 10, None))),)]
    assert result == expected


def test_intersect_nan_long():

    old_chunks = (tuple([float('nan')] * 4), (10,))
    new_chunks = (tuple([float('nan')] * 4), (5, 5))
    result = list(intersect_chunks(old_chunks, new_chunks))
    expected = [
        (((0, slice(0, None, None)), (0, slice(0, 5, None))),),
        (((0, slice(0, None, None)), (0, slice(5, 10, None))),),
        (((1, slice(0, None, None)), (0, slice(0, 5, None))),),
        (((1, slice(0, None, None)), (0, slice(5, 10, None))),),
        (((2, slice(0, None, None)), (0, slice(0, 5, None))),),
        (((2, slice(0, None, None)), (0, slice(5, 10, None))),),
        (((3, slice(0, None, None)), (0, slice(0, 5, None))),),
        (((3, slice(0, None, None)), (0, slice(5, 10, None))),)
    ]
    assert result == expected


def test_rechunk_unknown_from_pandas():
    dd = pytest.importorskip('dask.dataframe')
    pd = pytest.importorskip('pandas')

    arr = np.random.randn(50, 10)
    x = dd.from_pandas(pd.DataFrame(arr), 2).values
    result = x.rechunk((None, (5, 5)))
    assert np.isnan(x.chunks[0]).all()
    assert np.isnan(result.chunks[0]).all()
    assert result.chunks[1] == (5, 5)
    expected = da.from_array(arr, chunks=((25, 25), (10,))).rechunk((None, (5, 5)))
    assert_eq(result, expected)


def test_rechunk_unknown_from_array():
    dd = pytest.importorskip('dask.dataframe')
    # pd = pytest.importorskip('pandas')
    x = dd.from_array(da.ones(shape=(4, 4), chunks=(2, 2))).values
    # result = x.rechunk({1: 5})
    result = x.rechunk((None, 4))
    assert np.isnan(x.chunks[0]).all()
    assert np.isnan(result.chunks[0]).all()
    assert x.chunks[1] == (4,)
    assert_eq(x, result)


@pytest.mark.parametrize('x, chunks', [
    (da.ones(shape=(50, 10), chunks=(25, 10)), (None, 5)),
    (da.ones(shape=(50, 10), chunks=(25, 10)), {1: 5}),
    (da.ones(shape=(50, 10), chunks=(25, 10)), (None, (5, 5))),

    (da.ones(shape=(1000, 10), chunks=(5, 10)), (None, 5)),
    (da.ones(shape=(1000, 10), chunks=(5, 10)), {1: 5}),
    (da.ones(shape=(1000, 10), chunks=(5, 10)), (None, (5, 5))),

    (da.ones(shape=(10, 10), chunks=(10, 10)), (None, 5)),
    (da.ones(shape=(10, 10), chunks=(10, 10)), {1: 5}),
    (da.ones(shape=(10, 10), chunks=(10, 10)), (None, (5, 5))),

    (da.ones(shape=(10, 10), chunks=(10, 2)), (None, 5)),
    (da.ones(shape=(10, 10), chunks=(10, 2)), {1: 5}),
    (da.ones(shape=(10, 10), chunks=(10, 2)), (None, (5, 5))),
])
def test_rechunk_unknown(x, chunks):
    dd = pytest.importorskip('dask.dataframe')
    y = dd.from_array(x).values
    result = y.rechunk(chunks)
    expected = x.rechunk(chunks)

    assert_chunks_match(result.chunks, expected.chunks)
    assert_eq(result, expected)


def test_rechunk_unknown_explicit():
    dd = pytest.importorskip('dask.dataframe')
    x = da.ones(shape=(10, 10), chunks=(5, 2))
    y = dd.from_array(x).values
    result = y.rechunk(((float('nan'), float('nan')), (5, 5)))
    expected = x.rechunk((None, (5, 5)))
    assert_chunks_match(result.chunks, expected.chunks)
    assert_eq(result, expected)


def assert_chunks_match(left, right):
    for x, y in zip(left, right):
        if np.isnan(x).any():
            assert np.isnan(x).all()
        else:
            assert x == y


def test_rechunk_unknown_raises():
    dd = pytest.importorskip('dask.dataframe')

    x = dd.from_array(da.ones(shape=(10, 10), chunks=(5, 5))).values
    with pytest.raises(ValueError):
        x.rechunk((None, (5, 5, 5)))


def test_old_to_new_single():
    old = ((float('nan'), float('nan')), (8,))
    new = ((float('nan'), float('nan')), (4, 4))
    result = _old_to_new(old, new)

    expected = [[[(0, slice(0, None, None))], [(1, slice(0, None, None))]],
                [[(0, slice(0, 4, None))], [(0, slice(4, 8, None))]]]

    assert result == expected


def test_old_to_new():
    old = ((float('nan'),), (10,))
    new = ((float('nan'),), (5, 5))
    result = _old_to_new(old, new)
    expected = [[[(0, slice(0, None, None))]],
                [[(0, slice(0, 5, None))], [(0, slice(5, 10, None))]]]

    assert result == expected


def test_old_to_new_large():
    old = (tuple([float('nan')] * 4), (10,))
    new = (tuple([float('nan')] * 4), (5, 5))

    result = _old_to_new(old, new)
    expected = [[[(0, slice(0, None, None))],
                [(1, slice(0, None, None))],
                [(2, slice(0, None, None))],
                [(3, slice(0, None, None))]],
                [[(0, slice(0, 5, None))], [(0, slice(5, 10, None))]]]
    assert result == expected


def test_changing_raises():
    nan = float('nan')
    with pytest.raises(ValueError) as record:
        _old_to_new(((nan, nan), (4, 4)), ((nan, nan, nan), (4, 4)))

    assert 'unchanging' in str(record.value)


def test_old_to_new_known():
    old = ((10, 10, 10, 10, 10), )
    new = ((25, 5, 20), )
    result = _old_to_new(old, new)
    expected = [[[(0, slice(0, 10, None)), (1, slice(0, 10, None)), (2, slice(0, 5, None))],
                [(2, slice(5, 10, None))],
                [(3, slice(0, 10, None)), (4, slice(0, 10, None))]]]
    assert result == expected


def test_rechunk_zero_dim():
    da = pytest.importorskip('dask.array')

    x = da.ones((0, 10, 100), chunks=(0, 10, 10)).rechunk((0, 10, 50))
    assert len(x.compute()) == 0
