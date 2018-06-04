import pytest
pytest.importorskip('numpy')

import numpy as np
import dask
import dask.array as da
from dask.optimization import fuse
from dask.utils import SerializableLock
from dask.array.core import getter, getter_nofancy
from dask.array.optimization import (getitem, optimize, optimize_slices,
                                     fuse_slice)
from dask.array.utils import assert_eq


def test_fuse_getitem():
    pairs = [((getter, (getter, 'x', slice(1000, 2000)), slice(15, 20)),
              (getter, 'x', slice(1015, 1020))),

             ((getitem, (getter, 'x', (slice(1000, 2000), slice(100, 200))),
                        (slice(15, 20), slice(50, 60))),
              (getter, 'x', (slice(1015, 1020), slice(150, 160)))),

             ((getitem, (getter_nofancy, 'x', (slice(1000, 2000), slice(100, 200))),
                        (slice(15, 20), slice(50, 60))),
              (getter_nofancy, 'x', (slice(1015, 1020), slice(150, 160)))),

             ((getter, (getter, 'x', slice(1000, 2000)), 10),
              (getter, 'x', 1010)),

             ((getitem, (getter, 'x', (slice(1000, 2000), 10)),
                        (slice(15, 20),)),
              (getter, 'x', (slice(1015, 1020), 10))),

             ((getitem, (getter_nofancy, 'x', (slice(1000, 2000), 10)),
                        (slice(15, 20),)),
              (getter_nofancy, 'x', (slice(1015, 1020), 10))),

             ((getter, (getter, 'x', (10, slice(1000, 2000))),
               (slice(15, 20), )),
              (getter, 'x', (10, slice(1015, 1020)))),

             ((getter, (getter, 'x', (slice(1000, 2000), slice(100, 200))),
               (slice(None, None), slice(50, 60))),
              (getter, 'x', (slice(1000, 2000), slice(150, 160)))),

             ((getter, (getter, 'x', (None, slice(None, None))),
               (slice(None, None), 5)),
              (getter, 'x', (None, 5))),

             ((getter, (getter, 'x', (slice(1000, 2000), slice(10, 20))),
               (slice(5, 10),)),
              (getter, 'x', (slice(1005, 1010), slice(10, 20)))),

             ((getitem, (getitem, 'x', (slice(1000, 2000),)),
               (slice(5, 10), slice(10, 20))),
              (getitem, 'x', (slice(1005, 1010), slice(10, 20)))),

             ((getter, (getter, 'x', slice(1000, 2000), False, False), slice(15, 20)),
              (getter, 'x', slice(1015, 1020))),

             ((getter, (getter, 'x', slice(1000, 2000)), slice(15, 20), False, False),
              (getter, 'x', slice(1015, 1020))),

             ((getter, (getter_nofancy, 'x', slice(1000, 2000), False, False),
               slice(15, 20), False, False),
              (getter_nofancy, 'x', slice(1015, 1020), False, False)),
             ]

    for inp, expected in pairs:
        result = optimize_slices({'y': inp})
        assert result == {'y': expected}


def test_fuse_getitem_lock():
    lock1 = SerializableLock()
    lock2 = SerializableLock()

    pairs = [((getter, (getter, 'x', slice(1000, 2000), True, lock1), slice(15, 20)),
              (getter, 'x', slice(1015, 1020), True, lock1)),

             ((getitem, (getter, 'x', (slice(1000, 2000), slice(100, 200)), True, lock1),
                        (slice(15, 20), slice(50, 60))),
              (getter, 'x', (slice(1015, 1020), slice(150, 160)), True, lock1)),

             ((getitem, (getter_nofancy, 'x', (slice(1000, 2000), slice(100, 200)), True, lock1),
                        (slice(15, 20), slice(50, 60))),
              (getter_nofancy, 'x', (slice(1015, 1020), slice(150, 160)), True, lock1)),

             ((getter, (getter, 'x', slice(1000, 2000), True, lock1), slice(15, 20), True, lock2),
              (getter, (getter, 'x', slice(1000, 2000), True, lock1), slice(15, 20), True, lock2))]

    for inp, expected in pairs:
        result = optimize_slices({'y': inp})
        assert result == {'y': expected}


def test_optimize_with_getitem_fusion():
    dsk = {'a': 'some-array',
           'b': (getter, 'a', (slice(10, 20), slice(100, 200))),
           'c': (getter, 'b', (5, slice(50, 60)))}

    result = optimize(dsk, ['c'])
    expected_task = (getter, 'some-array', (15, slice(150, 160)))
    assert any(v == expected_task for v in result.values())
    assert len(result) < len(dsk)


def test_optimize_slicing():
    dsk = {'a': (range, 10),
           'b': (getter, 'a', (slice(None, None, None),)),
           'c': (getter, 'b', (slice(None, None, None),)),
           'd': (getter, 'c', (slice(0, 5, None),)),
           'e': (getter, 'd', (slice(None, None, None),))}

    expected = {'e': (getter, (range, 10), (slice(0, 5, None),))}
    result = optimize_slices(fuse(dsk, [], rename_keys=False)[0])
    assert result == expected

    # protect output keys
    expected = {'c': (getter, (range, 10), (slice(0, None, None),)),
                'd': (getter, 'c', (slice(0, 5, None),)),
                'e': (getter, 'd', (slice(None, None, None),))}
    result = optimize_slices(fuse(dsk, ['c', 'd', 'e'], rename_keys=False)[0])

    assert result == expected


def test_fuse_slice():
    assert fuse_slice(slice(10, 15), slice(0, 5, 2)) == slice(10, 15, 2)

    assert (fuse_slice((slice(100, 200),), (None, slice(10, 20))) ==
            (None, slice(110, 120)))
    assert (fuse_slice((slice(100, 200),), (slice(10, 20), None)) ==
            (slice(110, 120), None))
    assert (fuse_slice((1,), (None,)) ==
            (1, None))
    assert (fuse_slice((1, slice(10, 20)), (None, None, 3, None)) ==
            (1, None, None, 13, None))

    with pytest.raises(NotImplementedError):
        fuse_slice(slice(10, 15, 2), -1)
    # Regression test for #3076
    with pytest.raises(NotImplementedError):
        fuse_slice(None, np.array([0, 0]))


def test_fuse_slice_with_lists():
    assert fuse_slice(slice(10, 20, 2), [1, 2, 3]) == [12, 14, 16]
    assert fuse_slice([10, 20, 30, 40, 50], [3, 1, 2]) == [40, 20, 30]
    assert fuse_slice([10, 20, 30, 40, 50], 3) == 40
    assert fuse_slice([10, 20, 30, 40, 50], -1) == 50
    assert fuse_slice([10, 20, 30, 40, 50], slice(1, None, 2)) == [20, 40]
    assert fuse_slice((slice(None), slice(0, 10), [1, 2, 3]),
                      (slice(None), slice(1, 5), slice(None))) == (slice(0, None), slice(1, 5), [1, 2, 3])
    assert fuse_slice((slice(None), slice(None), [1, 2, 3]),
                      (slice(None), slice(1, 5), 1)) == (slice(0, None), slice(1, 5), 2)


def test_nonfusible_fancy_indexing():
    nil = slice(None)
    cases = [# x[:, list, :][int, :, :]
             ((nil, [1, 2, 3], nil), (0, nil, nil)),
             # x[int, :, :][:, list, :]
             ((0, nil, nil), (nil, [1, 2, 3], nil)),
             # x[:, list, :, :][:, :, :, int]
             ((nil, [1, 2], nil, nil), (nil, nil, nil, 0))]

    for a, b in cases:
        with pytest.raises(NotImplementedError):
            fuse_slice(a, b)


def test_hard_fuse_slice_cases():
    dsk = {'x': (getter, (getter, 'x', (None, slice(None, None))),
                 (slice(None, None), 5))}
    assert optimize_slices(dsk) == {'x': (getter, 'x', (None, 5))}


def test_dont_fuse_numpy_arrays():
    x = np.ones(10)
    for chunks in [(5,), (10,)]:
        y = da.from_array(x, chunks=(10,))

        dsk = y.__dask_optimize__(y.dask, y.__dask_keys__())
        assert sum(isinstance(v, np.ndarray) for v in dsk.values()) == 1


def test_minimize_data_transfer():
    x = np.ones(100)
    y = da.from_array(x, chunks=25)
    z = y + 1
    dsk = z.__dask_optimize__(z.dask, z.__dask_keys__())

    keys = list(dsk)
    results = dask.get(dsk, keys)
    big_key = [k for k, r in zip(keys, results) if r is x][0]
    dependencies, dependents = dask.core.get_deps(dsk)
    deps = dependents[big_key]

    assert len(deps) == 4
    for dep in deps:
        assert dsk[dep][0] in (getitem, getter)
        assert dsk[dep][1] == big_key


def test_fuse_slices_with_alias():
    dsk = {'x': np.arange(16).reshape((4, 4)),
           ('dx', 0, 0): (getter, 'x', (slice(0, 4), slice(0, 4))),
           ('alias', 0, 0): ('dx', 0, 0),
           ('dx2', 0): (getitem, ('alias', 0, 0), (slice(None), 0))}
    keys = [('dx2', 0)]
    dsk2 = optimize(dsk, keys)
    assert len(dsk2) == 3
    fused_key = set(dsk2).difference(['x', ('dx2', 0)]).pop()
    assert dsk2[fused_key] == (getter, 'x', (slice(0, 4), 0))


def test_dont_fuse_fancy_indexing_in_getter_nofancy():
    dsk = {'a': (getitem, (getter_nofancy, 'x', (slice(10, 20, None), slice(100, 200, None))),
                 ([1, 3], slice(50, 60, None)))}
    assert optimize_slices(dsk) == dsk

    dsk = {'a': (getitem, (getter_nofancy, 'x', [1, 2, 3]), 0)}
    assert optimize_slices(dsk) == dsk


@pytest.mark.parametrize('chunks', [10, 5, 3])
def test_fuse_getter_with_asarray(chunks):
    x = np.ones(10) * 1234567890
    y = da.ones(10, chunks=chunks)
    z = x + y
    dsk = z.__dask_optimize__(z.dask, z.__dask_keys__())
    assert any(v is x for v in dsk.values())
    for v in dsk.values():
        s = str(v)
        assert s.count('getitem') + s.count('getter') <= 1
        if v is not x:
            assert '1234567890' not in s
    n_getters = len([v for v in dsk.values() if v[0] in (getitem, getter)])
    if y.npartitions > 1:
        assert n_getters == y.npartitions
    else:
        assert n_getters == 0

    assert_eq(z, x + 1)


@pytest.mark.parametrize('get,remove',
                         [(getter, False), (getter_nofancy, False), (getitem, True)])
def test_remove_no_op_slices_if_get_is_not_getter_or_getter_nofancy(get, remove):
    # Test that no-op slices are removed as long as get is not getter or
    # getter_nofancy. This ensures that `get` calls are always made in all
    # tasks created by `from_array`, even after optimization
    null = slice(0,None)
    opts = [((get, 'x', null, False, False),
             'x' if remove else (get, 'x', null, False, False)),
            ((getitem, (get, 'x', null, False, False), null),
             'x' if remove else (get, 'x', null, False, False)),
            ((getitem, (get, 'x', (null, null), False, False), ()),
             'x' if remove else (get, 'x', (null, null), False, False))]
    for orig, final in opts:
        assert optimize_slices({'a': orig}) == {'a': final}


def test_turn_off_fusion():
    x = da.ones(10, chunks=(5,))
    y = da.sum(x + 1 + 2 + 3)

    a = y.__dask_optimize__(y.dask, y.__dask_keys__())

    with dask.set_options(fuse_ave_width=0):
        b = y.__dask_optimize__(y.dask, y.__dask_keys__())

    assert dask.get(a, y.__dask_keys__()) == dask.get(b, y.__dask_keys__())
    assert len(a) < len(b)
