import pytest
pytest.importorskip('numpy')

from dask.array.wrap import ones
import dask.array as da
import numpy as np
import dask


def test_ones():
    a = ones((10, 10), dtype='i4', chunks=(4, 4))
    x = np.array(a)
    assert (x == np.ones((10, 10), 'i4')).all()


def test_size_as_list():
    a = ones([10, 10], dtype='i4', chunks=(4, 4))
    x = np.array(a)
    assert (x == np.ones((10, 10), dtype='i4')).all()


def test_singleton_size():
    a = ones(10, dtype='i4', chunks=(4,))
    x = np.array(a)
    assert (x == np.ones(10, dtype='i4')).all()


def test_kwargs():
    a = ones(10, dtype='i4', chunks=(4,))
    x = np.array(a)
    assert (x == np.ones(10, dtype='i4')).all()


def test_full():
    a = da.full((3, 3), 100, chunks=(2, 2), dtype='i8')

    assert (a.compute() == 100).all()
    assert a.dtype == a.compute(get=dask.get).dtype == 'i8'


def test_can_make_really_big_array_of_ones():
    ones((1000000, 1000000), chunks=(100000, 100000))
    ones(shape=(1000000, 1000000), chunks=(100000, 100000))


def test_wrap_consistent_names():
    assert (sorted(ones(10, dtype='i4', chunks=(4,)).dask) ==
            sorted(ones(10, dtype='i4', chunks=(4,)).dask))
    assert (sorted(ones(10, dtype='i4', chunks=(4,)).dask) !=
            sorted(ones(10, chunks=(4,)).dask))
    assert (sorted(da.full((3, 3), 100, chunks=(2, 2), dtype='f8').dask) ==
            sorted(da.full((3, 3), 100, chunks=(2, 2), dtype='f8').dask))
    assert (sorted(da.full((3, 3), 100, chunks=(2, 2), dtype='f8').dask) !=
            sorted(da.full((3, 3), 100, chunks=(2, 2)).dask))
