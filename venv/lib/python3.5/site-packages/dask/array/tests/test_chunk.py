from __future__ import absolute_import, division, print_function

import pytest
pytest.importorskip('numpy')

import numpy as np

from dask.array.chunk import coarsen, keepdims_wrapper
import dask.array as da


def test_keepdims_wrapper_no_axis():
    def summer(a, axis=None):
        return a.sum(axis=axis)

    summer_wrapped = keepdims_wrapper(summer)

    assert summer_wrapped != summer

    a = np.arange(24).reshape(1, 2, 3, 4)

    r = summer(a)
    rw = summer_wrapped(a, keepdims=True)
    rwf = summer_wrapped(a, keepdims=False)

    assert r.ndim == 0
    assert r.shape == tuple()
    assert r == 276

    assert rw.ndim == 4
    assert rw.shape == (1, 1, 1, 1)
    assert (rw == 276).all()

    assert rwf.ndim == 0
    assert rwf.shape == tuple()
    assert rwf == 276


def test_keepdims_wrapper_one_axis():
    def summer(a, axis=None):
        return a.sum(axis=axis)

    summer_wrapped = keepdims_wrapper(summer)

    assert summer_wrapped != summer

    a = np.arange(24).reshape(1, 2, 3, 4)

    r = summer(a, axis=2)
    rw = summer_wrapped(a, axis=2, keepdims=True)
    rwf = summer_wrapped(a, axis=2, keepdims=False)

    assert r.ndim == 3
    assert r.shape == (1, 2, 4)
    assert (r == np.array([[[12, 15, 18, 21], [48, 51, 54, 57]]])).all()

    assert rw.ndim == 4
    assert rw.shape == (1, 2, 1, 4)
    assert (rw == np.array([[[[12, 15, 18, 21]], [[48, 51, 54, 57]]]])).all()

    assert rwf.ndim == 3
    assert rwf.shape == (1, 2, 4)
    assert (rwf == np.array([[[12, 15, 18, 21], [48, 51, 54, 57]]])).all()


def test_keepdims_wrapper_two_axes():
    def summer(a, axis=None):
        return a.sum(axis=axis)

    summer_wrapped = keepdims_wrapper(summer)

    assert summer_wrapped != summer

    a = np.arange(24).reshape(1, 2, 3, 4)

    r = summer(a, axis=(1, 3))
    rw = summer_wrapped(a, axis=(1, 3), keepdims=True)
    rwf = summer_wrapped(a, axis=(1, 3), keepdims=False)

    assert r.ndim == 2
    assert r.shape == (1, 3)
    assert (r == np.array([[60, 92, 124]])).all()

    assert rw.ndim == 4
    assert rw.shape == (1, 1, 3, 1)
    assert (rw == np.array([[[[60], [92], [124]]]])).all()

    assert rwf.ndim == 2
    assert rwf.shape == (1, 3)
    assert (rwf == np.array([[60, 92, 124]])).all()


def test_coarsen():
    x = np.random.randint(10, size=(24, 24))
    y = coarsen(np.sum, x, {0: 2, 1: 4})
    assert y.shape == (12, 6)
    assert y[0, 0] == np.sum(x[:2, :4])


"""
def test_coarsen_on_uneven_shape():
    x = np.random.randint(10, size=(23, 24))
    y = coarsen(np.sum, x, {0: 2, 1: 4})
    assert y.shape == (12, 6)
    assert y[0, 0] == np.sum(x[:2, :4])
    assert eq(y[11, :], x[23, :])
"""


def test_integer_input():
    assert da.zeros((4, 6), chunks=2).rechunk(3).chunks == ((3, 1), (3, 3))
