from itertools import combinations_with_replacement

import numpy as np

import pytest

import dask.array as da
import dask.array.fft
from dask.array.fft import fft_wrap
from dask.array.utils import assert_eq, same_keys

from dask.array.core import (
    normalize_chunks as _normalize_chunks,
)


all_1d_funcnames = [
    "fft",
    "ifft",
    "rfft",
    "irfft",
    "hfft",
    "ihfft",
]

all_nd_funcnames = [
    "fft2",
    "ifft2",
    "fftn",
    "ifftn",
    "rfft2",
    "irfft2",
    "rfftn",
    "irfftn",
]

nparr = np.arange(100).reshape(10, 10)
darr = da.from_array(nparr, chunks=(1, 10))
darr2 = da.from_array(nparr, chunks=(10, 1))
darr3 = da.from_array(nparr, chunks=(10, 10))


@pytest.mark.parametrize("funcname", all_1d_funcnames)
def test_cant_fft_chunked_axis(funcname):
    da_fft = getattr(da.fft, funcname)

    bad_darr = da.from_array(nparr, chunks=(5, 5))
    for i in range(bad_darr.ndim):
        with pytest.raises(ValueError):
            da_fft(bad_darr, axis=i)


@pytest.mark.parametrize("funcname", all_1d_funcnames)
def test_fft(funcname):
    da_fft = getattr(da.fft, funcname)
    np_fft = getattr(np.fft, funcname)

    assert_eq(da_fft(darr),
              np_fft(nparr))


@pytest.mark.parametrize("funcname", all_nd_funcnames)
def test_fft2n_shapes(funcname):
    da_fft = getattr(dask.array.fft, funcname)
    np_fft = getattr(np.fft, funcname)
    assert_eq(da_fft(darr3),
              np_fft(nparr))
    assert_eq(da_fft(darr3, (8, 9)),
              np_fft(nparr, (8, 9)))
    assert_eq(da_fft(darr3, (8, 9), axes=(1, 0)),
              np_fft(nparr, (8, 9), axes=(1, 0)))
    assert_eq(da_fft(darr3, (12, 11), axes=(1, 0)),
              np_fft(nparr, (12, 11), axes=(1, 0)))


@pytest.mark.parametrize("funcname", all_1d_funcnames)
def test_fft_n_kwarg(funcname):
    da_fft = getattr(da.fft, funcname)
    np_fft = getattr(np.fft, funcname)

    assert_eq(da_fft(darr, 5),
              np_fft(nparr, 5))
    assert_eq(da_fft(darr, 13),
              np_fft(nparr, 13))
    assert_eq(da_fft(darr2, axis=0),
              np_fft(nparr, axis=0))
    assert_eq(da_fft(darr2, 5, axis=0),
              np_fft(nparr, 5, axis=0))
    assert_eq(da_fft(darr2, 13, axis=0),
              np_fft(nparr, 13, axis=0))
    assert_eq(da_fft(darr2, 12, axis=0),
              np_fft(nparr, 12, axis=0))


@pytest.mark.parametrize("funcname", all_1d_funcnames)
def test_fft_consistent_names(funcname):
    da_fft = getattr(da.fft, funcname)

    assert same_keys(da_fft(darr, 5), da_fft(darr, 5))
    assert same_keys(da_fft(darr2, 5, axis=0), da_fft(darr2, 5, axis=0))
    assert not same_keys(da_fft(darr, 5), da_fft(darr, 13))


def test_wrap_bad_kind():
    with pytest.raises(ValueError):
        fft_wrap(np.ones)


@pytest.mark.parametrize("funcname", all_nd_funcnames)
@pytest.mark.parametrize("dtype", ["float32", "float64"])
def test_nd_ffts_axes(funcname, dtype):
    np_fft = getattr(np.fft, funcname)
    da_fft = getattr(da.fft, funcname)

    shape = (7, 8, 9)
    chunk_size = (3, 3, 3)
    a = np.arange(np.prod(shape), dtype=dtype).reshape(shape)
    d = da.from_array(a, chunks=chunk_size)

    for num_axes in range(1, d.ndim):
        for axes in combinations_with_replacement(range(d.ndim), num_axes):
            cs = list(chunk_size)
            for i in axes:
                cs[i] = shape[i]
            d2 = d.rechunk(cs)
            if len(set(axes)) < len(axes):
                with pytest.raises(ValueError):
                    da_fft(d2, axes=axes)
            else:
                r = da_fft(d2, axes=axes)
                er = np_fft(a, axes=axes)
                assert r.dtype == er.dtype
                assert r.shape == er.shape
                assert_eq(r, er)


@pytest.mark.parametrize("modname", ["numpy.fft", "scipy.fftpack"])
@pytest.mark.parametrize("funcname", all_1d_funcnames)
@pytest.mark.parametrize("dtype", ["float32", "float64"])
def test_wrap_ffts(modname, funcname, dtype):
    fft_mod = pytest.importorskip(modname)
    try:
        func = getattr(fft_mod, funcname)
    except AttributeError:
        pytest.skip("`%s` missing function `%s`." % (modname, funcname))

    darrc = darr.astype(dtype)
    darr2c = darr2.astype(dtype)
    nparrc = nparr.astype(dtype)

    if modname == "scipy.fftpack" and "rfft" in funcname:
        with pytest.raises(ValueError):
            fft_wrap(func)
    else:
        wfunc = fft_wrap(func)
        assert wfunc(darrc).dtype == func(nparrc).dtype
        assert wfunc(darrc).shape == func(nparrc).shape
        assert_eq(wfunc(darrc), func(nparrc))
        assert_eq(wfunc(darrc, axis=1), func(nparrc, axis=1))
        assert_eq(wfunc(darr2c, axis=0), func(nparrc, axis=0))
        assert_eq(wfunc(darrc, n=len(darrc) - 1),
                  func(nparrc, n=len(darrc) - 1))
        assert_eq(wfunc(darrc, axis=1, n=darrc.shape[1] - 1),
                  func(nparrc, n=darrc.shape[1] - 1))
        assert_eq(wfunc(darr2c, axis=0, n=darr2c.shape[0] - 1),
                  func(nparrc, axis=0, n=darr2c.shape[0] - 1))


@pytest.mark.parametrize("modname", ["numpy.fft", "scipy.fftpack"])
@pytest.mark.parametrize("funcname", all_nd_funcnames)
@pytest.mark.parametrize("dtype", ["float32", "float64"])
def test_wrap_fftns(modname, funcname, dtype):
    fft_mod = pytest.importorskip(modname)
    try:
        func = getattr(fft_mod, funcname)
    except AttributeError:
        pytest.skip("`%s` missing function `%s`." % (modname, funcname))

    darrc = darr.astype(dtype).rechunk(darr.shape)
    darr2c = darr2.astype(dtype).rechunk(darr2.shape)
    nparrc = nparr.astype(dtype)

    wfunc = fft_wrap(func)
    assert wfunc(darrc).dtype == func(nparrc).dtype
    assert wfunc(darrc).shape == func(nparrc).shape
    assert_eq(wfunc(darrc), func(nparrc))
    assert_eq(wfunc(darrc, axes=(1, 0)), func(nparrc, axes=(1, 0)))
    assert_eq(wfunc(darr2c, axes=(0, 1)), func(nparrc, axes=(0, 1)))
    assert_eq(
        wfunc(darr2c, (darr2c.shape[0] - 1, darr2c.shape[1] - 1), (0, 1)),
        func(nparrc, (nparrc.shape[0] - 1, nparrc.shape[1] - 1), (0, 1))
    )


@pytest.mark.parametrize("n", [1, 2, 3, 6, 7])
@pytest.mark.parametrize("d", [1.0, 0.5, 2 * np.pi])
@pytest.mark.parametrize("c", [lambda m: m, lambda m: (1, m - 1)])
def test_fftfreq(n, d, c):
    c = c(n)

    r1 = np.fft.fftfreq(n, d)
    r2 = da.fft.fftfreq(n, d, chunks=c)

    assert _normalize_chunks(c, r2.shape) == r2.chunks

    assert_eq(r1, r2)


@pytest.mark.parametrize("n", [1, 2, 3, 6, 7])
@pytest.mark.parametrize("d", [1.0, 0.5, 2 * np.pi])
@pytest.mark.parametrize("c", [lambda m: m // 2 + 1, lambda m: (1, m // 2)])
def test_rfftfreq(n, d, c):
    c = c(n)

    r1 = np.fft.rfftfreq(n, d)
    r2 = da.fft.rfftfreq(n, d, chunks=c)

    assert _normalize_chunks(c, r2.shape) == r2.chunks

    assert_eq(r1, r2)


@pytest.mark.parametrize("funcname", ["fftshift", "ifftshift"])
@pytest.mark.parametrize("axes", [
    None,
    0,
    1,
    2,
    (0, 1),
    (1, 2),
    (0, 2),
    (0, 1, 2),
])
@pytest.mark.parametrize("shape, chunks", [
    [(5, 6, 7), (2, 3, 4)],
    [(5, 6, 7), (2, 6, 4)],
    [(5, 6, 7), (5, 6, 7)],
])
def test_fftshift(funcname, shape, chunks, axes):
    np_func = getattr(np.fft, funcname)
    da_func = getattr(da.fft, funcname)

    a = np.arange(np.prod(shape)).reshape(shape)
    d = da.from_array(a, chunks=chunks)

    a_r = np_func(a, axes)
    d_r = da_func(d, axes)

    for each_d_chunks, each_d_r_chunks in zip(d.chunks, d_r.chunks):
        if len(each_d_chunks) == 1:
            assert len(each_d_r_chunks) == 1
            assert each_d_r_chunks == each_d_chunks
        else:
            assert len(each_d_r_chunks) != 1

    assert_eq(d_r, a_r)


@pytest.mark.parametrize("funcname1, funcname2", [
    ("fftshift", "ifftshift"),
    ("ifftshift", "fftshift"),
])
@pytest.mark.parametrize("axes", [
    None,
    0,
    1,
    2,
    (0, 1),
    (1, 2),
    (0, 2),
    (0, 1, 2),
])
@pytest.mark.parametrize("shape, chunks", [
    [(5, 6, 7), (2, 3, 4)],
    [(5, 6, 7), (2, 6, 4)],
    [(5, 6, 7), (5, 6, 7)],
])
def test_fftshift_identity(funcname1, funcname2, shape, chunks, axes):
    da_func1 = getattr(da.fft, funcname1)
    da_func2 = getattr(da.fft, funcname2)

    a = np.arange(np.prod(shape)).reshape(shape)
    d = da.from_array(a, chunks=chunks)

    d_r = da_func1(da_func2(d, axes), axes)

    for each_d_chunks, each_d_r_chunks in zip(d.chunks, d_r.chunks):
        if len(each_d_chunks) == 1:
            assert len(each_d_r_chunks) == 1
            assert each_d_r_chunks == each_d_chunks
        else:
            assert len(each_d_r_chunks) != 1

    assert_eq(d_r, d)
