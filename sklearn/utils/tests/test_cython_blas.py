import pytest
import cython

import numpy as np

from sklearn.utils.testing import assert_allclose

from sklearn.utils._cython_blas import _dot_memview
from sklearn.utils._cython_blas import _asum_memview
from sklearn.utils._cython_blas import _axpy_memview
from sklearn.utils._cython_blas import _nrm2_memview
from sklearn.utils._cython_blas import _copy_memview
from sklearn.utils._cython_blas import _scal_memview
from sklearn.utils._cython_blas import _gemv_memview
from sklearn.utils._cython_blas import _ger_memview
from sklearn.utils._cython_blas import _gemm_memview
from sklearn.utils._cython_blas import RowMajor, ColMajor
from sklearn.utils._cython_blas import Trans, NoTrans


NUMPY_TO_CYTHON = {np.float32: cython.float, np.float64: cython.double}
RTOL = {np.float32: 1e-6, np.float64: 1e-12}
ORDER = {RowMajor: 'C', ColMajor: 'F'}


def _no_op(x):
    return x


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_dot(dtype):
    dot = _dot_memview[NUMPY_TO_CYTHON[dtype]]

    rng = np.random.RandomState(0)
    x = rng.random_sample(10).astype(dtype, copy=False)
    y = rng.random_sample(10).astype(dtype, copy=False)

    expected = x.dot(y)
    actual = dot(x, y)

    assert_allclose(actual, expected, rtol=RTOL[dtype])


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_asum(dtype):
    asum = _asum_memview[NUMPY_TO_CYTHON[dtype]]

    rng = np.random.RandomState(0)
    x = rng.random_sample(10).astype(dtype, copy=False)

    expected = np.abs(x).sum()
    actual = asum(x)

    assert_allclose(actual, expected, rtol=RTOL[dtype])


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_axpy(dtype):
    axpy = _axpy_memview[NUMPY_TO_CYTHON[dtype]]

    rng = np.random.RandomState(0)
    x = rng.random_sample(10).astype(dtype, copy=False)
    y = rng.random_sample(10).astype(dtype, copy=False)
    alpha = 2.5

    expected = alpha * x + y
    axpy(alpha, x, y)

    assert_allclose(y, expected, rtol=RTOL[dtype])


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_nrm2(dtype):
    nrm2 = _nrm2_memview[NUMPY_TO_CYTHON[dtype]]

    rng = np.random.RandomState(0)
    x = rng.random_sample(10).astype(dtype, copy=False)

    expected = np.linalg.norm(x)
    actual = nrm2(x)

    assert_allclose(actual, expected, rtol=RTOL[dtype])


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_copy(dtype):
    copy = _copy_memview[NUMPY_TO_CYTHON[dtype]]

    rng = np.random.RandomState(0)
    x = rng.random_sample(10).astype(dtype, copy=False)
    y = np.empty_like(x)

    expected = x.copy()
    copy(x, y)

    assert_allclose(y, expected, rtol=RTOL[dtype])


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_scal(dtype):
    scal = _scal_memview[NUMPY_TO_CYTHON[dtype]]

    rng = np.random.RandomState(0)
    x = rng.random_sample(10).astype(dtype, copy=False)
    alpha = 2.5

    expected = alpha * x
    scal(alpha, x)

    assert_allclose(x, expected, rtol=RTOL[dtype])


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("opA, transA",
                         [(_no_op, NoTrans), (np.transpose, Trans)],
                         ids=["NoTrans", "Trans"])
@pytest.mark.parametrize("order", [RowMajor, ColMajor],
                         ids=["RowMajor", "ColMajor"])
def test_gemv(dtype, opA, transA, order):
    gemv = _gemv_memview[NUMPY_TO_CYTHON[dtype]]

    rng = np.random.RandomState(0)
    A = np.asarray(opA(rng.random_sample((20, 10)).astype(dtype, copy=False)),
                   order=ORDER[order])
    x = rng.random_sample(10).astype(dtype, copy=False)
    y = rng.random_sample(20).astype(dtype, copy=False)
    alpha, beta = 2.5, -0.5

    expected = alpha * opA(A).dot(x) + beta * y
    gemv(transA, alpha, A, x, beta, y)

    assert_allclose(y, expected, rtol=RTOL[dtype])


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("order", [RowMajor, ColMajor],
                         ids=["RowMajor", "ColMajor"])
def test_ger(dtype, order):
    ger = _ger_memview[NUMPY_TO_CYTHON[dtype]]

    rng = np.random.RandomState(0)
    x = rng.random_sample(10).astype(dtype, copy=False)
    y = rng.random_sample(20).astype(dtype, copy=False)
    A = np.asarray(rng.random_sample((10, 20)).astype(dtype, copy=False),
                   order=ORDER[order])
    alpha = 2.5

    expected = alpha * np.outer(x, y) + A
    ger(alpha, x, y, A)

    assert_allclose(A, expected, rtol=RTOL[dtype])


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("opB, transB",
                         [(_no_op, NoTrans), (np.transpose, Trans)],
                         ids=["NoTrans", "Trans"])
@pytest.mark.parametrize("opA, transA",
                         [(_no_op, NoTrans), (np.transpose, Trans)],
                         ids=["NoTrans", "Trans"])
@pytest.mark.parametrize("order", [RowMajor, ColMajor],
                         ids=["RowMajor", "ColMajor"])
def test_gemm(dtype, opA, transA, opB, transB, order):
    gemm = _gemm_memview[NUMPY_TO_CYTHON[dtype]]

    rng = np.random.RandomState(0)
    A = np.asarray(opA(rng.random_sample((30, 10)).astype(dtype, copy=False)),
                   order=ORDER[order])
    B = np.asarray(opB(rng.random_sample((10, 20)).astype(dtype, copy=False)),
                   order=ORDER[order])
    C = np.asarray(rng.random_sample((30, 20)).astype(dtype, copy=False),
                   order=ORDER[order])
    alpha, beta = 2.5, -0.5

    expected = alpha * opA(A).dot(opB(B)) + beta * C
    gemm(transA, transB, alpha, A, B, beta, C)

    assert_allclose(C, expected, rtol=RTOL[dtype])
