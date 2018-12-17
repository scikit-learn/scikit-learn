import pytest
import cython

import numpy as np

from sklearn.utils.testing import assert_allclose

from sklearn.utils._cython_blas import _xdot_memview
from sklearn.utils._cython_blas import _xasum_memview
from sklearn.utils._cython_blas import _xaxpy_memview
from sklearn.utils._cython_blas import _xnrm2_memview
from sklearn.utils._cython_blas import _xcopy_memview
from sklearn.utils._cython_blas import _xscal_memview
from sklearn.utils._cython_blas import _xgemv_memview
from sklearn.utils._cython_blas import _xger_memview
from sklearn.utils._cython_blas import _xgemm_memview


NUMPY_TO_CYTHON = {np.float32: cython.float, np.float64: cython.double}
RTOL = {np.float32: 1e-3, np.float64: 1e-13}


def _no_op(x):
    return x


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_dot(dtype):
    dot = _xdot_memview[NUMPY_TO_CYTHON[dtype]]

    rng = np.random.RandomState(0)
    x = rng.random_sample(10).astype(dtype, copy=False)
    y = rng.random_sample(10).astype(dtype, copy=False)

    expected = x.dot(y)
    actual = dot(x, y)

    assert_allclose(actual, expected)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_asum(dtype):
    asum = _xasum_memview[NUMPY_TO_CYTHON[dtype]]

    rng = np.random.RandomState(0)
    x = rng.random_sample(10).astype(dtype, copy=False)

    expected = np.abs(x).sum()
    actual = asum(x)

    assert_allclose(actual, expected)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_axpy(dtype):
    axpy = _xaxpy_memview[NUMPY_TO_CYTHON[dtype]]

    rng = np.random.RandomState(0)
    x = rng.random_sample(10).astype(dtype, copy=False)
    y = rng.random_sample(10).astype(dtype, copy=False)
    alpha = 1.23

    expected = alpha * x + y
    axpy(alpha, x, y)

    assert_allclose(y, expected)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_nrm2(dtype):
    nrm2 = _xnrm2_memview[NUMPY_TO_CYTHON[dtype]]

    rng = np.random.RandomState(0)
    x = rng.random_sample(10).astype(dtype, copy=False)

    expected = np.linalg.norm(x)
    actual = nrm2(x)

    assert_allclose(actual, expected)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_copy(dtype):
    copy = _xcopy_memview[NUMPY_TO_CYTHON[dtype]]

    rng = np.random.RandomState(0)
    x = rng.random_sample(10).astype(dtype, copy=False)
    y = np.empty_like(x)

    expected = x.copy()
    copy(x, y)

    assert_allclose(y, expected)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_scal(dtype):
    scal = _xscal_memview[NUMPY_TO_CYTHON[dtype]]

    rng = np.random.RandomState(0)
    x = rng.random_sample(10).astype(dtype, copy=False)
    alpha = 1.23

    expected = alpha * x
    scal(alpha, x)

    assert_allclose(x, expected)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("opA, transA",
                         [(_no_op, 'n'), (np.transpose, 't')],
                         ids=["A", "A.T"])
@pytest.mark.parametrize("layout", ['C', 'F'])
def test_gemv(dtype, opA, transA, layout):
    gemv = _xgemv_memview[NUMPY_TO_CYTHON[dtype]]

    rng = np.random.RandomState(0)
    A = np.asarray(opA(rng.random_sample((20, 10)).astype(dtype, copy=False)),
                   order=layout)
    x = rng.random_sample(10).astype(dtype, copy=False)
    y = rng.random_sample(20).astype(dtype, copy=False)
    alpha, beta = 1.23, -3.21

    expected = alpha * opA(A).dot(x) + beta * y
    gemv(layout, transA, alpha, A, x, beta, y)

    assert_allclose(y, expected, rtol=RTOL[dtype])


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("layout", ['C', 'F'])
def test_ger(dtype, layout):
    ger = _xger_memview[NUMPY_TO_CYTHON[dtype]]

    rng = np.random.RandomState(0)
    x = rng.random_sample(10).astype(dtype, copy=False)
    y = rng.random_sample(20).astype(dtype, copy=False)
    A = np.asarray(rng.random_sample((10, 20)).astype(dtype, copy=False),
                   order=layout)
    alpha = 1.23

    expected = alpha * np.outer(x, y) + A
    ger(layout, alpha, x, y, A)

    assert_allclose(A, expected, rtol=RTOL[dtype])


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("opB, transB",
                         [(_no_op, 'n'), (np.transpose, 't')],
                         ids=["B", "B.T"])
@pytest.mark.parametrize("opA, transA",
                         [(_no_op, 'n'), (np.transpose, 't')],
                         ids=["A", "A.T"])
@pytest.mark.parametrize("layout", ['C', 'F'])
def test_gemm(dtype, opA, transA, opB, transB, layout):
    gemm = _xgemm_memview[NUMPY_TO_CYTHON[dtype]]

    rng = np.random.RandomState(0)
    A = np.asarray(opA(rng.random_sample((30, 10)).astype(dtype, copy=False)),
                   order=layout)
    B = np.asarray(opB(rng.random_sample((10, 20)).astype(dtype, copy=False)),
                   order=layout)
    C = np.asarray(rng.random_sample((30, 20)).astype(dtype, copy=False),
                   order=layout)
    alpha, beta = 1.23, -3.21

    expected = alpha * opA(A).dot(opB(B)) + beta * C
    gemm(layout, transA, transB, alpha, A, B, beta, C)

    assert_allclose(C, expected, rtol=RTOL[dtype])
