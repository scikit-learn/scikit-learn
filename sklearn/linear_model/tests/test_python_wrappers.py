import numpy as np
import scipy
from sklearn.utils.testing import (assert_equal, assert_array_equal,
                                   assert_almost_equal)
from sklearn.utils import check_random_state
from sklearn.linear_model.python_wrappers import (
    _py_fused_dotc, _py_fused_dotu, _py_fused_asum, _py_fused_nrm2,
    _py_fused_nrm2_squared, _py_abs_max, _py_fused_copy, _py_diff_abs_max,
    _py_real_max)
from .prox_slow import nrm2
from .test_coordescendant import as_complex

IMAG_UNIT = scipy.sqrt(-1)


def test_fused_dot(random_state=0):
    rng = check_random_state(random_state)
    x = np.array([3. - 4 * IMAG_UNIT])
    assert_equal(_py_fused_dotc(x, x), 25.)

    for dtype in [np.float32, np.float64, np.complex64, np.complex128]:
        x = dtype(as_complex(*rng.randn(2, 10)))
        y = dtype(as_complex(*rng.randn(2, 10)))
    assert_equal(_py_fused_dotc(x, y), _py_fused_dotu(x.conj(), y))


def test_fused_asum(random_state=1):
    rng = check_random_state(random_state)
    for real in [True, False]:
        x = rng.randn(10)
        if not real:
            x = as_complex(x, rng.randn(*x.shape))
        for dtype in [np.float32, np.float64, np.complex64, np.complex128]:
            x = dtype(x)
            if dtype in [np.float32, np.complex64]:
                decimal = 6
            else:
                decimal = 13
            assert_almost_equal(_py_fused_asum(x.copy()), np.abs(x).sum(),
                                decimal=decimal)


def test_fused_nrm2(random_state=0):
    rng = check_random_state(random_state)
    for real in [True, False]:
        x = rng.randn(10)
        if not real:
            x = as_complex(x, rng.randn(*x.shape))
        for dtype in [np.float32, np.float64, np.complex64, np.complex128]:
            x = dtype(x)
            if dtype in [np.float32, np.complex64]:
                decimal = 6
            else:
                decimal = 13
            assert_almost_equal(_py_fused_nrm2(x), nrm2(x), decimal=decimal)


def test_fused_nrm2_squared(random_state=0):
    rng = check_random_state(random_state)
    for real in [True, False]:
        x = rng.randn(100)
        if not real:
            x = as_complex(x, rng.randn(*x.shape))
        for dtype in [np.float32, np.float64, np.complex64, np.complex128]:
            x = dtype(x)
            if dtype in [np.float32, np.complex64]:
                decimal = 5
            else:
                decimal = 13
            assert_almost_equal(_py_fused_nrm2_squared(x),
                                nrm2(x, squared=True), decimal=decimal)


def test_abs_max(random_state=0):
    rng = check_random_state(random_state)
    x = rng.randn(10)
    for real in [True, False]:
        x = rng.randn(10)
        if not real:
            x = as_complex(x, rng.randn(*x.shape))
        for dtype in [np.float32, np.float64, np.complex64, np.complex128]:
            x = dtype(x)
            assert_almost_equal(_py_abs_max(x), np.abs(x).max(),
                                decimal=13)


def test_real_max(random_state=0):
    rng = check_random_state(random_state)
    x = rng.randn(10)
    for dtype in [np.float32, np.float64]:
        x = dtype(x)
        assert_almost_equal(_py_real_max(x), np.max(x), decimal=13)


def test_fused_copy(random_state=0):
    rng = check_random_state(random_state)
    for real in [True, False]:
        x = rng.randn(10)
        if not real:
            x = as_complex(x, rng.randn(*x.shape))
        for dtype in [np.float32, np.float64, np.complex64, np.complex128]:
            x = dtype(x)
            y = np.zeros_like(x)
            _py_fused_copy(x, y)
            assert_array_equal(x, y)


def test_diff_abs_max(random_state=0):
    x = np.array([0.2524121, 0.15189841, -0.39494169, 0.57979797,
                  -0.00546893, 0.31447354, 0.2821645, 0.0450028,
                  0.01339086, 0.14667659, 0.16338791, -0.30872306,
                  0.00219627, 0.38066833 -0.0285587, -0.05008872,
                  -0.23163258, 0.59517036, 0.13490514, -0.03699011])
    y = np.zeros_like(x)
    assert_almost_equal(_py_diff_abs_max(x, y), np.abs(x - y).max())

    rng = check_random_state(random_state)
    for real in [True, False]:
        for dtype in [np.float32, np.float64]:
            x, y = map(dtype, rng.randn(2, 10))
            if not real:
                x = as_complex(x, dtype(rng.randn(*x.shape)))
                y = as_complex(y, dtype(rng.randn(*x.shape)))
            assert_almost_equal(_py_diff_abs_max(x, y), np.abs(x - y).max())
