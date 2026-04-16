import numpy as np
from numpy.testing import assert_allclose, assert_equal
import pytest
from scipy.special._ufuncs import _gen_harmonic, _normalized_gen_harmonic


#
# In the following tests, reference values were computed with mpmath.
#

@pytest.mark.parametrize('typ', [np.int32, np.int64, np.float64])
@pytest.mark.parametrize(
    'n, a, ref',
    [(8, 9.0, 1.0020083884212339),
     (1000, 2.5, 1.3414661912046497),
     (10, 1.5, 1.9953364933456017),
     (10000, 1.25, 4.1951168257387765),
     (10000,1.00001, 9.787182620770265),
     (80, 1.000002, 4.965460167788836),
     (75, 1 + 1e-12, 4.901355630543771),
     (100, 1 + 1e-14, 5.187377517639515),
     (100, 1 + 8e-16, 5.187377517639611),
     (100, 1.0, 5.187377517639621),
     (7, 1.0, 2.592857142857143),
     (8000, 1.0, 9.564474984261423),
     (5, 1 - 1e-12, 2.2833333333347143),
     (25000, 1 - 1e-12, 10.703866768669737),
     (1000, 0.995, 7.6058022857089975),
     (1000, 0.75, 19.055178975831392),
     (10000, 0.25, 1332.5700547197382),
     (5, 1e-8, 4.999999952125083),
     (15, 1e-16, 14.999999999999996),
     (100, 0.0, 100.0),
     (4, -1.0, 10.0),
     (75, -1.5, 19811.38815892374)]
)
def test_gen_harmonic(typ, n, a, ref):
    h = _gen_harmonic(typ(n), a)
    assert_allclose(h, ref, rtol=5e-15)


@pytest.mark.parametrize('typ', [np.int32, np.int64, np.float64])
@pytest.mark.parametrize(
    'n, a, ref',
    [(10, np.inf, 1.0),
     (1, np.nan, 1.0),
     (1, -np.inf, 1.0),
     (3, np.nan, np.nan),
     (-3, 1.0, np.nan)]
)
def test_gen_harmonic_exact_cases(typ, n, a, ref):
    h = _gen_harmonic(typ(n), a)
    assert_equal(h, ref)


def test_gen_harmonic_n_nan():
    h = _gen_harmonic(np.nan, 0.75)
    assert_equal(h, np.nan)


@pytest.mark.parametrize('typ', [np.int32, np.int64, np.float64])
@pytest.mark.parametrize(
    'j, k, n, a, ref',
    [(400, 5000, 5000, 10.0, 4.2821759663214485e-25),
     (400, 5000, 5000, 3.5, 1.11086549102426e-07),
     (1, 2, 3, 1.5, 0.8755176866163012),
     (300, 500, 500, 1 + 1e-14, 0.07559343891632035),
     (1500, 2500, 3000, 1 - 1e-12, 0.05957291246371843),
     (10, 12, 16, 0.5, 0.13601665344521513),
     (16, 16, 20, 0.125, 0.04583107002260924),
     (10, 12, 16, -0.5, 0.22359306724308234),
     (1, 8000, 10000, -1.5, 0.5724512895513029)]
)
def test_normalized_gen_harmonic(typ, j, k, n, a, ref):
    h = _normalized_gen_harmonic(typ(j), typ(k), typ(n), a)
    assert_allclose(h, ref, 5e-15)


@pytest.mark.parametrize('typ', [np.int32, np.int64, np.float64])
@pytest.mark.parametrize(
    'j, k, n, a, ref',
    [(1, 1, 1, 0.5, 1.0),
     (1, 1, 1, np.nan, 1.0),
     (1, 2, 5, np.nan, np.nan),
     (1, 2, 1, 1.25, np.nan),
     (1, 2, 3, np.inf, 1.0),
     (2, 3, 4, np.inf, 0.0),
     (1, 1, 10, -np.inf, 0.0),
     (2, 3, 4, -np.inf, np.nan),
     (3, 6, 8, 0.0, 0.5)]
)
def test_normalized_gen_harmonic_exact_cases(typ, j, k, n, a, ref):
    h = _normalized_gen_harmonic(typ(j), typ(k), typ(n), a)
    assert_equal(h, ref)


def test_normalized_gen_harmonic_input_nan():
    h = _normalized_gen_harmonic(1.0, np.nan, 10.0, 1.05)
    assert_equal(h, np.nan)
