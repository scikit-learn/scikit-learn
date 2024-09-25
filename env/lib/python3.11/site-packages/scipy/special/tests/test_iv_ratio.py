# This file contains unit tests for the iv_ratio() function.

import pytest
import numpy as np
from numpy.testing import assert_equal, assert_allclose
from scipy.special._ufuncs import _iv_ratio as iv_ratio  # type: ignore[attr-defined]


class TestIvRatio:

    @pytest.mark.parametrize('v,x,r', [
        (1, 0.3380952380952381, 0.1666773049170313),
        (1, 0.7083333333333333, 0.33366443586989925),
        (1, 1.1666666666666667, 0.5023355231537423),
        (1, 1.8666666666666665, 0.674616572252164),
        (1, 3.560606060606061, 0.844207659503163),
        (2.34, 0.7975238095238094, 0.16704903081553285),
        (2.34, 1.7133333333333334, 0.3360215931268845),
        (2.34, 2.953333333333333, 0.50681909317803),
        (2.34, 5.0826666666666656, 0.6755252698800679),
        (2.34, 10.869696969696973, 0.8379351104498762),
        (56.789, 19.46575238095238, 0.1667020505391409),
        (56.789, 42.55008333333333, 0.33353809996933026),
        (56.789, 75.552, 0.5003932381177826),
        (56.789, 135.76026666666667, 0.6670528221946127),
        (56.789, 307.8642424242425, 0.8334999441460798),
    ])
    def test_against_reference_values(self, v, x, r):
        """The reference values are computed using mpmath as follows.

        from mpmath import mp
        mp.dps = 100

        def iv_ratio_mp(v, x):
            return mp.besseli(v, x) / mp.besseli(v - 1, x)

        def _sample(n, *, v):
            '''Return n positive real numbers x such that iv_ratio(v, x) are
            roughly evenly spaced over (0, 1).  The formula is taken from [1].

            [1] Banerjee A., Dhillon, I. S., Ghosh, J., Sra, S. (2005).
                "Clustering on the Unit Hypersphere using von Mises-Fisher
                Distributions."  Journal of Machine Learning Research,
                6(46):1345-1382.
            '''
            r = np.arange(1, n+1) / (n+1)
            return r * (2*v-r*r) / (1-r*r)

        for v in (1, 2.34, 56.789):
            xs = _sample(5, v=v)
            for x in xs:
                print(f"({v}, {x}, {float(iv_ratio_mp_float(v,x))}),")
        """
        assert_allclose(iv_ratio(v, x), r, rtol=4e-16, atol=0)

    @pytest.mark.parametrize('v,x,r', [
        (1, np.inf, 1),
        (np.inf, 1, 0),
    ])
    def test_inf(self, v, x, r):
        """If exactly one of v or x is inf and the other is within domain,
        should return 0 or 1 accordingly.

        Also check that the function
        never returns -0.0."""
        assert_equal(iv_ratio(v, x), r)

    @pytest.mark.parametrize('v', [np.nextafter(1, 0), -np.inf, np.nan, np.inf])
    @pytest.mark.parametrize('x', [-np.finfo(float).smallest_normal,
                                   -np.finfo(float).smallest_subnormal,
                                   -np.inf, np.nan, np.inf])
    def test_nan(self, v, x):
        """If at least one argument is out of domain, or if v = x = inf,
        the function should return nan."""
        assert_equal(iv_ratio(v, x), np.nan)

    @pytest.mark.parametrize('v', [1, np.finfo(float).max, np.inf])
    def test_zero_x(self, v):
        """If x is +/-0.0, return x to agree with the limiting behavior."""
        assert_equal(iv_ratio(v, 0.0), 0.0)
        assert_equal(iv_ratio(v, -0.0), -0.0)

    @pytest.mark.parametrize('v,x', [
        (1, np.finfo(float).smallest_normal),
        (1, np.finfo(float).smallest_subnormal),
        (1, np.finfo(float).smallest_subnormal*2),
        (1e20, 123),
        (np.finfo(float).max, 1),
        (np.finfo(float).max, np.sqrt(np.finfo(float).max)),
    ])
    def test_tiny_x(self, v, x):
        """If x is much less than v, the bounds

                    x                                 x
        --------------------------- <= R <= -----------------------
        v-0.5+sqrt(x**2+(v+0.5)**2)         v-1+sqrt(x**2+(v+1)**2)

        collapses to R ~= x/2v.  Test against this asymptotic expression.
        """
        assert_equal(iv_ratio(v, x), (0.5*x)/v)

    @pytest.mark.parametrize('v,x', [
        (1, 1e16),
        (1e20, 1e40),
        (np.sqrt(np.finfo(float).max), np.finfo(float).max),
    ])
    def test_huge_x(self, v, x):
        """If x is much greater than v, the bounds

                    x                                 x
        --------------------------- <= R <= -----------------------
        v-0.5+sqrt(x**2+(v+0.5)**2)         v-1+sqrt(x**2+(v+1)**2)

        collapses to R ~= 1.  Test against this asymptotic expression.
        """
        assert_equal(iv_ratio(v, x), 1.0)

    @pytest.mark.parametrize('v,x', [
        (np.finfo(float).max, np.finfo(float).max),
        (np.finfo(float).max / 3, np.finfo(float).max),
        (np.finfo(float).max, np.finfo(float).max / 3),
    ])
    def test_huge_v_x(self, v, x):
        """If both x and v are very large, the bounds

                    x                                 x
        --------------------------- <= R <= -----------------------
        v-0.5+sqrt(x**2+(v+0.5)**2)         v-1+sqrt(x**2+(v+1)**2)

        collapses to R ~= x/(v+sqrt(x**2+v**2).  Test against this asymptotic
        expression, and in particular that no numerical overflow occurs during
        intermediate calculations.
        """
        t = x / v
        expected = t / (1 + np.hypot(1, t))
        assert_allclose(iv_ratio(v, x), expected, rtol=4e-16, atol=0)
