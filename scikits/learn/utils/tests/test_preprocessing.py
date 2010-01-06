from numpy.testing import *
import numpy as N

from ..preprocessing import scale, Scaler, nanscale, NanScaler

DEPS = N.finfo(N.float).eps

class test_scale:
    def __init__(self, *args, **kw):
        N.random.seed(0)

    def test_simple(self):
        a = 5 * N.random.randn(10, 4)
        s, t = scale(a)
        assert N.all(a <= 1. + DEPS) and N.all(a >= -1. - DEPS)

    def test_right(self):
        a = 5 * N.random.randn(10, 4)
        s, t = scale(a, 'right')
        assert N.all(a <= 1. + DEPS) and N.all(a >= 0. - DEPS)

    def test_nan_simple(self):
        a = 5 * N.random.randn(10, 4)
        a[9, 0] = N.nan
        # We keep a copy around so that we can check that nanscale does not
        # propagate Nan everywhere.
        at = a.copy()
        s, t = nanscale(at)
        assert N.all(at[N.isfinite(a)] <= 1. + DEPS) and \
               N.all(at[N.isfinite(a)] >= -1. - DEPS)

    def test_nan_right(self):
        a = 5 * N.random.randn(10, 4)
        a[9, 0] = N.nan
        # We keep a copy around so that we can check that nanscale does not
        # propagate Nan everywhere.
        at = a.copy()
        s, t = nanscale(at, 'right')
        assert N.all(at[N.isfinite(a)] <= 1. + DEPS) and \
               N.all(at[N.isfinite(a)] >= 0. - DEPS)

class test_Scaler:
    def __init__(self, *args, **kw):
        N.random.seed(0)

    def _generate_simple(self):
        a = 5 * N.random.randn(10, 4)
        b = 5 * N.random.randn(10, 4)
        return a, b

    def _simple_test(self, a, b, mode):
        at = a.copy()
        bt = b.copy()
        sc = Scaler(a, mode)

        # Test whether preprocess -> unprocess gives back the data
        sc.scale(a)
        if mode == 'sym':
            assert N.all(a <= 1. + DEPS) and N.all(a >= -1. - DEPS)
        elif mode == 'right':
            assert N.all(a <= 1. + DEPS) and N.all(a >= -0. - DEPS)
        else:
            raise ValueError("unexpected mode %s" % str(mode))
        sc.unscale(a)
        assert_array_almost_equal(a, at)

        sc.scale(b)
        sc.unscale(b)
        assert_array_almost_equal(b, bt)

    def _nan_test(self, a, b, mode):
        a[9, 0] = N.nan
        at = a.copy()
        bt = b.copy()
        sc = NanScaler(a, mode)

        # Test whether preprocess -> unprocess gives back the data
        sc.scale(at)
        if mode == 'sym':
            assert N.all(at[N.isfinite(a)] <= 1. + DEPS) and \
                   N.all(at[N.isfinite(a)] >= -1. - DEPS)
        elif mode == 'right':
            assert N.all(at[N.isfinite(a)] <= 1. + DEPS) and \
                   N.all(at[N.isfinite(a)] >= 0. - DEPS)
        else:
            raise ValueError("unexpected mode %s" % str(mode))
        sc.unscale(at)
        assert_array_almost_equal(a[N.isfinite(a)], 
                                  at[N.isfinite(a)])

        sc.scale(b)
        sc.unscale(b)
        assert_array_almost_equal(b, bt)

    def test_simple(self):
        a, b = self._generate_simple()
        self._simple_test(a, b, 'sym')

    def test_right(self):
        a, b = self._generate_simple()
        self._simple_test(a, b, 'right')

    def test_nan(self):
        a, b = self._generate_simple()
        self._nan_test(a, b, 'sym')

    def test_nan_right(self):
        a, b = self._generate_simple()
        self._nan_test(a, b, 'right')

    def test_feature_full_nan(self):
        """Test that instancing a scaler from data whose feature is full of Nan
        does not work."""
        a = N.random.randn(10, 5)
        a[:, 0] = N.nan
        try:
            sc = NanScaler(a)
            assert(0 == 1, "Creation of scaler from data with full Nan "\
                           "should not succeed")
        except ValueError:
            pass

