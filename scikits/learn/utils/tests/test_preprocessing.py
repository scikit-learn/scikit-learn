from numpy.testing import *
import numpy as N

set_package_path()
from preprocessing import scale, Scaler
restore_path()

DEPS = N.finfo(N.float).eps

class test_scale(NumpyTestCase):
    def __init__(self, *args, **kw):
        NumpyTestCase.__init__(self, *args, **kw)
        N.random.seed(0)

    def test_simple(self):
        a = 5 * N.random.randn(10, 4)
        s, t = scale(a)
        assert N.all(a <= 1 + DEPS) and N.all(a >= -1. - DEPS)

    def test_right(self):
        a = 5 * N.random.randn(10, 4)
        s, t = scale(a, 'right')
        assert N.all(a <= 1 + DEPS) and N.all(a >= 0. - DEPS)

class test_Scaler(NumpyTestCase):
    def __init__(self, *args, **kw):
        NumpyTestCase.__init__(self, *args, **kw)
        N.random.seed(0)

    def _common_simple(self, mode):
        a = 5 * N.random.randn(10, 4)
        b = 5 * N.random.randn(10, 4)
        at = a.copy()
        bt = b.copy()
        sc = Scaler(a, mode)

        # Test whether preprocess -> unprocess gives back the data
        sc.preprocess(a)
        if mode == 'sym':
            assert N.all(a <= 1 + DEPS) and N.all(a >= -1. - DEPS)
        elif mode == 'right':
            assert N.all(a <= 1 + DEPS) and N.all(a >= -0. - DEPS)
        else:
            raise ValueError("unepexted mode %s" % str(mode))
        sc.unprocess(a)
        assert_array_almost_equal(a, at)

        sc.preprocess(b)
        sc.unprocess(b)
        assert_array_almost_equal(b, bt)

    def test_simple(self):
        self._common_simple('sym')

    def test_right(self):
        self._common_simple('right')

if __name__ == "__main__":
    NumpyTest('scikits.learn.utils.preprocessing').run()
