import numpy as np
from numpy.testing import assert_, assert_equal
from skimage._shared import testing
import unittest

from skimage import data
from skimage.filters import LPIFilter2D, inverse, wiener


class TestLPIFilter2D(unittest.TestCase):
    img = data.camera()[:50, :50]

    def filt_func(self, r, c):
        return np.exp(-np.hypot(r, c) / 1)

    @testing.fixture(autouse=True)
    def setUp(self):
        self.f = LPIFilter2D(self.filt_func)

    def tst_shape(self, x):
        X = self.f(x)
        assert_equal(X.shape, x.shape)

    def test_ip_shape(self):
        rows, columns = self.img.shape[:2]

        for c_slice in [slice(0, columns), slice(0, columns - 5),
                        slice(0, columns - 20)]:
            yield (self.tst_shape, self.img[:, c_slice])

    def test_inverse(self):
        F = self.f(self.img)
        g = inverse(F, predefined_filter=self.f)
        assert_equal(g.shape, self.img.shape)

        g1 = inverse(F[::-1, ::-1], predefined_filter=self.f)
        assert_((g - g1[::-1, ::-1]).sum() < 55)

        # test cache
        g1 = inverse(F[::-1, ::-1], predefined_filter=self.f)
        assert_((g - g1[::-1, ::-1]).sum() < 55)

        g1 = inverse(F[::-1, ::-1], self.filt_func)
        assert_((g - g1[::-1, ::-1]).sum() < 55)

    def test_wiener(self):
        F = self.f(self.img)
        g = wiener(F, predefined_filter=self.f)
        assert_equal(g.shape, self.img.shape)

        g1 = wiener(F[::-1, ::-1], predefined_filter=self.f)
        assert_((g - g1[::-1, ::-1]).sum() < 1)

        g1 = wiener(F[::-1, ::-1], self.filt_func)
        assert_((g - g1[::-1, ::-1]).sum() < 1)

    def test_non_callable(self):
        with testing.raises(ValueError):
            LPIFilter2D(None)
