import unittest
import numpy as np
from skimage._shared.testing import assert_equal
from scipy.ndimage import binary_dilation, binary_erosion
import skimage.feature as F
from skimage import data, img_as_float


class TestCanny(unittest.TestCase):
    def test_00_00_zeros(self):
        '''Test that the Canny filter finds no points for a blank field'''
        result = F.canny(np.zeros((20, 20)), 4, 0, 0, np.ones((20, 20), bool))
        self.assertFalse(np.any(result))

    def test_00_01_zeros_mask(self):
        '''Test that the Canny filter finds no points in a masked image'''
        result = (F.canny(np.random.uniform(size=(20, 20)), 4, 0, 0,
                          np.zeros((20, 20), bool)))
        self.assertFalse(np.any(result))

    def test_01_01_circle(self):
        '''Test that the Canny filter finds the outlines of a circle'''
        i, j = np.mgrid[-200:200, -200:200].astype(float) / 200
        c = np.abs(np.sqrt(i * i + j * j) - .5) < .02
        result = F.canny(c.astype(float), 4, 0, 0, np.ones(c.shape, bool))
        #
        # erode and dilate the circle to get rings that should contain the
        # outlines
        #
        cd = binary_dilation(c, iterations=3)
        ce = binary_erosion(c, iterations=3)
        cde = np.logical_and(cd, np.logical_not(ce))
        self.assertTrue(np.all(cde[result]))
        #
        # The circle has a radius of 100. There are two rings here, one
        # for the inside edge and one for the outside. So that's
        # 100 * 2 * 2 * 3 for those places where pi is still 3.
        # The edge contains both pixels if there's a tie, so we
        # bump the count a little.
        point_count = np.sum(result)
        self.assertTrue(point_count > 1200)
        self.assertTrue(point_count < 1600)

    def test_01_02_circle_with_noise(self):
        '''Test that the Canny filter finds the circle outlines
         in a noisy image'''
        np.random.seed(0)
        i, j = np.mgrid[-200:200, -200:200].astype(float) / 200
        c = np.abs(np.sqrt(i * i + j * j) - .5) < .02
        cf = c.astype(float) * .5 + np.random.uniform(size=c.shape) * .5
        result = F.canny(cf, 4, .1, .2, np.ones(c.shape, bool))
        #
        # erode and dilate the circle to get rings that should contain the
        # outlines
        #
        cd = binary_dilation(c, iterations=4)
        ce = binary_erosion(c, iterations=4)
        cde = np.logical_and(cd, np.logical_not(ce))
        self.assertTrue(np.all(cde[result]))
        point_count = np.sum(result)
        self.assertTrue(point_count > 1200)
        self.assertTrue(point_count < 1600)

    def test_image_shape(self):
        self.assertRaises(ValueError, F.canny, np.zeros((20, 20, 20)), 4, 0, 0)

    def test_mask_none(self):
        result1 = F.canny(np.zeros((20, 20)), 4, 0, 0, np.ones((20, 20), bool))
        result2 = F.canny(np.zeros((20, 20)), 4, 0, 0)
        self.assertTrue(np.all(result1 == result2))

    def test_use_quantiles(self):
        image = img_as_float(data.camera()[::50, ::50])

        # Correct output produced manually with quantiles
        # of 0.8 and 0.6 for high and low respectively
        correct_output = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0],
           [0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0],
           [0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0],
           [0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0],
           [0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0],
           [0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0],
           [0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0],
           [0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]], dtype=bool)

        result = F.canny(image, low_threshold=0.6, high_threshold=0.8, use_quantiles=True)

        assert_equal(result, correct_output)

    def test_invalid_use_quantiles(self):
        image = img_as_float(data.camera()[::50, ::50])

        self.assertRaises(ValueError, F.canny, image, use_quantiles=True,
                          low_threshold=0.5, high_threshold=3.6)

        self.assertRaises(ValueError, F.canny, image, use_quantiles=True,
                          low_threshold=-5, high_threshold=0.5)

        self.assertRaises(ValueError, F.canny, image, use_quantiles=True,
                          low_threshold=99, high_threshold=0.9)

        self.assertRaises(ValueError, F.canny, image, use_quantiles=True,
                          low_threshold=0.5, high_threshold=-100)
