import math
import unittest

import numpy as np

from skimage.morphology import extrema
from scipy import ndimage as ndi

eps = 1e-12


def diff(a, b):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    t = ((a - b)**2).sum()
    return math.sqrt(t)


class TestExtrema(unittest.TestCase):

    def test_saturated_arithmetic(self):
        "Adding/subtracting a constant and clipping"
        # Test for unsigned integer
        data = np.array([[250, 251, 5, 5],
                         [100, 200, 253, 252],
                         [4, 10, 1, 3]],
                        dtype=np.uint8)
        # adding the constant
        img_constant_added = extrema._add_constant_clip(data, 4)
        expected = np.array([[254, 255, 9, 9],
                             [104, 204, 255, 255],
                             [8, 14, 5, 7]],
                            dtype=np.uint8)
        error = diff(img_constant_added, expected)
        assert error < eps
        img_constant_subtracted = extrema._subtract_constant_clip(data, 4)
        expected = np.array([[246, 247, 1, 1],
                             [96, 196, 249, 248],
                             [0, 6, 0, 0]],
                            dtype=np.uint8)
        error = diff(img_constant_subtracted, expected)
        assert error < eps

        # Test for signed integer
        data = np.array([[32767, 32766],
                         [-32768, -32767]],
                        dtype=np.int16)
        img_constant_added = extrema._add_constant_clip(data, 1)
        expected = np.array([[32767, 32767],
                             [-32767, -32766]],
                            dtype=np.int16)
        error = diff(img_constant_added, expected)
        assert error < eps
        img_constant_subtracted = extrema._subtract_constant_clip(data, 1)
        expected = np.array([[32766, 32765],
                             [-32768, -32768]],
                            dtype=np.int16)
        error = diff(img_constant_subtracted, expected)
        assert error < eps

    def test_local_maxima(self):
        "local maxima for various data types"
        data = np.array([[10,  11,  13,  14,  14,  15,  14,  14,  13,  11],
                         [11,  13,  15,  16,  16,  16,  16,  16,  15,  13],
                         [13,  15,  40,  40,  18,  18,  18,  60,  60,  15],
                         [14,  16,  40,  40,  19,  19,  19,  60,  60,  16],
                         [14,  16,  18,  19,  19,  19,  19,  19,  18,  16],
                         [15,  16,  18,  19,  19,  20,  19,  19,  18,  16],
                         [14,  16,  18,  19,  19,  19,  19,  19,  18,  16],
                         [14,  16,  80,  80,  19,  19,  19, 100, 100,  16],
                         [13,  15,  80,  80,  18,  18,  18, 100, 100,  15],
                         [11,  13,  15,  16,  16,  16,  16,  16,  15,  13]],
                        dtype=np.uint8)
        expected_result = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 1, 1, 0, 0, 0, 1, 1, 0],
                                    [0, 0, 1, 1, 0, 0, 0, 1, 1, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 1, 1, 0, 0, 0, 1, 1, 0],
                                    [0, 0, 1, 1, 0, 0, 0, 1, 1, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
                                   dtype=np.uint8)
        for dtype in [np.uint8, np.uint64, np.int8, np.int64]:

            test_data = data.astype(dtype)
            out = extrema.local_maxima(test_data)

            error = diff(expected_result, out)
            assert error < eps
            assert out.dtype == expected_result.dtype

    def test_local_minima(self):
        "local minima for various data types"

        data = np.array([[10,  11,  13,  14,  14,  15,  14,  14,  13,  11],
                         [11,  13,  15,  16,  16,  16,  16,  16,  15,  13],
                         [13,  15,  40,  40,  18,  18,  18,  60,  60,  15],
                         [14,  16,  40,  40,  19,  19,  19,  60,  60,  16],
                         [14,  16,  18,  19,  19,  19,  19,  19,  18,  16],
                         [15,  16,  18,  19,  19,  20,  19,  19,  18,  16],
                         [14,  16,  18,  19,  19,  19,  19,  19,  18,  16],
                         [14,  16,  80,  80,  19,  19,  19, 100, 100,  16],
                         [13,  15,  80,  80,  18,  18,  18, 100, 100,  15],
                         [11,  13,  15,  16,  16,  16,  16,  16,  15,  13]],
                        dtype=np.uint8)
        data = 100 - data
        expected_result = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 1, 1, 0, 0, 0, 1, 1, 0],
                                    [0, 0, 1, 1, 0, 0, 0, 1, 1, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 1, 1, 0, 0, 0, 1, 1, 0],
                                    [0, 0, 1, 1, 0, 0, 0, 1, 1, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
                                   dtype=np.uint8)
        for dtype in [np.uint8, np.uint64, np.int8, np.int64]:
            data = data.astype(dtype)
            out = extrema.local_minima(data)

            error = diff(expected_result, out)
            assert error < eps
            assert out.dtype == expected_result.dtype

    def test_h_maxima(self):
        "h-maxima for various data types"

        data = np.array([[10,  11,  13,  14,  14,  15,  14,  14,  13,  11],
                         [11,  13,  15,  16,  16,  16,  16,  16,  15,  13],
                         [13,  15,  40,  40,  18,  18,  18,  60,  60,  15],
                         [14,  16,  40,  40,  19,  19,  19,  60,  60,  16],
                         [14,  16,  18,  19,  19,  19,  19,  19,  18,  16],
                         [15,  16,  18,  19,  19,  20,  19,  19,  18,  16],
                         [14,  16,  18,  19,  19,  19,  19,  19,  18,  16],
                         [14,  16,  80,  80,  19,  19,  19, 100, 100,  16],
                         [13,  15,  80,  80,  18,  18,  18, 100, 100,  15],
                         [11,  13,  15,  16,  16,  16,  16,  16,  15,  13]],
                        dtype=np.uint8)

        expected_result = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 1, 1, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 1, 1, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 1, 1, 0, 0, 0, 1, 1, 0],
                                    [0, 0, 1, 1, 0, 0, 0, 1, 1, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
                                   dtype=np.uint8)
        for dtype in [np.uint8, np.uint64, np.int8, np.int64]:
            data = data.astype(dtype)
            out = extrema.h_maxima(data, 40)

            error = diff(expected_result, out)
            assert error < eps

    def test_h_minima(self):
        "h-minima for various data types"

        data = np.array([[10,  11,  13,  14,  14,  15,  14,  14,  13,  11],
                         [11,  13,  15,  16,  16,  16,  16,  16,  15,  13],
                         [13,  15,  40,  40,  18,  18,  18,  60,  60,  15],
                         [14,  16,  40,  40,  19,  19,  19,  60,  60,  16],
                         [14,  16,  18,  19,  19,  19,  19,  19,  18,  16],
                         [15,  16,  18,  19,  19,  20,  19,  19,  18,  16],
                         [14,  16,  18,  19,  19,  19,  19,  19,  18,  16],
                         [14,  16,  80,  80,  19,  19,  19, 100, 100,  16],
                         [13,  15,  80,  80,  18,  18,  18, 100, 100,  15],
                         [11,  13,  15,  16,  16,  16,  16,  16,  15,  13]],
                        dtype=np.uint8)
        data = 100 - data
        expected_result = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 1, 1, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 1, 1, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 1, 1, 0, 0, 0, 1, 1, 0],
                                    [0, 0, 1, 1, 0, 0, 0, 1, 1, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
                                   dtype=np.uint8)
        for dtype in [np.uint8, np.uint64, np.int8, np.int64]:
            data = data.astype(dtype)
            out = extrema.h_minima(data, 40)

            error = diff(expected_result, out)
            assert error < eps
            assert out.dtype == expected_result.dtype

    def test_local_extrema_uniform(self):
        "local extrema tests for uniform arrays with various data types"

        data = np.full((7, 6), 42, dtype=np.uint8)

        expected_result = np.zeros((7, 6), dtype=np.uint8)

        for dtype in [np.uint8, np.uint64, np.int8, np.int64]:
            data = data.astype(dtype)

            # test for local maxima
            out = extrema.local_maxima(data)
            error = diff(expected_result, out)
            assert error < eps
            assert out.dtype == expected_result.dtype

            # test for local minima
            out = extrema.local_minima(data)
            error = diff(expected_result, out)
            assert error < eps
            assert out.dtype == expected_result.dtype

    def test_extrema_float(self):
        "specific tests for float type"
        data = np.array([[0.10, 0.11, 0.13, 0.14, 0.14, 0.15, 0.14,
                          0.14, 0.13, 0.11],
                         [0.11, 0.13, 0.15, 0.16, 0.16, 0.16, 0.16,
                          0.16, 0.15, 0.13],
                         [0.13, 0.15, 0.40, 0.40, 0.18, 0.18, 0.18,
                          0.60, 0.60, 0.15],
                         [0.14, 0.16, 0.40, 0.40, 0.19, 0.19, 0.19,
                          0.60, 0.60, 0.16],
                         [0.14, 0.16, 0.18, 0.19, 0.19, 0.19, 0.19,
                          0.19, 0.18, 0.16],
                         [0.15, 0.182, 0.18, 0.19, 0.204, 0.20, 0.19,
                          0.19, 0.18, 0.16],
                         [0.14, 0.16, 0.18, 0.19, 0.19, 0.19, 0.19,
                          0.19, 0.18, 0.16],
                         [0.14, 0.16, 0.80, 0.80, 0.19, 0.19, 0.19,
                          1.0,  1.0, 0.16],
                         [0.13, 0.15, 0.80, 0.80, 0.18, 0.18, 0.18,
                          1.0, 1.0, 0.15],
                         [0.11, 0.13, 0.15, 0.16, 0.16, 0.16, 0.16,
                          0.16, 0.15, 0.13]],
                        dtype=np.float32)
        inverted_data = 1.0 - data

        expected_result = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 1, 1, 0, 0, 0, 1, 1, 0],
                                    [0, 0, 1, 1, 0, 0, 0, 1, 1, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 1, 0, 0, 1, 0, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 1, 1, 0, 0, 0, 1, 1, 0],
                                    [0, 0, 1, 1, 0, 0, 0, 1, 1, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
                                   dtype=np.uint8)

        # test for local maxima with automatic step calculation
        out = extrema.local_maxima(data)
        error = diff(expected_result, out)
        assert error < eps

        # test for local minima with automatic step calculation
        out = extrema.local_minima(inverted_data)
        error = diff(expected_result, out)
        assert error < eps

        out = extrema.h_maxima(data, 0.003)
        expected_result = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 1, 1, 0, 0, 0, 1, 1, 0],
                                    [0, 0, 1, 1, 0, 0, 0, 1, 1, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 1, 1, 0, 0, 0, 1, 1, 0],
                                    [0, 0, 1, 1, 0, 0, 0, 1, 1, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
                                   dtype=np.uint8)

        error = diff(expected_result, out)
        assert error < eps

        out = extrema.h_minima(inverted_data, 0.003)
        error = diff(expected_result, out)
        assert error < eps

    def test_3d(self):
        """tests the detection of maxima in 3D."""
        img = np.zeros((8, 8, 8), dtype=np.uint8)
        local_maxima = np.zeros((8, 8, 8), dtype=np.uint8)

        # first maximum: only one pixel
        img[1, 1:3, 1:3] = 100
        img[2, 2, 2] = 200
        img[3, 1:3, 1:3] = 100
        local_maxima[2, 2, 2] = 1

        # second maximum: three pixels in z-direction
        img[5:8, 1, 1] = 200
        local_maxima[5:8, 1, 1] = 1

        # third: two maxima in 0 and 3.
        img[0, 5:8, 5:8] = 200
        img[1, 6, 6] = 100
        img[2, 5:7, 5:7] = 200
        img[0:3, 5:8, 5:8] += 50
        local_maxima[0, 5:8, 5:8] = 1
        local_maxima[2, 5:7, 5:7] = 1

        # four : one maximum in the corner of the square
        img[6:8, 6:8, 6:8] = 200
        img[7, 7, 7] = 255
        local_maxima[7, 7, 7] = 1

        se = ndi.generate_binary_structure(3, 1)
        out = extrema.local_maxima(img, se)

        error = diff(local_maxima, out)
        assert error < eps


if __name__ == "__main__":
    np.testing.run_module_suite()
