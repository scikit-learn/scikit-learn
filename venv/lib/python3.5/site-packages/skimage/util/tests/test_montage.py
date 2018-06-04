from skimage._shared import testing
from skimage._shared.testing import (assert_equal, assert_array_equal,
                                     expected_warnings)

import numpy as np
from skimage.util import montage
from skimage.util import montage2d as montage2d_deprecated


def montage2d(*args, **kwargs):
    with expected_warnings(['deprecated']):
        return montage2d_deprecated(*args, **kwargs)


def test_montage2d_simple():
    n_images = 3
    height, width = 2, 3,
    arr_in = np.arange(n_images * height * width, dtype='float')
    arr_in = arr_in.reshape(n_images, height, width)

    arr_out = montage2d(arr_in)

    gt = np.array(
        [[  0. ,   1. ,   2. ,   6. ,   7. ,   8. ],
         [  3. ,   4. ,   5. ,   9. ,  10. ,  11. ],
         [ 12. ,  13. ,  14. ,   8.5,   8.5,   8.5],
         [ 15. ,  16. ,  17. ,   8.5,   8.5,   8.5]]
    )

    assert_array_equal(arr_out, gt)


def test_montage2d_fill():
    n_images = 3
    height, width = 2, 3,
    arr_in = np.arange(n_images * height * width)
    arr_in = arr_in.reshape(n_images, height, width)

    arr_out = montage2d(arr_in, fill=0)

    gt = np.array(
        [[  0. ,   1. ,   2. ,   6. ,   7. ,   8. ],
         [  3. ,   4. ,   5. ,   9. ,  10. ,  11. ],
         [ 12. ,  13. ,  14. ,   0. ,   0. ,   0. ],
         [ 15. ,  16. ,  17. ,   0. ,   0. ,   0. ]]
    )

    assert_array_equal(arr_out, gt)


def test_montage2d_shape():
    n_images = 15
    height, width = 11, 7
    arr_in = np.arange(n_images * height * width)
    arr_in = arr_in.reshape(n_images, height, width)

    alpha = int(np.ceil(np.sqrt(n_images)))

    arr_out = montage2d(arr_in)
    assert_equal(arr_out.shape, (alpha * height, alpha * width))


def test_montage2d_grid_shape():
    n_images = 6
    height, width = 2, 2
    arr_in = np.arange(n_images * height * width, dtype=np.float32)
    arr_in = arr_in.reshape(n_images, height, width)
    arr_out = montage2d(arr_in, grid_shape=(3, 2))
    correct_arr_out = np.array(
	[[  0.,   1.,   4.,   5.],
	 [  2.,   3.,   6.,   7.],
	 [  8.,   9.,  12.,  13.],
	 [ 10.,  11.,  14.,  15.],
	 [ 16.,  17.,  20.,  21.],
	 [ 18.,  19.,  22.,  23.]]
    )
    assert_array_equal(arr_out, correct_arr_out)


def test_montage2d_rescale_intensity():
    n_images = 4
    height, width = 3, 3
    arr_in = np.arange(n_images * height * width, dtype=np.float32)
    arr_in = arr_in.reshape(n_images, height, width)

    arr_out = montage2d(arr_in, rescale_intensity=True)

    gt = np.array(
        [[ 0.   ,  0.125,  0.25 ,  0.   ,  0.125,  0.25 ],
         [ 0.375,  0.5  ,  0.625,  0.375,  0.5  ,  0.625],
         [ 0.75 ,  0.875,  1.   ,  0.75 ,  0.875,  1.   ],
         [ 0.   ,  0.125,  0.25 ,  0.   ,  0.125,  0.25 ],
         [ 0.375,  0.5  ,  0.625,  0.375,  0.5  ,  0.625],
         [ 0.75 ,  0.875,  1.   ,  0.75 ,  0.875,  1.   ]]
        )

    assert_equal(arr_out.min(), 0.0)
    assert_equal(arr_out.max(), 1.0)
    assert_array_equal(arr_out, gt)


def test_montage2d_simple_padding():
    n_images = 2
    height, width = 2, 2,
    arr_in = np.arange(n_images * height * width)
    arr_in = arr_in.reshape(n_images, height, width)

    arr_out = montage2d(arr_in, padding_width=1)

    gt = np.array(
        [[0, 1, 0, 4, 5, 0],
         [2, 3, 0, 6, 7, 0],
         [0, 0, 0, 0, 0, 0],
         [3, 3, 3, 3, 3, 3],
         [3, 3, 3, 3, 3, 3],
         [3, 3, 3, 3, 3, 3]]
    )

    assert_array_equal(arr_out, gt)


def test_montage2d_error_ndim():
    arr_error = np.random.randn(1, 2, 3, 4)
    with testing.raises(AssertionError):
        montage2d(arr_error)


def test_montage_simple_gray():
    n_images, n_rows, n_cols = 3, 2, 3
    arr_in = np.arange(n_images * n_rows * n_cols, dtype=np.float)
    arr_in = arr_in.reshape(n_images, n_rows, n_cols)

    arr_out = montage(arr_in)
    arr_ref = np.array(
        [[  0. ,   1. ,   2. ,   6. ,   7. ,   8. ],
         [  3. ,   4. ,   5. ,   9. ,  10. ,  11. ],
         [ 12. ,  13. ,  14. ,   8.5,   8.5,   8.5],
         [ 15. ,  16. ,  17. ,   8.5,   8.5,   8.5]]
    )
    assert_array_equal(arr_out, arr_ref)


def test_montage_simple_rgb():
    n_images, n_rows, n_cols, n_channels = 2, 2, 2, 2
    arr_in = np.arange(n_images * n_rows * n_cols * n_channels, dtype=np.float)
    arr_in = arr_in.reshape(n_images, n_rows, n_cols, n_channels)

    arr_out = montage(arr_in, multichannel=True)
    arr_ref = np.array(
        [[[ 0,  1],
          [ 2,  3],
          [ 8,  9],
          [10, 11]],
         [[ 4,  5],
          [ 6,  7],
          [12, 13],
          [14, 15]],
         [[ 7,  8],
          [ 7,  8],
          [ 7,  8],
          [ 7,  8]],
         [[ 7,  8],
          [ 7,  8],
          [ 7,  8],
          [ 7,  8]]]
        )
    assert_array_equal(arr_out, arr_ref)


def test_montage_fill_gray():
    n_images, n_rows, n_cols = 3, 2, 3
    arr_in = np.arange(n_images*n_rows*n_cols, dtype=np.float)
    arr_in = arr_in.reshape(n_images, n_rows, n_cols)

    arr_out = montage(arr_in, fill=0)
    arr_ref = np.array(
        [[  0. ,   1. ,   2. ,   6. ,   7. ,   8. ],
         [  3. ,   4. ,   5. ,   9. ,  10. ,  11. ],
         [ 12. ,  13. ,  14. ,   0. ,   0. ,   0. ],
         [ 15. ,  16. ,  17. ,   0. ,   0. ,   0. ]]
    )
    assert_array_equal(arr_out, arr_ref)


def test_montage_grid_default_gray():
    n_images, n_rows, n_cols = 15, 11, 7
    arr_in = np.arange(n_images * n_rows * n_cols, dtype=np.float)
    arr_in = arr_in.reshape(n_images, n_rows, n_cols)

    n_tiles = int(np.ceil(np.sqrt(n_images)))
    arr_out = montage(arr_in)
    assert_equal(arr_out.shape, (n_tiles * n_rows, n_tiles * n_cols))


def test_montage_grid_custom_gray():
    n_images, n_rows, n_cols = 6, 2, 2
    arr_in = np.arange(n_images * n_rows * n_cols, dtype=np.float32)
    arr_in = arr_in.reshape(n_images, n_rows, n_cols)

    arr_out = montage(arr_in, grid_shape=(3, 2))
    arr_ref = np.array(
	[[  0.,   1.,   4.,   5.],
         [  2.,   3.,   6.,   7.],
         [  8.,   9.,  12.,  13.],
         [ 10.,  11.,  14.,  15.],
         [ 16.,  17.,  20.,  21.],
         [ 18.,  19.,  22.,  23.]]
    )
    assert_array_equal(arr_out, arr_ref)


def test_montage_rescale_intensity_gray():
    n_images, n_rows, n_cols = 4, 3, 3
    arr_in = np.arange(n_images * n_rows * n_cols, dtype=np.float32)
    arr_in = arr_in.reshape(n_images, n_rows, n_cols)

    arr_out = montage(arr_in, rescale_intensity=True)
    arr_ref = np.array(
        [[ 0.   ,  0.125,  0.25 ,  0.   ,  0.125,  0.25 ],
         [ 0.375,  0.5  ,  0.625,  0.375,  0.5  ,  0.625],
         [ 0.75 ,  0.875,  1.   ,  0.75 ,  0.875,  1.   ],
         [ 0.   ,  0.125,  0.25 ,  0.   ,  0.125,  0.25 ],
         [ 0.375,  0.5  ,  0.625,  0.375,  0.5  ,  0.625],
         [ 0.75 ,  0.875,  1.   ,  0.75 ,  0.875,  1.   ]]
    )
    assert_equal(arr_out.min(), 0.0)
    assert_equal(arr_out.max(), 1.0)
    assert_array_equal(arr_out, arr_ref)


def test_montage_simple_padding_gray():
    n_images, n_rows, n_cols = 2, 2, 2
    arr_in = np.arange(n_images * n_rows * n_cols)
    arr_in = arr_in.reshape(n_images, n_rows, n_cols)

    arr_out = montage(arr_in, padding_width=1)
    arr_ref = np.array(
        [[3, 3, 3, 3, 3, 3, 3],
         [3, 0, 1, 3, 4, 5, 3],
         [3, 2, 3, 3, 6, 7, 3],
         [3, 3, 3, 3, 3, 3, 3],
         [3, 3, 3, 3, 3, 3, 3],
         [3, 3, 3, 3, 3, 3, 3],
         [3, 3, 3, 3, 3, 3, 3]]
    )
    assert_array_equal(arr_out, arr_ref)


def test_error_ndim():
    arr_error = np.random.randn(1, 2)
    with testing.raises(ValueError):
        montage(arr_error)

    arr_error = np.random.randn(1, 2, 3, 4)
    with testing.raises(ValueError):
        montage(arr_error)

    arr_error = np.random.randn(1, 2, 3)
    with testing.raises(ValueError):
        montage(arr_error, multichannel=True)

    arr_error = np.random.randn(1, 2, 3, 4, 5)
    with testing.raises(ValueError):
        montage(arr_error, multichannel=True)
