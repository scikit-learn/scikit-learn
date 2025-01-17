from skimage._shared import testing
from skimage._shared.testing import assert_equal, assert_array_equal

import numpy as np
from skimage.util import montage

# TODO: when minimum numpy dependency is 1.25 use:
# np..exceptions.AxisError instead of AxisError
# and remove this try-except
try:
    from numpy import AxisError
except ImportError:
    from numpy.exceptions import AxisError


def test_montage_simple_gray():
    n_images, n_rows, n_cols = 3, 2, 3
    arr_in = np.arange(n_images * n_rows * n_cols, dtype=float)
    arr_in = arr_in.reshape(n_images, n_rows, n_cols)

    arr_out = montage(arr_in)
    arr_ref = np.array(
        [
            [0.0, 1.0, 2.0, 6.0, 7.0, 8.0],
            [3.0, 4.0, 5.0, 9.0, 10.0, 11.0],
            [12.0, 13.0, 14.0, 8.5, 8.5, 8.5],
            [15.0, 16.0, 17.0, 8.5, 8.5, 8.5],
        ]
    )
    assert_array_equal(arr_out, arr_ref)


def test_montage_simple_rgb():
    n_images, n_rows, n_cols, n_channels = 2, 2, 2, 2
    arr_in = np.arange(
        n_images * n_rows * n_cols * n_channels,
        dtype=float,
    )
    arr_in = arr_in.reshape(n_images, n_rows, n_cols, n_channels)

    arr_out = montage(arr_in, channel_axis=-1)
    arr_ref = np.array(
        [
            [[0, 1], [2, 3], [8, 9], [10, 11]],
            [[4, 5], [6, 7], [12, 13], [14, 15]],
            [[7, 8], [7, 8], [7, 8], [7, 8]],
            [[7, 8], [7, 8], [7, 8], [7, 8]],
        ]
    )
    assert_array_equal(arr_out, arr_ref)


@testing.parametrize('channel_axis', (0, 1, 2, 3, -1, -2, -3, -4))
def test_montage_simple_rgb_channel_axes(channel_axis):
    n_images, n_rows, n_cols, n_channels = 2, 2, 2, 2
    arr_in = np.arange(n_images * n_rows * n_cols * n_channels, dtype=float)
    arr_in = arr_in.reshape(n_images, n_rows, n_cols, n_channels)

    # place channels at the desired location
    arr_in = np.moveaxis(arr_in, -1, channel_axis)

    arr_out = montage(arr_in, channel_axis=channel_axis)
    arr_ref = np.array(
        [
            [[0, 1], [2, 3], [8, 9], [10, 11]],
            [[4, 5], [6, 7], [12, 13], [14, 15]],
            [[7, 8], [7, 8], [7, 8], [7, 8]],
            [[7, 8], [7, 8], [7, 8], [7, 8]],
        ]
    )
    assert_array_equal(arr_out, arr_ref)


@testing.parametrize('channel_axis', (4, -5))
def test_montage_invalid_channel_axes(channel_axis):
    arr_in = np.arange(16, dtype=float).reshape(2, 2, 2, 2)
    with testing.raises(AxisError):
        montage(arr_in, channel_axis=channel_axis)


def test_montage_fill_gray():
    n_images, n_rows, n_cols = 3, 2, 3
    arr_in = np.arange(n_images * n_rows * n_cols, dtype=float)
    arr_in = arr_in.reshape(n_images, n_rows, n_cols)

    arr_out = montage(arr_in, fill=0)
    arr_ref = np.array(
        [
            [0.0, 1.0, 2.0, 6.0, 7.0, 8.0],
            [3.0, 4.0, 5.0, 9.0, 10.0, 11.0],
            [12.0, 13.0, 14.0, 0.0, 0.0, 0.0],
            [15.0, 16.0, 17.0, 0.0, 0.0, 0.0],
        ]
    )
    assert_array_equal(arr_out, arr_ref)


def test_montage_grid_default_gray():
    n_images, n_rows, n_cols = 15, 11, 7
    arr_in = np.arange(n_images * n_rows * n_cols, dtype=float)
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
        [
            [0.0, 1.0, 4.0, 5.0],
            [2.0, 3.0, 6.0, 7.0],
            [8.0, 9.0, 12.0, 13.0],
            [10.0, 11.0, 14.0, 15.0],
            [16.0, 17.0, 20.0, 21.0],
            [18.0, 19.0, 22.0, 23.0],
        ]
    )
    assert_array_equal(arr_out, arr_ref)


def test_montage_rescale_intensity_gray():
    n_images, n_rows, n_cols = 4, 3, 3
    arr_in = np.arange(n_images * n_rows * n_cols, dtype=np.float32)
    arr_in = arr_in.reshape(n_images, n_rows, n_cols)

    arr_out = montage(arr_in, rescale_intensity=True)
    arr_ref = np.array(
        [
            [0.0, 0.125, 0.25, 0.0, 0.125, 0.25],
            [0.375, 0.5, 0.625, 0.375, 0.5, 0.625],
            [0.75, 0.875, 1.0, 0.75, 0.875, 1.0],
            [0.0, 0.125, 0.25, 0.0, 0.125, 0.25],
            [0.375, 0.5, 0.625, 0.375, 0.5, 0.625],
            [0.75, 0.875, 1.0, 0.75, 0.875, 1.0],
        ]
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
        [
            [3, 3, 3, 3, 3, 3, 3],
            [3, 0, 1, 3, 4, 5, 3],
            [3, 2, 3, 3, 6, 7, 3],
            [3, 3, 3, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 3, 3],
        ]
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
        montage(arr_error, channel_axis=-1)

    arr_error = np.random.randn(1, 2, 3, 4, 5)
    with testing.raises(ValueError):
        montage(arr_error, channel_axis=-1)
