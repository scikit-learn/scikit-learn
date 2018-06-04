import os
import tempfile

import numpy as np
from skimage import novice
from skimage.novice._novice import (array_to_xy_origin, xy_to_array_origin,
                                    rgb_transpose)
from skimage import data_dir

from skimage._shared import testing
from skimage._shared.testing import (assert_equal, assert_allclose,
                                     expected_warnings)
from skimage._shared.utils import all_warnings


IMAGE_PATH = os.path.join(data_dir, "chelsea.png")
SMALL_IMAGE_PATH = os.path.join(data_dir, "block.png")


def _array_2d_to_RGBA(array):
    return np.tile(array[:, :, np.newaxis], (1, 1, 4))


def test_xy_to_array_origin():
    h, w = 3, 5
    array = np.arange(h * w).reshape(h, w, 1)
    out = xy_to_array_origin(array_to_xy_origin(array.copy()))
    assert np.allclose(out, array)


def test_pic_info():
    pic = novice.open(IMAGE_PATH)
    assert_equal(pic.format, "png")
    assert_equal(pic.path, os.path.abspath(IMAGE_PATH))
    assert_equal(pic.size, (451, 300))
    assert_equal(pic.width, 451)
    assert_equal(pic.height, 300)
    assert not pic.modified


def test_pixel_iteration():
    pic = novice.open(SMALL_IMAGE_PATH)
    num_pixels = sum(1 for p in pic)
    assert_equal(num_pixels, pic.width * pic.height)


def test_modify():
    pic = novice.open(SMALL_IMAGE_PATH)
    assert_equal(pic.modified, False)

    for p in pic:
        if p.x < (pic.width / 2):
            p.red /= 2
            p.green /= 2
            p.blue /= 2

    for p in pic:
        if p.x < (pic.width / 2):
            assert p.red <= 128
            assert p.green <= 128
            assert p.blue <= 128

    s = pic.size
    with all_warnings():  # precision loss
        pic.size = (pic.width / 2, pic.height / 2)
    assert_equal(pic.size, (int(s[0] / 2), int(s[1] / 2)))

    assert pic.modified
    assert pic.path is None


def test_pixel_rgb():
    pic = novice.Picture.from_size((3, 3), color=(10, 10, 10))
    pixel = pic[0, 0]
    pixel.rgb = np.arange(3)

    assert_equal(pixel.rgb, np.arange(3))
    for i, channel in enumerate((pixel.red, pixel.green, pixel.blue)):
        assert_equal(channel, i)

    pixel.red = 3
    pixel.green = 4
    pixel.blue = 5
    assert_equal(pixel.rgb, np.arange(3) + 3)

    for i, channel in enumerate((pixel.red, pixel.green, pixel.blue)):
        assert_equal(channel, i + 3)

    pixel.rgb = np.arange(4)
    assert_equal(pixel.rgb, np.arange(3))

    assert pic.array.dtype == np.uint8


def test_pixel_rgba():
    pic = novice.Picture.from_size((3, 3), color=(10, 10, 10))
    pixel = pic[0, 0]
    pixel.rgba = np.arange(4)

    assert_equal(pixel.rgba, np.arange(4))
    pixel_channels = (pixel.red, pixel.green, pixel.blue, pixel.alpha)
    for i, channel in enumerate(pixel_channels):
        assert_equal(channel, i)

    pixel.red = 3
    pixel.green = 4
    pixel.blue = 5
    pixel.alpha = 6
    assert_equal(pixel.rgba, np.arange(4) + 3)

    pixel_channels = (pixel.red, pixel.green, pixel.blue, pixel.alpha)
    for i, channel in enumerate(pixel_channels):
        assert_equal(channel, i + 3)


def test_pixel_rgb_float():
    pixel = novice.Picture.from_size((1, 1))[0, 0]
    pixel.rgb = (1.1, 1.1, 1.1)
    assert_equal(pixel.rgb, (1, 1, 1))


def test_pixel_rgba_float():
    pixel = novice.Picture.from_size((1, 1))[0, 0]
    pixel.rgba = (1.1, 1.1, 1.1, 1.1)
    assert_equal(pixel.rgba, (1, 1, 1, 1))


def test_modified_on_set():
    pic = novice.Picture(SMALL_IMAGE_PATH)
    pic[0, 0] = (1, 1, 1)
    assert pic.modified
    assert pic.path is None


def test_modified_on_set_pixel():
    data = np.zeros(shape=(10, 5, 3), dtype=np.uint8)
    pic = novice.Picture(array=data)

    pixel = pic[0, 0]
    pixel.green = 1
    assert pic.modified


def test_reset():
    pic = novice.Picture(SMALL_IMAGE_PATH)
    v = pic[0, 0]
    pic[0, 0] = (1, 1, 1)
    pic.reset()
    assert pic[0, 0] == v


def test_update_on_save():
    pic = novice.Picture(array=np.zeros((3, 3, 3)))
    # prevent attempting to save low-contrast image
    pic[0, 0] = (255, 255, 255)

    with all_warnings():  # precision loss
        pic.size = (6, 6)
    assert pic.modified
    assert pic.path is None

    fd, filename = tempfile.mkstemp(suffix=".png")
    os.close(fd)
    try:
        pic.save(filename)

        assert not pic.modified
        assert_equal(pic.path, os.path.abspath(filename))
        assert_equal(pic.format, "png")
    finally:
        os.unlink(filename)


def test_save_with_alpha_channel():
    # create an image with an alpha channel
    pic = novice.Picture(array=np.zeros((3, 3, 4)))

    fd, filename = tempfile.mkstemp(suffix=".png")
    os.close(fd)
    with expected_warnings(['is a low contrast']):
        pic.save(filename)
    os.unlink(filename)


def test_indexing():
    array = 128 * np.ones((10, 10, 3), dtype=np.uint8)
    pic = novice.Picture(array=array)

    pic[0:5, 0:5] = (0, 0, 0)
    for p in pic:
        if (p.x < 5) and (p.y < 5):
            assert_equal(p.rgb, (0, 0, 0))
            assert_equal(p.red, 0)
            assert_equal(p.green, 0)
            assert_equal(p.blue, 0)

    pic[:5, :5] = (255, 255, 255)
    for p in pic:
        if (p.x < 5) and (p.y < 5):
            assert_equal(p.rgb, (255, 255, 255))
            assert_equal(p.red, 255)
            assert_equal(p.green, 255)
            assert_equal(p.blue, 255)

    pic[5:pic.width, 5:pic.height] = (255, 0, 255)
    for p in pic:
        if (p.x >= 5) and (p.y >= 5):
            assert_equal(p.rgb, (255, 0, 255))
            assert_equal(p.red, 255)
            assert_equal(p.green, 0)
            assert_equal(p.blue, 255)

    pic[5:, 5:] = (0, 0, 255)
    for p in pic:
        if (p.x >= 5) and (p.y >= 5):
            assert_equal(p.rgb, (0, 0, 255))
            assert_equal(p.red, 0)
            assert_equal(p.green, 0)
            assert_equal(p.blue, 255)


def test_picture_slice():
    array = _array_2d_to_RGBA(np.arange(0, 10)[np.newaxis, :])
    pic = novice.Picture(array=array)

    x_slice = slice(3, 8)
    subpic = pic[:, x_slice]
    assert_allclose(subpic.array, array[x_slice, :])


def test_move_slice():
    h, w = 3, 12
    array = _array_2d_to_RGBA(np.linspace(0, 255, h * w).reshape(h, w))
    array = array.astype(np.uint8)

    pic = novice.Picture(array=array)
    pic_orig = novice.Picture(array=array.copy())

    # Move left cut of image to the right side.
    cut = 5
    rest = pic.width - cut
    temp = pic[:cut, :]
    temp.array = temp.array.copy()
    pic[:rest, :] = pic[cut:, :]
    pic[rest:, :] = temp

    assert pic[rest:, :] == pic_orig[:cut, :]
    assert pic[:rest, :] == pic_orig[cut:, :]


def test_negative_index():
    n = 10
    array = _array_2d_to_RGBA(np.arange(0, n)[np.newaxis, :])
    # Test both x and y indices.
    pic = novice.Picture(array=array)
    assert pic[-1, 0] == pic[n - 1, 0]
    pic = novice.Picture(array=rgb_transpose(array))
    assert pic[0, -1] == pic[0, n - 1]


def test_negative_slice():
    n = 10
    array = _array_2d_to_RGBA(np.arange(0, n)[np.newaxis, :])
    # Test both x and y slices.
    pic = novice.Picture(array=array)
    assert pic[-3:, 0] == pic[n - 3:, 0]
    pic = novice.Picture(array=rgb_transpose(array))
    assert pic[0, -3:] == pic[0, n - 3:]


def test_getitem_with_step():
    h, w = 5, 5
    array = _array_2d_to_RGBA(np.linspace(0, 255, h * w).reshape(h, w))
    pic = novice.Picture(array=array)
    sliced_pic = pic[::2, ::2]
    assert sliced_pic == novice.Picture(array=array[::2, ::2])


def test_1d_getitem_raises():
    pic = novice.Picture.from_size((1, 1))
    with testing.raises(IndexError):
        pic[1]


def test_3d_getitem_raises():
    pic = novice.Picture.from_size((1, 1))
    with testing.raises(IndexError):
        pic[1, 2, 3]


def test_1d_setitem_raises():
    pic = novice.Picture.from_size((1, 1))
    with testing.raises(IndexError):
        pic[1] = 0


def test_3d_setitem_raises():
    pic = novice.Picture.from_size((1, 1))
    with testing.raises(IndexError):
        pic[1, 2, 3] = 0


def test_out_of_bounds_indexing():
    pic = novice.open(SMALL_IMAGE_PATH)
    with testing.raises(IndexError):
        pic[pic.width, pic.height]


def test_pixel_rgb_raises():
    pixel = novice.Picture.from_size((1, 1))[0, 0]
    with testing.raises(ValueError):
        pixel.rgb = (-1, -1, -1)


def test_pixel_red_raises():
    pixel = novice.Picture.from_size((1, 1))[0, 0]
    with testing.raises(ValueError):
        pixel.red = 256


def test_pixel_green_raises():
    pixel = novice.Picture.from_size((1, 1))[0, 0]
    with testing.raises(ValueError):
        pixel.green = 256


def test_pixel_blue_raises():
    pixel = novice.Picture.from_size((1, 1))[0, 0]
    with testing.raises(ValueError):
        pixel.blue = 256


def test_pixel_alpha_raises():
    pixel = novice.Picture.from_size((1, 1))[0, 0]
    with testing.raises(ValueError):
        pixel.alpha = 256
