import numpy as np
from numpy import testing

from skimage import data, color
from skimage.util import img_as_bool
from skimage.morphology import binary, grey, selem
from scipy import ndimage as ndi

import pytest

img = color.rgb2gray(data.astronaut())
bw_img = img > 100 / 255.


def test_non_square_image():
    strel = selem.square(3)
    binary_res = binary.binary_erosion(bw_img[:100, :200], strel)
    grey_res = img_as_bool(grey.erosion(bw_img[:100, :200], strel))
    testing.assert_array_equal(binary_res, grey_res)


def test_binary_erosion():
    strel = selem.square(3)
    binary_res = binary.binary_erosion(bw_img, strel)
    grey_res = img_as_bool(grey.erosion(bw_img, strel))
    testing.assert_array_equal(binary_res, grey_res)


def test_binary_dilation():
    strel = selem.square(3)
    binary_res = binary.binary_dilation(bw_img, strel)
    grey_res = img_as_bool(grey.dilation(bw_img, strel))
    testing.assert_array_equal(binary_res, grey_res)


def test_binary_closing():
    strel = selem.square(3)
    binary_res = binary.binary_closing(bw_img, strel)
    grey_res = img_as_bool(grey.closing(bw_img, strel))
    testing.assert_array_equal(binary_res, grey_res)


def test_binary_opening():
    strel = selem.square(3)
    binary_res = binary.binary_opening(bw_img, strel)
    grey_res = img_as_bool(grey.opening(bw_img, strel))
    testing.assert_array_equal(binary_res, grey_res)


def test_selem_overflow():
    strel = np.ones((17, 17), dtype=np.uint8)
    img = np.zeros((20, 20), dtype=bool)
    img[2:19, 2:19] = True
    binary_res = binary.binary_erosion(img, strel)
    grey_res = img_as_bool(grey.erosion(img, strel))
    testing.assert_array_equal(binary_res, grey_res)


def test_out_argument():
    for func in (binary.binary_erosion, binary.binary_dilation):
        strel = np.ones((3, 3), dtype=np.uint8)
        img = np.ones((10, 10))
        out = np.zeros_like(img)
        out_saved = out.copy()
        func(img, strel, out=out)
        testing.assert_(np.any(out != out_saved))
        testing.assert_array_equal(out, func(img, strel))


binary_functions = [binary.binary_erosion, binary.binary_dilation,
                    binary.binary_opening, binary.binary_closing]


@pytest.mark.parametrize("function", binary_functions)
def test_default_selem(function):
    strel = selem.diamond(radius=1)
    image = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                      [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
                      [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
                      [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
                      [0, 0, 1, 1, 1, 0, 0, 1, 0, 0],
                      [0, 0, 1, 1, 1, 0, 0, 1, 0, 0],
                      [0, 0, 1, 1, 1, 0, 0, 1, 0, 0],
                      [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
                      [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
                      [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
                      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]], np.uint8)
    im_expected = function(image, strel)
    im_test = function(image)
    testing.assert_array_equal(im_expected, im_test)

def test_3d_fallback_default_selem():
    # 3x3x3 cube inside a 7x7x7 image:
    image = np.zeros((7, 7, 7), np.bool)
    image[2:-2, 2:-2, 2:-2] = 1

    opened = binary.binary_opening(image)

    # expect a "hyper-cross" centered in the 5x5x5:
    image_expected = np.zeros((7, 7, 7), dtype=bool)
    image_expected[2:5, 2:5, 2:5] = ndi.generate_binary_structure(3, 1)
    testing.assert_array_equal(opened, image_expected)


binary_3d_fallback_functions = [binary.binary_opening, binary.binary_closing]


@pytest.mark.parametrize("function", binary_3d_fallback_functions)
def test_3d_fallback_cube_selem(function):
    # 3x3x3 cube inside a 7x7x7 image:
    image = np.zeros((7, 7, 7), np.bool)
    image[2:-2, 2:-2, 2:-2] = 1

    cube = np.ones((3, 3, 3), dtype=np.uint8)

    new_image = function(image, cube)
    testing.assert_array_equal(new_image, image)

def test_2d_ndimage_equivalence():
    image = np.zeros((9, 9), np.uint16)
    image[2:-2, 2:-2] = 2**14
    image[3:-3, 3:-3] = 2**15
    image[4, 4] = 2**16-1

    bin_opened = binary.binary_opening(image)
    bin_closed = binary.binary_closing(image)

    selem = ndi.generate_binary_structure(2, 1)
    ndimage_opened = ndi.binary_opening(image, structure=selem)
    ndimage_closed = ndi.binary_closing(image, structure=selem)

    testing.assert_array_equal(bin_opened, ndimage_opened)
    testing.assert_array_equal(bin_closed, ndimage_closed)

def test_binary_output_2d():
    image = np.zeros((9, 9), np.uint16)
    image[2:-2, 2:-2] = 2**14
    image[3:-3, 3:-3] = 2**15
    image[4, 4] = 2**16-1

    bin_opened = binary.binary_opening(image)
    bin_closed = binary.binary_closing(image)

    int_opened = np.empty_like(image, dtype=np.uint8)
    int_closed = np.empty_like(image, dtype=np.uint8)
    binary.binary_opening(image, out=int_opened)
    binary.binary_closing(image, out=int_closed)

    testing.assert_equal(bin_opened.dtype, np.bool)
    testing.assert_equal(bin_closed.dtype, np.bool)

    testing.assert_equal(int_opened.dtype, np.uint8)
    testing.assert_equal(int_closed.dtype, np.uint8)

def test_binary_output_3d():
    image = np.zeros((9, 9, 9), np.uint16)
    image[2:-2, 2:-2, 2:-2] = 2**14
    image[3:-3, 3:-3, 3:-3] = 2**15
    image[4, 4, 4] = 2**16-1

    bin_opened = binary.binary_opening(image)
    bin_closed = binary.binary_closing(image)

    int_opened = np.empty_like(image, dtype=np.uint8)
    int_closed = np.empty_like(image, dtype=np.uint8)
    binary.binary_opening(image, out=int_opened)
    binary.binary_closing(image, out=int_closed)

    testing.assert_equal(bin_opened.dtype, np.bool)
    testing.assert_equal(bin_closed.dtype, np.bool)

    testing.assert_equal(int_opened.dtype, np.uint8)
    testing.assert_equal(int_closed.dtype, np.uint8)

if __name__ == '__main__':
    testing.run_module_suite()
