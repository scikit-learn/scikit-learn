import math
import numpy as np
from skimage import data
from skimage.transform import pyramids

from skimage._shared import testing
from skimage._shared.testing import assert_array_equal, assert_, assert_equal
from skimage._shared._warnings import expected_warnings


image = data.astronaut()
image_gray = image[..., 0]


def test_pyramid_reduce_rgb():
    rows, cols, dim = image.shape
    with expected_warnings(['The default multichannel']):
        out = pyramids.pyramid_reduce(image, downscale=2)
    assert_array_equal(out.shape, (rows / 2, cols / 2, dim))


def test_pyramid_reduce_gray():
    rows, cols = image_gray.shape
    with expected_warnings(['The default multichannel']):
        out = pyramids.pyramid_reduce(image_gray, downscale=2)
    assert_array_equal(out.shape, (rows / 2, cols / 2))


def test_pyramid_reduce_nd():
    for ndim in [1, 2, 3, 4]:
        img = np.random.randn(*((8, ) * ndim))
        out = pyramids.pyramid_reduce(img, downscale=2,
                                      multichannel=False)
        expected_shape = np.asarray(img.shape) / 2
        assert_array_equal(out.shape, expected_shape)


def test_pyramid_expand_rgb():
    rows, cols, dim = image.shape
    with expected_warnings(['The default multichannel']):
        out = pyramids.pyramid_expand(image, upscale=2)
    assert_array_equal(out.shape, (rows * 2, cols * 2, dim))


def test_pyramid_expand_gray():
    rows, cols = image_gray.shape
    with expected_warnings(['The default multichannel']):
        out = pyramids.pyramid_expand(image_gray, upscale=2)
    assert_array_equal(out.shape, (rows * 2, cols * 2))


def test_pyramid_expand_nd():
    for ndim in [1, 2, 3, 4]:
        img = np.random.randn(*((4, ) * ndim))
        out = pyramids.pyramid_expand(img, upscale=2,
                                      multichannel=False)
        expected_shape = np.asarray(img.shape) * 2
        assert_array_equal(out.shape, expected_shape)


def test_build_gaussian_pyramid_rgb():
    rows, cols, dim = image.shape
    with expected_warnings(['The default multichannel']):
        pyramid = pyramids.pyramid_gaussian(image, downscale=2)
        for layer, out in enumerate(pyramid):
            layer_shape = (rows / 2 ** layer, cols / 2 ** layer, dim)
            assert_array_equal(out.shape, layer_shape)


def test_build_gaussian_pyramid_gray():
    rows, cols = image_gray.shape
    with expected_warnings(['The default multichannel']):
        pyramid = pyramids.pyramid_gaussian(image_gray, downscale=2)
        for layer, out in enumerate(pyramid):
            layer_shape = (rows / 2 ** layer, cols / 2 ** layer)
            assert_array_equal(out.shape, layer_shape)


def test_build_gaussian_pyramid_nd():
    for ndim in [1, 2, 3, 4]:
        img = np.random.randn(*((8, ) * ndim))
        original_shape = np.asarray(img.shape)
        pyramid = pyramids.pyramid_gaussian(img, downscale=2,
                                            multichannel=False)
        for layer, out in enumerate(pyramid):
            layer_shape = original_shape / 2 ** layer
            assert_array_equal(out.shape, layer_shape)


def test_build_laplacian_pyramid_rgb():
    rows, cols, dim = image.shape
    with expected_warnings(['The default multichannel']):
        pyramid = pyramids.pyramid_laplacian(image, downscale=2)
        for layer, out in enumerate(pyramid):
            layer_shape = (rows / 2 ** layer, cols / 2 ** layer, dim)
            assert_array_equal(out.shape, layer_shape)


def test_build_laplacian_pyramid_nd():
    for ndim in [1, 2, 3, 4]:
        img = np.random.randn(*(16, )*ndim)
        original_shape = np.asarray(img.shape)
        pyramid = pyramids.pyramid_laplacian(img, downscale=2,
                                             multichannel=False)
        for layer, out in enumerate(pyramid):
            print(out.shape)
            layer_shape = original_shape / 2 ** layer
            assert_array_equal(out.shape, layer_shape)


def test_laplacian_pyramid_max_layers():
    for downscale in [2, 3, 5, 7]:
        img = np.random.randn(32, 8)
        pyramid = pyramids.pyramid_laplacian(img, downscale=downscale,
                                             multichannel=False)
        max_layer = int(np.ceil(math.log(np.max(img.shape), downscale)))
        for layer, out in enumerate(pyramid):
            if layer < max_layer:
                # should not reach all axes as size 1 prior to final level
                assert_(np.max(out.shape) > 1)

        # total number of images is max_layer + 1
        assert_equal(max_layer, layer)

        # final layer should be size 1 on all axes
        assert_array_equal((out.shape), (1, 1))


def test_check_factor():
    with testing.raises(ValueError):
        pyramids._check_factor(0.99)
    with testing.raises(ValueError):
        pyramids._check_factor(- 2)
