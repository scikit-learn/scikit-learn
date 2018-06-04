import numpy as np
from skimage._shared.testing import assert_allclose

from skimage.feature.register_translation import (register_translation,
                                                  _upsampled_dft)
from skimage.data import camera, binary_blobs
from scipy.ndimage import fourier_shift
from skimage import img_as_float
from skimage._shared import testing


def test_correlation():
    reference_image = np.fft.fftn(camera())
    shift = (-7, 12)
    shifted_image = fourier_shift(reference_image, shift)

    # pixel precision
    result, error, diffphase = register_translation(reference_image,
                                                    shifted_image,
                                                    space="fourier")
    assert_allclose(result[:2], -np.array(shift))


def test_subpixel_precision():
    reference_image = np.fft.fftn(camera())
    subpixel_shift = (-2.4, 1.32)
    shifted_image = fourier_shift(reference_image, subpixel_shift)

    # subpixel precision
    result, error, diffphase = register_translation(reference_image,
                                                    shifted_image, 100,
                                                    space="fourier")
    assert_allclose(result[:2], -np.array(subpixel_shift), atol=0.05)


def test_real_input():
    reference_image = camera()
    subpixel_shift = (-2.4, 1.32)
    shifted_image = fourier_shift(np.fft.fftn(reference_image), subpixel_shift)
    shifted_image = np.fft.ifftn(shifted_image)

    # subpixel precision
    result, error, diffphase = register_translation(reference_image,
                                                    shifted_image, 100)
    assert_allclose(result[:2], -np.array(subpixel_shift), atol=0.05)


def test_size_one_dimension_input():
    # take a strip of the input image
    reference_image = np.fft.fftn(camera()[:, 15]).reshape((-1, 1))
    subpixel_shift = (-2.4, 4)
    shifted_image = fourier_shift(reference_image, subpixel_shift)

    # subpixel precision
    result, error, diffphase = register_translation(reference_image,
                                                    shifted_image, 100,
                                                    space="fourier")
    assert_allclose(result[:2], -np.array((-2.4, 0)), atol=0.05)


def test_3d_input():
    phantom = img_as_float(binary_blobs(length=32, n_dim=3))
    reference_image = np.fft.fftn(phantom)
    shift = (-2., 1., 5.)
    shifted_image = fourier_shift(reference_image, shift)

    result, error, diffphase = register_translation(reference_image,
                                                    shifted_image,
                                                    space="fourier")
    assert_allclose(result, -np.array(shift), atol=0.05)
    # subpixel precision not available for 3-D data
    subpixel_shift = (-2.3, 1., 5.)
    shifted_image = fourier_shift(reference_image, subpixel_shift)
    result, error, diffphase = register_translation(reference_image,
                                                    shifted_image,
                                                    space="fourier")
    assert_allclose(result, -np.array(shift), atol=0.5)
    with testing.raises(NotImplementedError):
        register_translation(
            reference_image,
            shifted_image, upsample_factor=100,
            space="fourier")


def test_unknown_space_input():
    image = np.ones((5, 5))
    with testing.raises(ValueError):
        register_translation(
            image, image,
            space="frank")


def test_wrong_input():
    # Dimensionality mismatch
    image = np.ones((5, 5, 1))
    template = np.ones((5, 5))
    with testing.raises(ValueError):
        register_translation(template, image)

    # Greater than 2 dimensions does not support subpixel precision
    #   (TODO: should support 3D at some point.)
    image = np.ones((5, 5, 5))
    template = np.ones((5, 5, 5))
    with testing.raises(NotImplementedError):
        register_translation(template, image, 2)

    # Size mismatch
    image = np.ones((5, 5))
    template = np.ones((4, 4))
    with testing.raises(ValueError):
        register_translation(template, image)


def test_mismatch_upsampled_region_size():
    with testing.raises(ValueError):
        _upsampled_dft(
            np.ones((4, 4)),
            upsampled_region_size=[3, 2, 1, 4])


def test_mismatch_offsets_size():
    with testing.raises(ValueError):
        _upsampled_dft(np.ones((4, 4)), 3,
                       axis_offsets=[3, 2, 1, 4])
