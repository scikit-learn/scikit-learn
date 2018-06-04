import warnings

import numpy as np
import skimage
from skimage import data
from skimage import exposure
from skimage.exposure.exposure import intensity_range
from skimage.color import rgb2gray
from skimage.util.dtype import dtype_range

from skimage._shared._warnings import expected_warnings
from skimage._shared import testing
from skimage._shared.testing import (assert_array_equal,
                                     assert_array_almost_equal,
                                     assert_almost_equal)


# Test integer histograms
# =======================

def test_negative_overflow():
    im = np.array([-1, 127], dtype=np.int8)
    frequencies, bin_centers = exposure.histogram(im)
    assert_array_equal(bin_centers, np.arange(-1, 128))
    assert frequencies[0] == 1
    assert frequencies[-1] == 1
    assert_array_equal(frequencies[1:-1], 0)


def test_all_negative_image():
    im = np.array([-128, -1], dtype=np.int8)
    frequencies, bin_centers = exposure.histogram(im)
    assert_array_equal(bin_centers, np.arange(-128, 0))
    assert frequencies[0] == 1
    assert frequencies[-1] == 1
    assert_array_equal(frequencies[1:-1], 0)


# Test histogram equalization
# ===========================

np.random.seed(0)

test_img_int = data.camera()
# squeeze image intensities to lower image contrast
test_img = skimage.img_as_float(test_img_int)
test_img = exposure.rescale_intensity(test_img / 5. + 100)


def test_equalize_uint8_approx():
    """Check integer bins used for uint8 images."""
    img_eq0 = exposure.equalize_hist(test_img_int)
    img_eq1 = exposure.equalize_hist(test_img_int, nbins=3)
    np.testing.assert_allclose(img_eq0, img_eq1)


def test_equalize_ubyte():
    with expected_warnings(['precision loss']):
        img = skimage.img_as_ubyte(test_img)
    img_eq = exposure.equalize_hist(img)

    cdf, bin_edges = exposure.cumulative_distribution(img_eq)
    check_cdf_slope(cdf)


def test_equalize_float():
    img = skimage.img_as_float(test_img)
    img_eq = exposure.equalize_hist(img)

    cdf, bin_edges = exposure.cumulative_distribution(img_eq)
    check_cdf_slope(cdf)


def test_equalize_masked():
    img = skimage.img_as_float(test_img)
    mask = np.zeros(test_img.shape)
    mask[50:150, 50:250] = 1
    img_mask_eq = exposure.equalize_hist(img, mask=mask)
    img_eq = exposure.equalize_hist(img)

    cdf, bin_edges = exposure.cumulative_distribution(img_mask_eq)
    check_cdf_slope(cdf)

    assert not (img_eq == img_mask_eq).all()


def check_cdf_slope(cdf):
    """Slope of cdf which should equal 1 for an equalized histogram."""
    norm_intensity = np.linspace(0, 1, len(cdf))
    slope, intercept = np.polyfit(norm_intensity, cdf, 1)
    assert 0.9 < slope < 1.1


# Test intensity range
# ====================


@testing.parametrize("test_input,expected", [
    ('image', [0, 1]),
    ('dtype', [0, 255]),
    ((10, 20), [10, 20])
])
def test_intensity_range_uint8(test_input, expected):
    image = np.array([0, 1], dtype=np.uint8)
    out = intensity_range(image, range_values=test_input)
    assert_array_equal(out, expected)


@testing.parametrize("test_input,expected", [
    ('image', [0.1, 0.2]),
    ('dtype', [-1, 1]),
    ((0.3, 0.4), [0.3, 0.4])
])
def test_intensity_range_float(test_input, expected):
    image = np.array([0.1, 0.2], dtype=np.float64)
    out = intensity_range(image, range_values=test_input)
    assert_array_equal(out, expected)


def test_intensity_range_clipped_float():
    image = np.array([0.1, 0.2], dtype=np.float64)
    out = intensity_range(image, range_values='dtype', clip_negative=True)
    assert_array_equal(out, (0, 1))


# Test rescale intensity
# ======================

uint10_max = 2**10 - 1
uint12_max = 2**12 - 1
uint14_max = 2**14 - 1
uint16_max = 2**16 - 1


def test_rescale_stretch():
    image = np.array([51, 102, 153], dtype=np.uint8)
    out = exposure.rescale_intensity(image)
    assert out.dtype == np.uint8
    assert_array_almost_equal(out, [0, 127, 255])


def test_rescale_shrink():
    image = np.array([51., 102., 153.])
    out = exposure.rescale_intensity(image)
    assert_array_almost_equal(out, [0, 0.5, 1])


def test_rescale_in_range():
    image = np.array([51., 102., 153.])
    out = exposure.rescale_intensity(image, in_range=(0, 255))
    assert_array_almost_equal(out, [0.2, 0.4, 0.6])


def test_rescale_in_range_clip():
    image = np.array([51., 102., 153.])
    out = exposure.rescale_intensity(image, in_range=(0, 102))
    assert_array_almost_equal(out, [0.5, 1, 1])


def test_rescale_out_range():
    image = np.array([-10, 0, 10], dtype=np.int8)
    out = exposure.rescale_intensity(image, out_range=(0, 127))
    assert out.dtype == np.int8
    assert_array_almost_equal(out, [0, 63, 127])


def test_rescale_named_in_range():
    image = np.array([0, uint10_max, uint10_max + 100], dtype=np.uint16)
    out = exposure.rescale_intensity(image, in_range='uint10')
    assert_array_almost_equal(out, [0, uint16_max, uint16_max])


def test_rescale_named_out_range():
    image = np.array([0, uint16_max], dtype=np.uint16)
    out = exposure.rescale_intensity(image, out_range='uint10')
    assert_array_almost_equal(out, [0, uint10_max])


def test_rescale_uint12_limits():
    image = np.array([0, uint16_max], dtype=np.uint16)
    out = exposure.rescale_intensity(image, out_range='uint12')
    assert_array_almost_equal(out, [0, uint12_max])


def test_rescale_uint14_limits():
    image = np.array([0, uint16_max], dtype=np.uint16)
    out = exposure.rescale_intensity(image, out_range='uint14')
    assert_array_almost_equal(out, [0, uint14_max])


# Test adaptive histogram equalization
# ====================================

def test_adapthist_grayscale():
    """Test a grayscale float image
    """
    img = skimage.img_as_float(data.astronaut())
    img = rgb2gray(img)
    img = np.dstack((img, img, img))
    with expected_warnings(['precision loss|non-contiguous input']):
        adapted = exposure.equalize_adapthist(img, kernel_size=(57, 51),
                                              clip_limit=0.01, nbins=128)
    assert img.shape == adapted.shape
    assert_almost_equal(peak_snr(img, adapted), 102.078, 3)
    assert_almost_equal(norm_brightness_err(img, adapted), 0.0529, 3)


def test_adapthist_color():
    """Test an RGB color uint16 image
    """
    img = skimage.img_as_uint(data.astronaut())
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        hist, bin_centers = exposure.histogram(img)
        assert len(w) > 0
    with expected_warnings(['precision loss']):
        adapted = exposure.equalize_adapthist(img, clip_limit=0.01)

    assert adapted.min() == 0
    assert adapted.max() == 1.0
    assert img.shape == adapted.shape
    full_scale = skimage.exposure.rescale_intensity(img)
    assert_almost_equal(peak_snr(full_scale, adapted), 109.393, 1)
    assert_almost_equal(norm_brightness_err(full_scale, adapted), 0.02, 2)
    return data, adapted


def test_adapthist_alpha():
    """Test an RGBA color image
    """
    img = skimage.img_as_float(data.astronaut())
    alpha = np.ones((img.shape[0], img.shape[1]), dtype=float)
    img = np.dstack((img, alpha))
    with expected_warnings(['precision loss']):
        adapted = exposure.equalize_adapthist(img)
    assert adapted.shape != img.shape
    img = img[:, :, :3]
    full_scale = skimage.exposure.rescale_intensity(img)
    assert img.shape == adapted.shape
    assert_almost_equal(peak_snr(full_scale, adapted), 109.393, 2)
    assert_almost_equal(norm_brightness_err(full_scale, adapted), 0.0248, 3)


def peak_snr(img1, img2):
    """Peak signal to noise ratio of two images

    Parameters
    ----------
    img1 : array-like
    img2 : array-like

    Returns
    -------
    peak_snr : float
        Peak signal to noise ratio
    """
    if img1.ndim == 3:
        img1, img2 = rgb2gray(img1.copy()), rgb2gray(img2.copy())
    img1 = skimage.img_as_float(img1)
    img2 = skimage.img_as_float(img2)
    mse = 1. / img1.size * np.square(img1 - img2).sum()
    _, max_ = dtype_range[img1.dtype.type]
    return 20 * np.log(max_ / mse)


def norm_brightness_err(img1, img2):
    """Normalized Absolute Mean Brightness Error between two images

    Parameters
    ----------
    img1 : array-like
    img2 : array-like

    Returns
    -------
    norm_brightness_error : float
        Normalized absolute mean brightness error
    """
    if img1.ndim == 3:
        img1, img2 = rgb2gray(img1), rgb2gray(img2)
    ambe = np.abs(img1.mean() - img2.mean())
    nbe = ambe / dtype_range[img1.dtype.type][1]
    return nbe


# Test Gamma Correction
# =====================


def test_adjust_gamma_one():
    """Same image should be returned for gamma equal to one"""
    image = np.random.uniform(0, 255, (8, 8))
    result = exposure.adjust_gamma(image, 1)
    assert_array_equal(result, image)


def test_adjust_gamma_zero():
    """White image should be returned for gamma equal to zero"""
    image = np.random.uniform(0, 255, (8, 8))
    result = exposure.adjust_gamma(image, 0)
    dtype = image.dtype.type
    assert_array_equal(result, dtype_range[dtype][1])


def test_adjust_gamma_less_one():
    """Verifying the output with expected results for gamma
    correction with gamma equal to half"""
    image = np.arange(0, 255, 4, np.uint8).reshape((8, 8))
    expected = np.array([
        [  0,  31,  45,  55,  63,  71,  78,  84],
        [ 90,  95, 100, 105, 110, 115, 119, 123],
        [127, 131, 135, 139, 142, 146, 149, 153],
        [156, 159, 162, 165, 168, 171, 174, 177],
        [180, 183, 186, 188, 191, 194, 196, 199],
        [201, 204, 206, 209, 211, 214, 216, 218],
        [221, 223, 225, 228, 230, 232, 234, 236],
        [238, 241, 243, 245, 247, 249, 251, 253]], dtype=np.uint8)

    result = exposure.adjust_gamma(image, 0.5)
    assert_array_equal(result, expected)


def test_adjust_gamma_greater_one():
    """Verifying the output with expected results for gamma
    correction with gamma equal to two"""
    image = np.arange(0, 255, 4, np.uint8).reshape((8, 8))
    expected = np.array([
        [  0,   0,   0,   0,   1,   1,   2,   3],
        [  4,   5,   6,   7,   9,  10,  12,  14],
        [ 16,  18,  20,  22,  25,  27,  30,  33],
        [ 36,  39,  42,  45,  49,  52,  56,  60],
        [ 64,  68,  72,  76,  81,  85,  90,  95],
        [100, 105, 110, 116, 121, 127, 132, 138],
        [144, 150, 156, 163, 169, 176, 182, 189],
        [196, 203, 211, 218, 225, 233, 241, 249]], dtype=np.uint8)

    result = exposure.adjust_gamma(image, 2)
    assert_array_equal(result, expected)


def test_adjust_gamma_neggative():
    image = np.arange(0, 255, 4, np.uint8).reshape((8, 8))
    with testing.raises(ValueError):
        exposure.adjust_gamma(image, -1)


# Test Logarithmic Correction
# ===========================

def test_adjust_log():
    """Verifying the output with expected results for logarithmic
    correction with multiplier constant multiplier equal to unity"""
    image = np.arange(0, 255, 4, np.uint8).reshape((8, 8))
    expected = np.array([
        [  0,   5,  11,  16,  22,  27,  33,  38],
        [ 43,  48,  53,  58,  63,  68,  73,  77],
        [ 82,  86,  91,  95, 100, 104, 109, 113],
        [117, 121, 125, 129, 133, 137, 141, 145],
        [149, 153, 157, 160, 164, 168, 172, 175],
        [179, 182, 186, 189, 193, 196, 199, 203],
        [206, 209, 213, 216, 219, 222, 225, 228],
        [231, 234, 238, 241, 244, 246, 249, 252]], dtype=np.uint8)

    result = exposure.adjust_log(image, 1)
    assert_array_equal(result, expected)


def test_adjust_inv_log():
    """Verifying the output with expected results for inverse logarithmic
    correction with multiplier constant multiplier equal to unity"""
    image = np.arange(0, 255, 4, np.uint8).reshape((8, 8))
    expected = np.array([
        [  0,   2,   5,   8,  11,  14,  17,  20],
        [ 23,  26,  29,  32,  35,  38,  41,  45],
        [ 48,  51,  55,  58,  61,  65,  68,  72],
        [ 76,  79,  83,  87,  90,  94,  98, 102],
        [106, 110, 114, 118, 122, 126, 130, 134],
        [138, 143, 147, 151, 156, 160, 165, 170],
        [174, 179, 184, 188, 193, 198, 203, 208],
        [213, 218, 224, 229, 234, 239, 245, 250]], dtype=np.uint8)

    result = exposure.adjust_log(image, 1, True)
    assert_array_equal(result, expected)


# Test Sigmoid Correction
# =======================

def test_adjust_sigmoid_cutoff_one():
    """Verifying the output with expected results for sigmoid correction
    with cutoff equal to one and gain of 5"""
    image = np.arange(0, 255, 4, np.uint8).reshape((8, 8))
    expected = np.array([
        [  1,   1,   1,   2,   2,   2,   2,   2],
        [  3,   3,   3,   4,   4,   4,   5,   5],
        [  5,   6,   6,   7,   7,   8,   9,  10],
        [ 10,  11,  12,  13,  14,  15,  16,  18],
        [ 19,  20,  22,  24,  25,  27,  29,  32],
        [ 34,  36,  39,  41,  44,  47,  50,  54],
        [ 57,  61,  64,  68,  72,  76,  80,  85],
        [ 89,  94,  99, 104, 108, 113, 118, 123]], dtype=np.uint8)

    result = exposure.adjust_sigmoid(image, 1, 5)
    assert_array_equal(result, expected)


def test_adjust_sigmoid_cutoff_zero():
    """Verifying the output with expected results for sigmoid correction
    with cutoff equal to zero and gain of 10"""
    image = np.arange(0, 255, 4, np.uint8).reshape((8, 8))
    expected = np.array([
        [127, 137, 147, 156, 166, 175, 183, 191],
        [198, 205, 211, 216, 221, 225, 229, 232],
        [235, 238, 240, 242, 244, 245, 247, 248],
        [249, 250, 250, 251, 251, 252, 252, 253],
        [253, 253, 253, 253, 254, 254, 254, 254],
        [254, 254, 254, 254, 254, 254, 254, 254],
        [254, 254, 254, 254, 254, 254, 254, 254],
        [254, 254, 254, 254, 254, 254, 254, 254]], dtype=np.uint8)

    result = exposure.adjust_sigmoid(image, 0, 10)
    assert_array_equal(result, expected)


def test_adjust_sigmoid_cutoff_half():
    """Verifying the output with expected results for sigmoid correction
    with cutoff equal to half and gain of 10"""
    image = np.arange(0, 255, 4, np.uint8).reshape((8, 8))
    expected = np.array([
        [  1,   1,   2,   2,   3,   3,   4,   5],
        [  5,   6,   7,   9,  10,  12,  14,  16],
        [ 19,  22,  25,  29,  34,  39,  44,  50],
        [ 57,  64,  72,  80,  89,  99, 108, 118],
        [128, 138, 148, 158, 167, 176, 184, 192],
        [199, 205, 211, 217, 221, 226, 229, 233],
        [236, 238, 240, 242, 244, 246, 247, 248],
        [249, 250, 250, 251, 251, 252, 252, 253]], dtype=np.uint8)

    result = exposure.adjust_sigmoid(image, 0.5, 10)
    assert_array_equal(result, expected)


def test_adjust_inv_sigmoid_cutoff_half():
    """Verifying the output with expected results for inverse sigmoid
    correction with cutoff equal to half and gain of 10"""
    image = np.arange(0, 255, 4, np.uint8).reshape((8, 8))
    expected = np.array([
        [253, 253, 252, 252, 251, 251, 250, 249],
        [249, 248, 247, 245, 244, 242, 240, 238],
        [235, 232, 229, 225, 220, 215, 210, 204],
        [197, 190, 182, 174, 165, 155, 146, 136],
        [126, 116, 106,  96,  87,  78,  70,  62],
        [ 55,  49,  43,  37,  33,  28,  25,  21],
        [ 18,  16,  14,  12,  10,   8,   7,   6],
        [  5,   4,   4,   3,   3,   2,   2,   1]], dtype=np.uint8)

    result = exposure.adjust_sigmoid(image, 0.5, 10, True)
    assert_array_equal(result, expected)


def test_negative():
    image = np.arange(-10, 245, 4).reshape((8, 8)).astype(np.double)
    with testing.raises(ValueError):
        exposure.adjust_gamma(image)


def test_is_low_contrast():
    image = np.linspace(0, 0.04, 100)
    assert exposure.is_low_contrast(image)
    image[-1] = 1
    assert exposure.is_low_contrast(image)
    assert not exposure.is_low_contrast(image, upper_percentile=100)

    image = (image * 255).astype(np.uint8)
    assert exposure.is_low_contrast(image)
    assert not exposure.is_low_contrast(image, upper_percentile=100)

    image = (image.astype(np.uint16)) * 2**8
    assert exposure.is_low_contrast(image)
    assert not exposure.is_low_contrast(image, upper_percentile=100)
