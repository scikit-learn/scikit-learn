import numpy as np
import pytest
from numpy.testing import (
    assert_almost_equal,
    assert_array_almost_equal,
    assert_array_equal,
    assert_array_less,
    assert_equal,
)
from scipy.ndimage import fourier_shift, shift as real_shift
import scipy.fft as fft

from skimage._shared.testing import fetch
from skimage._shared.utils import _supported_float_type
from skimage.data import camera, brain


from skimage.io import imread
from skimage.registration._masked_phase_cross_correlation import (
    _masked_phase_cross_correlation as masked_register_translation,
    cross_correlate_masked,
)
from skimage.registration import phase_cross_correlation


def test_masked_registration_vs_phase_cross_correlation():
    """masked_register_translation should give the same results as
    phase_cross_correlation in the case of trivial masks."""
    reference_image = camera()
    shift = (-7, 12)
    shifted = np.real(fft.ifft2(fourier_shift(fft.fft2(reference_image), shift)))
    trivial_mask = np.ones_like(reference_image)

    nonmasked_result, *_ = phase_cross_correlation(reference_image, shifted)
    masked_result = masked_register_translation(
        reference_image, shifted, reference_mask=trivial_mask, overlap_ratio=1 / 10
    )

    assert_equal(nonmasked_result, masked_result)


def test_masked_registration_random_masks():
    """masked_register_translation should be able to register translations
    between images even with random masks."""
    # See random number generator for reproducible results
    np.random.seed(23)

    reference_image = camera()
    shift = (-7, 12)
    shifted = np.real(fft.ifft2(fourier_shift(fft.fft2(reference_image), shift)))

    # Random masks with 75% of pixels being valid
    ref_mask = np.random.choice([True, False], reference_image.shape, p=[3 / 4, 1 / 4])
    shifted_mask = np.random.choice([True, False], shifted.shape, p=[3 / 4, 1 / 4])

    measured_shift = masked_register_translation(
        reference_image, shifted, reference_mask=ref_mask, moving_mask=shifted_mask
    )
    assert_equal(measured_shift, -np.array(shift))


def test_masked_registration_3d_contiguous_mask():
    """masked_register_translation should be able to register translations
    between volumes with contiguous masks."""
    ref_vol = brain()[:, ::2, ::2]

    offset = (1, -5, 10)

    # create square mask
    ref_mask = np.zeros_like(ref_vol, dtype=bool)
    ref_mask[:-2, 75:100, 75:100] = True
    ref_shifted = real_shift(ref_vol, offset)

    measured_offset = masked_register_translation(
        ref_vol, ref_shifted, reference_mask=ref_mask, moving_mask=ref_mask
    )

    assert_equal(offset, -np.array(measured_offset))


def test_masked_registration_random_masks_non_equal_sizes():
    """masked_register_translation should be able to register
    translations between images that are not the same size even
    with random masks."""
    # See random number generator for reproducible results
    np.random.seed(23)

    reference_image = camera()
    shift = (-7, 12)
    shifted = np.real(fft.ifft2(fourier_shift(fft.fft2(reference_image), shift)))

    # Crop the shifted image
    shifted = shifted[64:-64, 64:-64]

    # Random masks with 75% of pixels being valid
    ref_mask = np.random.choice([True, False], reference_image.shape, p=[3 / 4, 1 / 4])
    shifted_mask = np.random.choice([True, False], shifted.shape, p=[3 / 4, 1 / 4])

    measured_shift = masked_register_translation(
        reference_image,
        shifted,
        reference_mask=np.ones_like(ref_mask),
        moving_mask=np.ones_like(shifted_mask),
    )
    assert_equal(measured_shift, -np.array(shift))


def test_masked_registration_padfield_data():
    """Masked translation registration should behave like in the original
    publication"""
    # Test translated from MATLABimplementation `MaskedFFTRegistrationTest`
    # file. You can find the source code here:
    # http://www.dirkpadfield.com/Home/MaskedFFTRegistrationCode.zip

    shifts = [(75, 75), (-130, 130), (130, 130)]
    for xi, yi in shifts:
        fixed_image = imread(fetch(f'registration/tests/data/OriginalX{xi}Y{yi}.png'))
        moving_image = imread(
            fetch(f'registration/tests/data/TransformedX{xi}Y{yi}.png')
        )

        # Valid pixels are 1
        fixed_mask = fixed_image != 0
        moving_mask = moving_image != 0

        # Note that shifts in x and y and shifts in cols and rows
        shift_y, shift_x = masked_register_translation(
            fixed_image,
            moving_image,
            reference_mask=fixed_mask,
            moving_mask=moving_mask,
            overlap_ratio=0.1,
        )
        # Note: by looking at the test code from Padfield's
        # MaskedFFTRegistrationCode repository, the
        # shifts were not xi and yi, but xi and -yi
        assert_equal((shift_x, shift_y), (-xi, yi))


@pytest.mark.parametrize('dtype', [np.float16, np.float32, np.float64])
def test_cross_correlate_masked_output_shape(dtype):
    """Masked normalized cross-correlation should return a shape
    of N + M + 1 for each transform axis."""
    shape1 = (15, 4, 5)
    shape2 = (6, 12, 7)
    expected_full_shape = tuple(np.array(shape1) + np.array(shape2) - 1)
    expected_same_shape = shape1

    arr1 = np.zeros(shape1, dtype=dtype)
    arr2 = np.zeros(shape2, dtype=dtype)
    # Trivial masks
    m1 = np.ones_like(arr1)
    m2 = np.ones_like(arr2)

    float_dtype = _supported_float_type(dtype)

    full_xcorr = cross_correlate_masked(arr1, arr2, m1, m2, axes=(0, 1, 2), mode='full')
    assert_equal(full_xcorr.shape, expected_full_shape)
    assert full_xcorr.dtype == float_dtype

    same_xcorr = cross_correlate_masked(arr1, arr2, m1, m2, axes=(0, 1, 2), mode='same')
    assert_equal(same_xcorr.shape, expected_same_shape)
    assert same_xcorr.dtype == float_dtype


def test_cross_correlate_masked_test_against_mismatched_dimensions():
    """Masked normalized cross-correlation should raise an error if array
    dimensions along non-transformation axes are mismatched."""
    shape1 = (23, 1, 1)
    shape2 = (6, 2, 2)

    arr1 = np.zeros(shape1)
    arr2 = np.zeros(shape2)

    # Trivial masks
    m1 = np.ones_like(arr1)
    m2 = np.ones_like(arr2)

    with pytest.raises(ValueError):
        cross_correlate_masked(arr1, arr2, m1, m2, axes=(1, 2))


def test_cross_correlate_masked_output_range():
    """Masked normalized cross-correlation should return between 1 and -1."""
    # See random number generator for reproducible results
    np.random.seed(23)

    # Array dimensions must match along non-transformation axes, in
    # this case
    # axis 0
    shape1 = (15, 4, 5)
    shape2 = (15, 12, 7)

    # Initial array ranges between -5 and 5
    arr1 = 10 * np.random.random(shape1) - 5
    arr2 = 10 * np.random.random(shape2) - 5

    # random masks
    m1 = np.random.choice([True, False], arr1.shape)
    m2 = np.random.choice([True, False], arr2.shape)

    xcorr = cross_correlate_masked(arr1, arr2, m1, m2, axes=(1, 2))

    # No assert array less or equal, so we add an eps
    # Also could not find an `assert_array_greater`, Use (-xcorr) instead
    eps = np.finfo(float).eps
    assert_array_less(xcorr, 1 + eps)
    assert_array_less(-xcorr, 1 + eps)


def test_cross_correlate_masked_side_effects():
    """Masked normalized cross-correlation should not modify the inputs."""
    shape1 = (2, 2, 2)
    shape2 = (2, 2, 2)

    arr1 = np.zeros(shape1)
    arr2 = np.zeros(shape2)

    # Trivial masks
    m1 = np.ones_like(arr1)
    m2 = np.ones_like(arr2)

    for arr in (arr1, arr2, m1, m2):
        arr.setflags(write=False)

    cross_correlate_masked(arr1, arr2, m1, m2)


def test_cross_correlate_masked_over_axes():
    """Masked normalized cross-correlation over axes should be
    equivalent to a loop over non-transform axes."""
    # See random number generator for reproducible results
    np.random.seed(23)

    arr1 = np.random.random((8, 8, 5))
    arr2 = np.random.random((8, 8, 5))

    m1 = np.random.choice([True, False], arr1.shape)
    m2 = np.random.choice([True, False], arr2.shape)

    # Loop over last axis
    with_loop = np.empty_like(arr1, dtype=complex)
    for index in range(arr1.shape[-1]):
        with_loop[:, :, index] = cross_correlate_masked(
            arr1[:, :, index],
            arr2[:, :, index],
            m1[:, :, index],
            m2[:, :, index],
            axes=(0, 1),
            mode='same',
        )

    over_axes = cross_correlate_masked(arr1, arr2, m1, m2, axes=(0, 1), mode='same')

    assert_array_almost_equal(with_loop, over_axes)


def test_cross_correlate_masked_autocorrelation_trivial_masks():
    """Masked normalized cross-correlation between identical arrays
    should reduce to an autocorrelation even with random masks."""
    # See random number generator for reproducible results
    np.random.seed(23)

    arr1 = camera()

    # Random masks with 75% of pixels being valid
    m1 = np.random.choice([True, False], arr1.shape, p=[3 / 4, 1 / 4])
    m2 = np.random.choice([True, False], arr1.shape, p=[3 / 4, 1 / 4])

    xcorr = cross_correlate_masked(
        arr1, arr1, m1, m2, axes=(0, 1), mode='same', overlap_ratio=0
    ).real
    max_index = np.unravel_index(np.argmax(xcorr), xcorr.shape)

    # Autocorrelation should have maximum in center of array
    # uint8 inputs will be processed in float32, so reduce decimal to 5
    assert_almost_equal(xcorr.max(), 1, decimal=5)
    assert_array_equal(max_index, np.array(arr1.shape) / 2)
