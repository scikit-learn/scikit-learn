import numpy as np
from skimage._shared.testing import assert_array_equal
from skimage import data
from skimage.feature import BRIEF, corner_peaks, corner_harris
from skimage._shared import testing
from skimage._shared.testing import test_parallel


def test_color_image_unsupported_error():
    """Brief descriptors can be evaluated on gray-scale images only."""
    img = np.zeros((20, 20, 3))
    keypoints = np.asarray([[7, 5], [11, 13]])
    with testing.raises(ValueError):
        BRIEF().extract(img, keypoints)


def test_normal_mode():
    """Verify the computed BRIEF descriptors with expected for normal mode."""
    img = data.coins()

    keypoints = corner_peaks(corner_harris(img), min_distance=5,
                             threshold_abs=0, threshold_rel=0.1)

    extractor = BRIEF(descriptor_size=8, sigma=2)

    extractor.extract(img, keypoints[:8])

    expected = np.array([[False,  True, False, False,  True, False,  True, False],
                         [ True, False,  True,  True, False,  True, False, False],
                         [ True, False, False,  True, False,  True, False,  True],
                         [ True,  True,  True,  True, False,  True, False,  True],
                         [ True,  True,  True, False, False,  True,  True,  True],
                         [False, False, False, False,  True, False, False, False],
                         [False,  True, False, False,  True, False,  True, False],
                         [False, False, False, False, False, False, False, False]], dtype=bool)

    assert_array_equal(extractor.descriptors, expected)


def test_uniform_mode():
    """Verify the computed BRIEF descriptors with expected for uniform mode."""
    img = data.coins()

    keypoints = corner_peaks(corner_harris(img), min_distance=5,
                             threshold_abs=0, threshold_rel=0.1)

    extractor = BRIEF(descriptor_size=8, sigma=2, mode='uniform')

    extractor.extract(img, keypoints[:8])

    expected = np.array([[False, False, False,  True,  True,  True, False, False],
                         [ True,  True,  True, False,  True, False, False,  True],
                         [ True,  True,  True, False,  True,  True, False,  True],
                         [ True,  True,  True,  True, False,  True, False,  True],
                         [ True,  True,  True,  True,  True,  True, False, False],
                         [ True,  True,  True,  True,  True,  True,  True,  True],
                         [False, False, False,  True,  True,  True,  True,  True],
                         [False,  True, False,  True, False,  True,  True,  True]], dtype=bool)

    assert_array_equal(extractor.descriptors, expected)


def test_unsupported_mode():
    with testing.raises(ValueError):
        BRIEF(mode='foobar')


def test_border():
    img = np.zeros((100, 100))
    keypoints = np.array([[1, 1], [20, 20], [50, 50], [80, 80]])

    extractor = BRIEF(patch_size=41)
    extractor.extract(img, keypoints)

    assert extractor.descriptors.shape[0] == 3
    assert_array_equal(extractor.mask, (False, True, True, True))
