import numpy as np
from skimage._shared.testing import assert_equal
from skimage import data
from skimage import transform as tf
from skimage.color import rgb2gray
from skimage.feature import (BRIEF, match_descriptors,
                             corner_peaks, corner_harris)
from skimage._shared import testing


def test_binary_descriptors_unequal_descriptor_sizes_error():
    """Sizes of descriptors of keypoints to be matched should be equal."""
    descs1 = np.array([[True, True, False, True],
                       [False, True, False, True]])
    descs2 = np.array([[True, False, False, True, False],
                       [False, True, True, True, False]])
    with testing.raises(ValueError):
        match_descriptors(descs1, descs2)


def test_binary_descriptors():
    descs1 = np.array([[True, True, False, True, True],
                       [False, True, False, True, True]])
    descs2 = np.array([[True, False, False, True, False],
                       [False, False, True, True, True]])
    matches = match_descriptors(descs1, descs2)
    assert_equal(matches, [[0, 0], [1, 1]])


def test_binary_descriptors_rotation_crosscheck_false():
    """Verify matched keypoints and their corresponding masks results between
    image and its rotated version with the expected keypoint pairs with
    cross_check disabled."""
    img = data.astronaut()
    img = rgb2gray(img)
    tform = tf.SimilarityTransform(scale=1, rotation=0.15, translation=(0, 0))
    rotated_img = tf.warp(img, tform, clip=False)

    extractor = BRIEF(descriptor_size=512)

    keypoints1 = corner_peaks(corner_harris(img), min_distance=5,
                              threshold_abs=0, threshold_rel=0.1)
    extractor.extract(img, keypoints1)
    descriptors1 = extractor.descriptors

    keypoints2 = corner_peaks(corner_harris(rotated_img), min_distance=5,
                              threshold_abs=0, threshold_rel=0.1)
    extractor.extract(rotated_img, keypoints2)
    descriptors2 = extractor.descriptors

    matches = match_descriptors(descriptors1, descriptors2, cross_check=False)

    exp_matches1 = np.array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11,
                             12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                             24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35,
                             36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46])
    exp_matches2 = np.array([ 0, 31,  2,  3,  1,  4,  6,  4, 38,  5, 27,  7,
                             13, 10,  9, 27,  7, 11, 15,  8, 23, 14, 12, 16,
                             10, 25, 18, 19, 21, 20, 41, 24, 25, 26, 28, 27,
                             22, 23, 29, 30, 31, 32, 35, 33, 34, 30, 36])
    assert_equal(matches[:, 0], exp_matches1)
    assert_equal(matches[:, 1], exp_matches2)


def test_binary_descriptors_rotation_crosscheck_true():
    """Verify matched keypoints and their corresponding masks results between
    image and its rotated version with the expected keypoint pairs with
    cross_check enabled."""
    img = data.astronaut()
    img = rgb2gray(img)
    tform = tf.SimilarityTransform(scale=1, rotation=0.15, translation=(0, 0))
    rotated_img = tf.warp(img, tform, clip=False)

    extractor = BRIEF(descriptor_size=512)

    keypoints1 = corner_peaks(corner_harris(img), min_distance=5,
                              threshold_abs=0, threshold_rel=0.1)
    extractor.extract(img, keypoints1)
    descriptors1 = extractor.descriptors

    keypoints2 = corner_peaks(corner_harris(rotated_img), min_distance=5,
                              threshold_abs=0, threshold_rel=0.1)
    extractor.extract(rotated_img, keypoints2)
    descriptors2 = extractor.descriptors

    matches = match_descriptors(descriptors1, descriptors2, cross_check=True)

    exp_matches1 = np.array([ 0,  2,  3,  4,  5,  6,  9, 11, 12, 13, 14, 17,
                             18, 19, 21, 22, 23, 26, 27, 28, 29, 31, 32, 33,
                             34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 46])
    exp_matches2 = np.array([ 0,  2,  3,  1,  4,  6,  5,  7, 13, 10,  9, 11,
                             15,  8, 14, 12, 16, 18, 19, 21, 20, 24, 25, 26,
                             28, 27, 22, 23, 29, 30, 31, 32, 35, 33, 34, 36])
    assert_equal(matches[:, 0], exp_matches1)
    assert_equal(matches[:, 1], exp_matches2)


def test_max_distance():
    descs1 = np.zeros((10, 128))
    descs2 = np.zeros((15, 128))

    descs1[0, :] = 1

    matches =  match_descriptors(descs1, descs2, metric='euclidean',
                                 max_distance=0.1, cross_check=False)
    assert len(matches) == 9

    matches =  match_descriptors(descs1, descs2, metric='euclidean',
                                 max_distance=np.sqrt(128.1),
                                 cross_check=False)
    assert len(matches) == 10

    matches =  match_descriptors(descs1, descs2, metric='euclidean',
                                 max_distance=0.1,
                                 cross_check=True)
    assert_equal(matches, [[1, 0]])

    matches =  match_descriptors(descs1, descs2, metric='euclidean',
                                 max_distance=np.sqrt(128.1),
                                 cross_check=True)
    assert_equal(matches, [[1, 0]])


def test_max_ratio():
    descs1 = 10 * np.arange(10)[:, None].astype(np.float32)
    descs2 = 10 * np.arange(15)[:, None].astype(np.float32)

    descs2[0] = 5.0

    matches = match_descriptors(descs1, descs2, metric='euclidean',
                                max_ratio=1.0, cross_check=False)
    assert_equal(len(matches), 10)

    matches = match_descriptors(descs1, descs2, metric='euclidean',
                                max_ratio=0.6, cross_check=False)
    assert_equal(len(matches), 10)

    matches = match_descriptors(descs1, descs2, metric='euclidean',
                                max_ratio=0.5, cross_check=False)
    assert_equal(len(matches), 9)

    descs1[0] = 7.5

    matches = match_descriptors(descs1, descs2, metric='euclidean',
                                max_ratio=0.5, cross_check=False)
    assert_equal(len(matches), 9)

    descs2 = 10 * np.arange(1)[:, None].astype(np.float32)

    matches = match_descriptors(descs1, descs2, metric='euclidean',
                                max_ratio=1.0, cross_check=False)
    assert_equal(len(matches), 10)

    matches = match_descriptors(descs1, descs2, metric='euclidean',
                                max_ratio=0.5, cross_check=False)
    assert_equal(len(matches), 10)

    descs1 = 10 * np.arange(1)[:, None].astype(np.float32)

    matches = match_descriptors(descs1, descs2, metric='euclidean',
                                max_ratio=1.0, cross_check=False)
    assert_equal(len(matches), 1)

    matches = match_descriptors(descs1, descs2, metric='euclidean',
                                max_ratio=0.5, cross_check=False)
    assert_equal(len(matches), 1)
