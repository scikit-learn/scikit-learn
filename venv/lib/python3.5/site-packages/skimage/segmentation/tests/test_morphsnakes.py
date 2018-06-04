import numpy as np
from skimage.segmentation import (morphological_chan_vese,
                                  morphological_geodesic_active_contour,
                                  inverse_gaussian_gradient,
                                  circle_level_set)

from skimage._shared import testing
from skimage._shared.testing import assert_array_equal


def gaussian_blob():
    coords = np.mgrid[-5:6, -5:6]
    sqrdistances = (coords ** 2).sum(0)
    return np.exp(-sqrdistances / 10)


def test_morphsnakes_incorrect_image_shape():
    img = np.zeros((10, 10, 3))
    ls = np.zeros((10, 9))

    with testing.raises(ValueError):
        morphological_chan_vese(img, iterations=1, init_level_set=ls)
    with testing.raises(ValueError):
        morphological_geodesic_active_contour(img, iterations=1,
                                              init_level_set=ls)


def test_morphsnakes_incorrect_ndim():
    img = np.zeros((4, 4, 4, 4))
    ls = np.zeros((4, 4, 4, 4))

    with testing.raises(ValueError):
        morphological_chan_vese(img, iterations=1, init_level_set=ls)
    with testing.raises(ValueError):
        morphological_geodesic_active_contour(img, iterations=1,
                                              init_level_set=ls)


def test_morphsnakes_black():
    img = np.zeros((11, 11))
    ls = circle_level_set(img.shape, (5, 5), 3)

    ref_zeros = np.zeros(img.shape, dtype=np.int8)
    ref_ones = np.ones(img.shape, dtype=np.int8)

    acwe_ls = morphological_chan_vese(img, iterations=6, init_level_set=ls)
    assert_array_equal(acwe_ls, ref_zeros)

    gac_ls = morphological_geodesic_active_contour(img, iterations=6,
                                                   init_level_set=ls)
    assert_array_equal(gac_ls, ref_zeros)

    gac_ls2 = morphological_geodesic_active_contour(img, iterations=6,
                                                    init_level_set=ls,
                                                    balloon=1, threshold=-1,
                                                    smoothing=0)
    assert_array_equal(gac_ls2, ref_ones)

    assert acwe_ls.dtype == gac_ls.dtype == gac_ls2.dtype == np.int8


def test_morphsnakes_simple_shape_chan_vese():
    img = gaussian_blob()
    ls1 = circle_level_set(img.shape, (5, 5), 3)
    ls2 = circle_level_set(img.shape, (5, 5), 6)

    acwe_ls1 = morphological_chan_vese(img, iterations=10, init_level_set=ls1)
    acwe_ls2 = morphological_chan_vese(img, iterations=10, init_level_set=ls2)

    assert_array_equal(acwe_ls1, acwe_ls2)

    assert acwe_ls1.dtype == acwe_ls2.dtype == np.int8


def test_morphsnakes_simple_shape_geodesic_active_contour():
    img = np.float_(circle_level_set((11, 11), (5, 5), 3.5))
    gimg = inverse_gaussian_gradient(img, alpha=10.0, sigma=1.0)
    ls = circle_level_set(img.shape, (5, 5), 6)

    ref = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0],
                    [0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0],
                    [0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0],
                    [0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0],
                    [0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0],
                    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
                   dtype=np.int8)

    gac_ls = morphological_geodesic_active_contour(gimg, iterations=10,
                                                   init_level_set=ls,
                                                   balloon=-1)
    assert_array_equal(gac_ls, ref)
    assert gac_ls.dtype == np.int8


def test_init_level_sets():
    image = np.zeros((6, 6))
    checkerboard_ls = morphological_chan_vese(image, 0, 'checkerboard')
    checkerboard_ref = np.array([[0, 0, 0, 0, 0, 1],
                                 [0, 0, 0, 0, 0, 1],
                                 [0, 0, 0, 0, 0, 1],
                                 [0, 0, 0, 0, 0, 1],
                                 [0, 0, 0, 0, 0, 1],
                                 [1, 1, 1, 1, 1, 0]], dtype=np.int8)

    circle_ls = morphological_geodesic_active_contour(image, 0, 'circle')
    circle_ref = np.array([[0, 0, 0, 0, 0, 0],
                           [0, 0, 1, 1, 1, 0],
                           [0, 1, 1, 1, 1, 1],
                           [0, 1, 1, 1, 1, 1],
                           [0, 1, 1, 1, 1, 1],
                           [0, 0, 1, 1, 1, 0]], dtype=np.int8)

    assert_array_equal(checkerboard_ls, checkerboard_ref)
    assert_array_equal(circle_ls, circle_ref)


def test_morphsnakes_3d():
    image = np.zeros((7, 7, 7))

    evolution = []

    def callback(x):
        evolution.append(x.sum())

    ls = morphological_chan_vese(image, 5, 'circle',
                                 iter_callback=callback)

    # Check that the initial circle level set is correct
    assert evolution[0] == 81

    # Check that the final level set is correct
    assert ls.sum() == 0

    # Check that the contour is shrinking at every iteration
    for v1, v2 in zip(evolution[:-1], evolution[1:]):
        assert v1 >= v2
