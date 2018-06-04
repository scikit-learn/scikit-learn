import numpy as np
from skimage import filters
from skimage.filters.edges import _mask_filter_result

from skimage._shared import testing
from skimage._shared.testing import (assert_array_almost_equal,
                                     assert_, assert_allclose)


def test_roberts_zeros():
    """Roberts' filter on an array of all zeros."""
    result = filters.roberts(np.zeros((10, 10)), np.ones((10, 10), bool))
    assert (np.all(result == 0))


def test_roberts_diagonal1():
    """Roberts' filter on a diagonal edge should be a diagonal line."""
    image = np.tri(10, 10, 0)
    expected = ~(np.tri(10, 10, -1).astype(bool) |
                 np.tri(10, 10, -2).astype(bool).transpose())
    expected = _mask_filter_result(expected, None)
    result = filters.roberts(image).astype(bool)
    assert_array_almost_equal(result, expected)


def test_roberts_diagonal2():
    """Roberts' filter on a diagonal edge should be a diagonal line."""
    image = np.rot90(np.tri(10, 10, 0), 3)
    expected = ~np.rot90(np.tri(10, 10, -1).astype(bool) |
                         np.tri(10, 10, -2).astype(bool).transpose())
    expected = _mask_filter_result(expected, None)
    result = filters.roberts(image).astype(bool)
    assert_array_almost_equal(result, expected)


def test_sobel_zeros():
    """Sobel on an array of all zeros."""
    result = filters.sobel(np.zeros((10, 10)), np.ones((10, 10), bool))
    assert (np.all(result == 0))


def test_sobel_mask():
    """Sobel on a masked array should be zero."""
    np.random.seed(0)
    result = filters.sobel(np.random.uniform(size=(10, 10)),
                           np.zeros((10, 10), bool))
    assert (np.all(result == 0))


def test_sobel_horizontal():
    """Sobel on a horizontal edge should be a horizontal line."""
    i, j = np.mgrid[-5:6, -5:6]
    image = (i >= 0).astype(float)
    result = filters.sobel(image) * np.sqrt(2)
    # Fudge the eroded points
    i[np.abs(j) == 5] = 10000
    assert_allclose(result[i == 0], 1)
    assert (np.all(result[np.abs(i) > 1] == 0))


def test_sobel_vertical():
    """Sobel on a vertical edge should be a vertical line."""
    i, j = np.mgrid[-5:6, -5:6]
    image = (j >= 0).astype(float)
    result = filters.sobel(image) * np.sqrt(2)
    j[np.abs(i) == 5] = 10000
    assert (np.all(result[j == 0] == 1))
    assert (np.all(result[np.abs(j) > 1] == 0))


def test_sobel_h_zeros():
    """Horizontal sobel on an array of all zeros."""
    result = filters.sobel_h(np.zeros((10, 10)), np.ones((10, 10), bool))
    assert (np.all(result == 0))


def test_sobel_h_mask():
    """Horizontal Sobel on a masked array should be zero."""
    np.random.seed(0)
    result = filters.sobel_h(np.random.uniform(size=(10, 10)),
                             np.zeros((10, 10), bool))
    assert (np.all(result == 0))


def test_sobel_h_horizontal():
    """Horizontal Sobel on an edge should be a horizontal line."""
    i, j = np.mgrid[-5:6, -5:6]
    image = (i >= 0).astype(float)
    result = filters.sobel_h(image)
    # Fudge the eroded points
    i[np.abs(j) == 5] = 10000
    assert (np.all(result[i == 0] == 1))
    assert (np.all(result[np.abs(i) > 1] == 0))


def test_sobel_h_vertical():
    """Horizontal Sobel on a vertical edge should be zero."""
    i, j = np.mgrid[-5:6, -5:6]
    image = (j >= 0).astype(float) * np.sqrt(2)
    result = filters.sobel_h(image)
    assert_allclose(result, 0, atol=1e-10)


def test_sobel_v_zeros():
    """Vertical sobel on an array of all zeros."""
    result = filters.sobel_v(np.zeros((10, 10)), np.ones((10, 10), bool))
    assert_allclose(result, 0)


def test_sobel_v_mask():
    """Vertical Sobel on a masked array should be zero."""
    np.random.seed(0)
    result = filters.sobel_v(np.random.uniform(size=(10, 10)),
                             np.zeros((10, 10), bool))
    assert_allclose(result, 0)


def test_sobel_v_vertical():
    """Vertical Sobel on an edge should be a vertical line."""
    i, j = np.mgrid[-5:6, -5:6]
    image = (j >= 0).astype(float)
    result = filters.sobel_v(image)
    # Fudge the eroded points
    j[np.abs(i) == 5] = 10000
    assert (np.all(result[j == 0] == 1))
    assert (np.all(result[np.abs(j) > 1] == 0))


def test_sobel_v_horizontal():
    """vertical Sobel on a horizontal edge should be zero."""
    i, j = np.mgrid[-5:6, -5:6]
    image = (i >= 0).astype(float)
    result = filters.sobel_v(image)
    assert_allclose(result, 0)


def test_scharr_zeros():
    """Scharr on an array of all zeros."""
    result = filters.scharr(np.zeros((10, 10)), np.ones((10, 10), bool))
    assert (np.all(result < 1e-16))


def test_scharr_mask():
    """Scharr on a masked array should be zero."""
    np.random.seed(0)
    result = filters.scharr(np.random.uniform(size=(10, 10)),
                            np.zeros((10, 10), bool))
    assert_allclose(result, 0)


def test_scharr_horizontal():
    """Scharr on an edge should be a horizontal line."""
    i, j = np.mgrid[-5:6, -5:6]
    image = (i >= 0).astype(float)
    result = filters.scharr(image) * np.sqrt(2)
    # Fudge the eroded points
    i[np.abs(j) == 5] = 10000
    assert_allclose(result[i == 0], 1)
    assert (np.all(result[np.abs(i) > 1] == 0))


def test_scharr_vertical():
    """Scharr on a vertical edge should be a vertical line."""
    i, j = np.mgrid[-5:6, -5:6]
    image = (j >= 0).astype(float)
    result = filters.scharr(image) * np.sqrt(2)
    j[np.abs(i) == 5] = 10000
    assert_allclose(result[j == 0], 1)
    assert (np.all(result[np.abs(j) > 1] == 0))


def test_scharr_h_zeros():
    """Horizontal Scharr on an array of all zeros."""
    result = filters.scharr_h(np.zeros((10, 10)), np.ones((10, 10), bool))
    assert_allclose(result, 0)


def test_scharr_h_mask():
    """Horizontal Scharr on a masked array should be zero."""
    np.random.seed(0)
    result = filters.scharr_h(np.random.uniform(size=(10, 10)),
                              np.zeros((10, 10), bool))
    assert_allclose(result, 0)


def test_scharr_h_horizontal():
    """Horizontal Scharr on an edge should be a horizontal line."""
    i, j = np.mgrid[-5:6, -5:6]
    image = (i >= 0).astype(float)
    result = filters.scharr_h(image)
    # Fudge the eroded points
    i[np.abs(j) == 5] = 10000
    assert (np.all(result[i == 0] == 1))
    assert (np.all(result[np.abs(i) > 1] == 0))


def test_scharr_h_vertical():
    """Horizontal Scharr on a vertical edge should be zero."""
    i, j = np.mgrid[-5:6, -5:6]
    image = (j >= 0).astype(float)
    result = filters.scharr_h(image)
    assert_allclose(result, 0)


def test_scharr_v_zeros():
    """Vertical Scharr on an array of all zeros."""
    result = filters.scharr_v(np.zeros((10, 10)), np.ones((10, 10), bool))
    assert_allclose(result, 0)


def test_scharr_v_mask():
    """Vertical Scharr on a masked array should be zero."""
    np.random.seed(0)
    result = filters.scharr_v(np.random.uniform(size=(10, 10)),
                              np.zeros((10, 10), bool))
    assert_allclose(result, 0)


def test_scharr_v_vertical():
    """Vertical Scharr on an edge should be a vertical line."""
    i, j = np.mgrid[-5:6, -5:6]
    image = (j >= 0).astype(float)
    result = filters.scharr_v(image)
    # Fudge the eroded points
    j[np.abs(i) == 5] = 10000
    assert (np.all(result[j == 0] == 1))
    assert (np.all(result[np.abs(j) > 1] == 0))


def test_scharr_v_horizontal():
    """vertical Scharr on a horizontal edge should be zero."""
    i, j = np.mgrid[-5:6, -5:6]
    image = (i >= 0).astype(float)
    result = filters.scharr_v(image)
    assert_allclose(result, 0)


def test_prewitt_zeros():
    """Prewitt on an array of all zeros."""
    result = filters.prewitt(np.zeros((10, 10)), np.ones((10, 10), bool))
    assert_allclose(result, 0)


def test_prewitt_mask():
    """Prewitt on a masked array should be zero."""
    np.random.seed(0)
    result = filters.prewitt(np.random.uniform(size=(10, 10)),
                             np.zeros((10, 10), bool))
    assert_allclose(np.abs(result), 0)


def test_prewitt_horizontal():
    """Prewitt on an edge should be a horizontal line."""
    i, j = np.mgrid[-5:6, -5:6]
    image = (i >= 0).astype(float)
    result = filters.prewitt(image) * np.sqrt(2)
    # Fudge the eroded points
    i[np.abs(j) == 5] = 10000
    assert (np.all(result[i == 0] == 1))
    assert_allclose(result[np.abs(i) > 1], 0, atol=1e-10)


def test_prewitt_vertical():
    """Prewitt on a vertical edge should be a vertical line."""
    i, j = np.mgrid[-5:6, -5:6]
    image = (j >= 0).astype(float)
    result = filters.prewitt(image) * np.sqrt(2)
    j[np.abs(i) == 5] = 10000
    assert_allclose(result[j == 0], 1)
    assert_allclose(result[np.abs(j) > 1], 0, atol=1e-10)


def test_prewitt_h_zeros():
    """Horizontal prewitt on an array of all zeros."""
    result = filters.prewitt_h(np.zeros((10, 10)), np.ones((10, 10), bool))
    assert_allclose(result, 0)


def test_prewitt_h_mask():
    """Horizontal prewitt on a masked array should be zero."""
    np.random.seed(0)
    result = filters.prewitt_h(np.random.uniform(size=(10, 10)),
                               np.zeros((10, 10), bool))
    assert_allclose(result, 0)


def test_prewitt_h_horizontal():
    """Horizontal prewitt on an edge should be a horizontal line."""
    i, j = np.mgrid[-5:6, -5:6]
    image = (i >= 0).astype(float)
    result = filters.prewitt_h(image)
    # Fudge the eroded points
    i[np.abs(j) == 5] = 10000
    assert (np.all(result[i == 0] == 1))
    assert_allclose(result[np.abs(i) > 1], 0, atol=1e-10)


def test_prewitt_h_vertical():
    """Horizontal prewitt on a vertical edge should be zero."""
    i, j = np.mgrid[-5:6, -5:6]
    image = (j >= 0).astype(float)
    result = filters.prewitt_h(image)
    assert_allclose(result, 0, atol=1e-10)


def test_prewitt_v_zeros():
    """Vertical prewitt on an array of all zeros."""
    result = filters.prewitt_v(np.zeros((10, 10)), np.ones((10, 10), bool))
    assert_allclose(result, 0)


def test_prewitt_v_mask():
    """Vertical prewitt on a masked array should be zero."""
    np.random.seed(0)
    result = filters.prewitt_v(np.random.uniform(size=(10, 10)),
                               np.zeros((10, 10), bool))
    assert_allclose(result, 0)


def test_prewitt_v_vertical():
    """Vertical prewitt on an edge should be a vertical line."""
    i, j = np.mgrid[-5:6, -5:6]
    image = (j >= 0).astype(float)
    result = filters.prewitt_v(image)
    # Fudge the eroded points
    j[np.abs(i) == 5] = 10000
    assert (np.all(result[j == 0] == 1))
    assert_allclose(result[np.abs(j) > 1], 0, atol=1e-10)


def test_prewitt_v_horizontal():
    """Vertical prewitt on a horizontal edge should be zero."""
    i, j = np.mgrid[-5:6, -5:6]
    image = (i >= 0).astype(float)
    result = filters.prewitt_v(image)
    assert_allclose(result, 0)


def test_laplace_zeros():
    """Laplace on a square image."""
    # Create a synthetic 2D image
    image = np.zeros((9, 9))
    image[3:-3, 3:-3] = 1
    result = filters.laplace(image)
    res_chk = np.array([[ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                        [ 0.,  0.,  0., -1., -1., -1.,  0.,  0.,  0.],
                        [ 0.,  0., -1.,  2.,  1.,  2., -1.,  0.,  0.],
                        [ 0.,  0., -1.,  1.,  0.,  1., -1.,  0.,  0.],
                        [ 0.,  0., -1.,  2.,  1.,  2., -1.,  0.,  0.],
                        [ 0.,  0.,  0., -1., -1., -1.,  0.,  0.,  0.],
                        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.]])
    assert_allclose(result, res_chk)


def test_laplace_mask():
    """Laplace on a masked array should be zero."""
    # Create a synthetic 2D image
    image = np.zeros((9, 9))
    image[3:-3, 3:-3] = 1
    # Define the mask
    result = filters.laplace(image, ksize=3, mask=np.zeros((9, 9), bool))
    assert (np.all(result == 0))


@testing.parametrize("grad_func", (
    filters.prewitt_h, filters.sobel_h, filters.scharr_h
))
def test_horizontal_mask_line(grad_func):
    """Horizontal edge filters mask pixels surrounding input mask."""
    vgrad, _ = np.mgrid[:1:11j, :1:11j]  # vertical gradient with spacing 0.1
    vgrad[5, :] = 1                      # bad horizontal line

    mask = np.ones_like(vgrad)
    mask[5, :] = 0                       # mask bad line

    expected = np.zeros_like(vgrad)
    expected[1:-1, 1:-1] = 0.2           # constant gradient for most of image,
    expected[4:7, 1:-1] = 0              # but line and neighbors masked

    result = grad_func(vgrad, mask)
    assert_allclose(result, expected)


@testing.parametrize("grad_func", (
    filters.prewitt_v, filters.sobel_v, filters.scharr_v
))
def test_vertical_mask_line(grad_func):
    """Vertical edge filters mask pixels surrounding input mask."""
    _, hgrad = np.mgrid[:1:11j, :1:11j]  # horizontal gradient with spacing 0.1
    hgrad[:, 5] = 1                      # bad vertical line

    mask = np.ones_like(hgrad)
    mask[:, 5] = 0                       # mask bad line

    expected = np.zeros_like(hgrad)
    expected[1:-1, 1:-1] = 0.2           # constant gradient for most of image,
    expected[1:-1, 4:7] = 0              # but line and neighbors masked

    result = grad_func(hgrad, mask)
    assert_allclose(result, expected)


def test_range():
    """Output of edge detection should be in [0, 1]"""
    image = np.random.random((100, 100))
    for detector in (filters.sobel, filters.scharr,
                     filters.prewitt, filters.roberts):
        out = detector(image)
        assert_(out.min() >= 0,
                "Minimum of `{0}` is smaller than zero".format(
                    detector.__name__)
                )
        assert_(out.max() <= 1,
                "Maximum of `{0}` is larger than 1".format(
                    detector.__name__)
                )
