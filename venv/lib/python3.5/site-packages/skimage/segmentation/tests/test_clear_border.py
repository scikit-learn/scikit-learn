import numpy as np
from skimage.segmentation import clear_border

from skimage._shared.testing import assert_array_equal, assert_


def test_clear_border():
    image = np.array(
        [[0, 0, 0, 0, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 1, 0, 0, 0, 0],
         [1, 0, 0, 1, 0, 1, 0, 0, 0],
         [0, 0, 1, 1, 1, 1, 1, 0, 0],
         [0, 1, 1, 1, 1, 1, 1, 1, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0]])

    # test default case
    result = clear_border(image.copy())
    ref = image.copy()
    ref[2, 0] = 0
    ref[0, -2] = 0
    assert_array_equal(result, ref)

    # test buffer
    result = clear_border(image.copy(), 1)
    assert_array_equal(result, np.zeros(result.shape))

    # test background value
    result = clear_border(image.copy(), buffer_size=1, bgval=2)
    assert_array_equal(result, 2 * np.ones_like(image))


def test_clear_border_3d():
    image = np.array([
        [[0, 0, 0, 0],
         [0, 0, 0, 0],
         [0, 0, 0, 0],
         [1, 0, 0, 0]],
        [[0, 0, 0, 0],
         [0, 1, 1, 0],
         [0, 0, 1, 0],
         [0, 0, 0, 0]],
        [[0, 0, 0, 0],
         [0, 0, 0, 0],
         [0, 0, 0, 0],
         [0, 0, 0, 0]],
        ])
    # test default case
    result = clear_border(image.copy())
    ref = image.copy()
    ref[0, 3, 0] = 0
    assert_array_equal(result, ref)

    # test buffer
    result = clear_border(image.copy(), 1)
    assert_array_equal(result, np.zeros(result.shape))

    # test background value
    result = clear_border(image.copy(), buffer_size=1, bgval=2)
    assert_array_equal(result, 2 * np.ones_like(image))


def test_clear_border_non_binary():
    image = np.array([[1, 2, 3, 1, 2],
                      [3, 3, 5, 4, 2],
                      [3, 4, 5, 4, 2],
                      [3, 3, 2, 1, 2]])

    result = clear_border(image)
    expected = np.array([[0, 0, 0, 0, 0],
                         [0, 0, 5, 4, 0],
                         [0, 4, 5, 4, 0],
                         [0, 0, 0, 0, 0]])

    assert_array_equal(result, expected)
    assert_(not np.all(image == result))


def test_clear_border_non_binary_3d():
    image3d = np.array(
        [[[1, 2, 3, 1, 2],
          [3, 3, 3, 4, 2],
          [3, 4, 3, 4, 2],
          [3, 3, 2, 1, 2]],
         [[1, 2, 3, 1, 2],
          [3, 3, 5, 4, 2],
          [3, 4, 5, 4, 2],
          [3, 3, 2, 1, 2]],
         [[1, 2, 3, 1, 2],
          [3, 3, 3, 4, 2],
          [3, 4, 3, 4, 2],
          [3, 3, 2, 1, 2]],
         ])

    result = clear_border(image3d)
    expected = np.array(
        [[[0, 0, 0, 0, 0],
          [0, 0, 0, 0, 0],
          [0, 0, 0, 0, 0],
          [0, 0, 0, 0, 0]],
         [[0, 0, 0, 0, 0],
          [0, 0, 5, 0, 0],
          [0, 0, 5, 0, 0],
          [0, 0, 0, 0, 0]],
         [[0, 0, 0, 0, 0],
          [0, 0, 0, 0, 0],
          [0, 0, 0, 0, 0],
          [0, 0, 0, 0, 0]],
         ])

    assert_array_equal(result, expected)
    assert_(not np.all(image3d == result))


def test_clear_border_non_binary_inplace():
    image = np.array([[1, 2, 3, 1, 2],
                      [3, 3, 5, 4, 2],
                      [3, 4, 5, 4, 2],
                      [3, 3, 2, 1, 2]])

    result = clear_border(image, in_place=True)
    expected = np.array([[0, 0, 0, 0, 0],
                         [0, 0, 5, 4, 0],
                         [0, 4, 5, 4, 0],
                         [0, 0, 0, 0, 0]])

    assert_array_equal(result, expected)
    assert_array_equal(image, result)


def test_clear_border_non_binary_inplace_3d():
    image3d = np.array(
        [[[1, 2, 3, 1, 2],
          [3, 3, 3, 4, 2],
          [3, 4, 3, 4, 2],
          [3, 3, 2, 1, 2]],
         [[1, 2, 3, 1, 2],
          [3, 3, 5, 4, 2],
          [3, 4, 5, 4, 2],
          [3, 3, 2, 1, 2]],
         [[1, 2, 3, 1, 2],
          [3, 3, 3, 4, 2],
          [3, 4, 3, 4, 2],
          [3, 3, 2, 1, 2]],
         ])

    result = clear_border(image3d, in_place=True)
    expected = np.array(
        [[[0, 0, 0, 0, 0],
          [0, 0, 0, 0, 0],
          [0, 0, 0, 0, 0],
          [0, 0, 0, 0, 0]],
         [[0, 0, 0, 0, 0],
          [0, 0, 5, 0, 0],
          [0, 0, 5, 0, 0],
          [0, 0, 0, 0, 0]],
         [[0, 0, 0, 0, 0],
          [0, 0, 0, 0, 0],
          [0, 0, 0, 0, 0],
          [0, 0, 0, 0, 0]],
         ])

    assert_array_equal(result, expected)
    assert_array_equal(image3d, result)

