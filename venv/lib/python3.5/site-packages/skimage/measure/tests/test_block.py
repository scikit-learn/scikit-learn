import numpy as np
from skimage.measure import block_reduce

from skimage._shared import testing
from skimage._shared.testing import assert_equal


def test_block_reduce_sum():
    image1 = np.arange(4 * 6).reshape(4, 6)
    out1 = block_reduce(image1, (2, 3))
    expected1 = np.array([[ 24,  42],
                          [ 96, 114]])
    assert_equal(expected1, out1)

    image2 = np.arange(5 * 8).reshape(5, 8)
    out2 = block_reduce(image2, (3, 3))
    expected2 = np.array([[ 81, 108,  87],
                          [174, 192, 138]])
    assert_equal(expected2, out2)


def test_block_reduce_mean():
    image1 = np.arange(4 * 6).reshape(4, 6)
    out1 = block_reduce(image1, (2, 3), func=np.mean)
    expected1 = np.array([[  4.,   7.],
                          [ 16.,  19.]])
    assert_equal(expected1, out1)

    image2 = np.arange(5 * 8).reshape(5, 8)
    out2 = block_reduce(image2, (4, 5), func=np.mean)
    expected2 = np.array([[14. , 10.8],
                          [ 8.5,  5.7]])
    assert_equal(expected2, out2)


def test_block_reduce_median():
    image1 = np.arange(4 * 6).reshape(4, 6)
    out1 = block_reduce(image1, (2, 3), func=np.median)
    expected1 = np.array([[  4.,   7.],
                          [ 16.,  19.]])
    assert_equal(expected1, out1)

    image2 = np.arange(5 * 8).reshape(5, 8)
    out2 = block_reduce(image2, (4, 5), func=np.median)
    expected2 = np.array([[ 14.,  17.],
                          [  0.,   0.]])
    assert_equal(expected2, out2)

    image3 = np.array([[1, 5, 5, 5], [5, 5, 5, 1000]])
    out3 = block_reduce(image3, (2, 4), func=np.median)
    assert_equal(5, out3)


def test_block_reduce_min():
    image1 = np.arange(4 * 6).reshape(4, 6)
    out1 = block_reduce(image1, (2, 3), func=np.min)
    expected1 = np.array([[ 0, 3],
                          [12, 15]])
    assert_equal(expected1, out1)

    image2 = np.arange(5 * 8).reshape(5, 8)
    out2 = block_reduce(image2, (4, 5), func=np.min)
    expected2 = np.array([[0, 0],
                          [0, 0]])
    assert_equal(expected2, out2)


def test_block_reduce_max():
    image1 = np.arange(4 * 6).reshape(4, 6)
    out1 = block_reduce(image1, (2, 3), func=np.max)
    expected1 = np.array([[ 8, 11],
                          [20, 23]])
    assert_equal(expected1, out1)

    image2 = np.arange(5 * 8).reshape(5, 8)
    out2 = block_reduce(image2, (4, 5), func=np.max)
    expected2 = np.array([[28, 31],
                          [36, 39]])
    assert_equal(expected2, out2)


def test_invalid_block_size():
    image = np.arange(4 * 6).reshape(4, 6)

    with testing.raises(ValueError):
        block_reduce(image, [1, 2, 3])
    with testing.raises(ValueError):
        block_reduce(image, [1, 0.5])
