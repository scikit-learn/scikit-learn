import numpy as np
from skimage.segmentation import chan_vese

from skimage._shared import testing
from skimage._shared.testing import assert_array_equal


def test_chan_vese_flat_level_set():
    # because the algorithm evolves the level set around the
    # zero-level, it the level-set has no zero level, the algorithm
    # will not produce results in theory. However, since a continuous
    # approximation of the delta function is used, the algorithm
    # still affects the entirety of the level-set. Therefore with
    # infinite time, the segmentation will still converge.
    img = np.zeros((10, 10))
    img[3:6, 3:6] = np.ones((3, 3))
    ls = np.ones((10, 10)) * 1000
    result = chan_vese(img, mu=0.0, tol=1e-3, init_level_set=ls)
    assert_array_equal(result.astype(np.float), np.ones((10, 10)))
    result = chan_vese(img, mu=0.0, tol=1e-3, init_level_set=-ls)
    assert_array_equal(result.astype(np.float), np.zeros((10, 10)))


def test_chan_vese_small_disk_level_set():
    img = np.zeros((10, 10))
    img[3:6, 3:6] = np.ones((3, 3))
    result = chan_vese(img, mu=0.0, tol=1e-3, init_level_set="small disk")
    assert_array_equal(result.astype(np.float), img)


def test_chan_vese_simple_shape():
    img = np.zeros((10, 10))
    img[3:6, 3:6] = np.ones((3, 3))
    result = chan_vese(img, mu=0.0, tol=1e-8).astype(np.float)
    assert_array_equal(result, img)


def test_chan_vese_extended_output():
    img = np.zeros((10, 10))
    img[3:6, 3:6] = np.ones((3, 3))
    result = chan_vese(img, mu=0.0, tol=1e-8, extended_output=True)
    assert_array_equal(len(result), 3)


def test_chan_vese_remove_noise():
    ref = np.zeros((10, 10))
    ref[1:6, 1:6] = np.array([[0, 1, 1, 1, 0],
                              [1, 1, 1, 1, 1],
                              [1, 1, 1, 1, 1],
                              [1, 1, 1, 1, 1],
                              [0, 1, 1, 1, 0]])
    img = ref.copy()
    img[8, 3] = 1
    result = chan_vese(img, mu=0.3, tol=1e-3, max_iter=100, dt=10,
                       init_level_set="disk").astype(np.float)
    assert_array_equal(result, ref)


def test_chan_vese_incorrect_image_type():
    img = np.zeros((10, 10, 3))
    ls = np.zeros((10, 9))
    with testing.raises(ValueError):
        chan_vese(img, mu=0.0, init_level_set=ls)


def test_chan_vese_gap_closing():
    ref = np.zeros((20, 20))
    ref[8:15, :] = np.ones((7, 20))
    img = ref.copy()
    img[:, 6] = np.zeros((20))
    result = chan_vese(img, mu=0.7, tol=1e-3, max_iter=1000, dt=1000,
                       init_level_set="disk").astype(np.float)
    assert_array_equal(result, ref)


def test_chan_vese_incorrect_level_set():
    img = np.zeros((10, 10))
    ls = np.zeros((10, 9))
    with testing.raises(ValueError):
        chan_vese(img, mu=0.0, init_level_set=ls)
    with testing.raises(ValueError):
        chan_vese(img, mu=0.0, init_level_set="a")


def test_chan_vese_blank_image():
    img = np.zeros((10, 10))
    level_set = np.random.rand(10, 10)
    ref = level_set > 0
    result = chan_vese(img, mu=0.0, tol=0.0, init_level_set=level_set)
    assert_array_equal(result, ref)
