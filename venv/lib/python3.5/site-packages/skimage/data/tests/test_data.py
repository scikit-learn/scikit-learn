import numpy as np
import skimage.data as data
from skimage._shared.testing import assert_equal, assert_almost_equal


def test_astronaut():
    """ Test that "astronaut" image can be loaded. """
    astronaut = data.astronaut()
    assert_equal(astronaut.shape, (512, 512, 3))


def test_camera():
    """ Test that "camera" image can be loaded. """
    cameraman = data.camera()
    assert_equal(cameraman.ndim, 2)


def test_checkerboard():
    """ Test that "checkerboard" image can be loaded. """
    data.checkerboard()


def test_chelsea():
    """ Test that "chelsea" image can be loaded. """
    data.chelsea()


def test_clock():
    """ Test that "clock" image can be loaded. """
    data.clock()


def test_coffee():
    """ Test that "coffee" image can be loaded. """
    data.coffee()


def test_horse():
    """ Test that "horse" image can be loaded. """
    horse = data.horse()
    assert_equal(horse.ndim, 2)
    assert_equal(horse.dtype, np.dtype('bool'))


def test_hubble():
    """ Test that "Hubble" image can be loaded. """
    data.hubble_deep_field()


def test_immunohistochemistry():
    """ Test that "immunohistochemistry" image can be loaded. """
    data.immunohistochemistry()


def test_logo():
    """ Test that "logo" image can be loaded. """
    logo = data.logo()
    assert_equal(logo.ndim, 3)
    assert_equal(logo.shape[2], 4)


def test_moon():
    """ Test that "moon" image can be loaded. """
    data.moon()


def test_page():
    """ Test that "page" image can be loaded. """
    data.page()


def test_rocket():
    """ Test that "rocket" image can be loaded. """
    data.rocket()


def test_text():
    """ Test that "text" image can be loaded. """
    data.text()


def test_stereo_motorcycle():
    """ Test that "stereo_motorcycle" image can be loaded. """
    data.stereo_motorcycle()


def test_binary_blobs():
    blobs = data.binary_blobs(length=128)
    assert_almost_equal(blobs.mean(), 0.5, decimal=1)
    blobs = data.binary_blobs(length=128, volume_fraction=0.25)
    assert_almost_equal(blobs.mean(), 0.25, decimal=1)
    blobs = data.binary_blobs(length=32, volume_fraction=0.25, n_dim=3)
    assert_almost_equal(blobs.mean(), 0.25, decimal=1)
    other_realization = data.binary_blobs(length=32, volume_fraction=0.25,
                                          n_dim=3)
    assert not np.all(blobs == other_realization)


def test_lfw_subset():
    """ Test that "lfw_subset" can be loaded."""
    data.lfw_subset()
