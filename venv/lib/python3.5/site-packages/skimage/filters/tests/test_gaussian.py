import numpy as np
from skimage.filters._gaussian import gaussian
from skimage._shared import testing
from skimage._shared._warnings import expected_warnings


def test_negative_sigma():
    a = np.zeros((3, 3))
    a[1, 1] = 1.
    with testing.raises(ValueError):
        gaussian(a, sigma=-1.0)
    with testing.raises(ValueError):
        gaussian(a, sigma=[-1.0, 1.0])
    with testing.raises(ValueError):
        gaussian(a,
                 sigma=np.asarray([-1.0, 1.0]))


def test_null_sigma():
    a = np.zeros((3, 3))
    a[1, 1] = 1.
    assert np.all(gaussian(a, 0) == a)


def test_default_sigma():
    a = np.zeros((3, 3))
    a[1, 1] = 1.
    assert np.all(gaussian(a) == gaussian(a, sigma=1))


def test_energy_decrease():
    a = np.zeros((3, 3))
    a[1, 1] = 1.
    gaussian_a = gaussian(a, sigma=1, mode='reflect')
    assert gaussian_a.std() < a.std()


def test_multichannel():
    a = np.zeros((5, 5, 3))
    a[1, 1] = np.arange(1, 4)
    gaussian_rgb_a = gaussian(a, sigma=1, mode='reflect',
                              multichannel=True)
    # Check that the mean value is conserved in each channel
    # (color channels are not mixed together)
    assert np.allclose([a[..., i].mean() for i in range(3)],
                       [gaussian_rgb_a[..., i].mean() for i in range(3)])
    # Test multichannel = None
    with expected_warnings(['multichannel']):
        gaussian_rgb_a = gaussian(a, sigma=1, mode='reflect')
    # Check that the mean value is conserved in each channel
    # (color channels are not mixed together)
    assert np.allclose([a[..., i].mean() for i in range(3)],
                       [gaussian_rgb_a[..., i].mean() for i in range(3)])
    # Iterable sigma
    gaussian_rgb_a = gaussian(a, sigma=[1, 2], mode='reflect',
                              multichannel=True)
    assert np.allclose([a[..., i].mean() for i in range(3)],
                       [gaussian_rgb_a[..., i].mean() for i in range(3)])


def test_preserve_range():
    img = np.array([[10.0, -10.0], [-4, 3]], dtype=np.float32)
    gaussian(img, 1, preserve_range=True)
