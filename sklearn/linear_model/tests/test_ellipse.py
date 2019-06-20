import pytest

import numpy as np
from ..ellipse import LsqEllipse

from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_allclose


def make_dataset(center, width, height, phi, epsilon, num_points):
    """Generate Elliptical data with noise

    Returns
    -------
    X : numpy.ndarray
    Y : numpy.ndarray
    """

    t = np.linspace(0, 2 * np.pi, num_points)
    x_noise, y_noise = epsilon * np.random.rand(2, len(t))

    X = (center[0]
         + width * np.cos(t) * np.cos(phi)
         - height * np.sin(t) * np.sin(phi)
         + x_noise / 2.)
    y = (center[1]
         + width * np.cos(t) * np.sin(phi)
         + height * np.sin(t) * np.cos(phi)
         + y_noise / 2.)

    return X, y


# phi needs to be < (1/4 * pi)
# and width cannot = height
# or test is degenerate
@pytest.mark.parametrize('center', [[1, 1], [0, 1]])
@pytest.mark.parametrize('width', [.4, 10])
@pytest.mark.parametrize('height', [.2, 3])
@pytest.mark.parametrize('phi', [np.pi/5, np.pi/13])
def test_ellipse_fit(center, width, height, phi):
    x, y = make_dataset(
        center=center,
        width=width,
        height=height,
        phi=phi,
        epsilon=0,
        num_points=10
    )
    elp = LsqEllipse()
    elp.fit(x, y)
    _center, _width, _height, _phi = elp.get_parameters()

    assert_array_almost_equal(_center, center)
    assert_almost_equal(_width, width)
    assert_almost_equal(_height, height)
    assert_almost_equal(_phi, phi)


def test_minimum_data_points():
    x, y = make_dataset(
        center=[0, 0],
        width=1,
        height=.5,
        phi=0,
        epsilon=0,
        num_points=5
    )
    elp = LsqEllipse()
    elp.fit(x, y)
    _center, _width, _height, _phi = elp.get_parameters()

    assert_array_almost_equal(_center, [0, 0])
    assert_almost_equal(_width, 1)
    assert_almost_equal(_height, .5)
    assert_almost_equal(_phi, 0)


def test_less_than_minimum_data_points_raises_err():
    x, y = make_dataset(
        center=[0, 0],
        width=1,
        height=.5,
        phi=0,
        epsilon=0,
        num_points=4
    )
    elp = LsqEllipse()
    with pytest.raises(ValueError):
        elp.fit(x, y)
