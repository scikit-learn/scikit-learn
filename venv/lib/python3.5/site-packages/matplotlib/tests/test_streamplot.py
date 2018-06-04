from __future__ import absolute_import, division, print_function

import sys

import numpy as np
from numpy.testing import assert_array_almost_equal
import matplotlib.pyplot as plt
from matplotlib.testing.decorators import image_comparison
import matplotlib.transforms as mtransforms


on_win = (sys.platform == 'win32')


def velocity_field():
    Y, X = np.mgrid[-3:3:100j, -3:3:100j]
    U = -1 - X**2 + Y
    V = 1 + X - Y**2
    return X, Y, U, V


def swirl_velocity_field():
    x = np.linspace(-3., 3., 100)
    y = np.linspace(-3., 3., 100)
    X, Y = np.meshgrid(x, y)
    a = 0.1
    U = np.cos(a) * (-Y) - np.sin(a) * X
    V = np.sin(a) * (-Y) + np.cos(a) * X
    return x, y, U, V


@image_comparison(baseline_images=['streamplot_startpoints'])
def test_startpoints():
    X, Y, U, V = velocity_field()
    start_x = np.linspace(X.min(), X.max(), 10)
    start_y = np.linspace(Y.min(), Y.max(), 10)
    start_points = np.column_stack([start_x, start_y])
    plt.streamplot(X, Y, U, V, start_points=start_points)
    plt.plot(start_x, start_y, 'ok')


@image_comparison(baseline_images=['streamplot_colormap'],
                  tol=.02)
def test_colormap():
    X, Y, U, V = velocity_field()
    plt.streamplot(X, Y, U, V, color=U, density=0.6, linewidth=2,
                   cmap=plt.cm.autumn)
    plt.colorbar()


@image_comparison(baseline_images=['streamplot_linewidth'])
def test_linewidth():
    X, Y, U, V = velocity_field()
    speed = np.sqrt(U*U + V*V)
    lw = 5*speed/speed.max()
    df = 25. / 30.   # Compatibility factor for old test image
    plt.streamplot(X, Y, U, V, density=[0.5 * df, 1. * df], color='k',
                   linewidth=lw)


@image_comparison(baseline_images=['streamplot_masks_and_nans'],
                  tol=0.04 if on_win else 0)
def test_masks_and_nans():
    X, Y, U, V = velocity_field()
    mask = np.zeros(U.shape, dtype=bool)
    mask[40:60, 40:60] = 1
    U[:20, :20] = np.nan
    U = np.ma.array(U, mask=mask)
    with np.errstate(invalid='ignore'):
        plt.streamplot(X, Y, U, V, color=U, cmap=plt.cm.Blues)


@image_comparison(baseline_images=['streamplot_maxlength'],
                  extensions=['png'])
def test_maxlength():
    x, y, U, V = swirl_velocity_field()
    plt.streamplot(x, y, U, V, maxlength=10., start_points=[[0., 1.5]],
                   linewidth=2, density=2)


@image_comparison(baseline_images=['streamplot_direction'],
                  extensions=['png'])
def test_direction():
    x, y, U, V = swirl_velocity_field()
    plt.streamplot(x, y, U, V, integration_direction='backward',
                   maxlength=1.5, start_points=[[1.5, 0.]],
                   linewidth=2, density=2)


def test_streamplot_limits():
    ax = plt.axes()
    x = np.linspace(-5, 10, 20)
    y = np.linspace(-2, 4, 10)
    y, x = np.meshgrid(y, x)
    trans = mtransforms.Affine2D().translate(25, 32) + ax.transData
    plt.barbs(x, y, np.sin(x), np.cos(y), transform=trans)
    # The calculated bounds are approximately the bounds of the original data,
    # this is because the entire path is taken into account when updating the
    # datalim.
    assert_array_almost_equal(ax.dataLim.bounds, (20, 30, 15, 6),
                              decimal=1)
