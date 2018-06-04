from __future__ import print_function
import warnings
import numpy as np
import pytest
import sys
from matplotlib import pyplot as plt
from matplotlib.testing.decorators import image_comparison


def draw_quiver(ax, **kw):
    X, Y = np.meshgrid(np.arange(0, 2 * np.pi, 1),
                       np.arange(0, 2 * np.pi, 1))
    U = np.cos(X)
    V = np.sin(Y)

    Q = ax.quiver(U, V, **kw)
    return Q


def test_quiver_memory_leak():
    fig, ax = plt.subplots()

    Q = draw_quiver(ax)
    ttX = Q.X
    Q.remove()

    del Q

    assert sys.getrefcount(ttX) == 2


def test_quiver_key_memory_leak():
    fig, ax = plt.subplots()

    Q = draw_quiver(ax)

    qk = ax.quiverkey(Q, 0.5, 0.92, 2, r'$2 \frac{m}{s}$',
                      labelpos='W',
                      fontproperties={'weight': 'bold'})
    assert sys.getrefcount(qk) == 3
    qk.remove()
    assert sys.getrefcount(qk) == 2


def test_no_warnings():
    fig, ax = plt.subplots()

    X, Y = np.meshgrid(np.arange(15), np.arange(10))
    U = V = np.ones_like(X)

    phi = (np.random.rand(15, 10) - .5) * 150
    with warnings.catch_warnings(record=True) as w:
        ax.quiver(X, Y, U, V, angles=phi)
        fig.canvas.draw()
    assert len(w) == 0


def test_zero_headlength():
    # Based on report by Doug McNeil:
    # http://matplotlib.1069221.n5.nabble.com/quiver-warnings-td28107.html
    fig, ax = plt.subplots()
    X, Y = np.meshgrid(np.arange(10), np.arange(10))
    U, V = np.cos(X), np.sin(Y)
    with warnings.catch_warnings(record=True) as w:
        ax.quiver(U, V, headlength=0, headaxislength=0)
        fig.canvas.draw()
    assert len(w) == 0


@image_comparison(baseline_images=['quiver_animated_test_image'],
                  extensions=['png'])
def test_quiver_animate():
    # Tests fix for #2616
    fig, ax = plt.subplots()

    Q = draw_quiver(ax, animated=True)

    qk = ax.quiverkey(Q, 0.5, 0.92, 2, r'$2 \frac{m}{s}$',
                      labelpos='W',
                      fontproperties={'weight': 'bold'})


@image_comparison(baseline_images=['quiver_with_key_test_image'],
                  extensions=['png'])
def test_quiver_with_key():
    fig, ax = plt.subplots()
    ax.margins(0.1)

    Q = draw_quiver(ax)

    qk = ax.quiverkey(Q, 0.5, 0.95, 2,
                      r'$2\, \mathrm{m}\, \mathrm{s}^{-1}$',
                      angle=-10,
                      coordinates='figure',
                      labelpos='W',
                      fontproperties={'weight': 'bold',
                                      'size': 'large'})


@image_comparison(baseline_images=['quiver_single_test_image'],
                  extensions=['png'], remove_text=True)
def test_quiver_single():
    fig, ax = plt.subplots()
    ax.margins(0.1)

    ax.quiver([1], [1], [2], [2])


def test_quiver_copy():
    fig, ax = plt.subplots()
    uv = dict(u=np.array([1.1]), v=np.array([2.0]))
    q0 = ax.quiver([1], [1], uv['u'], uv['v'])
    uv['v'][0] = 0
    assert q0.V[0] == 2.0


@image_comparison(baseline_images=['quiver_key_pivot'],
                  extensions=['png'], remove_text=True)
def test_quiver_key_pivot():
    fig, ax = plt.subplots()

    u, v = np.mgrid[0:2*np.pi:10j, 0:2*np.pi:10j]

    q = ax.quiver(np.sin(u), np.cos(v))
    ax.set_xlim(-2, 11)
    ax.set_ylim(-2, 11)
    ax.quiverkey(q, 0.5, 1, 1, 'N', labelpos='N')
    ax.quiverkey(q, 1, 0.5, 1, 'E', labelpos='E')
    ax.quiverkey(q, 0.5, 0, 1, 'S', labelpos='S')
    ax.quiverkey(q, 0, 0.5, 1, 'W', labelpos='W')


@image_comparison(baseline_images=['barbs_test_image'],
                  extensions=['png'], remove_text=True)
def test_barbs():
    x = np.linspace(-5, 5, 5)
    X, Y = np.meshgrid(x, x)
    U, V = 12*X, 12*Y
    fig, ax = plt.subplots()
    ax.barbs(X, Y, U, V, np.sqrt(U*U + V*V), fill_empty=True, rounding=False,
             sizes=dict(emptybarb=0.25, spacing=0.2, height=0.3),
             cmap='viridis')


@image_comparison(baseline_images=['barbs_pivot_test_image'],
                  extensions=['png'], remove_text=True)
def test_barbs_pivot():
    x = np.linspace(-5, 5, 5)
    X, Y = np.meshgrid(x, x)
    U, V = 12*X, 12*Y
    fig, ax = plt.subplots()
    ax.barbs(X, Y, U, V, fill_empty=True, rounding=False, pivot=1.7,
             sizes=dict(emptybarb=0.25, spacing=0.2, height=0.3))
    ax.scatter(X, Y, s=49, c='black')


def test_bad_masked_sizes():
    'Test error handling when given differing sized masked arrays'
    x = np.arange(3)
    y = np.arange(3)
    u = np.ma.array(15. * np.ones((4,)))
    v = np.ma.array(15. * np.ones_like(u))
    u[1] = np.ma.masked
    v[1] = np.ma.masked
    fig, ax = plt.subplots()
    with pytest.raises(ValueError):
        ax.barbs(x, y, u, v)


def test_angles_and_scale():
    # angles array + scale_units kwarg
    fig, ax = plt.subplots()
    X, Y = np.meshgrid(np.arange(15), np.arange(10))
    U = V = np.ones_like(X)
    phi = (np.random.rand(15, 10) - .5) * 150
    ax.quiver(X, Y, U, V, angles=phi, scale_units='xy')


@image_comparison(baseline_images=['quiver_xy'],
                  extensions=['png'], remove_text=True)
def test_quiver_xy():
    # simple arrow pointing from SW to NE
    fig, ax = plt.subplots(subplot_kw=dict(aspect='equal'))
    ax.quiver(0, 0, 1, 1, angles='xy', scale_units='xy', scale=1)
    ax.set_xlim(0, 1.1)
    ax.set_ylim(0, 1.1)
    ax.grid()


def test_quiverkey_angles():
    # Check that only a single arrow is plotted for a quiverkey when an array
    # of angles is given to the original quiver plot
    fig, ax = plt.subplots()

    X, Y = np.meshgrid(np.arange(2), np.arange(2))
    U = V = angles = np.ones_like(X)

    q = ax.quiver(X, Y, U, V, angles=angles)
    qk = ax.quiverkey(q, 1, 1, 2, 'Label')
    # The arrows are only created when the key is drawn
    fig.canvas.draw()
    assert len(qk.verts) == 1
