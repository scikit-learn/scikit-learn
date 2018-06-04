"""
Tests specific to the lines module.
"""
from __future__ import absolute_import, division, print_function

import itertools
import matplotlib.lines as mlines
import pytest
from timeit import repeat
import numpy as np
from cycler import cycler

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.testing.decorators import image_comparison


# Runtimes on a loaded system are inherently flaky. Not so much that a rerun
# won't help, hopefully.
@pytest.mark.flaky(reruns=3)
def test_invisible_Line_rendering():
    """
    Github issue #1256 identified a bug in Line.draw method

    Despite visibility attribute set to False, the draw method was not
    returning early enough and some pre-rendering code was executed
    though not necessary.

    Consequence was an excessive draw time for invisible Line instances
    holding a large number of points (Npts> 10**6)
    """
    # Creates big x and y data:
    N = 10**7
    x = np.linspace(0,1,N)
    y = np.random.normal(size=N)

    # Create a plot figure:
    fig = plt.figure()
    ax = plt.subplot(111)

    # Create a "big" Line instance:
    l = mlines.Line2D(x,y)
    l.set_visible(False)
    # but don't add it to the Axis instance `ax`

    # [here Interactive panning and zooming is pretty responsive]
    # Time the canvas drawing:
    t_no_line = min(repeat(fig.canvas.draw, number=1, repeat=3))
    # (gives about 25 ms)

    # Add the big invisible Line:
    ax.add_line(l)

    # [Now interactive panning and zooming is very slow]
    # Time the canvas drawing:
    t_unvisible_line = min(repeat(fig.canvas.draw, number=1, repeat=3))
    # gives about 290 ms for N = 10**7 pts

    slowdown_factor = (t_unvisible_line/t_no_line)
    slowdown_threshold = 2 # trying to avoid false positive failures
    assert slowdown_factor < slowdown_threshold


def test_set_line_coll_dash():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    np.random.seed(0)
    # Testing setting linestyles for line collections.
    # This should not produce an error.
    cs = ax.contour(np.random.randn(20, 30), linestyles=[(0, (3, 3))])

    assert True


@image_comparison(baseline_images=['line_dashes'], remove_text=True)
def test_line_dashes():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    ax.plot(range(10), linestyle=(0, (3, 3)), lw=5)


def test_line_colors():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(range(10), color='none')
    ax.plot(range(10), color='r')
    ax.plot(range(10), color='.3')
    ax.plot(range(10), color=(1, 0, 0, 1))
    ax.plot(range(10), color=(1, 0, 0))
    fig.canvas.draw()
    assert True


def test_linestyle_variants():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    for ls in ["-", "solid", "--", "dashed",
               "-.", "dashdot", ":", "dotted"]:
        ax.plot(range(10), linestyle=ls)

    fig.canvas.draw()
    assert True


def test_valid_linestyles():
    line = mlines.Line2D([], [])
    with pytest.raises(ValueError):
        line.set_linestyle('aardvark')


@image_comparison(baseline_images=['drawstyle_variants'], remove_text=True,
                  extensions=["png"])
def test_drawstyle_variants():
    fig, axs = plt.subplots(6)
    dss = ["default", "steps-mid", "steps-pre", "steps-post", "steps", None]
    # We want to check that drawstyles are properly handled even for very long
    # lines (for which the subslice optimization is on); however, we need
    # to zoom in so that the difference between the drawstyles is actually
    # visible.
    for ax, ds in zip(axs.flat, dss):
        ax.plot(range(2000), drawstyle=ds)
        ax.set(xlim=(0, 2), ylim=(0, 2))


def test_valid_drawstyles():
    line = mlines.Line2D([], [])
    with pytest.raises(ValueError):
        line.set_drawstyle('foobar')


def test_set_drawstyle():
    x = np.linspace(0, 2*np.pi, 10)
    y = np.sin(x)

    fig, ax = plt.subplots()
    line, = ax.plot(x, y)
    line.set_drawstyle("steps-pre")
    assert len(line.get_path().vertices) == 2*len(x)-1

    line.set_drawstyle("default")
    assert len(line.get_path().vertices) == len(x)


@image_comparison(baseline_images=['line_collection_dashes'], remove_text=True)
def test_set_line_coll_dash_image():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    np.random.seed(0)
    cs = ax.contour(np.random.randn(20, 30), linestyles=[(0, (3, 3))])


@image_comparison(baseline_images=['marker_fill_styles'], remove_text=True,
                  extensions=['png'])
def test_marker_fill_styles():
    colors = itertools.cycle([[0, 0, 1], 'g', '#ff0000', 'c', 'm', 'y',
                              np.array([0, 0, 0])])
    altcolor = 'lightgreen'

    y = np.array([1, 1])
    x = np.array([0, 9])
    fig, ax = plt.subplots()

    for j, marker in enumerate(mlines.Line2D.filled_markers):
        for i, fs in enumerate(mlines.Line2D.fillStyles):
            color = next(colors)
            ax.plot(j * 10 + x, y + i + .5 * (j % 2),
                    marker=marker,
                    markersize=20,
                    markerfacecoloralt=altcolor,
                    fillstyle=fs,
                    label=fs,
                    linewidth=5,
                    color=color,
                    markeredgecolor=color,
                    markeredgewidth=2)

    ax.set_ylim([0, 7.5])
    ax.set_xlim([-5, 155])


@image_comparison(baseline_images=['scaled_lines'], style='default')
def test_lw_scaling():
    th = np.linspace(0, 32)
    fig, ax = plt.subplots()
    lins_styles = ['dashed', 'dotted', 'dashdot']
    cy = cycler(matplotlib.rcParams['axes.prop_cycle'])
    for j, (ls, sty) in enumerate(zip(lins_styles, cy)):
        for lw in np.linspace(.5, 10, 10):
            ax.plot(th, j*np.ones(50) + .1 * lw, linestyle=ls, lw=lw, **sty)


def test_nan_is_sorted():
    line = mlines.Line2D([], [])
    assert line._is_sorted(np.array([1, 2, 3]))
    assert line._is_sorted(np.array([1, np.nan, 3]))
    assert not line._is_sorted([3, 5] + [np.nan] * 100 + [0, 2])
