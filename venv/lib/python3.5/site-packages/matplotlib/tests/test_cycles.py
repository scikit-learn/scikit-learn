import warnings

from matplotlib.testing.decorators import image_comparison
from matplotlib.cbook import MatplotlibDeprecationWarning
import matplotlib.pyplot as plt
import numpy as np
import pytest

from cycler import cycler


@image_comparison(baseline_images=['color_cycle_basic'], remove_text=True,
                  extensions=['png'])
def test_colorcycle_basic():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_prop_cycle(cycler('color', ['r', 'g', 'y']))
    xs = np.arange(10)
    ys = 0.25 * xs + 2
    ax.plot(xs, ys, label='red', lw=4)
    ys = 0.45 * xs + 3
    ax.plot(xs, ys, label='green', lw=4)
    ys = 0.65 * xs + 4
    ax.plot(xs, ys, label='yellow', lw=4)
    ys = 0.85 * xs + 5
    ax.plot(xs, ys, label='red2', lw=4)
    ax.legend(loc='upper left')


@image_comparison(baseline_images=['marker_cycle', 'marker_cycle'],
                  remove_text=True, extensions=['png'])
def test_marker_cycle():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_prop_cycle(cycler('c', ['r', 'g', 'y']) +
                      cycler('marker', ['.', '*', 'x']))
    xs = np.arange(10)
    ys = 0.25 * xs + 2
    ax.plot(xs, ys, label='red dot', lw=4, ms=16)
    ys = 0.45 * xs + 3
    ax.plot(xs, ys, label='green star', lw=4, ms=16)
    ys = 0.65 * xs + 4
    ax.plot(xs, ys, label='yellow x', lw=4, ms=16)
    ys = 0.85 * xs + 5
    ax.plot(xs, ys, label='red2 dot', lw=4, ms=16)
    ax.legend(loc='upper left')

    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Test keyword arguments, numpy arrays, and generic iterators
    ax.set_prop_cycle(c=np.array(['r', 'g', 'y']),
                      marker=iter(['.', '*', 'x']))
    xs = np.arange(10)
    ys = 0.25 * xs + 2
    ax.plot(xs, ys, label='red dot', lw=4, ms=16)
    ys = 0.45 * xs + 3
    ax.plot(xs, ys, label='green star', lw=4, ms=16)
    ys = 0.65 * xs + 4
    ax.plot(xs, ys, label='yellow x', lw=4, ms=16)
    ys = 0.85 * xs + 5
    ax.plot(xs, ys, label='red2 dot', lw=4, ms=16)
    ax.legend(loc='upper left')


@image_comparison(baseline_images=['lineprop_cycle_basic'], remove_text=True,
                  extensions=['png'])
def test_linestylecycle_basic():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_prop_cycle(cycler('ls', ['-', '--', ':']))
    xs = np.arange(10)
    ys = 0.25 * xs + 2
    ax.plot(xs, ys, label='solid', lw=4, color='k')
    ys = 0.45 * xs + 3
    ax.plot(xs, ys, label='dashed', lw=4, color='k')
    ys = 0.65 * xs + 4
    ax.plot(xs, ys, label='dotted', lw=4, color='k')
    ys = 0.85 * xs + 5
    ax.plot(xs, ys, label='solid2', lw=4, color='k')
    ax.legend(loc='upper left')


@image_comparison(baseline_images=['fill_cycle_basic'], remove_text=True,
                  extensions=['png'])
def test_fillcycle_basic():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_prop_cycle(cycler('c',  ['r', 'g', 'y']) +
                      cycler('hatch', ['xx', 'O', '|-']) +
                      cycler('linestyle', ['-', '--', ':']))
    xs = np.arange(10)
    ys = 0.25 * xs**.5 + 2
    ax.fill(xs, ys, label='red, xx', linewidth=3)
    ys = 0.45 * xs**.5 + 3
    ax.fill(xs, ys, label='green, circle', linewidth=3)
    ys = 0.65 * xs**.5 + 4
    ax.fill(xs, ys, label='yellow, cross', linewidth=3)
    ys = 0.85 * xs**.5 + 5
    ax.fill(xs, ys, label='red2, xx', linewidth=3)
    ax.legend(loc='upper left')


@image_comparison(baseline_images=['fill_cycle_ignore'], remove_text=True,
                  extensions=['png'])
def test_fillcycle_ignore():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_prop_cycle(cycler('color',  ['r', 'g', 'y']) +
                      cycler('hatch', ['xx', 'O', '|-']) +
                      cycler('marker', ['.', '*', 'D']))
    xs = np.arange(10)
    ys = 0.25 * xs**.5 + 2
    # Should not advance the cycler, even though there is an
    # unspecified property in the cycler "marker".
    # "marker" is not a Polygon property, and should be ignored.
    ax.fill(xs, ys, 'r', hatch='xx', label='red, xx')
    ys = 0.45 * xs**.5 + 3
    # Allow the cycler to advance, but specify some properties
    ax.fill(xs, ys, hatch='O', label='red, circle')
    ys = 0.65 * xs**.5 + 4
    ax.fill(xs, ys, label='green, circle')
    ys = 0.85 * xs**.5 + 5
    ax.fill(xs, ys, label='yellow, cross')
    ax.legend(loc='upper left')


@image_comparison(baseline_images=['property_collision_plot'],
                  remove_text=True, extensions=['png'])
def test_property_collision_plot():
    fig, ax = plt.subplots()
    ax.set_prop_cycle('linewidth', [2, 4])
    for c in range(1, 4):
        ax.plot(np.arange(10), c * np.arange(10), lw=0.1, color='k')
    ax.plot(np.arange(10), 4 * np.arange(10), color='k')
    ax.plot(np.arange(10), 5 * np.arange(10), color='k')


@image_comparison(baseline_images=['property_collision_fill'],
                  remove_text=True, extensions=['png'])
def test_property_collision_fill():
    fig, ax = plt.subplots()
    xs = np.arange(10)
    ys = 0.25 * xs**.5 + 2
    ax.set_prop_cycle(linewidth=[2, 3, 4, 5, 6], facecolor='bgcmy')
    for c in range(1, 4):
        ax.fill(xs, c * ys, lw=0.1)
    ax.fill(xs, 4 * ys)
    ax.fill(xs, 5 * ys)


def test_valid_input_forms():
    fig, ax = plt.subplots()
    # These should not raise an error.
    ax.set_prop_cycle(None)
    ax.set_prop_cycle(cycler('linewidth', [1, 2]))
    ax.set_prop_cycle('color', 'rgywkbcm')
    ax.set_prop_cycle('lw', (1, 2))
    ax.set_prop_cycle('linewidth', [1, 2])
    ax.set_prop_cycle('linewidth', iter([1, 2]))
    ax.set_prop_cycle('linewidth', np.array([1, 2]))
    ax.set_prop_cycle('color', np.array([[1, 0, 0],
                                         [0, 1, 0],
                                         [0, 0, 1]]))
    ax.set_prop_cycle('dashes', [[], [13, 2], [8, 3, 1, 3], [None, None]])
    ax.set_prop_cycle(lw=[1, 2], color=['k', 'w'], ls=['-', '--'])
    ax.set_prop_cycle(lw=np.array([1, 2]),
                      color=np.array(['k', 'w']),
                      ls=np.array(['-', '--']))
    assert True


def test_cycle_reset():
    fig, ax = plt.subplots()

    # Can't really test a reset because only a cycle object is stored
    # but we can test the first item of the cycle.
    prop = next(ax._get_lines.prop_cycler)
    ax.set_prop_cycle(linewidth=[10, 9, 4])
    assert prop != next(ax._get_lines.prop_cycler)
    ax.set_prop_cycle(None)
    got = next(ax._get_lines.prop_cycler)
    assert prop == got

    fig, ax = plt.subplots()
    # Need to double-check the old set/get_color_cycle(), too
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", MatplotlibDeprecationWarning)
        prop = next(ax._get_lines.prop_cycler)
        ax.set_color_cycle(['c', 'm', 'y', 'k'])
        assert prop != next(ax._get_lines.prop_cycler)
        ax.set_color_cycle(None)
        got = next(ax._get_lines.prop_cycler)
        assert prop == got


def test_invalid_input_forms():
    fig, ax = plt.subplots()

    with pytest.raises((TypeError, ValueError)):
        ax.set_prop_cycle(1)
    with pytest.raises((TypeError, ValueError)):
        ax.set_prop_cycle([1, 2])

    with pytest.raises((TypeError, ValueError)):
        ax.set_prop_cycle('color', 'fish')

    with pytest.raises((TypeError, ValueError)):
        ax.set_prop_cycle('linewidth', 1)
    with pytest.raises((TypeError, ValueError)):
        ax.set_prop_cycle('linewidth', {'1': 1, '2': 2})
    with pytest.raises((TypeError, ValueError)):
        ax.set_prop_cycle(linewidth=1, color='r')

    with pytest.raises((TypeError, ValueError)):
        ax.set_prop_cycle('foobar', [1, 2])
    with pytest.raises((TypeError, ValueError)):
        ax.set_prop_cycle(foobar=[1, 2])

    with pytest.raises((TypeError, ValueError)):
        ax.set_prop_cycle(cycler(foobar=[1, 2]))
    with pytest.raises(ValueError):
        ax.set_prop_cycle(cycler(color='rgb', c='cmy'))
