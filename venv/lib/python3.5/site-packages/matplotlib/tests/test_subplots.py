from __future__ import absolute_import, division, print_function

import warnings

import numpy
import matplotlib.pyplot as plt
from matplotlib.testing.decorators import image_comparison

import pytest


def check_shared(axs, x_shared, y_shared):
    """
    x_shared and y_shared are n x n boolean matrices; entry (i, j) indicates
    whether the x (or y) axes of subplots i and j should be shared.
    """
    shared = [axs[0]._shared_x_axes, axs[0]._shared_y_axes]
    for (i1, ax1), (i2, ax2), (i3, (name, shared)) in zip(
            enumerate(axs),
            enumerate(axs),
            enumerate(zip("xy", [x_shared, y_shared]))):
        if i2 <= i1:
            continue
        assert shared[i3].joined(ax1, ax2) == shared[i1, i2], \
            "axes %i and %i incorrectly %ssharing %s axis" % (
                i1, i2, "not " if shared[i1, i2] else "", name)


def check_visible(axs, x_visible, y_visible):
    def tostr(v):
        return "invisible" if v else "visible"

    for (ax, vx, vy) in zip(axs, x_visible, y_visible):
        for l in ax.get_xticklabels() + [ax.get_xaxis().offsetText]:
            assert l.get_visible() == vx, \
                    "X axis was incorrectly %s" % (tostr(vx))
        for l in ax.get_yticklabels() + [ax.get_yaxis().offsetText]:
            assert l.get_visible() == vy, \
                    "Y axis was incorrectly %s" % (tostr(vy))


def test_shared():
    rdim = (4, 4, 2)
    share = {
            'all': numpy.ones(rdim[:2], dtype=bool),
            'none': numpy.zeros(rdim[:2], dtype=bool),
            'row': numpy.array([
                [False, True, False, False],
                [True, False, False, False],
                [False, False, False, True],
                [False, False, True, False]]),
            'col': numpy.array([
                [False, False, True, False],
                [False, False, False, True],
                [True, False, False, False],
                [False, True, False, False]]),
            }
    visible = {
            'x': {
                'all': [False, False, True, True],
                'col': [False, False, True, True],
                'row': [True] * 4,
                'none': [True] * 4,
                False: [True] * 4,
                True: [False, False, True, True],
                },
            'y': {
                'all': [True, False, True, False],
                'col': [True] * 4,
                'row': [True, False, True, False],
                'none': [True] * 4,
                False: [True] * 4,
                True: [True, False, True, False],
                },
            }
    share[False] = share['none']
    share[True] = share['all']

    # test default
    f, ((a1, a2), (a3, a4)) = plt.subplots(2, 2)
    axs = [a1, a2, a3, a4]
    check_shared(axs, share['none'], share['none'])
    plt.close(f)

    # test all option combinations
    ops = [False, True, 'all', 'none', 'row', 'col']
    for xo in ops:
        for yo in ops:
            f, ((a1, a2), (a3, a4)) = plt.subplots(2, 2, sharex=xo, sharey=yo)
            axs = [a1, a2, a3, a4]
            check_shared(axs, share[xo], share[yo])
            check_visible(axs, visible['x'][xo], visible['y'][yo])
            plt.close(f)

    # test label_outer
    f, ((a1, a2), (a3, a4)) = plt.subplots(2, 2, sharex=True, sharey=True)
    axs = [a1, a2, a3, a4]
    for ax in axs:
        ax.label_outer()
    check_visible(axs, [False, False, True, True], [True, False, True, False])


def test_shared_and_moved():
    # test if sharey is on, but then tick_left is called that labels don't
    # re-appear.  Seaborn does this just to be sure yaxis is on left...
    f, (a1, a2) = plt.subplots(1, 2, sharey=True)
    check_visible([a2], [True], [False])
    a2.yaxis.tick_left()
    check_visible([a2], [True], [False])

    f, (a1, a2) = plt.subplots(2, 1, sharex=True)
    check_visible([a1], [False], [True])
    a2.xaxis.tick_bottom()
    check_visible([a1], [False], [True])


def test_exceptions():
    # TODO should this test more options?
    with pytest.raises(ValueError):
        plt.subplots(2, 2, sharex='blah')
    with pytest.raises(ValueError):
        plt.subplots(2, 2, sharey='blah')
    # We filter warnings in this test which are genuine since
    # the point of this test is to ensure that this raises.
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',
                                message='.*sharex argument to subplots',
                                category=UserWarning)
        with pytest.raises(ValueError):
            plt.subplots(2, 2, -1)
        with pytest.raises(ValueError):
            plt.subplots(2, 2, 0)
        with pytest.raises(ValueError):
            plt.subplots(2, 2, 5)


@image_comparison(baseline_images=['subplots_offset_text'], remove_text=False)
def test_subplots_offsettext():
    x = numpy.arange(0, 1e10, 1e9)
    y = numpy.arange(0, 100, 10)+1e4
    fig, axes = plt.subplots(2, 2, sharex='col', sharey='all')
    axes[0, 0].plot(x, x)
    axes[1, 0].plot(x, x)
    axes[0, 1].plot(y, x)
    axes[1, 1].plot(y, x)
