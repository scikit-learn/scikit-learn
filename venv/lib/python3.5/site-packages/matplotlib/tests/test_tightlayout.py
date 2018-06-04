from __future__ import absolute_import, division, print_function

import six
import warnings

import numpy as np

from matplotlib.testing.decorators import image_comparison
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredOffsetbox, DrawingArea
from matplotlib.patches import Rectangle


def example_plot(ax, fontsize=12):
    ax.plot([1, 2])
    ax.locator_params(nbins=3)
    ax.set_xlabel('x-label', fontsize=fontsize)
    ax.set_ylabel('y-label', fontsize=fontsize)
    ax.set_title('Title', fontsize=fontsize)


@image_comparison(baseline_images=['tight_layout1'])
def test_tight_layout1():
    'Test tight_layout for a single subplot'
    fig = plt.figure()
    ax = fig.add_subplot(111)
    example_plot(ax, fontsize=24)
    plt.tight_layout()


@image_comparison(baseline_images=['tight_layout2'])
def test_tight_layout2():
    'Test tight_layout for multiple subplots'
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)
    example_plot(ax1)
    example_plot(ax2)
    example_plot(ax3)
    example_plot(ax4)
    plt.tight_layout()


@image_comparison(baseline_images=['tight_layout3'])
def test_tight_layout3():
    'Test tight_layout for multiple subplots'

    fig = plt.figure()

    ax1 = plt.subplot(221)
    ax2 = plt.subplot(223)
    ax3 = plt.subplot(122)

    example_plot(ax1)
    example_plot(ax2)
    example_plot(ax3)

    plt.tight_layout()


@image_comparison(baseline_images=['tight_layout4'],
                  freetype_version=('2.5.5', '2.6.1'))
def test_tight_layout4():
    'Test tight_layout for subplot2grid'

    fig = plt.figure()

    ax1 = plt.subplot2grid((3, 3), (0, 0))
    ax2 = plt.subplot2grid((3, 3), (0, 1), colspan=2)
    ax3 = plt.subplot2grid((3, 3), (1, 0), colspan=2, rowspan=2)
    ax4 = plt.subplot2grid((3, 3), (1, 2), rowspan=2)

    example_plot(ax1)
    example_plot(ax2)
    example_plot(ax3)
    example_plot(ax4)

    plt.tight_layout()


@image_comparison(baseline_images=['tight_layout5'])
def test_tight_layout5():
    'Test tight_layout for image'

    fig = plt.figure()

    ax = plt.subplot(111)
    arr = np.arange(100).reshape((10, 10))
    ax.imshow(arr, interpolation="none")

    plt.tight_layout()


@image_comparison(baseline_images=['tight_layout6'])
def test_tight_layout6():
    'Test tight_layout for gridspec'

    # This raises warnings since tight layout cannot
    # do this fully automatically. But the test is
    # correct since the layout is manually edited
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig = plt.figure()

        import matplotlib.gridspec as gridspec

        gs1 = gridspec.GridSpec(2, 1)
        ax1 = fig.add_subplot(gs1[0])
        ax2 = fig.add_subplot(gs1[1])

        example_plot(ax1)
        example_plot(ax2)

        gs1.tight_layout(fig, rect=[0, 0, 0.5, 1])

        gs2 = gridspec.GridSpec(3, 1)

        for ss in gs2:
            ax = fig.add_subplot(ss)
            example_plot(ax)
            ax.set_title("")
            ax.set_xlabel("")

        ax.set_xlabel("x-label", fontsize=12)

        gs2.tight_layout(fig, rect=[0.5, 0, 1, 1], h_pad=0.45)

        top = min(gs1.top, gs2.top)
        bottom = max(gs1.bottom, gs2.bottom)

        gs1.tight_layout(fig, rect=[None, 0 + (bottom-gs1.bottom),
                                    0.5, 1 - (gs1.top-top)])
        gs2.tight_layout(fig, rect=[0.5, 0 + (bottom-gs2.bottom),
                                    None, 1 - (gs2.top-top)],
                         h_pad=0.45)


@image_comparison(baseline_images=['tight_layout7'])
def test_tight_layout7():
    # tight layout with left and right titles
    fig = plt.figure()
    fontsize = 24
    ax = fig.add_subplot(111)
    ax.plot([1, 2])
    ax.locator_params(nbins=3)
    ax.set_xlabel('x-label', fontsize=fontsize)
    ax.set_ylabel('y-label', fontsize=fontsize)
    ax.set_title('Left Title', loc='left', fontsize=fontsize)
    ax.set_title('Right Title', loc='right', fontsize=fontsize)
    plt.tight_layout()


@image_comparison(baseline_images=['tight_layout8'])
def test_tight_layout8():
    'Test automatic use of tight_layout'
    fig = plt.figure()
    fig.set_tight_layout({'pad': .1})
    ax = fig.add_subplot(111)
    example_plot(ax, fontsize=24)


@image_comparison(baseline_images=['tight_layout9'])
def test_tight_layout9():
    # Test tight_layout for non-visible suplots
    # GH 8244
    f, axarr = plt.subplots(2, 2)
    axarr[1][1].set_visible(False)
    plt.tight_layout()


# The following test is misleading when the text is removed.
@image_comparison(baseline_images=['outward_ticks'], remove_text=False)
def test_outward_ticks():
    'Test automatic use of tight_layout'
    fig = plt.figure()
    ax = fig.add_subplot(221)
    ax.xaxis.set_tick_params(tickdir='out', length=16, width=3)
    ax.yaxis.set_tick_params(tickdir='out', length=16, width=3)
    ax.xaxis.set_tick_params(
        tickdir='out', length=32, width=3, tick1On=True, which='minor')
    ax.yaxis.set_tick_params(
        tickdir='out', length=32, width=3, tick1On=True, which='minor')
    # The following minor ticks are not labelled, and they
    # are drawn over the major ticks and labels--ugly!
    ax.xaxis.set_ticks([0], minor=True)
    ax.yaxis.set_ticks([0], minor=True)
    ax = fig.add_subplot(222)
    ax.xaxis.set_tick_params(tickdir='in', length=32, width=3)
    ax.yaxis.set_tick_params(tickdir='in', length=32, width=3)
    ax = fig.add_subplot(223)
    ax.xaxis.set_tick_params(tickdir='inout', length=32, width=3)
    ax.yaxis.set_tick_params(tickdir='inout', length=32, width=3)
    ax = fig.add_subplot(224)
    ax.xaxis.set_tick_params(tickdir='out', length=32, width=3)
    ax.yaxis.set_tick_params(tickdir='out', length=32, width=3)
    plt.tight_layout()


def add_offsetboxes(ax, size=10, margin=.1, color='black'):
    """
    Surround ax with OffsetBoxes
    """
    m, mp = margin, 1+margin
    anchor_points = [(-m, -m), (-m, .5), (-m, mp),
                     (mp, .5), (.5, mp), (mp, mp),
                     (.5, -m), (mp, -m), (.5, -m)]
    for point in anchor_points:
        da = DrawingArea(size, size)
        background = Rectangle((0, 0), width=size,
                               height=size,
                               facecolor=color,
                               edgecolor='None',
                               linewidth=0,
                               antialiased=False)
        da.add_artist(background)

        anchored_box = AnchoredOffsetbox(
            loc=10,
            child=da,
            pad=0.,
            frameon=False,
            bbox_to_anchor=point,
            bbox_transform=ax.transAxes,
            borderpad=0.)
        ax.add_artist(anchored_box)
    return anchored_box


@image_comparison(baseline_images=['tight_layout_offsetboxes1',
                                   'tight_layout_offsetboxes2'])
def test_tight_layout_offsetboxes():
    # 1.
    # - Create 4 subplots
    # - Plot a diagonal line on them
    # - Surround each plot with 7 boxes
    # - Use tight_layout
    # - See that the squares are included in the tight_layout
    #   and that the squares in the middle do not overlap
    #
    # 2.
    # - Make the squares around the right side axes invisible
    # - See that the invisible squares do not affect the
    #   tight_layout
    rows = cols = 2
    colors = ['red', 'blue', 'green', 'yellow']
    x = y = [0, 1]

    def _subplots():
        _, axs = plt.subplots(rows, cols)
        axs = axs.flat
        for ax, color in zip(axs, colors):
            ax.plot(x, y, color=color)
            add_offsetboxes(ax, 20, color=color)
        return axs

    # 1.
    axs = _subplots()
    plt.tight_layout()

    # 2.
    axs = _subplots()
    for ax in (axs[cols-1::rows]):
        for child in ax.get_children():
            if isinstance(child, AnchoredOffsetbox):
                child.set_visible(False)

    plt.tight_layout()


def test_empty_layout():
    """Tests that tight layout doesn't cause an error when there are
    no axes.
    """

    fig = plt.gcf()
    fig.tight_layout()
