from __future__ import absolute_import, division, print_function

import sys
import warnings

from matplotlib import rcParams
from matplotlib.testing.decorators import image_comparison
from matplotlib.axes import Axes
from matplotlib.ticker import AutoMinorLocator, FixedFormatter
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
import numpy as np
import pytest


@image_comparison(baseline_images=['figure_align_labels'])
def test_align_labels():
    # Check the figure.align_labels() command
    fig = plt.figure(tight_layout=True)
    gs = gridspec.GridSpec(3, 3)

    ax = fig.add_subplot(gs[0, :2])
    ax.plot(np.arange(0, 1e6, 1000))
    ax.set_ylabel('Ylabel0 0')
    ax = fig.add_subplot(gs[0, -1])
    ax.plot(np.arange(0, 1e4, 100))

    for i in range(3):
        ax = fig.add_subplot(gs[1, i])
        ax.set_ylabel('YLabel1 %d' % i)
        ax.set_xlabel('XLabel1 %d' % i)
        if i in [0, 2]:
            ax.xaxis.set_label_position("top")
            ax.xaxis.tick_top()
        if i == 0:
            for tick in ax.get_xticklabels():
                tick.set_rotation(90)
        if i == 2:
            ax.yaxis.set_label_position("right")
            ax.yaxis.tick_right()

    for i in range(3):
        ax = fig.add_subplot(gs[2, i])
        ax.set_xlabel('XLabel2 %d' % (i))
        ax.set_ylabel('YLabel2 %d' % (i))

        if i == 2:
            ax.plot(np.arange(0, 1e4, 10))
            ax.yaxis.set_label_position("right")
            ax.yaxis.tick_right()
            for tick in ax.get_xticklabels():
                tick.set_rotation(90)

    fig.align_labels()


def test_figure_label():
    # pyplot figure creation, selection and closing with figure label and
    # number
    plt.close('all')
    plt.figure('today')
    plt.figure(3)
    plt.figure('tomorrow')
    plt.figure()
    plt.figure(0)
    plt.figure(1)
    plt.figure(3)
    assert plt.get_fignums() == [0, 1, 3, 4, 5]
    assert plt.get_figlabels() == ['', 'today', '', 'tomorrow', '']
    plt.close(10)
    plt.close()
    plt.close(5)
    plt.close('tomorrow')
    assert plt.get_fignums() == [0, 1]
    assert plt.get_figlabels() == ['', 'today']


def test_fignum_exists():
    # pyplot figure creation, selection and closing with fignum_exists
    plt.figure('one')
    plt.figure(2)
    plt.figure('three')
    plt.figure()
    assert plt.fignum_exists('one')
    assert plt.fignum_exists(2)
    assert plt.fignum_exists('three')
    assert plt.fignum_exists(4)
    plt.close('one')
    plt.close(4)
    assert not plt.fignum_exists('one')
    assert not plt.fignum_exists(4)


def test_clf_keyword():
    # test if existing figure is cleared with figure() and subplots()
    text1 = 'A fancy plot'
    text2 = 'Really fancy!'

    fig0 = plt.figure(num=1)
    fig0.suptitle(text1)
    assert [t.get_text() for t in fig0.texts] == [text1]

    fig1 = plt.figure(num=1, clear=False)
    fig1.text(0.5, 0.5, text2)
    assert fig0 is fig1
    assert [t.get_text() for t in fig1.texts] == [text1, text2]

    fig2, ax2 = plt.subplots(2, 1, num=1, clear=True)
    assert fig0 is fig2
    assert [t.get_text() for t in fig2.texts] == []


@image_comparison(baseline_images=['figure_today'])
def test_figure():
    # named figure support
    fig = plt.figure('today')
    ax = fig.add_subplot(111)
    ax.set_title(fig.get_label())
    ax.plot(np.arange(5))
    # plot red line in a different figure.
    plt.figure('tomorrow')
    plt.plot([0, 1], [1, 0], 'r')
    # Return to the original; make sure the red line is not there.
    plt.figure('today')
    plt.close('tomorrow')


@image_comparison(baseline_images=['figure_legend'])
def test_figure_legend():
    fig, axes = plt.subplots(2)
    axes[0].plot([0, 1], [1, 0], label='x', color='g')
    axes[0].plot([0, 1], [0, 1], label='y', color='r')
    axes[0].plot([0, 1], [0.5, 0.5], label='y', color='k')

    axes[1].plot([0, 1], [1, 0], label='_y', color='r')
    axes[1].plot([0, 1], [0, 1], label='z', color='b')
    fig.legend()


def test_gca():
    fig = plt.figure()

    ax1 = fig.add_axes([0, 0, 1, 1])
    assert fig.gca(projection='rectilinear') is ax1
    assert fig.gca() is ax1

    ax2 = fig.add_subplot(121, projection='polar')
    assert fig.gca() is ax2
    assert fig.gca(polar=True)is ax2

    ax3 = fig.add_subplot(122)
    assert fig.gca() is ax3

    # the final request for a polar axes will end up creating one
    # with a spec of 111.
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        # Changing the projection will throw a warning
        assert fig.gca(polar=True) is not ax3
        assert len(w) == 1
    assert fig.gca(polar=True) is not ax2
    assert fig.gca().get_geometry() == (1, 1, 1)

    fig.sca(ax1)
    assert fig.gca(projection='rectilinear') is ax1
    assert fig.gca() is ax1


@image_comparison(baseline_images=['figure_suptitle'])
def test_suptitle():
    fig, _ = plt.subplots()
    fig.suptitle('hello', color='r')
    fig.suptitle('title', color='g', rotation='30')


def test_suptitle_fontproperties():
    from matplotlib.font_manager import FontProperties
    fig, ax = plt.subplots()
    fps = FontProperties(size='large', weight='bold')
    txt = fig.suptitle('fontprops title', fontproperties=fps)
    assert txt.get_fontsize() == fps.get_size_in_points()
    assert txt.get_weight() == fps.get_weight()


@image_comparison(baseline_images=['alpha_background'],
                  # only test png and svg. The PDF output appears correct,
                  # but Ghostscript does not preserve the background color.
                  extensions=['png', 'svg'],
                  savefig_kwarg={'facecolor': (0, 1, 0.4),
                                 'edgecolor': 'none'})
def test_alpha():
    # We want an image which has a background color and an
    # alpha of 0.4.
    fig = plt.figure(figsize=[2, 1])
    fig.set_facecolor((0, 1, 0.4))
    fig.patch.set_alpha(0.4)

    import matplotlib.patches as mpatches
    fig.patches.append(mpatches.CirclePolygon([20, 20],
                                              radius=15,
                                              alpha=0.6,
                                              facecolor='red'))


def test_too_many_figures():
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        for i in range(rcParams['figure.max_open_warning'] + 1):
            plt.figure()
        assert len(w) == 1


def test_iterability_axes_argument():

    # This is a regression test for matplotlib/matplotlib#3196. If one of the
    # arguments returned by _as_mpl_axes defines __getitem__ but is not
    # iterable, this would raise an execption. This is because we check
    # whether the arguments are iterable, and if so we try and convert them
    # to a tuple. However, the ``iterable`` function returns True if
    # __getitem__ is present, but some classes can define __getitem__ without
    # being iterable. The tuple conversion is now done in a try...except in
    # case it fails.

    class MyAxes(Axes):
        def __init__(self, *args, **kwargs):
            kwargs.pop('myclass', None)
            return Axes.__init__(self, *args, **kwargs)

    class MyClass(object):

        def __getitem__(self, item):
            if item != 'a':
                raise ValueError("item should be a")

        def _as_mpl_axes(self):
            return MyAxes, {'myclass': self}

    fig = plt.figure()
    fig.add_subplot(1, 1, 1, projection=MyClass())
    plt.close(fig)


def test_set_fig_size():
    fig = plt.figure()

    # check figwidth
    fig.set_figwidth(5)
    assert fig.get_figwidth() == 5

    # check figheight
    fig.set_figheight(1)
    assert fig.get_figheight() == 1

    # check using set_size_inches
    fig.set_size_inches(2, 4)
    assert fig.get_figwidth() == 2
    assert fig.get_figheight() == 4

    # check using tuple to first argument
    fig.set_size_inches((1, 3))
    assert fig.get_figwidth() == 1
    assert fig.get_figheight() == 3


def test_axes_remove():
    fig, axes = plt.subplots(2, 2)
    axes[-1, -1].remove()
    for ax in axes.ravel()[:-1]:
        assert ax in fig.axes
    assert axes[-1, -1] not in fig.axes
    assert len(fig.axes) == 3


def test_figaspect():
    w, h = plt.figaspect(np.float64(2) / np.float64(1))
    assert h / w == 2
    w, h = plt.figaspect(2)
    assert h / w == 2
    w, h = plt.figaspect(np.zeros((1, 2)))
    assert h / w == 0.5
    w, h = plt.figaspect(np.zeros((2, 2)))
    assert h / w == 1


@pytest.mark.parametrize('which', [None, 'both', 'major', 'minor'])
def test_autofmt_xdate(which):
    date = ['3 Jan 2013', '4 Jan 2013', '5 Jan 2013', '6 Jan 2013',
            '7 Jan 2013', '8 Jan 2013', '9 Jan 2013', '10 Jan 2013',
            '11 Jan 2013', '12 Jan 2013', '13 Jan 2013', '14 Jan 2013']

    time = ['16:44:00', '16:45:00', '16:46:00', '16:47:00', '16:48:00',
            '16:49:00', '16:51:00', '16:52:00', '16:53:00', '16:55:00',
            '16:56:00', '16:57:00']

    angle = 60
    minors = [1, 2, 3, 4, 5, 6, 7]

    x = mdates.datestr2num(date)
    y = mdates.datestr2num(time)

    fig, ax = plt.subplots()

    ax.plot(x, y)
    ax.yaxis_date()
    ax.xaxis_date()

    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.xaxis.set_minor_formatter(FixedFormatter(minors))

    fig.autofmt_xdate(0.2, angle, 'right', which)

    if which in ('both', 'major', None):
        for label in fig.axes[0].get_xticklabels(False, 'major'):
            assert int(label.get_rotation()) == angle

    if which in ('both', 'minor'):
        for label in fig.axes[0].get_xticklabels(True, 'minor'):
            assert int(label.get_rotation()) == angle


@pytest.mark.style('default')
def test_change_dpi():
    fig = plt.figure(figsize=(4, 4))
    fig.canvas.draw()
    assert fig.canvas.renderer.height == 400
    assert fig.canvas.renderer.width == 400
    fig.dpi = 50
    fig.canvas.draw()
    assert fig.canvas.renderer.height == 200
    assert fig.canvas.renderer.width == 200


def test_invalid_figure_size():
    with pytest.raises(ValueError):
        plt.figure(figsize=(1, np.nan))

    fig = plt.figure()
    with pytest.raises(ValueError):
        fig.set_size_inches(1, np.nan)

    with pytest.raises(ValueError):
        fig.add_axes((.1, .1, .5, np.nan))


def test_subplots_shareax_loglabels():
    fig, ax_arr = plt.subplots(2, 2, sharex=True, sharey=True, squeeze=False)
    for ax in ax_arr.flatten():
        ax.plot([10, 20, 30], [10, 20, 30])

    ax.set_yscale("log")
    ax.set_xscale("log")

    for ax in ax_arr[0, :]:
        assert 0 == len(ax.xaxis.get_ticklabels(which='both'))

    for ax in ax_arr[1, :]:
        assert 0 < len(ax.xaxis.get_ticklabels(which='both'))

    for ax in ax_arr[:, 1]:
        assert 0 == len(ax.yaxis.get_ticklabels(which='both'))

    for ax in ax_arr[:, 0]:
        assert 0 < len(ax.yaxis.get_ticklabels(which='both'))


def test_savefig():
    fig = plt.figure()
    msg = "savefig() takes 2 positional arguments but 3 were given"
    with pytest.raises(TypeError, message=msg):
        fig.savefig("fname1.png", "fname2.png")


def test_figure_repr():
    fig = plt.figure(figsize=(10, 20), dpi=10)
    assert repr(fig) == "<Figure size 100x200 with 0 Axes>"


@pytest.mark.skipif(sys.version_info < (3, 6), reason="requires Python 3.6+")
@pytest.mark.parametrize("fmt", ["png", "pdf", "ps", "eps", "svg"])
def test_fspath(fmt, tmpdir):
    from pathlib import Path
    out = Path(tmpdir, "test.{}".format(fmt))
    plt.savefig(out)
    with out.open("rb") as file:
        # All the supported formats include the format name (case-insensitive)
        # in the first 100 bytes.
        assert fmt.encode("ascii") in file.read(100).lower()
