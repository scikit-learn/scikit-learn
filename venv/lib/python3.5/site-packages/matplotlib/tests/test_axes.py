from __future__ import absolute_import, division, print_function

import six
from six.moves import xrange
from itertools import chain, product
from distutils.version import LooseVersion
import io

import datetime

import pytz

import numpy as np
from numpy import ma
from cycler import cycler
import pytest

import warnings

import matplotlib
from matplotlib.testing.decorators import image_comparison
import matplotlib.pyplot as plt
import matplotlib.markers as mmarkers
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from numpy.testing import assert_allclose, assert_array_equal
from matplotlib.cbook import (
    IgnoredKeywordWarning, MatplotlibDeprecationWarning)
from matplotlib.cbook._backports import broadcast_to

# Note: Some test cases are run twice: once normally and once with labeled data
#       These two must be defined in the same test function or need to have
#       different baseline images to prevent race conditions when pytest runs
#       the tests with multiple threads.


def test_get_labels():
    fig, ax = plt.subplots()
    ax.set_xlabel('x label')
    ax.set_ylabel('y label')
    assert ax.get_xlabel() == 'x label'
    assert ax.get_ylabel() == 'y label'


@image_comparison(baseline_images=['acorr'], extensions=['png'], style='mpl20')
def test_acorr():
    np.random.seed(19680801)
    n = 512
    x = np.random.normal(0, 1, n).cumsum()

    fig, ax = plt.subplots()
    ax.acorr(x, maxlags=n - 1, label='acorr')
    ax.legend()


@image_comparison(baseline_images=['spy'], extensions=['png'], style='mpl20')
def test_spy():
    np.random.seed(19680801)
    a = np.ones(32 * 32)
    a[:16 * 32] = 0
    np.random.shuffle(a)
    a = np.reshape(a, (32, 32))

    fig, ax = plt.subplots()
    ax.spy(a)


@image_comparison(baseline_images=['matshow'],
                  extensions=['png'], style='mpl20')
def test_matshow():
    np.random.seed(19680801)
    a = np.random.rand(32, 32)

    fig, ax = plt.subplots()
    ax.matshow(a)


@image_comparison(baseline_images=['formatter_ticker_001',
                                   'formatter_ticker_002',
                                   'formatter_ticker_003',
                                   'formatter_ticker_004',
                                   'formatter_ticker_005',
                                   ])
def test_formatter_ticker():
    import matplotlib.testing.jpl_units as units
    units.register()

    # This should affect the tick size.  (Tests issue #543)
    matplotlib.rcParams['lines.markeredgewidth'] = 30

    # This essentially test to see if user specified labels get overwritten
    # by the auto labeler functionality of the axes.
    xdata = [x*units.sec for x in range(10)]
    ydata1 = [(1.5*y - 0.5)*units.km for y in range(10)]
    ydata2 = [(1.75*y - 1.0)*units.km for y in range(10)]

    fig = plt.figure()
    ax = plt.subplot(111)
    ax.set_xlabel("x-label 001")

    fig = plt.figure()
    ax = plt.subplot(111)
    ax.set_xlabel("x-label 001")
    ax.plot(xdata, ydata1, color='blue', xunits="sec")

    fig = plt.figure()
    ax = plt.subplot(111)
    ax.set_xlabel("x-label 001")
    ax.plot(xdata, ydata1, color='blue', xunits="sec")
    ax.set_xlabel("x-label 003")

    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(xdata, ydata1, color='blue', xunits="sec")
    ax.plot(xdata, ydata2, color='green', xunits="hour")
    ax.set_xlabel("x-label 004")

    # See SF bug 2846058
    # https://sourceforge.net/tracker/?func=detail&aid=2846058&group_id=80706&atid=560720
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(xdata, ydata1, color='blue', xunits="sec")
    ax.plot(xdata, ydata2, color='green', xunits="hour")
    ax.set_xlabel("x-label 005")
    ax.autoscale_view()


@image_comparison(baseline_images=["formatter_large_small"])
def test_formatter_large_small():
    # github issue #617, pull #619
    if LooseVersion(np.__version__) >= LooseVersion('1.11.0'):
        pytest.skip("Fall out from a fixed numpy bug")
    fig, ax = plt.subplots(1)
    x = [0.500000001, 0.500000002]
    y = [1e64, 1.1e64]
    ax.plot(x, y)


@image_comparison(baseline_images=["twin_axis_locaters_formatters"])
def test_twin_axis_locaters_formatters():
    vals = np.linspace(0, 1, num=5, endpoint=True)
    locs = np.sin(np.pi * vals / 2.0)

    majl = plt.FixedLocator(locs)
    minl = plt.FixedLocator([0.1, 0.2, 0.3])

    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.plot([0.1, 100], [0, 1])
    ax1.yaxis.set_major_locator(majl)
    ax1.yaxis.set_minor_locator(minl)
    ax1.yaxis.set_major_formatter(plt.FormatStrFormatter('%08.2lf'))
    ax1.yaxis.set_minor_formatter(plt.FixedFormatter(['tricks', 'mind',
                                                      'jedi']))

    ax1.xaxis.set_major_locator(plt.LinearLocator())
    ax1.xaxis.set_minor_locator(plt.FixedLocator([15, 35, 55, 75]))
    ax1.xaxis.set_major_formatter(plt.FormatStrFormatter('%05.2lf'))
    ax1.xaxis.set_minor_formatter(plt.FixedFormatter(['c', '3', 'p', 'o']))
    ax2 = ax1.twiny()
    ax3 = ax1.twinx()


def test_twinx_cla():
    fig, ax = plt.subplots()
    ax2 = ax.twinx()
    ax3 = ax2.twiny()
    plt.draw()
    assert not ax2.xaxis.get_visible()
    assert not ax2.patch.get_visible()
    ax2.cla()
    ax3.cla()

    assert not ax2.xaxis.get_visible()
    assert not ax2.patch.get_visible()
    assert ax2.yaxis.get_visible()

    assert ax3.xaxis.get_visible()
    assert not ax3.patch.get_visible()
    assert not ax3.yaxis.get_visible()

    assert ax.xaxis.get_visible()
    assert ax.patch.get_visible()
    assert ax.yaxis.get_visible()


@image_comparison(baseline_images=['twin_autoscale'], extensions=['png'])
def test_twinx_axis_scales():
    x = np.array([0, 0.5, 1])
    y = 0.5 * x
    x2 = np.array([0, 1, 2])
    y2 = 2 * x2

    fig = plt.figure()
    ax = fig.add_axes((0, 0, 1, 1), autoscalex_on=False, autoscaley_on=False)
    ax.plot(x, y, color='blue', lw=10)

    ax2 = plt.twinx(ax)
    ax2.plot(x2, y2, 'r--', lw=5)

    ax.margins(0, 0)
    ax2.margins(0, 0)


def test_twin_inherit_autoscale_setting():
    fig, ax = plt.subplots()
    ax_x_on = ax.twinx()
    ax.set_autoscalex_on(False)
    ax_x_off = ax.twinx()

    assert ax_x_on.get_autoscalex_on()
    assert not ax_x_off.get_autoscalex_on()

    ax_y_on = ax.twiny()
    ax.set_autoscaley_on(False)
    ax_y_off = ax.twiny()

    assert ax_y_on.get_autoscaley_on()
    assert not ax_y_off.get_autoscaley_on()


def test_inverted_cla():
    # Github PR #5450. Setting autoscale should reset
    # axes to be non-inverted.
    # plotting an image, then 1d graph, axis is now down
    fig = plt.figure(0)
    ax = fig.gca()
    # 1. test that a new axis is not inverted per default
    assert not ax.xaxis_inverted()
    assert not ax.yaxis_inverted()
    img = np.random.random((100, 100))
    ax.imshow(img)
    # 2. test that a image axis is inverted
    assert not ax.xaxis_inverted()
    assert ax.yaxis_inverted()
    # 3. test that clearing and plotting a line, axes are
    # not inverted
    ax.cla()
    x = np.linspace(0, 2*np.pi, 100)
    ax.plot(x, np.cos(x))
    assert not ax.xaxis_inverted()
    assert not ax.yaxis_inverted()

    # 4. autoscaling should not bring back axes to normal
    ax.cla()
    ax.imshow(img)
    plt.autoscale()
    assert not(ax.xaxis_inverted())
    assert ax.yaxis_inverted()

    # 5. two shared axes. Clearing the master axis should bring axes in shared
    # axes back to normal
    ax0 = plt.subplot(211)
    ax1 = plt.subplot(212, sharey=ax0)
    ax0.imshow(img)
    ax1.plot(x, np.cos(x))
    ax0.cla()
    assert not(ax1.yaxis_inverted())
    ax1.cla()
    # 6. clearing the nonmaster should not touch limits
    ax0.imshow(img)
    ax1.plot(x, np.cos(x))
    ax1.cla()
    assert ax.yaxis_inverted()

    # clean up
    plt.close(fig)


@image_comparison(baseline_images=["minorticks_on_rcParams_both"],
                  extensions=['png'])
def test_minorticks_on_rcParams_both():
    fig = plt.figure()
    matplotlib.rcParams['xtick.minor.visible'] = True
    matplotlib.rcParams['ytick.minor.visible'] = True

    plt.plot([0, 1], [0, 1])
    plt.axis([0, 1, 0, 1])


@image_comparison(baseline_images=["autoscale_tiny_range"], remove_text=True)
def test_autoscale_tiny_range():
    # github pull #904
    fig, ax = plt.subplots(2, 2)
    ax = ax.flatten()
    for i in xrange(4):
        y1 = 10**(-11 - i)
        ax[i].plot([0, 1], [1, 1 + y1])


@pytest.mark.style('default')
def test_autoscale_tight():
    fig, ax = plt.subplots(1, 1)
    ax.plot([1, 2, 3, 4])
    ax.autoscale(enable=True, axis='x', tight=False)
    ax.autoscale(enable=True, axis='y', tight=True)
    assert_allclose(ax.get_xlim(), (-0.15, 3.15))
    assert_allclose(ax.get_ylim(), (1.0, 4.0))


@pytest.mark.style('default')
def test_autoscale_log_shared():
    # related to github #7587
    # array starts at zero to trigger _minpos handling
    x = np.arange(100, dtype=float)
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    ax1.loglog(x, x)
    ax2.semilogx(x, x)
    ax1.autoscale(tight=True)
    ax2.autoscale(tight=True)
    plt.draw()
    lims = (x[1], x[-1])
    assert_allclose(ax1.get_xlim(), lims)
    assert_allclose(ax1.get_ylim(), lims)
    assert_allclose(ax2.get_xlim(), lims)
    assert_allclose(ax2.get_ylim(), (x[0], x[-1]))


@pytest.mark.style('default')
def test_use_sticky_edges():
    fig, ax = plt.subplots()
    ax.imshow([[0, 1], [2, 3]], origin='lower')
    assert_allclose(ax.get_xlim(), (-0.5, 1.5))
    assert_allclose(ax.get_ylim(), (-0.5, 1.5))
    ax.use_sticky_edges = False
    ax.autoscale()
    xlim = (-0.5 - 2 * ax._xmargin, 1.5 + 2 * ax._xmargin)
    ylim = (-0.5 - 2 * ax._ymargin, 1.5 + 2 * ax._ymargin)
    assert_allclose(ax.get_xlim(), xlim)
    assert_allclose(ax.get_ylim(), ylim)
    # Make sure it is reversible:
    ax.use_sticky_edges = True
    ax.autoscale()
    assert_allclose(ax.get_xlim(), (-0.5, 1.5))
    assert_allclose(ax.get_ylim(), (-0.5, 1.5))


@image_comparison(baseline_images=['offset_points'],
                  remove_text=True)
def test_basic_annotate():
    # Setup some data
    t = np.arange(0.0, 5.0, 0.01)
    s = np.cos(2.0*np.pi * t)

    # Offset Points

    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(-1, 5), ylim=(-3, 5))
    line, = ax.plot(t, s, lw=3, color='purple')

    ax.annotate('local max', xy=(3, 1), xycoords='data',
                xytext=(3, 3), textcoords='offset points')


@image_comparison(baseline_images=['arrow_simple'],
                  extensions=['png'], remove_text=True)
def test_arrow_simple():
    # Simple image test for ax.arrow
    # kwargs that take discrete values
    length_includes_head = (True, False)
    shape = ('full', 'left', 'right')
    head_starts_at_zero = (True, False)
    # Create outer product of values
    kwargs = list(product(length_includes_head, shape, head_starts_at_zero))

    fig, axs = plt.subplots(3, 4)
    for i, (ax, kwarg) in enumerate(zip(axs.flatten(), kwargs)):
        ax.set_xlim(-2, 2)
        ax.set_ylim(-2, 2)
        # Unpack kwargs
        (length_includes_head, shape, head_starts_at_zero) = kwarg
        theta = 2 * np.pi * i / 12
        # Draw arrow
        ax.arrow(0, 0, np.sin(theta), np.cos(theta),
                 width=theta/100,
                 length_includes_head=length_includes_head,
                 shape=shape,
                 head_starts_at_zero=head_starts_at_zero,
                 head_width=theta / 10,
                 head_length=theta / 10)


def test_annotate_default_arrow():
    # Check that we can make an annotation arrow with only default properties.
    fig, ax = plt.subplots()
    ann = ax.annotate("foo", (0, 1), xytext=(2, 3))
    assert ann.arrow_patch is None
    ann = ax.annotate("foo", (0, 1), xytext=(2, 3), arrowprops={})
    assert ann.arrow_patch is not None


@image_comparison(baseline_images=['polar_axes'], style='default')
def test_polar_annotations():
    # you can specify the xypoint and the xytext in different
    # positions and coordinate systems, and optionally turn on a
    # connecting line and mark the point with a marker.  Annotations
    # work on polar axes too.  In the example below, the xy point is
    # in native coordinates (xycoords defaults to 'data').  For a
    # polar axes, this is in (theta, radius) space.  The text in this
    # example is placed in the fractional figure coordinate system.
    # Text keyword args like horizontal and vertical alignment are
    # respected

    # Setup some data
    r = np.arange(0.0, 1.0, 0.001)
    theta = 2.0 * 2.0 * np.pi * r

    fig = plt.figure()
    ax = fig.add_subplot(111, polar=True)
    line, = ax.plot(theta, r, color='#ee8d18', lw=3)
    line, = ax.plot((0, 0), (0, 1), color="#0000ff", lw=1)

    ind = 800
    thisr, thistheta = r[ind], theta[ind]
    ax.plot([thistheta], [thisr], 'o')
    ax.annotate('a polar annotation',
                xy=(thistheta, thisr),  # theta, radius
                xytext=(0.05, 0.05),    # fraction, fraction
                textcoords='figure fraction',
                arrowprops=dict(facecolor='black', shrink=0.05),
                horizontalalignment='left',
                verticalalignment='baseline',
                )

    ax.tick_params(axis='x', tick1On=True, tick2On=True, direction='out')


@image_comparison(baseline_images=['polar_coords'], style='default',
                  remove_text=True)
def test_polar_coord_annotations():
    # You can also use polar notation on a catesian axes.  Here the
    # native coordinate system ('data') is cartesian, so you need to
    # specify the xycoords and textcoords as 'polar' if you want to
    # use (theta, radius)
    from matplotlib.patches import Ellipse
    el = Ellipse((0, 0), 10, 20, facecolor='r', alpha=0.5)

    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')

    ax.add_artist(el)
    el.set_clip_box(ax.bbox)

    ax.annotate('the top',
                xy=(np.pi/2., 10.),      # theta, radius
                xytext=(np.pi/3, 20.),   # theta, radius
                xycoords='polar',
                textcoords='polar',
                arrowprops=dict(facecolor='black', shrink=0.05),
                horizontalalignment='left',
                verticalalignment='baseline',
                clip_on=True,  # clip to the axes bounding box
                )

    ax.set_xlim(-20, 20)
    ax.set_ylim(-20, 20)


@image_comparison(baseline_images=['fill_units'], extensions=['png'],
                  savefig_kwarg={'dpi': 60})
def test_fill_units():
    from datetime import datetime
    import matplotlib.testing.jpl_units as units
    units.register()

    # generate some data
    t = units.Epoch("ET", dt=datetime(2009, 4, 27))
    value = 10.0 * units.deg
    day = units.Duration("ET", 24.0 * 60.0 * 60.0)

    fig = plt.figure()

    # Top-Left
    ax1 = fig.add_subplot(221)
    ax1.plot([t], [value], yunits='deg', color='red')
    ax1.fill([733525.0, 733525.0, 733526.0, 733526.0],
             [0.0, 0.0, 90.0, 0.0], 'b')

    # Top-Right
    ax2 = fig.add_subplot(222)
    ax2.plot([t], [value], yunits='deg', color='red')
    ax2.fill([t, t, t + day, t + day],
             [0.0, 0.0, 90.0, 0.0], 'b')

    # Bottom-Left
    ax3 = fig.add_subplot(223)
    ax3.plot([t], [value], yunits='deg', color='red')
    ax3.fill([733525.0, 733525.0, 733526.0, 733526.0],
             [0 * units.deg, 0 * units.deg, 90 * units.deg, 0 * units.deg],
             'b')

    # Bottom-Right
    ax4 = fig.add_subplot(224)
    ax4.plot([t], [value], yunits='deg', color='red')
    ax4.fill([t, t, t + day, t + day],
             [0 * units.deg, 0 * units.deg, 90 * units.deg, 0 * units.deg],
             facecolor="blue")

    fig.autofmt_xdate()


@image_comparison(baseline_images=['single_point', 'single_point'])
def test_single_point():
    # Issue #1796: don't let lines.marker affect the grid
    matplotlib.rcParams['lines.marker'] = 'o'
    matplotlib.rcParams['axes.grid'] = True

    fig = plt.figure()
    plt.subplot(211)
    plt.plot([0], [0], 'o')

    plt.subplot(212)
    plt.plot([1], [1], 'o')

    # Reuse testcase from above for a labeled data test
    data = {'a': [0], 'b': [1]}

    fig = plt.figure()
    plt.subplot(211)
    plt.plot('a', 'a', 'o', data=data)

    plt.subplot(212)
    plt.plot('b', 'b', 'o', data=data)


@image_comparison(baseline_images=['single_date'])
def test_single_date():
    time1 = [721964.0]
    data1 = [-65.54]

    fig = plt.figure()
    plt.subplot(211)
    plt.plot_date(time1, data1, 'o', color='r')

    plt.subplot(212)
    plt.plot(time1, data1, 'o', color='r')


@image_comparison(baseline_images=['shaped_data'])
def test_shaped_data():
    xdata = np.array([[0.53295185, 0.23052951, 0.19057629, 0.66724975,
                       0.96577916, 0.73136095, 0.60823287, 0.01792100,
                       0.29744742, 0.27164665],
                      [0.27980120, 0.25814229, 0.02818193, 0.12966456,
                       0.57446277, 0.58167607, 0.71028245, 0.69112737,
                       0.89923072, 0.99072476],
                      [0.81218578, 0.80464528, 0.76071809, 0.85616314,
                       0.12757994, 0.94324936, 0.73078663, 0.09658102,
                       0.60703967, 0.77664978],
                      [0.28332265, 0.81479711, 0.86985333, 0.43797066,
                       0.32540082, 0.43819229, 0.92230363, 0.49414252,
                       0.68168256, 0.05922372],
                      [0.10721335, 0.93904142, 0.79163075, 0.73232848,
                       0.90283839, 0.68408046, 0.25502302, 0.95976614,
                       0.59214115, 0.13663711],
                      [0.28087456, 0.33127607, 0.15530412, 0.76558121,
                       0.83389773, 0.03735974, 0.98717738, 0.71432229,
                       0.54881366, 0.86893953],
                      [0.77995937, 0.99555600, 0.29688434, 0.15646162,
                       0.05184800, 0.37161935, 0.12998491, 0.09377296,
                       0.36882507, 0.36583435],
                      [0.37851836, 0.05315792, 0.63144617, 0.25003433,
                       0.69586032, 0.11393988, 0.92362096, 0.88045438,
                       0.93530252, 0.68275072],
                      [0.86486596, 0.83236675, 0.82960664, 0.57796630,
                       0.25724233, 0.84841095, 0.90862812, 0.64414887,
                       0.35652720, 0.71026066],
                      [0.01383268, 0.34060930, 0.76084285, 0.70800694,
                       0.87634056, 0.08213693, 0.54655021, 0.98123181,
                       0.44080053, 0.86815815]])

    y1 = np.arange(10).reshape((1, -1))
    y2 = np.arange(10).reshape((-1, 1))

    fig = plt.figure()
    plt.subplot(411)
    plt.plot(y1)
    plt.subplot(412)
    plt.plot(y2)

    plt.subplot(413)
    with pytest.raises(ValueError):
        plt.plot((y1, y2))

    plt.subplot(414)
    plt.plot(xdata[:, 1], xdata[1, :], 'o')


@image_comparison(baseline_images=['const_xy'])
def test_const_xy():
    fig = plt.figure()

    plt.subplot(311)
    plt.plot(np.arange(10), np.ones((10,)))

    plt.subplot(312)
    plt.plot(np.ones((10,)), np.arange(10))

    plt.subplot(313)
    plt.plot(np.ones((10,)), np.ones((10,)), 'o')


@image_comparison(baseline_images=['polar_wrap_180', 'polar_wrap_360'],
                  style='default')
def test_polar_wrap():
    fig = plt.figure()
    plt.subplot(111, polar=True)
    plt.polar(np.deg2rad([179, -179]), [0.2, 0.1], "b.-")
    plt.polar(np.deg2rad([179,  181]), [0.2, 0.1], "g.-")
    plt.rgrids([0.05, 0.1, 0.15, 0.2, 0.25, 0.3])
    assert len(fig.axes) == 1, 'More than one polar axes created.'

    fig = plt.figure()
    plt.subplot(111, polar=True)
    plt.polar(np.deg2rad([2, -2]), [0.2, 0.1], "b.-")
    plt.polar(np.deg2rad([2, 358]), [0.2, 0.1], "g.-")
    plt.polar(np.deg2rad([358, 2]), [0.2, 0.1], "r.-")
    plt.rgrids([0.05, 0.1, 0.15, 0.2, 0.25, 0.3])


@image_comparison(baseline_images=['polar_units', 'polar_units_2'],
                  style='default')
def test_polar_units():
    import matplotlib.testing.jpl_units as units
    units.register()

    pi = np.pi
    deg = units.deg
    km = units.km

    x1 = [pi/6.0, pi/4.0, pi/3.0, pi/2.0]
    x2 = [30.0*deg, 45.0*deg, 60.0*deg, 90.0*deg]

    y1 = [1.0, 2.0, 3.0, 4.0]
    y2 = [4.0, 3.0, 2.0, 1.0]

    fig = plt.figure()

    plt.polar(x2, y1, color="blue")

    # polar(x2, y1, color = "red", xunits="rad")
    # polar(x2, y2, color = "green")

    fig = plt.figure()

    # make sure runits and theta units work
    y1 = [y*km for y in y1]
    plt.polar(x2, y1, color="blue", thetaunits="rad", runits="km")
    assert isinstance(plt.gca().get_xaxis().get_major_formatter(),
                      units.UnitDblFormatter)


@image_comparison(baseline_images=['polar_rmin'], style='default')
def test_polar_rmin():
    r = np.arange(0, 3.0, 0.01)
    theta = 2*np.pi*r

    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
    ax.plot(theta, r)
    ax.set_rmax(2.0)
    ax.set_rmin(0.5)


@image_comparison(baseline_images=['polar_negative_rmin'], style='default')
def test_polar_negative_rmin():
    r = np.arange(-3.0, 0.0, 0.01)
    theta = 2*np.pi*r

    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
    ax.plot(theta, r)
    ax.set_rmax(0.0)
    ax.set_rmin(-3.0)


@image_comparison(baseline_images=['polar_rorigin'], style='default')
def test_polar_rorigin():
    r = np.arange(0, 3.0, 0.01)
    theta = 2*np.pi*r

    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
    ax.plot(theta, r)
    ax.set_rmax(2.0)
    ax.set_rmin(0.5)
    ax.set_rorigin(0.0)


@image_comparison(baseline_images=['polar_theta_position'], style='default')
def test_polar_theta_position():
    r = np.arange(0, 3.0, 0.01)
    theta = 2*np.pi*r

    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
    ax.plot(theta, r)
    ax.set_theta_zero_location("NW", 30)
    ax.set_theta_direction('clockwise')


@image_comparison(baseline_images=['polar_rlabel_position'], style='default')
def test_polar_rlabel_position():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='polar')
    ax.set_rlabel_position(315)
    ax.tick_params(rotation='auto')


@image_comparison(baseline_images=['polar_theta_wedge'], style='default',
                  tol=0.01 if six.PY2 else 0)
def test_polar_theta_limits():
    r = np.arange(0, 3.0, 0.01)
    theta = 2*np.pi*r

    theta_mins = np.arange(15.0, 361.0, 90.0)
    theta_maxs = np.arange(50.0, 361.0, 90.0)
    DIRECTIONS = ('out', 'in', 'inout')

    fig, axes = plt.subplots(len(theta_mins), len(theta_maxs),
                             subplot_kw={'polar': True},
                             figsize=(8, 6))

    for i, start in enumerate(theta_mins):
        for j, end in enumerate(theta_maxs):
            ax = axes[i, j]
            ax.plot(theta, r)
            if start < end:
                ax.set_thetamin(start)
                ax.set_thetamax(end)
            else:
                # Plot with clockwise orientation instead.
                ax.set_thetamin(end)
                ax.set_thetamax(start)
                ax.set_theta_direction('clockwise')
            ax.tick_params(tick1On=True, tick2On=True,
                           direction=DIRECTIONS[i % len(DIRECTIONS)],
                           rotation='auto')
            ax.yaxis.set_tick_params(label2On=True, rotation='auto')


@image_comparison(baseline_images=['axvspan_epoch'])
def test_axvspan_epoch():
    from datetime import datetime
    import matplotlib.testing.jpl_units as units
    units.register()

    # generate some data
    t0 = units.Epoch("ET", dt=datetime(2009, 1, 20))
    tf = units.Epoch("ET", dt=datetime(2009, 1, 21))

    dt = units.Duration("ET", units.day.convert("sec"))

    fig = plt.figure()

    plt.axvspan(t0, tf, facecolor="blue", alpha=0.25)

    ax = plt.gca()
    ax.set_xlim(t0 - 5.0*dt, tf + 5.0*dt)


@image_comparison(baseline_images=['axhspan_epoch'])
def test_axhspan_epoch():
    from datetime import datetime
    import matplotlib.testing.jpl_units as units
    units.register()

    # generate some data
    t0 = units.Epoch("ET", dt=datetime(2009, 1, 20))
    tf = units.Epoch("ET", dt=datetime(2009, 1, 21))

    dt = units.Duration("ET", units.day.convert("sec"))

    fig = plt.figure()

    plt.axhspan(t0, tf, facecolor="blue", alpha=0.25)

    ax = plt.gca()
    ax.set_ylim(t0 - 5.0*dt, tf + 5.0*dt)


@image_comparison(baseline_images=['hexbin_extent', 'hexbin_extent'],
                  remove_text=True, extensions=['png'])
def test_hexbin_extent():
    # this test exposes sf bug 2856228
    fig = plt.figure()

    ax = fig.add_subplot(111)
    data = (np.arange(2000) / 2000).reshape((2, 1000))
    x, y = data

    ax.hexbin(x, y, extent=[.1, .3, .6, .7])

    # Reuse testcase from above for a labeled data test
    data = {"x": x, "y": y}

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hexbin("x", "y", extent=[.1, .3, .6, .7], data=data)


@image_comparison(baseline_images=['hexbin_empty'], remove_text=True,
                  extensions=['png'])
def test_hexbin_empty():
    # From #3886: creating hexbin from empty dataset raises ValueError
    ax = plt.gca()
    ax.hexbin([], [])


def test_hexbin_pickable():
    # From #1973: Test that picking a hexbin collection works
    class FauxMouseEvent:
        def __init__(self, x, y):
            self.x = x
            self.y = y

    fig = plt.figure()

    ax = fig.add_subplot(111)
    data = (np.arange(200) / 200).reshape((2, 100))
    x, y = data
    hb = ax.hexbin(x, y, extent=[.1, .3, .6, .7], picker=-1)

    assert hb.contains(FauxMouseEvent(400, 300))[0]


@image_comparison(baseline_images=['hexbin_log'],
                  remove_text=True,
                  extensions=['png'])
def test_hexbin_log():
    # Issue #1636
    fig = plt.figure()

    np.random.seed(0)
    n = 100000
    x = np.random.standard_normal(n)
    y = 2.0 + 3.0 * x + 4.0 * np.random.standard_normal(n)
    y = np.power(2, y * 0.5)
    ax = fig.add_subplot(111)
    ax.hexbin(x, y, yscale='log')


def test_inverted_limits():
    # Test gh:1553
    # Calling invert_xaxis prior to plotting should not disable autoscaling
    # while still maintaining the inverted direction
    fig = plt.figure()
    ax = fig.gca()
    ax.invert_xaxis()
    ax.plot([-5, -3, 2, 4], [1, 2, -3, 5])

    assert ax.get_xlim() == (4, -5)
    assert ax.get_ylim() == (-3, 5)
    plt.close()

    fig = plt.figure()
    ax = fig.gca()
    ax.invert_yaxis()
    ax.plot([-5, -3, 2, 4], [1, 2, -3, 5])

    assert ax.get_xlim() == (-5, 4)
    assert ax.get_ylim() == (5, -3)
    plt.close()


@image_comparison(baseline_images=['nonfinite_limits'])
def test_nonfinite_limits():
    x = np.arange(0., np.e, 0.01)
    # silence divide by zero warning from log(0)
    olderr = np.seterr(divide='ignore')
    try:
        y = np.log(x)
    finally:
        np.seterr(**olderr)
    x[len(x)//2] = np.nan
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y)


@image_comparison(baseline_images=['imshow', 'imshow'],
                  remove_text=True, style='mpl20')
def test_imshow():
    # Create a NxN image
    N = 100
    (x, y) = np.indices((N, N))
    x -= N//2
    y -= N//2
    r = np.sqrt(x**2+y**2-x*y)

    # Create a contour plot at N/4 and extract both the clip path and transform
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.imshow(r)

    # Reuse testcase from above for a labeled data test
    data = {"r": r}
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow("r", data=data)


@image_comparison(baseline_images=['imshow_clip'], style='mpl20')
def test_imshow_clip():
    # As originally reported by Gellule Xg <gellule.xg@free.fr>

    # Create a NxN image
    N = 100
    (x, y) = np.indices((N, N))
    x -= N//2
    y -= N//2
    r = np.sqrt(x**2+y**2-x*y)

    # Create a contour plot at N/4 and extract both the clip path and transform
    fig = plt.figure()
    ax = fig.add_subplot(111)

    c = ax.contour(r, [N/4])
    x = c.collections[0]
    clipPath = x.get_paths()[0]
    clipTransform = x.get_transform()

    from matplotlib.transforms import TransformedPath
    clip_path = TransformedPath(clipPath, clipTransform)

    # Plot the image clipped by the contour
    ax.imshow(r, clip_path=clip_path)


@image_comparison(baseline_images=['polycollection_joinstyle'],
                  remove_text=True)
def test_polycollection_joinstyle():
    # Bug #2890979 reported by Matthew West

    from matplotlib import collections as mcoll

    fig = plt.figure()
    ax = fig.add_subplot(111)
    verts = np.array([[1, 1], [1, 2], [2, 2], [2, 1]])
    c = mcoll.PolyCollection([verts], linewidths=40)
    ax.add_collection(c)
    ax.set_xbound(0, 3)
    ax.set_ybound(0, 3)


@pytest.mark.parametrize(
    'x, y1, y2', [
        (np.zeros((2, 2)), 3, 3),
        (np.arange(0.0, 2, 0.02), np.zeros((2, 2)), 3),
        (np.arange(0.0, 2, 0.02), 3, np.zeros((2, 2)))
    ], ids=[
        '2d_x_input',
        '2d_y1_input',
        '2d_y2_input'
    ]
)
def test_fill_between_input(x, y1, y2):
    fig = plt.figure()
    ax = fig.add_subplot(211)
    with pytest.raises(ValueError):
        ax.fill_between(x, y1, y2)


@pytest.mark.parametrize(
    'y, x1, x2', [
        (np.zeros((2, 2)), 3, 3),
        (np.arange(0.0, 2, 0.02), np.zeros((2, 2)), 3),
        (np.arange(0.0, 2, 0.02), 3, np.zeros((2, 2)))
    ], ids=[
        '2d_y_input',
        '2d_x1_input',
        '2d_x2_input'
    ]
)
def test_fill_betweenx_input(y, x1, x2):
    fig = plt.figure()
    ax = fig.add_subplot(211)
    with pytest.raises(ValueError):
        ax.fill_betweenx(y, x1, x2)


@image_comparison(baseline_images=['fill_between_interpolate'],
                  remove_text=True)
def test_fill_between_interpolate():
    x = np.arange(0.0, 2, 0.02)
    y1 = np.sin(2*np.pi*x)
    y2 = 1.2*np.sin(4*np.pi*x)

    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax.plot(x, y1, x, y2, color='black')
    ax.fill_between(x, y1, y2, where=y2 >= y1, facecolor='white', hatch='/',
                    interpolate=True)
    ax.fill_between(x, y1, y2, where=y2 <= y1, facecolor='red',
                    interpolate=True)

    # Test support for masked arrays.
    y2 = np.ma.masked_greater(y2, 1.0)
    # Test that plotting works for masked arrays with the first element masked
    y2[0] = np.ma.masked
    ax1 = fig.add_subplot(212, sharex=ax)
    ax1.plot(x, y1, x, y2, color='black')
    ax1.fill_between(x, y1, y2, where=y2 >= y1, facecolor='green',
                     interpolate=True)
    ax1.fill_between(x, y1, y2, where=y2 <= y1, facecolor='red',
                     interpolate=True)


@image_comparison(baseline_images=['fill_between_interpolate_decreasing'],
                  style='mpl20', remove_text=True)
def test_fill_between_interpolate_decreasing():
    p = np.array([724.3, 700, 655])
    t = np.array([9.4, 7, 2.2])
    prof = np.array([7.9, 6.6, 3.8])

    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_subplot(1, 1, 1)

    ax.plot(t, p, 'tab:red')
    ax.plot(prof, p, 'k')

    ax.fill_betweenx(p, t, prof, where=prof < t,
                     facecolor='blue', interpolate=True, alpha=0.4)
    ax.fill_betweenx(p, t, prof, where=prof > t,
                     facecolor='red', interpolate=True, alpha=0.4)

    ax.set_xlim(0, 30)
    ax.set_ylim(800, 600)


@image_comparison(baseline_images=['symlog'])
def test_symlog():
    x = np.array([0, 1, 2, 4, 6, 9, 12, 24])
    y = np.array([1000000, 500000, 100000, 100, 5, 0, 0, 0])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y)
    ax.set_yscale('symlog')
    ax.set_xscale('linear')
    ax.set_ylim(-1, 10000000)


@image_comparison(baseline_images=['symlog2'],
                  remove_text=True)
def test_symlog2():
    # Numbers from -50 to 50, with 0.1 as step
    x = np.arange(-50, 50, 0.001)

    fig = plt.figure()
    ax = fig.add_subplot(511)
    # Plots a simple linear function 'f(x) = x'
    ax.plot(x, x)
    ax.set_xscale('symlog', linthreshx=20.0)
    ax.grid(True)

    ax = fig.add_subplot(512)
    # Plots a simple linear function 'f(x) = x'
    ax.plot(x, x)
    ax.set_xscale('symlog', linthreshx=2.0)
    ax.grid(True)

    ax = fig.add_subplot(513)
    # Plots a simple linear function 'f(x) = x'
    ax.plot(x, x)
    ax.set_xscale('symlog', linthreshx=1.0)
    ax.grid(True)

    ax = fig.add_subplot(514)
    # Plots a simple linear function 'f(x) = x'
    ax.plot(x, x)
    ax.set_xscale('symlog', linthreshx=0.1)
    ax.grid(True)

    ax = fig.add_subplot(515)
    # Plots a simple linear function 'f(x) = x'
    ax.plot(x, x)
    ax.set_xscale('symlog', linthreshx=0.01)
    ax.grid(True)
    ax.set_ylim(-0.1, 0.1)


def test_pcolorargs_5205():
    # Smoketest to catch issue found in gh:5205
    x = [-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5]
    y = [-1.5, -1.25, -1.0, -0.75, -0.5, -0.25, 0,
         0.25, 0.5, 0.75, 1.0, 1.25, 1.5]
    X, Y = np.meshgrid(x, y)
    Z = np.hypot(X, Y)

    plt.pcolor(Z)
    plt.pcolor(list(Z))
    plt.pcolor(x, y, Z)
    plt.pcolor(X, Y, list(Z))


@image_comparison(baseline_images=['pcolormesh'], remove_text=True)
def test_pcolormesh():
    n = 12
    x = np.linspace(-1.5, 1.5, n)
    y = np.linspace(-1.5, 1.5, n*2)
    X, Y = np.meshgrid(x, y)
    Qx = np.cos(Y) - np.cos(X)
    Qz = np.sin(Y) + np.sin(X)
    Qx = (Qx + 1.1)
    Z = np.hypot(X, Y) / 5
    Z = (Z - Z.min()) / Z.ptp()

    # The color array can include masked values:
    Zm = ma.masked_where(np.abs(Qz) < 0.5 * np.max(Qz), Z)

    fig = plt.figure()
    ax = fig.add_subplot(131)
    ax.pcolormesh(Qx, Qz, Z, lw=0.5, edgecolors='k')

    ax = fig.add_subplot(132)
    ax.pcolormesh(Qx, Qz, Z, lw=2, edgecolors=['b', 'w'])

    ax = fig.add_subplot(133)
    ax.pcolormesh(Qx, Qz, Z, shading="gouraud")


@image_comparison(baseline_images=['pcolormesh_datetime_axis'],
                  extensions=['png'], remove_text=False)
def test_pcolormesh_datetime_axis():
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.4, top=0.98, bottom=.15)
    base = datetime.datetime(2013, 1, 1)
    x = np.array([base + datetime.timedelta(days=d) for d in range(21)])
    y = np.arange(21)
    z1, z2 = np.meshgrid(np.arange(20), np.arange(20))
    z = z1 * z2
    plt.subplot(221)
    plt.pcolormesh(x[:-1], y[:-1], z)
    plt.subplot(222)
    plt.pcolormesh(x, y, z)
    x = np.repeat(x[np.newaxis], 21, axis=0)
    y = np.repeat(y[:, np.newaxis], 21, axis=1)
    plt.subplot(223)
    plt.pcolormesh(x[:-1, :-1], y[:-1, :-1], z)
    plt.subplot(224)
    plt.pcolormesh(x, y, z)
    for ax in fig.get_axes():
        for label in ax.get_xticklabels():
            label.set_ha('right')
            label.set_rotation(30)


@image_comparison(baseline_images=['pcolor_datetime_axis'],
                  extensions=['png'], remove_text=False)
def test_pcolor_datetime_axis():
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.4, top=0.98, bottom=.15)
    base = datetime.datetime(2013, 1, 1)
    x = np.array([base + datetime.timedelta(days=d) for d in range(21)])
    y = np.arange(21)
    z1, z2 = np.meshgrid(np.arange(20), np.arange(20))
    z = z1 * z2
    plt.subplot(221)
    plt.pcolor(x[:-1], y[:-1], z)
    plt.subplot(222)
    plt.pcolor(x, y, z)
    x = np.repeat(x[np.newaxis], 21, axis=0)
    y = np.repeat(y[:, np.newaxis], 21, axis=1)
    plt.subplot(223)
    plt.pcolor(x[:-1, :-1], y[:-1, :-1], z)
    plt.subplot(224)
    plt.pcolor(x, y, z)
    for ax in fig.get_axes():
        for label in ax.get_xticklabels():
            label.set_ha('right')
            label.set_rotation(30)


def test_pcolorargs():
    n = 12
    x = np.linspace(-1.5, 1.5, n)
    y = np.linspace(-1.5, 1.5, n*2)
    X, Y = np.meshgrid(x, y)
    Z = np.sqrt(X**2 + Y**2)/5

    _, ax = plt.subplots()
    with pytest.raises(TypeError):
        ax.pcolormesh(y, x, Z)
    with pytest.raises(TypeError):
        ax.pcolormesh(X, Y, Z.T)
    with pytest.raises(TypeError):
        ax.pcolormesh(x, y, Z[:-1, :-1], shading="gouraud")
    with pytest.raises(TypeError):
        ax.pcolormesh(X, Y, Z[:-1, :-1], shading="gouraud")
    x[0] = np.NaN
    with pytest.raises(ValueError):
        ax.pcolormesh(x, y, Z[:-1, :-1])
    with np.errstate(invalid='ignore'):
        x = np.ma.array(x, mask=(x < 0))
    with pytest.raises(ValueError):
        ax.pcolormesh(x, y, Z[:-1, :-1])


@image_comparison(baseline_images=['canonical'])
def test_canonical():
    fig, ax = plt.subplots()
    ax.plot([1, 2, 3])


@image_comparison(baseline_images=['arc_angles'], remove_text=True,
                  style='default', extensions=['png'])
def test_arc_angles():
    from matplotlib import patches
    # Ellipse parameters
    w = 2
    h = 1
    centre = (0.2, 0.5)
    scale = 2

    fig, axs = plt.subplots(3, 3)
    for i, ax in enumerate(axs.flat):
        theta2 = i * 360 / 9
        theta1 = theta2 - 45

        ax.add_patch(patches.Ellipse(centre, w, h, alpha=0.3))
        ax.add_patch(patches.Arc(centre, w, h, theta1=theta1, theta2=theta2))
        # Straight lines intersecting start and end of arc
        ax.plot([scale * np.cos(np.deg2rad(theta1)) + centre[0],
                 centre[0],
                 scale * np.cos(np.deg2rad(theta2)) + centre[0]],
                [scale * np.sin(np.deg2rad(theta1)) + centre[1],
                 centre[1],
                 scale * np.sin(np.deg2rad(theta2)) + centre[1]])

        ax.set_xlim(-scale, scale)
        ax.set_ylim(-scale, scale)

        # This looks the same, but it triggers a different code path when it
        # gets large enough.
        w *= 10
        h *= 10
        centre = (centre[0] * 10, centre[1] * 10)
        scale *= 10


@image_comparison(baseline_images=['arc_ellipse'],
                  remove_text=True)
def test_arc_ellipse():
    from matplotlib import patches
    xcenter, ycenter = 0.38, 0.52
    width, height = 1e-1, 3e-1
    angle = -30

    theta = np.deg2rad(np.arange(360))
    x = width / 2. * np.cos(theta)
    y = height / 2. * np.sin(theta)

    rtheta = np.deg2rad(angle)
    R = np.array([
        [np.cos(rtheta), -np.sin(rtheta)],
        [np.sin(rtheta), np.cos(rtheta)]])

    x, y = np.dot(R, np.array([x, y]))
    x += xcenter
    y += ycenter

    fig = plt.figure()
    ax = fig.add_subplot(211, aspect='auto')
    ax.fill(x, y, alpha=0.2, facecolor='yellow', edgecolor='yellow',
            linewidth=1, zorder=1)

    e1 = patches.Arc((xcenter, ycenter), width, height,
                     angle=angle, linewidth=2, fill=False, zorder=2)

    ax.add_patch(e1)

    ax = fig.add_subplot(212, aspect='equal')
    ax.fill(x, y, alpha=0.2, facecolor='green', edgecolor='green', zorder=1)
    e2 = patches.Arc((xcenter, ycenter), width, height,
                     angle=angle, linewidth=2, fill=False, zorder=2)

    ax.add_patch(e2)


@image_comparison(baseline_images=['markevery'],
                  remove_text=True)
def test_markevery():
    x = np.linspace(0, 10, 100)
    y = np.sin(x) * np.sqrt(x/10 + 0.5)

    # check marker only plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y, 'o', label='default')
    ax.plot(x, y, 'd', markevery=None, label='mark all')
    ax.plot(x, y, 's', markevery=10, label='mark every 10')
    ax.plot(x, y, '+', markevery=(5, 20), label='mark every 5 starting at 10')
    ax.legend()


@image_comparison(baseline_images=['markevery_line'],
                  remove_text=True)
def test_markevery_line():
    x = np.linspace(0, 10, 100)
    y = np.sin(x) * np.sqrt(x/10 + 0.5)

    # check line/marker combos
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y, '-o', label='default')
    ax.plot(x, y, '-d', markevery=None, label='mark all')
    ax.plot(x, y, '-s', markevery=10, label='mark every 10')
    ax.plot(x, y, '-+', markevery=(5, 20), label='mark every 5 starting at 10')
    ax.legend()


@image_comparison(baseline_images=['markevery_linear_scales'],
                  remove_text=True)
def test_markevery_linear_scales():
    cases = [None,
             8,
             (30, 8),
             [16, 24, 30], [0, -1],
             slice(100, 200, 3),
             0.1, 0.3, 1.5,
             (0.0, 0.1), (0.45, 0.1)]

    cols = 3
    gs = matplotlib.gridspec.GridSpec(len(cases) // cols + 1, cols)

    delta = 0.11
    x = np.linspace(0, 10 - 2 * delta, 200) + delta
    y = np.sin(x) + 1.0 + delta

    for i, case in enumerate(cases):
        row = (i // cols)
        col = i % cols
        plt.subplot(gs[row, col])
        plt.title('markevery=%s' % str(case))
        plt.plot(x, y, 'o', ls='-', ms=4,  markevery=case)


@image_comparison(baseline_images=['markevery_linear_scales_zoomed'],
                  remove_text=True)
def test_markevery_linear_scales_zoomed():
    cases = [None,
             8,
             (30, 8),
             [16, 24, 30], [0, -1],
             slice(100, 200, 3),
             0.1, 0.3, 1.5,
             (0.0, 0.1), (0.45, 0.1)]

    cols = 3
    gs = matplotlib.gridspec.GridSpec(len(cases) // cols + 1, cols)

    delta = 0.11
    x = np.linspace(0, 10 - 2 * delta, 200) + delta
    y = np.sin(x) + 1.0 + delta

    for i, case in enumerate(cases):
        row = (i // cols)
        col = i % cols
        plt.subplot(gs[row, col])
        plt.title('markevery=%s' % str(case))
        plt.plot(x, y, 'o', ls='-', ms=4,  markevery=case)
        plt.xlim((6, 6.7))
        plt.ylim((1.1, 1.7))


@image_comparison(baseline_images=['markevery_log_scales'],
                  remove_text=True)
def test_markevery_log_scales():
    cases = [None,
             8,
             (30, 8),
             [16, 24, 30], [0, -1],
             slice(100, 200, 3),
             0.1, 0.3, 1.5,
             (0.0, 0.1), (0.45, 0.1)]

    cols = 3
    gs = matplotlib.gridspec.GridSpec(len(cases) // cols + 1, cols)

    delta = 0.11
    x = np.linspace(0, 10 - 2 * delta, 200) + delta
    y = np.sin(x) + 1.0 + delta

    for i, case in enumerate(cases):
        row = (i // cols)
        col = i % cols
        plt.subplot(gs[row, col])
        plt.title('markevery=%s' % str(case))
        plt.xscale('log')
        plt.yscale('log')
        plt.plot(x, y, 'o', ls='-', ms=4,  markevery=case)


@image_comparison(baseline_images=['markevery_polar'], style='default',
                  remove_text=True)
def test_markevery_polar():
    cases = [None,
             8,
             (30, 8),
             [16, 24, 30], [0, -1],
             slice(100, 200, 3),
             0.1, 0.3, 1.5,
             (0.0, 0.1), (0.45, 0.1)]

    cols = 3
    gs = matplotlib.gridspec.GridSpec(len(cases) // cols + 1, cols)

    r = np.linspace(0, 3.0, 200)
    theta = 2 * np.pi * r

    for i, case in enumerate(cases):
        row = (i // cols)
        col = i % cols
        plt.subplot(gs[row, col], polar=True)
        plt.title('markevery=%s' % str(case))
        plt.plot(theta, r, 'o', ls='-', ms=4,  markevery=case)


@image_comparison(baseline_images=['marker_edges'],
                  remove_text=True)
def test_marker_edges():
    x = np.linspace(0, 1, 10)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, np.sin(x), 'y.', ms=30.0, mew=0, mec='r')
    ax.plot(x+0.1, np.sin(x), 'y.', ms=30.0, mew=1, mec='r')
    ax.plot(x+0.2, np.sin(x), 'y.', ms=30.0, mew=2, mec='b')


@image_comparison(baseline_images=['bar_tick_label_single',
                                   'bar_tick_label_single'],
                  extensions=['png'])
def test_bar_tick_label_single():
    # From 2516: plot bar with array of string labels for x axis
    ax = plt.gca()
    ax.bar(0, 1, tick_label='0')

    # Reuse testcase from above for a labeled data test
    data = {"a": 0, "b": 1}
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax = plt.gca()
    ax.bar("a", "b", tick_label='0', data=data)


def test_bar_ticklabel_fail():
    fig, ax = plt.subplots()
    ax.bar([], [])


@image_comparison(baseline_images=['bar_tick_label_multiple'],
                  extensions=['png'])
def test_bar_tick_label_multiple():
    # From 2516: plot bar with array of string labels for x axis
    ax = plt.gca()
    ax.bar([1, 2.5], [1, 2], width=[0.2, 0.5], tick_label=['a', 'b'],
           align='center')


@image_comparison(
    baseline_images=['bar_tick_label_multiple_old_label_alignment'],
    extensions=['png'])
def test_bar_tick_label_multiple_old_alignment():
    # Test that the algnment for class is backward compatible
    matplotlib.rcParams["ytick.alignment"] = "center"
    ax = plt.gca()
    ax.bar([1, 2.5], [1, 2], width=[0.2, 0.5], tick_label=['a', 'b'],
           align='center')


@image_comparison(baseline_images=['barh_tick_label'],
                  extensions=['png'])
def test_barh_tick_label():
    # From 2516: plot barh with array of string labels for y axis
    ax = plt.gca()
    ax.barh([1, 2.5], [1, 2], height=[0.2, 0.5], tick_label=['a', 'b'],
            align='center')


@image_comparison(baseline_images=['hist_log'],
                  remove_text=True)
def test_hist_log():
    data0 = np.linspace(0, 1, 200)**3
    data = np.r_[1-data0, 1+data0]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(data, fill=False, log=True)


@image_comparison(baseline_images=['hist_bar_empty'], remove_text=True,
                  extensions=['png'])
def test_hist_bar_empty():
    # From #3886: creating hist from empty dataset raises ValueError
    ax = plt.gca()
    ax.hist([], histtype='bar')


@image_comparison(baseline_images=['hist_step_empty'], remove_text=True,
                  extensions=['png'])
def test_hist_step_empty():
    # From #3886: creating hist from empty dataset raises ValueError
    ax = plt.gca()
    ax.hist([], histtype='step')


@image_comparison(baseline_images=['hist_steplog'], remove_text=True, tol=0.1)
def test_hist_steplog():
    np.random.seed(0)
    data = np.random.standard_normal(2000)
    data += -2.0 - np.min(data)
    data_pos = data + 2.1
    data_big = data_pos + 30
    weights = np.ones_like(data) * 1.e-5

    ax = plt.subplot(4, 1, 1)
    plt.hist(data, 100, histtype='stepfilled', log=True)

    ax = plt.subplot(4, 1, 2)
    plt.hist(data_pos, 100, histtype='stepfilled', log=True)

    ax = plt.subplot(4, 1, 3)
    plt.hist(data, 100, weights=weights, histtype='stepfilled', log=True)

    ax = plt.subplot(4, 1, 4)
    plt.hist(data_big, 100, histtype='stepfilled', log=True,
             orientation='horizontal')


@image_comparison(baseline_images=['hist_step_filled'], remove_text=True,
                  extensions=['png'])
def test_hist_step_filled():
    np.random.seed(0)
    x = np.random.randn(1000, 3)
    n_bins = 10

    kwargs = [{'fill': True}, {'fill': False}, {'fill': None}, {}]*2
    types = ['step']*4+['stepfilled']*4
    fig, axes = plt.subplots(nrows=2, ncols=4)
    axes = axes.flatten()

    for kg, _type, ax in zip(kwargs, types, axes):
        ax.hist(x, n_bins, histtype=_type, stacked=True, **kg)
        ax.set_title('%s/%s' % (kg, _type))
        ax.set_ylim(ymin=-50)

    patches = axes[0].patches
    assert all([p.get_facecolor() == p.get_edgecolor() for p in patches])


@image_comparison(baseline_images=['hist_density'], extensions=['png'])
def test_hist_density():
    np.random.seed(19680801)
    data = np.random.standard_normal(2000)
    fig, ax = plt.subplots()
    ax.hist(data, density=True)


@image_comparison(baseline_images=['hist_step_log_bottom'],
                  remove_text=True, extensions=['png'])
def test_hist_step_log_bottom():
    # check that bottom doesn't get overwritten by the 'minimum' on a
    # log scale histogram (https://github.com/matplotlib/matplotlib/pull/4608)
    np.random.seed(0)
    data = np.random.standard_normal(2000)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # normal hist (should clip minimum to 1/base)
    ax.hist(data, bins=10, log=True, histtype='stepfilled',
            alpha=0.5, color='b')
    # manual bottom < 1/base (previously buggy, see #4608)
    ax.hist(data, bins=10, log=True, histtype='stepfilled',
            alpha=0.5, color='g', bottom=1e-2)
    # manual bottom > 1/base
    ax.hist(data, bins=10, log=True, histtype='stepfilled',
            alpha=0.5, color='r', bottom=0.5)
    # array bottom with some less than 1/base (should clip to 1/base)
    ax.hist(data, bins=10, log=True, histtype='stepfilled',
            alpha=0.5, color='y', bottom=np.arange(10))
    ax.set_ylim(9e-3, 1e3)


def test_hist_unequal_bins_density():
    # Test correct behavior of normalized histogram with unequal bins
    # https://github.com/matplotlib/matplotlib/issues/9557
    rng = np.random.RandomState(57483)
    t = rng.randn(100)
    bins = [-3, -1, -0.5, 0, 1, 5]
    mpl_heights, _, _ = plt.hist(t, bins=bins, density=True)
    np_heights, _ = np.histogram(t, bins=bins, density=True)
    assert_allclose(mpl_heights, np_heights)


def test_hist_datetime_datasets():
    data = [[datetime.datetime(2017, 1, 1), datetime.datetime(2017, 1, 1)],
            [datetime.datetime(2017, 1, 1), datetime.datetime(2017, 1, 2)]]
    fig, ax = plt.subplots()
    ax.hist(data, stacked=True)
    ax.hist(data, stacked=False)


def contour_dat():
    x = np.linspace(-3, 5, 150)
    y = np.linspace(-3, 5, 120)
    z = np.cos(x) + np.sin(y[:, np.newaxis])
    return x, y, z


@image_comparison(baseline_images=['contour_hatching'])
def test_contour_hatching():
    x, y, z = contour_dat()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    cs = ax.contourf(x, y, z, hatches=['-', '/', '\\', '//'],
                     cmap=plt.get_cmap('gray'),
                     extend='both', alpha=0.5)


@image_comparison(baseline_images=['contour_colorbar'])
def test_contour_colorbar():
    x, y, z = contour_dat()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    cs = ax.contourf(x, y, z, levels=np.arange(-1.8, 1.801, 0.2),
                     cmap=plt.get_cmap('RdBu'),
                     vmin=-0.6,
                     vmax=0.6,
                     extend='both')
    cs1 = ax.contour(x, y, z, levels=np.arange(-2.2, -0.599, 0.2),
                     colors=['y'],
                     linestyles='solid',
                     linewidths=2)
    cs2 = ax.contour(x, y, z, levels=np.arange(0.6, 2.2, 0.2),
                     colors=['c'],
                     linewidths=2)
    cbar = fig.colorbar(cs, ax=ax)
    cbar.add_lines(cs1)
    cbar.add_lines(cs2, erase=False)


@image_comparison(baseline_images=['hist2d', 'hist2d'])
def test_hist2d():
    np.random.seed(0)
    # make it not symmetric in case we switch x and y axis
    x = np.random.randn(100)*2+5
    y = np.random.randn(100)-2
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist2d(x, y, bins=10)

    # Reuse testcase from above for a labeled data test
    data = {"x": x, "y": y}
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist2d("x", "y", bins=10, data=data)


@image_comparison(baseline_images=['hist2d_transpose'])
def test_hist2d_transpose():
    np.random.seed(0)
    # make sure the output from np.histogram is transposed before
    # passing to pcolorfast
    x = np.array([5]*100)
    y = np.random.randn(100)-2
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist2d(x, y, bins=10)


@image_comparison(baseline_images=['scatter', 'scatter'])
def test_scatter_plot():
    fig, ax = plt.subplots()
    data = {"x": [3, 4, 2, 6], "y": [2, 5, 2, 3], "c": ['r', 'y', 'b', 'lime'],
            "s": [24, 15, 19, 29]}

    ax.scatter(data["x"], data["y"], c=data["c"], s=data["s"])

    # Reuse testcase from above for a labeled data test
    fig, ax = plt.subplots()
    ax.scatter("x", "y", c="c", s="s", data=data)


@image_comparison(baseline_images=['scatter_marker'], remove_text=True,
                  extensions=['png'])
def test_scatter_marker():
    fig, (ax0, ax1, ax2) = plt.subplots(ncols=3)
    ax0.scatter([3, 4, 2, 6], [2, 5, 2, 3],
                c=[(1, 0, 0), 'y', 'b', 'lime'],
                s=[60, 50, 40, 30],
                edgecolors=['k', 'r', 'g', 'b'],
                marker='s')
    ax1.scatter([3, 4, 2, 6], [2, 5, 2, 3],
                c=[(1, 0, 0), 'y', 'b', 'lime'],
                s=[60, 50, 40, 30],
                edgecolors=['k', 'r', 'g', 'b'],
                marker=mmarkers.MarkerStyle('o', fillstyle='top'))
    # unit area ellipse
    rx, ry = 3, 1
    area = rx * ry * np.pi
    theta = np.linspace(0, 2 * np.pi, 21)
    verts = np.column_stack([np.cos(theta) * rx / area,
                             np.sin(theta) * ry / area])
    ax2.scatter([3, 4, 2, 6], [2, 5, 2, 3],
                c=[(1, 0, 0), 'y', 'b', 'lime'],
                s=[60, 50, 40, 30],
                edgecolors=['k', 'r', 'g', 'b'],
                verts=verts)


@image_comparison(baseline_images=['scatter_2D'], remove_text=True,
                  extensions=['png'])
def test_scatter_2D():
    x = np.arange(3)
    y = np.arange(2)
    x, y = np.meshgrid(x, y)
    z = x + y
    fig, ax = plt.subplots()
    ax.scatter(x, y, c=z, s=200, edgecolors='face')


def test_scatter_color():
    # Try to catch cases where 'c' kwarg should have been used.
    with pytest.raises(ValueError):
        plt.scatter([1, 2], [1, 2], color=[0.1, 0.2])
    with pytest.raises(ValueError):
        plt.scatter([1, 2, 3], [1, 2, 3], color=[1, 2, 3])


def test_as_mpl_axes_api():
    # tests the _as_mpl_axes api
    from matplotlib.projections.polar import PolarAxes
    import matplotlib.axes as maxes

    class Polar(object):
        def __init__(self):
            self.theta_offset = 0

        def _as_mpl_axes(self):
            # implement the matplotlib axes interface
            return PolarAxes, {'theta_offset': self.theta_offset}
    prj = Polar()
    prj2 = Polar()
    prj2.theta_offset = np.pi
    prj3 = Polar()

    # testing axes creation with plt.axes
    ax = plt.axes([0, 0, 1, 1], projection=prj)
    assert type(ax) == PolarAxes
    ax_via_gca = plt.gca(projection=prj)
    assert ax_via_gca is ax
    plt.close()

    # testing axes creation with gca
    ax = plt.gca(projection=prj)
    assert type(ax) == maxes._subplots._subplot_classes[PolarAxes]
    ax_via_gca = plt.gca(projection=prj)
    assert ax_via_gca is ax
    # try getting the axes given a different polar projection
    with pytest.warns(UserWarning) as rec:
        ax_via_gca = plt.gca(projection=prj2)
        assert len(rec) == 1
        assert 'Requested projection is different' in str(rec[0].message)
    assert ax_via_gca is not ax
    assert ax.get_theta_offset() == 0
    assert ax_via_gca.get_theta_offset() == np.pi
    # try getting the axes given an == (not is) polar projection
    with pytest.warns(UserWarning):
        ax_via_gca = plt.gca(projection=prj3)
        assert len(rec) == 1
        assert 'Requested projection is different' in str(rec[0].message)
    assert ax_via_gca is ax
    plt.close()

    # testing axes creation with subplot
    ax = plt.subplot(121, projection=prj)
    assert type(ax) == maxes._subplots._subplot_classes[PolarAxes]
    plt.close()


def test_pyplot_axes():
    # test focusing of Axes in other Figure
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    assert ax1 is plt.axes(ax1)
    assert ax1 is plt.gca()
    assert fig1 is plt.gcf()
    plt.close(fig1)
    plt.close(fig2)


@image_comparison(baseline_images=['log_scales'])
def test_log_scales():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(np.log(np.linspace(0.1, 100)))
    ax.set_yscale('log', basey=5.5)
    ax.invert_yaxis()
    ax.set_xscale('log', basex=9.0)


@image_comparison(baseline_images=['stackplot_test_image',
                                   'stackplot_test_image'])
def test_stackplot():
    fig = plt.figure()
    x = np.linspace(0, 10, 10)
    y1 = 1.0 * x
    y2 = 2.0 * x + 1
    y3 = 3.0 * x + 2
    ax = fig.add_subplot(1, 1, 1)
    ax.stackplot(x, y1, y2, y3)
    ax.set_xlim((0, 10))
    ax.set_ylim((0, 70))

    # Reuse testcase from above for a labeled data test
    data = {"x": x, "y1": y1, "y2": y2, "y3": y3}
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.stackplot("x", "y1", "y2", "y3", data=data)
    ax.set_xlim((0, 10))
    ax.set_ylim((0, 70))


@image_comparison(baseline_images=['stackplot_test_baseline'],
                  remove_text=True)
def test_stackplot_baseline():
    np.random.seed(0)

    def layers(n, m):
        a = np.zeros((m, n))
        for i in range(n):
            for j in range(5):
                x = 1 / (.1 + np.random.random())
                y = 2 * np.random.random() - .5
                z = 10 / (.1 + np.random.random())
                a[:, i] += x * np.exp(-((np.arange(m) / m - y) * z) ** 2)
        return a

    d = layers(3, 100)
    d[50, :] = 0  # test for fixed weighted wiggle (issue #6313)

    fig, axs = plt.subplots(2, 2)

    axs[0, 0].stackplot(range(100), d.T, baseline='zero')
    axs[0, 1].stackplot(range(100), d.T, baseline='sym')
    axs[1, 0].stackplot(range(100), d.T, baseline='wiggle')
    axs[1, 1].stackplot(range(100), d.T, baseline='weighted_wiggle')


@image_comparison(baseline_images=['bxp_baseline'],
                  extensions=['png'],
                  savefig_kwarg={'dpi': 40},
                  style='default')
def test_bxp_baseline():
    np.random.seed(937)
    logstats = matplotlib.cbook.boxplot_stats(
        np.random.lognormal(mean=1.25, sigma=1., size=(37, 4))
    )

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.bxp(logstats)


@image_comparison(baseline_images=['bxp_rangewhis'],
                  extensions=['png'],
                  savefig_kwarg={'dpi': 40},
                  style='default')
def test_bxp_rangewhis():
    np.random.seed(937)
    logstats = matplotlib.cbook.boxplot_stats(
        np.random.lognormal(mean=1.25, sigma=1., size=(37, 4)),
        whis='range'
    )

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.bxp(logstats)


@image_comparison(baseline_images=['bxp_precentilewhis'],
                  extensions=['png'],
                  savefig_kwarg={'dpi': 40},
                  style='default')
def test_bxp_precentilewhis():
    np.random.seed(937)
    logstats = matplotlib.cbook.boxplot_stats(
        np.random.lognormal(mean=1.25, sigma=1., size=(37, 4)),
        whis=[5, 95]
    )

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.bxp(logstats)


@image_comparison(baseline_images=['bxp_with_xlabels'],
                  extensions=['png'],
                  savefig_kwarg={'dpi': 40},
                  style='default')
def test_bxp_with_xlabels():
    np.random.seed(937)
    logstats = matplotlib.cbook.boxplot_stats(
        np.random.lognormal(mean=1.25, sigma=1., size=(37, 4))
    )
    for stats, label in zip(logstats, list('ABCD')):
        stats['label'] = label

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.bxp(logstats)


@image_comparison(baseline_images=['bxp_horizontal'],
                  remove_text=True, extensions=['png'],
                  savefig_kwarg={'dpi': 40},
                  style='default',
                  tol=0.1)
def test_bxp_horizontal():
    np.random.seed(937)
    logstats = matplotlib.cbook.boxplot_stats(
        np.random.lognormal(mean=1.25, sigma=1., size=(37, 4))
    )

    fig, ax = plt.subplots()
    ax.set_xscale('log')
    ax.bxp(logstats, vert=False)


@image_comparison(baseline_images=['bxp_with_ylabels'],
                  extensions=['png'],
                  savefig_kwarg={'dpi': 40},
                  style='default',
                  tol=0.1,)
def test_bxp_with_ylabels():
    np.random.seed(937)
    logstats = matplotlib.cbook.boxplot_stats(
        np.random.lognormal(mean=1.25, sigma=1., size=(37, 4))
    )
    for stats, label in zip(logstats, list('ABCD')):
        stats['label'] = label

    fig, ax = plt.subplots()
    ax.set_xscale('log')
    ax.bxp(logstats, vert=False)


@image_comparison(baseline_images=['bxp_patchartist'],
                  remove_text=True, extensions=['png'],
                  savefig_kwarg={'dpi': 40},
                  style='default')
def test_bxp_patchartist():
    np.random.seed(937)
    logstats = matplotlib.cbook.boxplot_stats(
        np.random.lognormal(mean=1.25, sigma=1., size=(37, 4))
    )

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.bxp(logstats, patch_artist=True)


@image_comparison(baseline_images=['bxp_custompatchartist'],
                  remove_text=True, extensions=['png'],
                  savefig_kwarg={'dpi': 100},
                  style='default')
def test_bxp_custompatchartist():
    np.random.seed(937)
    logstats = matplotlib.cbook.boxplot_stats(
        np.random.lognormal(mean=1.25, sigma=1., size=(37, 4))
    )

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    boxprops = dict(facecolor='yellow', edgecolor='green', linestyle='dotted')
    ax.bxp(logstats, patch_artist=True, boxprops=boxprops)


@image_comparison(baseline_images=['bxp_customoutlier'],
                  remove_text=True, extensions=['png'],
                  savefig_kwarg={'dpi': 40},
                  style='default')
def test_bxp_customoutlier():
    np.random.seed(937)
    logstats = matplotlib.cbook.boxplot_stats(
        np.random.lognormal(mean=1.25, sigma=1., size=(37, 4))
    )

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    flierprops = dict(linestyle='none', marker='d', markerfacecolor='g')
    ax.bxp(logstats, flierprops=flierprops)


@image_comparison(baseline_images=['bxp_withmean_custompoint'],
                  remove_text=True, extensions=['png'],
                  savefig_kwarg={'dpi': 40},
                  style='default')
def test_bxp_showcustommean():
    np.random.seed(937)
    logstats = matplotlib.cbook.boxplot_stats(
        np.random.lognormal(mean=1.25, sigma=1., size=(37, 4))
    )

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    meanprops = dict(linestyle='none', marker='d', markerfacecolor='green')
    ax.bxp(logstats, showmeans=True, meanprops=meanprops)


@image_comparison(baseline_images=['bxp_custombox'],
                  remove_text=True, extensions=['png'],
                  savefig_kwarg={'dpi': 40},
                  style='default')
def test_bxp_custombox():
    np.random.seed(937)
    logstats = matplotlib.cbook.boxplot_stats(
        np.random.lognormal(mean=1.25, sigma=1., size=(37, 4))
    )

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    boxprops = dict(linestyle='--', color='b', linewidth=3)
    ax.bxp(logstats, boxprops=boxprops)


@image_comparison(baseline_images=['bxp_custommedian'],
                  remove_text=True, extensions=['png'],
                  savefig_kwarg={'dpi': 40},
                  style='default')
def test_bxp_custommedian():
    np.random.seed(937)
    logstats = matplotlib.cbook.boxplot_stats(
        np.random.lognormal(mean=1.25, sigma=1., size=(37, 4))
    )

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    medianprops = dict(linestyle='--', color='b', linewidth=3)
    ax.bxp(logstats, medianprops=medianprops)


@image_comparison(baseline_images=['bxp_customcap'],
                  remove_text=True, extensions=['png'],
                  savefig_kwarg={'dpi': 40},
                  style='default')
def test_bxp_customcap():
    np.random.seed(937)
    logstats = matplotlib.cbook.boxplot_stats(
        np.random.lognormal(mean=1.25, sigma=1., size=(37, 4))
    )

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    capprops = dict(linestyle='--', color='g', linewidth=3)
    ax.bxp(logstats, capprops=capprops)


@image_comparison(baseline_images=['bxp_customwhisker'],
                  remove_text=True, extensions=['png'],
                  savefig_kwarg={'dpi': 40},
                  style='default')
def test_bxp_customwhisker():
    np.random.seed(937)
    logstats = matplotlib.cbook.boxplot_stats(
        np.random.lognormal(mean=1.25, sigma=1., size=(37, 4))
    )

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    whiskerprops = dict(linestyle='-', color='m', linewidth=3)
    ax.bxp(logstats, whiskerprops=whiskerprops)


@image_comparison(baseline_images=['bxp_withnotch'],
                  remove_text=True, extensions=['png'],
                  savefig_kwarg={'dpi': 40},
                  style='default')
def test_bxp_shownotches():
    np.random.seed(937)
    logstats = matplotlib.cbook.boxplot_stats(
        np.random.lognormal(mean=1.25, sigma=1., size=(37, 4))
    )

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.bxp(logstats, shownotches=True)


@image_comparison(baseline_images=['bxp_nocaps'],
                  remove_text=True, extensions=['png'],
                  savefig_kwarg={'dpi': 40},
                  style='default')
def test_bxp_nocaps():
    np.random.seed(937)
    logstats = matplotlib.cbook.boxplot_stats(
        np.random.lognormal(mean=1.25, sigma=1., size=(37, 4))
    )

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.bxp(logstats, showcaps=False)


@image_comparison(baseline_images=['bxp_nobox'],
                  remove_text=True, extensions=['png'],
                  savefig_kwarg={'dpi': 40},
                  style='default')
def test_bxp_nobox():
    np.random.seed(937)
    logstats = matplotlib.cbook.boxplot_stats(
        np.random.lognormal(mean=1.25, sigma=1., size=(37, 4))
    )

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.bxp(logstats, showbox=False)


@image_comparison(baseline_images=['bxp_no_flier_stats'],
                  remove_text=True, extensions=['png'],
                  savefig_kwarg={'dpi': 40},
                  style='default')
def test_bxp_no_flier_stats():
    np.random.seed(937)
    logstats = matplotlib.cbook.boxplot_stats(
        np.random.lognormal(mean=1.25, sigma=1., size=(37, 4))
    )
    for ls in logstats:
        ls.pop('fliers', None)

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.bxp(logstats, showfliers=False)


@image_comparison(baseline_images=['bxp_withmean_point'],
                  remove_text=True, extensions=['png'],
                  savefig_kwarg={'dpi': 40},
                  style='default')
def test_bxp_showmean():
    np.random.seed(937)
    logstats = matplotlib.cbook.boxplot_stats(
        np.random.lognormal(mean=1.25, sigma=1., size=(37, 4))
    )

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.bxp(logstats, showmeans=True, meanline=False)


@image_comparison(baseline_images=['bxp_withmean_line'],
                  remove_text=True, extensions=['png'],
                  savefig_kwarg={'dpi': 40},
                  style='default')
def test_bxp_showmeanasline():
    np.random.seed(937)
    logstats = matplotlib.cbook.boxplot_stats(
        np.random.lognormal(mean=1.25, sigma=1., size=(37, 4))
    )

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.bxp(logstats, showmeans=True, meanline=True)


@image_comparison(baseline_images=['bxp_scalarwidth'],
                  remove_text=True, extensions=['png'],
                  savefig_kwarg={'dpi': 40},
                  style='default')
def test_bxp_scalarwidth():
    np.random.seed(937)
    logstats = matplotlib.cbook.boxplot_stats(
        np.random.lognormal(mean=1.25, sigma=1., size=(37, 4))
    )

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.bxp(logstats, widths=0.25)


@image_comparison(baseline_images=['bxp_customwidths'],
                  remove_text=True, extensions=['png'],
                  savefig_kwarg={'dpi': 40},
                  style='default')
def test_bxp_customwidths():
    np.random.seed(937)
    logstats = matplotlib.cbook.boxplot_stats(
        np.random.lognormal(mean=1.25, sigma=1., size=(37, 4))
    )

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.bxp(logstats, widths=[0.10, 0.25, 0.65, 0.85])


@image_comparison(baseline_images=['bxp_custompositions'],
                  remove_text=True, extensions=['png'],
                  savefig_kwarg={'dpi': 40},
                  style='default')
def test_bxp_custompositions():
    np.random.seed(937)
    logstats = matplotlib.cbook.boxplot_stats(
        np.random.lognormal(mean=1.25, sigma=1., size=(37, 4))
    )

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.bxp(logstats, positions=[1, 5, 6, 7])


def test_bxp_bad_widths():
    np.random.seed(937)
    logstats = matplotlib.cbook.boxplot_stats(
        np.random.lognormal(mean=1.25, sigma=1., size=(37, 4))
    )

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    with pytest.raises(ValueError):
        ax.bxp(logstats, widths=[1])


def test_bxp_bad_positions():
    np.random.seed(937)
    logstats = matplotlib.cbook.boxplot_stats(
        np.random.lognormal(mean=1.25, sigma=1., size=(37, 4))
    )

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    with pytest.raises(ValueError):
        ax.bxp(logstats, positions=[2, 3])


@image_comparison(baseline_images=['boxplot', 'boxplot'],
                  tol=1,
                  style='default')
def test_boxplot():
    # Randomness used for bootstrapping.
    np.random.seed(937)

    x = np.linspace(-7, 7, 140)
    x = np.hstack([-25, x, 25])
    fig, ax = plt.subplots()

    ax.boxplot([x, x], bootstrap=10000, notch=1)
    ax.set_ylim((-30, 30))

    # Reuse testcase from above for a labeled data test
    data = {"x": [x, x]}
    fig, ax = plt.subplots()
    ax.boxplot("x", bootstrap=10000, notch=1, data=data)
    ax.set_ylim((-30, 30))


@image_comparison(baseline_images=['boxplot_sym2'],
                  remove_text=True, extensions=['png'],
                  style='default')
def test_boxplot_sym2():
    # Randomness used for bootstrapping.
    np.random.seed(937)

    x = np.linspace(-7, 7, 140)
    x = np.hstack([-25, x, 25])
    fig, [ax1, ax2] = plt.subplots(1, 2)

    ax1.boxplot([x, x], bootstrap=10000, sym='^')
    ax1.set_ylim((-30, 30))

    ax2.boxplot([x, x], bootstrap=10000, sym='g')
    ax2.set_ylim((-30, 30))


@image_comparison(baseline_images=['boxplot_sym'],
                  remove_text=True, extensions=['png'],
                  savefig_kwarg={'dpi': 40},
                  style='default')
def test_boxplot_sym():
    x = np.linspace(-7, 7, 140)
    x = np.hstack([-25, x, 25])
    fig, ax = plt.subplots()

    ax.boxplot([x, x], sym='gs')
    ax.set_ylim((-30, 30))


@image_comparison(
    baseline_images=['boxplot_autorange_false_whiskers',
                     'boxplot_autorange_true_whiskers'],
    extensions=['png'],
    style='default'
)
def test_boxplot_autorange_whiskers():
    # Randomness used for bootstrapping.
    np.random.seed(937)

    x = np.ones(140)
    x = np.hstack([0, x, 2])

    fig1, ax1 = plt.subplots()
    ax1.boxplot([x, x], bootstrap=10000, notch=1)
    ax1.set_ylim((-5, 5))

    fig2, ax2 = plt.subplots()
    ax2.boxplot([x, x], bootstrap=10000, notch=1, autorange=True)
    ax2.set_ylim((-5, 5))


def _rc_test_bxp_helper(ax, rc_dict):
    x = np.linspace(-7, 7, 140)
    x = np.hstack([-25, x, 25])
    with matplotlib.rc_context(rc_dict):
        ax.boxplot([x, x])
    return ax


@image_comparison(baseline_images=['boxplot_rc_parameters'],
                  savefig_kwarg={'dpi': 100}, remove_text=True,
                  tol=1, style='default')
def test_boxplot_rc_parameters():
    # Randomness used for bootstrapping.
    np.random.seed(937)

    fig, ax = plt.subplots(3)

    rc_axis0 = {
        'boxplot.notch': True,
        'boxplot.whiskers': [5, 95],
        'boxplot.bootstrap': 10000,

        'boxplot.flierprops.color': 'b',
        'boxplot.flierprops.marker': 'o',
        'boxplot.flierprops.markerfacecolor': 'g',
        'boxplot.flierprops.markeredgecolor': 'b',
        'boxplot.flierprops.markersize': 5,
        'boxplot.flierprops.linestyle': '--',
        'boxplot.flierprops.linewidth': 2.0,

        'boxplot.boxprops.color': 'r',
        'boxplot.boxprops.linewidth': 2.0,
        'boxplot.boxprops.linestyle': '--',

        'boxplot.capprops.color': 'c',
        'boxplot.capprops.linewidth': 2.0,
        'boxplot.capprops.linestyle': '--',

        'boxplot.medianprops.color': 'k',
        'boxplot.medianprops.linewidth': 2.0,
        'boxplot.medianprops.linestyle': '--',
    }

    rc_axis1 = {
        'boxplot.vertical': False,
        'boxplot.whiskers': 'range',
        'boxplot.patchartist': True,
    }

    rc_axis2 = {
        'boxplot.whiskers': 2.0,
        'boxplot.showcaps': False,
        'boxplot.showbox': False,
        'boxplot.showfliers': False,
        'boxplot.showmeans': True,
        'boxplot.meanline': True,

        'boxplot.meanprops.color': 'c',
        'boxplot.meanprops.linewidth': 2.0,
        'boxplot.meanprops.linestyle': '--',

        'boxplot.whiskerprops.color': 'r',
        'boxplot.whiskerprops.linewidth': 2.0,
        'boxplot.whiskerprops.linestyle': '-.',
    }
    dict_list = [rc_axis0, rc_axis1, rc_axis2]
    for axis, rc_axis in zip(ax, dict_list):
        _rc_test_bxp_helper(axis, rc_axis)

    assert (matplotlib.patches.PathPatch in
            [type(t) for t in ax[1].get_children()])


@image_comparison(baseline_images=['boxplot_with_CIarray'],
                  remove_text=True, extensions=['png'],
                  savefig_kwarg={'dpi': 40}, style='default')
def test_boxplot_with_CIarray():
    # Randomness used for bootstrapping.
    np.random.seed(937)

    x = np.linspace(-7, 7, 140)
    x = np.hstack([-25, x, 25])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    CIs = np.array([[-1.5, 3.], [-1., 3.5]])

    # show 1 boxplot with mpl medians/conf. interfals, 1 with manual values
    ax.boxplot([x, x], bootstrap=10000, usermedians=[None, 1.0],
               conf_intervals=CIs, notch=1)
    ax.set_ylim((-30, 30))


@image_comparison(baseline_images=['boxplot_no_inverted_whisker'],
                  remove_text=True, extensions=['png'],
                  savefig_kwarg={'dpi': 40}, style='default')
def test_boxplot_no_weird_whisker():
    x = np.array([3, 9000, 150, 88, 350, 200000, 1400, 960],
                 dtype=np.float64)
    ax1 = plt.axes()
    ax1.boxplot(x)
    ax1.set_yscale('log')
    ax1.yaxis.grid(False, which='minor')
    ax1.xaxis.grid(False)


def test_boxplot_bad_medians_1():
    x = np.linspace(-7, 7, 140)
    x = np.hstack([-25, x, 25])
    fig, ax = plt.subplots()
    with pytest.raises(ValueError):
        ax.boxplot(x, usermedians=[1, 2])


def test_boxplot_bad_medians_2():
    x = np.linspace(-7, 7, 140)
    x = np.hstack([-25, x, 25])
    fig, ax = plt.subplots()
    with pytest.raises(ValueError):
        ax.boxplot([x, x], usermedians=[[1, 2], [1, 2]])


def test_boxplot_bad_ci_1():
    x = np.linspace(-7, 7, 140)
    x = np.hstack([-25, x, 25])
    fig, ax = plt.subplots()
    with pytest.raises(ValueError):
        ax.boxplot([x, x], conf_intervals=[[1, 2]])


def test_boxplot_zorder():
    x = np.arange(10)
    fix, ax = plt.subplots()
    assert ax.boxplot(x)['boxes'][0].get_zorder() == 2
    assert ax.boxplot(x, zorder=10)['boxes'][0].get_zorder() == 10


def test_boxplot_bad_ci_2():
    x = np.linspace(-7, 7, 140)
    x = np.hstack([-25, x, 25])
    fig, ax = plt.subplots()
    with pytest.raises(ValueError):
        ax.boxplot([x, x], conf_intervals=[[1, 2], [1]])


@image_comparison(baseline_images=['boxplot_mod_artists_after_plotting'],
                  remove_text=True, extensions=['png'],
                  savefig_kwarg={'dpi': 40}, style='default')
def test_boxplot_mod_artist_after_plotting():
    x = [0.15, 0.11, 0.06, 0.06, 0.12, 0.56, -0.56]
    fig, ax = plt.subplots()
    bp = ax.boxplot(x, sym="o")
    for key in bp:
        for obj in bp[key]:
            obj.set_color('green')


@image_comparison(baseline_images=['violinplot_vert_baseline',
                                   'violinplot_vert_baseline'],
                  extensions=['png'])
def test_vert_violinplot_baseline():
    # First 9 digits of frac(sqrt(2))
    np.random.seed(414213562)
    data = [np.random.normal(size=100) for i in range(4)]
    ax = plt.axes()
    ax.violinplot(data, positions=range(4), showmeans=0, showextrema=0,
                  showmedians=0)

    # Reuse testcase from above for a labeled data test
    data = {"d": data}
    fig, ax = plt.subplots()
    ax = plt.axes()
    ax.violinplot("d", positions=range(4), showmeans=0, showextrema=0,
                  showmedians=0, data=data)


@image_comparison(baseline_images=['violinplot_vert_showmeans'],
                  extensions=['png'])
def test_vert_violinplot_showmeans():
    ax = plt.axes()
    # First 9 digits of frac(sqrt(3))
    np.random.seed(732050807)
    data = [np.random.normal(size=100) for i in range(4)]
    ax.violinplot(data, positions=range(4), showmeans=1, showextrema=0,
                  showmedians=0)


@image_comparison(baseline_images=['violinplot_vert_showextrema'],
                  extensions=['png'])
def test_vert_violinplot_showextrema():
    ax = plt.axes()
    # First 9 digits of frac(sqrt(5))
    np.random.seed(236067977)
    data = [np.random.normal(size=100) for i in range(4)]
    ax.violinplot(data, positions=range(4), showmeans=0, showextrema=1,
                  showmedians=0)


@image_comparison(baseline_images=['violinplot_vert_showmedians'],
                  extensions=['png'])
def test_vert_violinplot_showmedians():
    ax = plt.axes()
    # First 9 digits of frac(sqrt(7))
    np.random.seed(645751311)
    data = [np.random.normal(size=100) for i in range(4)]
    ax.violinplot(data, positions=range(4), showmeans=0, showextrema=0,
                  showmedians=1)


@image_comparison(baseline_images=['violinplot_vert_showall'],
                  extensions=['png'])
def test_vert_violinplot_showall():
    ax = plt.axes()
    # First 9 digits of frac(sqrt(11))
    np.random.seed(316624790)
    data = [np.random.normal(size=100) for i in range(4)]
    ax.violinplot(data, positions=range(4), showmeans=1, showextrema=1,
                  showmedians=1)


@image_comparison(baseline_images=['violinplot_vert_custompoints_10'],
                  extensions=['png'])
def test_vert_violinplot_custompoints_10():
    ax = plt.axes()
    # First 9 digits of frac(sqrt(13))
    np.random.seed(605551275)
    data = [np.random.normal(size=100) for i in range(4)]
    ax.violinplot(data, positions=range(4), showmeans=0, showextrema=0,
                  showmedians=0, points=10)


@image_comparison(baseline_images=['violinplot_vert_custompoints_200'],
                  extensions=['png'])
def test_vert_violinplot_custompoints_200():
    ax = plt.axes()
    # First 9 digits of frac(sqrt(17))
    np.random.seed(123105625)
    data = [np.random.normal(size=100) for i in range(4)]
    ax.violinplot(data, positions=range(4), showmeans=0, showextrema=0,
                  showmedians=0, points=200)


@image_comparison(baseline_images=['violinplot_horiz_baseline'],
                  extensions=['png'])
def test_horiz_violinplot_baseline():
    ax = plt.axes()
    # First 9 digits of frac(sqrt(19))
    np.random.seed(358898943)
    data = [np.random.normal(size=100) for i in range(4)]
    ax.violinplot(data, positions=range(4), vert=False, showmeans=0,
                  showextrema=0, showmedians=0)


@image_comparison(baseline_images=['violinplot_horiz_showmedians'],
                  extensions=['png'])
def test_horiz_violinplot_showmedians():
    ax = plt.axes()
    # First 9 digits of frac(sqrt(23))
    np.random.seed(795831523)
    data = [np.random.normal(size=100) for i in range(4)]
    ax.violinplot(data, positions=range(4), vert=False, showmeans=0,
                  showextrema=0, showmedians=1)


@image_comparison(baseline_images=['violinplot_horiz_showmeans'],
                  extensions=['png'])
def test_horiz_violinplot_showmeans():
    ax = plt.axes()
    # First 9 digits of frac(sqrt(29))
    np.random.seed(385164807)
    data = [np.random.normal(size=100) for i in range(4)]
    ax.violinplot(data, positions=range(4), vert=False, showmeans=1,
                  showextrema=0, showmedians=0)


@image_comparison(baseline_images=['violinplot_horiz_showextrema'],
                  extensions=['png'])
def test_horiz_violinplot_showextrema():
    ax = plt.axes()
    # First 9 digits of frac(sqrt(31))
    np.random.seed(567764362)
    data = [np.random.normal(size=100) for i in range(4)]
    ax.violinplot(data, positions=range(4), vert=False, showmeans=0,
                  showextrema=1, showmedians=0)


@image_comparison(baseline_images=['violinplot_horiz_showall'],
                  extensions=['png'])
def test_horiz_violinplot_showall():
    ax = plt.axes()
    # First 9 digits of frac(sqrt(37))
    np.random.seed(82762530)
    data = [np.random.normal(size=100) for i in range(4)]
    ax.violinplot(data, positions=range(4), vert=False, showmeans=1,
                  showextrema=1, showmedians=1)


@image_comparison(baseline_images=['violinplot_horiz_custompoints_10'],
                  extensions=['png'])
def test_horiz_violinplot_custompoints_10():
    ax = plt.axes()
    # First 9 digits of frac(sqrt(41))
    np.random.seed(403124237)
    data = [np.random.normal(size=100) for i in range(4)]
    ax.violinplot(data, positions=range(4), vert=False, showmeans=0,
                  showextrema=0, showmedians=0, points=10)


@image_comparison(baseline_images=['violinplot_horiz_custompoints_200'],
                  extensions=['png'])
def test_horiz_violinplot_custompoints_200():
    ax = plt.axes()
    # First 9 digits of frac(sqrt(43))
    np.random.seed(557438524)
    data = [np.random.normal(size=100) for i in range(4)]
    ax.violinplot(data, positions=range(4), vert=False, showmeans=0,
                  showextrema=0, showmedians=0, points=200)


def test_violinplot_bad_positions():
    ax = plt.axes()
    # First 9 digits of frac(sqrt(47))
    np.random.seed(855654600)
    data = [np.random.normal(size=100) for i in range(4)]
    with pytest.raises(ValueError):
        ax.violinplot(data, positions=range(5))


def test_violinplot_bad_widths():
    ax = plt.axes()
    # First 9 digits of frac(sqrt(53))
    np.random.seed(280109889)
    data = [np.random.normal(size=100) for i in range(4)]
    with pytest.raises(ValueError):
        ax.violinplot(data, positions=range(4), widths=[1, 2, 3])


def test_manage_xticks():
    _, ax = plt.subplots()
    ax.set_xlim(0, 4)
    old_xlim = ax.get_xlim()
    np.random.seed(0)
    y1 = np.random.normal(10, 3, 20)
    y2 = np.random.normal(3, 1, 20)
    ax.boxplot([y1, y2], positions=[1, 2],
               manage_xticks=False)
    new_xlim = ax.get_xlim()
    assert_array_equal(old_xlim, new_xlim)


def test_tick_space_size_0():
    # allow font size to be zero, which affects ticks when there is
    # no other text in the figure.
    plt.plot([0, 1], [0, 1])
    matplotlib.rcParams.update({'font.size': 0})
    b = io.BytesIO()
    plt.savefig(b, dpi=80, format='raw')


@image_comparison(baseline_images=['errorbar_basic', 'errorbar_mixed',
                                   'errorbar_basic'])
def test_errorbar():
    x = np.arange(0.1, 4, 0.5)
    y = np.exp(-x)

    yerr = 0.1 + 0.2*np.sqrt(x)
    xerr = 0.1 + yerr

    # First illustrate basic pyplot interface, using defaults where possible.
    fig = plt.figure()
    ax = fig.gca()
    ax.errorbar(x, y, xerr=0.2, yerr=0.4)
    ax.set_title("Simplest errorbars, 0.2 in x, 0.4 in y")

    # Now switch to a more OO interface to exercise more features.
    fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True)
    ax = axs[0, 0]
    # Try a Nx1 shaped error just to check
    ax.errorbar(x, y, yerr=np.reshape(yerr, (len(y), 1)), fmt='o')
    ax.set_title('Vert. symmetric')

    # With 4 subplots, reduce the number of axis ticks to avoid crowding.
    ax.locator_params(nbins=4)

    ax = axs[0, 1]
    ax.errorbar(x, y, xerr=xerr, fmt='o', alpha=0.4)
    ax.set_title('Hor. symmetric w/ alpha')

    ax = axs[1, 0]
    ax.errorbar(x, y, yerr=[yerr, 2*yerr], xerr=[xerr, 2*xerr], fmt='--o')
    ax.set_title('H, V asymmetric')

    ax = axs[1, 1]
    ax.set_yscale('log')
    # Here we have to be careful to keep all y values positive:
    ylower = np.maximum(1e-2, y - yerr)
    yerr_lower = y - ylower

    ax.errorbar(x, y, yerr=[yerr_lower, 2*yerr], xerr=xerr,
                fmt='o', ecolor='g', capthick=2)
    ax.set_title('Mixed sym., log y')

    fig.suptitle('Variable errorbars')

    # Reuse the first testcase from above for a labeled data test
    data = {"x": x, "y": y}
    fig = plt.figure()
    ax = fig.gca()
    ax.errorbar("x", "y", xerr=0.2, yerr=0.4, data=data)
    ax.set_title("Simplest errorbars, 0.2 in x, 0.4 in y")


def test_errorbar_colorcycle():

    f, ax = plt.subplots()
    x = np.arange(10)
    y = 2*x

    e1, _, _ = ax.errorbar(x, y, c=None)
    e2, _, _ = ax.errorbar(x, 2*y, c=None)
    ln1, = ax.plot(x, 4*y)

    assert mcolors.to_rgba(e1.get_color()) == mcolors.to_rgba('C0')
    assert mcolors.to_rgba(e2.get_color()) == mcolors.to_rgba('C1')
    assert mcolors.to_rgba(ln1.get_color()) == mcolors.to_rgba('C2')


def test_errorbar_shape():
    fig = plt.figure()
    ax = fig.gca()

    x = np.arange(0.1, 4, 0.5)
    y = np.exp(-x)
    yerr1 = 0.1 + 0.2*np.sqrt(x)
    yerr = np.vstack((yerr1, 2*yerr1)).T
    xerr = 0.1 + yerr

    with pytest.raises(ValueError):
        ax.errorbar(x, y, yerr=yerr, fmt='o')
    with pytest.raises(ValueError):
        ax.errorbar(x, y, xerr=xerr, fmt='o')
    with pytest.raises(ValueError):
        ax.errorbar(x, y, yerr=yerr, xerr=xerr, fmt='o')


@image_comparison(baseline_images=['errorbar_limits'])
def test_errorbar_limits():
    x = np.arange(0.5, 5.5, 0.5)
    y = np.exp(-x)
    xerr = 0.1
    yerr = 0.2
    ls = 'dotted'

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    # standard error bars
    plt.errorbar(x, y, xerr=xerr, yerr=yerr, ls=ls, color='blue')

    # including upper limits
    uplims = np.zeros_like(x)
    uplims[[1, 5, 9]] = True
    plt.errorbar(x, y+0.5, xerr=xerr, yerr=yerr, uplims=uplims, ls=ls,
                 color='green')

    # including lower limits
    lolims = np.zeros_like(x)
    lolims[[2, 4, 8]] = True
    plt.errorbar(x, y+1.0, xerr=xerr, yerr=yerr, lolims=lolims, ls=ls,
                 color='red')

    # including upper and lower limits
    plt.errorbar(x, y+1.5, marker='o', ms=8, xerr=xerr, yerr=yerr,
                 lolims=lolims, uplims=uplims, ls=ls, color='magenta')

    # including xlower and xupper limits
    xerr = 0.2
    yerr = np.zeros_like(x) + 0.2
    yerr[[3, 6]] = 0.3
    xlolims = lolims
    xuplims = uplims
    lolims = np.zeros_like(x)
    uplims = np.zeros_like(x)
    lolims[[6]] = True
    uplims[[3]] = True
    plt.errorbar(x, y+2.1, marker='o', ms=8, xerr=xerr, yerr=yerr,
                 xlolims=xlolims, xuplims=xuplims, uplims=uplims,
                 lolims=lolims, ls='none', mec='blue', capsize=0,
                 color='cyan')
    ax.set_xlim((0, 5.5))
    ax.set_title('Errorbar upper and lower limits')


def test_errobar_nonefmt():
    # Check that passing 'none' as a format still plots errorbars
    x = np.arange(5)
    y = np.arange(5)

    plotline, _, barlines = plt.errorbar(x, y, xerr=1, yerr=1, fmt='none')
    assert plotline is None
    for errbar in barlines:
        assert np.all(errbar.get_color() == mcolors.to_rgba('C0'))


@image_comparison(baseline_images=['errorbar_with_prop_cycle'],
                  extensions=['png'], style='mpl20', remove_text=True)
def test_errorbar_with_prop_cycle():
    _cycle = cycler(ls=['--', ':'], marker=['s', 's'], mfc=['k', 'w'])
    plt.rc("axes", prop_cycle=_cycle)
    fig, ax = plt.subplots()
    ax.errorbar(x=[2, 4, 10], y=[3, 2, 4], yerr=0.5)
    ax.errorbar(x=[2, 4, 10], y=[6, 4, 2], yerr=0.5)


@image_comparison(baseline_images=['hist_stacked_stepfilled',
                                   'hist_stacked_stepfilled'])
def test_hist_stacked_stepfilled():
    # make some data
    d1 = np.linspace(1, 3, 20)
    d2 = np.linspace(0, 10, 50)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist((d1, d2), histtype="stepfilled", stacked=True)

    # Reuse testcase from above for a labeled data test
    data = {"x": (d1, d2)}
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist("x", histtype="stepfilled", stacked=True, data=data)


@image_comparison(baseline_images=['hist_offset'])
def test_hist_offset():
    # make some data
    d1 = np.linspace(0, 10, 50)
    d2 = np.linspace(1, 3, 20)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(d1, bottom=5)
    ax.hist(d2, bottom=15)


@image_comparison(baseline_images=['hist_step'], extensions=['png'],
                  remove_text=True)
def test_hist_step():
    # make some data
    d1 = np.linspace(1, 3, 20)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(d1, histtype="step")
    ax.set_ylim(0, 10)
    ax.set_xlim(-1, 5)


@image_comparison(baseline_images=['hist_step_horiz'], extensions=['png'])
def test_hist_step_horiz():
    # make some data
    d1 = np.linspace(0, 10, 50)
    d2 = np.linspace(1, 3, 20)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist((d1, d2), histtype="step", orientation="horizontal")


@image_comparison(baseline_images=['hist_stacked_weights'])
def test_hist_stacked_weighted():
    # make some data
    d1 = np.linspace(0, 10, 50)
    d2 = np.linspace(1, 3, 20)
    w1 = np.linspace(0.01, 3.5, 50)
    w2 = np.linspace(0.05, 2., 20)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist((d1, d2), weights=(w1, w2), histtype="stepfilled", stacked=True)


def test_stem_args():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    x = list(xrange(10))
    y = list(xrange(10))

    # Test the call signatures
    ax.stem(y)
    ax.stem(x, y)
    ax.stem(x, y, 'r--')
    ax.stem(x, y, 'r--', basefmt='b--')


def test_stem_dates():
    fig, ax = plt.subplots(1, 1)
    from dateutil import parser
    x = parser.parse("2013-9-28 11:00:00")
    y = 100

    x1 = parser.parse("2013-9-28 12:00:00")
    y1 = 200

    ax.stem([x, x1], [y, y1], "*-")


@image_comparison(baseline_images=['hist_stacked_stepfilled_alpha'])
def test_hist_stacked_stepfilled_alpha():
    # make some data
    d1 = np.linspace(1, 3, 20)
    d2 = np.linspace(0, 10, 50)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist((d1, d2), histtype="stepfilled", stacked=True, alpha=0.5)


@image_comparison(baseline_images=['hist_stacked_step'])
def test_hist_stacked_step():
    # make some data
    d1 = np.linspace(1, 3, 20)
    d2 = np.linspace(0, 10, 50)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist((d1, d2), histtype="step", stacked=True)


@image_comparison(baseline_images=['hist_stacked_normed'])
def test_hist_stacked_normed():
    # make some data
    d1 = np.linspace(1, 3, 20)
    d2 = np.linspace(0, 10, 50)
    fig, ax = plt.subplots()
    with pytest.warns(UserWarning):
        ax.hist((d1, d2), stacked=True, normed=True)


@image_comparison(baseline_images=['hist_stacked_normed'], extensions=['png'])
def test_hist_stacked_density():
    # make some data
    d1 = np.linspace(1, 3, 20)
    d2 = np.linspace(0, 10, 50)
    fig, ax = plt.subplots()
    ax.hist((d1, d2), stacked=True, density=True)


@pytest.mark.parametrize('normed', [False, True])
@pytest.mark.parametrize('density', [False, True])
def test_hist_normed_density(normed, density):
    # Normed and density should not be used simultaneously
    d1 = np.linspace(1, 3, 20)
    d2 = np.linspace(0, 10, 50)
    fig, ax = plt.subplots()
    # test that kwargs normed and density cannot be set both.
    with pytest.raises(Exception):
        ax.hist((d1, d2), stacked=True, normed=normed, density=density)


@image_comparison(baseline_images=['hist_step_bottom'], extensions=['png'],
                  remove_text=True)
def test_hist_step_bottom():
    # make some data
    d1 = np.linspace(1, 3, 20)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(d1, bottom=np.arange(10), histtype="stepfilled")


@image_comparison(baseline_images=['hist_stacked_bar'])
def test_hist_stacked_bar():
    # make some data
    d = [[100, 100, 100, 100, 200, 320, 450, 80, 20, 600, 310, 800],
         [20, 23, 50, 11, 100, 420], [120, 120, 120, 140, 140, 150, 180],
         [60, 60, 60, 60, 300, 300, 5, 5, 5, 5, 10, 300],
         [555, 555, 555, 30, 30, 30, 30, 30, 100, 100, 100, 100, 30, 30],
         [30, 30, 30, 30, 400, 400, 400, 400, 400, 400, 400, 400]]
    colors = [(0.5759849696758961, 1.0, 0.0), (0.0, 1.0, 0.350624650815206),
              (0.0, 1.0, 0.6549834156005998), (0.0, 0.6569064625276622, 1.0),
              (0.28302699607823545, 0.0, 1.0), (0.6849123462299822, 0.0, 1.0)]
    labels = ['green', 'orange', ' yellow', 'magenta', 'black']
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(d, bins=10, histtype='barstacked', align='mid', color=colors,
            label=labels)
    ax.legend(loc='upper right', bbox_to_anchor=(1.0, 1.0), ncol=1)


def test_hist_emptydata():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist([[], range(10), range(10)], histtype="step")


@image_comparison(baseline_images=['transparent_markers'], remove_text=True)
def test_transparent_markers():
    np.random.seed(0)
    data = np.random.random(50)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(data, 'D', mfc='none', markersize=100)


@image_comparison(baseline_images=['rgba_markers'], remove_text=True)
def test_rgba_markers():
    fig, axs = plt.subplots(ncols=2)
    rcolors = [(1, 0, 0, 1), (1, 0, 0, 0.5)]
    bcolors = [(0, 0, 1, 1), (0, 0, 1, 0.5)]
    alphas = [None, 0.2]
    kw = dict(ms=100, mew=20)
    for i, alpha in enumerate(alphas):
        for j, rcolor in enumerate(rcolors):
            for k, bcolor in enumerate(bcolors):
                axs[i].plot(j+1, k+1, 'o', mfc=bcolor, mec=rcolor,
                            alpha=alpha, **kw)
                axs[i].plot(j+1, k+3, 'x', mec=rcolor, alpha=alpha, **kw)
    for ax in axs:
        ax.axis([-1, 4, 0, 5])


@image_comparison(baseline_images=['mollweide_grid'], remove_text=True)
def test_mollweide_grid():
    # test that both horizontal and vertical gridlines appear on the Mollweide
    # projection
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='mollweide')
    ax.grid()


def test_mollweide_forward_inverse_closure():
    # test that the round-trip Mollweide forward->inverse transformation is an
    # approximate identity
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='mollweide')

    # set up 1-degree grid in longitude, latitude
    lon = np.linspace(-np.pi, np.pi, 360)
    lat = np.linspace(-np.pi / 2.0, np.pi / 2.0, 180)
    lon, lat = np.meshgrid(lon, lat)
    ll = np.vstack((lon.flatten(), lat.flatten())).T

    # perform forward transform
    xy = ax.transProjection.transform(ll)

    # perform inverse transform
    ll2 = ax.transProjection.inverted().transform(xy)

    # compare
    np.testing.assert_array_almost_equal(ll, ll2, 3)


def test_mollweide_inverse_forward_closure():
    # test that the round-trip Mollweide inverse->forward transformation is an
    # approximate identity
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='mollweide')

    # set up grid in x, y
    x = np.linspace(0, 1, 500)
    x, y = np.meshgrid(x, x)
    xy = np.vstack((x.flatten(), y.flatten())).T

    # perform inverse transform
    ll = ax.transProjection.inverted().transform(xy)

    # perform forward transform
    xy2 = ax.transProjection.transform(ll)

    # compare
    np.testing.assert_array_almost_equal(xy, xy2, 3)


@image_comparison(baseline_images=['test_alpha'], remove_text=True)
def test_alpha():
    np.random.seed(0)
    data = np.random.random(50)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # alpha=.5 markers, solid line
    ax.plot(data, '-D', color=[1, 0, 0], mfc=[1, 0, 0, .5],
            markersize=20, lw=10)

    # everything solid by kwarg
    ax.plot(data + 2, '-D', color=[1, 0, 0, .5], mfc=[1, 0, 0, .5],
            markersize=20, lw=10,
            alpha=1)

    # everything alpha=.5 by kwarg
    ax.plot(data + 4, '-D', color=[1, 0, 0], mfc=[1, 0, 0],
            markersize=20, lw=10,
            alpha=.5)

    # everything alpha=.5 by colors
    ax.plot(data + 6, '-D', color=[1, 0, 0, .5], mfc=[1, 0, 0, .5],
            markersize=20, lw=10)

    # alpha=.5 line, solid markers
    ax.plot(data + 8, '-D', color=[1, 0, 0, .5], mfc=[1, 0, 0],
            markersize=20, lw=10)


@image_comparison(baseline_images=['eventplot', 'eventplot'], remove_text=True)
def test_eventplot():
    '''
    test that eventplot produces the correct output
    '''
    np.random.seed(0)

    data1 = np.random.random([32, 20]).tolist()
    data2 = np.random.random([6, 20]).tolist()
    data = data1 + data2
    num_datasets = len(data)

    colors1 = [[0, 1, .7]] * len(data1)
    colors2 = [[1, 0, 0],
               [0, 1, 0],
               [0, 0, 1],
               [1, .75, 0],
               [1, 0, 1],
               [0, 1, 1]]
    colors = colors1 + colors2

    lineoffsets1 = 12 + np.arange(0, len(data1)) * .33
    lineoffsets2 = [-15, -3, 1, 1.5, 6, 10]
    lineoffsets = lineoffsets1.tolist() + lineoffsets2

    linelengths1 = [.33] * len(data1)
    linelengths2 = [5, 2, 1, 1, 3, 1.5]
    linelengths = linelengths1 + linelengths2

    fig = plt.figure()
    axobj = fig.add_subplot(111)
    colls = axobj.eventplot(data, colors=colors, lineoffsets=lineoffsets,
                            linelengths=linelengths)

    num_collections = len(colls)
    assert num_collections == num_datasets

    # Reuse testcase from above for a labeled data test
    data = {"pos": data, "c": colors, "lo": lineoffsets, "ll": linelengths}
    fig = plt.figure()
    axobj = fig.add_subplot(111)
    colls = axobj.eventplot("pos", colors="c", lineoffsets="lo",
                            linelengths="ll", data=data)
    num_collections = len(colls)
    assert num_collections == num_datasets


@image_comparison(baseline_images=['test_eventplot_defaults'],
                  extensions=['png'], remove_text=True)
def test_eventplot_defaults():
    '''
    test that eventplot produces the correct output given the default params
    (see bug #3728)
    '''
    np.random.seed(0)

    data1 = np.random.random([32, 20]).tolist()
    data2 = np.random.random([6, 20]).tolist()
    data = data1 + data2

    fig = plt.figure()
    axobj = fig.add_subplot(111)
    colls = axobj.eventplot(data)


@pytest.mark.parametrize(('colors'), [
    ('0.5',),  # string color with multiple characters: not OK before #8193 fix
    ('tab:orange', 'tab:pink', 'tab:cyan', 'bLacK'),  # case-insensitive
    ('red', (0, 1, 0), None, (1, 0, 1, 0.5)),  # a tricky case mixing types
    ('rgbk',)  # len('rgbk') == len(data) and each character is a valid color
])
def test_eventplot_colors(colors):
    '''Test the *colors* parameter of eventplot. Inspired by the issue #8193.
    '''
    data = [[i] for i in range(4)]  # 4 successive events of different nature

    # Build the list of the expected colors
    expected = [c if c is not None else 'C0' for c in colors]
    # Convert the list into an array of RGBA values
    # NB: ['rgbk'] is not a valid argument for to_rgba_array, while 'rgbk' is.
    if len(expected) == 1:
        expected = expected[0]
    expected = broadcast_to(mcolors.to_rgba_array(expected), (len(data), 4))

    fig, ax = plt.subplots()
    if len(colors) == 1:  # tuple with a single string (like '0.5' or 'rgbk')
        colors = colors[0]
    collections = ax.eventplot(data, colors=colors)

    for coll, color in zip(collections, expected):
        assert_allclose(coll.get_color(), color)


@image_comparison(baseline_images=['test_eventplot_problem_kwargs'],
                  extensions=['png'], remove_text=True)
def test_eventplot_problem_kwargs():
    '''
    test that 'singular' versions of LineCollection props raise an
    IgnoredKeywordWarning rather than overriding the 'plural' versions (e.g.
    to prevent 'color' from overriding 'colors', see issue #4297)
    '''
    np.random.seed(0)

    data1 = np.random.random([20]).tolist()
    data2 = np.random.random([10]).tolist()
    data = [data1, data2]

    fig = plt.figure()
    axobj = fig.add_subplot(111)

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        colls = axobj.eventplot(data,
                                colors=['r', 'b'],
                                color=['c', 'm'],
                                linewidths=[2, 1],
                                linewidth=[1, 2],
                                linestyles=['solid', 'dashed'],
                                linestyle=['dashdot', 'dotted'])

        # check that three IgnoredKeywordWarnings were raised
        assert len(w) == 3
        assert all(issubclass(wi.category, IgnoredKeywordWarning) for wi in w)


def test_empty_eventplot():
    fig, ax = plt.subplots(1, 1)
    ax.eventplot([[]], colors=[(0.0, 0.0, 0.0, 0.0)])
    plt.draw()


@pytest.mark.parametrize('data, orientation', product(
    ([[]], [[], [0, 1]], [[0, 1], []]),
    ('_empty', 'vertical', 'horizontal', None, 'none')))
def test_eventplot_orientation(data, orientation):
    """Introduced when fixing issue #6412. """
    opts = {} if orientation == "_empty" else {'orientation': orientation}
    fig, ax = plt.subplots(1, 1)
    ax.eventplot(data, **opts)
    plt.draw()


@image_comparison(baseline_images=['marker_styles'], extensions=['png'],
                  remove_text=True)
def test_marker_styles():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for y, marker in enumerate(sorted(matplotlib.markers.MarkerStyle.markers,
                                      key=lambda x: str(type(x))+str(x))):
        ax.plot((y % 2)*5 + np.arange(10)*10, np.ones(10)*10*y, linestyle='',
                marker=marker, markersize=10+y/5, label=marker)


@image_comparison(baseline_images=['rc_markerfill'], extensions=['png'])
def test_markers_fillstyle_rcparams():
    fig, ax = plt.subplots()
    x = np.arange(7)
    for idx, (style, marker) in enumerate(
            [('top', 's'), ('bottom', 'o'), ('none', '^')]):
        matplotlib.rcParams['markers.fillstyle'] = style
        ax.plot(x+idx, marker=marker)


@image_comparison(baseline_images=['vertex_markers'], extensions=['png'],
                  remove_text=True)
def test_vertex_markers():
    data = list(xrange(10))
    marker_as_tuple = ((-1, -1), (1, -1), (1, 1), (-1, 1))
    marker_as_list = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(data, linestyle='', marker=marker_as_tuple, mfc='k')
    ax.plot(data[::-1], linestyle='', marker=marker_as_list, mfc='b')
    ax.set_xlim([-1, 10])
    ax.set_ylim([-1, 10])


@image_comparison(baseline_images=['vline_hline_zorder',
                                   'errorbar_zorder'])
def test_eb_line_zorder():
    x = list(xrange(10))

    # First illustrate basic pyplot interface, using defaults where possible.
    fig = plt.figure()
    ax = fig.gca()
    ax.plot(x, lw=10, zorder=5)
    ax.axhline(1, color='red', lw=10, zorder=1)
    ax.axhline(5, color='green', lw=10, zorder=10)
    ax.axvline(7, color='m', lw=10, zorder=7)
    ax.axvline(2, color='k', lw=10, zorder=3)

    ax.set_title("axvline and axhline zorder test")

    # Now switch to a more OO interface to exercise more features.
    fig = plt.figure()
    ax = fig.gca()
    x = list(xrange(10))
    y = np.zeros(10)
    yerr = list(xrange(10))
    ax.errorbar(x, y, yerr=yerr, zorder=5, lw=5, color='r')
    for j in range(10):
        ax.axhline(j, lw=5, color='k', zorder=j)
        ax.axhline(-j, lw=5, color='k', zorder=j)

    ax.set_title("errorbar zorder test")


@image_comparison(
    baseline_images=['vlines_basic', 'vlines_with_nan', 'vlines_masked'],
    extensions=['png']
)
def test_vlines():
    # normal
    x1 = [2, 3, 4, 5, 7]
    y1 = [2, -6, 3, 8, 2]
    fig1, ax1 = plt.subplots()
    ax1.vlines(x1, 0, y1, colors='g', linewidth=5)

    # GH #7406
    x2 = [2, 3, 4, 5, 6, 7]
    y2 = [2, -6, 3, 8, np.nan, 2]
    fig2, (ax2, ax3, ax4) = plt.subplots(nrows=3, figsize=(4, 8))
    ax2.vlines(x2, 0, y2, colors='g', linewidth=5)

    x3 = [2, 3, 4, 5, 6, 7]
    y3 = [np.nan, 2, -6, 3, 8, 2]
    ax3.vlines(x3, 0, y3, colors='r', linewidth=3, linestyle='--')

    x4 = [2, 3, 4, 5, 6, 7]
    y4 = [np.nan, 2, -6, 3, 8, np.nan]
    ax4.vlines(x4, 0, y4, colors='k', linewidth=2)

    # tweak the x-axis so we can see the lines better
    for ax in [ax1, ax2, ax3, ax4]:
        ax.set_xlim(0, 10)

    # check that the y-lims are all automatically the same
    assert ax1.get_ylim() == ax2.get_ylim()
    assert ax1.get_ylim() == ax3.get_ylim()
    assert ax1.get_ylim() == ax4.get_ylim()

    fig3, ax5 = plt.subplots()
    x5 = np.ma.masked_equal([2, 4, 6, 8, 10, 12], 8)
    ymin5 = np.ma.masked_equal([0, 1, -1, 0, 2, 1], 2)
    ymax5 = np.ma.masked_equal([13, 14, 15, 16, 17, 18], 18)
    ax5.vlines(x5, ymin5, ymax5, colors='k', linewidth=2)
    ax5.set_xlim(0, 15)


@image_comparison(
    baseline_images=['hlines_basic', 'hlines_with_nan', 'hlines_masked'],
    extensions=['png']
)
def test_hlines():
    # normal
    y1 = [2, 3, 4, 5, 7]
    x1 = [2, -6, 3, 8, 2]
    fig1, ax1 = plt.subplots()
    ax1.hlines(y1, 0, x1, colors='g', linewidth=5)

    # GH #7406
    y2 = [2, 3, 4, 5, 6, 7]
    x2 = [2, -6, 3, 8, np.nan, 2]
    fig2, (ax2, ax3, ax4) = plt.subplots(nrows=3, figsize=(4, 8))
    ax2.hlines(y2, 0, x2, colors='g', linewidth=5)

    y3 = [2, 3, 4, 5, 6, 7]
    x3 = [np.nan, 2, -6, 3, 8, 2]
    ax3.hlines(y3, 0, x3, colors='r', linewidth=3, linestyle='--')

    y4 = [2, 3, 4, 5, 6, 7]
    x4 = [np.nan, 2, -6, 3, 8, np.nan]
    ax4.hlines(y4, 0, x4, colors='k', linewidth=2)

    # tweak the y-axis so we can see the lines better
    for ax in [ax1, ax2, ax3, ax4]:
        ax.set_ylim(0, 10)

    # check that the x-lims are all automatically the same
    assert ax1.get_xlim() == ax2.get_xlim()
    assert ax1.get_xlim() == ax3.get_xlim()
    assert ax1.get_xlim() == ax4.get_xlim()

    fig3, ax5 = plt.subplots()
    y5 = np.ma.masked_equal([2, 4, 6, 8, 10, 12], 8)
    xmin5 = np.ma.masked_equal([0, 1, -1, 0, 2, 1], 2)
    xmax5 = np.ma.masked_equal([13, 14, 15, 16, 17, 18], 18)
    ax5.hlines(y5, xmin5, xmax5, colors='k', linewidth=2)
    ax5.set_ylim(0, 15)


@image_comparison(baseline_images=['step_linestyle', 'step_linestyle'],
                  remove_text=True)
def test_step_linestyle():
    x = y = np.arange(10)

    # First illustrate basic pyplot interface, using defaults where possible.
    fig, ax_lst = plt.subplots(2, 2)
    ax_lst = ax_lst.flatten()

    ln_styles = ['-', '--', '-.', ':']

    for ax, ls in zip(ax_lst, ln_styles):
        ax.step(x, y, lw=5, linestyle=ls, where='pre')
        ax.step(x, y + 1, lw=5, linestyle=ls, where='mid')
        ax.step(x, y + 2, lw=5, linestyle=ls, where='post')
        ax.set_xlim([-1, 5])
        ax.set_ylim([-1, 7])

    # Reuse testcase from above for a labeled data test
    data = {"x": x, "y": y, "y1": y+1, "y2": y+2}
    fig, ax_lst = plt.subplots(2, 2)
    ax_lst = ax_lst.flatten()
    ln_styles = ['-', '--', '-.', ':']
    for ax, ls in zip(ax_lst, ln_styles):
        ax.step("x", "y", lw=5, linestyle=ls, where='pre', data=data)
        ax.step("x", "y1", lw=5, linestyle=ls, where='mid', data=data)
        ax.step("x", "y2", lw=5, linestyle=ls, where='post', data=data)
        ax.set_xlim([-1, 5])
        ax.set_ylim([-1, 7])


@image_comparison(baseline_images=['mixed_collection'], remove_text=True)
def test_mixed_collection():
    from matplotlib import patches
    from matplotlib import collections

    x = list(xrange(10))

    # First illustrate basic pyplot interface, using defaults where possible.
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    c = patches.Circle((8, 8), radius=4, facecolor='none', edgecolor='green')

    # PDF can optimize this one
    p1 = collections.PatchCollection([c], match_original=True)
    p1.set_offsets([[0, 0], [24, 24]])
    p1.set_linewidths([1, 5])

    # PDF can't optimize this one, because the alpha of the edge changes
    p2 = collections.PatchCollection([c], match_original=True)
    p2.set_offsets([[48, 0], [-32, -16]])
    p2.set_linewidths([1, 5])
    p2.set_edgecolors([[0, 0, 0.1, 1.0], [0, 0, 0.1, 0.5]])

    ax.patch.set_color('0.5')
    ax.add_collection(p1)
    ax.add_collection(p2)

    ax.set_xlim(0, 16)
    ax.set_ylim(0, 16)


def test_subplot_key_hash():
    ax = plt.subplot(np.float64(5.5), np.int64(1), np.float64(1.2))
    ax.twinx()
    assert ax.get_subplotspec().get_geometry() == (5, 1, 0, 0)


@image_comparison(baseline_images=['specgram_freqs',
                                   'specgram_freqs_linear'],
                  remove_text=True, extensions=['png'], tol=0.07,
                  style='default')
def test_specgram_freqs():
    '''test axes.specgram in default (psd) mode with sinusoidal stimuli'''
    n = 1000
    Fs = 10.

    fstims1 = [Fs/4, Fs/5, Fs/11]
    fstims2 = [Fs/4.7, Fs/5.6, Fs/11.9]

    NFFT = int(10 * Fs / min(fstims1 + fstims2))
    noverlap = int(NFFT / 2)
    pad_to = int(2 ** np.ceil(np.log2(NFFT)))

    x = np.arange(0, n, 1/Fs)

    y1 = np.zeros(x.size)
    y2 = np.zeros(x.size)
    for fstim1, fstim2 in zip(fstims1, fstims2):
        y1 += np.sin(fstim1 * x * np.pi * 2)
        y2 += np.sin(fstim2 * x * np.pi * 2)
    y = np.hstack([y1, y2])

    fig1 = plt.figure()
    fig2 = plt.figure()

    ax11 = fig1.add_subplot(3, 1, 1)
    ax12 = fig1.add_subplot(3, 1, 2)
    ax13 = fig1.add_subplot(3, 1, 3)

    ax21 = fig2.add_subplot(3, 1, 1)
    ax22 = fig2.add_subplot(3, 1, 2)
    ax23 = fig2.add_subplot(3, 1, 3)

    spec11 = ax11.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='default')
    spec12 = ax12.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='onesided')
    spec13 = ax13.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='twosided')

    spec21 = ax21.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='default',
                           scale='linear', norm=matplotlib.colors.LogNorm())
    spec22 = ax22.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='onesided',
                           scale='linear', norm=matplotlib.colors.LogNorm())
    spec23 = ax23.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='twosided',
                           scale='linear', norm=matplotlib.colors.LogNorm())


@image_comparison(baseline_images=['specgram_noise',
                                   'specgram_noise_linear'],
                  remove_text=True, extensions=['png'], tol=0.01,
                  style='default')
def test_specgram_noise():
    '''test axes.specgram in default (psd) mode with noise stimuli'''
    np.random.seed(0)

    n = 1000
    Fs = 10.

    NFFT = int(10 * Fs / 11)
    noverlap = int(NFFT / 2)
    pad_to = int(2 ** np.ceil(np.log2(NFFT)))

    y1 = np.random.standard_normal(n)
    y2 = np.random.rand(n)
    y = np.hstack([y1, y2])

    fig1 = plt.figure()
    fig2 = plt.figure()

    ax11 = fig1.add_subplot(3, 1, 1)
    ax12 = fig1.add_subplot(3, 1, 2)
    ax13 = fig1.add_subplot(3, 1, 3)

    ax21 = fig2.add_subplot(3, 1, 1)
    ax22 = fig2.add_subplot(3, 1, 2)
    ax23 = fig2.add_subplot(3, 1, 3)

    spec11 = ax11.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='default')
    spec12 = ax12.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='onesided')
    spec13 = ax13.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='twosided')

    spec21 = ax21.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='default',
                           scale='linear', norm=matplotlib.colors.LogNorm())
    spec22 = ax22.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='onesided',
                           scale='linear', norm=matplotlib.colors.LogNorm())
    spec23 = ax23.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='twosided',
                           scale='linear', norm=matplotlib.colors.LogNorm())


@image_comparison(baseline_images=['specgram_magnitude_freqs',
                                   'specgram_magnitude_freqs_linear'],
                  remove_text=True, extensions=['png'], tol=0.07,
                  style='default')
def test_specgram_magnitude_freqs():
    '''test axes.specgram in magnitude mode with sinusoidal stimuli'''
    n = 1000
    Fs = 10.

    fstims1 = [Fs/4, Fs/5, Fs/11]
    fstims2 = [Fs/4.7, Fs/5.6, Fs/11.9]

    NFFT = int(100 * Fs / min(fstims1 + fstims2))
    noverlap = int(NFFT / 2)
    pad_to = int(2 ** np.ceil(np.log2(NFFT)))

    x = np.arange(0, n, 1/Fs)

    y1 = np.zeros(x.size)
    y2 = np.zeros(x.size)
    for i, (fstim1, fstim2) in enumerate(zip(fstims1, fstims2)):
        y1 += np.sin(fstim1 * x * np.pi * 2)
        y2 += np.sin(fstim2 * x * np.pi * 2)
        y1[-1] = y1[-1]/y1[-1]
        y2[-1] = y2[-1]/y2[-1]
    y = np.hstack([y1, y2])

    fig1 = plt.figure()
    fig2 = plt.figure()

    ax11 = fig1.add_subplot(3, 1, 1)
    ax12 = fig1.add_subplot(3, 1, 2)
    ax13 = fig1.add_subplot(3, 1, 3)

    ax21 = fig2.add_subplot(3, 1, 1)
    ax22 = fig2.add_subplot(3, 1, 2)
    ax23 = fig2.add_subplot(3, 1, 3)

    spec11 = ax11.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='default', mode='magnitude')
    spec12 = ax12.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='onesided', mode='magnitude')
    spec13 = ax13.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='twosided', mode='magnitude')

    spec21 = ax21.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='default', mode='magnitude',
                           scale='linear', norm=matplotlib.colors.LogNorm())
    spec22 = ax22.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='onesided', mode='magnitude',
                           scale='linear', norm=matplotlib.colors.LogNorm())
    spec23 = ax23.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='twosided', mode='magnitude',
                           scale='linear', norm=matplotlib.colors.LogNorm())


@image_comparison(baseline_images=['specgram_magnitude_noise',
                                   'specgram_magnitude_noise_linear'],
                  remove_text=True, extensions=['png'],
                  style='default')
def test_specgram_magnitude_noise():
    '''test axes.specgram in magnitude mode with noise stimuli'''
    np.random.seed(0)

    n = 1000
    Fs = 10.

    NFFT = int(10 * Fs / 11)
    noverlap = int(NFFT / 2)
    pad_to = int(2 ** np.ceil(np.log2(NFFT)))

    y1 = np.random.standard_normal(n)
    y2 = np.random.rand(n)
    y = np.hstack([y1, y2])

    fig1 = plt.figure()
    fig2 = plt.figure()

    ax11 = fig1.add_subplot(3, 1, 1)
    ax12 = fig1.add_subplot(3, 1, 2)
    ax13 = fig1.add_subplot(3, 1, 3)

    ax21 = fig2.add_subplot(3, 1, 1)
    ax22 = fig2.add_subplot(3, 1, 2)
    ax23 = fig2.add_subplot(3, 1, 3)

    spec11 = ax11.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='default', mode='magnitude')
    spec12 = ax12.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='onesided', mode='magnitude')
    spec13 = ax13.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='twosided', mode='magnitude')

    spec21 = ax21.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='default', mode='magnitude',
                           scale='linear', norm=matplotlib.colors.LogNorm())
    spec22 = ax22.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='onesided', mode='magnitude',
                           scale='linear', norm=matplotlib.colors.LogNorm())
    spec23 = ax23.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='twosided', mode='magnitude',
                           scale='linear', norm=matplotlib.colors.LogNorm())


@image_comparison(baseline_images=['specgram_angle_freqs'],
                  remove_text=True, extensions=['png'], tol=0.007,
                  style='default')
def test_specgram_angle_freqs():
    '''test axes.specgram in angle mode with sinusoidal stimuli'''
    n = 1000
    Fs = 10.

    fstims1 = [Fs/4, Fs/5, Fs/11]
    fstims2 = [Fs/4.7, Fs/5.6, Fs/11.9]

    NFFT = int(10 * Fs / min(fstims1 + fstims2))
    noverlap = int(NFFT / 2)
    pad_to = int(2 ** np.ceil(np.log2(NFFT)))

    x = np.arange(0, n, 1/Fs)

    y1 = np.zeros(x.size)
    y2 = np.zeros(x.size)
    for i, (fstim1, fstim2) in enumerate(zip(fstims1, fstims2)):
        y1 += np.sin(fstim1 * x * np.pi * 2)
        y2 += np.sin(fstim2 * x * np.pi * 2)
        y1[-1] = y1[-1]/y1[-1]
        y2[-1] = y2[-1]/y2[-1]
    y = np.hstack([y1, y2])

    fig1 = plt.figure()

    ax11 = fig1.add_subplot(3, 1, 1)
    ax12 = fig1.add_subplot(3, 1, 2)
    ax13 = fig1.add_subplot(3, 1, 3)

    spec11 = ax11.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='default', mode='angle')
    spec12 = ax12.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='onesided', mode='angle')
    spec13 = ax13.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='twosided', mode='angle')

    with pytest.raises(ValueError):
        ax11.specgram(y, NFFT=NFFT, Fs=Fs,
                      noverlap=noverlap, pad_to=pad_to, sides='default',
                      mode='phase', scale='dB')

    with pytest.raises(ValueError):
        ax12.specgram(y, NFFT=NFFT, Fs=Fs,
                      noverlap=noverlap, pad_to=pad_to, sides='onesided',
                      mode='phase', scale='dB')

    with pytest.raises(ValueError):
        ax13.specgram(y, NFFT=NFFT, Fs=Fs,
                      noverlap=noverlap, pad_to=pad_to, sides='twosided',
                      mode='phase', scale='dB')


@image_comparison(baseline_images=['specgram_angle_noise'],
                  remove_text=True, extensions=['png'],
                  style='default')
def test_specgram_noise_angle():
    '''test axes.specgram in angle mode with noise stimuli'''
    np.random.seed(0)

    n = 1000
    Fs = 10.

    NFFT = int(10 * Fs / 11)
    noverlap = int(NFFT / 2)
    pad_to = int(2 ** np.ceil(np.log2(NFFT)))

    y1 = np.random.standard_normal(n)
    y2 = np.random.rand(n)
    y = np.hstack([y1, y2])

    fig1 = plt.figure()

    ax11 = fig1.add_subplot(3, 1, 1)
    ax12 = fig1.add_subplot(3, 1, 2)
    ax13 = fig1.add_subplot(3, 1, 3)

    spec11 = ax11.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='default', mode='angle')
    spec12 = ax12.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='onesided', mode='angle')
    spec13 = ax13.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='twosided', mode='angle')

    with pytest.raises(ValueError):
        ax11.specgram(y, NFFT=NFFT, Fs=Fs,
                      noverlap=noverlap, pad_to=pad_to, sides='default',
                      mode='phase', scale='dB')

    with pytest.raises(ValueError):
        ax12.specgram(y, NFFT=NFFT, Fs=Fs,
                      noverlap=noverlap, pad_to=pad_to, sides='onesided',
                      mode='phase', scale='dB')

    with pytest.raises(ValueError):
        ax13.specgram(y, NFFT=NFFT, Fs=Fs,
                      noverlap=noverlap, pad_to=pad_to, sides='twosided',
                      mode='phase', scale='dB')


@image_comparison(baseline_images=['specgram_phase_freqs'],
                  remove_text=True, extensions=['png'],
                  style='default')
def test_specgram_freqs_phase():
    '''test axes.specgram in phase mode with sinusoidal stimuli'''
    n = 1000
    Fs = 10.

    fstims1 = [Fs/4, Fs/5, Fs/11]
    fstims2 = [Fs/4.7, Fs/5.6, Fs/11.9]

    NFFT = int(10 * Fs / min(fstims1 + fstims2))
    noverlap = int(NFFT / 2)
    pad_to = int(2 ** np.ceil(np.log2(NFFT)))

    x = np.arange(0, n, 1/Fs)

    y1 = np.zeros(x.size)
    y2 = np.zeros(x.size)
    for i, (fstim1, fstim2) in enumerate(zip(fstims1, fstims2)):
        y1 += np.sin(fstim1 * x * np.pi * 2)
        y2 += np.sin(fstim2 * x * np.pi * 2)
        y1[-1] = y1[-1]/y1[-1]
        y2[-1] = y2[-1]/y2[-1]
    y = np.hstack([y1, y2])

    fig1 = plt.figure()

    ax11 = fig1.add_subplot(3, 1, 1)
    ax12 = fig1.add_subplot(3, 1, 2)
    ax13 = fig1.add_subplot(3, 1, 3)

    spec11 = ax11.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='default', mode='phase')
    spec12 = ax12.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='onesided', mode='phase')
    spec13 = ax13.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='twosided', mode='phase')

    with pytest.raises(ValueError):
        ax11.specgram(y, NFFT=NFFT, Fs=Fs,
                      noverlap=noverlap, pad_to=pad_to, sides='default',
                      mode='phase', scale='dB')

    with pytest.raises(ValueError):
        ax12.specgram(y, NFFT=NFFT, Fs=Fs,
                      noverlap=noverlap, pad_to=pad_to, sides='onesided',
                      mode='phase', scale='dB')

    with pytest.raises(ValueError):
        ax13.specgram(y, NFFT=NFFT, Fs=Fs,
                      noverlap=noverlap, pad_to=pad_to, sides='twosided',
                      mode='phase', scale='dB')


@image_comparison(baseline_images=['specgram_phase_noise'],
                  remove_text=True, extensions=['png'],
                  style='default')
def test_specgram_noise_phase():
    '''test axes.specgram in phase mode with noise stimuli'''
    np.random.seed(0)

    n = 1000
    Fs = 10.

    NFFT = int(10 * Fs / 11)
    noverlap = int(NFFT / 2)
    pad_to = int(2 ** np.ceil(np.log2(NFFT)))

    y1 = np.random.standard_normal(n)
    y2 = np.random.rand(n)
    y = np.hstack([y1, y2])

    fig1 = plt.figure()

    ax11 = fig1.add_subplot(3, 1, 1)
    ax12 = fig1.add_subplot(3, 1, 2)
    ax13 = fig1.add_subplot(3, 1, 3)

    spec11 = ax11.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='default',
                           mode='phase', )
    spec12 = ax12.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='onesided',
                           mode='phase', )
    spec13 = ax13.specgram(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='twosided',
                           mode='phase', )

    with pytest.raises(ValueError):
        ax11.specgram(y, NFFT=NFFT, Fs=Fs,
                      noverlap=noverlap, pad_to=pad_to, sides='default',
                      mode='phase', scale='dB')

    with pytest.raises(ValueError):
        ax12.specgram(y, NFFT=NFFT, Fs=Fs,
                      noverlap=noverlap, pad_to=pad_to, sides='onesided',
                      mode='phase', scale='dB')

    with pytest.raises(ValueError):
        ax13.specgram(y, NFFT=NFFT, Fs=Fs,
                      noverlap=noverlap, pad_to=pad_to, sides='twosided',
                      mode='phase', scale='dB')


@image_comparison(baseline_images=['psd_freqs'], remove_text=True,
                  extensions=['png'])
def test_psd_freqs():
    '''test axes.psd with sinusoidal stimuli'''
    n = 10000
    Fs = 100.

    fstims1 = [Fs/4, Fs/5, Fs/11]
    fstims2 = [Fs/4.7, Fs/5.6, Fs/11.9]

    NFFT = int(1000 * Fs / min(fstims1 + fstims2))
    noverlap = int(NFFT / 2)
    pad_to = int(2 ** np.ceil(np.log2(NFFT)))

    x = np.arange(0, n, 1/Fs)

    y1 = np.zeros(x.size)
    y2 = np.zeros(x.size)
    for fstim1, fstim2 in zip(fstims1, fstims2):
        y1 += np.sin(fstim1 * x * np.pi * 2)
        y2 += np.sin(fstim2 * x * np.pi * 2)
    y = np.hstack([y1, y2])

    fig = plt.figure()
    ax1 = fig.add_subplot(3, 1, 1)
    ax2 = fig.add_subplot(3, 1, 2)
    ax3 = fig.add_subplot(3, 1, 3)

    psd1, freqs1 = ax1.psd(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='default')
    psd2, freqs2 = ax2.psd(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='onesided',
                           return_line=False)
    psd3, freqs3, line3 = ax3.psd(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                                  pad_to=pad_to, sides='twosided',
                                  return_line=True)

    ax1.set_xlabel('')
    ax2.set_xlabel('')
    ax3.set_xlabel('')
    ax1.set_ylabel('')
    ax2.set_ylabel('')
    ax3.set_ylabel('')


@image_comparison(baseline_images=['psd_noise'], remove_text=True,
                  extensions=['png'])
def test_psd_noise():
    '''test axes.psd with noise stimuli'''
    np.random.seed(0)

    n = 10000
    Fs = 100.

    NFFT = int(1000 * Fs / 11)
    noverlap = int(NFFT / 2)
    pad_to = int(2 ** np.ceil(np.log2(NFFT)))

    y1 = np.random.standard_normal(n)
    y2 = np.random.rand(n)
    y = np.hstack([y1, y2])

    fig = plt.figure()
    ax1 = fig.add_subplot(3, 1, 1)
    ax2 = fig.add_subplot(3, 1, 2)
    ax3 = fig.add_subplot(3, 1, 3)

    psd1, freqs1 = ax1.psd(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='default')
    psd2, freqs2 = ax2.psd(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='onesided',
                           return_line=False)
    psd3, freqs3, line3 = ax3.psd(y, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                                  pad_to=pad_to, sides='twosided',
                                  return_line=True)

    ax1.set_xlabel('')
    ax2.set_xlabel('')
    ax3.set_xlabel('')
    ax1.set_ylabel('')
    ax2.set_ylabel('')
    ax3.set_ylabel('')


@image_comparison(baseline_images=['csd_freqs'], remove_text=True,
                  extensions=['png'])
def test_csd_freqs():
    '''test axes.csd with sinusoidal stimuli'''
    n = 10000
    Fs = 100.

    fstims1 = [Fs/4, Fs/5, Fs/11]
    fstims2 = [Fs/4.7, Fs/5.6, Fs/11.9]

    NFFT = int(1000 * Fs / min(fstims1 + fstims2))
    noverlap = int(NFFT / 2)
    pad_to = int(2 ** np.ceil(np.log2(NFFT)))

    x = np.arange(0, n, 1/Fs)

    y1 = np.zeros(x.size)
    y2 = np.zeros(x.size)
    for fstim1, fstim2 in zip(fstims1, fstims2):
        y1 += np.sin(fstim1 * x * np.pi * 2)
        y2 += np.sin(fstim2 * x * np.pi * 2)

    fig = plt.figure()
    ax1 = fig.add_subplot(3, 1, 1)
    ax2 = fig.add_subplot(3, 1, 2)
    ax3 = fig.add_subplot(3, 1, 3)

    csd1, freqs1 = ax1.csd(y1, y2, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='default')
    csd2, freqs2 = ax2.csd(y1, y2, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='onesided',
                           return_line=False)
    csd3, freqs3, line3 = ax3.csd(y1, y2, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                                  pad_to=pad_to, sides='twosided',
                                  return_line=True)

    ax1.set_xlabel('')
    ax2.set_xlabel('')
    ax3.set_xlabel('')
    ax1.set_ylabel('')
    ax2.set_ylabel('')
    ax3.set_ylabel('')


@image_comparison(baseline_images=['csd_noise'], remove_text=True,
                  extensions=['png'])
def test_csd_noise():
    '''test axes.csd with noise stimuli'''
    np.random.seed(0)

    n = 10000
    Fs = 100.

    NFFT = int(1000 * Fs / 11)
    noverlap = int(NFFT / 2)
    pad_to = int(2 ** np.ceil(np.log2(NFFT)))

    y1 = np.random.standard_normal(n)
    y2 = np.random.rand(n)

    fig = plt.figure()
    ax1 = fig.add_subplot(3, 1, 1)
    ax2 = fig.add_subplot(3, 1, 2)
    ax3 = fig.add_subplot(3, 1, 3)

    csd1, freqs1 = ax1.csd(y1, y2, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='default')
    csd2, freqs2 = ax2.csd(y1, y2, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                           pad_to=pad_to, sides='onesided',
                           return_line=False)
    csd3, freqs3, line3 = ax3.csd(y1, y2, NFFT=NFFT, Fs=Fs, noverlap=noverlap,
                                  pad_to=pad_to, sides='twosided',
                                  return_line=True)

    ax1.set_xlabel('')
    ax2.set_xlabel('')
    ax3.set_xlabel('')
    ax1.set_ylabel('')
    ax2.set_ylabel('')
    ax3.set_ylabel('')


@image_comparison(baseline_images=['magnitude_spectrum_freqs_linear',
                                   'magnitude_spectrum_freqs_dB'],
                  remove_text=True,
                  extensions=['png'])
def test_magnitude_spectrum_freqs():
    '''test axes.magnitude_spectrum with sinusoidal stimuli'''
    n = 10000
    Fs = 100.

    fstims1 = [Fs/4, Fs/5, Fs/11]

    NFFT = int(1000 * Fs / min(fstims1))
    pad_to = int(2 ** np.ceil(np.log2(NFFT)))

    x = np.arange(0, n, 1/Fs)

    y = np.zeros(x.size)
    for i, fstim1 in enumerate(fstims1):
        y += np.sin(fstim1 * x * np.pi * 2) * 10**i
    y = y

    fig1 = plt.figure()
    fig2 = plt.figure()

    ax11 = fig1.add_subplot(3, 1, 1)
    ax12 = fig1.add_subplot(3, 1, 2)
    ax13 = fig1.add_subplot(3, 1, 3)

    ax21 = fig2.add_subplot(3, 1, 1)
    ax22 = fig2.add_subplot(3, 1, 2)
    ax23 = fig2.add_subplot(3, 1, 3)

    spec11, freqs11, line11 = ax11.magnitude_spectrum(y, Fs=Fs, pad_to=pad_to,
                                                      sides='default')
    spec12, freqs12, line12 = ax12.magnitude_spectrum(y, Fs=Fs, pad_to=pad_to,
                                                      sides='onesided')
    spec13, freqs13, line13 = ax13.magnitude_spectrum(y, Fs=Fs, pad_to=pad_to,
                                                      sides='twosided')

    spec21, freqs21, line21 = ax21.magnitude_spectrum(y, Fs=Fs, pad_to=pad_to,
                                                      sides='default',
                                                      scale='dB')
    spec22, freqs22, line22 = ax22.magnitude_spectrum(y, Fs=Fs, pad_to=pad_to,
                                                      sides='onesided',
                                                      scale='dB')
    spec23, freqs23, line23 = ax23.magnitude_spectrum(y, Fs=Fs, pad_to=pad_to,
                                                      sides='twosided',
                                                      scale='dB')

    ax11.set_xlabel('')
    ax12.set_xlabel('')
    ax13.set_xlabel('')
    ax11.set_ylabel('')
    ax12.set_ylabel('')
    ax13.set_ylabel('')

    ax21.set_xlabel('')
    ax22.set_xlabel('')
    ax23.set_xlabel('')
    ax21.set_ylabel('')
    ax22.set_ylabel('')
    ax23.set_ylabel('')


@image_comparison(baseline_images=['magnitude_spectrum_noise_linear',
                                   'magnitude_spectrum_noise_dB'],
                  remove_text=True,
                  extensions=['png'])
def test_magnitude_spectrum_noise():
    '''test axes.magnitude_spectrum with noise stimuli'''
    np.random.seed(0)

    n = 10000
    Fs = 100.

    NFFT = int(1000 * Fs / 11)
    pad_to = int(2 ** np.ceil(np.log2(NFFT)))

    y1 = np.random.standard_normal(n)
    y2 = np.random.rand(n)
    y = np.hstack([y1, y2]) - .5

    fig1 = plt.figure()
    fig2 = plt.figure()

    ax11 = fig1.add_subplot(3, 1, 1)
    ax12 = fig1.add_subplot(3, 1, 2)
    ax13 = fig1.add_subplot(3, 1, 3)

    ax21 = fig2.add_subplot(3, 1, 1)
    ax22 = fig2.add_subplot(3, 1, 2)
    ax23 = fig2.add_subplot(3, 1, 3)

    spec11, freqs11, line11 = ax11.magnitude_spectrum(y, Fs=Fs, pad_to=pad_to,
                                                      sides='default')
    spec12, freqs12, line12 = ax12.magnitude_spectrum(y, Fs=Fs, pad_to=pad_to,
                                                      sides='onesided')
    spec13, freqs13, line13 = ax13.magnitude_spectrum(y, Fs=Fs, pad_to=pad_to,
                                                      sides='twosided')

    spec21, freqs21, line21 = ax21.magnitude_spectrum(y, Fs=Fs, pad_to=pad_to,
                                                      sides='default',
                                                      scale='dB')
    spec22, freqs22, line22 = ax22.magnitude_spectrum(y, Fs=Fs, pad_to=pad_to,
                                                      sides='onesided',
                                                      scale='dB')
    spec23, freqs23, line23 = ax23.magnitude_spectrum(y, Fs=Fs, pad_to=pad_to,
                                                      sides='twosided',
                                                      scale='dB')

    ax11.set_xlabel('')
    ax12.set_xlabel('')
    ax13.set_xlabel('')
    ax11.set_ylabel('')
    ax12.set_ylabel('')
    ax13.set_ylabel('')

    ax21.set_xlabel('')
    ax22.set_xlabel('')
    ax23.set_xlabel('')
    ax21.set_ylabel('')
    ax22.set_ylabel('')
    ax23.set_ylabel('')


@image_comparison(baseline_images=['angle_spectrum_freqs'],
                  remove_text=True,
                  extensions=['png'])
def test_angle_spectrum_freqs():
    '''test axes.angle_spectrum with sinusoidal stimuli'''
    n = 10000
    Fs = 100.

    fstims1 = [Fs/4, Fs/5, Fs/11]

    NFFT = int(1000 * Fs / min(fstims1))
    pad_to = int(2 ** np.ceil(np.log2(NFFT)))

    x = np.arange(0, n, 1/Fs)

    y = np.zeros(x.size)
    for i, fstim1 in enumerate(fstims1):
        y += np.sin(fstim1 * x * np.pi * 2) * 10**i
    y = y

    fig = plt.figure()
    ax1 = fig.add_subplot(3, 1, 1)
    ax2 = fig.add_subplot(3, 1, 2)
    ax3 = fig.add_subplot(3, 1, 3)

    spec1, freqs1, line1 = ax1.angle_spectrum(y, Fs=Fs, pad_to=pad_to,
                                              sides='default')
    spec2, freqs2, line2 = ax2.angle_spectrum(y, Fs=Fs, pad_to=pad_to,
                                              sides='onesided')
    spec3, freqs3, line3 = ax3.angle_spectrum(y, Fs=Fs, pad_to=pad_to,
                                              sides='twosided')

    ax1.set_xlabel('')
    ax2.set_xlabel('')
    ax3.set_xlabel('')
    ax1.set_ylabel('')
    ax2.set_ylabel('')
    ax3.set_ylabel('')


@image_comparison(baseline_images=['angle_spectrum_noise'],
                  remove_text=True,
                  extensions=['png'])
def test_angle_spectrum_noise():
    '''test axes.angle_spectrum with noise stimuli'''
    np.random.seed(0)

    n = 10000
    Fs = 100.

    NFFT = int(1000 * Fs / 11)
    pad_to = int(2 ** np.ceil(np.log2(NFFT)))

    y1 = np.random.standard_normal(n)
    y2 = np.random.rand(n)
    y = np.hstack([y1, y2]) - .5

    fig = plt.figure()
    ax1 = fig.add_subplot(3, 1, 1)
    ax2 = fig.add_subplot(3, 1, 2)
    ax3 = fig.add_subplot(3, 1, 3)

    spec1, freqs1, line1 = ax1.angle_spectrum(y, Fs=Fs, pad_to=pad_to,
                                              sides='default')
    spec2, freqs2, line2 = ax2.angle_spectrum(y, Fs=Fs, pad_to=pad_to,
                                              sides='onesided')
    spec3, freqs3, line3 = ax3.angle_spectrum(y, Fs=Fs, pad_to=pad_to,
                                              sides='twosided')

    ax1.set_xlabel('')
    ax2.set_xlabel('')
    ax3.set_xlabel('')
    ax1.set_ylabel('')
    ax2.set_ylabel('')
    ax3.set_ylabel('')


@image_comparison(baseline_images=['phase_spectrum_freqs'],
                  remove_text=True,
                  extensions=['png'])
def test_phase_spectrum_freqs():
    '''test axes.phase_spectrum with sinusoidal stimuli'''
    n = 10000
    Fs = 100.

    fstims1 = [Fs/4, Fs/5, Fs/11]

    NFFT = int(1000 * Fs / min(fstims1))
    pad_to = int(2 ** np.ceil(np.log2(NFFT)))

    x = np.arange(0, n, 1/Fs)

    y = np.zeros(x.size)
    for i, fstim1 in enumerate(fstims1):
        y += np.sin(fstim1 * x * np.pi * 2) * 10**i
    y = y

    fig = plt.figure()
    ax1 = fig.add_subplot(3, 1, 1)
    ax2 = fig.add_subplot(3, 1, 2)
    ax3 = fig.add_subplot(3, 1, 3)

    spec1, freqs1, line1 = ax1.phase_spectrum(y, Fs=Fs, pad_to=pad_to,
                                              sides='default')
    spec2, freqs2, line2 = ax2.phase_spectrum(y, Fs=Fs, pad_to=pad_to,
                                              sides='onesided')
    spec3, freqs3, line3 = ax3.phase_spectrum(y, Fs=Fs, pad_to=pad_to,
                                              sides='twosided')

    ax1.set_xlabel('')
    ax2.set_xlabel('')
    ax3.set_xlabel('')
    ax1.set_ylabel('')
    ax2.set_ylabel('')
    ax3.set_ylabel('')


@image_comparison(baseline_images=['phase_spectrum_noise'],
                  remove_text=True,
                  extensions=['png'])
def test_phase_spectrum_noise():
    '''test axes.phase_spectrum with noise stimuli'''
    np.random.seed(0)

    n = 10000
    Fs = 100.

    NFFT = int(1000 * Fs / 11)
    pad_to = int(2 ** np.ceil(np.log2(NFFT)))

    y1 = np.random.standard_normal(n)
    y2 = np.random.rand(n)
    y = np.hstack([y1, y2]) - .5

    fig = plt.figure()
    ax1 = fig.add_subplot(3, 1, 1)
    ax2 = fig.add_subplot(3, 1, 2)
    ax3 = fig.add_subplot(3, 1, 3)

    spec1, freqs1, line1 = ax1.phase_spectrum(y, Fs=Fs, pad_to=pad_to,
                                              sides='default')
    spec2, freqs2, line2 = ax2.phase_spectrum(y, Fs=Fs, pad_to=pad_to,
                                              sides='onesided')
    spec3, freqs3, line3 = ax3.phase_spectrum(y, Fs=Fs, pad_to=pad_to,
                                              sides='twosided')

    ax1.set_xlabel('')
    ax2.set_xlabel('')
    ax3.set_xlabel('')
    ax1.set_ylabel('')
    ax2.set_ylabel('')
    ax3.set_ylabel('')


@image_comparison(baseline_images=['twin_spines'], remove_text=True,
                  extensions=['png'])
def test_twin_spines():

    def make_patch_spines_invisible(ax):
        ax.set_frame_on(True)
        ax.patch.set_visible(False)
        for sp in six.itervalues(ax.spines):
            sp.set_visible(False)

    fig = plt.figure(figsize=(4, 3))
    fig.subplots_adjust(right=0.75)

    host = fig.add_subplot(111)
    par1 = host.twinx()
    par2 = host.twinx()

    # Offset the right spine of par2.  The ticks and label have already been
    # placed on the right by twinx above.
    par2.spines["right"].set_position(("axes", 1.2))
    # Having been created by twinx, par2 has its frame off, so the line of
    # its detached spine is invisible.  First, activate the frame but make
    # the patch and spines invisible.
    make_patch_spines_invisible(par2)
    # Second, show the right spine.
    par2.spines["right"].set_visible(True)

    p1, = host.plot([0, 1, 2], [0, 1, 2], "b-")
    p2, = par1.plot([0, 1, 2], [0, 3, 2], "r-")
    p3, = par2.plot([0, 1, 2], [50, 30, 15], "g-")

    host.set_xlim(0, 2)
    host.set_ylim(0, 2)
    par1.set_ylim(0, 4)
    par2.set_ylim(1, 65)

    host.yaxis.label.set_color(p1.get_color())
    par1.yaxis.label.set_color(p2.get_color())
    par2.yaxis.label.set_color(p3.get_color())

    tkw = dict(size=4, width=1.5)
    host.tick_params(axis='y', colors=p1.get_color(), **tkw)
    par1.tick_params(axis='y', colors=p2.get_color(), **tkw)
    par2.tick_params(axis='y', colors=p3.get_color(), **tkw)
    host.tick_params(axis='x', **tkw)


@image_comparison(baseline_images=['twin_spines_on_top', 'twin_spines_on_top'],
                  extensions=['png'], remove_text=True)
def test_twin_spines_on_top():
    matplotlib.rcParams['axes.linewidth'] = 48.0
    matplotlib.rcParams['lines.linewidth'] = 48.0

    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)

    data = np.array([[1000, 1100, 1200, 1250],
                     [310, 301, 360, 400]])

    ax2 = ax1.twinx()

    ax1.plot(data[0], data[1]/1E3, color='#BEAED4')
    ax1.fill_between(data[0], data[1]/1E3, color='#BEAED4', alpha=.8)

    ax2.plot(data[0], data[1]/1E3, color='#7FC97F')
    ax2.fill_between(data[0], data[1]/1E3, color='#7FC97F', alpha=.5)

    # Reuse testcase from above for a labeled data test
    data = {"i": data[0], "j": data[1]/1E3}
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)
    ax2 = ax1.twinx()
    ax1.plot("i", "j", color='#BEAED4', data=data)
    ax1.fill_between("i", "j", color='#BEAED4', alpha=.8, data=data)
    ax2.plot("i", "j", color='#7FC97F', data=data)
    ax2.fill_between("i", "j", color='#7FC97F', alpha=.5, data=data)


def test_rcparam_grid_minor():
    orig_grid = matplotlib.rcParams['axes.grid']
    orig_locator = matplotlib.rcParams['axes.grid.which']

    matplotlib.rcParams['axes.grid'] = True

    values = (
        (('both'), (True, True)),
        (('major'), (True, False)),
        (('minor'), (False, True))
        )

    for locator, result in values:
        matplotlib.rcParams['axes.grid.which'] = locator
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        assert (ax.xaxis._gridOnMajor, ax.xaxis._gridOnMinor) == result

    matplotlib.rcParams['axes.grid'] = orig_grid
    matplotlib.rcParams['axes.grid.which'] = orig_locator


def test_vline_limit():
    fig = plt.figure()
    ax = fig.gca()
    ax.axvline(0.5)
    ax.plot([-0.1, 0, 0.2, 0.1])
    (ymin, ymax) = ax.get_ylim()
    assert_allclose(ax.get_ylim(), (-.1, .2))


def test_empty_shared_subplots():
    # empty plots with shared axes inherit limits from populated plots
    fig, axes = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True)
    axes[0].plot([1, 2, 3], [2, 4, 6])
    x0, x1 = axes[1].get_xlim()
    y0, y1 = axes[1].get_ylim()
    assert x0 <= 1
    assert x1 >= 3
    assert y0 <= 2
    assert y1 >= 6


def test_shared_with_aspect_1():
    # allow sharing one axis
    for adjustable in ['box', 'datalim']:
        fig, axes = plt.subplots(nrows=2, sharex=True)
        axes[0].set_aspect(2, adjustable=adjustable, share=True)
        assert axes[1].get_aspect() == 2
        assert axes[1].get_adjustable() == adjustable

        fig, axes = plt.subplots(nrows=2, sharex=True)
        axes[0].set_aspect(2, adjustable=adjustable)
        assert axes[1].get_aspect() == 'auto'


def test_shared_with_aspect_2():
    # Share 2 axes only with 'box':
    fig, axes = plt.subplots(nrows=2, sharex=True, sharey=True)
    axes[0].set_aspect(2, share=True)
    axes[0].plot([1, 2], [3, 4])
    axes[1].plot([3, 4], [1, 2])
    plt.draw()  # Trigger apply_aspect().
    assert axes[0].get_xlim() == axes[1].get_xlim()
    assert axes[0].get_ylim() == axes[1].get_ylim()


def test_shared_with_aspect_3():
    # Different aspect ratios:
    for adjustable in ['box', 'datalim']:
        fig, axes = plt.subplots(nrows=2, sharey=True)
        axes[0].set_aspect(2, adjustable=adjustable)
        axes[1].set_aspect(0.5, adjustable=adjustable)
        axes[0].plot([1, 2], [3, 4])
        axes[1].plot([3, 4], [1, 2])
        plt.draw()  # Trigger apply_aspect().
        assert axes[0].get_xlim() != axes[1].get_xlim()
        assert axes[0].get_ylim() == axes[1].get_ylim()
        fig_aspect = fig.bbox_inches.height / fig.bbox_inches.width
        for ax in axes:
            p = ax.get_position()
            box_aspect = p.height / p.width
            lim_aspect = ax.viewLim.height / ax.viewLim.width
            expected = fig_aspect * box_aspect / lim_aspect
            assert round(expected, 4) == round(ax.get_aspect(), 4)


@pytest.mark.parametrize('twin', ('x', 'y'))
def test_twin_with_aspect(twin):
    fig, ax = plt.subplots()
    # test twinx or twiny
    ax_twin = getattr(ax, 'twin{}'.format(twin))()
    ax.set_aspect(5)
    ax_twin.set_aspect(2)
    assert_array_equal(ax.bbox.extents,
                       ax_twin.bbox.extents)


def test_relim_visible_only():
    x1 = (0., 10.)
    y1 = (0., 10.)
    x2 = (-10., 20.)
    y2 = (-10., 30.)

    fig = matplotlib.figure.Figure()
    ax = fig.add_subplot(111)
    ax.plot(x1, y1)
    assert ax.get_xlim() == x1
    assert ax.get_ylim() == y1
    l = ax.plot(x2, y2)
    assert ax.get_xlim() == x2
    assert ax.get_ylim() == y2
    l[0].set_visible(False)
    assert ax.get_xlim() == x2
    assert ax.get_ylim() == y2

    ax.relim(visible_only=True)
    ax.autoscale_view()

    assert ax.get_xlim() == x1
    assert ax.get_ylim() == y1


def test_text_labelsize():
    """
    tests for issue #1172
    """
    fig = plt.figure()
    ax = fig.gca()
    ax.tick_params(labelsize='large')
    ax.tick_params(direction='out')


@image_comparison(baseline_images=['pie_linewidth_0', 'pie_linewidth_0',
                                   'pie_linewidth_0'],
                  extensions=['png'])
def test_pie_linewidth_0():
    # The slices will be ordered and plotted counter-clockwise.
    labels = 'Frogs', 'Hogs', 'Dogs', 'Logs'
    sizes = [15, 30, 45, 10]
    colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral']
    explode = (0, 0.1, 0, 0)  # only "explode" the 2nd slice (i.e. 'Hogs')

    plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=90,
            wedgeprops={'linewidth': 0})
    # Set aspect ratio to be equal so that pie is drawn as a circle.
    plt.axis('equal')

    # Reuse testcase from above for a labeled data test
    data = {"l": labels, "s": sizes, "c": colors, "ex": explode}
    fig = plt.figure()
    ax = fig.gca()
    ax.pie("s", explode="ex", labels="l", colors="c",
           autopct='%1.1f%%', shadow=True, startangle=90,
           wedgeprops={'linewidth': 0}, data=data)
    ax.axis('equal')

    # And again to test the pyplot functions which should also be able to be
    # called with a data kwarg
    plt.figure()
    plt.pie("s", explode="ex", labels="l", colors="c",
            autopct='%1.1f%%', shadow=True, startangle=90,
            wedgeprops={'linewidth': 0}, data=data)
    plt.axis('equal')


@image_comparison(baseline_images=['pie_center_radius'], extensions=['png'])
def test_pie_center_radius():
    # The slices will be ordered and plotted counter-clockwise.
    labels = 'Frogs', 'Hogs', 'Dogs', 'Logs'
    sizes = [15, 30, 45, 10]
    colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral']
    explode = (0, 0.1, 0, 0)  # only "explode" the 2nd slice (i.e. 'Hogs')

    plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=90,
            wedgeprops={'linewidth': 0}, center=(1, 2), radius=1.5)

    plt.annotate("Center point", xy=(1, 2), xytext=(1, 1.5),
                 arrowprops=dict(arrowstyle="->",
                                 connectionstyle="arc3"))
    # Set aspect ratio to be equal so that pie is drawn as a circle.
    plt.axis('equal')


@image_comparison(baseline_images=['pie_linewidth_2'], extensions=['png'])
def test_pie_linewidth_2():
    # The slices will be ordered and plotted counter-clockwise.
    labels = 'Frogs', 'Hogs', 'Dogs', 'Logs'
    sizes = [15, 30, 45, 10]
    colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral']
    explode = (0, 0.1, 0, 0)  # only "explode" the 2nd slice (i.e. 'Hogs')

    plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=90,
            wedgeprops={'linewidth': 2})
    # Set aspect ratio to be equal so that pie is drawn as a circle.
    plt.axis('equal')


@image_comparison(baseline_images=['pie_ccw_true'], extensions=['png'])
def test_pie_ccw_true():
    # The slices will be ordered and plotted counter-clockwise.
    labels = 'Frogs', 'Hogs', 'Dogs', 'Logs'
    sizes = [15, 30, 45, 10]
    colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral']
    explode = (0, 0.1, 0, 0)  # only "explode" the 2nd slice (i.e. 'Hogs')

    plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=90,
            counterclock=True)
    # Set aspect ratio to be equal so that pie is drawn as a circle.
    plt.axis('equal')


@image_comparison(baseline_images=['pie_frame_grid'], extensions=['png'])
def test_pie_frame_grid():
    # The slices will be ordered and plotted counter-clockwise.
    labels = 'Frogs', 'Hogs', 'Dogs', 'Logs'
    sizes = [15, 30, 45, 10]
    colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral']
    # only "explode" the 2nd slice (i.e. 'Hogs')
    explode = (0, 0.1, 0, 0)

    plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=90,
            wedgeprops={'linewidth': 0},
            frame=True, center=(2, 2))

    plt.pie(sizes[::-1], explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=90,
            wedgeprops={'linewidth': 0},
            frame=True, center=(5, 2))

    plt.pie(sizes, explode=explode[::-1], labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=90,
            wedgeprops={'linewidth': 0},
            frame=True, center=(3, 5))
    # Set aspect ratio to be equal so that pie is drawn as a circle.
    plt.axis('equal')


@image_comparison(baseline_images=['pie_rotatelabels_true'],
                  extensions=['png'])
def test_pie_rotatelabels_true():
    # The slices will be ordered and plotted counter-clockwise.
    labels = 'Hogwarts', 'Frogs', 'Dogs', 'Logs'
    sizes = [15, 30, 45, 10]
    colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral']
    explode = (0, 0.1, 0, 0)  # only "explode" the 2nd slice (i.e. 'Hogs')

    plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=90,
            rotatelabels=True)
    # Set aspect ratio to be equal so that pie is drawn as a circle.
    plt.axis('equal')


@image_comparison(baseline_images=['set_get_ticklabels'], extensions=['png'])
def test_set_get_ticklabels():
    # test issue 2246
    fig, ax = plt.subplots(2)
    ha = ['normal', 'set_x/yticklabels']

    ax[0].plot(np.arange(10))
    ax[0].set_title(ha[0])

    ax[1].plot(np.arange(10))
    ax[1].set_title(ha[1])

    # set ticklabel to 1 plot in normal way
    ax[0].set_xticklabels(('a', 'b', 'c', 'd'))
    ax[0].set_yticklabels(('11', '12', '13', '14'))

    # set ticklabel to the other plot, expect the 2 plots have same label
    # setting pass get_ticklabels return value as ticklabels argument
    ax[1].set_xticklabels(ax[0].get_xticklabels())
    ax[1].set_yticklabels(ax[0].get_yticklabels())


def test_tick_label_update():
    # test issue 9397

    fig, ax = plt.subplots()

    # Set up a dummy formatter
    def formatter_func(x, pos):
        return "unit value" if x == 1 else ""
    ax.xaxis.set_major_formatter(plt.FuncFormatter(formatter_func))

    # Force some of the x-axis ticks to be outside of the drawn range
    ax.set_xticks([-1, 0, 1, 2, 3])
    ax.set_xlim(-0.5, 2.5)

    ax.figure.canvas.draw()
    tick_texts = [tick.get_text() for tick in ax.xaxis.get_ticklabels()]
    assert tick_texts == ["", "", "unit value", "", ""]


@image_comparison(baseline_images=['o_marker_path_snap'], extensions=['png'],
                  savefig_kwarg={'dpi': 72})
def test_o_marker_path_snap():
    fig, ax = plt.subplots()
    ax.margins(.1)
    for ms in range(1, 15):
        ax.plot([1, 2, ], np.ones(2) + ms, 'o', ms=ms)

    for ms in np.linspace(1, 10, 25):
        ax.plot([3, 4, ], np.ones(2) + ms, 'o', ms=ms)


def test_margins():
    # test all ways margins can be called
    data = [1, 10]
    xmin = 0.0
    xmax = len(data) - 1.0
    ymin = min(data)
    ymax = max(data)

    fig1, ax1 = plt.subplots(1, 1)
    ax1.plot(data)
    ax1.margins(1)
    assert ax1.margins() == (1, 1)
    assert ax1.get_xlim() == (xmin - (xmax - xmin) * 1,
                              xmax + (xmax - xmin) * 1)
    assert ax1.get_ylim() == (ymin - (ymax - ymin) * 1,
                              ymax + (ymax - ymin) * 1)

    fig2, ax2 = plt.subplots(1, 1)
    ax2.plot(data)
    ax2.margins(0.5, 2)
    assert ax2.margins() == (0.5, 2)
    assert ax2.get_xlim() == (xmin - (xmax - xmin) * 0.5,
                              xmax + (xmax - xmin) * 0.5)
    assert ax2.get_ylim() == (ymin - (ymax - ymin) * 2,
                              ymax + (ymax - ymin) * 2)

    fig3, ax3 = plt.subplots(1, 1)
    ax3.plot(data)
    ax3.margins(x=-0.2, y=0.5)
    assert ax3.margins() == (-0.2, 0.5)
    assert ax3.get_xlim() == (xmin - (xmax - xmin) * -0.2,
                              xmax + (xmax - xmin) * -0.2)
    assert ax3.get_ylim() == (ymin - (ymax - ymin) * 0.5,
                              ymax + (ymax - ymin) * 0.5)


def test_length_one_hist():
    fig, ax = plt.subplots()
    ax.hist(1)
    ax.hist([1])


def test_pathological_hexbin():
    # issue #2863
    out = io.BytesIO()

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        mylist = [10] * 100
        fig, ax = plt.subplots(1, 1)
        ax.hexbin(mylist, mylist)
        fig.savefig(out)
        assert len(w) == 0


def test_color_None():
    # issue 3855
    fig, ax = plt.subplots()
    ax.plot([1, 2], [1, 2], color=None)


def test_color_alias():
    # issues 4157 and 4162
    fig, ax = plt.subplots()
    line = ax.plot([0, 1], c='lime')[0]
    assert 'lime' == line.get_color()


def test_numerical_hist_label():
    fig, ax = plt.subplots()
    ax.hist([range(15)] * 5, label=range(5))
    ax.legend()


def test_unicode_hist_label():
    fig, ax = plt.subplots()
    a = (b'\xe5\xbe\x88\xe6\xbc\x82\xe4\xba\xae, ' +
         b'r\xc3\xb6m\xc3\xa4n ch\xc3\xa4r\xc3\xa1ct\xc3\xa8rs')
    b = b'\xd7\xa9\xd7\x9c\xd7\x95\xd7\x9d'
    labels = [a.decode('utf-8'),
              'hi aardvark',
              b.decode('utf-8'),
              ]

    ax.hist([range(15)] * 3, label=labels)
    ax.legend()


def test_move_offsetlabel():
    data = np.random.random(10) * 1e-22
    fig, ax = plt.subplots()
    ax.plot(data)
    ax.yaxis.tick_right()
    assert (1, 0.5) == ax.yaxis.offsetText.get_position()


@image_comparison(baseline_images=['rc_spines'], extensions=['png'],
                  savefig_kwarg={'dpi': 40})
def test_rc_spines():
    rc_dict = {
        'axes.spines.left': False,
        'axes.spines.right': False,
        'axes.spines.top': False,
        'axes.spines.bottom': False}
    with matplotlib.rc_context(rc_dict):
        fig, ax = plt.subplots()


@image_comparison(baseline_images=['rc_grid'], extensions=['png'],
                  savefig_kwarg={'dpi': 40})
def test_rc_grid():
    fig = plt.figure()
    rc_dict0 = {
        'axes.grid': True,
        'axes.grid.axis': 'both'
    }
    rc_dict1 = {
        'axes.grid': True,
        'axes.grid.axis': 'x'
    }
    rc_dict2 = {
        'axes.grid': True,
        'axes.grid.axis': 'y'
    }
    dict_list = [rc_dict0, rc_dict1, rc_dict2]

    i = 1
    for rc_dict in dict_list:
        with matplotlib.rc_context(rc_dict):
            fig.add_subplot(3, 1, i)
            i += 1


def test_rc_tick():
    d = {'xtick.bottom': False, 'xtick.top': True,
         'ytick.left': True, 'ytick.right': False}
    with plt.rc_context(rc=d):
        fig = plt.figure()
        ax1 = fig.add_subplot(1, 1, 1)
        xax = ax1.xaxis
        yax = ax1.yaxis
        # tick1On bottom/left
        assert not xax._major_tick_kw['tick1On']
        assert xax._major_tick_kw['tick2On']
        assert not xax._minor_tick_kw['tick1On']
        assert xax._minor_tick_kw['tick2On']

        assert yax._major_tick_kw['tick1On']
        assert not yax._major_tick_kw['tick2On']
        assert yax._minor_tick_kw['tick1On']
        assert not yax._minor_tick_kw['tick2On']


def test_rc_major_minor_tick():
    d = {'xtick.top': True, 'ytick.right': True,  # Enable all ticks
         'xtick.bottom': True, 'ytick.left': True,
         # Selectively disable
         'xtick.minor.bottom': False, 'xtick.major.bottom': False,
         'ytick.major.left': False, 'ytick.minor.left': False}
    with plt.rc_context(rc=d):
        fig = plt.figure()
        ax1 = fig.add_subplot(1, 1, 1)
        xax = ax1.xaxis
        yax = ax1.yaxis
        # tick1On bottom/left
        assert not xax._major_tick_kw['tick1On']
        assert xax._major_tick_kw['tick2On']
        assert not xax._minor_tick_kw['tick1On']
        assert xax._minor_tick_kw['tick2On']

        assert not yax._major_tick_kw['tick1On']
        assert yax._major_tick_kw['tick2On']
        assert not yax._minor_tick_kw['tick1On']
        assert yax._minor_tick_kw['tick2On']


def test_square_plot():
    x = np.arange(4)
    y = np.array([1., 3., 5., 7.])
    fig, ax = plt.subplots()
    ax.plot(x, y, 'mo')
    ax.axis('square')
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    assert np.diff(xlim) == np.diff(ylim)
    assert ax.get_aspect() == 'equal'


def test_no_None():
    fig, ax = plt.subplots()
    with pytest.raises(ValueError):
        plt.plot(None)
    with pytest.raises(ValueError):
        plt.plot(None, None)


def test_pcolor_fast_non_uniform():
    Z = np.arange(6).reshape((3, 2))
    X = np.array([0, 1, 2, 10])
    Y = np.array([0, 1, 2])

    plt.figure()
    ax = plt.subplot(111)
    ax.pcolorfast(X, Y, Z.T)


def test_shared_scale():
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True)

    axs[0, 0].set_xscale("log")
    axs[0, 0].set_yscale("log")

    for ax in axs.flat:
        assert ax.get_yscale() == 'log'
        assert ax.get_xscale() == 'log'

    axs[1, 1].set_xscale("linear")
    axs[1, 1].set_yscale("linear")

    for ax in axs.flat:
        assert ax.get_yscale() == 'linear'
        assert ax.get_xscale() == 'linear'


def test_violin_point_mass():
    """Violin plot should handle point mass pdf gracefully."""
    plt.violinplot(np.array([0, 0]))


def generate_errorbar_inputs():
    base_xy = cycler('x', [np.arange(5)]) + cycler('y', [np.ones((5, ))])
    err_cycler = cycler('err', [1,
                                [1, 1, 1, 1, 1],
                                [[1, 1, 1, 1, 1],
                                 [1, 1, 1, 1, 1]],
                                [[1]] * 5,
                                np.ones(5),
                                np.ones((2, 5)),
                                np.ones((5, 1)),
                                None
                                ])
    xerr_cy = cycler('xerr', err_cycler)
    yerr_cy = cycler('yerr', err_cycler)

    empty = ((cycler('x', [[]]) + cycler('y', [[]])) *
             cycler('xerr', [[], None]) * cycler('yerr', [[], None]))
    xerr_only = base_xy * xerr_cy
    yerr_only = base_xy * yerr_cy
    both_err = base_xy * yerr_cy * xerr_cy

    test_cyclers = chain(xerr_only, yerr_only, both_err, empty)

    return test_cyclers


@pytest.mark.parametrize('kwargs', generate_errorbar_inputs())
def test_errorbar_inputs_shotgun(kwargs):
    ax = plt.gca()
    eb = ax.errorbar(**kwargs)
    eb.remove()


@image_comparison(baseline_images=["dash_offset"], remove_text=True)
def test_dash_offset():
    fig, ax = plt.subplots()
    x = np.linspace(0, 10)
    y = np.ones_like(x)
    for j in range(0, 100, 2):
        ax.plot(x, j*y, ls=(j, (10, 10)), lw=5, color='k')


def test_title_pad():
    # check that title padding puts the title in the right
    # place...
    fig, ax = plt.subplots()
    ax.set_title('aardvark', pad=30.)
    m = ax.titleOffsetTrans.get_matrix()
    assert m[1, -1] == (30. / 72. * fig.dpi)
    ax.set_title('aardvark', pad=0.)
    m = ax.titleOffsetTrans.get_matrix()
    assert m[1, -1] == 0.
    # check that it is reverted...
    ax.set_title('aardvark', pad=None)
    m = ax.titleOffsetTrans.get_matrix()
    assert m[1, -1] == (matplotlib.rcParams['axes.titlepad'] / 72. * fig.dpi)


def test_title_location_roundtrip():
    fig, ax = plt.subplots()
    ax.set_title('aardvark')
    ax.set_title('left', loc='left')
    ax.set_title('right', loc='right')

    assert 'left' == ax.get_title(loc='left')
    assert 'right' == ax.get_title(loc='right')
    assert 'aardvark' == ax.get_title()

    with pytest.raises(ValueError):
        ax.get_title(loc='foo')
    with pytest.raises(ValueError):
        ax.set_title('fail', loc='foo')


@image_comparison(baseline_images=["loglog"], remove_text=True,
                  extensions=['png'])
def test_loglog():
    fig, ax = plt.subplots()
    x = np.arange(1, 11)
    ax.loglog(x, x**3, lw=5)
    ax.tick_params(length=25, width=2)
    ax.tick_params(length=15, width=2, which='minor')


@image_comparison(baseline_images=["test_loglog_nonpos"],
                  remove_text=True, extensions=['png'], style='mpl20')
def test_loglog_nonpos():
    fig, ax = plt.subplots(3, 3)
    x = np.arange(1, 11)
    y = x**3
    y[7] = -3.
    x[4] = -10
    for nn, mcx in enumerate(['mask', 'clip', '']):
        for mm, mcy in enumerate(['mask', 'clip', '']):
            kws = {}
            if mcx:
                kws['nonposx'] = mcx
            if mcy:
                kws['nonposy'] = mcy
            ax[mm, nn].loglog(x, y**3, lw=2, **kws)


@pytest.mark.style('default')
def test_axes_margins():
    fig, ax = plt.subplots()
    ax.plot([0, 1, 2, 3])
    assert ax.get_ybound()[0] != 0

    fig, ax = plt.subplots()
    ax.bar([0, 1, 2, 3], [1, 1, 1, 1])
    assert ax.get_ybound()[0] == 0

    fig, ax = plt.subplots()
    ax.barh([0, 1, 2, 3], [1, 1, 1, 1])
    assert ax.get_xbound()[0] == 0

    fig, ax = plt.subplots()
    ax.pcolor(np.zeros((10, 10)))
    assert ax.get_xbound() == (0, 10)
    assert ax.get_ybound() == (0, 10)

    fig, ax = plt.subplots()
    ax.pcolorfast(np.zeros((10, 10)))
    assert ax.get_xbound() == (0, 10)
    assert ax.get_ybound() == (0, 10)

    fig, ax = plt.subplots()
    ax.hist(np.arange(10))
    assert ax.get_ybound()[0] == 0

    fig, ax = plt.subplots()
    ax.imshow(np.zeros((10, 10)))
    assert ax.get_xbound() == (-0.5, 9.5)
    assert ax.get_ybound() == (-0.5, 9.5)


@pytest.fixture(params=['x', 'y'])
def shared_axis_remover(request):
    def _helper_x(ax):
        ax2 = ax.twinx()
        ax2.remove()
        ax.set_xlim(0, 15)
        r = ax.xaxis.get_major_locator()()
        assert r[-1] > 14

    def _helper_y(ax):
        ax2 = ax.twiny()
        ax2.remove()
        ax.set_ylim(0, 15)
        r = ax.yaxis.get_major_locator()()
        assert r[-1] > 14

    return {"x": _helper_x, "y": _helper_y}[request.param]


@pytest.fixture(params=['gca', 'subplots', 'subplots_shared', 'add_axes'])
def shared_axes_generator(request):
    # test all of the ways to get fig/ax sets
    if request.param == 'gca':
        fig = plt.figure()
        ax = fig.gca()
    elif request.param == 'subplots':
        fig, ax = plt.subplots()
    elif request.param == 'subplots_shared':
        fig, ax_lst = plt.subplots(2, 2, sharex='all', sharey='all')
        ax = ax_lst[0][0]
    elif request.param == 'add_axes':
        fig = plt.figure()
        ax = fig.add_axes([.1, .1, .8, .8])
    return fig, ax


def test_remove_shared_axes(shared_axes_generator, shared_axis_remover):
    # test all of the ways to get fig/ax sets
    fig, ax = shared_axes_generator
    shared_axis_remover(ax)


def test_remove_shared_axes_relim():
    fig, ax_lst = plt.subplots(2, 2, sharex='all', sharey='all')
    ax = ax_lst[0][0]
    orig_xlim = ax_lst[0][1].get_xlim()
    ax.remove()
    ax.set_xlim(0, 5)
    assert_array_equal(ax_lst[0][1].get_xlim(), orig_xlim)


def test_adjust_numtick_aspect():
    fig, ax = plt.subplots()
    ax.yaxis.get_major_locator().set_params(nbins='auto')
    ax.set_xlim(0, 1000)
    ax.set_aspect('equal')
    fig.canvas.draw()
    assert len(ax.yaxis.get_major_locator()()) == 2
    ax.set_ylim(0, 1000)
    fig.canvas.draw()
    assert len(ax.yaxis.get_major_locator()()) > 2


@image_comparison(baseline_images=["auto_numticks"], style='default',
                  extensions=['png'])
def test_auto_numticks():
    # Make tiny, empty subplots, verify that there are only 3 ticks.
    fig, axes = plt.subplots(4, 4)


@image_comparison(baseline_images=["auto_numticks_log"], style='default',
                  extensions=['png'])
def test_auto_numticks_log():
    # Verify that there are not too many ticks with a large log range.
    fig, ax = plt.subplots()
    matplotlib.rcParams['axes.autolimit_mode'] = 'round_numbers'
    ax.loglog([1e-20, 1e5], [1e-16, 10])


def test_broken_barh_empty():
    fig, ax = plt.subplots()
    ax.broken_barh([], (.1, .5))


def test_pandas_pcolormesh(pd):
    time = pd.date_range('2000-01-01', periods=10)
    depth = np.arange(20)
    data = np.random.rand(20, 10)

    fig, ax = plt.subplots()
    ax.pcolormesh(time, depth, data)


def test_pandas_indexing_dates(pd):
    dates = np.arange('2005-02', '2005-03', dtype='datetime64[D]')
    values = np.sin(np.array(range(len(dates))))
    df = pd.DataFrame({'dates': dates, 'values': values})

    ax = plt.gca()

    without_zero_index = df[np.array(df.index) % 2 == 1].copy()
    ax.plot('dates', 'values', data=without_zero_index)


def test_pandas_errorbar_indexing(pd):
    df = pd.DataFrame(np.random.uniform(size=(5, 4)),
                      columns=['x', 'y', 'xe', 'ye'],
                      index=[1, 2, 3, 4, 5])
    fig, ax = plt.subplots()
    ax.errorbar('x', 'y', xerr='xe', yerr='ye', data=df)


def test_pandas_indexing_hist(pd):
    ser_1 = pd.Series(data=[1, 2, 2, 3, 3, 4, 4, 4, 4, 5])
    ser_2 = ser_1.iloc[1:]
    fig, axes = plt.subplots()
    axes.hist(ser_2)


def test_pandas_bar_align_center(pd):
    # Tests fix for issue 8767
    df = pd.DataFrame({'a': range(2), 'b': range(2)})

    fig, ax = plt.subplots(1)

    ax.bar(df.loc[df['a'] == 1, 'b'],
           df.loc[df['a'] == 1, 'b'],
           align='center')

    fig.canvas.draw()


def test_axis_set_tick_params_labelsize_labelcolor():
    # Tests fix for issue 4346
    axis_1 = plt.subplot()
    axis_1.yaxis.set_tick_params(labelsize=30, labelcolor='red',
                                 direction='out')

    # Expected values after setting the ticks
    assert axis_1.yaxis.majorTicks[0]._size == 4.0
    assert axis_1.yaxis.majorTicks[0]._color == 'k'
    assert axis_1.yaxis.majorTicks[0]._labelsize == 30.0
    assert axis_1.yaxis.majorTicks[0]._labelcolor == 'red'


def test_axes_tick_params_gridlines():
    # Now treating grid params like other Tick params
    ax = plt.subplot()
    ax.tick_params(grid_color='b', grid_linewidth=5, grid_alpha=0.5,
                   grid_linestyle='dashdot')
    for axis in ax.xaxis, ax.yaxis:
        assert axis.majorTicks[0]._grid_color == 'b'
        assert axis.majorTicks[0]._grid_linewidth == 5
        assert axis.majorTicks[0]._grid_alpha == 0.5
        assert axis.majorTicks[0]._grid_linestyle == 'dashdot'


def test_axes_tick_params_ylabelside():
    # Tests fix for issue 10267
    ax = plt.subplot()
    ax.tick_params(labelleft=False, labelright=True,
                   which='major')
    ax.tick_params(labelleft=False, labelright=True,
                   which='minor')
    # expects left false, right true
    assert ax.yaxis.majorTicks[0].label1On is False
    assert ax.yaxis.majorTicks[0].label2On is True
    assert ax.yaxis.minorTicks[0].label1On is False
    assert ax.yaxis.minorTicks[0].label2On is True


def test_axes_tick_params_xlabelside():
    # Tests fix for issue 10267
    ax = plt.subplot()
    ax.tick_params(labeltop=True, labelbottom=False,
                   which='major')
    ax.tick_params(labeltop=True, labelbottom=False,
                   which='minor')
    # expects top True, bottom False
    # label1On mapped to labelbottom
    # label2On mapped to labeltop
    assert ax.xaxis.majorTicks[0].label1On is False
    assert ax.xaxis.majorTicks[0].label2On is True
    assert ax.xaxis.minorTicks[0].label1On is False
    assert ax.xaxis.minorTicks[0].label2On is True


def test_none_kwargs():
    fig, ax = plt.subplots()
    ln, = ax.plot(range(32), linestyle=None)
    assert ln.get_linestyle() == '-'


def test_ls_ds_conflict():
    with pytest.raises(ValueError):
        plt.plot(range(32), linestyle='steps-pre:', drawstyle='steps-post')


def test_bar_uint8():
    xs = [0, 1, 2, 3]
    b = plt.bar(np.array(xs, dtype=np.uint8), [2, 3, 4, 5])
    for (patch, x) in zip(b.patches, xs):
        assert patch.xy[0] == x


@image_comparison(baseline_images=['date_timezone_x'], extensions=['png'])
def test_date_timezone_x():
    # Tests issue 5575
    time_index = [pytz.timezone('Canada/Eastern').localize(datetime.datetime(
        year=2016, month=2, day=22, hour=x)) for x in range(3)]

    # Same Timezone
    fig = plt.figure(figsize=(20, 12))
    plt.subplot(2, 1, 1)
    plt.plot_date(time_index, [3] * 3, tz='Canada/Eastern')

    # Different Timezone
    plt.subplot(2, 1, 2)
    plt.plot_date(time_index, [3] * 3, tz='UTC')


@image_comparison(baseline_images=['date_timezone_y'],
                  extensions=['png'])
def test_date_timezone_y():
    # Tests issue 5575
    time_index = [pytz.timezone('Canada/Eastern').localize(datetime.datetime(
        year=2016, month=2, day=22, hour=x)) for x in range(3)]

    # Same Timezone
    fig = plt.figure(figsize=(20, 12))
    plt.subplot(2, 1, 1)
    plt.plot_date([3] * 3,
                  time_index, tz='Canada/Eastern', xdate=False, ydate=True)

    # Different Timezone
    plt.subplot(2, 1, 2)
    plt.plot_date([3] * 3, time_index, tz='UTC', xdate=False, ydate=True)


@image_comparison(baseline_images=['date_timezone_x_and_y'],
                  extensions=['png'])
def test_date_timezone_x_and_y():
    # Tests issue 5575
    time_index = [pytz.timezone('UTC').localize(datetime.datetime(
        year=2016, month=2, day=22, hour=x)) for x in range(3)]

    # Same Timezone
    fig = plt.figure(figsize=(20, 12))
    plt.subplot(2, 1, 1)
    plt.plot_date(time_index, time_index, tz='UTC', ydate=True)

    # Different Timezone
    plt.subplot(2, 1, 2)
    plt.plot_date(time_index, time_index, tz='US/Eastern', ydate=True)


@image_comparison(baseline_images=['axisbelow'],
                  extensions=['png'], remove_text=True)
def test_axisbelow():
    # Test 'line' setting added in 6287.
    # Show only grids, not frame or ticks, to make this test
    # independent of future change to drawing order of those elements.
    fig, axs = plt.subplots(ncols=3, sharex=True, sharey=True)
    settings = (False, 'line', True)

    for ax, setting in zip(axs, settings):
        ax.plot((0, 10), (0, 10), lw=10, color='m')
        circ = mpatches.Circle((3, 3), color='r')
        ax.add_patch(circ)
        ax.grid(color='c', linestyle='-', linewidth=3)
        ax.tick_params(top=False, bottom=False,
                       left=False, right=False)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_axisbelow(setting)


def test_offset_label_color():
    # Tests issue 6440
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot([1.01e9, 1.02e9, 1.03e9])
    ax.yaxis.set_tick_params(labelcolor='red')
    assert ax.yaxis.get_offset_text().get_color() == 'red'


def test_large_offset():
    fig, ax = plt.subplots()
    ax.plot((1 + np.array([0, 1.e-12])) * 1.e27)
    fig.canvas.draw()


def test_barb_units():
    fig, ax = plt.subplots()
    dates = [datetime.datetime(2017, 7, 15, 18, i) for i in range(0, 60, 10)]
    y = np.linspace(0, 5, len(dates))
    u = v = np.linspace(0, 50, len(dates))
    ax.barbs(dates, y, u, v)


def test_quiver_units():
    fig, ax = plt.subplots()
    dates = [datetime.datetime(2017, 7, 15, 18, i) for i in range(0, 60, 10)]
    y = np.linspace(0, 5, len(dates))
    u = v = np.linspace(0, 50, len(dates))
    ax.quiver(dates, y, u, v)


def test_bar_color_cycle():
    ccov = mcolors.colorConverter.to_rgb
    fig, ax = plt.subplots()
    for j in range(5):
        ln, = ax.plot(range(3))
        brs = ax.bar(range(3), range(3))
        for br in brs:
            assert ccov(ln.get_color()) == ccov(br.get_facecolor())


def test_tick_param_label_rotation():
    fix, (ax, ax2) = plt.subplots(1, 2)
    ax.plot([0, 1], [0, 1])
    ax2.plot([0, 1], [0, 1])
    ax.xaxis.set_tick_params(which='both', rotation=75)
    ax.yaxis.set_tick_params(which='both', rotation=90)
    for text in ax.get_xticklabels(which='both'):
        assert text.get_rotation() == 75
    for text in ax.get_yticklabels(which='both'):
        assert text.get_rotation() == 90

    ax2.tick_params(axis='x', labelrotation=53)
    ax2.tick_params(axis='y', rotation=35)
    for text in ax2.get_xticklabels(which='major'):
        assert text.get_rotation() == 53
    for text in ax2.get_yticklabels(which='major'):
        assert text.get_rotation() == 35


@pytest.mark.style('default')
def test_fillbetween_cycle():
    fig, ax = plt.subplots()

    for j in range(3):
        cc = ax.fill_between(range(3), range(3))
        target = mcolors.to_rgba('C{}'.format(j))
        assert tuple(cc.get_facecolors().squeeze()) == tuple(target)

    for j in range(3, 6):
        cc = ax.fill_betweenx(range(3), range(3))
        target = mcolors.to_rgba('C{}'.format(j))
        assert tuple(cc.get_facecolors().squeeze()) == tuple(target)

    target = mcolors.to_rgba('k')

    for al in ['facecolor', 'facecolors', 'color']:
        cc = ax.fill_between(range(3), range(3), **{al: 'k'})
        assert tuple(cc.get_facecolors().squeeze()) == tuple(target)

    edge_target = mcolors.to_rgba('k')
    for j, el in enumerate(['edgecolor', 'edgecolors'], start=6):
        cc = ax.fill_between(range(3), range(3), **{el: 'k'})
        face_target = mcolors.to_rgba('C{}'.format(j))
        assert tuple(cc.get_facecolors().squeeze()) == tuple(face_target)
        assert tuple(cc.get_edgecolors().squeeze()) == tuple(edge_target)


def test_log_margins():
    plt.rcParams['axes.autolimit_mode'] = 'data'
    fig, ax = plt.subplots()
    margin = 0.05
    ax.set_xmargin(margin)
    ax.semilogx([10, 100], [10, 100])
    xlim0, xlim1 = ax.get_xlim()
    transform = ax.xaxis.get_transform()
    xlim0t, xlim1t = transform.transform([xlim0, xlim1])
    x0t, x1t = transform.transform([10, 100])
    delta = (x1t - x0t) * margin
    assert_allclose([xlim0t + delta, xlim1t - delta], [x0t, x1t])


def test_color_length_mismatch():
    N = 5
    x, y = np.arange(N), np.arange(N)
    colors = np.arange(N+1)
    fig, ax = plt.subplots()
    with pytest.raises(ValueError):
        ax.scatter(x, y, c=colors)
    c_rgb = (0.5, 0.5, 0.5)
    ax.scatter(x, y, c=c_rgb)
    ax.scatter(x, y, c=[c_rgb] * N)


def test_scatter_color_masking():
    x = np.array([1, 2, 3])
    y = np.array([1, np.nan, 3])
    colors = np.array(['k', 'w', 'k'])
    linewidths = np.array([1, 2, 3])
    s = plt.scatter(x, y, color=colors, linewidths=linewidths)

    facecolors = s.get_facecolors()
    linecolors = s.get_edgecolors()
    linewidths = s.get_linewidths()
    assert_array_equal(facecolors[1], np.array([0, 0, 0, 1]))
    assert_array_equal(linecolors[1], np.array([0, 0, 0, 1]))
    assert linewidths[1] == 3


def test_eventplot_legend():
    plt.eventplot([1.0], label='Label')
    plt.legend()


def test_bar_broadcast_args():
    fig, ax = plt.subplots()
    # Check that a bar chart with a single height for all bars works.
    ax.bar(range(4), 1)
    # Check that a horizontal chart with one width works.
    ax.bar(0, 1, bottom=range(4), width=1, orientation='horizontal')
    # Check that edgecolor gets broadcasted.
    rect1, rect2 = ax.bar([0, 1], [0, 1], edgecolor=(.1, .2, .3, .4))
    assert rect1.get_edgecolor() == rect2.get_edgecolor() == (.1, .2, .3, .4)


def test_invalid_axis_limits():
    plt.plot([0, 1], [0, 1])
    with pytest.raises(ValueError):
        plt.xlim(np.nan)
    with pytest.raises(ValueError):
        plt.xlim(np.inf)
    with pytest.raises(ValueError):
        plt.ylim(np.nan)
    with pytest.raises(ValueError):
        plt.ylim(np.inf)


# Test all 4 combinations of logs/symlogs for minorticks_on()
@pytest.mark.parametrize('xscale', ['symlog', 'log'])
@pytest.mark.parametrize('yscale', ['symlog', 'log'])
def test_minorticks_on(xscale, yscale):
    ax = plt.subplot(111)
    ax.plot([1, 2, 3, 4])
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.minorticks_on()


def test_twinx_knows_limits():
    fig, ax = plt.subplots()

    ax.axvspan(1, 2)
    xtwin = ax.twinx()
    xtwin.plot([0, 0.5], [1, 2])
    # control axis
    fig2, ax2 = plt.subplots()

    ax2.axvspan(1, 2)
    ax2.plot([0, 0.5], [1, 2])

    assert_array_equal(xtwin.viewLim.intervalx, ax2.viewLim.intervalx)


@pytest.mark.style('mpl20')
@pytest.mark.parametrize('args, kwargs, warning_count',
                         [((1, 1), {'width': 1, 'bottom': 1}, 0),
                          ((1, ), {'height': 1, 'bottom': 1}, 0),
                          ((), {'x': 1, 'height': 1}, 0),
                          ((), {'left': 1, 'height': 1}, 1)])
def test_bar_signature(args, kwargs, warning_count):
    fig, ax = plt.subplots()
    with warnings.catch_warnings(record=True) as w:
        r, = ax.bar(*args, **kwargs)

    assert r.get_width() == kwargs.get('width', 0.8)
    assert r.get_y() == kwargs.get('bottom', 0)
    assert len(w) == warning_count


@pytest.mark.style('mpl20')
@pytest.mark.parametrize('args, kwargs, warning_count',
                         [((1, 1), {'height': 1, 'left': 1}, 0),
                          ((1, ), {'width': 1, 'left': 1}, 0),
                          ((), {'y': 1, 'width': 1}, 0),
                          ((), {'bottom': 1, 'width': 1}, 1)])
def test_barh_signature(args, kwargs, warning_count):
    fig, ax = plt.subplots()
    with warnings.catch_warnings(record=True) as w:
        r, = ax.barh(*args, **kwargs)

    assert r.get_height() == kwargs.get('height', 0.8)
    assert r.get_x() == kwargs.get('left', 0)
    assert len(w) == warning_count


def test_zero_linewidth():
    # Check that setting a zero linewidth doesn't error
    plt.plot([0, 1], [0, 1], ls='--', lw=0)


def test_patch_deprecations():
    fig, ax = plt.subplots()
    with warnings.catch_warnings(record=True) as w:
        assert ax.patch == ax.axesPatch
        assert fig.patch == fig.figurePatch

    assert len(w) == 2


def test_polar_gridlines():
    fig = plt.figure()
    ax = fig.add_subplot(111, polar=True)

    # make all major grid lines lighter, only x grid lines set in 2.1.0
    ax.grid(alpha=0.2)

    # hide y tick labels, no effect in 2.1.0
    plt.setp(ax.yaxis.get_ticklabels(), visible=False)

    fig.canvas.draw()

    assert ax.xaxis.majorTicks[0].gridline.get_alpha() == .2
    assert ax.yaxis.majorTicks[0].gridline.get_alpha() == .2


def test_empty_errorbar_legend():
    fig, ax = plt.subplots()
    ax.errorbar([], [], xerr=[], label='empty y')
    ax.errorbar([], [], yerr=[], label='empty x')
    ax.legend()


def test_plot_columns_cycle_deprecation():
    with pytest.warns(MatplotlibDeprecationWarning):
        plt.plot(np.zeros((2, 2)), np.zeros((2, 3)))
