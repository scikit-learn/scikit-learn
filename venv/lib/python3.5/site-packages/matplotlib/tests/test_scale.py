from __future__ import print_function, unicode_literals

from matplotlib.testing.decorators import image_comparison
import matplotlib.pyplot as plt
from matplotlib.scale import Log10Transform, InvertedLog10Transform
import numpy as np
import io
import pytest


@image_comparison(baseline_images=['log_scales'], remove_text=True)
def test_log_scales():
    ax = plt.figure().add_subplot(122, yscale='log', xscale='symlog')

    ax.axvline(24.1)
    ax.axhline(24.1)


@image_comparison(baseline_images=['logit_scales'], remove_text=True,
                  extensions=['png'])
def test_logit_scales():
    ax = plt.figure().add_subplot(111, xscale='logit')

    # Typical extinction curve for logit
    x = np.array([0.001, 0.003, 0.01, 0.03, 0.1, 0.2, 0.3, 0.4, 0.5,
                  0.6, 0.7, 0.8, 0.9, 0.97, 0.99, 0.997, 0.999])
    y = 1.0 / x

    ax.plot(x, y)
    ax.grid(True)


def test_log_scatter():
    """Issue #1799"""
    fig, ax = plt.subplots(1)

    x = np.arange(10)
    y = np.arange(10) - 1

    ax.scatter(x, y)

    buf = io.BytesIO()
    fig.savefig(buf, format='pdf')

    buf = io.BytesIO()
    fig.savefig(buf, format='eps')

    buf = io.BytesIO()
    fig.savefig(buf, format='svg')


def test_logscale_subs():
    fig, ax = plt.subplots()
    ax.set_yscale('log', subsy=np.array([2, 3, 4]))
    # force draw
    fig.canvas.draw()


@image_comparison(baseline_images=['logscale_mask'], remove_text=True,
                  extensions=['png'])
def test_logscale_mask():
    # Check that zero values are masked correctly on log scales.
    # See github issue 8045
    xs = np.linspace(0, 50, 1001)

    fig, ax = plt.subplots()
    ax.plot(np.exp(-xs**2))
    fig.canvas.draw()
    ax.set(yscale="log")


def test_extra_kwargs_raise():
    fig, ax = plt.subplots()
    with pytest.raises(ValueError):
        ax.set_yscale('log', nonpos='mask')


def test_logscale_invert_transform():
    fig, ax = plt.subplots()
    ax.set_yscale('log')
    # get transformation from data to axes
    tform = (ax.transAxes + ax.transData.inverted()).inverted()

    # direct test of log transform inversion
    assert isinstance(Log10Transform().inverted(), InvertedLog10Transform)


def test_logscale_transform_repr():
    # check that repr of log transform succeeds
    fig, ax = plt.subplots()
    ax.set_yscale('log')
    s = repr(ax.transData)

    # check that repr of log transform returns correct string
    s = repr(Log10Transform(nonpos='clip'))
    assert s == "Log10Transform({!r})".format('clip')


@image_comparison(baseline_images=['logscale_nonpos_values'], remove_text=True,
                  extensions=['png'], style='mpl20')
def test_logscale_nonpos_values():
    np.random.seed(19680801)
    xs = np.random.normal(size=int(1e3))
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    ax1.hist(xs, range=(-5, 5), bins=10)
    ax1.set_yscale('log')
    ax2.hist(xs, range=(-5, 5), bins=10)
    ax2.set_yscale('log', nonposy='mask')

    xdata = np.arange(0, 10, 0.01)
    ydata = np.exp(-xdata)
    edata = 0.2*(10-xdata)*np.cos(5*xdata)*np.exp(-xdata)

    ax3.fill_between(xdata, ydata - edata, ydata + edata)
    ax3.set_yscale('log')

    x = np.logspace(-1, 1)
    y = x ** 3
    yerr = x**2
    ax4.errorbar(x, y, yerr=yerr)

    ax4.set_yscale('log')
    ax4.set_xscale('log')
