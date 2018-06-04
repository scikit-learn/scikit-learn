from __future__ import absolute_import, division, print_function

import io

import numpy as np
from numpy.testing import assert_array_almost_equal
import pytest

from matplotlib import (
    collections, path, pyplot as plt, transforms as mtransforms, rcParams)
from matplotlib.image import imread
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.testing.decorators import image_comparison


def test_repeated_save_with_alpha():
    # We want an image which has a background color of bluish green, with an
    # alpha of 0.25.

    fig = Figure([1, 0.4])
    canvas = FigureCanvas(fig)
    fig.set_facecolor((0, 1, 0.4))
    fig.patch.set_alpha(0.25)

    # The target color is fig.patch.get_facecolor()

    buf = io.BytesIO()

    fig.savefig(buf,
                facecolor=fig.get_facecolor(),
                edgecolor='none')

    # Save the figure again to check that the
    # colors don't bleed from the previous renderer.
    buf.seek(0)
    fig.savefig(buf,
                facecolor=fig.get_facecolor(),
                edgecolor='none')

    # Check the first pixel has the desired color & alpha
    # (approx: 0, 1.0, 0.4, 0.25)
    buf.seek(0)
    assert_array_almost_equal(tuple(imread(buf)[0, 0]),
                              (0.0, 1.0, 0.4, 0.250),
                              decimal=3)


def test_large_single_path_collection():
    buff = io.BytesIO()

    # Generates a too-large single path in a path collection that
    # would cause a segfault if the draw_markers optimization is
    # applied.
    f, ax = plt.subplots()
    collection = collections.PathCollection(
        [path.Path([[-10, 5], [10, 5], [10, -5], [-10, -5], [-10, 5]])])
    ax.add_artist(collection)
    ax.set_xlim(10**-3, 1)
    plt.savefig(buff)


def test_marker_with_nan():
    # This creates a marker with nans in it, which was segfaulting the
    # Agg backend (see #3722)
    fig, ax = plt.subplots(1)
    steps = 1000
    data = np.arange(steps)
    ax.semilogx(data)
    ax.fill_between(data, data*0.8, data*1.2)
    buf = io.BytesIO()
    fig.savefig(buf, format='png')


def test_long_path():
    buff = io.BytesIO()

    fig, ax = plt.subplots()
    np.random.seed(0)
    points = np.random.rand(70000)
    ax.plot(points)
    fig.savefig(buff, format='png')


@image_comparison(baseline_images=['agg_filter'],
                  extensions=['png'], remove_text=True)
def test_agg_filter():
    def smooth1d(x, window_len):
        s = np.r_[2*x[0] - x[window_len:1:-1],
                  x,
                  2*x[-1] - x[-1:-window_len:-1]]
        w = np.hanning(window_len)
        y = np.convolve(w/w.sum(), s, mode='same')
        return y[window_len-1:-window_len+1]

    def smooth2d(A, sigma=3):
        window_len = max(int(sigma), 3)*2 + 1
        A1 = np.array([smooth1d(x, window_len) for x in np.asarray(A)])
        A2 = np.transpose(A1)
        A3 = np.array([smooth1d(x, window_len) for x in A2])
        A4 = np.transpose(A3)

        return A4

    class BaseFilter(object):
        def prepare_image(self, src_image, dpi, pad):
            ny, nx, depth = src_image.shape
            padded_src = np.zeros([pad*2 + ny, pad*2 + nx, depth], dtype="d")
            padded_src[pad:-pad, pad:-pad, :] = src_image[:, :, :]

            return padded_src  # , tgt_image

        def get_pad(self, dpi):
            return 0

        def __call__(self, im, dpi):
            pad = self.get_pad(dpi)
            padded_src = self.prepare_image(im, dpi, pad)
            tgt_image = self.process_image(padded_src, dpi)
            return tgt_image, -pad, -pad

    class OffsetFilter(BaseFilter):
        def __init__(self, offsets=None):
            if offsets is None:
                self.offsets = (0, 0)
            else:
                self.offsets = offsets

        def get_pad(self, dpi):
            return int(max(*self.offsets)/72.*dpi)

        def process_image(self, padded_src, dpi):
            ox, oy = self.offsets
            a1 = np.roll(padded_src, int(ox/72.*dpi), axis=1)
            a2 = np.roll(a1, -int(oy/72.*dpi), axis=0)
            return a2

    class GaussianFilter(BaseFilter):
        "simple gauss filter"

        def __init__(self, sigma, alpha=0.5, color=None):
            self.sigma = sigma
            self.alpha = alpha
            if color is None:
                self.color = (0, 0, 0)
            else:
                self.color = color

        def get_pad(self, dpi):
            return int(self.sigma*3/72.*dpi)

        def process_image(self, padded_src, dpi):
            tgt_image = np.zeros_like(padded_src)
            aa = smooth2d(padded_src[:, :, -1]*self.alpha,
                          self.sigma/72.*dpi)
            tgt_image[:, :, -1] = aa
            tgt_image[:, :, :-1] = self.color
            return tgt_image

    class DropShadowFilter(BaseFilter):
        def __init__(self, sigma, alpha=0.3, color=None, offsets=None):
            self.gauss_filter = GaussianFilter(sigma, alpha, color)
            self.offset_filter = OffsetFilter(offsets)

        def get_pad(self, dpi):
            return max(self.gauss_filter.get_pad(dpi),
                       self.offset_filter.get_pad(dpi))

        def process_image(self, padded_src, dpi):
            t1 = self.gauss_filter.process_image(padded_src, dpi)
            t2 = self.offset_filter.process_image(t1, dpi)
            return t2

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # draw lines
    l1, = ax.plot([0.1, 0.5, 0.9], [0.1, 0.9, 0.5], "bo-",
                  mec="b", mfc="w", lw=5, mew=3, ms=10, label="Line 1")
    l2, = ax.plot([0.1, 0.5, 0.9], [0.5, 0.2, 0.7], "ro-",
                  mec="r", mfc="w", lw=5, mew=3, ms=10, label="Line 1")

    gauss = DropShadowFilter(4)

    for l in [l1, l2]:

        # draw shadows with same lines with slight offset.

        xx = l.get_xdata()
        yy = l.get_ydata()
        shadow, = ax.plot(xx, yy)
        shadow.update_from(l)

        # offset transform
        ot = mtransforms.offset_copy(l.get_transform(), ax.figure,
                                     x=4.0, y=-6.0, units='points')

        shadow.set_transform(ot)

        # adjust zorder of the shadow lines so that it is drawn below the
        # original lines
        shadow.set_zorder(l.get_zorder() - 0.5)
        shadow.set_agg_filter(gauss)
        shadow.set_rasterized(True)  # to support mixed-mode renderers

    ax.set_xlim(0., 1.)
    ax.set_ylim(0., 1.)

    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)


def test_too_large_image():
    fig = plt.figure(figsize=(300, 1000))
    buff = io.BytesIO()
    with pytest.raises(ValueError):
        fig.savefig(buff)


def test_chunksize():
    x = range(200)

    # Test without chunksize
    fig, ax = plt.subplots()
    ax.plot(x, np.sin(x))
    fig.canvas.draw()

    # Test with chunksize
    fig, ax = plt.subplots()
    rcParams['agg.path.chunksize'] = 105
    ax.plot(x, np.sin(x))
    fig.canvas.draw()


@pytest.mark.backend('Agg')
def test_jpeg_dpi():
    Image = pytest.importorskip("PIL.Image")
    # Check that dpi is set correctly in jpg files.
    plt.plot([0, 1, 2], [0, 1, 0])
    buf = io.BytesIO()
    plt.savefig(buf, format="jpg", dpi=200)
    im = Image.open(buf)
    assert im.info['dpi'] == (200, 200)
