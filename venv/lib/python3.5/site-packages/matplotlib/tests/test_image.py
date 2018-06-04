from __future__ import absolute_import, division, print_function

import six

from copy import copy
import io
import os
import sys
import warnings

import numpy as np
from numpy import ma
from numpy.testing import assert_array_equal

from matplotlib import (
    colors, image as mimage, patches, pyplot as plt,
    rc_context, rcParams)
from matplotlib.image import (AxesImage, BboxImage, FigureImage,
                              NonUniformImage, PcolorImage)
from matplotlib.testing.decorators import image_comparison
from matplotlib.transforms import Bbox, Affine2D, TransformedBbox

import pytest


try:
    from PIL import Image
    HAS_PIL = True
except ImportError:
    HAS_PIL = False
needs_pillow = pytest.mark.xfail(not HAS_PIL, reason='Test requires Pillow')


@image_comparison(baseline_images=['image_interps'], style='mpl20')
def test_image_interps():
    'make the basic nearest, bilinear and bicubic interps'
    X = np.arange(100)
    X = X.reshape(5, 20)

    fig = plt.figure()
    ax1 = fig.add_subplot(311)
    ax1.imshow(X, interpolation='nearest')
    ax1.set_title('three interpolations')
    ax1.set_ylabel('nearest')

    ax2 = fig.add_subplot(312)
    ax2.imshow(X, interpolation='bilinear')
    ax2.set_ylabel('bilinear')

    ax3 = fig.add_subplot(313)
    ax3.imshow(X, interpolation='bicubic')
    ax3.set_ylabel('bicubic')


@image_comparison(baseline_images=['interp_nearest_vs_none'],
                  extensions=['pdf', 'svg'], remove_text=True)
def test_interp_nearest_vs_none():
    'Test the effect of "nearest" and "none" interpolation'
    # Setting dpi to something really small makes the difference very
    # visible. This works fine with pdf, since the dpi setting doesn't
    # affect anything but images, but the agg output becomes unusably
    # small.
    rcParams['savefig.dpi'] = 3
    X = np.array([[[218, 165, 32], [122, 103, 238]],
                  [[127, 255, 0], [255, 99, 71]]], dtype=np.uint8)
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax1.imshow(X, interpolation='none')
    ax1.set_title('interpolation none')
    ax2 = fig.add_subplot(122)
    ax2.imshow(X, interpolation='nearest')
    ax2.set_title('interpolation nearest')


def do_figimage(suppressComposite):
    """ Helper for the next two tests """
    fig = plt.figure(figsize=(2,2), dpi=100)
    fig.suppressComposite = suppressComposite
    x,y = np.ix_(np.arange(100.0)/100.0, np.arange(100.0)/100.0)
    z = np.sin(x**2 + y**2 - x*y)
    c = np.sin(20*x**2 + 50*y**2)
    img = z + c/5

    fig.figimage(img, xo=0, yo=0, origin='lower')
    fig.figimage(img[::-1,:], xo=0, yo=100, origin='lower')
    fig.figimage(img[:,::-1], xo=100, yo=0, origin='lower')
    fig.figimage(img[::-1,::-1], xo=100, yo=100, origin='lower')


@image_comparison(baseline_images=['figimage-0'],
                  extensions=['png','pdf'])
def test_figimage0():
    'test the figimage method'

    suppressComposite = False
    do_figimage(suppressComposite)


@image_comparison(baseline_images=['figimage-1'],
                  extensions=['png','pdf'])
def test_figimage1():
    'test the figimage method'
    suppressComposite = True
    do_figimage(suppressComposite)


def test_image_python_io():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot([1,2,3])
    buffer = io.BytesIO()
    fig.savefig(buffer)
    buffer.seek(0)
    plt.imread(buffer)


@needs_pillow
def test_imread_pil_uint16():
    img = plt.imread(os.path.join(os.path.dirname(__file__),
                     'baseline_images', 'test_image', 'uint16.tif'))
    assert (img.dtype == np.uint16)
    assert np.sum(img) == 134184960


@pytest.mark.skipif(sys.version_info < (3, 6), reason="requires Python 3.6+")
@needs_pillow
def test_imread_fspath():
    from pathlib import Path
    img = plt.imread(
        Path(__file__).parent / 'baseline_images/test_image/uint16.tif')
    assert img.dtype == np.uint16
    assert np.sum(img) == 134184960


def test_imsave():
    # The goal here is that the user can specify an output logical DPI
    # for the image, but this will not actually add any extra pixels
    # to the image, it will merely be used for metadata purposes.

    # So we do the traditional case (dpi == 1), and the new case (dpi
    # == 100) and read the resulting PNG files back in and make sure
    # the data is 100% identical.
    np.random.seed(1)
    data = np.random.rand(256, 128)

    buff_dpi1 = io.BytesIO()
    plt.imsave(buff_dpi1, data, dpi=1)

    buff_dpi100 = io.BytesIO()
    plt.imsave(buff_dpi100, data, dpi=100)

    buff_dpi1.seek(0)
    arr_dpi1 = plt.imread(buff_dpi1)

    buff_dpi100.seek(0)
    arr_dpi100 = plt.imread(buff_dpi100)

    assert arr_dpi1.shape == (256, 128, 4)
    assert arr_dpi100.shape == (256, 128, 4)

    assert_array_equal(arr_dpi1, arr_dpi100)


@pytest.mark.skipif(sys.version_info < (3, 6), reason="requires Python 3.6+")
@pytest.mark.parametrize("fmt", ["png", "pdf", "ps", "eps", "svg"])
def test_imsave_fspath(fmt):
    Path = pytest.importorskip("pathlib").Path
    plt.imsave(Path(os.devnull), np.array([[0, 1]]), format=fmt)


def test_imsave_color_alpha():
    # Test that imsave accept arrays with ndim=3 where the third dimension is
    # color and alpha without raising any exceptions, and that the data is
    # acceptably preserved through a save/read roundtrip.
    np.random.seed(1)

    for origin in ['lower', 'upper']:
        data = np.random.rand(16, 16, 4)
        buff = io.BytesIO()
        plt.imsave(buff, data, origin=origin, format="png")

        buff.seek(0)
        arr_buf = plt.imread(buff)

        # Recreate the float -> uint8 conversion of the data
        # We can only expect to be the same with 8 bits of precision,
        # since that's what the PNG file used.
        data = (255*data).astype('uint8')
        if origin == 'lower':
            data = data[::-1]
        arr_buf = (255*arr_buf).astype('uint8')

        assert_array_equal(data, arr_buf)


@image_comparison(baseline_images=['image_alpha'], remove_text=True)
def test_image_alpha():
    plt.figure()

    np.random.seed(0)
    Z = np.random.rand(6, 6)

    plt.subplot(131)
    plt.imshow(Z, alpha=1.0, interpolation='none')

    plt.subplot(132)
    plt.imshow(Z, alpha=0.5, interpolation='none')

    plt.subplot(133)
    plt.imshow(Z, alpha=0.5, interpolation='nearest')

def test_cursor_data():
    from matplotlib.backend_bases import MouseEvent

    fig, ax = plt.subplots()
    im = ax.imshow(np.arange(100).reshape(10, 10), origin='upper')

    x, y = 4, 4
    xdisp, ydisp = ax.transData.transform_point([x, y])

    event = MouseEvent('motion_notify_event', fig.canvas, xdisp, ydisp)
    assert im.get_cursor_data(event) == 44

    # Now try for a point outside the image
    # Tests issue #4957
    x, y = 10.1, 4
    xdisp, ydisp = ax.transData.transform_point([x, y])

    event = MouseEvent('motion_notify_event', fig.canvas, xdisp, ydisp)
    assert im.get_cursor_data(event) is None

    # Hmm, something is wrong here... I get 0, not None...
    # But, this works further down in the tests with extents flipped
    #x, y = 0.1, -0.1
    #xdisp, ydisp = ax.transData.transform_point([x, y])
    #event = MouseEvent('motion_notify_event', fig.canvas, xdisp, ydisp)
    #z = im.get_cursor_data(event)
    #assert z is None, "Did not get None, got %d" % z

    ax.clear()
    # Now try with the extents flipped.
    im = ax.imshow(np.arange(100).reshape(10, 10), origin='lower')

    x, y = 4, 4
    xdisp, ydisp = ax.transData.transform_point([x, y])

    event = MouseEvent('motion_notify_event', fig.canvas, xdisp, ydisp)
    assert im.get_cursor_data(event) == 44

    fig, ax = plt.subplots()
    im = ax.imshow(np.arange(100).reshape(10, 10), extent=[0, 0.5, 0, 0.5])

    x, y = 0.25, 0.25
    xdisp, ydisp = ax.transData.transform_point([x, y])

    event = MouseEvent('motion_notify_event', fig.canvas, xdisp, ydisp)
    assert im.get_cursor_data(event) == 55

    # Now try for a point outside the image
    # Tests issue #4957
    x, y = 0.75, 0.25
    xdisp, ydisp = ax.transData.transform_point([x, y])

    event = MouseEvent('motion_notify_event', fig.canvas, xdisp, ydisp)
    assert im.get_cursor_data(event) is None

    x, y = 0.01, -0.01
    xdisp, ydisp = ax.transData.transform_point([x, y])

    event = MouseEvent('motion_notify_event', fig.canvas, xdisp, ydisp)
    assert im.get_cursor_data(event) is None


@image_comparison(baseline_images=['image_clip'], style='mpl20')
def test_image_clip():
    d = [[1, 2], [3, 4]]

    fig, ax = plt.subplots()
    im = ax.imshow(d)
    patch = patches.Circle((0, 0), radius=1, transform=ax.transData)
    im.set_clip_path(patch)


@image_comparison(baseline_images=['image_cliprect'], style='mpl20')
def test_image_cliprect():
    import matplotlib.patches as patches

    fig = plt.figure()
    ax = fig.add_subplot(111)
    d = [[1,2],[3,4]]

    im = ax.imshow(d, extent=(0,5,0,5))

    rect = patches.Rectangle(xy=(1,1), width=2, height=2, transform=im.axes.transData)
    im.set_clip_path(rect)


@image_comparison(baseline_images=['imshow'], remove_text=True, style='mpl20')
def test_imshow():
    import numpy as np
    import matplotlib.pyplot as plt

    fig = plt.figure()
    arr = np.arange(100).reshape((10, 10))
    ax = fig.add_subplot(111)
    ax.imshow(arr, interpolation="bilinear", extent=(1,2,1,2))
    ax.set_xlim(0,3)
    ax.set_ylim(0,3)


@image_comparison(baseline_images=['no_interpolation_origin'],
                  remove_text=True)
def test_no_interpolation_origin():
    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax.imshow(np.arange(100).reshape((2, 50)), origin="lower",
              interpolation='none')

    ax = fig.add_subplot(212)
    ax.imshow(np.arange(100).reshape((2, 50)), interpolation='none')


@image_comparison(baseline_images=['image_shift'], remove_text=True,
                  extensions=['pdf', 'svg'])
def test_image_shift():
    from matplotlib.colors import LogNorm

    imgData = [[1.0/(x) + 1.0/(y) for x in range(1,100)] for y in range(1,100)]
    tMin=734717.945208
    tMax=734717.946366

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(imgData, norm=LogNorm(), interpolation='none',
              extent=(tMin, tMax, 1, 100))
    ax.set_aspect('auto')


def test_image_edges():
    f = plt.figure(figsize=[1, 1])
    ax = f.add_axes([0, 0, 1, 1], frameon=False)

    data = np.tile(np.arange(12), 15).reshape(20, 9)

    im = ax.imshow(data, origin='upper', extent=[-10, 10, -10, 10],
                   interpolation='none', cmap='gray')

    x = y = 2
    ax.set_xlim([-x, x])
    ax.set_ylim([-y, y])

    ax.set_xticks([])
    ax.set_yticks([])

    buf = io.BytesIO()
    f.savefig(buf, facecolor=(0, 1, 0))

    buf.seek(0)

    im = plt.imread(buf)
    r, g, b, a = sum(im[:, 0])
    r, g, b, a = sum(im[:, -1])

    assert g != 100, 'Expected a non-green edge - but sadly, it was.'


@image_comparison(baseline_images=['image_composite_background'],
                  remove_text=True,
                  style='mpl20')
def test_image_composite_background():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    arr = np.arange(12).reshape(4, 3)
    ax.imshow(arr, extent=[0, 2, 15, 0])
    ax.imshow(arr, extent=[4, 6, 15, 0])
    ax.set_facecolor((1, 0, 0, 0.5))
    ax.set_xlim([0, 12])


@image_comparison(baseline_images=['image_composite_alpha'],
                  remove_text=True)
def test_image_composite_alpha():
    """
    Tests that the alpha value is recognized and correctly applied in the
    process of compositing images together.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    arr = np.zeros((11, 21, 4))
    arr[:, :, 0] = 1
    arr[:, :, 3] = np.concatenate((np.arange(0, 1.1, 0.1), np.arange(0, 1, 0.1)[::-1]))
    arr2 = np.zeros((21, 11, 4))
    arr2[:, :, 0] = 1
    arr2[:, :, 1] = 1
    arr2[:, :, 3] = np.concatenate((np.arange(0, 1.1, 0.1), np.arange(0, 1, 0.1)[::-1]))[:, np.newaxis]
    ax.imshow(arr, extent=[1, 2, 5, 0], alpha=0.3)
    ax.imshow(arr, extent=[2, 3, 5, 0], alpha=0.6)
    ax.imshow(arr, extent=[3, 4, 5, 0])
    ax.imshow(arr2, extent=[0, 5, 1, 2])
    ax.imshow(arr2, extent=[0, 5, 2, 3], alpha=0.6)
    ax.imshow(arr2, extent=[0, 5, 3, 4], alpha=0.3)
    ax.set_facecolor((0, 0.5, 0, 1))
    ax.set_xlim([0, 5])
    ax.set_ylim([5, 0])


@image_comparison(baseline_images=['rasterize_10dpi'],
                  extensions=['pdf', 'svg'],
                  remove_text=True, style='mpl20')
def test_rasterize_dpi():
    # This test should check rasterized rendering with high output resolution.
    # It plots a rasterized line and a normal image with implot. So it will catch
    # when images end up in the wrong place in case of non-standard dpi setting.
    # Instead of high-res rasterization i use low-res.  Therefore the fact that the
    # resolution is non-standard is easily checked by image_comparison.
    import numpy as np
    import matplotlib.pyplot as plt

    img = np.asarray([[1, 2], [3, 4]])

    fig, axes = plt.subplots(1, 3, figsize = (3, 1))

    axes[0].imshow(img)

    axes[1].plot([0,1],[0,1], linewidth=20., rasterized=True)
    axes[1].set(xlim = (0,1), ylim = (-1, 2))

    axes[2].plot([0,1],[0,1], linewidth=20.)
    axes[2].set(xlim = (0,1), ylim = (-1, 2))

    # Low-dpi PDF rasterization errors prevent proper image comparison tests.
    # Hide detailed structures like the axes spines.
    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

    rcParams['savefig.dpi'] = 10


@image_comparison(baseline_images=['bbox_image_inverted'], remove_text=True,
                  style='mpl20')
def test_bbox_image_inverted():
    # This is just used to produce an image to feed to BboxImage
    image = np.arange(100).reshape((10, 10))

    ax = plt.subplot(111)
    bbox_im = BboxImage(
        TransformedBbox(Bbox([[100, 100], [0, 0]]), ax.transData))
    bbox_im.set_data(image)
    bbox_im.set_clip_on(False)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.add_artist(bbox_im)

    image = np.identity(10)

    bbox_im = BboxImage(
        TransformedBbox(Bbox([[0.1, 0.2], [0.3, 0.25]]), ax.figure.transFigure))
    bbox_im.set_data(image)
    bbox_im.set_clip_on(False)
    ax.add_artist(bbox_im)


def test_get_window_extent_for_AxisImage():
    # Create a figure of known size (1000x1000 pixels), place an image
    # object at a given location and check that get_window_extent()
    # returns the correct bounding box values (in pixels).

    im = np.array([[0.25, 0.75, 1.0, 0.75], [0.1, 0.65, 0.5, 0.4],
                   [0.6, 0.3, 0.0, 0.2], [0.7, 0.9, 0.4, 0.6]])
    fig = plt.figure(figsize=(10, 10), dpi=100)
    ax = plt.subplot()
    ax.set_position([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    im_obj = ax.imshow(im, extent=[0.4, 0.7, 0.2, 0.9], interpolation='nearest')

    fig.canvas.draw()
    renderer = fig.canvas.renderer
    im_bbox = im_obj.get_window_extent(renderer)

    assert_array_equal(im_bbox.get_points(), [[400, 200], [700, 900]])


@image_comparison(baseline_images=['zoom_and_clip_upper_origin'],
                  remove_text=True,
                  extensions=['png'],
                  style='mpl20')
def test_zoom_and_clip_upper_origin():
    image = np.arange(100)
    image = image.reshape((10, 10))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(image)
    ax.set_ylim(2.0, -0.5)
    ax.set_xlim(-0.5, 2.0)


def test_nonuniformimage_setcmap():
    ax = plt.gca()
    im = NonUniformImage(ax)
    im.set_cmap('Blues')


def test_nonuniformimage_setnorm():
    ax = plt.gca()
    im = NonUniformImage(ax)
    im.set_norm(plt.Normalize())


@needs_pillow
def test_jpeg_2d():
    # smoke test that mode-L pillow images work.
    imd = np.ones((10, 10), dtype='uint8')
    for i in range(10):
        imd[i, :] = np.linspace(0.0, 1.0, 10) * 255
    im = Image.new('L', (10, 10))
    im.putdata(imd.flatten())
    fig, ax = plt.subplots()
    ax.imshow(im)


@needs_pillow
def test_jpeg_alpha():
    plt.figure(figsize=(1, 1), dpi=300)
    # Create an image that is all black, with a gradient from 0-1 in
    # the alpha channel from left to right.
    im = np.zeros((300, 300, 4), dtype=float)
    im[..., 3] = np.linspace(0.0, 1.0, 300)

    plt.figimage(im)

    buff = io.BytesIO()
    with rc_context({'savefig.facecolor': 'red'}):
        plt.savefig(buff, transparent=True, format='jpg', dpi=300)

    buff.seek(0)
    image = Image.open(buff)

    # If this fails, there will be only one color (all black). If this
    # is working, we should have all 256 shades of grey represented.
    num_colors = len(image.getcolors(256))
    assert 175 <= num_colors <= 185
    # The fully transparent part should be red.
    corner_pixel = image.getpixel((0, 0))
    assert corner_pixel == (254, 0, 0)


def test_nonuniformimage_setdata():
    ax = plt.gca()
    im = NonUniformImage(ax)
    x = np.arange(3, dtype=float)
    y = np.arange(4, dtype=float)
    z = np.arange(12, dtype=float).reshape((4, 3))
    im.set_data(x, y, z)
    x[0] = y[0] = z[0, 0] = 9.9
    assert im._A[0, 0] == im._Ax[0] == im._Ay[0] == 0, 'value changed'


def test_axesimage_setdata():
    ax = plt.gca()
    im = AxesImage(ax)
    z = np.arange(12, dtype=float).reshape((4, 3))
    im.set_data(z)
    z[0, 0] = 9.9
    assert im._A[0, 0] == 0, 'value changed'


def test_figureimage_setdata():
    fig = plt.gcf()
    im = FigureImage(fig)
    z = np.arange(12, dtype=float).reshape((4, 3))
    im.set_data(z)
    z[0, 0] = 9.9
    assert im._A[0, 0] == 0, 'value changed'


def test_pcolorimage_setdata():
    ax = plt.gca()
    im = PcolorImage(ax)
    x = np.arange(3, dtype=float)
    y = np.arange(4, dtype=float)
    z = np.arange(6, dtype=float).reshape((3, 2))
    im.set_data(x, y, z)
    x[0] = y[0] = z[0, 0] = 9.9
    assert im._A[0, 0] == im._Ax[0] == im._Ay[0] == 0, 'value changed'


def test_pcolorimage_extent():
    im = plt.hist2d([1, 2, 3], [3, 5, 6],
                    bins=[[0, 3, 7], [1, 2, 3]])[-1]
    assert im.get_extent() == (0, 7, 1, 3)


def test_minimized_rasterized():
    # This ensures that the rasterized content in the colorbars is
    # only as thick as the colorbar, and doesn't extend to other parts
    # of the image.  See #5814.  While the original bug exists only
    # in Postscript, the best way to detect it is to generate SVG
    # and then parse the output to make sure the two colorbar images
    # are the same size.
    from xml.etree import ElementTree

    np.random.seed(0)
    data = np.random.rand(10, 10)

    fig, ax = plt.subplots(1, 2)
    p1 = ax[0].pcolormesh(data)
    p2 = ax[1].pcolormesh(data)

    plt.colorbar(p1, ax=ax[0])
    plt.colorbar(p2, ax=ax[1])

    buff = io.BytesIO()
    plt.savefig(buff, format='svg')

    buff = io.BytesIO(buff.getvalue())
    tree = ElementTree.parse(buff)
    width = None
    for image in tree.iter('image'):
        if width is None:
            width = image['width']
        else:
            if image['width'] != width:
                assert False


@pytest.mark.network
def test_load_from_url():
    req = six.moves.urllib.request.urlopen(
        "http://matplotlib.org/_static/logo_sidebar_horiz.png")
    plt.imread(req)


@image_comparison(baseline_images=['log_scale_image'],
                  remove_text=True)
# The recwarn fixture captures a warning in image_comparison.
def test_log_scale_image(recwarn):
    Z = np.zeros((10, 10))
    Z[::2] = 1

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.imshow(Z, extent=[1, 100, 1, 100], cmap='viridis',
              vmax=1, vmin=-1)
    ax.set_yscale('log')


@image_comparison(baseline_images=['rotate_image'],
                  remove_text=True)
def test_rotate_image():
    delta = 0.25
    x = y = np.arange(-3.0, 3.0, delta)
    X, Y = np.meshgrid(x, y)
    Z1 = np.exp(-(X**2 + Y**2) / 2) / (2 * np.pi)
    Z2 = (np.exp(-(((X - 1) / 1.5)**2 + ((Y - 1) / 0.5)**2) / 2) /
          (2 * np.pi * 0.5 * 1.5))
    Z = Z2 - Z1  # difference of Gaussians

    fig, ax1 = plt.subplots(1, 1)
    im1 = ax1.imshow(Z, interpolation='none', cmap='viridis',
                     origin='lower',
                     extent=[-2, 4, -3, 2], clip_on=True)

    trans_data2 = Affine2D().rotate_deg(30) + ax1.transData
    im1.set_transform(trans_data2)

    # display intended extent of the image
    x1, x2, y1, y2 = im1.get_extent()

    ax1.plot([x1, x2, x2, x1, x1], [y1, y1, y2, y2, y1], "r--", lw=3,
             transform=trans_data2)

    ax1.set_xlim(2, 5)
    ax1.set_ylim(0, 4)


def test_image_preserve_size():
    buff = io.BytesIO()

    im = np.zeros((481, 321))
    plt.imsave(buff, im, format="png")

    buff.seek(0)
    img = plt.imread(buff)

    assert img.shape[:2] == im.shape


def test_image_preserve_size2():
    n = 7
    data = np.identity(n, float)

    fig = plt.figure(figsize=(n, n), frameon=False)

    ax = plt.Axes(fig, [0.0, 0.0, 1.0, 1.0])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(data, interpolation='nearest', origin='lower',aspect='auto')
    buff = io.BytesIO()
    fig.savefig(buff, dpi=1)

    buff.seek(0)
    img = plt.imread(buff)

    assert img.shape == (7, 7, 4)

    assert_array_equal(np.asarray(img[:, :, 0], bool),
                       np.identity(n, bool)[::-1])


@image_comparison(baseline_images=['mask_image_over_under'],
                  remove_text=True, extensions=['png'])
def test_mask_image_over_under():
    delta = 0.025
    x = y = np.arange(-3.0, 3.0, delta)
    X, Y = np.meshgrid(x, y)
    Z1 = np.exp(-(X**2 + Y**2) / 2) / (2 * np.pi)
    Z2 = (np.exp(-(((X - 1) / 1.5)**2 + ((Y - 1) / 0.5)**2) / 2) /
          (2 * np.pi * 0.5 * 1.5))
    Z = 10*(Z2 - Z1)  # difference of Gaussians

    palette = copy(plt.cm.gray)
    palette.set_over('r', 1.0)
    palette.set_under('g', 1.0)
    palette.set_bad('b', 1.0)
    Zm = ma.masked_where(Z > 1.2, Z)
    fig, (ax1, ax2) = plt.subplots(1, 2)
    im = ax1.imshow(Zm, interpolation='bilinear',
                    cmap=palette,
                    norm=colors.Normalize(vmin=-1.0, vmax=1.0, clip=False),
                    origin='lower', extent=[-3, 3, -3, 3])
    ax1.set_title('Green=low, Red=high, Blue=bad')
    fig.colorbar(im, extend='both', orientation='horizontal',
                 ax=ax1, aspect=10)

    im = ax2.imshow(Zm, interpolation='nearest',
                    cmap=palette,
                    norm=colors.BoundaryNorm([-1, -0.5, -0.2, 0, 0.2, 0.5, 1],
                                             ncolors=256, clip=False),
                    origin='lower', extent=[-3, 3, -3, 3])
    ax2.set_title('With BoundaryNorm')
    fig.colorbar(im, extend='both', spacing='proportional',
                 orientation='horizontal', ax=ax2, aspect=10)


@image_comparison(baseline_images=['mask_image'],
                  remove_text=True)
def test_mask_image():
    # Test mask image two ways: Using nans and using a masked array.

    fig, (ax1, ax2) = plt.subplots(1, 2)

    A = np.ones((5, 5))
    A[1:2, 1:2] = np.nan

    ax1.imshow(A, interpolation='nearest')

    A = np.zeros((5, 5), dtype=bool)
    A[1:2, 1:2] = True
    A = np.ma.masked_array(np.ones((5, 5), dtype=np.uint16), A)

    ax2.imshow(A, interpolation='nearest')


@image_comparison(baseline_images=['imshow_endianess'],
                  remove_text=True, extensions=['png'])
def test_imshow_endianess():
    x = np.arange(10)
    X, Y = np.meshgrid(x, x)
    Z = ((X-5)**2 + (Y-5)**2)**0.5

    fig, (ax1, ax2) = plt.subplots(1, 2)

    kwargs = dict(origin="lower", interpolation='nearest',
                  cmap='viridis')

    ax1.imshow(Z.astype('<f8'), **kwargs)
    ax2.imshow(Z.astype('>f8'), **kwargs)


@image_comparison(baseline_images=['imshow_masked_interpolation'],
                  remove_text=True, style='mpl20')
def test_imshow_masked_interpolation():

    cm = copy(plt.get_cmap('viridis'))
    cm.set_over('r')
    cm.set_under('b')
    cm.set_bad('k')

    N = 20
    n = colors.Normalize(vmin=0, vmax=N*N-1)

    # data = np.random.random((N, N))*N*N
    data = np.arange(N*N, dtype='float').reshape(N, N)

    data[5, 5] = -1
    # This will cause crazy ringing for the higher-order
    # interpolations
    data[15, 5] = 1e5

    # data[3, 3] = np.nan

    data[15, 15] = np.inf

    mask = np.zeros_like(data).astype('bool')
    mask[5, 15] = True

    data = np.ma.masked_array(data, mask)

    fig, ax_grid = plt.subplots(3, 6)

    for interp, ax in zip(sorted(mimage._interpd_), ax_grid.ravel()):
        ax.set_title(interp)
        ax.imshow(data, norm=n, cmap=cm, interpolation=interp)
        ax.axis('off')


def test_imshow_no_warn_invalid():
    with warnings.catch_warnings(record=True) as warns:
        warnings.simplefilter("always")
        plt.imshow([[1, 2], [3, np.nan]])
    assert len(warns) == 0


@pytest.mark.parametrize(
    'dtype', [np.dtype(s) for s in 'u2 u4 i2 i4 i8 f4 f8'.split()])
def test_imshow_clips_rgb_to_valid_range(dtype):
    arr = np.arange(300, dtype=dtype).reshape((10, 10, 3))
    if dtype.kind != 'u':
        arr -= 10
    too_low = arr < 0
    too_high = arr > 255
    if dtype.kind == 'f':
        arr = arr / 255
    _, ax = plt.subplots()
    out = ax.imshow(arr).get_array()
    assert (out[too_low] == 0).all()
    if dtype.kind == 'f':
        assert (out[too_high] == 1).all()
        assert out.dtype.kind == 'f'
    else:
        assert (out[too_high] == 255).all()
        assert out.dtype == np.uint8


@image_comparison(baseline_images=['imshow_flatfield'],
                  remove_text=True, style='mpl20',
                  extensions=['png'])
def test_imshow_flatfield():
    fig, ax = plt.subplots()
    im = ax.imshow(np.ones((5, 5)))
    im.set_clim(.5, 1.5)


@image_comparison(baseline_images=['imshow_bignumbers'],
                  remove_text=True, style='mpl20',
                  extensions=['png'])
def test_imshow_bignumbers():
    # putting a big number in an array of integers shouldn't
    # ruin the dynamic range of the resolved bits.
    fig, ax = plt.subplots()
    img = np.array([[1, 2, 1e12],[3, 1, 4]], dtype=np.uint64)
    pc = ax.imshow(img)
    pc.set_clim(0, 5)


@image_comparison(baseline_images=['imshow_bignumbers_real'],
                  remove_text=True, style='mpl20',
                  extensions=['png'])
def test_imshow_bignumbers_real():
    # putting a big number in an array of integers shouldn't
    # ruin the dynamic range of the resolved bits.
    fig, ax = plt.subplots()
    img = np.array([[2., 1., 1.e22],[4., 1., 3.]])
    pc = ax.imshow(img)
    pc.set_clim(0, 5)


@pytest.mark.parametrize(
    "make_norm",
    [colors.Normalize,
     colors.LogNorm,
     lambda: colors.SymLogNorm(1),
     lambda: colors.PowerNorm(1)])
def test_empty_imshow(make_norm):
    fig, ax = plt.subplots()
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", "Attempting to set identical left==right")
        im = ax.imshow([[]], norm=make_norm())
    im.set_extent([-5, 5, -5, 5])
    fig.canvas.draw()

    with pytest.raises(RuntimeError):
        im.make_image(fig._cachedRenderer)


def test_imshow_float128():
    fig, ax = plt.subplots()
    ax.imshow(np.zeros((3, 3), dtype=np.longdouble))


def test_imshow_bool():
    fig, ax = plt.subplots()
    ax.imshow(np.array([[True, False], [False, True]], dtype=bool))


def test_imshow_deprecated_interd_warn():
    im = plt.imshow([[1, 2], [3, np.nan]])
    for k in ('_interpd', '_interpdr', 'iterpnames'):
        with warnings.catch_warnings(record=True) as warns:
            getattr(im, k)
        assert len(warns) == 1


def test_full_invalid():
    x = np.ones((10, 10))
    x[:] = np.nan

    f, ax = plt.subplots()
    ax.imshow(x)

    f.canvas.draw()


@pytest.mark.parametrize("fmt,counted",
                         [("ps", b" colorimage"), ("svg", b"<image")])
@pytest.mark.parametrize("composite_image,count", [(True, 1), (False, 2)])
def test_composite(fmt, counted, composite_image, count):
    # Test that figures can be saved with and without combining multiple images
    # (on a single set of axes) into a single composite image.
    X, Y = np.meshgrid(np.arange(-5, 5, 1), np.arange(-5, 5, 1))
    Z = np.sin(Y ** 2)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlim(0, 3)
    ax.imshow(Z, extent=[0, 1, 0, 1])
    ax.imshow(Z[::-1], extent=[2, 3, 0, 1])
    plt.rcParams['image.composite_image'] = composite_image
    buf = io.BytesIO()
    fig.savefig(buf, format=fmt)
    assert buf.getvalue().count(counted) == count


def test_relim():
    fig, ax = plt.subplots()
    ax.imshow([[0]], extent=(0, 1, 0, 1))
    ax.relim()
    ax.autoscale()
    assert ax.get_xlim() == ax.get_ylim() == (0, 1)
