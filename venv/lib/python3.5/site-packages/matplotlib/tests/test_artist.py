from __future__ import absolute_import, division, print_function

import io
import warnings
from itertools import chain

import numpy as np

import pytest

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.path as mpath
import matplotlib.transforms as mtransforms
import matplotlib.collections as mcollections
import matplotlib.artist as martist
from matplotlib.testing.decorators import image_comparison


def test_patch_transform_of_none():
    # tests the behaviour of patches added to an Axes with various transform
    # specifications

    ax = plt.axes()
    ax.set_xlim([1, 3])
    ax.set_ylim([1, 3])

    # Draw an ellipse over data coord (2,2) by specifying device coords.
    xy_data = (2, 2)
    xy_pix = ax.transData.transform_point(xy_data)

    # Not providing a transform of None puts the ellipse in data coordinates .
    e = mpatches.Ellipse(xy_data, width=1, height=1, fc='yellow', alpha=0.5)
    ax.add_patch(e)
    assert e._transform == ax.transData

    # Providing a transform of None puts the ellipse in device coordinates.
    e = mpatches.Ellipse(xy_pix, width=120, height=120, fc='coral',
                         transform=None, alpha=0.5)
    assert e.is_transform_set() is True
    ax.add_patch(e)
    assert isinstance(e._transform, mtransforms.IdentityTransform)

    # Providing an IdentityTransform puts the ellipse in device coordinates.
    e = mpatches.Ellipse(xy_pix, width=100, height=100,
                         transform=mtransforms.IdentityTransform(), alpha=0.5)
    ax.add_patch(e)
    assert isinstance(e._transform, mtransforms.IdentityTransform)

    # Not providing a transform, and then subsequently "get_transform" should
    # not mean that "is_transform_set".
    e = mpatches.Ellipse(xy_pix, width=120, height=120, fc='coral',
                         alpha=0.5)
    intermediate_transform = e.get_transform()
    assert e.is_transform_set() is False
    ax.add_patch(e)
    assert e.get_transform() != intermediate_transform
    assert e.is_transform_set() is True
    assert e._transform == ax.transData


def test_collection_transform_of_none():
    # tests the behaviour of collections added to an Axes with various
    # transform specifications

    ax = plt.axes()
    ax.set_xlim([1, 3])
    ax.set_ylim([1, 3])

    # draw an ellipse over data coord (2,2) by specifying device coords
    xy_data = (2, 2)
    xy_pix = ax.transData.transform_point(xy_data)

    # not providing a transform of None puts the ellipse in data coordinates
    e = mpatches.Ellipse(xy_data, width=1, height=1)
    c = mcollections.PatchCollection([e], facecolor='yellow', alpha=0.5)
    ax.add_collection(c)
    # the collection should be in data coordinates
    assert c.get_offset_transform() + c.get_transform() == ax.transData

    # providing a transform of None puts the ellipse in device coordinates
    e = mpatches.Ellipse(xy_pix, width=120, height=120)
    c = mcollections.PatchCollection([e], facecolor='coral',
                                     alpha=0.5)
    c.set_transform(None)
    ax.add_collection(c)
    assert isinstance(c.get_transform(), mtransforms.IdentityTransform)

    # providing an IdentityTransform puts the ellipse in device coordinates
    e = mpatches.Ellipse(xy_pix, width=100, height=100)
    c = mcollections.PatchCollection([e],
                                 transform=mtransforms.IdentityTransform(),
                                 alpha=0.5)
    ax.add_collection(c)
    assert isinstance(c._transOffset, mtransforms.IdentityTransform)


@image_comparison(baseline_images=["clip_path_clipping"], remove_text=True)
def test_clipping():
    exterior = mpath.Path.unit_rectangle().deepcopy()
    exterior.vertices *= 4
    exterior.vertices -= 2
    interior = mpath.Path.unit_circle().deepcopy()
    interior.vertices = interior.vertices[::-1]
    clip_path = mpath.Path(vertices=np.concatenate([exterior.vertices,
                                                    interior.vertices]),
                           codes=np.concatenate([exterior.codes,
                                                 interior.codes]))

    star = mpath.Path.unit_regular_star(6).deepcopy()
    star.vertices *= 2.6

    ax1 = plt.subplot(121)
    col = mcollections.PathCollection([star], lw=5, edgecolor='blue',
                                      facecolor='red', alpha=0.7, hatch='*')
    col.set_clip_path(clip_path, ax1.transData)
    ax1.add_collection(col)

    ax2 = plt.subplot(122, sharex=ax1, sharey=ax1)
    patch = mpatches.PathPatch(star, lw=5, edgecolor='blue', facecolor='red',
                               alpha=0.7, hatch='*')
    patch.set_clip_path(clip_path, ax2.transData)
    ax2.add_patch(patch)

    ax1.set_xlim([-3, 3])
    ax1.set_ylim([-3, 3])


def test_cull_markers():
    x = np.random.random(20000)
    y = np.random.random(20000)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y, 'k.')
    ax.set_xlim(2, 3)

    pdf = io.BytesIO()
    fig.savefig(pdf, format="pdf")
    assert len(pdf.getvalue()) < 8000

    svg = io.BytesIO()
    fig.savefig(svg, format="svg")
    assert len(svg.getvalue()) < 20000


@image_comparison(baseline_images=['hatching'], remove_text=True,
                  style='default')
def test_hatching():
    fig, ax = plt.subplots(1, 1)

    # Default hatch color.
    rect1 = mpatches.Rectangle((0, 0), 3, 4, hatch='/')
    ax.add_patch(rect1)

    rect2 = mcollections.RegularPolyCollection(4, sizes=[16000],
                                               offsets=[(1.5, 6.5)],
                                               transOffset=ax.transData,
                                               hatch='/')
    ax.add_collection(rect2)

    # Ensure edge color is not applied to hatching.
    rect3 = mpatches.Rectangle((4, 0), 3, 4, hatch='/', edgecolor='C1')
    ax.add_patch(rect3)

    rect4 = mcollections.RegularPolyCollection(4, sizes=[16000],
                                               offsets=[(5.5, 6.5)],
                                               transOffset=ax.transData,
                                               hatch='/', edgecolor='C1')
    ax.add_collection(rect4)

    ax.set_xlim(0, 7)
    ax.set_ylim(0, 9)


def test_remove():
    fig, ax = plt.subplots()
    im = ax.imshow(np.arange(36).reshape(6, 6))
    ln, = ax.plot(range(5))

    assert fig.stale
    assert ax.stale

    fig.canvas.draw()
    assert not fig.stale
    assert not ax.stale
    assert not ln.stale

    assert im in ax.mouseover_set
    assert ln not in ax.mouseover_set
    assert im.axes is ax

    im.remove()
    ln.remove()

    for art in [im, ln]:
        assert art.axes is None
        assert art.figure is None

    assert im not in ax.mouseover_set
    assert fig.stale
    assert ax.stale


@image_comparison(baseline_images=["default_edges"], remove_text=True,
                  extensions=['png'], style='default')
def test_default_edges():
    fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2, 2)

    ax1.plot(np.arange(10), np.arange(10), 'x',
             np.arange(10) + 1, np.arange(10), 'o')
    ax2.bar(np.arange(10), np.arange(10), align='edge')
    ax3.text(0, 0, "BOX", size=24, bbox=dict(boxstyle='sawtooth'))
    ax3.set_xlim((-1, 1))
    ax3.set_ylim((-1, 1))
    pp1 = mpatches.PathPatch(
        mpath.Path([(0, 0), (1, 0), (1, 1), (0, 0)],
                   [mpath.Path.MOVETO, mpath.Path.CURVE3,
                    mpath.Path.CURVE3, mpath.Path.CLOSEPOLY]),
        fc="none", transform=ax4.transData)
    ax4.add_patch(pp1)


def test_properties():
    ln = mlines.Line2D([], [])
    with warnings.catch_warnings(record=True) as w:
        # Cause all warnings to always be triggered.
        warnings.simplefilter("always")
        ln.properties()
        assert len(w) == 0


def test_setp():
    # Check empty list
    plt.setp([])
    plt.setp([[]])

    # Check arbitrary iterables
    fig, axes = plt.subplots()
    lines1 = axes.plot(range(3))
    lines2 = axes.plot(range(3))
    martist.setp(chain(lines1, lines2), 'lw', 5)
    plt.setp(axes.spines.values(), color='green')

    # Check `file` argument
    sio = io.StringIO()
    plt.setp(lines1, 'zorder', file=sio)
    assert sio.getvalue() == '  zorder: float \n'


def test_None_zorder():
    fig, ax = plt.subplots()
    ln, = ax.plot(range(5), zorder=None)
    assert ln.get_zorder() == mlines.Line2D.zorder
    ln.set_zorder(123456)
    assert ln.get_zorder() == 123456
    ln.set_zorder(None)
    assert ln.get_zorder() == mlines.Line2D.zorder


@pytest.mark.parametrize('accept_clause, expected', [
    ('', 'unknown'),
    ("ACCEPTS: [ '-' | '--' | '-.' ]", "[ '-' | '--' | '-.' ] "),
    ('ACCEPTS: Some description.', 'Some description. '),
    ('.. ACCEPTS: Some description.', 'Some description. '),
])
def test_artist_inspector_get_valid_values(accept_clause, expected):
    class TestArtist(martist.Artist):
        def set_f(self):
            pass

    func = TestArtist.set_f
    if hasattr(func, '__func__'):
        func = func.__func__  # python 2 must write via __func__.__doc__
    func.__doc__ = """
    Some text.

    %s
    """ % accept_clause
    valid_values = martist.ArtistInspector(TestArtist).get_valid_values('f')
    assert valid_values == expected
