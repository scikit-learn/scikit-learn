from __future__ import absolute_import, division, print_function
import copy

import numpy as np

from numpy.testing import assert_array_equal
import pytest

from matplotlib.path import Path
from matplotlib.patches import Polygon
from matplotlib.testing.decorators import image_comparison
import matplotlib.pyplot as plt
from matplotlib import transforms


def test_readonly_path():
    path = Path.unit_circle()

    def modify_vertices():
        path.vertices = path.vertices * 2.0

    with pytest.raises(AttributeError):
        modify_vertices()


def test_point_in_path():
    # Test #1787
    verts2 = [(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)]

    path = Path(verts2, closed=True)
    points = [(0.5, 0.5), (1.5, 0.5)]
    ret = path.contains_points(points)
    assert ret.dtype == 'bool'
    assert np.all(ret == [True, False])


def test_contains_points_negative_radius():
    path = Path.unit_circle()

    points = [(0.0, 0.0), (1.25, 0.0), (0.9, 0.9)]
    expected = [True, False, False]
    result = path.contains_points(points, radius=-0.5)

    assert np.all(result == expected)


def test_point_in_path_nan():
    box = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]])
    p = Path(box)
    test = np.array([[np.nan, 0.5]])
    contains = p.contains_points(test)
    assert len(contains) == 1
    assert not contains[0]


def test_nonlinear_containment():
    fig, ax = plt.subplots()
    ax.set(xscale="log", ylim=(0, 1))
    polygon = ax.axvspan(1, 10)
    assert polygon.get_path().contains_point(
        ax.transData.transform_point((5, .5)), ax.transData)
    assert not polygon.get_path().contains_point(
        ax.transData.transform_point((.5, .5)), ax.transData)
    assert not polygon.get_path().contains_point(
        ax.transData.transform_point((50, .5)), ax.transData)


@image_comparison(baseline_images=['path_clipping'],
                  extensions=['svg'], remove_text=True)
def test_path_clipping():
    fig = plt.figure(figsize=(6.0, 6.2))

    for i, xy in enumerate([
            [(200, 200), (200, 350), (400, 350), (400, 200)],
            [(200, 200), (200, 350), (400, 350), (400, 100)],
            [(200, 100), (200, 350), (400, 350), (400, 100)],
            [(200, 100), (200, 415), (400, 350), (400, 100)],
            [(200, 100), (200, 415), (400, 415), (400, 100)],
            [(200, 415), (400, 415), (400, 100), (200, 100)],
            [(400, 415), (400, 100), (200, 100), (200, 415)]]):
        ax = fig.add_subplot(4, 2, i+1)
        bbox = [0, 140, 640, 260]
        ax.set_xlim(bbox[0], bbox[0] + bbox[2])
        ax.set_ylim(bbox[1], bbox[1] + bbox[3])
        ax.add_patch(Polygon(
            xy, facecolor='none', edgecolor='red', closed=True))


@image_comparison(baseline_images=['semi_log_with_zero'], extensions=['png'],
                  style='mpl20')
def test_log_transform_with_zero():
    x = np.arange(-10, 10)
    y = (1.0 - 1.0/(x**2+1))**20

    fig, ax = plt.subplots()
    ax.semilogy(x, y, "-o", lw=15, markeredgecolor='k')
    ax.set_ylim(1e-7, 1)
    ax.grid(True)


def test_make_compound_path_empty():
    # We should be able to make a compound path with no arguments.
    # This makes it easier to write generic path based code.
    r = Path.make_compound_path()
    assert r.vertices.shape == (0, 2)


@image_comparison(baseline_images=['xkcd'], remove_text=True)
def test_xkcd():
    np.random.seed(0)

    x = np.linspace(0, 2 * np.pi, 100)
    y = np.sin(x)

    with plt.xkcd():
        fig, ax = plt.subplots()
        ax.plot(x, y)


@image_comparison(baseline_images=['marker_paths'], extensions=['pdf'],
                  remove_text=True)
def test_marker_paths_pdf():
    N = 7

    plt.errorbar(np.arange(N),
                 np.ones(N) + 4,
                 np.ones(N))
    plt.xlim(-1, N)
    plt.ylim(-1, 7)


@image_comparison(baseline_images=['nan_path'], style='default',
                  remove_text=True, extensions=['pdf', 'svg', 'eps', 'png'])
def test_nan_isolated_points():

    y0 = [0, np.nan, 2, np.nan, 4, 5, 6]
    y1 = [np.nan, 7, np.nan, 9, 10, np.nan, 12]

    fig, ax = plt.subplots()

    ax.plot(y0, '-o')
    ax.plot(y1, '-o')


def test_path_no_doubled_point_in_to_polygon():
    hand = np.array(
        [[1.64516129, 1.16145833],
         [1.64516129, 1.59375],
         [1.35080645, 1.921875],
         [1.375, 2.18229167],
         [1.68548387, 1.9375],
         [1.60887097, 2.55208333],
         [1.68548387, 2.69791667],
         [1.76209677, 2.56770833],
         [1.83064516, 1.97395833],
         [1.89516129, 2.75],
         [1.9516129, 2.84895833],
         [2.01209677, 2.76041667],
         [1.99193548, 1.99479167],
         [2.11290323, 2.63020833],
         [2.2016129, 2.734375],
         [2.25403226, 2.60416667],
         [2.14919355, 1.953125],
         [2.30645161, 2.36979167],
         [2.39112903, 2.36979167],
         [2.41532258, 2.1875],
         [2.1733871, 1.703125],
         [2.07782258, 1.16666667]])

    (r0, c0, r1, c1) = (1.0, 1.5, 2.1, 2.5)

    poly = Path(np.vstack((hand[:, 1], hand[:, 0])).T, closed=True)
    clip_rect = transforms.Bbox([[r0, c0], [r1, c1]])
    poly_clipped = poly.clip_to_bbox(clip_rect).to_polygons()[0]

    assert np.all(poly_clipped[-2] != poly_clipped[-1])
    assert np.all(poly_clipped[-1] == poly_clipped[0])


def test_path_to_polygons():
    data = [[10, 10], [20, 20]]
    p = Path(data)

    assert_array_equal(p.to_polygons(width=40, height=40), [])
    assert_array_equal(p.to_polygons(width=40, height=40, closed_only=False),
                       [data])
    assert_array_equal(p.to_polygons(), [])
    assert_array_equal(p.to_polygons(closed_only=False), [data])

    data = [[10, 10], [20, 20], [30, 30]]
    closed_data = [[10, 10], [20, 20], [30, 30], [10, 10]]
    p = Path(data)

    assert_array_equal(p.to_polygons(width=40, height=40), [closed_data])
    assert_array_equal(p.to_polygons(width=40, height=40, closed_only=False),
                       [data])
    assert_array_equal(p.to_polygons(), [closed_data])
    assert_array_equal(p.to_polygons(closed_only=False), [data])


def test_path_deepcopy():
    # Should not raise any error
    verts = [[0, 0], [1, 1]]
    codes = [Path.MOVETO, Path.LINETO]
    path1 = Path(verts)
    path2 = Path(verts, codes)
    copy.deepcopy(path1)
    copy.deepcopy(path2)


@pytest.mark.parametrize('offset', range(-720, 361, 45))
def test_full_arc(offset):
    low = offset
    high = 360 + offset

    path = Path.arc(low, high)
    mins = np.min(path.vertices, axis=0)
    maxs = np.max(path.vertices, axis=0)
    np.testing.assert_allclose(mins, -1)
    assert np.allclose(maxs, 1)
