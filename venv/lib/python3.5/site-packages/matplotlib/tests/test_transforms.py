from __future__ import absolute_import, division, print_function

from six.moves import zip

import unittest

import numpy as np
from numpy.testing import (assert_allclose, assert_almost_equal,
                           assert_array_equal, assert_array_almost_equal)
import pytest

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.transforms as mtransforms
from matplotlib.path import Path
from matplotlib.scale import LogScale
from matplotlib.testing.decorators import image_comparison


def test_non_affine_caching():
    class AssertingNonAffineTransform(mtransforms.Transform):
        """
        This transform raises an assertion error when called when it
        shouldn't be and self.raise_on_transform is True.

        """
        input_dims = output_dims = 2
        is_affine = False

        def __init__(self, *args, **kwargs):
            mtransforms.Transform.__init__(self, *args, **kwargs)
            self.raise_on_transform = False
            self.underlying_transform = mtransforms.Affine2D().scale(10, 10)

        def transform_path_non_affine(self, path):
            assert not self.raise_on_transform, \
                'Invalidated affine part of transform unnecessarily.'
            return self.underlying_transform.transform_path(path)
        transform_path = transform_path_non_affine

        def transform_non_affine(self, path):
            assert not self.raise_on_transform, \
                'Invalidated affine part of transform unnecessarily.'
            return self.underlying_transform.transform(path)
        transform = transform_non_affine

    my_trans = AssertingNonAffineTransform()
    ax = plt.axes()
    plt.plot(np.arange(10), transform=my_trans + ax.transData)
    plt.draw()
    # enable the transform to raise an exception if it's non-affine transform
    # method is triggered again.
    my_trans.raise_on_transform = True
    ax.transAxes.invalidate()
    plt.draw()


def test_external_transform_api():
    class ScaledBy(object):
        def __init__(self, scale_factor):
            self._scale_factor = scale_factor

        def _as_mpl_transform(self, axes):
            return (mtransforms.Affine2D().scale(self._scale_factor)
                    + axes.transData)

    ax = plt.axes()
    line, = plt.plot(np.arange(10), transform=ScaledBy(10))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    # assert that the top transform of the line is the scale transform.
    assert_allclose(line.get_transform()._a.get_matrix(),
                    mtransforms.Affine2D().scale(10).get_matrix())


@image_comparison(baseline_images=['pre_transform_data'],
                  tol=0.08)
def test_pre_transform_plotting():
    # a catch-all for as many as possible plot layouts which handle
    # pre-transforming the data NOTE: The axis range is important in this
    # plot. It should be x10 what the data suggests it should be
    ax = plt.axes()
    times10 = mtransforms.Affine2D().scale(10)

    ax.contourf(np.arange(48).reshape(6, 8), transform=times10 + ax.transData)

    ax.pcolormesh(np.linspace(0, 4, 7),
                  np.linspace(5.5, 8, 9),
                  np.arange(48).reshape(8, 6),
                  transform=times10 + ax.transData)

    ax.scatter(np.linspace(0, 10), np.linspace(10, 0),
               transform=times10 + ax.transData)

    x = np.linspace(8, 10, 20)
    y = np.linspace(1, 5, 20)
    u = 2*np.sin(x) + np.cos(y[:, np.newaxis])
    v = np.sin(x) - np.cos(y[:, np.newaxis])

    df = 25. / 30.   # Compatibility factor for old test image
    ax.streamplot(x, y, u, v, transform=times10 + ax.transData,
                  density=(df, df), linewidth=u**2 + v**2)

    # reduce the vector data down a bit for barb and quiver plotting
    x, y = x[::3], y[::3]
    u, v = u[::3, ::3], v[::3, ::3]

    ax.quiver(x, y + 5, u, v, transform=times10 + ax.transData)

    ax.barbs(x - 3, y + 5, u**2, v**2, transform=times10 + ax.transData)


def test_contour_pre_transform_limits():
    ax = plt.axes()
    xs, ys = np.meshgrid(np.linspace(15, 20, 15), np.linspace(12.4, 12.5, 20))
    ax.contourf(xs, ys, np.log(xs * ys),
                transform=mtransforms.Affine2D().scale(0.1) + ax.transData)

    expected = np.array([[1.5, 1.24],
                         [2., 1.25]])
    assert_almost_equal(expected, ax.dataLim.get_points())


def test_pcolor_pre_transform_limits():
    # Based on test_contour_pre_transform_limits()
    ax = plt.axes()
    xs, ys = np.meshgrid(np.linspace(15, 20, 15), np.linspace(12.4, 12.5, 20))
    ax.pcolor(xs, ys, np.log(xs * ys),
              transform=mtransforms.Affine2D().scale(0.1) + ax.transData)

    expected = np.array([[1.5, 1.24],
                         [2., 1.25]])
    assert_almost_equal(expected, ax.dataLim.get_points())


def test_pcolormesh_pre_transform_limits():
    # Based on test_contour_pre_transform_limits()
    ax = plt.axes()
    xs, ys = np.meshgrid(np.linspace(15, 20, 15), np.linspace(12.4, 12.5, 20))
    ax.pcolormesh(xs, ys, np.log(xs * ys),
                  transform=mtransforms.Affine2D().scale(0.1) + ax.transData)

    expected = np.array([[1.5, 1.24],
                         [2., 1.25]])
    assert_almost_equal(expected, ax.dataLim.get_points())


def test_Affine2D_from_values():
    points = np.array([[0, 0],
                       [10, 20],
                       [-1, 0],
                       ])

    t = mtransforms.Affine2D.from_values(1, 0, 0, 0, 0, 0)
    actual = t.transform(points)
    expected = np.array([[0, 0], [10, 0], [-1, 0]])
    assert_almost_equal(actual, expected)

    t = mtransforms.Affine2D.from_values(0, 2, 0, 0, 0, 0)
    actual = t.transform(points)
    expected = np.array([[0, 0], [0, 20], [0, -2]])
    assert_almost_equal(actual, expected)

    t = mtransforms.Affine2D.from_values(0, 0, 3, 0, 0, 0)
    actual = t.transform(points)
    expected = np.array([[0, 0], [60, 0], [0, 0]])
    assert_almost_equal(actual, expected)

    t = mtransforms.Affine2D.from_values(0, 0, 0, 4, 0, 0)
    actual = t.transform(points)
    expected = np.array([[0, 0], [0, 80], [0, 0]])
    assert_almost_equal(actual, expected)

    t = mtransforms.Affine2D.from_values(0, 0, 0, 0, 5, 0)
    actual = t.transform(points)
    expected = np.array([[5, 0], [5, 0], [5, 0]])
    assert_almost_equal(actual, expected)

    t = mtransforms.Affine2D.from_values(0, 0, 0, 0, 0, 6)
    actual = t.transform(points)
    expected = np.array([[0, 6], [0, 6], [0, 6]])
    assert_almost_equal(actual, expected)


def test_clipping_of_log():
    # issue 804
    M, L, C = Path.MOVETO, Path.LINETO, Path.CLOSEPOLY
    points = [(0.2, -99), (0.4, -99), (0.4, 20), (0.2, 20), (0.2, -99)]
    codes = [M, L, L, L, C]
    path = Path(points, codes)

    # something like this happens in plotting logarithmic histograms
    trans = mtransforms.BlendedGenericTransform(mtransforms.Affine2D(),
                                            LogScale.Log10Transform('clip'))
    tpath = trans.transform_path_non_affine(path)
    result = tpath.iter_segments(trans.get_affine(),
                                 clip=(0, 0, 100, 100),
                                 simplify=False)

    tpoints, tcodes = zip(*result)
    assert_allclose(tcodes, [M, L, L, L, C])


class NonAffineForTest(mtransforms.Transform):
    """
    A class which looks like a non affine transform, but does whatever
    the given transform does (even if it is affine). This is very useful
    for testing NonAffine behaviour with a simple Affine transform.

    """
    is_affine = False
    output_dims = 2
    input_dims = 2

    def __init__(self, real_trans, *args, **kwargs):
        self.real_trans = real_trans
        mtransforms.Transform.__init__(self, *args, **kwargs)

    def transform_non_affine(self, values):
        return self.real_trans.transform(values)

    def transform_path_non_affine(self, path):
        return self.real_trans.transform_path(path)


class BasicTransformTests(unittest.TestCase):
    def setUp(self):

        self.ta1 = mtransforms.Affine2D(shorthand_name='ta1').rotate(np.pi / 2)
        self.ta2 = mtransforms.Affine2D(shorthand_name='ta2').translate(10, 0)
        self.ta3 = mtransforms.Affine2D(shorthand_name='ta3').scale(1, 2)

        self.tn1 = NonAffineForTest(mtransforms.Affine2D().translate(1, 2),
                                    shorthand_name='tn1')
        self.tn2 = NonAffineForTest(mtransforms.Affine2D().translate(1, 2),
                                    shorthand_name='tn2')
        self.tn3 = NonAffineForTest(mtransforms.Affine2D().translate(1, 2),
                                    shorthand_name='tn3')

        # creates a transform stack which looks like ((A, (N, A)), A)
        self.stack1 = (self.ta1 + (self.tn1 + self.ta2)) + self.ta3
        # creates a transform stack which looks like (((A, N), A), A)
        self.stack2 = self.ta1 + self.tn1 + self.ta2 + self.ta3
        # creates a transform stack which is a subset of stack2
        self.stack2_subset = self.tn1 + self.ta2 + self.ta3

        # when in debug, the transform stacks can produce dot images:
#        self.stack1.write_graphviz(file('stack1.dot', 'w'))
#        self.stack2.write_graphviz(file('stack2.dot', 'w'))
#        self.stack2_subset.write_graphviz(file('stack2_subset.dot', 'w'))

    def test_transform_depth(self):
        assert self.stack1.depth == 4
        assert self.stack2.depth == 4
        assert self.stack2_subset.depth == 3

    def test_left_to_right_iteration(self):
        stack3 = (self.ta1 + (self.tn1 + (self.ta2 + self.tn2))) + self.ta3
#        stack3.write_graphviz(file('stack3.dot', 'w'))

        target_transforms = [stack3,
                             (self.tn1 + (self.ta2 + self.tn2)) + self.ta3,
                             (self.ta2 + self.tn2) + self.ta3,
                             self.tn2 + self.ta3,
                             self.ta3,
                             ]
        r = [rh for _, rh in stack3._iter_break_from_left_to_right()]
        assert len(r) == len(target_transforms)

        for target_stack, stack in zip(target_transforms, r):
            assert target_stack == stack

    def test_transform_shortcuts(self):
        assert self.stack1 - self.stack2_subset == self.ta1
        assert self.stack2 - self.stack2_subset == self.ta1

        assert self.stack2_subset - self.stack2 == self.ta1.inverted()
        assert (self.stack2_subset - self.stack2).depth == 1

        with pytest.raises(ValueError):
            self.stack1 - self.stack2

        aff1 = self.ta1 + (self.ta2 + self.ta3)
        aff2 = self.ta2 + self.ta3

        assert aff1 - aff2 == self.ta1
        assert aff1 - self.ta2 == aff1 + self.ta2.inverted()

        assert self.stack1 - self.ta3 == self.ta1 + (self.tn1 + self.ta2)
        assert self.stack2 - self.ta3 == self.ta1 + self.tn1 + self.ta2

        assert ((self.ta2 + self.ta3) - self.ta3 + self.ta3 ==
                self.ta2 + self.ta3)

    def test_contains_branch(self):
        r1 = (self.ta2 + self.ta1)
        r2 = (self.ta2 + self.ta1)
        assert r1 == r2
        assert r1 != self.ta1
        assert r1.contains_branch(r2)
        assert r1.contains_branch(self.ta1)
        assert not r1.contains_branch(self.ta2)
        assert not r1.contains_branch((self.ta2 + self.ta2))

        assert r1 == r2

        assert self.stack1.contains_branch(self.ta3)
        assert self.stack2.contains_branch(self.ta3)

        assert self.stack1.contains_branch(self.stack2_subset)
        assert self.stack2.contains_branch(self.stack2_subset)

        assert not self.stack2_subset.contains_branch(self.stack1)
        assert not self.stack2_subset.contains_branch(self.stack2)

        assert self.stack1.contains_branch((self.ta2 + self.ta3))
        assert self.stack2.contains_branch((self.ta2 + self.ta3))

        assert not self.stack1.contains_branch((self.tn1 + self.ta2))

    def test_affine_simplification(self):
        # tests that a transform stack only calls as much is absolutely
        # necessary "non-affine" allowing the best possible optimization with
        # complex transformation stacks.
        points = np.array([[0, 0], [10, 20], [np.nan, 1], [-1, 0]],
                          dtype=np.float64)
        na_pts = self.stack1.transform_non_affine(points)
        all_pts = self.stack1.transform(points)

        na_expected = np.array([[1., 2.], [-19., 12.],
                                [np.nan, np.nan], [1., 1.]], dtype=np.float64)
        all_expected = np.array([[11., 4.], [-9., 24.],
                                 [np.nan, np.nan], [11., 2.]],
                                dtype=np.float64)

        # check we have the expected results from doing the affine part only
        assert_array_almost_equal(na_pts, na_expected)
        # check we have the expected results from a full transformation
        assert_array_almost_equal(all_pts, all_expected)
        # check we have the expected results from doing the transformation in
        # two steps
        assert_array_almost_equal(self.stack1.transform_affine(na_pts),
                                  all_expected)
        # check that getting the affine transformation first, then fully
        # transforming using that yields the same result as before.
        assert_array_almost_equal(self.stack1.get_affine().transform(na_pts),
                                  all_expected)

        # check that the affine part of stack1 & stack2 are equivalent
        # (i.e. the optimization is working)
        expected_result = (self.ta2 + self.ta3).get_matrix()
        result = self.stack1.get_affine().get_matrix()
        assert_array_equal(expected_result, result)

        result = self.stack2.get_affine().get_matrix()
        assert_array_equal(expected_result, result)


class TestTransformPlotInterface(unittest.TestCase):
    def tearDown(self):
        plt.close()

    def test_line_extent_axes_coords(self):
        # a simple line in axes coordinates
        ax = plt.axes()
        ax.plot([0.1, 1.2, 0.8], [0.9, 0.5, 0.8], transform=ax.transAxes)
        assert_array_equal(ax.dataLim.get_points(),
                           np.array([[np.inf, np.inf],
                                     [-np.inf, -np.inf]]))

    def test_line_extent_data_coords(self):
        # a simple line in data coordinates
        ax = plt.axes()
        ax.plot([0.1, 1.2, 0.8], [0.9, 0.5, 0.8], transform=ax.transData)
        assert_array_equal(ax.dataLim.get_points(),
                           np.array([[0.1,  0.5], [1.2,  0.9]]))

    def test_line_extent_compound_coords1(self):
        # a simple line in data coordinates in the y component, and in axes
        # coordinates in the x
        ax = plt.axes()
        trans = mtransforms.blended_transform_factory(ax.transAxes,
                                                      ax.transData)
        ax.plot([0.1, 1.2, 0.8], [35, -5, 18], transform=trans)
        assert_array_equal(ax.dataLim.get_points(),
                           np.array([[np.inf, -5.],
                                     [-np.inf, 35.]]))
        plt.close()

    def test_line_extent_predata_transform_coords(self):
        # a simple line in (offset + data) coordinates
        ax = plt.axes()
        trans = mtransforms.Affine2D().scale(10) + ax.transData
        ax.plot([0.1, 1.2, 0.8], [35, -5, 18], transform=trans)
        assert_array_equal(ax.dataLim.get_points(),
                           np.array([[1., -50.], [12., 350.]]))
        plt.close()

    def test_line_extent_compound_coords2(self):
        # a simple line in (offset + data) coordinates in the y component, and
        # in axes coordinates in the x
        ax = plt.axes()
        trans = mtransforms.blended_transform_factory(ax.transAxes,
            mtransforms.Affine2D().scale(10) + ax.transData)
        ax.plot([0.1, 1.2, 0.8], [35, -5, 18], transform=trans)
        assert_array_equal(ax.dataLim.get_points(),
                           np.array([[np.inf, -50.], [-np.inf, 350.]]))
        plt.close()

    def test_line_extents_affine(self):
        ax = plt.axes()
        offset = mtransforms.Affine2D().translate(10, 10)
        plt.plot(np.arange(10), transform=offset + ax.transData)
        expected_data_lim = np.array([[0., 0.], [9.,  9.]]) + 10
        assert_array_almost_equal(ax.dataLim.get_points(), expected_data_lim)

    def test_line_extents_non_affine(self):
        ax = plt.axes()
        offset = mtransforms.Affine2D().translate(10, 10)
        na_offset = NonAffineForTest(mtransforms.Affine2D().translate(10, 10))
        plt.plot(np.arange(10), transform=offset + na_offset + ax.transData)
        expected_data_lim = np.array([[0., 0.], [9.,  9.]]) + 20
        assert_array_almost_equal(ax.dataLim.get_points(), expected_data_lim)

    def test_pathc_extents_non_affine(self):
        ax = plt.axes()
        offset = mtransforms.Affine2D().translate(10, 10)
        na_offset = NonAffineForTest(mtransforms.Affine2D().translate(10, 10))
        pth = Path(np.array([[0, 0], [0, 10], [10, 10], [10, 0]]))
        patch = mpatches.PathPatch(pth,
                                   transform=offset + na_offset + ax.transData)
        ax.add_patch(patch)
        expected_data_lim = np.array([[0., 0.], [10.,  10.]]) + 20
        assert_array_almost_equal(ax.dataLim.get_points(), expected_data_lim)

    def test_pathc_extents_affine(self):
        ax = plt.axes()
        offset = mtransforms.Affine2D().translate(10, 10)
        pth = Path(np.array([[0, 0], [0, 10], [10, 10], [10, 0]]))
        patch = mpatches.PathPatch(pth, transform=offset + ax.transData)
        ax.add_patch(patch)
        expected_data_lim = np.array([[0., 0.], [10.,  10.]]) + 10
        assert_array_almost_equal(ax.dataLim.get_points(), expected_data_lim)

    def test_line_extents_for_non_affine_transData(self):
        ax = plt.axes(projection='polar')
        # add 10 to the radius of the data
        offset = mtransforms.Affine2D().translate(0, 10)

        plt.plot(np.arange(10), transform=offset + ax.transData)
        # the data lim of a polar plot is stored in coordinates
        # before a transData transformation, hence the data limits
        # are not what is being shown on the actual plot.
        expected_data_lim = np.array([[0., 0.], [9.,  9.]]) + [0, 10]
        assert_array_almost_equal(ax.dataLim.get_points(), expected_data_lim)


def assert_bbox_eq(bbox1, bbox2):
    assert_array_equal(bbox1.bounds, bbox2.bounds)


def test_bbox_intersection():
    bbox_from_ext = mtransforms.Bbox.from_extents
    inter = mtransforms.Bbox.intersection

    r1 = bbox_from_ext(0, 0, 1, 1)
    r2 = bbox_from_ext(0.5, 0.5, 1.5, 1.5)
    r3 = bbox_from_ext(0.5, 0, 0.75, 0.75)
    r4 = bbox_from_ext(0.5, 1.5, 1, 2.5)
    r5 = bbox_from_ext(1, 1, 2, 2)

    # self intersection -> no change
    assert_bbox_eq(inter(r1, r1), r1)
    # simple intersection
    assert_bbox_eq(inter(r1, r2), bbox_from_ext(0.5, 0.5, 1, 1))
    # r3 contains r2
    assert_bbox_eq(inter(r1, r3), r3)
    # no intersection
    assert inter(r1, r4) is None
    # single point
    assert_bbox_eq(inter(r1, r5), bbox_from_ext(1, 1, 1, 1))


def test_bbox_as_strings():
    b = mtransforms.Bbox([[.5, 0], [.75, .75]])
    assert_bbox_eq(b, eval(repr(b), {'Bbox': mtransforms.Bbox}))
    asdict = eval(str(b), {'Bbox': dict})
    for k, v in asdict.items():
        assert getattr(b, k) == v
    fmt = '.1f'
    asdict = eval(format(b, fmt), {'Bbox': dict})
    for k, v in asdict.items():
        assert eval(format(getattr(b, k), fmt)) == v


def test_transform_single_point():
    t = mtransforms.Affine2D()
    r = t.transform_affine((1, 1))
    assert r.shape == (2,)


def test_log_transform():
    # Tests that the last line runs without exception (previously the
    # transform would fail if one of the axes was logarithmic).
    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.transData.transform((1, 1))


def test_nan_overlap():
    a = mtransforms.Bbox([[0, 0], [1, 1]])
    b = mtransforms.Bbox([[0, 0], [1, np.nan]])
    assert not a.overlaps(b)


def test_transform_angles():
    t = mtransforms.Affine2D()  # Identity transform
    angles = np.array([20, 45, 60])
    points = np.array([[0, 0], [1, 1], [2, 2]])

    # Identity transform does not change angles
    new_angles = t.transform_angles(angles, points)
    assert_array_almost_equal(angles, new_angles)

    # points missing a 2nd dimension
    with pytest.raises(ValueError):
        t.transform_angles(angles, points[0:2, 0:1])

    # Number of angles != Number of points
    with pytest.raises(ValueError):
        t.transform_angles(angles, points[0:2, :])


def test_nonsingular():
    # test for zero-expansion type cases; other cases may be added later
    zero_expansion = np.array([-0.001, 0.001])
    cases = [(0, np.nan), (0, 0), (0, 7.9e-317)]
    for args in cases:
        out = np.array(mtransforms.nonsingular(*args))
        assert_array_equal(out, zero_expansion)


def test_invalid_arguments():
    t = mtransforms.Affine2D()
    # There are two different exceptions, since the wrong number of
    # dimensions is caught when constructing an array_view, and that
    # raises a ValueError, and a wrong shape with a possible number
    # of dimensions is caught by our CALL_CPP macro, which always
    # raises the less precise RuntimeError.
    with pytest.raises(ValueError):
        t.transform(1)
    with pytest.raises(ValueError):
        t.transform([[[1]]])
    with pytest.raises(RuntimeError):
        t.transform([])
    with pytest.raises(RuntimeError):
        t.transform([1])
    with pytest.raises(RuntimeError):
        t.transform([[1]])
    with pytest.raises(RuntimeError):
        t.transform([[1, 2, 3]])


def test_transformed_path():
    points = [(0, 0), (1, 0), (1, 1), (0, 1)]
    codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
    path = Path(points, codes)

    trans = mtransforms.Affine2D()
    trans_path = mtransforms.TransformedPath(path, trans)
    assert_allclose(trans_path.get_fully_transformed_path().vertices, points)

    # Changing the transform should change the result.
    r2 = 1 / np.sqrt(2)
    trans.rotate(np.pi / 4)
    assert_allclose(trans_path.get_fully_transformed_path().vertices,
                    [(0, 0), (r2, r2), (0, 2 * r2), (-r2, r2)],
                    atol=1e-15)

    # Changing the path does not change the result (it's cached).
    path.points = [(0, 0)] * 4
    assert_allclose(trans_path.get_fully_transformed_path().vertices,
                    [(0, 0), (r2, r2), (0, 2 * r2), (-r2, r2)],
                    atol=1e-15)


def test_transformed_patch_path():
    trans = mtransforms.Affine2D()
    patch = mpatches.Wedge((0, 0), 1, 45, 135, transform=trans)

    tpatch = mtransforms.TransformedPatchPath(patch)
    points = tpatch.get_fully_transformed_path().vertices

    # Changing the transform should change the result.
    trans.scale(2)
    assert_allclose(tpatch.get_fully_transformed_path().vertices, points * 2)

    # Changing the path should change the result (and cancel out the scaling
    # from the transform).
    patch.set_radius(0.5)
    assert_allclose(tpatch.get_fully_transformed_path().vertices, points)


@pytest.mark.parametrize('locked_element', ['x0', 'y0', 'x1', 'y1'])
def test_lockable_bbox(locked_element):
    other_elements = ['x0', 'y0', 'x1', 'y1']
    other_elements.remove(locked_element)

    orig = mtransforms.Bbox.unit()
    locked = mtransforms.LockableBbox(orig, **{locked_element: 2})

    # LockableBbox should keep its locked element as specified in __init__.
    assert getattr(locked, locked_element) == 2
    assert getattr(locked, 'locked_' + locked_element) == 2
    for elem in other_elements:
        assert getattr(locked, elem) == getattr(orig, elem)

    # Changing underlying Bbox should update everything but locked element.
    orig.set_points(orig.get_points() + 10)
    assert getattr(locked, locked_element) == 2
    assert getattr(locked, 'locked_' + locked_element) == 2
    for elem in other_elements:
        assert getattr(locked, elem) == getattr(orig, elem)

    # Unlocking element should revert values back to the underlying Bbox.
    setattr(locked, 'locked_' + locked_element, None)
    assert getattr(locked, 'locked_' + locked_element) is None
    assert np.all(orig.get_points() == locked.get_points())

    # Relocking an element should change its value, but not others.
    setattr(locked, 'locked_' + locked_element, 3)
    assert getattr(locked, locked_element) == 3
    assert getattr(locked, 'locked_' + locked_element) == 3
    for elem in other_elements:
        assert getattr(locked, elem) == getattr(orig, elem)
