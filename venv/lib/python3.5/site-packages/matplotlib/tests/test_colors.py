from __future__ import absolute_import, division, print_function

import copy
import six
import itertools
import warnings
from distutils.version import LooseVersion as V

import numpy as np
import pytest

from numpy.testing.utils import assert_array_equal, assert_array_almost_equal

from matplotlib import cycler
import matplotlib
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.colorbar as mcolorbar
import matplotlib.cbook as cbook
import matplotlib.pyplot as plt
from matplotlib.testing.decorators import image_comparison


def test_resample():
    """
    Github issue #6025 pointed to incorrect ListedColormap._resample;
    here we test the method for LinearSegmentedColormap as well.
    """
    n = 101
    colorlist = np.empty((n, 4), float)
    colorlist[:, 0] = np.linspace(0, 1, n)
    colorlist[:, 1] = 0.2
    colorlist[:, 2] = np.linspace(1, 0, n)
    colorlist[:, 3] = 0.7
    lsc = mcolors.LinearSegmentedColormap.from_list('lsc', colorlist)
    lc = mcolors.ListedColormap(colorlist)
    lsc3 = lsc._resample(3)
    lc3 = lc._resample(3)
    expected = np.array([[0.0, 0.2, 1.0, 0.7],
                         [0.5, 0.2, 0.5, 0.7],
                         [1.0, 0.2, 0.0, 0.7]], float)
    assert_array_almost_equal(lsc3([0, 0.5, 1]), expected)
    assert_array_almost_equal(lc3([0, 0.5, 1]), expected)


def test_colormap_copy():
    cm = plt.cm.Reds
    cm_copy = copy.copy(cm)
    with np.errstate(invalid='ignore'):
        ret1 = cm_copy([-1, 0, .5, 1, np.nan, np.inf])
    cm2 = copy.copy(cm_copy)
    cm2.set_bad('g')
    with np.errstate(invalid='ignore'):
        ret2 = cm_copy([-1, 0, .5, 1, np.nan, np.inf])
    assert_array_equal(ret1, ret2)


def test_colormap_endian():
    """
    Github issue #1005: a bug in putmask caused erroneous
    mapping of 1.0 when input from a non-native-byteorder
    array.
    """
    cmap = cm.get_cmap("jet")
    # Test under, over, and invalid along with values 0 and 1.
    a = [-0.5, 0, 0.5, 1, 1.5, np.nan]
    for dt in ["f2", "f4", "f8"]:
        anative = np.ma.masked_invalid(np.array(a, dtype=dt))
        aforeign = anative.byteswap().newbyteorder()
        assert_array_equal(cmap(anative), cmap(aforeign))


def test_BoundaryNorm():
    """
    Github issue #1258: interpolation was failing with numpy
    1.7 pre-release.
    """

    boundaries = [0, 1.1, 2.2]
    vals = [-1, 0, 1, 2, 2.2, 4]

    # Without interpolation
    expected = [-1, 0, 0, 1, 2, 2]
    ncolors = len(boundaries) - 1
    bn = mcolors.BoundaryNorm(boundaries, ncolors)
    assert_array_equal(bn(vals), expected)

    # ncolors != len(boundaries) - 1 triggers interpolation
    expected = [-1, 0, 0, 2, 3, 3]
    ncolors = len(boundaries)
    bn = mcolors.BoundaryNorm(boundaries, ncolors)
    assert_array_equal(bn(vals), expected)

    # more boundaries for a third color
    boundaries = [0, 1, 2, 3]
    vals = [-1, 0.1, 1.1, 2.2, 4]
    ncolors = 5
    expected = [-1, 0, 2, 4, 5]
    bn = mcolors.BoundaryNorm(boundaries, ncolors)
    assert_array_equal(bn(vals), expected)

    # a scalar as input should not trigger an error and should return a scalar
    boundaries = [0, 1, 2]
    vals = [-1, 0.1, 1.1, 2.2]
    bn = mcolors.BoundaryNorm(boundaries, 2)
    expected = [-1, 0, 1, 2]
    for v, ex in zip(vals, expected):
        ret = bn(v)
        assert isinstance(ret, six.integer_types)
        assert_array_equal(ret, ex)
        assert_array_equal(bn([v]), ex)

    # same with interp
    bn = mcolors.BoundaryNorm(boundaries, 3)
    expected = [-1, 0, 2, 3]
    for v, ex in zip(vals, expected):
        ret = bn(v)
        assert isinstance(ret, six.integer_types)
        assert_array_equal(ret, ex)
        assert_array_equal(bn([v]), ex)

    # Clipping
    bn = mcolors.BoundaryNorm(boundaries, 3, clip=True)
    expected = [0, 0, 2, 2]
    for v, ex in zip(vals, expected):
        ret = bn(v)
        assert isinstance(ret, six.integer_types)
        assert_array_equal(ret, ex)
        assert_array_equal(bn([v]), ex)

    # Masked arrays
    boundaries = [0, 1.1, 2.2]
    vals = np.ma.masked_invalid([-1., np.NaN, 0, 1.4, 9])

    # Without interpolation
    ncolors = len(boundaries) - 1
    bn = mcolors.BoundaryNorm(boundaries, ncolors)
    expected = np.ma.masked_array([-1, -99, 0, 1, 2], mask=[0, 1, 0, 0, 0])
    assert_array_equal(bn(vals), expected)

    # With interpolation
    bn = mcolors.BoundaryNorm(boundaries, len(boundaries))
    expected = np.ma.masked_array([-1, -99, 0, 2, 3], mask=[0, 1, 0, 0, 0])
    assert_array_equal(bn(vals), expected)

    # Non-trivial masked arrays
    vals = np.ma.masked_invalid([np.Inf, np.NaN])
    assert np.all(bn(vals).mask)
    vals = np.ma.masked_invalid([np.Inf])
    assert np.all(bn(vals).mask)


def test_LogNorm():
    """
    LogNorm ignored clip, now it has the same
    behavior as Normalize, e.g., values > vmax are bigger than 1
    without clip, with clip they are 1.
    """
    ln = mcolors.LogNorm(clip=True, vmax=5)
    assert_array_equal(ln([1, 6]), [0, 1.0])


def test_PowerNorm():
    a = np.array([0, 0.5, 1, 1.5], dtype=float)
    pnorm = mcolors.PowerNorm(1)
    norm = mcolors.Normalize()
    assert_array_almost_equal(norm(a), pnorm(a))

    a = np.array([-0.5, 0, 2, 4, 8], dtype=float)
    expected = [0, 0, 1/16, 1/4, 1]
    pnorm = mcolors.PowerNorm(2, vmin=0, vmax=8)
    assert_array_almost_equal(pnorm(a), expected)
    assert pnorm(a[0]) == expected[0]
    assert pnorm(a[2]) == expected[2]
    assert_array_almost_equal(a[1:], pnorm.inverse(pnorm(a))[1:])

    # Clip = True
    a = np.array([-0.5, 0, 1, 8, 16], dtype=float)
    expected = [0, 0, 0, 1, 1]
    pnorm = mcolors.PowerNorm(2, vmin=2, vmax=8, clip=True)
    assert_array_almost_equal(pnorm(a), expected)
    assert pnorm(a[0]) == expected[0]
    assert pnorm(a[-1]) == expected[-1]

    # Clip = True at call time
    a = np.array([-0.5, 0, 1, 8, 16], dtype=float)
    expected = [0, 0, 0, 1, 1]
    pnorm = mcolors.PowerNorm(2, vmin=2, vmax=8, clip=False)
    assert_array_almost_equal(pnorm(a, clip=True), expected)
    assert pnorm(a[0], clip=True) == expected[0]
    assert pnorm(a[-1], clip=True) == expected[-1]


def test_Normalize():
    norm = mcolors.Normalize()
    vals = np.arange(-10, 10, 1, dtype=float)
    _inverse_tester(norm, vals)
    _scalar_tester(norm, vals)
    _mask_tester(norm, vals)

    # Handle integer input correctly (don't overflow when computing max-min,
    # i.e. 127-(-128) here).
    vals = np.array([-128, 127], dtype=np.int8)
    norm = mcolors.Normalize(vals.min(), vals.max())
    assert_array_equal(np.asarray(norm(vals)), [0, 1])

    # Don't lose precision on longdoubles (float128 on Linux):
    # for array inputs...
    vals = np.array([1.2345678901, 9.8765432109], dtype=np.longdouble)
    norm = mcolors.Normalize(vals.min(), vals.max())
    assert_array_equal(np.asarray(norm(vals)), [0, 1])
    # and for scalar ones.
    eps = np.finfo(np.longdouble).resolution
    norm = plt.Normalize(1, 1 + 100 * eps)
    # This returns exactly 0.5 when longdouble is extended precision (80-bit),
    # but only a value close to it when it is quadruple precision (128-bit).
    assert 0 < norm(1 + 50 * eps) < 1


def test_SymLogNorm():
    """
    Test SymLogNorm behavior
    """
    norm = mcolors.SymLogNorm(3, vmax=5, linscale=1.2)
    vals = np.array([-30, -1, 2, 6], dtype=float)
    normed_vals = norm(vals)
    expected = [0., 0.53980074, 0.826991, 1.02758204]
    assert_array_almost_equal(normed_vals, expected)
    _inverse_tester(norm, vals)
    _scalar_tester(norm, vals)
    _mask_tester(norm, vals)

    # Ensure that specifying vmin returns the same result as above
    norm = mcolors.SymLogNorm(3, vmin=-30, vmax=5, linscale=1.2)
    normed_vals = norm(vals)
    assert_array_almost_equal(normed_vals, expected)


def test_SymLogNorm_colorbar():
    """
    Test un-called SymLogNorm in a colorbar.
    """
    norm = mcolors.SymLogNorm(0.1, vmin=-1, vmax=1, linscale=1)
    fig = plt.figure()
    cbar = mcolorbar.ColorbarBase(fig.add_subplot(111), norm=norm)
    plt.close(fig)


def test_SymLogNorm_single_zero():
    """
    Test SymLogNorm to ensure it is not adding sub-ticks to zero label
    """
    fig = plt.figure()
    norm = mcolors.SymLogNorm(1e-5, vmin=-1, vmax=1)
    cbar = mcolorbar.ColorbarBase(fig.add_subplot(111), norm=norm)
    ticks = cbar.get_ticks()
    assert sum(ticks == 0) == 1
    plt.close(fig)


def _inverse_tester(norm_instance, vals):
    """
    Checks if the inverse of the given normalization is working.
    """
    assert_array_almost_equal(norm_instance.inverse(norm_instance(vals)), vals)


def _scalar_tester(norm_instance, vals):
    """
    Checks if scalars and arrays are handled the same way.
    Tests only for float.
    """
    scalar_result = [norm_instance(float(v)) for v in vals]
    assert_array_almost_equal(scalar_result, norm_instance(vals))


def _mask_tester(norm_instance, vals):
    """
    Checks mask handling
    """
    masked_array = np.ma.array(vals)
    masked_array[0] = np.ma.masked
    assert_array_equal(masked_array.mask, norm_instance(masked_array).mask)


@image_comparison(baseline_images=['levels_and_colors'],
                  extensions=['png'])
def test_cmap_and_norm_from_levels_and_colors():
    data = np.linspace(-2, 4, 49).reshape(7, 7)
    levels = [-1, 2, 2.5, 3]
    colors = ['red', 'green', 'blue', 'yellow', 'black']
    extend = 'both'
    cmap, norm = mcolors.from_levels_and_colors(levels, colors, extend=extend)

    ax = plt.axes()
    m = plt.pcolormesh(data, cmap=cmap, norm=norm)
    plt.colorbar(m)

    # Hide the axes labels (but not the colorbar ones, as they are useful)
    ax.tick_params(labelleft=False, labelbottom=False)


def test_cmap_and_norm_from_levels_and_colors2():
    levels = [-1, 2, 2.5, 3]
    colors = ['red', (0, 1, 0), 'blue', (0.5, 0.5, 0.5), (0.0, 0.0, 0.0, 1.0)]
    clr = mcolors.to_rgba_array(colors)
    bad = (0.1, 0.1, 0.1, 0.1)
    no_color = (0.0, 0.0, 0.0, 0.0)
    masked_value = 'masked_value'

    # Define the test values which are of interest.
    # Note: levels are lev[i] <= v < lev[i+1]
    tests = [('both', None, {-2: clr[0],
                             -1: clr[1],
                             2: clr[2],
                             2.25: clr[2],
                             3: clr[4],
                             3.5: clr[4],
                             masked_value: bad}),

             ('min', -1, {-2: clr[0],
                          -1: clr[1],
                          2: clr[2],
                          2.25: clr[2],
                          3: no_color,
                          3.5: no_color,
                          masked_value: bad}),

             ('max', -1, {-2: no_color,
                          -1: clr[0],
                          2: clr[1],
                          2.25: clr[1],
                          3: clr[3],
                          3.5: clr[3],
                          masked_value: bad}),

             ('neither', -2, {-2: no_color,
                              -1: clr[0],
                              2: clr[1],
                              2.25: clr[1],
                              3: no_color,
                              3.5: no_color,
                              masked_value: bad}),
             ]

    for extend, i1, cases in tests:
        cmap, norm = mcolors.from_levels_and_colors(levels, colors[0:i1],
                                                    extend=extend)
        cmap.set_bad(bad)
        for d_val, expected_color in cases.items():
            if d_val == masked_value:
                d_val = np.ma.array([1], mask=True)
            else:
                d_val = [d_val]
            assert_array_equal(expected_color, cmap(norm(d_val))[0],
                               'Wih extend={0!r} and data '
                               'value={1!r}'.format(extend, d_val))

    with pytest.raises(ValueError):
        mcolors.from_levels_and_colors(levels, colors)


def test_rgb_hsv_round_trip():
    for a_shape in [(500, 500, 3), (500, 3), (1, 3), (3,)]:
        np.random.seed(0)
        tt = np.random.random(a_shape)
        assert_array_almost_equal(tt,
            mcolors.hsv_to_rgb(mcolors.rgb_to_hsv(tt)))
        assert_array_almost_equal(tt,
            mcolors.rgb_to_hsv(mcolors.hsv_to_rgb(tt)))


def test_autoscale_masked():
    # Test for #2336. Previously fully masked data would trigger a ValueError.
    data = np.ma.masked_all((12, 20))
    plt.pcolor(data)
    plt.draw()


def test_colors_no_float():
    # Gray must be a string to distinguish 3-4 grays from RGB or RGBA.
    with pytest.raises(ValueError):
        mcolors.to_rgba(0.4)


@image_comparison(baseline_images=['light_source_shading_topo'],
                  extensions=['png'])
def test_light_source_topo_surface():
    """Shades a DEM using different v.e.'s and blend modes."""
    fname = cbook.get_sample_data('jacksboro_fault_dem.npz', asfileobj=False)
    dem = np.load(fname)
    elev = dem['elevation']
    # Get the true cellsize in meters for accurate vertical exaggeration
    #   Convert from decimal degrees to meters
    dx, dy = dem['dx'], dem['dy']
    dx = 111320.0 * dx * np.cos(dem['ymin'])
    dy = 111320.0 * dy
    dem.close()

    ls = mcolors.LightSource(315, 45)
    cmap = cm.gist_earth

    fig, axes = plt.subplots(nrows=3, ncols=3)
    for row, mode in zip(axes, ['hsv', 'overlay', 'soft']):
        for ax, ve in zip(row, [0.1, 1, 10]):
            rgb = ls.shade(elev, cmap, vert_exag=ve, dx=dx, dy=dy,
                           blend_mode=mode)
            ax.imshow(rgb)
            ax.set(xticks=[], yticks=[])


def test_light_source_shading_default():
    """Array comparison test for the default "hsv" blend mode. Ensure the
    default result doesn't change without warning."""
    y, x = np.mgrid[-1.2:1.2:8j, -1.2:1.2:8j]
    z = 10 * np.cos(x**2 + y**2)

    cmap = plt.cm.copper
    ls = mcolors.LightSource(315, 45)
    rgb = ls.shade(z, cmap)

    # Result stored transposed and rounded for more compact display...
    expect = np.array(
        [[[0.00, 0.45, 0.90, 0.90, 0.82, 0.62, 0.28, 0.00],
          [0.45, 0.94, 0.99, 1.00, 1.00, 0.96, 0.65, 0.17],
          [0.90, 0.99, 1.00, 1.00, 1.00, 1.00, 0.94, 0.35],
          [0.90, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.49],
          [0.82, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.41],
          [0.62, 0.96, 1.00, 1.00, 1.00, 1.00, 0.90, 0.07],
          [0.28, 0.65, 0.94, 1.00, 1.00, 0.90, 0.35, 0.01],
          [0.00, 0.17, 0.35, 0.49, 0.41, 0.07, 0.01, 0.00]],

         [[0.00, 0.28, 0.59, 0.72, 0.62, 0.40, 0.18, 0.00],
          [0.28, 0.78, 0.93, 0.92, 0.83, 0.66, 0.39, 0.11],
          [0.59, 0.93, 0.99, 1.00, 0.92, 0.75, 0.50, 0.21],
          [0.72, 0.92, 1.00, 0.99, 0.93, 0.76, 0.51, 0.18],
          [0.62, 0.83, 0.92, 0.93, 0.87, 0.68, 0.42, 0.08],
          [0.40, 0.66, 0.75, 0.76, 0.68, 0.52, 0.23, 0.02],
          [0.18, 0.39, 0.50, 0.51, 0.42, 0.23, 0.00, 0.00],
          [0.00, 0.11, 0.21, 0.18, 0.08, 0.02, 0.00, 0.00]],

         [[0.00, 0.18, 0.38, 0.46, 0.39, 0.26, 0.11, 0.00],
          [0.18, 0.50, 0.70, 0.75, 0.64, 0.44, 0.25, 0.07],
          [0.38, 0.70, 0.91, 0.98, 0.81, 0.51, 0.29, 0.13],
          [0.46, 0.75, 0.98, 0.96, 0.84, 0.48, 0.22, 0.12],
          [0.39, 0.64, 0.81, 0.84, 0.71, 0.31, 0.11, 0.05],
          [0.26, 0.44, 0.51, 0.48, 0.31, 0.10, 0.03, 0.01],
          [0.11, 0.25, 0.29, 0.22, 0.11, 0.03, 0.00, 0.00],
          [0.00, 0.07, 0.13, 0.12, 0.05, 0.01, 0.00, 0.00]],

         [[1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00],
          [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00],
          [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00],
          [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00],
          [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00],
          [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00],
          [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00],
          [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]]
        ]).T

    if (V(np.__version__) == V('1.9.0')):
        # Numpy 1.9.0 uses a 2. order algorithm on the edges by default
        # This was changed back again in 1.9.1
        expect = expect[1:-1, 1:-1, :]
        rgb = rgb[1:-1, 1:-1, :]

    assert_array_almost_equal(rgb, expect, decimal=2)


@pytest.mark.xfail(V('1.7.0') <= V(np.__version__) <= V('1.9.0'),
                   reason='NumPy version is not buggy')
# Numpy 1.9.1 fixed a bug in masked arrays which resulted in
# additional elements being masked when calculating the gradient thus
# the output is different with earlier numpy versions.
def test_light_source_masked_shading():
    """Array comparison test for a surface with a masked portion. Ensures that
    we don't wind up with "fringes" of odd colors around masked regions."""
    y, x = np.mgrid[-1.2:1.2:8j, -1.2:1.2:8j]
    z = 10 * np.cos(x**2 + y**2)

    z = np.ma.masked_greater(z, 9.9)

    cmap = plt.cm.copper
    ls = mcolors.LightSource(315, 45)
    rgb = ls.shade(z, cmap)

    # Result stored transposed and rounded for more compact display...
    expect = np.array(
        [[[0.00, 0.46, 0.91, 0.91, 0.84, 0.64, 0.29, 0.00],
          [0.46, 0.96, 1.00, 1.00, 1.00, 0.97, 0.67, 0.18],
          [0.91, 1.00, 1.00, 1.00, 1.00, 1.00, 0.96, 0.36],
          [0.91, 1.00, 1.00, 0.00, 0.00, 1.00, 1.00, 0.51],
          [0.84, 1.00, 1.00, 0.00, 0.00, 1.00, 1.00, 0.44],
          [0.64, 0.97, 1.00, 1.00, 1.00, 1.00, 0.94, 0.09],
          [0.29, 0.67, 0.96, 1.00, 1.00, 0.94, 0.38, 0.01],
          [0.00, 0.18, 0.36, 0.51, 0.44, 0.09, 0.01, 0.00]],

         [[0.00, 0.29, 0.61, 0.75, 0.64, 0.41, 0.18, 0.00],
          [0.29, 0.81, 0.95, 0.93, 0.85, 0.68, 0.40, 0.11],
          [0.61, 0.95, 1.00, 0.78, 0.78, 0.77, 0.52, 0.22],
          [0.75, 0.93, 0.78, 0.00, 0.00, 0.78, 0.54, 0.19],
          [0.64, 0.85, 0.78, 0.00, 0.00, 0.78, 0.45, 0.08],
          [0.41, 0.68, 0.77, 0.78, 0.78, 0.55, 0.25, 0.02],
          [0.18, 0.40, 0.52, 0.54, 0.45, 0.25, 0.00, 0.00],
          [0.00, 0.11, 0.22, 0.19, 0.08, 0.02, 0.00, 0.00]],

         [[0.00, 0.19, 0.39, 0.48, 0.41, 0.26, 0.12, 0.00],
          [0.19, 0.52, 0.73, 0.78, 0.66, 0.46, 0.26, 0.07],
          [0.39, 0.73, 0.95, 0.50, 0.50, 0.53, 0.30, 0.14],
          [0.48, 0.78, 0.50, 0.00, 0.00, 0.50, 0.23, 0.12],
          [0.41, 0.66, 0.50, 0.00, 0.00, 0.50, 0.11, 0.05],
          [0.26, 0.46, 0.53, 0.50, 0.50, 0.11, 0.03, 0.01],
          [0.12, 0.26, 0.30, 0.23, 0.11, 0.03, 0.00, 0.00],
          [0.00, 0.07, 0.14, 0.12, 0.05, 0.01, 0.00, 0.00]],

         [[1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00],
          [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00],
          [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00],
          [1.00, 1.00, 1.00, 0.00, 0.00, 1.00, 1.00, 1.00],
          [1.00, 1.00, 1.00, 0.00, 0.00, 1.00, 1.00, 1.00],
          [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00],
          [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00],
          [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]],
        ]).T

    assert_array_almost_equal(rgb, expect, decimal=2)


def test_light_source_hillshading():
    """Compare the current hillshading method against one that should be
    mathematically equivalent. Illuminates a cone from a range of angles."""

    def alternative_hillshade(azimuth, elev, z):
        illum = _sph2cart(*_azimuth2math(azimuth, elev))
        illum = np.array(illum)

        dy, dx = np.gradient(-z)
        dy = -dy
        dz = np.ones_like(dy)
        normals = np.dstack([dx, dy, dz])
        dividers = np.zeros_like(z)[..., None]
        for i, mat in enumerate(normals):
            for j, vec in enumerate(mat):
                dividers[i, j, 0] = np.linalg.norm(vec)
        normals /= dividers
        # once we drop support for numpy 1.7.x the above can be written as
        # normals /= np.linalg.norm(normals, axis=2)[..., None]
        # aviding the double loop.

        intensity = np.tensordot(normals, illum, axes=(2, 0))
        intensity -= intensity.min()
        intensity /= intensity.ptp()
        return intensity

    y, x = np.mgrid[5:0:-1, :5]
    z = -np.hypot(x - x.mean(), y - y.mean())

    for az, elev in itertools.product(range(0, 390, 30), range(0, 105, 15)):
        ls = mcolors.LightSource(az, elev)
        h1 = ls.hillshade(z)
        h2 = alternative_hillshade(az, elev, z)
        assert_array_almost_equal(h1, h2)


def test_light_source_planar_hillshading():
    """Ensure that the illumination intensity is correct for planar
    surfaces."""

    def plane(azimuth, elevation, x, y):
        """Create a plane whose normal vector is at the given azimuth and
        elevation."""
        theta, phi = _azimuth2math(azimuth, elevation)
        a, b, c = _sph2cart(theta, phi)
        z = -(a*x + b*y) / c
        return z

    def angled_plane(azimuth, elevation, angle, x, y):
        """Create a plane whose normal vector is at an angle from the given
        azimuth and elevation."""
        elevation = elevation + angle
        if elevation > 90:
            azimuth = (azimuth + 180) % 360
            elevation = (90 - elevation) % 90
        return plane(azimuth, elevation, x, y)

    y, x = np.mgrid[5:0:-1, :5]
    for az, elev in itertools.product(range(0, 390, 30), range(0, 105, 15)):
        ls = mcolors.LightSource(az, elev)

        # Make a plane at a range of angles to the illumination
        for angle in range(0, 105, 15):
            z = angled_plane(az, elev, angle, x, y)
            h = ls.hillshade(z)
            assert_array_almost_equal(h, np.cos(np.radians(angle)))


def test_color_names():
    assert mcolors.to_hex("blue") == "#0000ff"
    assert mcolors.to_hex("xkcd:blue") == "#0343df"
    assert mcolors.to_hex("tab:blue") == "#1f77b4"


def _sph2cart(theta, phi):
    x = np.cos(theta) * np.sin(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(phi)
    return x, y, z


def _azimuth2math(azimuth, elevation):
    """Converts from clockwise-from-north and up-from-horizontal to
    mathematical conventions."""
    theta = np.radians((90 - azimuth) % 360)
    phi = np.radians(90 - elevation)
    return theta, phi


def test_pandas_iterable(pd):
    # Using a list or series yields equivalent
    # color maps, i.e the series isn't seen as
    # a single color
    lst = ['red', 'blue', 'green']
    s = pd.Series(lst)
    cm1 = mcolors.ListedColormap(lst, N=5)
    cm2 = mcolors.ListedColormap(s, N=5)
    assert_array_equal(cm1.colors, cm2.colors)


@pytest.mark.parametrize('name', sorted(cm.cmap_d))
def test_colormap_reversing(name):
    """Check the generated _lut data of a colormap and corresponding
    reversed colormap if they are almost the same."""
    cmap = plt.get_cmap(name)
    cmap_r = cmap.reversed()
    if not cmap_r._isinit:
        cmap._init()
        cmap_r._init()
    assert_array_almost_equal(cmap._lut[:-3], cmap_r._lut[-4::-1])


def test_cn():
    matplotlib.rcParams['axes.prop_cycle'] = cycler('color',
                                                    ['blue', 'r'])
    assert mcolors.to_hex("C0") == '#0000ff'
    assert mcolors.to_hex("C1") == '#ff0000'

    matplotlib.rcParams['axes.prop_cycle'] = cycler('color',
                                                    ['xkcd:blue', 'r'])
    assert mcolors.to_hex("C0") == '#0343df'
    assert mcolors.to_hex("C1") == '#ff0000'

    matplotlib.rcParams['axes.prop_cycle'] = cycler('color', ['8e4585', 'r'])

    assert mcolors.to_hex("C0") == '#8e4585'
    # if '8e4585' gets parsed as a float before it gets detected as a hex
    # colour it will be interpreted as a very large number.
    # this mustn't happen.
    assert mcolors.to_rgb("C0")[0] != np.inf


def test_conversions():
    # to_rgba_array("none") returns a (0, 4) array.
    assert_array_equal(mcolors.to_rgba_array("none"), np.zeros((0, 4)))
    # a list of grayscale levels, not a single color.
    assert_array_equal(
        mcolors.to_rgba_array([".2", ".5", ".8"]),
        np.vstack([mcolors.to_rgba(c) for c in [".2", ".5", ".8"]]))
    # alpha is properly set.
    assert mcolors.to_rgba((1, 1, 1), .5) == (1, 1, 1, .5)
    assert mcolors.to_rgba(".1", .5) == (.1, .1, .1, .5)
    # builtin round differs between py2 and py3.
    assert mcolors.to_hex((.7, .7, .7)) == "#b2b2b2"
    # hex roundtrip.
    hex_color = "#1234abcd"
    assert mcolors.to_hex(mcolors.to_rgba(hex_color), keep_alpha=True) == \
        hex_color


def test_grey_gray():
    color_mapping = mcolors._colors_full_map
    for k in color_mapping.keys():
        if 'grey' in k:
            assert color_mapping[k] == color_mapping[k.replace('grey', 'gray')]
        if 'gray' in k:
            assert color_mapping[k] == color_mapping[k.replace('gray', 'grey')]


def test_tableau_order():
    dflt_cycle = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
                  '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
                  '#bcbd22', '#17becf']

    assert list(mcolors.TABLEAU_COLORS.values()) == dflt_cycle


def test_ndarray_subclass_norm(recwarn):
    # Emulate an ndarray subclass that handles units
    # which objects when adding or subtracting with other
    # arrays. See #6622 and #8696
    class MyArray(np.ndarray):
        def __isub__(self, other):
            raise RuntimeError

        def __add__(self, other):
            raise RuntimeError

    data = np.arange(-10, 10, 1, dtype=float)
    data.shape = (10, 2)
    mydata = data.view(MyArray)

    for norm in [mcolors.Normalize(), mcolors.LogNorm(),
                 mcolors.SymLogNorm(3, vmax=5, linscale=1),
                 mcolors.Normalize(vmin=mydata.min(), vmax=mydata.max()),
                 mcolors.SymLogNorm(3, vmin=mydata.min(), vmax=mydata.max()),
                 mcolors.PowerNorm(1)]:
        assert_array_equal(norm(mydata), norm(data))
        fig, ax = plt.subplots()
        ax.imshow(mydata, norm=norm)
        fig.canvas.draw()
        if isinstance(norm, mcolors.PowerNorm):
            assert len(recwarn) == 1
            warn = recwarn.pop(UserWarning)
            assert ('Power-law scaling on negative values is ill-defined'
                    in str(warn.message))
        else:
            assert len(recwarn) == 0
        recwarn.clear()


def test_same_color():
    assert mcolors.same_color('k', (0, 0, 0))
    assert not mcolors.same_color('w', (1, 1, 0))
