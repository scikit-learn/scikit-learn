"""Tests for seaborn utility functions."""
import tempfile
from urllib.request import urlopen
from http.client import HTTPException

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from cycler import cycler

import pytest
from numpy.testing import (
    assert_array_equal,
)
from pandas.testing import (
    assert_series_equal,
    assert_frame_equal,
)

from distutils.version import LooseVersion

from .. import utils, rcmod
from ..utils import (
    get_dataset_names,
    get_color_cycle,
    remove_na,
    load_dataset,
    _assign_default_kwargs,
    _draw_figure,
)


a_norm = np.random.randn(100)


def _network(t=None, url="https://github.com"):
    """
    Decorator that will skip a test if `url` is unreachable.

    Parameters
    ----------
    t : function, optional
    url : str, optional

    """
    if t is None:
        return lambda x: _network(x, url=url)

    def wrapper(*args, **kwargs):
        # attempt to connect
        try:
            f = urlopen(url)
        except (IOError, HTTPException):
            pytest.skip("No internet connection")
        else:
            f.close()
            return t(*args, **kwargs)
    return wrapper


def test_pmf_hist_basics():
    """Test the function to return barplot args for pmf hist."""
    with pytest.warns(FutureWarning):
        out = utils.pmf_hist(a_norm)
    assert len(out) == 3
    x, h, w = out
    assert len(x) == len(h)

    # Test simple case
    a = np.arange(10)
    with pytest.warns(FutureWarning):
        x, h, w = utils.pmf_hist(a, 10)
    assert np.all(h == h[0])

    # Test width
    with pytest.warns(FutureWarning):
        x, h, w = utils.pmf_hist(a_norm)
    assert x[1] - x[0] == w

    # Test normalization
    with pytest.warns(FutureWarning):
        x, h, w = utils.pmf_hist(a_norm)
    assert sum(h) == pytest.approx(1)
    assert h.max() <= 1

    # Test bins
    with pytest.warns(FutureWarning):
        x, h, w = utils.pmf_hist(a_norm, 20)
    assert len(x) == 20


def test_ci_to_errsize():
    """Test behavior of ci_to_errsize."""
    cis = [[.5, .5],
           [1.25, 1.5]]

    heights = [1, 1.5]

    actual_errsize = np.array([[.5, 1],
                               [.25, 0]])

    test_errsize = utils.ci_to_errsize(cis, heights)
    assert_array_equal(actual_errsize, test_errsize)


def test_desaturate():
    """Test color desaturation."""
    out1 = utils.desaturate("red", .5)
    assert out1 == (.75, .25, .25)

    out2 = utils.desaturate("#00FF00", .5)
    assert out2 == (.25, .75, .25)

    out3 = utils.desaturate((0, 0, 1), .5)
    assert out3 == (.25, .25, .75)

    out4 = utils.desaturate("red", .5)
    assert out4 == (.75, .25, .25)


def test_desaturation_prop():
    """Test that pct outside of [0, 1] raises exception."""
    with pytest.raises(ValueError):
        utils.desaturate("blue", 50)


def test_saturate():
    """Test performance of saturation function."""
    out = utils.saturate((.75, .25, .25))
    assert out == (1, 0, 0)


@pytest.mark.parametrize(
    "p,annot", [(.0001, "***"), (.001, "**"), (.01, "*"), (.09, "."), (1, "")]
)
def test_sig_stars(p, annot):
    """Test the sig stars function."""
    with pytest.warns(FutureWarning):
        stars = utils.sig_stars(p)
        assert stars == annot


def test_iqr():
    """Test the IQR function."""
    a = np.arange(5)
    with pytest.warns(FutureWarning):
        iqr = utils.iqr(a)
    assert iqr == 2


@pytest.mark.parametrize(
    "s,exp",
    [
        ("a", "a"),
        ("abc", "abc"),
        (b"a", "a"),
        (b"abc", "abc"),
        (bytearray("abc", "utf-8"), "abc"),
        (bytearray(), ""),
        (1, "1"),
        (0, "0"),
        ([], str([])),
    ],
)
def test_to_utf8(s, exp):
    """Test the to_utf8 function: object to string"""
    u = utils.to_utf8(s)
    assert type(u) == str
    assert u == exp


class TestSpineUtils(object):

    sides = ["left", "right", "bottom", "top"]
    outer_sides = ["top", "right"]
    inner_sides = ["left", "bottom"]

    offset = 10
    original_position = ("outward", 0)
    offset_position = ("outward", offset)

    def test_despine(self):
        f, ax = plt.subplots()
        for side in self.sides:
            assert ax.spines[side].get_visible()

        utils.despine()
        for side in self.outer_sides:
            assert ~ax.spines[side].get_visible()
        for side in self.inner_sides:
            assert ax.spines[side].get_visible()

        utils.despine(**dict(zip(self.sides, [True] * 4)))
        for side in self.sides:
            assert ~ax.spines[side].get_visible()

    def test_despine_specific_axes(self):
        f, (ax1, ax2) = plt.subplots(2, 1)

        utils.despine(ax=ax2)

        for side in self.sides:
            assert ax1.spines[side].get_visible()

        for side in self.outer_sides:
            assert ~ax2.spines[side].get_visible()
        for side in self.inner_sides:
            assert ax2.spines[side].get_visible()

    def test_despine_with_offset(self):
        f, ax = plt.subplots()

        for side in self.sides:
            pos = ax.spines[side].get_position()
            assert pos == self.original_position

        utils.despine(ax=ax, offset=self.offset)

        for side in self.sides:
            is_visible = ax.spines[side].get_visible()
            new_position = ax.spines[side].get_position()
            if is_visible:
                assert new_position == self.offset_position
            else:
                assert new_position == self.original_position

    def test_despine_side_specific_offset(self):

        f, ax = plt.subplots()
        utils.despine(ax=ax, offset=dict(left=self.offset))

        for side in self.sides:
            is_visible = ax.spines[side].get_visible()
            new_position = ax.spines[side].get_position()
            if is_visible and side == "left":
                assert new_position == self.offset_position
            else:
                assert new_position == self.original_position

    def test_despine_with_offset_specific_axes(self):
        f, (ax1, ax2) = plt.subplots(2, 1)

        utils.despine(offset=self.offset, ax=ax2)

        for side in self.sides:
            pos1 = ax1.spines[side].get_position()
            pos2 = ax2.spines[side].get_position()
            assert pos1 == self.original_position
            if ax2.spines[side].get_visible():
                assert pos2 == self.offset_position
            else:
                assert pos2 == self.original_position

    def test_despine_trim_spines(self):

        f, ax = plt.subplots()
        ax.plot([1, 2, 3], [1, 2, 3])
        ax.set_xlim(.75, 3.25)

        utils.despine(trim=True)
        for side in self.inner_sides:
            bounds = ax.spines[side].get_bounds()
            assert bounds == (1, 3)

    def test_despine_trim_inverted(self):

        f, ax = plt.subplots()
        ax.plot([1, 2, 3], [1, 2, 3])
        ax.set_ylim(.85, 3.15)
        ax.invert_yaxis()

        utils.despine(trim=True)
        for side in self.inner_sides:
            bounds = ax.spines[side].get_bounds()
            assert bounds == (1, 3)

    def test_despine_trim_noticks(self):

        f, ax = plt.subplots()
        ax.plot([1, 2, 3], [1, 2, 3])
        ax.set_yticks([])
        utils.despine(trim=True)
        assert ax.get_yticks().size == 0

    def test_despine_trim_categorical(self):

        f, ax = plt.subplots()
        ax.plot(["a", "b", "c"], [1, 2, 3])

        utils.despine(trim=True)

        bounds = ax.spines["left"].get_bounds()
        assert bounds == (1, 3)

        bounds = ax.spines["bottom"].get_bounds()
        assert bounds == (0, 2)

    def test_despine_moved_ticks(self):

        f, ax = plt.subplots()
        for t in ax.yaxis.majorTicks:
            t.tick1line.set_visible(True)
        utils.despine(ax=ax, left=True, right=False)
        for t in ax.yaxis.majorTicks:
            assert t.tick2line.get_visible()
        plt.close(f)

        f, ax = plt.subplots()
        for t in ax.yaxis.majorTicks:
            t.tick1line.set_visible(False)
        utils.despine(ax=ax, left=True, right=False)
        for t in ax.yaxis.majorTicks:
            assert not t.tick2line.get_visible()
        plt.close(f)

        f, ax = plt.subplots()
        for t in ax.xaxis.majorTicks:
            t.tick1line.set_visible(True)
        utils.despine(ax=ax, bottom=True, top=False)
        for t in ax.xaxis.majorTicks:
            assert t.tick2line.get_visible()
        plt.close(f)

        f, ax = plt.subplots()
        for t in ax.xaxis.majorTicks:
            t.tick1line.set_visible(False)
        utils.despine(ax=ax, bottom=True, top=False)
        for t in ax.xaxis.majorTicks:
            assert not t.tick2line.get_visible()
        plt.close(f)


def test_ticklabels_overlap():

    rcmod.set()
    f, ax = plt.subplots(figsize=(2, 2))
    f.tight_layout()  # This gets the Agg renderer working

    assert not utils.axis_ticklabels_overlap(ax.get_xticklabels())

    big_strings = "abcdefgh", "ijklmnop"
    ax.set_xlim(-.5, 1.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(big_strings)

    assert utils.axis_ticklabels_overlap(ax.get_xticklabels())

    x, y = utils.axes_ticklabels_overlap(ax)
    assert x
    assert not y


def test_locator_to_legend_entries():

    locator = mpl.ticker.MaxNLocator(nbins=3)
    limits = (0.09, 0.4)
    levels, str_levels = utils.locator_to_legend_entries(
        locator, limits, float
    )
    assert str_levels == ["0.15", "0.30"]

    limits = (0.8, 0.9)
    levels, str_levels = utils.locator_to_legend_entries(
        locator, limits, float
    )
    assert str_levels == ["0.80", "0.84", "0.88"]

    limits = (1, 6)
    levels, str_levels = utils.locator_to_legend_entries(locator, limits, int)
    assert str_levels == ["2", "4", "6"]

    locator = mpl.ticker.LogLocator(numticks=5)
    limits = (5, 1425)
    levels, str_levels = utils.locator_to_legend_entries(locator, limits, int)
    if LooseVersion(mpl.__version__) >= "3.1":
        assert str_levels == ['10', '100', '1000']

    limits = (0.00003, 0.02)
    levels, str_levels = utils.locator_to_legend_entries(
        locator, limits, float
    )
    if LooseVersion(mpl.__version__) >= "3.1":
        assert str_levels == ['1e-04', '1e-03', '1e-02']


def test_move_legend_matplotlib_objects():

    fig, ax = plt.subplots()

    colors = "C2", "C5"
    labels = "first label", "second label"
    title = "the legend"

    for color, label in zip(colors, labels):
        ax.plot([0, 1], color=color, label=label)
    ax.legend(loc="upper right", title=title)
    utils._draw_figure(fig)
    xfm = ax.transAxes.inverted().transform

    # --- Test axes legend

    old_pos = xfm(ax.legend_.legendPatch.get_extents())

    new_fontsize = 14
    utils.move_legend(ax, "lower left", title_fontsize=new_fontsize)
    utils._draw_figure(fig)
    new_pos = xfm(ax.legend_.legendPatch.get_extents())

    assert (new_pos < old_pos).all()
    assert ax.legend_.get_title().get_text() == title
    assert ax.legend_.get_title().get_size() == new_fontsize

    # --- Test title replacement

    new_title = "new title"
    utils.move_legend(ax, "lower left", title=new_title)
    utils._draw_figure(fig)
    assert ax.legend_.get_title().get_text() == new_title

    # --- Test figure legend

    fig.legend(loc="upper right", title=title)
    _draw_figure(fig)
    xfm = fig.transFigure.inverted().transform
    old_pos = xfm(fig.legends[0].legendPatch.get_extents())

    utils.move_legend(fig, "lower left", title=new_title)
    _draw_figure(fig)

    new_pos = xfm(fig.legends[0].legendPatch.get_extents())
    assert (new_pos < old_pos).all()
    assert fig.legends[0].get_title().get_text() == new_title


def test_move_legend_grid_object(long_df):

    from seaborn.axisgrid import FacetGrid

    hue_var = "a"
    g = FacetGrid(long_df, hue=hue_var)
    g.map(plt.plot, "x", "y")

    g.add_legend()
    _draw_figure(g.figure)

    xfm = g.figure.transFigure.inverted().transform
    old_pos = xfm(g.legend.legendPatch.get_extents())

    fontsize = 20
    utils.move_legend(g, "lower left", title_fontsize=fontsize)
    _draw_figure(g.figure)

    new_pos = xfm(g.legend.legendPatch.get_extents())
    assert (new_pos < old_pos).all()
    assert g.legend.get_title().get_text() == hue_var
    assert g.legend.get_title().get_size() == fontsize

    assert g.legend.legendHandles
    for i, h in enumerate(g.legend.legendHandles):
        assert mpl.colors.to_rgb(h.get_color()) == mpl.colors.to_rgb(f"C{i}")


def test_move_legend_input_checks():

    ax = plt.figure().subplots()
    with pytest.raises(TypeError):
        utils.move_legend(ax.xaxis, "best")

    with pytest.raises(ValueError):
        utils.move_legend(ax, "best")

    with pytest.raises(ValueError):
        utils.move_legend(ax.figure, "best")


def check_load_dataset(name):
    ds = load_dataset(name, cache=False)
    assert(isinstance(ds, pd.DataFrame))


def check_load_cached_dataset(name):
    # Test the cacheing using a temporary file.
    with tempfile.TemporaryDirectory() as tmpdir:
        # download and cache
        ds = load_dataset(name, cache=True, data_home=tmpdir)

        # use cached version
        ds2 = load_dataset(name, cache=True, data_home=tmpdir)
        assert_frame_equal(ds, ds2)


@_network(url="https://github.com/mwaskom/seaborn-data")
def test_get_dataset_names():
    names = get_dataset_names()
    assert names
    assert "tips" in names


@_network(url="https://github.com/mwaskom/seaborn-data")
def test_load_datasets():

    # Heavy test to verify that we can load all available datasets
    for name in get_dataset_names():
        # unfortunately @network somehow obscures this generator so it
        # does not get in effect, so we need to call explicitly
        # yield check_load_dataset, name
        check_load_dataset(name)


@_network(url="https://github.com/mwaskom/seaborn-data")
def test_load_dataset_string_error():

    name = "bad_name"
    err = f"'{name}' is not one of the example datasets."
    with pytest.raises(ValueError, match=err):
        load_dataset(name)


def test_load_dataset_passed_data_error():

    df = pd.DataFrame()
    err = "This function accepts only strings"
    with pytest.raises(TypeError, match=err):
        load_dataset(df)


@_network(url="https://github.com/mwaskom/seaborn-data")
def test_load_cached_datasets():

    # Heavy test to verify that we can load all available datasets
    for name in get_dataset_names():
        # unfortunately @network somehow obscures this generator so it
        # does not get in effect, so we need to call explicitly
        # yield check_load_dataset, name
        check_load_cached_dataset(name)


def test_relative_luminance():
    """Test relative luminance."""
    out1 = utils.relative_luminance("white")
    assert out1 == 1

    out2 = utils.relative_luminance("#000000")
    assert out2 == 0

    out3 = utils.relative_luminance((.25, .5, .75))
    assert out3 == pytest.approx(0.201624536)

    rgbs = mpl.cm.RdBu(np.linspace(0, 1, 10))
    lums1 = [utils.relative_luminance(rgb) for rgb in rgbs]
    lums2 = utils.relative_luminance(rgbs)

    for lum1, lum2 in zip(lums1, lums2):
        assert lum1 == pytest.approx(lum2)


@pytest.mark.parametrize(
    "cycler,result",
    [
        (cycler(color=["y"]), ["y"]),
        (cycler(color=["k"]), ["k"]),
        (cycler(color=["k", "y"]), ["k", "y"]),
        (cycler(color=["y", "k"]), ["y", "k"]),
        (cycler(color=["b", "r"]), ["b", "r"]),
        (cycler(color=["r", "b"]), ["r", "b"]),
        (cycler(lw=[1, 2]), [".15"]),  # no color in cycle
    ],
)
def test_get_color_cycle(cycler, result):
    with mpl.rc_context(rc={"axes.prop_cycle": cycler}):
        assert get_color_cycle() == result


def test_remove_na():

    a_array = np.array([1, 2, np.nan, 3])
    a_array_rm = remove_na(a_array)
    assert_array_equal(a_array_rm, np.array([1, 2, 3]))

    a_series = pd.Series([1, 2, np.nan, 3])
    a_series_rm = remove_na(a_series)
    assert_series_equal(a_series_rm, pd.Series([1., 2, 3], [0, 1, 3]))


def test_assign_default_kwargs():

    def f(a, b, c, d):
        pass

    def g(c=1, d=2):
        pass

    kws = {"c": 3}

    kws = _assign_default_kwargs(kws, f, g)
    assert kws == {"c": 3, "d": 2}


def test_draw_figure():

    f, ax = plt.subplots()
    ax.plot(["a", "b", "c"], [1, 2, 3])
    _draw_figure(f)
    assert not f.stale
    # ticklabels are not populated until a draw, but this may change
    assert ax.get_xticklabels()[0].get_text() == "a"
