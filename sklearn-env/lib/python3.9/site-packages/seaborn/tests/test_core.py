import itertools
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

import pytest
from numpy.testing import assert_array_equal
from pandas.testing import assert_frame_equal

from ..axisgrid import FacetGrid
from .._core import (
    SemanticMapping,
    HueMapping,
    SizeMapping,
    StyleMapping,
    VectorPlotter,
    variable_type,
    infer_orient,
    unique_dashes,
    unique_markers,
    categorical_order,
)

from ..palettes import color_palette


try:
    from pandas import NA as PD_NA
except ImportError:
    PD_NA = None


class TestSemanticMapping:

    def test_call_lookup(self):

        m = SemanticMapping(VectorPlotter())
        lookup_table = dict(zip("abc", (1, 2, 3)))
        m.lookup_table = lookup_table
        for key, val in lookup_table.items():
            assert m(key) == val


class TestHueMapping:

    def test_init_from_map(self, long_df):

        p_orig = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", hue="a")
        )
        palette = "Set2"
        p = HueMapping.map(p_orig, palette=palette)
        assert p is p_orig
        assert isinstance(p._hue_map, HueMapping)
        assert p._hue_map.palette == palette

    def test_plotter_default_init(self, long_df):

        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y"),
        )
        assert isinstance(p._hue_map, HueMapping)
        assert p._hue_map.map_type is None

        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", hue="a"),
        )
        assert isinstance(p._hue_map, HueMapping)
        assert p._hue_map.map_type == p.var_types["hue"]

    def test_plotter_reinit(self, long_df):

        p_orig = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", hue="a"),
        )
        palette = "muted"
        hue_order = ["b", "a", "c"]
        p = p_orig.map_hue(palette=palette, order=hue_order)
        assert p is p_orig
        assert p._hue_map.palette == palette
        assert p._hue_map.levels == hue_order

    def test_hue_map_null(self, flat_series, null_series):

        p = VectorPlotter(variables=dict(x=flat_series, hue=null_series))
        m = HueMapping(p)
        assert m.levels is None
        assert m.map_type is None
        assert m.palette is None
        assert m.cmap is None
        assert m.norm is None
        assert m.lookup_table is None

    def test_hue_map_categorical(self, wide_df, long_df):

        p = VectorPlotter(data=wide_df)
        m = HueMapping(p)
        assert m.levels == wide_df.columns.tolist()
        assert m.map_type == "categorical"
        assert m.cmap is None

        # Test named palette
        palette = "Blues"
        expected_colors = color_palette(palette, wide_df.shape[1])
        expected_lookup_table = dict(zip(wide_df.columns, expected_colors))
        m = HueMapping(p, palette=palette)
        assert m.palette == "Blues"
        assert m.lookup_table == expected_lookup_table

        # Test list palette
        palette = color_palette("Reds", wide_df.shape[1])
        expected_lookup_table = dict(zip(wide_df.columns, palette))
        m = HueMapping(p, palette=palette)
        assert m.palette == palette
        assert m.lookup_table == expected_lookup_table

        # Test dict palette
        colors = color_palette("Set1", 8)
        palette = dict(zip(wide_df.columns, colors))
        m = HueMapping(p, palette=palette)
        assert m.palette == palette
        assert m.lookup_table == palette

        # Test dict with missing keys
        palette = dict(zip(wide_df.columns[:-1], colors))
        with pytest.raises(ValueError):
            HueMapping(p, palette=palette)

        # Test dict with missing keys
        palette = dict(zip(wide_df.columns[:-1], colors))
        with pytest.raises(ValueError):
            HueMapping(p, palette=palette)

        # Test list with wrong number of colors
        palette = colors[:-1]
        with pytest.raises(ValueError):
            HueMapping(p, palette=palette)

        # Test hue order
        hue_order = ["a", "c", "d"]
        m = HueMapping(p, order=hue_order)
        assert m.levels == hue_order

        # Test long data
        p = VectorPlotter(data=long_df, variables=dict(x="x", y="y", hue="a"))
        m = HueMapping(p)
        assert m.levels == categorical_order(long_df["a"])
        assert m.map_type == "categorical"
        assert m.cmap is None

        # Test default palette
        m = HueMapping(p)
        hue_levels = categorical_order(long_df["a"])
        expected_colors = color_palette(n_colors=len(hue_levels))
        expected_lookup_table = dict(zip(hue_levels, expected_colors))
        assert m.lookup_table == expected_lookup_table

        # Test missing data
        m = HueMapping(p)
        assert m(np.nan) == (0, 0, 0, 0)

        # Test default palette with many levels
        x = y = np.arange(26)
        hue = pd.Series(list("abcdefghijklmnopqrstuvwxyz"))
        p = VectorPlotter(variables=dict(x=x, y=y, hue=hue))
        m = HueMapping(p)
        expected_colors = color_palette("husl", n_colors=len(hue))
        expected_lookup_table = dict(zip(hue, expected_colors))
        assert m.lookup_table == expected_lookup_table

        # Test binary data
        p = VectorPlotter(data=long_df, variables=dict(x="x", y="y", hue="c"))
        m = HueMapping(p)
        assert m.levels == [0, 1]
        assert m.map_type == "categorical"

        for val in [0, 1]:
            p = VectorPlotter(
                data=long_df[long_df["c"] == val],
                variables=dict(x="x", y="y", hue="c"),
            )
            m = HueMapping(p)
            assert m.levels == [val]
            assert m.map_type == "categorical"

        # Test Timestamp data
        p = VectorPlotter(data=long_df, variables=dict(x="x", y="y", hue="t"))
        m = HueMapping(p)
        assert m.levels == [pd.Timestamp(t) for t in long_df["t"].unique()]
        assert m.map_type == "datetime"

        # Test excplicit categories
        p = VectorPlotter(data=long_df, variables=dict(x="x", hue="a_cat"))
        m = HueMapping(p)
        assert m.levels == long_df["a_cat"].cat.categories.tolist()
        assert m.map_type == "categorical"

        # Test numeric data with category type
        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", hue="s_cat")
        )
        m = HueMapping(p)
        assert m.levels == categorical_order(long_df["s_cat"])
        assert m.map_type == "categorical"
        assert m.cmap is None

        # Test categorical palette specified for numeric data
        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", hue="s")
        )
        palette = "deep"
        levels = categorical_order(long_df["s"])
        expected_colors = color_palette(palette, n_colors=len(levels))
        expected_lookup_table = dict(zip(levels, expected_colors))
        m = HueMapping(p, palette=palette)
        assert m.lookup_table == expected_lookup_table
        assert m.map_type == "categorical"

    def test_hue_map_numeric(self, long_df):

        # Test default colormap
        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", hue="s")
        )
        hue_levels = list(np.sort(long_df["s"].unique()))
        m = HueMapping(p)
        assert m.levels == hue_levels
        assert m.map_type == "numeric"
        assert m.cmap.name == "seaborn_cubehelix"

        # Test named colormap
        palette = "Purples"
        m = HueMapping(p, palette=palette)
        assert m.cmap is mpl.cm.get_cmap(palette)

        # Test colormap object
        palette = mpl.cm.get_cmap("Greens")
        m = HueMapping(p, palette=palette)
        assert m.cmap is mpl.cm.get_cmap(palette)

        # Test cubehelix shorthand
        palette = "ch:2,0,light=.2"
        m = HueMapping(p, palette=palette)
        assert isinstance(m.cmap, mpl.colors.ListedColormap)

        # Test specified hue limits
        hue_norm = 1, 4
        m = HueMapping(p, norm=hue_norm)
        assert isinstance(m.norm, mpl.colors.Normalize)
        assert m.norm.vmin == hue_norm[0]
        assert m.norm.vmax == hue_norm[1]

        # Test Normalize object
        hue_norm = mpl.colors.PowerNorm(2, vmin=1, vmax=10)
        m = HueMapping(p, norm=hue_norm)
        assert m.norm is hue_norm

        # Test default colormap values
        hmin, hmax = p.plot_data["hue"].min(), p.plot_data["hue"].max()
        m = HueMapping(p)
        assert m.lookup_table[hmin] == pytest.approx(m.cmap(0.0))
        assert m.lookup_table[hmax] == pytest.approx(m.cmap(1.0))

        # Test specified colormap values
        hue_norm = hmin - 1, hmax - 1
        m = HueMapping(p, norm=hue_norm)
        norm_min = (hmin - hue_norm[0]) / (hue_norm[1] - hue_norm[0])
        assert m.lookup_table[hmin] == pytest.approx(m.cmap(norm_min))
        assert m.lookup_table[hmax] == pytest.approx(m.cmap(1.0))

        # Test list of colors
        hue_levels = list(np.sort(long_df["s"].unique()))
        palette = color_palette("Blues", len(hue_levels))
        m = HueMapping(p, palette=palette)
        assert m.lookup_table == dict(zip(hue_levels, palette))

        palette = color_palette("Blues", len(hue_levels) + 1)
        with pytest.raises(ValueError):
            HueMapping(p, palette=palette)

        # Test dictionary of colors
        palette = dict(zip(hue_levels, color_palette("Reds")))
        m = HueMapping(p, palette=palette)
        assert m.lookup_table == palette

        palette.pop(hue_levels[0])
        with pytest.raises(ValueError):
            HueMapping(p, palette=palette)

        # Test invalid palette
        with pytest.raises(ValueError):
            HueMapping(p, palette="not a valid palette")

        # Test bad norm argument
        with pytest.raises(ValueError):
            HueMapping(p, norm="not a norm")


class TestSizeMapping:

    def test_init_from_map(self, long_df):

        p_orig = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", size="a")
        )
        sizes = 1, 6
        p = SizeMapping.map(p_orig, sizes=sizes)
        assert p is p_orig
        assert isinstance(p._size_map, SizeMapping)
        assert min(p._size_map.lookup_table.values()) == sizes[0]
        assert max(p._size_map.lookup_table.values()) == sizes[1]

    def test_plotter_default_init(self, long_df):

        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y"),
        )
        assert isinstance(p._size_map, SizeMapping)
        assert p._size_map.map_type is None

        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", size="a"),
        )
        assert isinstance(p._size_map, SizeMapping)
        assert p._size_map.map_type == p.var_types["size"]

    def test_plotter_reinit(self, long_df):

        p_orig = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", size="a"),
        )
        sizes = [1, 4, 2]
        size_order = ["b", "a", "c"]
        p = p_orig.map_size(sizes=sizes, order=size_order)
        assert p is p_orig
        assert p._size_map.lookup_table == dict(zip(size_order, sizes))
        assert p._size_map.levels == size_order

    def test_size_map_null(self, flat_series, null_series):

        p = VectorPlotter(variables=dict(x=flat_series, size=null_series))
        m = HueMapping(p)
        assert m.levels is None
        assert m.map_type is None
        assert m.norm is None
        assert m.lookup_table is None

    def test_map_size_numeric(self, long_df):

        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", size="s"),
        )

        # Test default range of keys in the lookup table values
        m = SizeMapping(p)
        size_values = m.lookup_table.values()
        value_range = min(size_values), max(size_values)
        assert value_range == p._default_size_range

        # Test specified range of size values
        sizes = 1, 5
        m = SizeMapping(p, sizes=sizes)
        size_values = m.lookup_table.values()
        assert min(size_values), max(size_values) == sizes

        # Test size values with normalization range
        norm = 1, 10
        m = SizeMapping(p, sizes=sizes, norm=norm)
        normalize = mpl.colors.Normalize(*norm, clip=True)
        for key, val in m.lookup_table.items():
            assert val == sizes[0] + (sizes[1] - sizes[0]) * normalize(key)

        # Test size values with normalization object
        norm = mpl.colors.LogNorm(1, 10, clip=False)
        m = SizeMapping(p, sizes=sizes, norm=norm)
        assert m.norm.clip
        for key, val in m.lookup_table.items():
            assert val == sizes[0] + (sizes[1] - sizes[0]) * norm(key)

        # Test bad sizes argument
        with pytest.raises(ValueError):
            SizeMapping(p, sizes="bad_sizes")

        # Test bad sizes argument
        with pytest.raises(ValueError):
            SizeMapping(p, sizes=(1, 2, 3))

        # Test bad norm argument
        with pytest.raises(ValueError):
            SizeMapping(p, norm="bad_norm")

    def test_map_size_categorical(self, long_df):

        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", size="a"),
        )

        # Test specified size order
        levels = p.plot_data["size"].unique()
        sizes = [1, 4, 6]
        order = [levels[1], levels[2], levels[0]]
        m = SizeMapping(p, sizes=sizes, order=order)
        assert m.lookup_table == dict(zip(order, sizes))

        # Test list of sizes
        order = categorical_order(p.plot_data["size"])
        sizes = list(np.random.rand(len(levels)))
        m = SizeMapping(p, sizes=sizes)
        assert m.lookup_table == dict(zip(order, sizes))

        # Test dict of sizes
        sizes = dict(zip(levels, np.random.rand(len(levels))))
        m = SizeMapping(p, sizes=sizes)
        assert m.lookup_table == sizes

        # Test specified size range
        sizes = (2, 5)
        m = SizeMapping(p, sizes=sizes)
        values = np.linspace(*sizes, len(m.levels))[::-1]
        assert m.lookup_table == dict(zip(m.levels, values))

        # Test explicit categories
        p = VectorPlotter(data=long_df, variables=dict(x="x", size="a_cat"))
        m = SizeMapping(p)
        assert m.levels == long_df["a_cat"].cat.categories.tolist()
        assert m.map_type == "categorical"

        # Test sizes list with wrong length
        sizes = list(np.random.rand(len(levels) + 1))
        with pytest.raises(ValueError):
            SizeMapping(p, sizes=sizes)

        # Test sizes dict with missing levels
        sizes = dict(zip(levels, np.random.rand(len(levels) - 1)))
        with pytest.raises(ValueError):
            SizeMapping(p, sizes=sizes)

        # Test bad sizes argument
        with pytest.raises(ValueError):
            SizeMapping(p, sizes="bad_size")


class TestStyleMapping:

    def test_init_from_map(self, long_df):

        p_orig = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", style="a")
        )
        markers = ["s", "p", "h"]
        p = StyleMapping.map(p_orig, markers=markers)
        assert p is p_orig
        assert isinstance(p._style_map, StyleMapping)
        assert p._style_map(p._style_map.levels, "marker") == markers

    def test_plotter_default_init(self, long_df):

        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y"),
        )
        assert isinstance(p._style_map, StyleMapping)

        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", style="a"),
        )
        assert isinstance(p._style_map, StyleMapping)

    def test_plotter_reinit(self, long_df):

        p_orig = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", style="a"),
        )
        markers = ["s", "p", "h"]
        style_order = ["b", "a", "c"]
        p = p_orig.map_style(markers=markers, order=style_order)
        assert p is p_orig
        assert p._style_map.levels == style_order
        assert p._style_map(style_order, "marker") == markers

    def test_style_map_null(self, flat_series, null_series):

        p = VectorPlotter(variables=dict(x=flat_series, style=null_series))
        m = HueMapping(p)
        assert m.levels is None
        assert m.map_type is None
        assert m.lookup_table is None

    def test_map_style(self, long_df):

        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", style="a"),
        )

        # Test defaults
        m = StyleMapping(p, markers=True, dashes=True)

        n = len(m.levels)
        for key, dashes in zip(m.levels, unique_dashes(n)):
            assert m(key, "dashes") == dashes

        actual_marker_paths = {
            k: mpl.markers.MarkerStyle(m(k, "marker")).get_path()
            for k in m.levels
        }
        expected_marker_paths = {
            k: mpl.markers.MarkerStyle(m).get_path()
            for k, m in zip(m.levels, unique_markers(n))
        }
        assert actual_marker_paths == expected_marker_paths

        # Test lists
        markers, dashes = ["o", "s", "d"], [(1, 0), (1, 1), (2, 1, 3, 1)]
        m = StyleMapping(p, markers=markers, dashes=dashes)
        for key, mark, dash in zip(m.levels, markers, dashes):
            assert m(key, "marker") == mark
            assert m(key, "dashes") == dash

        # Test dicts
        markers = dict(zip(p.plot_data["style"].unique(), markers))
        dashes = dict(zip(p.plot_data["style"].unique(), dashes))
        m = StyleMapping(p, markers=markers, dashes=dashes)
        for key in m.levels:
            assert m(key, "marker") == markers[key]
            assert m(key, "dashes") == dashes[key]

        # Test excplicit categories
        p = VectorPlotter(data=long_df, variables=dict(x="x", style="a_cat"))
        m = StyleMapping(p)
        assert m.levels == long_df["a_cat"].cat.categories.tolist()

        # Test style order with defaults
        order = p.plot_data["style"].unique()[[1, 2, 0]]
        m = StyleMapping(p, markers=True, dashes=True, order=order)
        n = len(order)
        for key, mark, dash in zip(order, unique_markers(n), unique_dashes(n)):
            assert m(key, "dashes") == dash
            assert m(key, "marker") == mark
            obj = mpl.markers.MarkerStyle(mark)
            path = obj.get_path().transformed(obj.get_transform())
            assert_array_equal(m(key, "path").vertices, path.vertices)

        # Test too many levels with style lists
        with pytest.raises(ValueError):
            StyleMapping(p, markers=["o", "s"], dashes=False)

        with pytest.raises(ValueError):
            StyleMapping(p, markers=False, dashes=[(2, 1)])

        # Test too many levels with style dicts
        markers, dashes = {"a": "o", "b": "s"}, False
        with pytest.raises(ValueError):
            StyleMapping(p, markers=markers, dashes=dashes)

        markers, dashes = False, {"a": (1, 0), "b": (2, 1)}
        with pytest.raises(ValueError):
            StyleMapping(p, markers=markers, dashes=dashes)

        # Test mixture of filled and unfilled markers
        markers, dashes = ["o", "x", "s"], None
        with pytest.raises(ValueError):
            StyleMapping(p, markers=markers, dashes=dashes)


class TestVectorPlotter:

    def test_flat_variables(self, flat_data):

        p = VectorPlotter()
        p.assign_variables(data=flat_data)
        assert p.input_format == "wide"
        assert list(p.variables) == ["x", "y"]
        assert len(p.plot_data) == len(flat_data)

        try:
            expected_x = flat_data.index
            expected_x_name = flat_data.index.name
        except AttributeError:
            expected_x = np.arange(len(flat_data))
            expected_x_name = None

        x = p.plot_data["x"]
        assert_array_equal(x, expected_x)

        expected_y = flat_data
        expected_y_name = getattr(flat_data, "name", None)

        y = p.plot_data["y"]
        assert_array_equal(y, expected_y)

        assert p.variables["x"] == expected_x_name
        assert p.variables["y"] == expected_y_name

    # TODO note that most of the other tests that exercise the core
    # variable assignment code still live in test_relational

    @pytest.mark.parametrize("name", [3, 4.5])
    def test_long_numeric_name(self, long_df, name):

        long_df[name] = long_df["x"]
        p = VectorPlotter()
        p.assign_variables(data=long_df, variables={"x": name})
        assert_array_equal(p.plot_data["x"], long_df[name])
        assert p.variables["x"] == name

    def test_long_hierarchical_index(self, rng):

        cols = pd.MultiIndex.from_product([["a"], ["x", "y"]])
        data = rng.uniform(size=(50, 2))
        df = pd.DataFrame(data, columns=cols)

        name = ("a", "y")
        var = "y"

        p = VectorPlotter()
        p.assign_variables(data=df, variables={var: name})
        assert_array_equal(p.plot_data[var], df[name])
        assert p.variables[var] == name

    def test_long_scalar_and_data(self, long_df):

        val = 22
        p = VectorPlotter(data=long_df, variables={"x": "x", "y": val})
        assert (p.plot_data["y"] == val).all()
        assert p.variables["y"] is None

    def test_wide_semantic_error(self, wide_df):

        err = "The following variable cannot be assigned with wide-form data: `hue`"
        with pytest.raises(ValueError, match=err):
            VectorPlotter(data=wide_df, variables={"hue": "a"})

    def test_long_unknown_error(self, long_df):

        err = "Could not interpret value `what` for parameter `hue`"
        with pytest.raises(ValueError, match=err):
            VectorPlotter(data=long_df, variables={"x": "x", "hue": "what"})

    def test_long_unmatched_size_error(self, long_df, flat_array):

        err = "Length of ndarray vectors must match length of `data`"
        with pytest.raises(ValueError, match=err):
            VectorPlotter(data=long_df, variables={"x": "x", "hue": flat_array})

    def test_wide_categorical_columns(self, wide_df):

        wide_df.columns = pd.CategoricalIndex(wide_df.columns)
        p = VectorPlotter(data=wide_df)
        assert_array_equal(p.plot_data["hue"].unique(), ["a", "b", "c"])

    def test_iter_data_quantitites(self, long_df):

        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y"),
        )
        out = p.iter_data("hue")
        assert len(list(out)) == 1

        var = "a"
        n_subsets = len(long_df[var].unique())

        semantics = ["hue", "size", "style"]
        for semantic in semantics:

            p = VectorPlotter(
                data=long_df,
                variables={"x": "x", "y": "y", semantic: var},
            )
            out = p.iter_data(semantics)
            assert len(list(out)) == n_subsets

        var = "a"
        n_subsets = len(long_df[var].unique())

        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", hue=var, style=var),
        )
        out = p.iter_data(semantics)
        assert len(list(out)) == n_subsets

        # --

        out = p.iter_data(semantics, reverse=True)
        assert len(list(out)) == n_subsets

        # --

        var1, var2 = "a", "s"

        n_subsets = len(long_df[var1].unique())

        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", hue=var1, style=var2),
        )
        out = p.iter_data(["hue"])
        assert len(list(out)) == n_subsets

        n_subsets = len(set(list(map(tuple, long_df[[var1, var2]].values))))

        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", hue=var1, style=var2),
        )
        out = p.iter_data(semantics)
        assert len(list(out)) == n_subsets

        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", hue=var1, size=var2, style=var1),
        )
        out = p.iter_data(semantics)
        assert len(list(out)) == n_subsets

        # --

        var1, var2, var3 = "a", "s", "b"
        cols = [var1, var2, var3]
        n_subsets = len(set(list(map(tuple, long_df[cols].values))))

        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", hue=var1, size=var2, style=var3),
        )
        out = p.iter_data(semantics)
        assert len(list(out)) == n_subsets

    def test_iter_data_keys(self, long_df):

        semantics = ["hue", "size", "style"]

        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y"),
        )
        for sub_vars, _ in p.iter_data("hue"):
            assert sub_vars == {}

        # --

        var = "a"

        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", hue=var),
        )
        for sub_vars, _ in p.iter_data("hue"):
            assert list(sub_vars) == ["hue"]
            assert sub_vars["hue"] in long_df[var].values

        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", size=var),
        )
        for sub_vars, _ in p.iter_data("size"):
            assert list(sub_vars) == ["size"]
            assert sub_vars["size"] in long_df[var].values

        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", hue=var, style=var),
        )
        for sub_vars, _ in p.iter_data(semantics):
            assert list(sub_vars) == ["hue", "style"]
            assert sub_vars["hue"] in long_df[var].values
            assert sub_vars["style"] in long_df[var].values
            assert sub_vars["hue"] == sub_vars["style"]

        var1, var2 = "a", "s"

        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", hue=var1, size=var2),
        )
        for sub_vars, _ in p.iter_data(semantics):
            assert list(sub_vars) == ["hue", "size"]
            assert sub_vars["hue"] in long_df[var1].values
            assert sub_vars["size"] in long_df[var2].values

        semantics = ["hue", "col", "row"]
        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", hue=var1, col=var2),
        )
        for sub_vars, _ in p.iter_data("hue"):
            assert list(sub_vars) == ["hue", "col"]
            assert sub_vars["hue"] in long_df[var1].values
            assert sub_vars["col"] in long_df[var2].values

    def test_iter_data_values(self, long_df):

        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y"),
        )

        p.sort = True
        _, sub_data = next(p.iter_data("hue"))
        assert_frame_equal(sub_data, p.plot_data)

        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", hue="a"),
        )

        for sub_vars, sub_data in p.iter_data("hue"):
            rows = p.plot_data["hue"] == sub_vars["hue"]
            assert_frame_equal(sub_data, p.plot_data[rows])

        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", hue="a", size="s"),
        )
        for sub_vars, sub_data in p.iter_data(["hue", "size"]):
            rows = p.plot_data["hue"] == sub_vars["hue"]
            rows &= p.plot_data["size"] == sub_vars["size"]
            assert_frame_equal(sub_data, p.plot_data[rows])

    def test_iter_data_reverse(self, long_df):

        reversed_order = categorical_order(long_df["a"])[::-1]
        p = VectorPlotter(
            data=long_df,
            variables=dict(x="x", y="y", hue="a")
        )
        iterator = p.iter_data("hue", reverse=True)
        for i, (sub_vars, _) in enumerate(iterator):
            assert sub_vars["hue"] == reversed_order[i]

    def test_axis_labels(self, long_df):

        f, ax = plt.subplots()

        p = VectorPlotter(data=long_df, variables=dict(x="a"))

        p._add_axis_labels(ax)
        assert ax.get_xlabel() == "a"
        assert ax.get_ylabel() == ""
        ax.clear()

        p = VectorPlotter(data=long_df, variables=dict(y="a"))
        p._add_axis_labels(ax)
        assert ax.get_xlabel() == ""
        assert ax.get_ylabel() == "a"
        ax.clear()

        p = VectorPlotter(data=long_df, variables=dict(x="a"))

        p._add_axis_labels(ax, default_y="default")
        assert ax.get_xlabel() == "a"
        assert ax.get_ylabel() == "default"
        ax.clear()

        p = VectorPlotter(data=long_df, variables=dict(y="a"))
        p._add_axis_labels(ax, default_x="default", default_y="default")
        assert ax.get_xlabel() == "default"
        assert ax.get_ylabel() == "a"
        ax.clear()

        p = VectorPlotter(data=long_df, variables=dict(x="x", y="a"))
        ax.set(xlabel="existing", ylabel="also existing")
        p._add_axis_labels(ax)
        assert ax.get_xlabel() == "existing"
        assert ax.get_ylabel() == "also existing"

        f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
        p = VectorPlotter(data=long_df, variables=dict(x="x", y="y"))

        p._add_axis_labels(ax1)
        p._add_axis_labels(ax2)

        assert ax1.get_xlabel() == "x"
        assert ax1.get_ylabel() == "y"
        assert ax1.yaxis.label.get_visible()

        assert ax2.get_xlabel() == "x"
        assert ax2.get_ylabel() == "y"
        assert not ax2.yaxis.label.get_visible()

    @pytest.mark.parametrize(
        "variables",
        [
            dict(x="x", y="y"),
            dict(x="x"),
            dict(y="y"),
            dict(x="t", y="y"),
            dict(x="x", y="a"),
        ]
    )
    def test_attach_basics(self, long_df, variables):

        _, ax = plt.subplots()
        p = VectorPlotter(data=long_df, variables=variables)
        p._attach(ax)
        assert p.ax is ax

    def test_attach_disallowed(self, long_df):

        _, ax = plt.subplots()
        p = VectorPlotter(data=long_df, variables={"x": "a"})

        with pytest.raises(TypeError):
            p._attach(ax, allowed_types="numeric")

        with pytest.raises(TypeError):
            p._attach(ax, allowed_types=["datetime", "numeric"])

        _, ax = plt.subplots()
        p = VectorPlotter(data=long_df, variables={"x": "x"})

        with pytest.raises(TypeError):
            p._attach(ax, allowed_types="categorical")

        _, ax = plt.subplots()
        p = VectorPlotter(data=long_df, variables={"x": "x", "y": "t"})

        with pytest.raises(TypeError):
            p._attach(ax, allowed_types=["numeric", "categorical"])

    def test_attach_log_scale(self, long_df):

        _, ax = plt.subplots()
        p = VectorPlotter(data=long_df, variables={"x": "x"})
        p._attach(ax, log_scale=True)
        assert ax.xaxis.get_scale() == "log"
        assert ax.yaxis.get_scale() == "linear"
        assert p._log_scaled("x")
        assert not p._log_scaled("y")

        _, ax = plt.subplots()
        p = VectorPlotter(data=long_df, variables={"x": "x"})
        p._attach(ax, log_scale=2)
        assert ax.xaxis.get_scale() == "log"
        assert ax.yaxis.get_scale() == "linear"
        assert p._log_scaled("x")
        assert not p._log_scaled("y")

        _, ax = plt.subplots()
        p = VectorPlotter(data=long_df, variables={"y": "y"})
        p._attach(ax, log_scale=True)
        assert ax.xaxis.get_scale() == "linear"
        assert ax.yaxis.get_scale() == "log"
        assert not p._log_scaled("x")
        assert p._log_scaled("y")

        _, ax = plt.subplots()
        p = VectorPlotter(data=long_df, variables={"x": "x", "y": "y"})
        p._attach(ax, log_scale=True)
        assert ax.xaxis.get_scale() == "log"
        assert ax.yaxis.get_scale() == "log"
        assert p._log_scaled("x")
        assert p._log_scaled("y")

        _, ax = plt.subplots()
        p = VectorPlotter(data=long_df, variables={"x": "x", "y": "y"})
        p._attach(ax, log_scale=(True, False))
        assert ax.xaxis.get_scale() == "log"
        assert ax.yaxis.get_scale() == "linear"
        assert p._log_scaled("x")
        assert not p._log_scaled("y")

        _, ax = plt.subplots()
        p = VectorPlotter(data=long_df, variables={"x": "x", "y": "y"})
        p._attach(ax, log_scale=(False, 2))
        assert ax.xaxis.get_scale() == "linear"
        assert ax.yaxis.get_scale() == "log"
        assert not p._log_scaled("x")
        assert p._log_scaled("y")

    def test_attach_converters(self, long_df):

        _, ax = plt.subplots()
        p = VectorPlotter(data=long_df, variables={"x": "x", "y": "t"})
        p._attach(ax)
        assert ax.xaxis.converter is None
        assert isinstance(ax.yaxis.converter, mpl.dates.DateConverter)

        _, ax = plt.subplots()
        p = VectorPlotter(data=long_df, variables={"x": "a", "y": "y"})
        p._attach(ax)
        assert isinstance(ax.xaxis.converter, mpl.category.StrCategoryConverter)
        assert ax.yaxis.converter is None

    def test_attach_facets(self, long_df):

        g = FacetGrid(long_df, col="a")
        p = VectorPlotter(data=long_df, variables={"x": "x", "col": "a"})
        p._attach(g)
        assert p.ax is None
        assert p.facets == g

    def test_get_axes_single(self, long_df):

        ax = plt.figure().subplots()
        p = VectorPlotter(data=long_df, variables={"x": "x", "hue": "a"})
        p._attach(ax)
        assert p._get_axes({"hue": "a"}) is ax

    def test_get_axes_facets(self, long_df):

        g = FacetGrid(long_df, col="a")
        p = VectorPlotter(data=long_df, variables={"x": "x", "col": "a"})
        p._attach(g)
        assert p._get_axes({"col": "b"}) is g.axes_dict["b"]

        g = FacetGrid(long_df, col="a", row="c")
        p = VectorPlotter(
            data=long_df, variables={"x": "x", "col": "a", "row": "c"}
        )
        p._attach(g)
        assert p._get_axes({"row": 1, "col": "b"}) is g.axes_dict[(1, "b")]

    def test_comp_data(self, long_df):

        p = VectorPlotter(data=long_df, variables={"x": "x", "y": "t"})

        # We have disabled this check for now, while it remains part of
        # the internal API, because it will require updating a number of tests
        # with pytest.raises(AttributeError):
        #     p.comp_data

        _, ax = plt.subplots()
        p._attach(ax)

        assert_array_equal(p.comp_data["x"], p.plot_data["x"])
        assert_array_equal(
            p.comp_data["y"], ax.yaxis.convert_units(p.plot_data["y"])
        )

        p = VectorPlotter(data=long_df, variables={"x": "a"})

        _, ax = plt.subplots()
        p._attach(ax)

        assert_array_equal(
            p.comp_data["x"], ax.xaxis.convert_units(p.plot_data["x"])
        )

    def test_comp_data_log(self, long_df):

        p = VectorPlotter(data=long_df, variables={"x": "z", "y": "y"})
        _, ax = plt.subplots()
        p._attach(ax, log_scale=(True, False))

        assert_array_equal(
            p.comp_data["x"], np.log10(p.plot_data["x"])
        )
        assert_array_equal(p.comp_data["y"], p.plot_data["y"])

    def test_comp_data_category_order(self):

        s = (pd.Series(["a", "b", "c", "a"], dtype="category")
             .cat.set_categories(["b", "c", "a"], ordered=True))

        p = VectorPlotter(variables={"x": s})
        _, ax = plt.subplots()
        p._attach(ax)
        assert_array_equal(
            p.comp_data["x"],
            [2, 0, 1, 2],
        )

    @pytest.fixture(
        params=itertools.product(
            [None, np.nan, PD_NA],
            ["numeric", "category", "datetime"]
        )
    )
    @pytest.mark.parametrize(
        "NA,var_type",
    )
    def comp_data_missing_fixture(self, request):

        # This fixture holds the logic for parametrizing
        # the following test (test_comp_data_missing)

        NA, var_type = request.param

        if NA is None:
            pytest.skip("No pandas.NA available")

        comp_data = [0, 1, np.nan, 2, np.nan, 1]
        if var_type == "numeric":
            orig_data = [0, 1, NA, 2, np.inf, 1]
        elif var_type == "category":
            orig_data = ["a", "b", NA, "c", NA, "b"]
        elif var_type == "datetime":
            # Use 1-based numbers to avoid issue on matplotlib<3.2
            # Could simplify the test a bit when we roll off that version
            comp_data = [1, 2, np.nan, 3, np.nan, 2]
            numbers = [1, 2, 3, 2]

            orig_data = mpl.dates.num2date(numbers)
            orig_data.insert(2, NA)
            orig_data.insert(4, np.inf)

        return orig_data, comp_data

    def test_comp_data_missing(self, comp_data_missing_fixture):

        orig_data, comp_data = comp_data_missing_fixture
        p = VectorPlotter(variables={"x": orig_data})
        ax = plt.figure().subplots()
        p._attach(ax)
        assert_array_equal(p.comp_data["x"], comp_data)

    def test_var_order(self, long_df):

        order = ["c", "b", "a"]
        for var in ["hue", "size", "style"]:
            p = VectorPlotter(data=long_df, variables={"x": "x", var: "a"})

            mapper = getattr(p, f"map_{var}")
            mapper(order=order)

            assert p.var_levels[var] == order


class TestCoreFunc:

    def test_unique_dashes(self):

        n = 24
        dashes = unique_dashes(n)

        assert len(dashes) == n
        assert len(set(dashes)) == n
        assert dashes[0] == ""
        for spec in dashes[1:]:
            assert isinstance(spec, tuple)
            assert not len(spec) % 2

    def test_unique_markers(self):

        n = 24
        markers = unique_markers(n)

        assert len(markers) == n
        assert len(set(markers)) == n
        for m in markers:
            assert mpl.markers.MarkerStyle(m).is_filled()

    def test_variable_type(self):

        s = pd.Series([1., 2., 3.])
        assert variable_type(s) == "numeric"
        assert variable_type(s.astype(int)) == "numeric"
        assert variable_type(s.astype(object)) == "numeric"
        # assert variable_type(s.to_numpy()) == "numeric"
        assert variable_type(s.values) == "numeric"
        # assert variable_type(s.to_list()) == "numeric"
        assert variable_type(s.tolist()) == "numeric"

        s = pd.Series([1, 2, 3, np.nan], dtype=object)
        assert variable_type(s) == "numeric"

        s = pd.Series([np.nan, np.nan])
        # s = pd.Series([pd.NA, pd.NA])
        assert variable_type(s) == "numeric"

        s = pd.Series(["1", "2", "3"])
        assert variable_type(s) == "categorical"
        # assert variable_type(s.to_numpy()) == "categorical"
        assert variable_type(s.values) == "categorical"
        # assert variable_type(s.to_list()) == "categorical"
        assert variable_type(s.tolist()) == "categorical"

        s = pd.Series([True, False, False])
        assert variable_type(s) == "numeric"
        assert variable_type(s, boolean_type="categorical") == "categorical"
        s_cat = s.astype("category")
        assert variable_type(s_cat, boolean_type="categorical") == "categorical"
        assert variable_type(s_cat, boolean_type="numeric") == "categorical"

        s = pd.Series([pd.Timestamp(1), pd.Timestamp(2)])
        assert variable_type(s) == "datetime"
        assert variable_type(s.astype(object)) == "datetime"
        # assert variable_type(s.to_numpy()) == "datetime"
        assert variable_type(s.values) == "datetime"
        # assert variable_type(s.to_list()) == "datetime"
        assert variable_type(s.tolist()) == "datetime"

    def test_infer_orient(self):

        nums = pd.Series(np.arange(6))
        cats = pd.Series(["a", "b"] * 3)

        assert infer_orient(cats, nums) == "v"
        assert infer_orient(nums, cats) == "h"

        assert infer_orient(nums, None) == "h"
        with pytest.warns(UserWarning, match="Vertical .+ `x`"):
            assert infer_orient(nums, None, "v") == "h"

        assert infer_orient(None, nums) == "v"
        with pytest.warns(UserWarning, match="Horizontal .+ `y`"):
            assert infer_orient(None, nums, "h") == "v"

        infer_orient(cats, None, require_numeric=False) == "h"
        with pytest.raises(TypeError, match="Horizontal .+ `x`"):
            infer_orient(cats, None)

        infer_orient(cats, None, require_numeric=False) == "v"
        with pytest.raises(TypeError, match="Vertical .+ `y`"):
            infer_orient(None, cats)

        assert infer_orient(nums, nums, "vert") == "v"
        assert infer_orient(nums, nums, "hori") == "h"

        assert infer_orient(cats, cats, "h", require_numeric=False) == "h"
        assert infer_orient(cats, cats, "v", require_numeric=False) == "v"
        assert infer_orient(cats, cats, require_numeric=False) == "v"

        with pytest.raises(TypeError, match="Vertical .+ `y`"):
            infer_orient(cats, cats, "v")
        with pytest.raises(TypeError, match="Horizontal .+ `x`"):
            infer_orient(cats, cats, "h")
        with pytest.raises(TypeError, match="Neither"):
            infer_orient(cats, cats)

    def test_categorical_order(self):

        x = ["a", "c", "c", "b", "a", "d"]
        y = [3, 2, 5, 1, 4]
        order = ["a", "b", "c", "d"]

        out = categorical_order(x)
        assert out == ["a", "c", "b", "d"]

        out = categorical_order(x, order)
        assert out == order

        out = categorical_order(x, ["b", "a"])
        assert out == ["b", "a"]

        out = categorical_order(np.array(x))
        assert out == ["a", "c", "b", "d"]

        out = categorical_order(pd.Series(x))
        assert out == ["a", "c", "b", "d"]

        out = categorical_order(y)
        assert out == [1, 2, 3, 4, 5]

        out = categorical_order(np.array(y))
        assert out == [1, 2, 3, 4, 5]

        out = categorical_order(pd.Series(y))
        assert out == [1, 2, 3, 4, 5]

        x = pd.Categorical(x, order)
        out = categorical_order(x)
        assert out == list(x.categories)

        x = pd.Series(x)
        out = categorical_order(x)
        assert out == list(x.cat.categories)

        out = categorical_order(x, ["b", "a"])
        assert out == ["b", "a"]

        x = ["a", np.nan, "c", "c", "b", "a", "d"]
        out = categorical_order(x)
        assert out == ["a", "c", "b", "d"]
