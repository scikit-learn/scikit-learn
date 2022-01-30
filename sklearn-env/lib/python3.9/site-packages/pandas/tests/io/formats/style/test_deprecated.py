"""
modules collects tests for Styler methods which have been deprecated
"""
import numpy as np
import pytest

jinja2 = pytest.importorskip("jinja2")

from pandas import (
    DataFrame,
    IndexSlice,
    NaT,
    Timestamp,
)
import pandas._testing as tm


@pytest.fixture
def df():
    return DataFrame({"A": [0, 1], "B": np.random.randn(2)})


@pytest.mark.parametrize("axis", ["index", "columns"])
def test_hide_index_columns(df, axis):
    with tm.assert_produces_warning(FutureWarning):
        getattr(df.style, "hide_" + axis)()


def test_set_non_numeric_na():
    # GH 21527 28358
    df = DataFrame(
        {
            "object": [None, np.nan, "foo"],
            "datetime": [None, NaT, Timestamp("20120101")],
        }
    )

    with tm.assert_produces_warning(FutureWarning):
        ctx = df.style.set_na_rep("NA")._translate(True, True)
    assert ctx["body"][0][1]["display_value"] == "NA"
    assert ctx["body"][0][2]["display_value"] == "NA"
    assert ctx["body"][1][1]["display_value"] == "NA"
    assert ctx["body"][1][2]["display_value"] == "NA"


def test_where_with_one_style(df):
    # GH 17474
    def f(x):
        return x > 0.5

    style1 = "foo: bar"

    with tm.assert_produces_warning(FutureWarning):
        result = df.style.where(f, style1)._compute().ctx
    expected = {
        (r, c): [("foo", "bar")]
        for r, row in enumerate(df.index)
        for c, col in enumerate(df.columns)
        if f(df.loc[row, col])
    }
    assert result == expected


@pytest.mark.parametrize(
    "slice_",
    [
        IndexSlice[:],
        IndexSlice[:, ["A"]],
        IndexSlice[[1], :],
        IndexSlice[[1], ["A"]],
        IndexSlice[:2, ["A", "B"]],
    ],
)
def test_where_subset(df, slice_):
    # GH 17474
    def f(x):
        return x > 0.5

    style1 = "foo: bar"
    style2 = "baz: foo"

    with tm.assert_produces_warning(FutureWarning):
        res = df.style.where(f, style1, style2, subset=slice_)._compute().ctx
    expected = {
        (r, c): [("foo", "bar") if f(df.loc[row, col]) else ("baz", "foo")]
        for r, row in enumerate(df.index)
        for c, col in enumerate(df.columns)
        if row in df.loc[slice_].index and col in df.loc[slice_].columns
    }
    assert res == expected


def test_where_subset_compare_with_applymap(df):
    # GH 17474
    def f(x):
        return x > 0.5

    style1 = "foo: bar"
    style2 = "baz: foo"

    def g(x):
        return style1 if f(x) else style2

    slices = [
        IndexSlice[:],
        IndexSlice[:, ["A"]],
        IndexSlice[[1], :],
        IndexSlice[[1], ["A"]],
        IndexSlice[:2, ["A", "B"]],
    ]

    for slice_ in slices:
        with tm.assert_produces_warning(FutureWarning):
            result = df.style.where(f, style1, style2, subset=slice_)._compute().ctx
        expected = df.style.applymap(g, subset=slice_)._compute().ctx
        assert result == expected


def test_where_kwargs():
    df = DataFrame([[1, 2], [3, 4]])

    def f(x, val):
        return x > val

    with tm.assert_produces_warning(FutureWarning):
        res = df.style.where(f, "color:green;", "color:red;", val=2)._compute().ctx
    expected = {
        (0, 0): [("color", "red")],
        (0, 1): [("color", "red")],
        (1, 0): [("color", "green")],
        (1, 1): [("color", "green")],
    }
    assert res == expected


def test_set_na_rep():
    # GH 21527 28358
    df = DataFrame([[None, None], [1.1, 1.2]], columns=["A", "B"])

    with tm.assert_produces_warning(FutureWarning):
        ctx = df.style.set_na_rep("NA")._translate(True, True)
    assert ctx["body"][0][1]["display_value"] == "NA"
    assert ctx["body"][0][2]["display_value"] == "NA"

    with tm.assert_produces_warning(FutureWarning):
        ctx = (
            df.style.set_na_rep("NA")
            .format(None, na_rep="-", subset=["B"])
            ._translate(True, True)
        )
    assert ctx["body"][0][1]["display_value"] == "NA"
    assert ctx["body"][0][2]["display_value"] == "-"


def test_precision(df):
    styler = df.style
    with tm.assert_produces_warning(FutureWarning):
        s2 = styler.set_precision(1)
    assert styler is s2
    assert styler.precision == 1


def test_render(df):
    with tm.assert_produces_warning(FutureWarning):
        df.style.render()
