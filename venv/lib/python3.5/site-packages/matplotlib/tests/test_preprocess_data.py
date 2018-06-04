from __future__ import (absolute_import, division, print_function)

import re

import numpy as np
import pytest

from matplotlib import _preprocess_data


# Notes on testing the plotting functions itself
# *   the individual decorated plotting functions are tested in 'test_axes.py'
# *   that pyplot functions accept a data kwarg is only tested in
#     test_axes.test_pie_linewidth_0


# these two get used in multiple tests, so define them here
@_preprocess_data(replace_names=["x", "y"], label_namer="y")
def plot_func(ax, x, y, ls="x", label=None, w="xyz"):
    return ("x: %s, y: %s, ls: %s, w: %s, label: %s" % (
        list(x), list(y), ls, w, label))


@_preprocess_data(replace_names=["x", "y"], label_namer="y",
                  positional_parameter_names=["x", "y", "ls", "label", "w"])
def plot_func_varargs(ax, *args, **kwargs):
    all_args = [None, None, "x", None, "xyz"]
    for i, v in enumerate(args):
        all_args[i] = v
    for i, k in enumerate(["x", "y", "ls", "label", "w"]):
        if k in kwargs:
            all_args[i] = kwargs[k]
    x, y, ls, label, w = all_args
    return ("x: %s, y: %s, ls: %s, w: %s, label: %s" % (
        list(x), list(y), ls, w, label))


all_funcs = [plot_func, plot_func_varargs]
all_func_ids = ['plot_func', 'plot_func_varargs']


def test_compiletime_checks():
    """test decorator invocations -> no replacements"""

    def func(ax, x, y): pass

    def func_args(ax, x, y, *args): pass

    def func_kwargs(ax, x, y, **kwargs): pass

    def func_no_ax_args(*args, **kwargs): pass

    # this is ok
    _preprocess_data(replace_names=["x", "y"])(func)
    _preprocess_data(replace_names=["x", "y"])(func_kwargs)
    # this has "enough" information to do all the replaces
    _preprocess_data(replace_names=["x", "y"])(func_args)

    # no positional_parameter_names but needed due to replaces
    with pytest.raises(AssertionError):
        # z is unknown
        _preprocess_data(replace_names=["x", "y", "z"])(func_args)

    with pytest.raises(AssertionError):
        _preprocess_data(replace_names=["x", "y"])(func_no_ax_args)

    # no replacements at all -> all ok...
    _preprocess_data(replace_names=[], label_namer=None)(func)
    _preprocess_data(replace_names=[], label_namer=None)(func_args)
    _preprocess_data(replace_names=[], label_namer=None)(func_kwargs)
    _preprocess_data(replace_names=[], label_namer=None)(func_no_ax_args)

    # label namer is unknown
    with pytest.raises(AssertionError):
        _preprocess_data(label_namer="z")(func)

    with pytest.raises(AssertionError):
        _preprocess_data(label_namer="z")(func_args)

    # but "ok-ish", if func has kwargs -> will show up at runtime :-(
    _preprocess_data(label_namer="z")(func_kwargs)
    _preprocess_data(label_namer="z")(func_no_ax_args)


def test_label_problems_at_runtime():
    """Tests for behaviour which would actually be nice to get rid of."""

    @_preprocess_data(label_namer="z")
    def func(*args, **kwargs):
        pass

    # This is a programming mistake: the parameter which should add the
    # label is not present in the function call. Unfortunately this was masked
    # due to the **kwargs usage
    # This would be nice to handle as a compiletime check (see above...)
    with pytest.warns(RuntimeWarning):
        func(None, x="a", y="b")

    def real_func(x, y):
        pass

    @_preprocess_data(label_namer="x")
    def func(*args, **kwargs):
        real_func(**kwargs)

    # This sets a label although the function can't handle it.
    with pytest.raises(TypeError):
        func(None, x="a", y="b")


@pytest.mark.parametrize('func', all_funcs, ids=all_func_ids)
def test_function_call_without_data(func):
    """test without data -> no replacements"""
    assert (func(None, "x", "y") ==
            "x: ['x'], y: ['y'], ls: x, w: xyz, label: None")
    assert (func(None, x="x", y="y") ==
            "x: ['x'], y: ['y'], ls: x, w: xyz, label: None")
    assert (func(None, "x", "y", label="") ==
            "x: ['x'], y: ['y'], ls: x, w: xyz, label: ")
    assert (func(None, "x", "y", label="text") ==
            "x: ['x'], y: ['y'], ls: x, w: xyz, label: text")
    assert (func(None, x="x", y="y", label="") ==
            "x: ['x'], y: ['y'], ls: x, w: xyz, label: ")
    assert (func(None, x="x", y="y", label="text") ==
            "x: ['x'], y: ['y'], ls: x, w: xyz, label: text")


@pytest.mark.parametrize('func', all_funcs, ids=all_func_ids)
def test_function_call_with_dict_data(func):
    """Test with dict data -> label comes from the value of 'x' parameter """
    data = {"a": [1, 2], "b": [8, 9], "w": "NOT"}
    assert (func(None, "a", "b", data=data) ==
            "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: b")
    assert (func(None, x="a", y="b", data=data) ==
            "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: b")
    assert (func(None, "a", "b", label="", data=data) ==
            "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: ")
    assert (func(None, "a", "b", label="text", data=data) ==
            "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: text")
    assert (func(None, x="a", y="b", label="", data=data) ==
            "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: ")
    assert (func(None, x="a", y="b", label="text", data=data) ==
            "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: text")


@pytest.mark.parametrize('func', all_funcs, ids=all_func_ids)
def test_function_call_with_dict_data_not_in_data(func):
    "test for the case that one var is not in data -> half replaces, half kept"
    data = {"a": [1, 2], "w": "NOT"}
    assert (func(None, "a", "b", data=data) ==
            "x: [1, 2], y: ['b'], ls: x, w: xyz, label: b")
    assert (func(None, x="a", y="b", data=data) ==
            "x: [1, 2], y: ['b'], ls: x, w: xyz, label: b")
    assert (func(None, "a", "b", label="", data=data) ==
            "x: [1, 2], y: ['b'], ls: x, w: xyz, label: ")
    assert (func(None, "a", "b", label="text", data=data) ==
            "x: [1, 2], y: ['b'], ls: x, w: xyz, label: text")
    assert (func(None, x="a", y="b", label="", data=data) ==
            "x: [1, 2], y: ['b'], ls: x, w: xyz, label: ")
    assert (func(None, x="a", y="b", label="text", data=data) ==
            "x: [1, 2], y: ['b'], ls: x, w: xyz, label: text")


@pytest.mark.parametrize('func', all_funcs, ids=all_func_ids)
def test_function_call_with_pandas_data(func, pd):
    """test with pandas dataframe -> label comes from data["col"].name """
    data = pd.DataFrame({"a": np.array([1, 2], dtype=np.int32),
                         "b": np.array([8, 9], dtype=np.int32),
                         "w": ["NOT", "NOT"]})

    assert (func(None, "a", "b", data=data) ==
            "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: b")
    assert (func(None, x="a", y="b", data=data) ==
            "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: b")
    assert (func(None, "a", "b", label="", data=data) ==
            "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: ")
    assert (func(None, "a", "b", label="text", data=data) ==
            "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: text")
    assert (func(None, x="a", y="b", label="", data=data) ==
            "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: ")
    assert (func(None, x="a", y="b", label="text", data=data) ==
            "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: text")


def test_function_call_replace_all():
    """Test without a "replace_names" argument, all vars should be replaced"""
    data = {"a": [1, 2], "b": [8, 9], "x": "xyz"}

    @_preprocess_data(label_namer="y")
    def func_replace_all(ax, x, y, ls="x", label=None, w="NOT"):
        return "x: %s, y: %s, ls: %s, w: %s, label: %s" % (
            list(x), list(y), ls, w, label)

    assert (func_replace_all(None, "a", "b", w="x", data=data) ==
            "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: b")
    assert (func_replace_all(None, x="a", y="b", w="x", data=data) ==
            "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: b")
    assert (func_replace_all(None, "a", "b", w="x", label="", data=data) ==
            "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: ")
    assert (
        func_replace_all(None, "a", "b", w="x", label="text", data=data) ==
        "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: text")
    assert (
        func_replace_all(None, x="a", y="b", w="x", label="", data=data) ==
        "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: ")
    assert (
        func_replace_all(None, x="a", y="b", w="x", label="text", data=data) ==
        "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: text")

    @_preprocess_data(label_namer="y")
    def func_varags_replace_all(ax, *args, **kwargs):
        all_args = [None, None, "x", None, "xyz"]
        for i, v in enumerate(args):
            all_args[i] = v
        for i, k in enumerate(["x", "y", "ls", "label", "w"]):
            if k in kwargs:
                all_args[i] = kwargs[k]
        x, y, ls, label, w = all_args
        return "x: %s, y: %s, ls: %s, w: %s, label: %s" % (
            list(x), list(y), ls, w, label)

    # in the first case, we can't get a "y" argument,
    # as we don't know the names of the *args
    assert (func_varags_replace_all(None, x="a", y="b", w="x", data=data) ==
            "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: b")
    assert (
        func_varags_replace_all(None, "a", "b", w="x", label="", data=data) ==
        "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: ")
    assert (
        func_varags_replace_all(None, "a", "b", w="x", label="text",
                                data=data) ==
        "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: text")
    assert (
        func_varags_replace_all(None, x="a", y="b", w="x", label="",
                                data=data) ==
        "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: ")
    assert (
        func_varags_replace_all(None, x="a", y="b", w="x", label="text",
                                data=data) ==
        "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: text")

    with pytest.warns(RuntimeWarning):
        assert (func_varags_replace_all(None, "a", "b", w="x", data=data) ==
                "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: None")


def test_no_label_replacements():
    """Test with "label_namer=None" -> no label replacement at all"""

    @_preprocess_data(replace_names=["x", "y"], label_namer=None)
    def func_no_label(ax, x, y, ls="x", label=None, w="xyz"):
        return "x: %s, y: %s, ls: %s, w: %s, label: %s" % (
            list(x), list(y), ls, w, label)

    data = {"a": [1, 2], "b": [8, 9], "w": "NOT"}
    assert (func_no_label(None, "a", "b", data=data) ==
            "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: None")
    assert (func_no_label(None, x="a", y="b", data=data) ==
            "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: None")
    assert (func_no_label(None, "a", "b", label="", data=data) ==
            "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: ")
    assert (func_no_label(None, "a", "b", label="text", data=data) ==
            "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: text")


def test_more_args_than_pos_parameter():
    @_preprocess_data(replace_names=["x", "y"], label_namer="y")
    def func(ax, x, y, z=1):
        pass

    data = {"a": [1, 2], "b": [8, 9], "w": "NOT"}
    with pytest.raises(RuntimeError):
        func(None, "a", "b", "z", "z", data=data)


def test_function_call_with_replace_all_args():
    """Test with a "replace_all_args" argument, all *args should be replaced"""
    data = {"a": [1, 2], "b": [8, 9], "x": "xyz"}

    def funcy(ax, *args, **kwargs):
        all_args = [None, None, "x", None, "NOT"]
        for i, v in enumerate(args):
            all_args[i] = v
        for i, k in enumerate(["x", "y", "ls", "label", "w"]):
            if k in kwargs:
                all_args[i] = kwargs[k]
        x, y, ls, label, w = all_args
        return "x: %s, y: %s, ls: %s, w: %s, label: %s" % (
            list(x), list(y), ls, w, label)

    func = _preprocess_data(replace_all_args=True, replace_names=["w"],
                            label_namer="y")(funcy)

    assert (func(None, "a", "b", w="x", label="", data=data) ==
            "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: ")
    assert (func(None, "a", "b", w="x", label="text", data=data) ==
            "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: text")

    func2 = _preprocess_data(replace_all_args=True, replace_names=["w"],
                             label_namer="y",
                             positional_parameter_names=["x", "y", "ls",
                                                         "label", "w"])(funcy)

    assert (func2(None, "a", "b", w="x", data=data) ==
            "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: b")
    assert (func2(None, "a", "b", w="x", label="", data=data) ==
            "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: ")
    assert (func2(None, "a", "b", w="x", label="text", data=data) ==
            "x: [1, 2], y: [8, 9], ls: x, w: xyz, label: text")


def test_docstring_addition():
    @_preprocess_data()
    def funcy(ax, *args, **kwargs):
        """Funcy does nothing"""
        pass

    assert re.search(r".*All positional and all keyword arguments\.",
                     funcy.__doc__)
    assert not re.search(r".*All positional arguments\.", funcy.__doc__)
    assert not re.search(r".*All arguments with the following names: .*",
                         funcy.__doc__)

    @_preprocess_data(replace_all_args=True, replace_names=[])
    def funcy(ax, x, y, z, bar=None):
        """Funcy does nothing"""
        pass

    assert re.search(r".*All positional arguments\.",
                     funcy.__doc__)
    assert not re.search(r".*All positional and all keyword arguments\.",
                         funcy.__doc__)
    assert not re.search(r".*All arguments with the following names: .*",
                         funcy.__doc__)

    @_preprocess_data(replace_all_args=True, replace_names=["bar"])
    def funcy(ax, x, y, z, bar=None):
        """Funcy does nothing"""
        pass

    assert re.search(r".*All positional arguments\.", funcy.__doc__)
    assert re.search(r".*All arguments with the following names: 'bar'\.",
                     funcy.__doc__)
    assert not re.search(r".*All positional and all keyword arguments\.",
                         funcy.__doc__)

    @_preprocess_data(replace_names=["x", "bar"])
    def funcy(ax, x, y, z, bar=None):
        """Funcy does nothing"""
        pass

    # lists can print in any order, so test for both x,bar and bar,x
    assert re.search(r".*All arguments with the following names: '.*', '.*'\.",
                     funcy.__doc__)
    assert re.search(r".*'x'.*", funcy.__doc__)
    assert re.search(r".*'bar'.*", funcy.__doc__)
    assert not re.search(r".*All positional and all keyword arguments\.",
                         funcy.__doc__)
    assert not re.search(r".*All positional arguments\.",
                         funcy.__doc__)


def test_positional_parameter_names_as_function():
    # Also test the _plot_arg_replacer for plot...
    from matplotlib.axes._axes import _plot_args_replacer

    @_preprocess_data(replace_names=["x", "y"],
                      positional_parameter_names=_plot_args_replacer)
    def funcy(ax, *args, **kwargs):
        return "{args} | {kwargs}".format(args=args, kwargs=kwargs)

    # the normal case...
    data = {"x": "X", "hy1": "Y"}
    assert funcy(None, "x", "hy1", data=data) == "('X', 'Y') | {}"
    assert funcy(None, "x", "hy1", "c", data=data) == "('X', 'Y', 'c') | {}"

    # no arbitrary long args with data
    with pytest.raises(ValueError):
        assert (funcy(None, "x", "y", "c", "x", "y", "x", "y", data=data) ==
                "('X', 'Y', 'c', 'X', 'Y', 'X', 'Y') | {}")

    # In the two arg case, if a valid color spec is in data, we warn but use
    # it as data...
    data = {"x": "X", "y": "Y", "ro": "!!"}
    with pytest.warns(RuntimeWarning):
        assert funcy(None, "y", "ro", data=data) == "('Y', '!!') | {}"
