from distutils.version import LooseVersion

import pytest
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy.testing as npt

from .. import rcmod, palettes, utils
from ..conftest import has_verdana


class RCParamTester:

    def flatten_list(self, orig_list):

        iter_list = map(np.atleast_1d, orig_list)
        flat_list = [item for sublist in iter_list for item in sublist]
        return flat_list

    def assert_rc_params(self, params):

        for k, v in params.items():
            # Various subtle issues in matplotlib lead to unexpected
            # values for the backend rcParam, which isn't relevant here
            if k == "backend":
                continue
            if isinstance(v, np.ndarray):
                npt.assert_array_equal(mpl.rcParams[k], v)
            else:
                assert mpl.rcParams[k] == v

    def assert_rc_params_equal(self, params1, params2):

        for key, v1 in params1.items():
            # Various subtle issues in matplotlib lead to unexpected
            # values for the backend rcParam, which isn't relevant here
            if key == "backend":
                continue

            v2 = params2[key]
            if isinstance(v1, np.ndarray):
                npt.assert_array_equal(v1, v2)
            else:
                assert v1 == v2


class TestAxesStyle(RCParamTester):

    styles = ["white", "dark", "whitegrid", "darkgrid", "ticks"]

    def test_default_return(self):

        current = rcmod.axes_style()
        self.assert_rc_params(current)

    def test_key_usage(self):

        _style_keys = set(rcmod._style_keys)
        for style in self.styles:
            assert not set(rcmod.axes_style(style)) ^ _style_keys

    def test_bad_style(self):

        with pytest.raises(ValueError):
            rcmod.axes_style("i_am_not_a_style")

    def test_rc_override(self):

        rc = {"axes.facecolor": "blue", "foo.notaparam": "bar"}
        out = rcmod.axes_style("darkgrid", rc)
        assert out["axes.facecolor"] == "blue"
        assert "foo.notaparam" not in out

    def test_set_style(self):

        for style in self.styles:

            style_dict = rcmod.axes_style(style)
            rcmod.set_style(style)
            self.assert_rc_params(style_dict)

    def test_style_context_manager(self):

        rcmod.set_style("darkgrid")
        orig_params = rcmod.axes_style()
        context_params = rcmod.axes_style("whitegrid")

        with rcmod.axes_style("whitegrid"):
            self.assert_rc_params(context_params)
        self.assert_rc_params(orig_params)

        @rcmod.axes_style("whitegrid")
        def func():
            self.assert_rc_params(context_params)
        func()
        self.assert_rc_params(orig_params)

    def test_style_context_independence(self):

        assert set(rcmod._style_keys) ^ set(rcmod._context_keys)

    def test_set_rc(self):

        rcmod.set_theme(rc={"lines.linewidth": 4})
        assert mpl.rcParams["lines.linewidth"] == 4
        rcmod.set_theme()

    def test_set_with_palette(self):

        rcmod.reset_orig()

        rcmod.set_theme(palette="deep")
        assert utils.get_color_cycle() == palettes.color_palette("deep", 10)
        rcmod.reset_orig()

        rcmod.set_theme(palette="deep", color_codes=False)
        assert utils.get_color_cycle() == palettes.color_palette("deep", 10)
        rcmod.reset_orig()

        pal = palettes.color_palette("deep")
        rcmod.set_theme(palette=pal)
        assert utils.get_color_cycle() == palettes.color_palette("deep", 10)
        rcmod.reset_orig()

        rcmod.set_theme(palette=pal, color_codes=False)
        assert utils.get_color_cycle() == palettes.color_palette("deep", 10)
        rcmod.reset_orig()

        rcmod.set_theme()

    def test_reset_defaults(self):

        rcmod.reset_defaults()
        self.assert_rc_params(mpl.rcParamsDefault)
        rcmod.set_theme()

    def test_reset_orig(self):

        rcmod.reset_orig()
        self.assert_rc_params(mpl.rcParamsOrig)
        rcmod.set_theme()

    def test_set_is_alias(self):

        rcmod.set_theme(context="paper", style="white")
        params1 = mpl.rcParams.copy()
        rcmod.reset_orig()

        rcmod.set_theme(context="paper", style="white")
        params2 = mpl.rcParams.copy()

        self.assert_rc_params_equal(params1, params2)

        rcmod.set_theme()


class TestPlottingContext(RCParamTester):

    contexts = ["paper", "notebook", "talk", "poster"]

    def test_default_return(self):

        current = rcmod.plotting_context()
        self.assert_rc_params(current)

    def test_key_usage(self):

        _context_keys = set(rcmod._context_keys)
        for context in self.contexts:
            missing = set(rcmod.plotting_context(context)) ^ _context_keys
            assert not missing

    def test_bad_context(self):

        with pytest.raises(ValueError):
            rcmod.plotting_context("i_am_not_a_context")

    def test_font_scale(self):

        notebook_ref = rcmod.plotting_context("notebook")
        notebook_big = rcmod.plotting_context("notebook", 2)

        font_keys = ["axes.labelsize", "axes.titlesize", "legend.fontsize",
                     "xtick.labelsize", "ytick.labelsize", "font.size"]

        if LooseVersion(mpl.__version__) >= "3.0":
            font_keys.append("legend.title_fontsize")

        for k in font_keys:
            assert notebook_ref[k] * 2 == notebook_big[k]

    def test_rc_override(self):

        key, val = "grid.linewidth", 5
        rc = {key: val, "foo": "bar"}
        out = rcmod.plotting_context("talk", rc=rc)
        assert out[key] == val
        assert "foo" not in out

    def test_set_context(self):

        for context in self.contexts:

            context_dict = rcmod.plotting_context(context)
            rcmod.set_context(context)
            self.assert_rc_params(context_dict)

    def test_context_context_manager(self):

        rcmod.set_context("notebook")
        orig_params = rcmod.plotting_context()
        context_params = rcmod.plotting_context("paper")

        with rcmod.plotting_context("paper"):
            self.assert_rc_params(context_params)
        self.assert_rc_params(orig_params)

        @rcmod.plotting_context("paper")
        def func():
            self.assert_rc_params(context_params)
        func()
        self.assert_rc_params(orig_params)


class TestPalette:

    def test_set_palette(self):

        rcmod.set_palette("deep")
        assert utils.get_color_cycle() == palettes.color_palette("deep", 10)

        rcmod.set_palette("pastel6")
        assert utils.get_color_cycle() == palettes.color_palette("pastel6", 6)

        rcmod.set_palette("dark", 4)
        assert utils.get_color_cycle() == palettes.color_palette("dark", 4)

        rcmod.set_palette("Set2", color_codes=True)
        assert utils.get_color_cycle() == palettes.color_palette("Set2", 8)


class TestFonts:

    _no_verdana = not has_verdana()

    @pytest.mark.skipif(_no_verdana, reason="Verdana font is not present")
    def test_set_font(self):

        rcmod.set_theme(font="Verdana")

        _, ax = plt.subplots()
        ax.set_xlabel("foo")

        assert ax.xaxis.label.get_fontname() == "Verdana"

        rcmod.set_theme()

    def test_set_serif_font(self):

        rcmod.set_theme(font="serif")

        _, ax = plt.subplots()
        ax.set_xlabel("foo")

        assert ax.xaxis.label.get_fontname() in mpl.rcParams["font.serif"]

        rcmod.set_theme()

    @pytest.mark.skipif(_no_verdana, reason="Verdana font is not present")
    def test_different_sans_serif(self):

        rcmod.set_theme()
        rcmod.set_style(rc={"font.sans-serif": ["Verdana"]})

        _, ax = plt.subplots()
        ax.set_xlabel("foo")

        assert ax.xaxis.label.get_fontname() == "Verdana"

        rcmod.set_theme()
