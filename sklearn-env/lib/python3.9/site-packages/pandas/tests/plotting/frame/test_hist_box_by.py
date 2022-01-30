import re

import numpy as np
import pytest

import pandas.util._test_decorators as td

from pandas import DataFrame
import pandas._testing as tm
from pandas.tests.plotting.common import (
    TestPlotBase,
    _check_plot_works,
)


def _create_hist_box_with_by_df():
    np.random.seed(0)
    df = DataFrame(np.random.randn(30, 2), columns=["A", "B"])
    df["C"] = np.random.choice(["a", "b", "c"], 30)
    df["D"] = np.random.choice(["a", "b", "c"], 30)
    return df


@td.skip_if_no_mpl
class TestHistWithBy(TestPlotBase):
    def setup_method(self, method):
        TestPlotBase.setup_method(self, method)
        import matplotlib as mpl

        mpl.rcdefaults()
        self.hist_df = _create_hist_box_with_by_df()

    @pytest.mark.parametrize(
        "by, column, titles, legends",
        [
            ("C", "A", ["a", "b", "c"], [["A"]] * 3),
            ("C", ["A", "B"], ["a", "b", "c"], [["A", "B"]] * 3),
            ("C", None, ["a", "b", "c"], [["A", "B"]] * 3),
            (
                ["C", "D"],
                "A",
                [
                    "(a, a)",
                    "(a, b)",
                    "(a, c)",
                    "(b, a)",
                    "(b, b)",
                    "(b, c)",
                    "(c, a)",
                    "(c, b)",
                    "(c, c)",
                ],
                [["A"]] * 9,
            ),
            (
                ["C", "D"],
                ["A", "B"],
                [
                    "(a, a)",
                    "(a, b)",
                    "(a, c)",
                    "(b, a)",
                    "(b, b)",
                    "(b, c)",
                    "(c, a)",
                    "(c, b)",
                    "(c, c)",
                ],
                [["A", "B"]] * 9,
            ),
            (
                ["C", "D"],
                None,
                [
                    "(a, a)",
                    "(a, b)",
                    "(a, c)",
                    "(b, a)",
                    "(b, b)",
                    "(b, c)",
                    "(c, a)",
                    "(c, b)",
                    "(c, c)",
                ],
                [["A", "B"]] * 9,
            ),
        ],
    )
    def test_hist_plot_by_argument(self, by, column, titles, legends):
        # GH 15079
        axes = _check_plot_works(self.hist_df.plot.hist, column=column, by=by)
        result_titles = [ax.get_title() for ax in axes]
        result_legends = [
            [legend.get_text() for legend in ax.get_legend().texts] for ax in axes
        ]

        assert result_legends == legends
        assert result_titles == titles

    @pytest.mark.parametrize(
        "by, column, titles, legends",
        [
            (0, "A", ["a", "b", "c"], [["A"]] * 3),
            (0, None, ["a", "b", "c"], [["A", "B"]] * 3),
            (
                [0, "D"],
                "A",
                [
                    "(a, a)",
                    "(a, b)",
                    "(a, c)",
                    "(b, a)",
                    "(b, b)",
                    "(b, c)",
                    "(c, a)",
                    "(c, b)",
                    "(c, c)",
                ],
                [["A"]] * 9,
            ),
        ],
    )
    def test_hist_plot_by_0(self, by, column, titles, legends):
        # GH 15079
        df = self.hist_df.copy()
        df = df.rename(columns={"C": 0})

        axes = _check_plot_works(df.plot.hist, column=column, by=by)
        result_titles = [ax.get_title() for ax in axes]
        result_legends = [
            [legend.get_text() for legend in ax.get_legend().texts] for ax in axes
        ]

        assert result_legends == legends
        assert result_titles == titles

    @pytest.mark.parametrize(
        "by, column",
        [
            ([], ["A"]),
            ([], ["A", "B"]),
            ((), None),
            ((), ["A", "B"]),
        ],
    )
    def test_hist_plot_empty_list_string_tuple_by(self, by, column):
        # GH 15079
        msg = "No group keys passed"
        with pytest.raises(ValueError, match=msg):
            _check_plot_works(self.hist_df.plot.hist, column=column, by=by)

    @pytest.mark.slow
    @pytest.mark.parametrize(
        "by, column, layout, axes_num",
        [
            (["C"], "A", (2, 2), 3),
            ("C", "A", (2, 2), 3),
            (["C"], ["A"], (1, 3), 3),
            ("C", None, (3, 1), 3),
            ("C", ["A", "B"], (3, 1), 3),
            (["C", "D"], "A", (9, 1), 9),
            (["C", "D"], "A", (3, 3), 9),
            (["C", "D"], ["A"], (5, 2), 9),
            (["C", "D"], ["A", "B"], (9, 1), 9),
            (["C", "D"], None, (9, 1), 9),
            (["C", "D"], ["A", "B"], (5, 2), 9),
        ],
    )
    def test_hist_plot_layout_with_by(self, by, column, layout, axes_num):
        # GH 15079
        # _check_plot_works adds an ax so catch warning. see GH #13188
        with tm.assert_produces_warning(UserWarning):
            axes = _check_plot_works(
                self.hist_df.plot.hist, column=column, by=by, layout=layout
            )
        self._check_axes_shape(axes, axes_num=axes_num, layout=layout)

    @pytest.mark.parametrize(
        "msg, by, layout",
        [
            ("larger than required size", ["C", "D"], (1, 1)),
            (re.escape("Layout must be a tuple of (rows, columns)"), "C", (1,)),
            ("At least one dimension of layout must be positive", "C", (-1, -1)),
        ],
    )
    def test_hist_plot_invalid_layout_with_by_raises(self, msg, by, layout):
        # GH 15079, test if error is raised when invalid layout is given

        with pytest.raises(ValueError, match=msg):
            self.hist_df.plot.hist(column=["A", "B"], by=by, layout=layout)

    @pytest.mark.slow
    def test_axis_share_x_with_by(self):
        # GH 15079
        ax1, ax2, ax3 = self.hist_df.plot.hist(column="A", by="C", sharex=True)

        # share x
        assert self.get_x_axis(ax1).joined(ax1, ax2)
        assert self.get_x_axis(ax2).joined(ax1, ax2)
        assert self.get_x_axis(ax3).joined(ax1, ax3)
        assert self.get_x_axis(ax3).joined(ax2, ax3)

        # don't share y
        assert not self.get_y_axis(ax1).joined(ax1, ax2)
        assert not self.get_y_axis(ax2).joined(ax1, ax2)
        assert not self.get_y_axis(ax3).joined(ax1, ax3)
        assert not self.get_y_axis(ax3).joined(ax2, ax3)

    @pytest.mark.slow
    def test_axis_share_y_with_by(self):
        # GH 15079
        ax1, ax2, ax3 = self.hist_df.plot.hist(column="A", by="C", sharey=True)

        # share y
        assert self.get_y_axis(ax1).joined(ax1, ax2)
        assert self.get_y_axis(ax2).joined(ax1, ax2)
        assert self.get_y_axis(ax3).joined(ax1, ax3)
        assert self.get_y_axis(ax3).joined(ax2, ax3)

        # don't share x
        assert not self.get_x_axis(ax1).joined(ax1, ax2)
        assert not self.get_x_axis(ax2).joined(ax1, ax2)
        assert not self.get_x_axis(ax3).joined(ax1, ax3)
        assert not self.get_x_axis(ax3).joined(ax2, ax3)

    @pytest.mark.parametrize("figsize", [(12, 8), (20, 10)])
    def test_figure_shape_hist_with_by(self, figsize):
        # GH 15079
        axes = self.hist_df.plot.hist(column="A", by="C", figsize=figsize)
        self._check_axes_shape(axes, axes_num=3, figsize=figsize)


@td.skip_if_no_mpl
class TestBoxWithBy(TestPlotBase):
    def setup_method(self, method):
        TestPlotBase.setup_method(self, method)
        import matplotlib as mpl

        mpl.rcdefaults()
        self.box_df = _create_hist_box_with_by_df()

    @pytest.mark.parametrize(
        "by, column, titles, xticklabels",
        [
            ("C", "A", ["A"], [["a", "b", "c"]]),
            (
                ["C", "D"],
                "A",
                ["A"],
                [
                    [
                        "(a, a)",
                        "(a, b)",
                        "(a, c)",
                        "(b, a)",
                        "(b, b)",
                        "(b, c)",
                        "(c, a)",
                        "(c, b)",
                        "(c, c)",
                    ]
                ],
            ),
            ("C", ["A", "B"], ["A", "B"], [["a", "b", "c"]] * 2),
            (
                ["C", "D"],
                ["A", "B"],
                ["A", "B"],
                [
                    [
                        "(a, a)",
                        "(a, b)",
                        "(a, c)",
                        "(b, a)",
                        "(b, b)",
                        "(b, c)",
                        "(c, a)",
                        "(c, b)",
                        "(c, c)",
                    ]
                ]
                * 2,
            ),
            (["C"], None, ["A", "B"], [["a", "b", "c"]] * 2),
        ],
    )
    def test_box_plot_by_argument(self, by, column, titles, xticklabels):
        # GH 15079
        axes = _check_plot_works(self.box_df.plot.box, column=column, by=by)
        result_titles = [ax.get_title() for ax in axes]
        result_xticklabels = [
            [label.get_text() for label in ax.get_xticklabels()] for ax in axes
        ]

        assert result_xticklabels == xticklabels
        assert result_titles == titles

    @pytest.mark.parametrize(
        "by, column, titles, xticklabels",
        [
            (0, "A", ["A"], [["a", "b", "c"]]),
            (
                [0, "D"],
                "A",
                ["A"],
                [
                    [
                        "(a, a)",
                        "(a, b)",
                        "(a, c)",
                        "(b, a)",
                        "(b, b)",
                        "(b, c)",
                        "(c, a)",
                        "(c, b)",
                        "(c, c)",
                    ]
                ],
            ),
            (0, None, ["A", "B"], [["a", "b", "c"]] * 2),
        ],
    )
    def test_box_plot_by_0(self, by, column, titles, xticklabels):
        # GH 15079
        df = self.box_df.copy()
        df = df.rename(columns={"C": 0})

        axes = _check_plot_works(df.plot.box, column=column, by=by)
        result_titles = [ax.get_title() for ax in axes]
        result_xticklabels = [
            [label.get_text() for label in ax.get_xticklabels()] for ax in axes
        ]

        assert result_xticklabels == xticklabels
        assert result_titles == titles

    @pytest.mark.parametrize(
        "by, column",
        [
            ([], ["A"]),
            ((), "A"),
            ([], None),
            ((), ["A", "B"]),
        ],
    )
    def test_box_plot_with_none_empty_list_by(self, by, column):
        # GH 15079
        msg = "No group keys passed"
        with pytest.raises(ValueError, match=msg):
            _check_plot_works(self.box_df.plot.box, column=column, by=by)

    @pytest.mark.slow
    @pytest.mark.parametrize(
        "by, column, layout, axes_num",
        [
            (["C"], "A", (1, 1), 1),
            ("C", "A", (1, 1), 1),
            ("C", None, (2, 1), 2),
            ("C", ["A", "B"], (1, 2), 2),
            (["C", "D"], "A", (1, 1), 1),
            (["C", "D"], None, (1, 2), 2),
        ],
    )
    def test_box_plot_layout_with_by(self, by, column, layout, axes_num):
        # GH 15079
        axes = _check_plot_works(
            self.box_df.plot.box, column=column, by=by, layout=layout
        )
        self._check_axes_shape(axes, axes_num=axes_num, layout=layout)

    @pytest.mark.parametrize(
        "msg, by, layout",
        [
            ("larger than required size", ["C", "D"], (1, 1)),
            (re.escape("Layout must be a tuple of (rows, columns)"), "C", (1,)),
            ("At least one dimension of layout must be positive", "C", (-1, -1)),
        ],
    )
    def test_box_plot_invalid_layout_with_by_raises(self, msg, by, layout):
        # GH 15079, test if error is raised when invalid layout is given

        with pytest.raises(ValueError, match=msg):
            self.box_df.plot.box(column=["A", "B"], by=by, layout=layout)

    @pytest.mark.parametrize("figsize", [(12, 8), (20, 10)])
    def test_figure_shape_hist_with_by(self, figsize):
        # GH 15079
        axes = self.box_df.plot.box(column="A", by="C", figsize=figsize)
        self._check_axes_shape(axes, axes_num=1, figsize=figsize)
