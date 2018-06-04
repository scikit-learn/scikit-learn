# coding: utf-8

""" Test cases for misc plot functions """

import pytest

from pandas import DataFrame
from pandas.compat import lmap
import pandas.util.testing as tm
import pandas.util._test_decorators as td

import numpy as np
from numpy import random
from numpy.random import randn

import pandas.plotting as plotting
from pandas.tests.plotting.common import TestPlotBase, _check_plot_works


@td.skip_if_mpl
def test_import_error_message():
    # GH-19810
    df = DataFrame({"A": [1, 2]})

    with tm.assert_raises_regex(ImportError, 'matplotlib is required'):
        df.plot()


@td.skip_if_no_mpl
class TestSeriesPlots(TestPlotBase):

    def setup_method(self, method):
        TestPlotBase.setup_method(self, method)
        import matplotlib as mpl
        mpl.rcdefaults()

        self.ts = tm.makeTimeSeries()
        self.ts.name = 'ts'

    @pytest.mark.slow
    def test_autocorrelation_plot(self):
        from pandas.plotting import autocorrelation_plot
        _check_plot_works(autocorrelation_plot, series=self.ts)
        _check_plot_works(autocorrelation_plot, series=self.ts.values)

        ax = autocorrelation_plot(self.ts, label='Test')
        self._check_legend_labels(ax, labels=['Test'])

    @pytest.mark.slow
    def test_lag_plot(self):
        from pandas.plotting import lag_plot
        _check_plot_works(lag_plot, series=self.ts)
        _check_plot_works(lag_plot, series=self.ts, lag=5)

    @pytest.mark.slow
    def test_bootstrap_plot(self):
        from pandas.plotting import bootstrap_plot
        _check_plot_works(bootstrap_plot, series=self.ts, size=10)


@td.skip_if_no_mpl
class TestDataFramePlots(TestPlotBase):

    @td.xfail_if_mpl_2_2
    @td.skip_if_no_scipy
    def test_scatter_matrix_axis(self):
        scatter_matrix = plotting.scatter_matrix

        with tm.RNGContext(42):
            df = DataFrame(randn(100, 3))

        # we are plotting multiples on a sub-plot
        with tm.assert_produces_warning(UserWarning):
            axes = _check_plot_works(scatter_matrix, filterwarnings='always',
                                     frame=df, range_padding=.1)
        axes0_labels = axes[0][0].yaxis.get_majorticklabels()

        # GH 5662
        if self.mpl_ge_2_0_0:
            expected = ['-2', '0', '2']
        else:
            expected = ['-2', '-1', '0', '1', '2']
        self._check_text_labels(axes0_labels, expected)
        self._check_ticks_props(
            axes, xlabelsize=8, xrot=90, ylabelsize=8, yrot=0)

        df[0] = ((df[0] - 2) / 3)

        # we are plotting multiples on a sub-plot
        with tm.assert_produces_warning(UserWarning):
            axes = _check_plot_works(scatter_matrix, filterwarnings='always',
                                     frame=df, range_padding=.1)
        axes0_labels = axes[0][0].yaxis.get_majorticklabels()
        if self.mpl_ge_2_0_0:
            expected = ['-1.0', '-0.5', '0.0']
        else:
            expected = ['-1.2', '-1.0', '-0.8', '-0.6', '-0.4', '-0.2', '0.0']
        self._check_text_labels(axes0_labels, expected)
        self._check_ticks_props(
            axes, xlabelsize=8, xrot=90, ylabelsize=8, yrot=0)

    @pytest.mark.slow
    def test_andrews_curves(self):
        from pandas.plotting import andrews_curves
        from matplotlib import cm

        df = self.iris

        _check_plot_works(andrews_curves, frame=df, class_column='Name')

        rgba = ('#556270', '#4ECDC4', '#C7F464')
        ax = _check_plot_works(andrews_curves, frame=df,
                               class_column='Name', color=rgba)
        self._check_colors(
            ax.get_lines()[:10], linecolors=rgba, mapping=df['Name'][:10])

        cnames = ['dodgerblue', 'aquamarine', 'seagreen']
        ax = _check_plot_works(andrews_curves, frame=df,
                               class_column='Name', color=cnames)
        self._check_colors(
            ax.get_lines()[:10], linecolors=cnames, mapping=df['Name'][:10])

        ax = _check_plot_works(andrews_curves, frame=df,
                               class_column='Name', colormap=cm.jet)
        cmaps = lmap(cm.jet, np.linspace(0, 1, df['Name'].nunique()))
        self._check_colors(
            ax.get_lines()[:10], linecolors=cmaps, mapping=df['Name'][:10])

        length = 10
        df = DataFrame({"A": random.rand(length),
                        "B": random.rand(length),
                        "C": random.rand(length),
                        "Name": ["A"] * length})

        _check_plot_works(andrews_curves, frame=df, class_column='Name')

        rgba = ('#556270', '#4ECDC4', '#C7F464')
        ax = _check_plot_works(andrews_curves, frame=df,
                               class_column='Name', color=rgba)
        self._check_colors(
            ax.get_lines()[:10], linecolors=rgba, mapping=df['Name'][:10])

        cnames = ['dodgerblue', 'aquamarine', 'seagreen']
        ax = _check_plot_works(andrews_curves, frame=df,
                               class_column='Name', color=cnames)
        self._check_colors(
            ax.get_lines()[:10], linecolors=cnames, mapping=df['Name'][:10])

        ax = _check_plot_works(andrews_curves, frame=df,
                               class_column='Name', colormap=cm.jet)
        cmaps = lmap(cm.jet, np.linspace(0, 1, df['Name'].nunique()))
        self._check_colors(
            ax.get_lines()[:10], linecolors=cmaps, mapping=df['Name'][:10])

        colors = ['b', 'g', 'r']
        df = DataFrame({"A": [1, 2, 3],
                        "B": [1, 2, 3],
                        "C": [1, 2, 3],
                        "Name": colors})
        ax = andrews_curves(df, 'Name', color=colors)
        handles, labels = ax.get_legend_handles_labels()
        self._check_colors(handles, linecolors=colors)

        with tm.assert_produces_warning(FutureWarning):
            andrews_curves(data=df, class_column='Name')

    @pytest.mark.slow
    def test_parallel_coordinates(self):
        from pandas.plotting import parallel_coordinates
        from matplotlib import cm

        df = self.iris

        ax = _check_plot_works(parallel_coordinates,
                               frame=df, class_column='Name')
        nlines = len(ax.get_lines())
        nxticks = len(ax.xaxis.get_ticklabels())

        rgba = ('#556270', '#4ECDC4', '#C7F464')
        ax = _check_plot_works(parallel_coordinates,
                               frame=df, class_column='Name', color=rgba)
        self._check_colors(
            ax.get_lines()[:10], linecolors=rgba, mapping=df['Name'][:10])

        cnames = ['dodgerblue', 'aquamarine', 'seagreen']
        ax = _check_plot_works(parallel_coordinates,
                               frame=df, class_column='Name', color=cnames)
        self._check_colors(
            ax.get_lines()[:10], linecolors=cnames, mapping=df['Name'][:10])

        ax = _check_plot_works(parallel_coordinates,
                               frame=df, class_column='Name', colormap=cm.jet)
        cmaps = lmap(cm.jet, np.linspace(0, 1, df['Name'].nunique()))
        self._check_colors(
            ax.get_lines()[:10], linecolors=cmaps, mapping=df['Name'][:10])

        ax = _check_plot_works(parallel_coordinates,
                               frame=df, class_column='Name', axvlines=False)
        assert len(ax.get_lines()) == (nlines - nxticks)

        colors = ['b', 'g', 'r']
        df = DataFrame({"A": [1, 2, 3],
                        "B": [1, 2, 3],
                        "C": [1, 2, 3],
                        "Name": colors})
        ax = parallel_coordinates(df, 'Name', color=colors)
        handles, labels = ax.get_legend_handles_labels()
        self._check_colors(handles, linecolors=colors)

        with tm.assert_produces_warning(FutureWarning):
            parallel_coordinates(data=df, class_column='Name')
        with tm.assert_produces_warning(FutureWarning):
            parallel_coordinates(df, 'Name', colors=colors)

    @pytest.mark.xfail(reason="unreliable test")
    def test_parallel_coordinates_with_sorted_labels(self):
        """ For #15908 """
        from pandas.plotting import parallel_coordinates

        df = DataFrame({"feat": [i for i in range(30)],
                        "class": [2 for _ in range(10)] +
                        [3 for _ in range(10)] +
                        [1 for _ in range(10)]})
        ax = parallel_coordinates(df, 'class', sort_labels=True)
        polylines, labels = ax.get_legend_handles_labels()
        color_label_tuples = \
            zip([polyline.get_color() for polyline in polylines], labels)
        ordered_color_label_tuples = sorted(color_label_tuples,
                                            key=lambda x: x[1])
        prev_next_tupels = zip([i for i in ordered_color_label_tuples[0:-1]],
                               [i for i in ordered_color_label_tuples[1:]])
        for prev, nxt in prev_next_tupels:
            # labels and colors are ordered strictly increasing
            assert prev[1] < nxt[1] and prev[0] < nxt[0]

    @pytest.mark.slow
    def test_radviz(self):
        from pandas.plotting import radviz
        from matplotlib import cm

        df = self.iris
        _check_plot_works(radviz, frame=df, class_column='Name')

        rgba = ('#556270', '#4ECDC4', '#C7F464')
        ax = _check_plot_works(
            radviz, frame=df, class_column='Name', color=rgba)
        # skip Circle drawn as ticks
        patches = [p for p in ax.patches[:20] if p.get_label() != '']
        self._check_colors(
            patches[:10], facecolors=rgba, mapping=df['Name'][:10])

        cnames = ['dodgerblue', 'aquamarine', 'seagreen']
        _check_plot_works(radviz, frame=df, class_column='Name', color=cnames)
        patches = [p for p in ax.patches[:20] if p.get_label() != '']
        self._check_colors(patches, facecolors=cnames, mapping=df['Name'][:10])

        _check_plot_works(radviz, frame=df,
                          class_column='Name', colormap=cm.jet)
        cmaps = lmap(cm.jet, np.linspace(0, 1, df['Name'].nunique()))
        patches = [p for p in ax.patches[:20] if p.get_label() != '']
        self._check_colors(patches, facecolors=cmaps, mapping=df['Name'][:10])

        colors = [[0., 0., 1., 1.],
                  [0., 0.5, 1., 1.],
                  [1., 0., 0., 1.]]
        df = DataFrame({"A": [1, 2, 3],
                        "B": [2, 1, 3],
                        "C": [3, 2, 1],
                        "Name": ['b', 'g', 'r']})
        ax = radviz(df, 'Name', color=colors)
        handles, labels = ax.get_legend_handles_labels()
        self._check_colors(handles, facecolors=colors)

    @pytest.mark.slow
    def test_subplot_titles(self):
        df = self.iris.drop('Name', axis=1).head()
        # Use the column names as the subplot titles
        title = list(df.columns)

        # Case len(title) == len(df)
        plot = df.plot(subplots=True, title=title)
        assert [p.get_title() for p in plot] == title

        # Case len(title) > len(df)
        pytest.raises(ValueError, df.plot, subplots=True,
                      title=title + ["kittens > puppies"])

        # Case len(title) < len(df)
        pytest.raises(ValueError, df.plot, subplots=True, title=title[:2])

        # Case subplots=False and title is of type list
        pytest.raises(ValueError, df.plot, subplots=False, title=title)

        # Case df with 3 numeric columns but layout of (2,2)
        plot = df.drop('SepalWidth', axis=1).plot(subplots=True, layout=(2, 2),
                                                  title=title[:-1])
        title_list = [ax.get_title() for sublist in plot for ax in sublist]
        assert title_list == title[:3] + ['']

    def test_get_standard_colors_random_seed(self):
        # GH17525
        df = DataFrame(np.zeros((10, 10)))

        # Make sure that the random seed isn't reset by _get_standard_colors
        plotting.parallel_coordinates(df, 0)
        rand1 = random.random()
        plotting.parallel_coordinates(df, 0)
        rand2 = random.random()
        assert rand1 != rand2

        # Make sure it produces the same colors every time it's called
        from pandas.plotting._style import _get_standard_colors
        color1 = _get_standard_colors(1, color_type='random')
        color2 = _get_standard_colors(1, color_type='random')
        assert color1 == color2
