# coding: utf-8

import string

import pandas as pd
import pandas.util.testing as tm
import pandas.util._test_decorators as td
import pytest

from numpy.random import randn

import pandas.tools.plotting as plotting

from pandas.tests.plotting.common import TestPlotBase


"""
Test cases for plot functions imported from deprecated
pandas.tools.plotting
"""


@td.skip_if_no_mpl
class TestDeprecatedNameSpace(TestPlotBase):

    @pytest.mark.slow
    @td.skip_if_no_scipy
    def test_scatter_plot_legacy(self):
        df = pd.DataFrame(randn(100, 2))

        with tm.assert_produces_warning(FutureWarning):
            plotting.scatter_matrix(df)

        with tm.assert_produces_warning(FutureWarning):
            pd.scatter_matrix(df)

    @pytest.mark.slow
    def test_boxplot_deprecated(self):
        df = pd.DataFrame(randn(6, 4),
                          index=list(string.ascii_letters[:6]),
                          columns=['one', 'two', 'three', 'four'])
        df['indic'] = ['foo', 'bar'] * 3

        with tm.assert_produces_warning(FutureWarning):
            plotting.boxplot(df, column=['one', 'two'],
                             by='indic')

    @pytest.mark.slow
    def test_radviz_deprecated(self):
        df = self.iris
        with tm.assert_produces_warning(FutureWarning):
            plotting.radviz(frame=df, class_column='Name')

    @pytest.mark.slow
    def test_plot_params(self):

        with tm.assert_produces_warning(FutureWarning):
            pd.plot_params['xaxis.compat'] = True
