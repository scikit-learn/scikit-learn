import string

import pytest
import numpy as np

from pandas.api.types import CategoricalDtype
from pandas import Categorical
from pandas.tests.extension import base


def make_data():
    return np.random.choice(list(string.ascii_letters), size=100)


@pytest.fixture
def dtype():
    return CategoricalDtype()


@pytest.fixture
def data():
    """Length-100 PeriodArray for semantics test."""
    return Categorical(make_data())


@pytest.fixture
def data_missing():
    """Length 2 array with [NA, Valid]"""
    return Categorical([np.nan, 'A'])


@pytest.fixture
def data_for_sorting():
    return Categorical(['A', 'B', 'C'], categories=['C', 'A', 'B'],
                       ordered=True)


@pytest.fixture
def data_missing_for_sorting():
    return Categorical(['A', None, 'B'], categories=['B', 'A'],
                       ordered=True)


@pytest.fixture
def na_value():
    return np.nan


@pytest.fixture
def data_for_grouping():
    return Categorical(['a', 'a', None, None, 'b', 'b', 'a', 'c'])


class TestDtype(base.BaseDtypeTests):
    pass


class TestInterface(base.BaseInterfaceTests):
    @pytest.mark.skip(reason="Memory usage doesn't match")
    def test_memory_usage(self):
        # Is this deliberate?
        pass


class TestConstructors(base.BaseConstructorsTests):
    pass


class TestReshaping(base.BaseReshapingTests):
    @pytest.mark.skip(reason="Unobserved categories preseved in concat.")
    def test_concat_columns(self, data, na_value):
        pass

    @pytest.mark.skip(reason="Unobserved categories preseved in concat.")
    def test_align(self, data, na_value):
        pass

    @pytest.mark.skip(reason="Unobserved categories preseved in concat.")
    def test_align_frame(self, data, na_value):
        pass

    @pytest.mark.skip(reason="Unobserved categories preseved in concat.")
    def test_merge(self, data, na_value):
        pass


class TestGetitem(base.BaseGetitemTests):
    skip_take = pytest.mark.skip(reason="GH-20664.")

    @pytest.mark.skip(reason="Backwards compatibility")
    def test_getitem_scalar(self):
        # CategoricalDtype.type isn't "correct" since it should
        # be a parent of the elements (object). But don't want
        # to break things by changing.
        pass

    @skip_take
    def test_take(self):
        # TODO remove this once Categorical.take is fixed
        pass

    @skip_take
    def test_take_negative(self):
        pass

    @skip_take
    def test_take_pandas_style_negative_raises(self):
        pass

    @skip_take
    def test_take_non_na_fill_value(self):
        pass

    @skip_take
    def test_take_out_of_bounds_raises(self):
        pass

    @pytest.mark.skip(reason="GH-20747. Unobserved categories.")
    def test_take_series(self):
        pass

    @skip_take
    def test_reindex_non_na_fill_value(self):
        pass

    @pytest.mark.xfail(reason="Categorical.take buggy")
    def test_take_empty(self):
        pass

    @pytest.mark.xfail(reason="test not written correctly for categorical")
    def test_reindex(self):
        pass


class TestSetitem(base.BaseSetitemTests):
    pass


class TestMissing(base.BaseMissingTests):

    @pytest.mark.skip(reason="Not implemented")
    def test_fillna_limit_pad(self):
        pass

    @pytest.mark.skip(reason="Not implemented")
    def test_fillna_limit_backfill(self):
        pass


class TestMethods(base.BaseMethodsTests):
    pass

    @pytest.mark.skip(reason="Unobserved categories included")
    def test_value_counts(self, all_data, dropna):
        pass


class TestCasting(base.BaseCastingTests):
    pass
