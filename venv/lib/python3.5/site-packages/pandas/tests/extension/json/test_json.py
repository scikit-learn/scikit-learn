import operator
import collections

import pytest

import pandas as pd
import pandas.util.testing as tm
from pandas.compat import PY2, PY36
from pandas.tests.extension import base

from .array import JSONArray, JSONDtype, make_data

pytestmark = pytest.mark.skipif(PY2, reason="Py2 doesn't have a UserDict")


@pytest.fixture
def dtype():
    return JSONDtype()


@pytest.fixture
def data():
    """Length-100 PeriodArray for semantics test."""
    data = make_data()

    # Why the while loop? NumPy is unable to construct an ndarray from
    # equal-length ndarrays. Many of our operations involve coercing the
    # EA to an ndarray of objects. To avoid random test failures, we ensure
    # that our data is coercable to an ndarray. Several tests deal with only
    # the first two elements, so that's what we'll check.

    while len(data[0]) == len(data[1]):
        data = make_data()

    return JSONArray(data)


@pytest.fixture
def data_missing():
    """Length 2 array with [NA, Valid]"""
    return JSONArray([{}, {'a': 10}])


@pytest.fixture
def data_for_sorting():
    return JSONArray([{'b': 1}, {'c': 4}, {'a': 2, 'c': 3}])


@pytest.fixture
def data_missing_for_sorting():
    return JSONArray([{'b': 1}, {}, {'a': 4}])


@pytest.fixture
def na_value(dtype):
    return dtype.na_value


@pytest.fixture
def na_cmp():
    return operator.eq


@pytest.fixture
def data_for_grouping():
    return JSONArray([
        {'b': 1}, {'b': 1},
        {}, {},
        {'a': 0, 'c': 2}, {'a': 0, 'c': 2},
        {'b': 1},
        {'c': 2},
    ])


class BaseJSON(object):
    # NumPy doesn't handle an array of equal-length UserDicts.
    # The default assert_series_equal eventually does a
    # Series.values, which raises. We work around it by
    # converting the UserDicts to dicts.
    def assert_series_equal(self, left, right, **kwargs):
        if left.dtype.name == 'json':
            assert left.dtype == right.dtype
            left = pd.Series(JSONArray(left.values.astype(object)),
                             index=left.index, name=left.name)
            right = pd.Series(JSONArray(right.values.astype(object)),
                              index=right.index, name=right.name)
        tm.assert_series_equal(left, right, **kwargs)

    def assert_frame_equal(self, left, right, *args, **kwargs):
        tm.assert_index_equal(
            left.columns, right.columns,
            exact=kwargs.get('check_column_type', 'equiv'),
            check_names=kwargs.get('check_names', True),
            check_exact=kwargs.get('check_exact', False),
            check_categorical=kwargs.get('check_categorical', True),
            obj='{obj}.columns'.format(obj=kwargs.get('obj', 'DataFrame')))

        jsons = (left.dtypes == 'json').index

        for col in jsons:
            self.assert_series_equal(left[col], right[col],
                                     *args, **kwargs)

        left = left.drop(columns=jsons)
        right = right.drop(columns=jsons)
        tm.assert_frame_equal(left, right, *args, **kwargs)


class TestDtype(BaseJSON, base.BaseDtypeTests):
    pass


class TestInterface(BaseJSON, base.BaseInterfaceTests):
    def test_custom_asserts(self):
        # This would always trigger the KeyError from trying to put
        # an array of equal-length UserDicts inside an ndarray.
        data = JSONArray([collections.UserDict({'a': 1}),
                          collections.UserDict({'b': 2}),
                          collections.UserDict({'c': 3})])
        a = pd.Series(data)
        self.assert_series_equal(a, a)
        self.assert_frame_equal(a.to_frame(), a.to_frame())

        b = pd.Series(data.take([0, 0, 1]))
        with pytest.raises(AssertionError):
            self.assert_series_equal(a, b)

        with pytest.raises(AssertionError):
            self.assert_frame_equal(a.to_frame(), b.to_frame())


class TestConstructors(BaseJSON, base.BaseConstructorsTests):
    pass


class TestReshaping(BaseJSON, base.BaseReshapingTests):
    pass


class TestGetitem(BaseJSON, base.BaseGetitemTests):
    pass


class TestMissing(BaseJSON, base.BaseMissingTests):
    @pytest.mark.xfail(reason="Setting a dict as a scalar")
    def test_fillna_series(self):
        """We treat dictionaries as a mapping in fillna, not a scalar."""

    @pytest.mark.xfail(reason="Setting a dict as a scalar")
    def test_fillna_frame(self):
        """We treat dictionaries as a mapping in fillna, not a scalar."""


unhashable = pytest.mark.skip(reason="Unhashable")
unstable = pytest.mark.skipif(not PY36,  # 3.6 or higher
                              reason="Dictionary order unstable")


class TestMethods(BaseJSON, base.BaseMethodsTests):
    @unhashable
    def test_value_counts(self, all_data, dropna):
        pass

    @unhashable
    def test_sort_values_frame(self):
        # TODO (EA.factorize): see if _values_for_factorize allows this.
        pass

    @unstable
    def test_argsort(self, data_for_sorting):
        super(TestMethods, self).test_argsort(data_for_sorting)

    @unstable
    def test_argsort_missing(self, data_missing_for_sorting):
        super(TestMethods, self).test_argsort_missing(
            data_missing_for_sorting)

    @unstable
    @pytest.mark.parametrize('ascending', [True, False])
    def test_sort_values(self, data_for_sorting, ascending):
        super(TestMethods, self).test_sort_values(
            data_for_sorting, ascending)

    @unstable
    @pytest.mark.parametrize('ascending', [True, False])
    def test_sort_values_missing(self, data_missing_for_sorting, ascending):
        super(TestMethods, self).test_sort_values_missing(
            data_missing_for_sorting, ascending)


class TestCasting(BaseJSON, base.BaseCastingTests):
    @pytest.mark.xfail
    def test_astype_str(self):
        """This currently fails in NumPy on np.array(self, dtype=str) with

        *** ValueError: setting an array element with a sequence
        """


# We intentionally don't run base.BaseSetitemTests because pandas'
# internals has trouble setting sequences of values into scalar positions.


class TestGroupby(BaseJSON, base.BaseGroupbyTests):

    @unhashable
    def test_groupby_extension_transform(self):
        """
        This currently fails in Series.name.setter, since the
        name must be hashable, but the value is a dictionary.
        I think this is what we want, i.e. `.name` should be the original
        values, and not the values for factorization.
        """

    @unhashable
    def test_groupby_extension_apply(self):
        """
        This fails in Index._do_unique_check with

        >   hash(val)
        E   TypeError: unhashable type: 'UserDict' with

        I suspect that once we support Index[ExtensionArray],
        we'll be able to dispatch unique.
        """

    @unstable
    @pytest.mark.parametrize('as_index', [True, False])
    def test_groupby_extension_agg(self, as_index, data_for_grouping):
        super(TestGroupby, self).test_groupby_extension_agg(
            as_index, data_for_grouping
        )
