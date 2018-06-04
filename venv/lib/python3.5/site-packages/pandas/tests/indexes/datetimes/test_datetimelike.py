""" generic tests from the Datetimelike class """

from pandas.util import testing as tm
from pandas import DatetimeIndex, date_range

from ..datetimelike import DatetimeLike


class TestDatetimeIndex(DatetimeLike):
    _holder = DatetimeIndex

    def setup_method(self, method):
        self.indices = dict(index=tm.makeDateIndex(10),
                            index_dec=date_range('20130110', periods=10,
                                                 freq='-1D'))
        self.setup_indices()

    def create_index(self):
        return date_range('20130101', periods=5)

    def test_shift(self):
        pass  # handled in test_ops

    def test_pickle_compat_construction(self):
        pass

    def test_intersection(self):
        pass  # handled in test_setops

    def test_union(self):
        pass  # handled in test_setops
