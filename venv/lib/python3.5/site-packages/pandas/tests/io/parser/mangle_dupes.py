# -*- coding: utf-8 -*-

"""
Tests that duplicate columns are handled appropriately when parsed by the
CSV engine. In general, the expected result is that they are either thoroughly
de-duplicated (if mangling requested) or ignored otherwise.
"""

from pandas.compat import StringIO
from pandas import DataFrame

import pandas.util.testing as tm


class DupeColumnTests(object):
    def test_basic(self):
        # TODO: add test for condition "mangle_dupe_cols=False"
        # once it is actually supported (gh-12935)
        data = "a,a,b,b,b\n1,2,3,4,5"

        for method in ("read_csv", "read_table"):
            # Check default behavior.
            expected = ["a", "a.1", "b", "b.1", "b.2"]
            df = getattr(self, method)(StringIO(data), sep=",")
            assert list(df.columns) == expected

            df = getattr(self, method)(StringIO(data), sep=",",
                                       mangle_dupe_cols=True)
            assert list(df.columns) == expected

    def test_basic_names(self):
        # See gh-7160
        data = "a,b,a\n0,1,2\n3,4,5"
        expected = DataFrame([[0, 1, 2], [3, 4, 5]],
                             columns=["a", "b", "a.1"])

        df = self.read_csv(StringIO(data))
        tm.assert_frame_equal(df, expected)

        with tm.assert_produces_warning(UserWarning, check_stacklevel=False):
            data = "0,1,2\n3,4,5"
            df = self.read_csv(StringIO(data),
                               names=["a", "b", "a"])
            tm.assert_frame_equal(df, expected)

    def test_thorough_mangle_columns(self):
        # see gh-17060
        data = "a,a,a.1\n1,2,3"
        df = self.read_csv(StringIO(data), sep=",", mangle_dupe_cols=True)
        assert list(df.columns) == ["a", "a.1", "a.1.1"]

        data = "a,a,a.1,a.1.1,a.1.1.1,a.1.1.1.1\n1,2,3,4,5,6"
        df = self.read_csv(StringIO(data), sep=",", mangle_dupe_cols=True)
        assert list(df.columns) == ["a", "a.1", "a.1.1", "a.1.1.1",
                                    "a.1.1.1.1", "a.1.1.1.1.1"]

        data = "a,a,a.3,a.1,a.2,a,a\n1,2,3,4,5,6,7"
        df = self.read_csv(StringIO(data), sep=",", mangle_dupe_cols=True)
        assert list(df.columns) == ["a", "a.1", "a.3", "a.1.1",
                                    "a.2", "a.2.1", "a.3.1"]

    def test_thorough_mangle_names(self):
        # see gh-17095
        data = "a,b,b\n1,2,3"
        names = ["a.1", "a.1", "a.1.1"]

        with tm.assert_produces_warning(UserWarning, check_stacklevel=False):
            df = self.read_csv(StringIO(data), sep=",", names=names,
                               mangle_dupe_cols=True)
            assert list(df.columns) == ["a.1", "a.1.1", "a.1.1.1"]

        data = "a,b,c,d,e,f\n1,2,3,4,5,6"
        names = ["a", "a", "a.1", "a.1.1", "a.1.1.1", "a.1.1.1.1"]

        with tm.assert_produces_warning(UserWarning, check_stacklevel=False):
            df = self.read_csv(StringIO(data), sep=",", names=names,
                               mangle_dupe_cols=True)
            assert list(df.columns) == ["a", "a.1", "a.1.1", "a.1.1.1",
                                        "a.1.1.1.1", "a.1.1.1.1.1"]

        data = "a,b,c,d,e,f,g\n1,2,3,4,5,6,7"
        names = ["a", "a", "a.3", "a.1", "a.2", "a", "a"]

        with tm.assert_produces_warning(UserWarning, check_stacklevel=False):
            df = self.read_csv(StringIO(data), sep=",", names=names,
                               mangle_dupe_cols=True)
            assert list(df.columns) == ["a", "a.1", "a.3", "a.1.1",
                                        "a.2", "a.2.1", "a.3.1"]
