# -*- coding: utf-8 -*-

import pytest

from warnings import catch_warnings
import numpy as np
from pandas import (Series, DataFrame, Index, Float64Index, Int64Index,
                    RangeIndex)
from pandas.util.testing import assert_series_equal, assert_almost_equal
import pandas.util.testing as tm


class TestFloatIndexers(object):

    def check(self, result, original, indexer, getitem):
        """
        comparator for results
        we need to take care if we are indexing on a
        Series or a frame
        """
        if isinstance(original, Series):
            expected = original.iloc[indexer]
        else:
            if getitem:
                expected = original.iloc[:, indexer]
            else:
                expected = original.iloc[indexer]

        assert_almost_equal(result, expected)

    def test_scalar_error(self):

        # GH 4892
        # float_indexers should raise exceptions
        # on appropriate Index types & accessors
        # this duplicates the code below
        # but is spefically testing for the error
        # message

        for index in [tm.makeStringIndex, tm.makeUnicodeIndex,
                      tm.makeCategoricalIndex,
                      tm.makeDateIndex, tm.makeTimedeltaIndex,
                      tm.makePeriodIndex, tm.makeIntIndex,
                      tm.makeRangeIndex]:

            i = index(5)

            s = Series(np.arange(len(i)), index=i)

            def f():
                s.iloc[3.0]
            tm.assert_raises_regex(TypeError,
                                   'cannot do positional indexing',
                                   f)

            def f():
                s.iloc[3.0] = 0
            pytest.raises(TypeError, f)

    def test_scalar_non_numeric(self):

        # GH 4892
        # float_indexers should raise exceptions
        # on appropriate Index types & accessors

        for index in [tm.makeStringIndex, tm.makeUnicodeIndex,
                      tm.makeCategoricalIndex,
                      tm.makeDateIndex, tm.makeTimedeltaIndex,
                      tm.makePeriodIndex]:

            i = index(5)

            for s in [Series(
                    np.arange(len(i)), index=i), DataFrame(
                        np.random.randn(
                            len(i), len(i)), index=i, columns=i)]:

                # getting
                for idxr, getitem in [(lambda x: x.ix, False),
                                      (lambda x: x.iloc, False),
                                      (lambda x: x, True)]:

                    def f():
                        with catch_warnings(record=True):
                            idxr(s)[3.0]

                    # gettitem on a DataFrame is a KeyError as it is indexing
                    # via labels on the columns
                    if getitem and isinstance(s, DataFrame):
                        error = KeyError
                    else:
                        error = TypeError
                    pytest.raises(error, f)

                # label based can be a TypeError or KeyError
                def f():
                    s.loc[3.0]

                if s.index.inferred_type in ['string', 'unicode', 'mixed']:
                    error = KeyError
                else:
                    error = TypeError
                pytest.raises(error, f)

                # contains
                assert 3.0 not in s

                # setting with a float fails with iloc
                def f():
                    s.iloc[3.0] = 0
                pytest.raises(TypeError, f)

                # setting with an indexer
                if s.index.inferred_type in ['categorical']:
                    # Value or Type Error
                    pass
                elif s.index.inferred_type in ['datetime64', 'timedelta64',
                                               'period']:

                    # these should prob work
                    # and are inconsisten between series/dataframe ATM
                    # for idxr in [lambda x: x.ix,
                    #             lambda x: x]:
                    #    s2 = s.copy()
                    #    def f():
                    #        idxr(s2)[3.0] = 0
                    #    pytest.raises(TypeError, f)
                    pass

                else:

                    s2 = s.copy()
                    s2.loc[3.0] = 10
                    assert s2.index.is_object()

                    for idxr in [lambda x: x.ix,
                                 lambda x: x]:
                        s2 = s.copy()
                        with catch_warnings(record=True):
                            idxr(s2)[3.0] = 0
                        assert s2.index.is_object()

            # fallsback to position selection, series only
            s = Series(np.arange(len(i)), index=i)
            s[3]
            pytest.raises(TypeError, lambda: s[3.0])

    def test_scalar_with_mixed(self):

        s2 = Series([1, 2, 3], index=['a', 'b', 'c'])
        s3 = Series([1, 2, 3], index=['a', 'b', 1.5])

        # lookup in a pure string index
        # with an invalid indexer
        for idxr in [lambda x: x.ix,
                     lambda x: x,
                     lambda x: x.iloc]:

            def f():
                with catch_warnings(record=True):
                    idxr(s2)[1.0]

            pytest.raises(TypeError, f)

        pytest.raises(KeyError, lambda: s2.loc[1.0])

        result = s2.loc['b']
        expected = 2
        assert result == expected

        # mixed index so we have label
        # indexing
        for idxr in [lambda x: x]:

            def f():
                idxr(s3)[1.0]

            pytest.raises(TypeError, f)

            result = idxr(s3)[1]
            expected = 2
            assert result == expected

        # mixed index so we have label
        # indexing
        for idxr in [lambda x: x.ix]:
            with catch_warnings(record=True):

                def f():
                    idxr(s3)[1.0]

                pytest.raises(TypeError, f)

                result = idxr(s3)[1]
                expected = 2
                assert result == expected

        pytest.raises(TypeError, lambda: s3.iloc[1.0])
        pytest.raises(KeyError, lambda: s3.loc[1.0])

        result = s3.loc[1.5]
        expected = 3
        assert result == expected

    def test_scalar_integer(self):

        # test how scalar float indexers work on int indexes

        # integer index
        for i in [Int64Index(range(5)), RangeIndex(5)]:

            for s in [Series(np.arange(len(i))),
                      DataFrame(np.random.randn(len(i), len(i)),
                                index=i, columns=i)]:

                # coerce to equal int
                for idxr, getitem in [(lambda x: x.ix, False),
                                      (lambda x: x.loc, False),
                                      (lambda x: x, True)]:

                    with catch_warnings(record=True):
                        result = idxr(s)[3.0]
                    self.check(result, s, 3, getitem)

                # coerce to equal int
                for idxr, getitem in [(lambda x: x.ix, False),
                                      (lambda x: x.loc, False),
                                      (lambda x: x, True)]:

                    if isinstance(s, Series):
                        def compare(x, y):
                            assert x == y
                        expected = 100
                    else:
                        compare = tm.assert_series_equal
                        if getitem:
                            expected = Series(100,
                                              index=range(len(s)), name=3)
                        else:
                            expected = Series(100.,
                                              index=range(len(s)), name=3)

                    s2 = s.copy()
                    with catch_warnings(record=True):
                        idxr(s2)[3.0] = 100

                        result = idxr(s2)[3.0]
                        compare(result, expected)

                        result = idxr(s2)[3]
                        compare(result, expected)

                # contains
                # coerce to equal int
                assert 3.0 in s

    def test_scalar_float(self):

        # scalar float indexers work on a float index
        index = Index(np.arange(5.))
        for s in [Series(np.arange(len(index)), index=index),
                  DataFrame(np.random.randn(len(index), len(index)),
                            index=index, columns=index)]:

            # assert all operations except for iloc are ok
            indexer = index[3]
            for idxr, getitem in [(lambda x: x.ix, False),
                                  (lambda x: x.loc, False),
                                  (lambda x: x, True)]:

                # getting
                with catch_warnings(record=True):
                    result = idxr(s)[indexer]
                self.check(result, s, 3, getitem)

                # setting
                s2 = s.copy()

                def f():
                    with catch_warnings(record=True):
                        idxr(s2)[indexer] = expected
                with catch_warnings(record=True):
                    result = idxr(s2)[indexer]
                self.check(result, s, 3, getitem)

                # random integer is a KeyError
                with catch_warnings(record=True):
                    pytest.raises(KeyError, lambda: idxr(s)[3.5])

            # contains
            assert 3.0 in s

            # iloc succeeds with an integer
            expected = s.iloc[3]
            s2 = s.copy()

            s2.iloc[3] = expected
            result = s2.iloc[3]
            self.check(result, s, 3, False)

            # iloc raises with a float
            pytest.raises(TypeError, lambda: s.iloc[3.0])

            def g():
                s2.iloc[3.0] = 0
            pytest.raises(TypeError, g)

    def test_slice_non_numeric(self):

        # GH 4892
        # float_indexers should raise exceptions
        # on appropriate Index types & accessors

        for index in [tm.makeStringIndex, tm.makeUnicodeIndex,
                      tm.makeDateIndex, tm.makeTimedeltaIndex,
                      tm.makePeriodIndex]:

            index = index(5)
            for s in [Series(range(5), index=index),
                      DataFrame(np.random.randn(5, 2), index=index)]:

                # getitem
                for l in [slice(3.0, 4),
                          slice(3, 4.0),
                          slice(3.0, 4.0)]:

                    def f():
                        s.iloc[l]
                    pytest.raises(TypeError, f)

                    for idxr in [lambda x: x.ix,
                                 lambda x: x.loc,
                                 lambda x: x.iloc,
                                 lambda x: x]:

                        def f():
                            with catch_warnings(record=True):
                                idxr(s)[l]
                        pytest.raises(TypeError, f)

                # setitem
                for l in [slice(3.0, 4),
                          slice(3, 4.0),
                          slice(3.0, 4.0)]:

                    def f():
                        s.iloc[l] = 0
                    pytest.raises(TypeError, f)

                    for idxr in [lambda x: x.ix,
                                 lambda x: x.loc,
                                 lambda x: x.iloc,
                                 lambda x: x]:
                        def f():
                            with catch_warnings(record=True):
                                idxr(s)[l] = 0
                        pytest.raises(TypeError, f)

    def test_slice_integer(self):

        # same as above, but for Integer based indexes
        # these coerce to a like integer
        # oob indicates if we are out of bounds
        # of positional indexing
        for index, oob in [(Int64Index(range(5)), False),
                           (RangeIndex(5), False),
                           (Int64Index(range(5)) + 10, True)]:

            # s is an in-range index
            s = Series(range(5), index=index)

            # getitem
            for l in [slice(3.0, 4),
                      slice(3, 4.0),
                      slice(3.0, 4.0)]:

                for idxr in [lambda x: x.loc,
                             lambda x: x.ix]:

                    with catch_warnings(record=True):
                        result = idxr(s)[l]

                    # these are all label indexing
                    # except getitem which is positional
                    # empty
                    if oob:
                        indexer = slice(0, 0)
                    else:
                        indexer = slice(3, 5)
                    self.check(result, s, indexer, False)

                # positional indexing
                def f():
                    s[l]

                pytest.raises(TypeError, f)

            # getitem out-of-bounds
            for l in [slice(-6, 6),
                      slice(-6.0, 6.0)]:

                for idxr in [lambda x: x.loc,
                             lambda x: x.ix]:
                    with catch_warnings(record=True):
                        result = idxr(s)[l]

                    # these are all label indexing
                    # except getitem which is positional
                    # empty
                    if oob:
                        indexer = slice(0, 0)
                    else:
                        indexer = slice(-6, 6)
                    self.check(result, s, indexer, False)

            # positional indexing
            def f():
                s[slice(-6.0, 6.0)]

            pytest.raises(TypeError, f)

            # getitem odd floats
            for l, res1 in [(slice(2.5, 4), slice(3, 5)),
                            (slice(2, 3.5), slice(2, 4)),
                            (slice(2.5, 3.5), slice(3, 4))]:

                for idxr in [lambda x: x.loc,
                             lambda x: x.ix]:

                    with catch_warnings(record=True):
                        result = idxr(s)[l]
                    if oob:
                        res = slice(0, 0)
                    else:
                        res = res1

                    self.check(result, s, res, False)

                # positional indexing
                def f():
                    s[l]

                pytest.raises(TypeError, f)

            # setitem
            for l in [slice(3.0, 4),
                      slice(3, 4.0),
                      slice(3.0, 4.0)]:

                for idxr in [lambda x: x.loc,
                             lambda x: x.ix]:
                    sc = s.copy()
                    with catch_warnings(record=True):
                        idxr(sc)[l] = 0
                        result = idxr(sc)[l].values.ravel()
                    assert (result == 0).all()

                # positional indexing
                def f():
                    s[l] = 0

                pytest.raises(TypeError, f)

    def test_integer_positional_indexing(self):
        """ make sure that we are raising on positional indexing
        w.r.t. an integer index """

        s = Series(range(2, 6), index=range(2, 6))

        result = s[2:4]
        expected = s.iloc[2:4]
        assert_series_equal(result, expected)

        for idxr in [lambda x: x,
                     lambda x: x.iloc]:

            for l in [slice(2, 4.0),
                      slice(2.0, 4),
                      slice(2.0, 4.0)]:

                def f():
                    idxr(s)[l]

                pytest.raises(TypeError, f)

    def test_slice_integer_frame_getitem(self):

        # similar to above, but on the getitem dim (of a DataFrame)
        for index in [Int64Index(range(5)), RangeIndex(5)]:

            s = DataFrame(np.random.randn(5, 2), index=index)

            def f(idxr):

                # getitem
                for l in [slice(0.0, 1),
                          slice(0, 1.0),
                          slice(0.0, 1.0)]:

                    result = idxr(s)[l]
                    indexer = slice(0, 2)
                    self.check(result, s, indexer, False)

                    # positional indexing
                    def f():
                        s[l]

                    pytest.raises(TypeError, f)

                # getitem out-of-bounds
                for l in [slice(-10, 10),
                          slice(-10.0, 10.0)]:

                    result = idxr(s)[l]
                    self.check(result, s, slice(-10, 10), True)

                # positional indexing
                def f():
                    s[slice(-10.0, 10.0)]

                pytest.raises(TypeError, f)

                # getitem odd floats
                for l, res in [(slice(0.5, 1), slice(1, 2)),
                               (slice(0, 0.5), slice(0, 1)),
                               (slice(0.5, 1.5), slice(1, 2))]:

                    result = idxr(s)[l]
                    self.check(result, s, res, False)

                    # positional indexing
                    def f():
                        s[l]

                    pytest.raises(TypeError, f)

                # setitem
                for l in [slice(3.0, 4),
                          slice(3, 4.0),
                          slice(3.0, 4.0)]:

                    sc = s.copy()
                    idxr(sc)[l] = 0
                    result = idxr(sc)[l].values.ravel()
                    assert (result == 0).all()

                    # positional indexing
                    def f():
                        s[l] = 0

                    pytest.raises(TypeError, f)

            f(lambda x: x.loc)
            with catch_warnings(record=True):
                f(lambda x: x.ix)

    def test_slice_float(self):

        # same as above, but for floats
        index = Index(np.arange(5.)) + 0.1
        for s in [Series(range(5), index=index),
                  DataFrame(np.random.randn(5, 2), index=index)]:

            for l in [slice(3.0, 4),
                      slice(3, 4.0),
                      slice(3.0, 4.0)]:

                expected = s.iloc[3:4]
                for idxr in [lambda x: x.ix,
                             lambda x: x.loc,
                             lambda x: x]:

                    # getitem
                    with catch_warnings(record=True):
                        result = idxr(s)[l]
                    if isinstance(s, Series):
                        tm.assert_series_equal(result, expected)
                    else:
                        tm.assert_frame_equal(result, expected)
                    # setitem
                    s2 = s.copy()
                    with catch_warnings(record=True):
                        idxr(s2)[l] = 0
                        result = idxr(s2)[l].values.ravel()
                    assert (result == 0).all()

    def test_floating_index_doc_example(self):

        index = Index([1.5, 2, 3, 4.5, 5])
        s = Series(range(5), index=index)
        assert s[3] == 2
        assert s.loc[3] == 2
        assert s.loc[3] == 2
        assert s.iloc[3] == 3

    def test_floating_misc(self):

        # related 236
        # scalar/slicing of a float index
        s = Series(np.arange(5), index=np.arange(5) * 2.5, dtype=np.int64)

        # label based slicing
        result1 = s[1.0:3.0]
        result2 = s.loc[1.0:3.0]
        result3 = s.loc[1.0:3.0]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)

        # exact indexing when found
        result1 = s[5.0]
        result2 = s.loc[5.0]
        result3 = s.loc[5.0]
        assert result1 == result2
        assert result1 == result3

        result1 = s[5]
        result2 = s.loc[5]
        result3 = s.loc[5]
        assert result1 == result2
        assert result1 == result3

        assert s[5.0] == s[5]

        # value not found (and no fallbacking at all)

        # scalar integers
        pytest.raises(KeyError, lambda: s.loc[4])
        pytest.raises(KeyError, lambda: s.loc[4])
        pytest.raises(KeyError, lambda: s[4])

        # fancy floats/integers create the correct entry (as nan)
        # fancy tests
        expected = Series([2, 0], index=Float64Index([5.0, 0.0]))
        for fancy_idx in [[5.0, 0.0], np.array([5.0, 0.0])]:  # float
            assert_series_equal(s[fancy_idx], expected)
            assert_series_equal(s.loc[fancy_idx], expected)
            assert_series_equal(s.loc[fancy_idx], expected)

        expected = Series([2, 0], index=Index([5, 0], dtype='int64'))
        for fancy_idx in [[5, 0], np.array([5, 0])]:  # int
            assert_series_equal(s[fancy_idx], expected)
            assert_series_equal(s.loc[fancy_idx], expected)
            assert_series_equal(s.loc[fancy_idx], expected)

        # all should return the same as we are slicing 'the same'
        result1 = s.loc[2:5]
        result2 = s.loc[2.0:5.0]
        result3 = s.loc[2.0:5]
        result4 = s.loc[2.1:5]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)
        assert_series_equal(result1, result4)

        # previously this did fallback indexing
        result1 = s[2:5]
        result2 = s[2.0:5.0]
        result3 = s[2.0:5]
        result4 = s[2.1:5]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)
        assert_series_equal(result1, result4)

        result1 = s.loc[2:5]
        result2 = s.loc[2.0:5.0]
        result3 = s.loc[2.0:5]
        result4 = s.loc[2.1:5]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)
        assert_series_equal(result1, result4)

        # combined test
        result1 = s.loc[2:5]
        result2 = s.loc[2:5]
        result3 = s[2:5]

        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)

        # list selection
        result1 = s[[0.0, 5, 10]]
        result2 = s.loc[[0.0, 5, 10]]
        result3 = s.loc[[0.0, 5, 10]]
        result4 = s.iloc[[0, 2, 4]]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)
        assert_series_equal(result1, result4)

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            result1 = s[[1.6, 5, 10]]
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            result2 = s.loc[[1.6, 5, 10]]
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            result3 = s.loc[[1.6, 5, 10]]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)
        assert_series_equal(result1, Series(
            [np.nan, 2, 4], index=[1.6, 5, 10]))

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            result1 = s[[0, 1, 2]]
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            result2 = s.loc[[0, 1, 2]]
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            result3 = s.loc[[0, 1, 2]]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)
        assert_series_equal(result1, Series(
            [0.0, np.nan, np.nan], index=[0, 1, 2]))

        result1 = s.loc[[2.5, 5]]
        result2 = s.loc[[2.5, 5]]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, Series([1, 2], index=[2.5, 5.0]))

        result1 = s[[2.5]]
        result2 = s.loc[[2.5]]
        result3 = s.loc[[2.5]]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)
        assert_series_equal(result1, Series([1], index=[2.5]))

    def test_floating_tuples(self):
        # see gh-13509
        s = Series([(1, 1), (2, 2), (3, 3)], index=[0.0, 0.1, 0.2], name='foo')

        result = s[0.0]
        assert result == (1, 1)

        expected = Series([(1, 1), (2, 2)], index=[0.0, 0.0], name='foo')
        s = Series([(1, 1), (2, 2), (3, 3)], index=[0.0, 0.0, 0.2], name='foo')

        result = s[0.0]
        tm.assert_series_equal(result, expected)

    def test_float64index_slicing_bug(self):
        # GH 5557, related to slicing a float index
        ser = {256: 2321.0,
               1: 78.0,
               2: 2716.0,
               3: 0.0,
               4: 369.0,
               5: 0.0,
               6: 269.0,
               7: 0.0,
               8: 0.0,
               9: 0.0,
               10: 3536.0,
               11: 0.0,
               12: 24.0,
               13: 0.0,
               14: 931.0,
               15: 0.0,
               16: 101.0,
               17: 78.0,
               18: 9643.0,
               19: 0.0,
               20: 0.0,
               21: 0.0,
               22: 63761.0,
               23: 0.0,
               24: 446.0,
               25: 0.0,
               26: 34773.0,
               27: 0.0,
               28: 729.0,
               29: 78.0,
               30: 0.0,
               31: 0.0,
               32: 3374.0,
               33: 0.0,
               34: 1391.0,
               35: 0.0,
               36: 361.0,
               37: 0.0,
               38: 61808.0,
               39: 0.0,
               40: 0.0,
               41: 0.0,
               42: 6677.0,
               43: 0.0,
               44: 802.0,
               45: 0.0,
               46: 2691.0,
               47: 0.0,
               48: 3582.0,
               49: 0.0,
               50: 734.0,
               51: 0.0,
               52: 627.0,
               53: 70.0,
               54: 2584.0,
               55: 0.0,
               56: 324.0,
               57: 0.0,
               58: 605.0,
               59: 0.0,
               60: 0.0,
               61: 0.0,
               62: 3989.0,
               63: 10.0,
               64: 42.0,
               65: 0.0,
               66: 904.0,
               67: 0.0,
               68: 88.0,
               69: 70.0,
               70: 8172.0,
               71: 0.0,
               72: 0.0,
               73: 0.0,
               74: 64902.0,
               75: 0.0,
               76: 347.0,
               77: 0.0,
               78: 36605.0,
               79: 0.0,
               80: 379.0,
               81: 70.0,
               82: 0.0,
               83: 0.0,
               84: 3001.0,
               85: 0.0,
               86: 1630.0,
               87: 7.0,
               88: 364.0,
               89: 0.0,
               90: 67404.0,
               91: 9.0,
               92: 0.0,
               93: 0.0,
               94: 7685.0,
               95: 0.0,
               96: 1017.0,
               97: 0.0,
               98: 2831.0,
               99: 0.0,
               100: 2963.0,
               101: 0.0,
               102: 854.0,
               103: 0.0,
               104: 0.0,
               105: 0.0,
               106: 0.0,
               107: 0.0,
               108: 0.0,
               109: 0.0,
               110: 0.0,
               111: 0.0,
               112: 0.0,
               113: 0.0,
               114: 0.0,
               115: 0.0,
               116: 0.0,
               117: 0.0,
               118: 0.0,
               119: 0.0,
               120: 0.0,
               121: 0.0,
               122: 0.0,
               123: 0.0,
               124: 0.0,
               125: 0.0,
               126: 67744.0,
               127: 22.0,
               128: 264.0,
               129: 0.0,
               260: 197.0,
               268: 0.0,
               265: 0.0,
               269: 0.0,
               261: 0.0,
               266: 1198.0,
               267: 0.0,
               262: 2629.0,
               258: 775.0,
               257: 0.0,
               263: 0.0,
               259: 0.0,
               264: 163.0,
               250: 10326.0,
               251: 0.0,
               252: 1228.0,
               253: 0.0,
               254: 2769.0,
               255: 0.0}

        # smoke test for the repr
        s = Series(ser)
        result = s.value_counts()
        str(result)
