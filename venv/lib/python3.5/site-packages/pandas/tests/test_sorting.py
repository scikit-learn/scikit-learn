import pytest
from itertools import product
from collections import defaultdict
import warnings
from datetime import datetime

import numpy as np
from numpy import nan
from pandas.core import common as com
from pandas import (DataFrame, MultiIndex, merge, concat, Series, compat,
                    _np_version_under1p10)
from pandas.util import testing as tm
from pandas.util.testing import assert_frame_equal, assert_series_equal
from pandas.core.sorting import (is_int64_overflow_possible,
                                 decons_group_index,
                                 get_group_index,
                                 nargsort,
                                 lexsort_indexer,
                                 safe_sort)


class TestSorting(object):

    @pytest.mark.slow
    def test_int64_overflow(self):

        B = np.concatenate((np.arange(1000), np.arange(1000), np.arange(500)))
        A = np.arange(2500)
        df = DataFrame({'A': A,
                        'B': B,
                        'C': A,
                        'D': B,
                        'E': A,
                        'F': B,
                        'G': A,
                        'H': B,
                        'values': np.random.randn(2500)})

        lg = df.groupby(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'])
        rg = df.groupby(['H', 'G', 'F', 'E', 'D', 'C', 'B', 'A'])

        left = lg.sum()['values']
        right = rg.sum()['values']

        exp_index, _ = left.index.sortlevel()
        tm.assert_index_equal(left.index, exp_index)

        exp_index, _ = right.index.sortlevel(0)
        tm.assert_index_equal(right.index, exp_index)

        tups = list(map(tuple, df[['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'
                                   ]].values))
        tups = com._asarray_tuplesafe(tups)

        expected = df.groupby(tups).sum()['values']

        for k, v in compat.iteritems(expected):
            assert left[k] == right[k[::-1]]
            assert left[k] == v
        assert len(left) == len(right)

    def test_int64_overflow_moar(self):

        # GH9096
        values = range(55109)
        data = DataFrame.from_dict(
            {'a': values, 'b': values, 'c': values, 'd': values})
        grouped = data.groupby(['a', 'b', 'c', 'd'])
        assert len(grouped) == len(values)

        arr = np.random.randint(-1 << 12, 1 << 12, (1 << 15, 5))
        i = np.random.choice(len(arr), len(arr) * 4)
        arr = np.vstack((arr, arr[i]))  # add sume duplicate rows

        i = np.random.permutation(len(arr))
        arr = arr[i]  # shuffle rows

        df = DataFrame(arr, columns=list('abcde'))
        df['jim'], df['joe'] = np.random.randn(2, len(df)) * 10
        gr = df.groupby(list('abcde'))

        # verify this is testing what it is supposed to test!
        assert is_int64_overflow_possible(gr.grouper.shape)

        # manually compute groupings
        jim, joe = defaultdict(list), defaultdict(list)
        for key, a, b in zip(map(tuple, arr), df['jim'], df['joe']):
            jim[key].append(a)
            joe[key].append(b)

        assert len(gr) == len(jim)
        mi = MultiIndex.from_tuples(jim.keys(), names=list('abcde'))

        def aggr(func):
            f = lambda a: np.fromiter(map(func, a), dtype='f8')
            arr = np.vstack((f(jim.values()), f(joe.values()))).T
            res = DataFrame(arr, columns=['jim', 'joe'], index=mi)
            return res.sort_index()

        assert_frame_equal(gr.mean(), aggr(np.mean))
        assert_frame_equal(gr.median(), aggr(np.median))

    def test_lexsort_indexer(self):
        keys = [[nan] * 5 + list(range(100)) + [nan] * 5]
        # orders=True, na_position='last'
        result = lexsort_indexer(keys, orders=True, na_position='last')
        exp = list(range(5, 105)) + list(range(5)) + list(range(105, 110))
        tm.assert_numpy_array_equal(result, np.array(exp, dtype=np.intp))

        # orders=True, na_position='first'
        result = lexsort_indexer(keys, orders=True, na_position='first')
        exp = list(range(5)) + list(range(105, 110)) + list(range(5, 105))
        tm.assert_numpy_array_equal(result, np.array(exp, dtype=np.intp))

        # orders=False, na_position='last'
        result = lexsort_indexer(keys, orders=False, na_position='last')
        exp = list(range(104, 4, -1)) + list(range(5)) + list(range(105, 110))
        tm.assert_numpy_array_equal(result, np.array(exp, dtype=np.intp))

        # orders=False, na_position='first'
        result = lexsort_indexer(keys, orders=False, na_position='first')
        exp = list(range(5)) + list(range(105, 110)) + list(range(104, 4, -1))
        tm.assert_numpy_array_equal(result, np.array(exp, dtype=np.intp))

    def test_nargsort(self):
        # np.argsort(items) places NaNs last
        items = [nan] * 5 + list(range(100)) + [nan] * 5
        # np.argsort(items2) may not place NaNs first
        items2 = np.array(items, dtype='O')

        try:
            # GH 2785; due to a regression in NumPy1.6.2
            np.argsort(np.array([[1, 2], [1, 3], [1, 2]], dtype='i'))
            np.argsort(items2, kind='mergesort')
        except TypeError:
            pytest.skip('requested sort not available for type')

        # mergesort is the most difficult to get right because we want it to be
        # stable.

        # According to numpy/core/tests/test_multiarray, """The number of
        # sorted items must be greater than ~50 to check the actual algorithm
        # because quick and merge sort fall over to insertion sort for small
        # arrays."""

        # mergesort, ascending=True, na_position='last'
        result = nargsort(items, kind='mergesort', ascending=True,
                          na_position='last')
        exp = list(range(5, 105)) + list(range(5)) + list(range(105, 110))
        tm.assert_numpy_array_equal(result, np.array(exp), check_dtype=False)

        # mergesort, ascending=True, na_position='first'
        result = nargsort(items, kind='mergesort', ascending=True,
                          na_position='first')
        exp = list(range(5)) + list(range(105, 110)) + list(range(5, 105))
        tm.assert_numpy_array_equal(result, np.array(exp), check_dtype=False)

        # mergesort, ascending=False, na_position='last'
        result = nargsort(items, kind='mergesort', ascending=False,
                          na_position='last')
        exp = list(range(104, 4, -1)) + list(range(5)) + list(range(105, 110))
        tm.assert_numpy_array_equal(result, np.array(exp), check_dtype=False)

        # mergesort, ascending=False, na_position='first'
        result = nargsort(items, kind='mergesort', ascending=False,
                          na_position='first')
        exp = list(range(5)) + list(range(105, 110)) + list(range(104, 4, -1))
        tm.assert_numpy_array_equal(result, np.array(exp), check_dtype=False)

        # mergesort, ascending=True, na_position='last'
        result = nargsort(items2, kind='mergesort', ascending=True,
                          na_position='last')
        exp = list(range(5, 105)) + list(range(5)) + list(range(105, 110))
        tm.assert_numpy_array_equal(result, np.array(exp), check_dtype=False)

        # mergesort, ascending=True, na_position='first'
        result = nargsort(items2, kind='mergesort', ascending=True,
                          na_position='first')
        exp = list(range(5)) + list(range(105, 110)) + list(range(5, 105))
        tm.assert_numpy_array_equal(result, np.array(exp), check_dtype=False)

        # mergesort, ascending=False, na_position='last'
        result = nargsort(items2, kind='mergesort', ascending=False,
                          na_position='last')
        exp = list(range(104, 4, -1)) + list(range(5)) + list(range(105, 110))
        tm.assert_numpy_array_equal(result, np.array(exp), check_dtype=False)

        # mergesort, ascending=False, na_position='first'
        result = nargsort(items2, kind='mergesort', ascending=False,
                          na_position='first')
        exp = list(range(5)) + list(range(105, 110)) + list(range(104, 4, -1))
        tm.assert_numpy_array_equal(result, np.array(exp), check_dtype=False)


class TestMerge(object):

    @pytest.mark.slow
    def test_int64_overflow_issues(self):

        # #2690, combinatorial explosion
        df1 = DataFrame(np.random.randn(1000, 7),
                        columns=list('ABCDEF') + ['G1'])
        df2 = DataFrame(np.random.randn(1000, 7),
                        columns=list('ABCDEF') + ['G2'])

        # it works!
        result = merge(df1, df2, how='outer')
        assert len(result) == 2000

        low, high, n = -1 << 10, 1 << 10, 1 << 20
        left = DataFrame(np.random.randint(low, high, (n, 7)),
                         columns=list('ABCDEFG'))
        left['left'] = left.sum(axis=1)

        # one-2-one match
        i = np.random.permutation(len(left))
        right = left.iloc[i].copy()
        right.columns = right.columns[:-1].tolist() + ['right']
        right.index = np.arange(len(right))
        right['right'] *= -1

        out = merge(left, right, how='outer')
        assert len(out) == len(left)
        assert_series_equal(out['left'], - out['right'], check_names=False)
        result = out.iloc[:, :-2].sum(axis=1)
        assert_series_equal(out['left'], result, check_names=False)
        assert result.name is None

        out.sort_values(out.columns.tolist(), inplace=True)
        out.index = np.arange(len(out))
        for how in ['left', 'right', 'outer', 'inner']:
            assert_frame_equal(out, merge(left, right, how=how, sort=True))

        # check that left merge w/ sort=False maintains left frame order
        out = merge(left, right, how='left', sort=False)
        assert_frame_equal(left, out[left.columns.tolist()])

        out = merge(right, left, how='left', sort=False)
        assert_frame_equal(right, out[right.columns.tolist()])

        # one-2-many/none match
        n = 1 << 11
        left = DataFrame(np.random.randint(low, high, (n, 7)).astype('int64'),
                         columns=list('ABCDEFG'))

        # confirm that this is checking what it is supposed to check
        shape = left.apply(Series.nunique).values
        assert is_int64_overflow_possible(shape)

        # add duplicates to left frame
        left = concat([left, left], ignore_index=True)

        right = DataFrame(np.random.randint(low, high, (n // 2, 7))
                          .astype('int64'),
                          columns=list('ABCDEFG'))

        # add duplicates & overlap with left to the right frame
        i = np.random.choice(len(left), n)
        right = concat([right, right, left.iloc[i]], ignore_index=True)

        left['left'] = np.random.randn(len(left))
        right['right'] = np.random.randn(len(right))

        # shuffle left & right frames
        i = np.random.permutation(len(left))
        left = left.iloc[i].copy()
        left.index = np.arange(len(left))

        i = np.random.permutation(len(right))
        right = right.iloc[i].copy()
        right.index = np.arange(len(right))

        # manually compute outer merge
        ldict, rdict = defaultdict(list), defaultdict(list)

        for idx, row in left.set_index(list('ABCDEFG')).iterrows():
            ldict[idx].append(row['left'])

        for idx, row in right.set_index(list('ABCDEFG')).iterrows():
            rdict[idx].append(row['right'])

        vals = []
        for k, lval in ldict.items():
            rval = rdict.get(k, [np.nan])
            for lv, rv in product(lval, rval):
                vals.append(k + tuple([lv, rv]))

        for k, rval in rdict.items():
            if k not in ldict:
                for rv in rval:
                    vals.append(k + tuple([np.nan, rv]))

        def align(df):
            df = df.sort_values(df.columns.tolist())
            df.index = np.arange(len(df))
            return df

        def verify_order(df):
            kcols = list('ABCDEFG')
            assert_frame_equal(df[kcols].copy(),
                               df[kcols].sort_values(kcols, kind='mergesort'))

        out = DataFrame(vals, columns=list('ABCDEFG') + ['left', 'right'])
        out = align(out)

        jmask = {'left': out['left'].notna(),
                 'right': out['right'].notna(),
                 'inner': out['left'].notna() & out['right'].notna(),
                 'outer': np.ones(len(out), dtype='bool')}

        for how in 'left', 'right', 'outer', 'inner':
            mask = jmask[how]
            frame = align(out[mask].copy())
            assert mask.all() ^ mask.any() or how == 'outer'

            for sort in [False, True]:
                res = merge(left, right, how=how, sort=sort)
                if sort:
                    verify_order(res)

                # as in GH9092 dtypes break with outer/right join
                assert_frame_equal(frame, align(res),
                                   check_dtype=how not in ('right', 'outer'))


def test_decons():

    def testit(label_list, shape):
        group_index = get_group_index(label_list, shape, sort=True, xnull=True)
        label_list2 = decons_group_index(group_index, shape)

        for a, b in zip(label_list, label_list2):
            tm.assert_numpy_array_equal(a, b)

    shape = (4, 5, 6)
    label_list = [np.tile([0, 1, 2, 3, 0, 1, 2, 3], 100).astype(np.int64),
                  np.tile([0, 2, 4, 3, 0, 1, 2, 3], 100).astype(np.int64),
                  np.tile([5, 1, 0, 2, 3, 0, 5, 4], 100).astype(np.int64)]
    testit(label_list, shape)

    shape = (10000, 10000)
    label_list = [np.tile(np.arange(10000, dtype=np.int64), 5),
                  np.tile(np.arange(10000, dtype=np.int64), 5)]
    testit(label_list, shape)


class TestSafeSort(object):

    def test_basic_sort(self):
        values = [3, 1, 2, 0, 4]
        result = safe_sort(values)
        expected = np.array([0, 1, 2, 3, 4])
        tm.assert_numpy_array_equal(result, expected)

        values = list("baaacb")
        result = safe_sort(values)
        expected = np.array(list("aaabbc"), dtype='object')
        tm.assert_numpy_array_equal(result, expected)

        values = []
        result = safe_sort(values)
        expected = np.array([])
        tm.assert_numpy_array_equal(result, expected)

    def test_labels(self):
        values = [3, 1, 2, 0, 4]
        expected = np.array([0, 1, 2, 3, 4])

        labels = [0, 1, 1, 2, 3, 0, -1, 4]
        result, result_labels = safe_sort(values, labels)
        expected_labels = np.array([3, 1, 1, 2, 0, 3, -1, 4], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)
        tm.assert_numpy_array_equal(result_labels, expected_labels)

        # na_sentinel
        labels = [0, 1, 1, 2, 3, 0, 99, 4]
        result, result_labels = safe_sort(values, labels,
                                          na_sentinel=99)
        expected_labels = np.array([3, 1, 1, 2, 0, 3, 99, 4], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)
        tm.assert_numpy_array_equal(result_labels, expected_labels)

        # out of bound indices
        labels = [0, 101, 102, 2, 3, 0, 99, 4]
        result, result_labels = safe_sort(values, labels)
        expected_labels = np.array([3, -1, -1, 2, 0, 3, -1, 4], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)
        tm.assert_numpy_array_equal(result_labels, expected_labels)

        labels = []
        result, result_labels = safe_sort(values, labels)
        expected_labels = np.array([], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)
        tm.assert_numpy_array_equal(result_labels, expected_labels)

    def test_mixed_integer(self):
        values = np.array(['b', 1, 0, 'a', 0, 'b'], dtype=object)
        result = safe_sort(values)
        expected = np.array([0, 0, 1, 'a', 'b', 'b'], dtype=object)
        tm.assert_numpy_array_equal(result, expected)

        values = np.array(['b', 1, 0, 'a'], dtype=object)
        labels = [0, 1, 2, 3, 0, -1, 1]
        result, result_labels = safe_sort(values, labels)
        expected = np.array([0, 1, 'a', 'b'], dtype=object)
        expected_labels = np.array([3, 1, 0, 2, 3, -1, 1], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)
        tm.assert_numpy_array_equal(result_labels, expected_labels)

    def test_mixed_integer_from_list(self):
        values = ['b', 1, 0, 'a', 0, 'b']
        result = safe_sort(values)
        expected = np.array([0, 0, 1, 'a', 'b', 'b'], dtype=object)
        tm.assert_numpy_array_equal(result, expected)

    def test_unsortable(self):
        # GH 13714
        arr = np.array([1, 2, datetime.now(), 0, 3], dtype=object)
        if compat.PY2 and not _np_version_under1p10:
            # RuntimeWarning: tp_compare didn't return -1 or -2 for exception
            with warnings.catch_warnings():
                pytest.raises(TypeError, safe_sort, arr)
        else:
            pytest.raises(TypeError, safe_sort, arr)

    def test_exceptions(self):
        with tm.assert_raises_regex(TypeError,
                                    "Only list-like objects are allowed"):
            safe_sort(values=1)

        with tm.assert_raises_regex(TypeError,
                                    "Only list-like objects or None"):
            safe_sort(values=[0, 1, 2], labels=1)

        with tm.assert_raises_regex(ValueError,
                                    "values should be unique"):
            safe_sort(values=[0, 1, 2, 1], labels=[0, 1])
