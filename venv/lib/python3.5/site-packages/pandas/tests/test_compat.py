# -*- coding: utf-8 -*-
"""
Testing that functions from compat work as expected
"""

import pytest
from pandas.compat import (range, zip, map, filter, lrange, lzip, lmap,
                           lfilter, builtins, iterkeys, itervalues, iteritems,
                           next, get_range_parameters, PY2)


class TestBuiltinIterators(object):

    @classmethod
    def check_result(cls, actual, expected, lengths):
        for (iter_res, list_res), exp, length in zip(actual, expected,
                                                     lengths):
            assert not isinstance(iter_res, list)
            assert isinstance(list_res, list)

            iter_res = list(iter_res)

            assert len(list_res) == length
            assert len(iter_res) == length
            assert iter_res == exp
            assert list_res == exp

    def test_range(self):
        actual1 = range(10)
        actual2 = lrange(10)
        actual = [actual1, actual2],
        expected = list(builtins.range(10)),
        lengths = 10,

        actual1 = range(1, 10, 2)
        actual2 = lrange(1, 10, 2)
        actual += [actual1, actual2],
        lengths += 5,
        expected += list(builtins.range(1, 10, 2)),
        self.check_result(actual, expected, lengths)

    def test_map(self):
        func = lambda x, y, z: x + y + z
        lst = [builtins.range(10), builtins.range(10), builtins.range(10)]
        actual1 = map(func, *lst)
        actual2 = lmap(func, *lst)
        actual = [actual1, actual2],
        expected = list(builtins.map(func, *lst)),
        lengths = 10,
        self.check_result(actual, expected, lengths)

    def test_filter(self):
        func = lambda x: x
        lst = list(builtins.range(10))
        actual1 = filter(func, lst)
        actual2 = lfilter(func, lst)
        actual = [actual1, actual2],
        lengths = 9,
        expected = list(builtins.filter(func, lst)),
        self.check_result(actual, expected, lengths)

    def test_zip(self):
        lst = [builtins.range(10), builtins.range(10), builtins.range(10)]
        actual = [zip(*lst), lzip(*lst)],
        expected = list(builtins.zip(*lst)),
        lengths = 10,
        self.check_result(actual, expected, lengths)

    def test_dict_iterators(self):
        assert next(itervalues({1: 2})) == 2
        assert next(iterkeys({1: 2})) == 1
        assert next(iteritems({1: 2})) == (1, 2)


class TestCompatFunctions(object):

    @pytest.mark.parametrize(
        'start,stop,step', [(0, 10, 2), (11, -2, -1), (0, -5, 1), (2, 4, 8)])
    def test_get_range_parameters(self, start, stop, step):
        rng = range(start, stop, step)
        if PY2 and len(rng) == 0:
            start_expected, stop_expected, step_expected = 0, 0, 1
        elif PY2 and len(rng) == 1:
            start_expected, stop_expected, step_expected = start, start + 1, 1
        else:
            start_expected, stop_expected, step_expected = start, stop, step

        start_result, stop_result, step_result = get_range_parameters(rng)
        assert start_result == start_expected
        assert stop_result == stop_expected
        assert step_result == step_expected
