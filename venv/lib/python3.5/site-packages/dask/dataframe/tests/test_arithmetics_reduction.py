from datetime import datetime

import pytest
import numpy as np
import pandas as pd

import dask.dataframe as dd
from dask.dataframe.utils import assert_eq, assert_dask_graph, make_meta


@pytest.mark.slow
def test_arithmetics():
    dsk = {('x', 0): pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]},
                                  index=[0, 1, 3]),
           ('x', 1): pd.DataFrame({'a': [4, 5, 6], 'b': [3, 2, 1]},
                                  index=[5, 6, 8]),
           ('x', 2): pd.DataFrame({'a': [7, 8, 9], 'b': [0, 0, 0]},
                                  index=[9, 9, 9])}
    meta = make_meta({'a': 'i8', 'b': 'i8'}, index=pd.Index([], 'i8'))
    ddf1 = dd.DataFrame(dsk, 'x', meta, [0, 4, 9, 9])
    pdf1 = ddf1.compute()

    pdf2 = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6, 7, 8],
                         'b': [5, 6, 7, 8, 1, 2, 3, 4]})
    pdf3 = pd.DataFrame({'a': [5, 6, 7, 8, 4, 3, 2, 1],
                         'b': [2, 4, 5, 3, 4, 2, 1, 0]})
    ddf2 = dd.from_pandas(pdf2, 3)
    ddf3 = dd.from_pandas(pdf3, 2)

    dsk4 = {('y', 0): pd.DataFrame({'a': [3, 2, 1], 'b': [7, 8, 9]},
                                   index=[0, 1, 3]),
            ('y', 1): pd.DataFrame({'a': [5, 2, 8], 'b': [4, 2, 3]},
                                   index=[5, 6, 8]),
            ('y', 2): pd.DataFrame({'a': [1, 4, 10], 'b': [1, 0, 5]},
                                   index=[9, 9, 9])}
    ddf4 = dd.DataFrame(dsk4, 'y', meta, [0, 4, 9, 9])
    pdf4 = ddf4.compute()

    # Arithmetics
    cases = [(ddf1, ddf1, pdf1, pdf1),
             (ddf1, ddf1.repartition([0, 1, 3, 6, 9]), pdf1, pdf1),
             (ddf2, ddf3, pdf2, pdf3),
             (ddf2.repartition([0, 3, 6, 7]), ddf3.repartition([0, 7]),
              pdf2, pdf3),
             (ddf2.repartition([0, 7]), ddf3.repartition([0, 2, 4, 5, 7]),
              pdf2, pdf3),
             (ddf1, ddf4, pdf1, pdf4),
             (ddf1, ddf4.repartition([0, 9]), pdf1, pdf4),
             (ddf1.repartition([0, 3, 9]), ddf4.repartition([0, 5, 9]),
              pdf1, pdf4),
             # dask + pandas
             (ddf1, pdf4, pdf1, pdf4), (ddf2, pdf3, pdf2, pdf3)]

    for (l, r, el, er) in cases:
        check_series_arithmetics(l.a, r.b, el.a, er.b)
        check_frame_arithmetics(l, r, el, er)

    # different index, pandas raises ValueError in comparison ops

    pdf5 = pd.DataFrame({'a': [3, 2, 1, 5, 2, 8, 1, 4, 10],
                         'b': [7, 8, 9, 4, 2, 3, 1, 0, 5]},
                        index=[0, 1, 3, 5, 6, 8, 9, 9, 9])
    ddf5 = dd.from_pandas(pdf5, 2)

    pdf6 = pd.DataFrame({'a': [3, 2, 1, 5, 2, 8, 1, 4, 10],
                         'b': [7, 8, 9, 5, 7, 8, 4, 2, 5]},
                        index=[0, 1, 2, 3, 4, 5, 6, 7, 9])
    ddf6 = dd.from_pandas(pdf6, 4)

    pdf7 = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6, 7, 8],
                         'b': [5, 6, 7, 8, 1, 2, 3, 4]},
                        index=list('aaabcdeh'))
    pdf8 = pd.DataFrame({'a': [5, 6, 7, 8, 4, 3, 2, 1],
                         'b': [2, 4, 5, 3, 4, 2, 1, 0]},
                        index=list('abcdefgh'))
    ddf7 = dd.from_pandas(pdf7, 3)
    ddf8 = dd.from_pandas(pdf8, 4)

    pdf9 = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6, 7, 8],
                         'b': [5, 6, 7, 8, 1, 2, 3, 4],
                         'c': [5, 6, 7, 8, 1, 2, 3, 4]},
                        index=list('aaabcdeh'))
    pdf10 = pd.DataFrame({'b': [5, 6, 7, 8, 4, 3, 2, 1],
                          'c': [2, 4, 5, 3, 4, 2, 1, 0],
                          'd': [2, 4, 5, 3, 4, 2, 1, 0]},
                         index=list('abcdefgh'))
    ddf9 = dd.from_pandas(pdf9, 3)
    ddf10 = dd.from_pandas(pdf10, 4)

    # Arithmetics with different index
    cases = [(ddf5, ddf6, pdf5, pdf6),
             (ddf5.repartition([0, 9]), ddf6, pdf5, pdf6),
             (ddf5.repartition([0, 5, 9]), ddf6.repartition([0, 7, 9]),
              pdf5, pdf6),
             (ddf7, ddf8, pdf7, pdf8),
             (ddf7.repartition(['a', 'c', 'h']), ddf8.repartition(['a', 'h']),
              pdf7, pdf8),
             (ddf7.repartition(['a', 'b', 'e', 'h']),
              ddf8.repartition(['a', 'e', 'h']), pdf7, pdf8),
             (ddf9, ddf10, pdf9, pdf10),
             (ddf9.repartition(['a', 'c', 'h']), ddf10.repartition(['a', 'h']),
              pdf9, pdf10),
             # dask + pandas
             (ddf5, pdf6, pdf5, pdf6), (ddf7, pdf8, pdf7, pdf8),
             (ddf9, pdf10, pdf9, pdf10)]

    for (l, r, el, er) in cases:
        check_series_arithmetics(l.a, r.b, el.a, er.b,
                                 allow_comparison_ops=False)
        check_frame_arithmetics(l, r, el, er,
                                allow_comparison_ops=False)


def test_deterministic_arithmetic_names():
    df = pd.DataFrame({'x': [1, 2, 3, 4], 'y': [5, 6, 7, 8]})
    a = dd.from_pandas(df, npartitions=2)

    assert sorted((a.x + a.y ** 2).dask) == sorted((a.x + a.y ** 2).dask)
    assert sorted((a.x + a.y ** 2).dask) != sorted((a.x + a.y ** 3).dask)
    assert sorted((a.x + a.y ** 2).dask) != sorted((a.x - a.y ** 2).dask)


@pytest.mark.slow
def test_arithmetics_different_index():
    # index are different, but overwraps
    pdf1 = pd.DataFrame({'a': [1, 2, 3, 4, 5], 'b': [3, 5, 2, 5, 7]},
                        index=[1, 2, 3, 4, 5])
    ddf1 = dd.from_pandas(pdf1, 2)
    pdf2 = pd.DataFrame({'a': [3, 2, 6, 7, 8], 'b': [9, 4, 2, 6, 2]},
                        index=[3, 4, 5, 6, 7])
    ddf2 = dd.from_pandas(pdf2, 2)

    # index are not overwrapped
    pdf3 = pd.DataFrame({'a': [1, 2, 3, 4, 5], 'b': [3, 5, 2, 5, 7]},
                        index=[1, 2, 3, 4, 5])
    ddf3 = dd.from_pandas(pdf3, 2)
    pdf4 = pd.DataFrame({'a': [3, 2, 6, 7, 8], 'b': [9, 4, 2, 6, 2]},
                        index=[10, 11, 12, 13, 14])
    ddf4 = dd.from_pandas(pdf4, 2)

    # index is included in another
    pdf5 = pd.DataFrame({'a': [1, 2, 3, 4, 5], 'b': [3, 5, 2, 5, 7]},
                        index=[1, 3, 5, 7, 9])
    ddf5 = dd.from_pandas(pdf5, 2)
    pdf6 = pd.DataFrame({'a': [3, 2, 6, 7, 8], 'b': [9, 4, 2, 6, 2]},
                        index=[2, 3, 4, 5, 6])
    ddf6 = dd.from_pandas(pdf6, 2)

    cases = [(ddf1, ddf2, pdf1, pdf2),
             (ddf2, ddf1, pdf2, pdf1),
             (ddf1.repartition([1, 3, 5]), ddf2.repartition([3, 4, 7]),
              pdf1, pdf2),
             (ddf2.repartition([3, 4, 5, 7]), ddf1.repartition([1, 2, 4, 5]),
              pdf2, pdf1),
             (ddf3, ddf4, pdf3, pdf4),
             (ddf4, ddf3, pdf4, pdf3),
             (ddf3.repartition([1, 2, 3, 4, 5]),
              ddf4.repartition([10, 11, 12, 13, 14]), pdf3, pdf4),
             (ddf4.repartition([10, 14]), ddf3.repartition([1, 3, 4, 5]),
              pdf4, pdf3),
             (ddf5, ddf6, pdf5, pdf6),
             (ddf6, ddf5, pdf6, pdf5),
             (ddf5.repartition([1, 7, 8, 9]), ddf6.repartition([2, 3, 4, 6]),
              pdf5, pdf6),
             (ddf6.repartition([2, 6]), ddf5.repartition([1, 3, 7, 9]),
              pdf6, pdf5),
             # dask + pandas
             (ddf1, pdf2, pdf1, pdf2), (ddf2, pdf1, pdf2, pdf1),
             (ddf3, pdf4, pdf3, pdf4), (ddf4, pdf3, pdf4, pdf3),
             (ddf5, pdf6, pdf5, pdf6), (ddf6, pdf5, pdf6, pdf5)]

    for (l, r, el, er) in cases:
        check_series_arithmetics(l.a, r.b, el.a, er.b,
                                 allow_comparison_ops=False)
        check_frame_arithmetics(l, r, el, er,
                                allow_comparison_ops=False)

    pdf7 = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6, 7, 8],
                         'b': [5, 6, 7, 8, 1, 2, 3, 4]},
                        index=[0, 2, 4, 8, 9, 10, 11, 13])
    pdf8 = pd.DataFrame({'a': [5, 6, 7, 8, 4, 3, 2, 1],
                         'b': [2, 4, 5, 3, 4, 2, 1, 0]},
                        index=[1, 3, 4, 8, 9, 11, 12, 13])
    ddf7 = dd.from_pandas(pdf7, 3)
    ddf8 = dd.from_pandas(pdf8, 2)

    pdf9 = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6, 7, 8],
                         'b': [5, 6, 7, 8, 1, 2, 3, 4]},
                        index=[0, 2, 4, 8, 9, 10, 11, 13])
    pdf10 = pd.DataFrame({'a': [5, 6, 7, 8, 4, 3, 2, 1],
                          'b': [2, 4, 5, 3, 4, 2, 1, 0]},
                         index=[0, 3, 4, 8, 9, 11, 12, 13])
    ddf9 = dd.from_pandas(pdf9, 3)
    ddf10 = dd.from_pandas(pdf10, 2)

    cases = [(ddf7, ddf8, pdf7, pdf8),
             (ddf8, ddf7, pdf8, pdf7),
             # (ddf7.repartition([0, 13]),
             #  ddf8.repartition([0, 4, 11, 14], force=True),
             #  pdf7, pdf8),
             (ddf8.repartition([-5, 10, 15], force=True),
              ddf7.repartition([-1, 4, 11, 14], force=True), pdf8, pdf7),
             (ddf7.repartition([0, 8, 12, 13]),
              ddf8.repartition([0, 2, 8, 12, 13], force=True), pdf7, pdf8),
             (ddf8.repartition([-5, 0, 10, 20], force=True),
              ddf7.repartition([-1, 4, 11, 13], force=True), pdf8, pdf7),
             (ddf9, ddf10, pdf9, pdf10),
             (ddf10, ddf9, pdf10, pdf9),
             # dask + pandas
             (ddf7, pdf8, pdf7, pdf8), (ddf8, pdf7, pdf8, pdf7),
             (ddf9, pdf10, pdf9, pdf10), (ddf10, pdf9, pdf10, pdf9)]

    for (l, r, el, er) in cases:
        check_series_arithmetics(l.a, r.b, el.a, er.b,
                                 allow_comparison_ops=False)
        check_frame_arithmetics(l, r, el, er,
                                allow_comparison_ops=False)


def check_series_arithmetics(l, r, el, er, allow_comparison_ops=True):
    assert isinstance(l, dd.Series)
    assert isinstance(r, (dd.Series, pd.Series))
    assert isinstance(el, pd.Series)
    assert isinstance(er, pd.Series)

    # l, r may be repartitioned, test whether repartition keeps original data
    assert_eq(l, el)
    assert_eq(r, er)

    assert_eq(l + r, el + er)
    assert_eq(l * r, el * er)
    assert_eq(l - r, el - er)
    assert_eq(l / r, el / er)
    assert_eq(l // r, el // er)
    assert_eq(l ** r, el ** er)
    assert_eq(l % r, el % er)

    if allow_comparison_ops:
        # comparison is allowed if data have same index
        assert_eq(l & r, el & er)
        assert_eq(l | r, el | er)
        assert_eq(l ^ r, el ^ er)
        assert_eq(l > r, el > er)
        assert_eq(l < r, el < er)
        assert_eq(l >= r, el >= er)
        assert_eq(l <= r, el <= er)
        assert_eq(l == r, el == er)
        assert_eq(l != r, el != er)
        assert_eq(l.lt(r), el.lt(er))
        assert_eq(l.gt(r), el.gt(er))
        assert_eq(l.le(r), el.le(er))
        assert_eq(l.ge(r), el.ge(er))
        assert_eq(l.ne(r), el.ne(er))
        assert_eq(l.eq(r), el.eq(er))

    assert_eq(l + 2, el + 2)
    assert_eq(l * 2, el * 2)
    assert_eq(l - 2, el - 2)
    assert_eq(l / 2, el / 2)
    assert_eq(l & True, el & True)
    assert_eq(l | True, el | True)
    assert_eq(l ^ True, el ^ True)
    assert_eq(l // 2, el // 2)
    assert_eq(l ** 2, el ** 2)
    assert_eq(l % 2, el % 2)
    assert_eq(l > 2, el > 2)
    assert_eq(l < 2, el < 2)
    assert_eq(l >= 2, el >= 2)
    assert_eq(l <= 2, el <= 2)
    assert_eq(l == 2, el == 2)
    assert_eq(l != 2, el != 2)

    assert_eq(2 + r, 2 + er)
    assert_eq(2 * r, 2 * er)
    assert_eq(2 - r, 2 - er)
    assert_eq(2 / r, 2 / er)
    assert_eq(True & r, True & er)
    assert_eq(True | r, True | er)
    assert_eq(True ^ r, True ^ er)
    assert_eq(2 // r, 2 // er)
    assert_eq(2 ** r, 2 ** er)
    assert_eq(2 % r, 2 % er)
    assert_eq(2 > r, 2 > er)
    assert_eq(2 < r, 2 < er)
    assert_eq(2 >= r, 2 >= er)
    assert_eq(2 <= r, 2 <= er)
    assert_eq(2 == r, 2 == er)
    assert_eq(2 != r, 2 != er)

    assert_eq(l.lt(2), el.lt(2))
    assert_eq(l.gt(2), el.gt(2))
    assert_eq(l.le(2), el.le(2))
    assert_eq(l.ge(2), el.ge(2))
    assert_eq(l.ne(2), el.ne(2))
    assert_eq(l.eq(2), el.eq(2))

    assert_eq(-l, -el)
    assert_eq(abs(l), abs(el))

    if allow_comparison_ops:
        # comparison is allowed if data have same index
        assert_eq(~(l == r), ~(el == er))


def check_frame_arithmetics(l, r, el, er, allow_comparison_ops=True):
    assert isinstance(l, dd.DataFrame)
    assert isinstance(r, (dd.DataFrame, pd.DataFrame))
    assert isinstance(el, pd.DataFrame)
    assert isinstance(er, pd.DataFrame)
    # l, r may be repartitioned, test whether repartition keeps original data
    assert_eq(l, el)
    assert_eq(r, er)

    assert_eq(l + r, el + er)
    assert_eq(l * r, el * er)
    assert_eq(l - r, el - er)
    assert_eq(l / r, el / er)
    assert_eq(l // r, el // er)
    assert_eq(l ** r, el ** er)
    assert_eq(l % r, el % er)

    if allow_comparison_ops:
        # comparison is allowed if data have same index
        assert_eq(l & r, el & er)
        assert_eq(l | r, el | er)
        assert_eq(l ^ r, el ^ er)
        assert_eq(l > r, el > er)
        assert_eq(l < r, el < er)
        assert_eq(l >= r, el >= er)
        assert_eq(l <= r, el <= er)
        assert_eq(l == r, el == er)
        assert_eq(l != r, el != er)
        assert_eq(l.lt(r), el.lt(er))
        assert_eq(l.gt(r), el.gt(er))
        assert_eq(l.le(r), el.le(er))
        assert_eq(l.ge(r), el.ge(er))
        assert_eq(l.ne(r), el.ne(er))
        assert_eq(l.eq(r), el.eq(er))

    assert_eq(l + 2, el + 2)
    assert_eq(l * 2, el * 2)
    assert_eq(l - 2, el - 2)
    assert_eq(l / 2, el / 2)
    assert_eq(l & True, el & True)
    assert_eq(l | True, el | True)
    assert_eq(l ^ True, el ^ True)
    assert_eq(l // 2, el // 2)
    assert_eq(l ** 2, el ** 2)
    assert_eq(l % 2, el % 2)
    assert_eq(l > 2, el > 2)
    assert_eq(l < 2, el < 2)
    assert_eq(l >= 2, el >= 2)
    assert_eq(l <= 2, el <= 2)
    assert_eq(l == 2, el == 2)
    assert_eq(l != 2, el != 2)

    assert_eq(2 + l, 2 + el)
    assert_eq(2 * l, 2 * el)
    assert_eq(2 - l, 2 - el)
    assert_eq(2 / l, 2 / el)
    assert_eq(True & l, True & el)
    assert_eq(True | l, True | el)
    assert_eq(True ^ l, True ^ el)
    assert_eq(2 // l, 2 // el)
    assert_eq(2 ** l, 2 ** el)
    assert_eq(2 % l, 2 % el)
    assert_eq(2 > l, 2 > el)
    assert_eq(2 < l, 2 < el)
    assert_eq(2 >= l, 2 >= el)
    assert_eq(2 <= l, 2 <= el)
    assert_eq(2 == l, 2 == el)
    assert_eq(2 != l, 2 != el)

    assert_eq(l.lt(2), el.lt(2))
    assert_eq(l.gt(2), el.gt(2))
    assert_eq(l.le(2), el.le(2))
    assert_eq(l.ge(2), el.ge(2))
    assert_eq(l.ne(2), el.ne(2))
    assert_eq(l.eq(2), el.eq(2))

    assert_eq(-l, -el)
    assert_eq(abs(l), abs(el))

    if allow_comparison_ops:
        # comparison is allowed if data have same index
        assert_eq(~(l == r), ~(el == er))


def test_scalar_arithmetics():
    el = np.int64(10)
    er = np.int64(4)
    l = dd.core.Scalar({('l', 0): el}, 'l', 'i8')
    r = dd.core.Scalar({('r', 0): er}, 'r', 'i8')

    assert isinstance(l, dd.core.Scalar)
    assert isinstance(r, dd.core.Scalar)

    assert_eq(l, el)
    assert_eq(r, er)

    assert_eq(l + r, el + er)
    assert_eq(l * r, el * er)
    assert_eq(l - r, el - er)
    assert_eq(l / r, el / er)
    assert_eq(l // r, el // er)
    assert_eq(l ** r, el ** er)
    assert_eq(l % r, el % er)

    assert_eq(l & r, el & er)
    assert_eq(l | r, el | er)
    assert_eq(l ^ r, el ^ er)
    assert_eq(l > r, el > er)
    assert_eq(l < r, el < er)
    assert_eq(l >= r, el >= er)
    assert_eq(l <= r, el <= er)
    assert_eq(l == r, el == er)
    assert_eq(l != r, el != er)

    assert_eq(l + 2, el + 2)
    assert_eq(l * 2, el * 2)
    assert_eq(l - 2, el - 2)
    assert_eq(l / 2, el / 2)
    assert_eq(l & True, el & True)
    assert_eq(l | True, el | True)
    assert_eq(l ^ True, el ^ True)
    assert_eq(l // 2, el // 2)
    assert_eq(l ** 2, el ** 2)
    assert_eq(l % 2, el % 2)
    assert_eq(l > 2, el > 2)
    assert_eq(l < 2, el < 2)
    assert_eq(l >= 2, el >= 2)
    assert_eq(l <= 2, el <= 2)
    assert_eq(l == 2, el == 2)
    assert_eq(l != 2, el != 2)

    assert_eq(2 + r, 2 + er)
    assert_eq(2 * r, 2 * er)
    assert_eq(2 - r, 2 - er)
    assert_eq(2 / r, 2 / er)
    assert_eq(True & r, True & er)
    assert_eq(True | r, True | er)
    assert_eq(True ^ r, True ^ er)
    assert_eq(2 // r, 2 // er)
    assert_eq(2 ** r, 2 ** er)
    assert_eq(2 % r, 2 % er)
    assert_eq(2 > r, 2 > er)
    assert_eq(2 < r, 2 < er)
    assert_eq(2 >= r, 2 >= er)
    assert_eq(2 <= r, 2 <= er)
    assert_eq(2 == r, 2 == er)
    assert_eq(2 != r, 2 != er)

    assert_eq(-l, -el)
    assert_eq(abs(l), abs(el))

    assert_eq(~(l == r), ~(el == er))


def test_scalar_arithmetics_with_dask_instances():
    s = dd.core.Scalar({('s', 0): 10}, 's', 'i8')
    e = 10

    pds = pd.Series([1, 2, 3, 4, 5, 6, 7])
    dds = dd.from_pandas(pds, 2)

    pdf = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6, 7],
                        'b': [7, 6, 5, 4, 3, 2, 1]})
    ddf = dd.from_pandas(pdf, 2)

    # pandas Series
    result = pds + s   # this result pd.Series (automatically computed)
    assert isinstance(result, pd.Series)
    assert_eq(result, pds + e)

    result = s + pds   # this result dd.Series
    assert isinstance(result, dd.Series)
    assert_eq(result, pds + e)

    # dask Series
    result = dds + s   # this result dd.Series
    assert isinstance(result, dd.Series)
    assert_eq(result, pds + e)

    result = s + dds   # this result dd.Series
    assert isinstance(result, dd.Series)
    assert_eq(result, pds + e)

    # pandas DataFrame
    result = pdf + s   # this result pd.DataFrame (automatically computed)
    assert isinstance(result, pd.DataFrame)
    assert_eq(result, pdf + e)

    result = s + pdf   # this result dd.DataFrame
    assert isinstance(result, dd.DataFrame)
    assert_eq(result, pdf + e)

    # dask DataFrame
    result = ddf + s   # this result dd.DataFrame
    assert isinstance(result, dd.DataFrame)
    assert_eq(result, pdf + e)

    result = s + ddf   # this result dd.DataFrame
    assert isinstance(result, dd.DataFrame)
    assert_eq(result, pdf + e)


def test_frame_series_arithmetic_methods():
    pdf1 = pd.DataFrame({'A': np.arange(10),
                         'B': [np.nan, 1, 2, 3, 4] * 2,
                         'C': [np.nan] * 10,
                         'D': np.arange(10)},
                        index=list('abcdefghij'), columns=list('ABCD'))
    pdf2 = pd.DataFrame(np.random.randn(10, 4),
                        index=list('abcdefghjk'), columns=list('ABCX'))
    ps1 = pdf1.A
    ps2 = pdf2.A
    ps3 = pd.Series(np.random.randn(10), index=list('ABCDXabcde'))

    ddf1 = dd.from_pandas(pdf1, 2)
    ddf2 = dd.from_pandas(pdf2, 2)
    ds1 = ddf1.A
    ds2 = ddf2.A

    s = dd.core.Scalar({('s', 0): 4}, 's', 'i8')

    for l, r, el, er in [(ddf1, ddf2, pdf1, pdf2), (ds1, ds2, ps1, ps2),
                         (ddf1.repartition(['a', 'f', 'j']), ddf2, pdf1, pdf2),
                         (ds1.repartition(['a', 'b', 'f', 'j']), ds2, ps1, ps2),
                         (ddf1, ddf2.repartition(['a', 'k']), pdf1, pdf2),
                         (ds1, ds2.repartition(['a', 'b', 'd', 'h', 'k']), ps1, ps2),
                         (ddf1, 3, pdf1, 3), (ds1, 3, ps1, 3),
                         (ddf1, s, pdf1, 4), (ds1, s, ps1, 4)]:
        # l, r may be repartitioned, test whether repartition keeps original data
        assert_eq(l, el)
        assert_eq(r, er)

        assert_eq(l.add(r, fill_value=0), el.add(er, fill_value=0))
        assert_eq(l.sub(r, fill_value=0), el.sub(er, fill_value=0))
        assert_eq(l.mul(r, fill_value=0), el.mul(er, fill_value=0))
        assert_eq(l.div(r, fill_value=0), el.div(er, fill_value=0))
        assert_eq(l.truediv(r, fill_value=0), el.truediv(er, fill_value=0))
        assert_eq(l.floordiv(r, fill_value=1), el.floordiv(er, fill_value=1))
        assert_eq(l.mod(r, fill_value=0), el.mod(er, fill_value=0))
        assert_eq(l.pow(r, fill_value=0), el.pow(er, fill_value=0))

        assert_eq(l.radd(r, fill_value=0), el.radd(er, fill_value=0))
        assert_eq(l.rsub(r, fill_value=0), el.rsub(er, fill_value=0))
        assert_eq(l.rmul(r, fill_value=0), el.rmul(er, fill_value=0))
        assert_eq(l.rdiv(r, fill_value=0), el.rdiv(er, fill_value=0))
        assert_eq(l.rtruediv(r, fill_value=0), el.rtruediv(er, fill_value=0))
        assert_eq(l.rfloordiv(r, fill_value=1), el.rfloordiv(er, fill_value=1))
        assert_eq(l.rmod(r, fill_value=0), el.rmod(er, fill_value=0))
        assert_eq(l.rpow(r, fill_value=0), el.rpow(er, fill_value=0))

    for l, r, el, er in [(ddf1, ds2, pdf1, ps2), (ddf1, ddf2.X, pdf1, pdf2.X)]:
        assert_eq(l, el)
        assert_eq(r, er)

        # must specify axis=0 to add Series to each column
        # axis=1 is not supported (add to each row)
        assert_eq(l.add(r, axis=0), el.add(er, axis=0))
        assert_eq(l.sub(r, axis=0), el.sub(er, axis=0))
        assert_eq(l.mul(r, axis=0), el.mul(er, axis=0))
        assert_eq(l.div(r, axis=0), el.div(er, axis=0))
        assert_eq(l.truediv(r, axis=0), el.truediv(er, axis=0))
        assert_eq(l.floordiv(r, axis=0), el.floordiv(er, axis=0))
        assert_eq(l.mod(r, axis=0), el.mod(er, axis=0))
        assert_eq(l.pow(r, axis=0), el.pow(er, axis=0))

        assert_eq(l.radd(r, axis=0), el.radd(er, axis=0))
        assert_eq(l.rsub(r, axis=0), el.rsub(er, axis=0))
        assert_eq(l.rmul(r, axis=0), el.rmul(er, axis=0))
        assert_eq(l.rdiv(r, axis=0), el.rdiv(er, axis=0))
        assert_eq(l.rtruediv(r, axis=0), el.rtruediv(er, axis=0))
        assert_eq(l.rfloordiv(r, axis=0), el.rfloordiv(er, axis=0))
        assert_eq(l.rmod(r, axis=0), el.rmod(er, axis=0))
        assert_eq(l.rpow(r, axis=0), el.rpow(er, axis=0))

        pytest.raises(ValueError, lambda: l.add(r, axis=1))

    for l, r, el, er in [(ddf1, pdf2, pdf1, pdf2), (ddf1, ps3, pdf1, ps3)]:
        assert_eq(l, el)
        assert_eq(r, er)

        for axis in [0, 1, 'index', 'columns']:
            assert_eq(l.add(r, axis=axis), el.add(er, axis=axis))
            assert_eq(l.sub(r, axis=axis), el.sub(er, axis=axis))
            assert_eq(l.mul(r, axis=axis), el.mul(er, axis=axis))
            assert_eq(l.div(r, axis=axis), el.div(er, axis=axis))
            assert_eq(l.truediv(r, axis=axis), el.truediv(er, axis=axis))
            assert_eq(l.floordiv(r, axis=axis), el.floordiv(er, axis=axis))
            assert_eq(l.mod(r, axis=axis), el.mod(er, axis=axis))
            assert_eq(l.pow(r, axis=axis), el.pow(er, axis=axis))

            assert_eq(l.radd(r, axis=axis), el.radd(er, axis=axis))
            assert_eq(l.rsub(r, axis=axis), el.rsub(er, axis=axis))
            assert_eq(l.rmul(r, axis=axis), el.rmul(er, axis=axis))
            assert_eq(l.rdiv(r, axis=axis), el.rdiv(er, axis=axis))
            assert_eq(l.rtruediv(r, axis=axis), el.rtruediv(er, axis=axis))
            assert_eq(l.rfloordiv(r, axis=axis), el.rfloordiv(er, axis=axis))
            assert_eq(l.rmod(r, axis=axis), el.rmod(er, axis=axis))
            assert_eq(l.rpow(r, axis=axis), el.rpow(er, axis=axis))


@pytest.mark.parametrize('split_every', [False, 2])
def test_reductions(split_every):
    dsk = {('x', 0): pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]},
                                  index=[0, 1, 3]),
           ('x', 1): pd.DataFrame({'a': [4, 5, 6], 'b': [3, 2, 1]},
                                  index=[5, 6, 8]),
           ('x', 2): pd.DataFrame({'a': [7, 8, 9], 'b': [0, 0, 0]},
                                  index=[9, 9, 9])}
    meta = make_meta({'a': 'i8', 'b': 'i8'}, index=pd.Index([], 'i8'))
    ddf1 = dd.DataFrame(dsk, 'x', meta, [0, 4, 9, 9])
    pdf1 = ddf1.compute()

    nans1 = pd.Series([1] + [np.nan] * 4 + [2] + [np.nan] * 3)
    nands1 = dd.from_pandas(nans1, 2)
    nans2 = pd.Series([1] + [np.nan] * 8)
    nands2 = dd.from_pandas(nans2, 2)
    nans3 = pd.Series([np.nan] * 9)
    nands3 = dd.from_pandas(nans3, 2)

    bools = pd.Series([True, False, True, False, True], dtype=bool)
    boolds = dd.from_pandas(bools, 2)

    for dds, pds in [(ddf1.b, pdf1.b), (ddf1.a, pdf1.a),
                     (ddf1['a'], pdf1['a']), (ddf1['b'], pdf1['b']),
                     (nands1, nans1), (nands2, nans2), (nands3, nans3),
                     (boolds, bools)]:
        assert isinstance(dds, dd.Series)
        assert isinstance(pds, pd.Series)

        assert_eq(dds.sum(split_every=split_every), pds.sum())
        assert_eq(dds.prod(split_every=split_every), pds.prod())
        assert_eq(dds.min(split_every=split_every), pds.min())
        assert_eq(dds.max(split_every=split_every), pds.max())
        assert_eq(dds.count(split_every=split_every), pds.count())
        with pytest.warns(None):
            # runtime warnings; https://github.com/dask/dask/issues/2381
            assert_eq(dds.std(split_every=split_every), pds.std())
        with pytest.warns(None):
            # runtime warnings; https://github.com/dask/dask/issues/2381
            assert_eq(dds.var(split_every=split_every), pds.var())
        with pytest.warns(None):
            # runtime warnings; https://github.com/dask/dask/issues/2381
            assert_eq(dds.sem(split_every=split_every), pds.sem())
        assert_eq(dds.std(ddof=0, split_every=split_every), pds.std(ddof=0))
        assert_eq(dds.var(ddof=0, split_every=split_every), pds.var(ddof=0))
        assert_eq(dds.sem(ddof=0, split_every=split_every), pds.sem(ddof=0))
        assert_eq(dds.mean(split_every=split_every), pds.mean())
        assert_eq(dds.nunique(split_every=split_every), pds.nunique())

        assert_eq(dds.sum(skipna=False, split_every=split_every),
                  pds.sum(skipna=False))
        assert_eq(dds.prod(skipna=False, split_every=split_every),
                  pds.prod(skipna=False))
        assert_eq(dds.min(skipna=False, split_every=split_every),
                  pds.min(skipna=False))
        assert_eq(dds.max(skipna=False, split_every=split_every),
                  pds.max(skipna=False))
        assert_eq(dds.std(skipna=False, split_every=split_every),
                  pds.std(skipna=False))
        assert_eq(dds.var(skipna=False, split_every=split_every),
                  pds.var(skipna=False))
        assert_eq(dds.sem(skipna=False, split_every=split_every),
                  pds.sem(skipna=False))
        assert_eq(dds.std(skipna=False, ddof=0, split_every=split_every),
                  pds.std(skipna=False, ddof=0))
        assert_eq(dds.var(skipna=False, ddof=0, split_every=split_every),
                  pds.var(skipna=False, ddof=0))
        assert_eq(dds.sem(skipna=False, ddof=0, split_every=split_every),
                  pds.sem(skipna=False, ddof=0))
        assert_eq(dds.mean(skipna=False, split_every=split_every),
                  pds.mean(skipna=False))

    assert_dask_graph(ddf1.b.sum(split_every=split_every), 'series-sum')
    assert_dask_graph(ddf1.b.prod(split_every=split_every), 'series-prod')
    assert_dask_graph(ddf1.b.min(split_every=split_every), 'series-min')
    assert_dask_graph(ddf1.b.max(split_every=split_every), 'series-max')
    assert_dask_graph(ddf1.b.count(split_every=split_every), 'series-count')
    assert_dask_graph(ddf1.b.std(split_every=split_every), 'series-std')
    assert_dask_graph(ddf1.b.var(split_every=split_every), 'series-var')
    assert_dask_graph(ddf1.b.sem(split_every=split_every), 'series-sem')
    assert_dask_graph(ddf1.b.std(ddof=0, split_every=split_every), 'series-std')
    assert_dask_graph(ddf1.b.var(ddof=0, split_every=split_every), 'series-var')
    assert_dask_graph(ddf1.b.sem(ddof=0, split_every=split_every), 'series-sem')
    assert_dask_graph(ddf1.b.mean(split_every=split_every), 'series-mean')
    # nunique is performed using drop-duplicates
    assert_dask_graph(ddf1.b.nunique(split_every=split_every), 'drop-duplicates')

    # testing index
    assert_eq(ddf1.index.min(split_every=split_every), pdf1.index.min())
    assert_eq(ddf1.index.max(split_every=split_every), pdf1.index.max())
    assert_eq(ddf1.index.count(split_every=split_every), pd.notnull(pdf1.index).sum())


@pytest.mark.parametrize('frame,axis,out',
                         [(pd.DataFrame({'a': [1, 2, 3],
                                         'b': [4, 5, 6]},
                                        index=[0, 1, 3]),
                          0, pd.Series([])),
                          (pd.DataFrame({'a': [1, 2, 3],
                                        'b': [4, 5, 6]},
                                        index=[0, 1, 3]),
                           1, pd.Series([])),
                          (pd.Series([1, 2.5, 6]), None, None)])
@pytest.mark.parametrize('redfunc', ['sum', 'prod', 'min', 'max', 'mean', 'var', 'std'])
def test_reductions_out(frame, axis, out, redfunc):
    dsk_in = dd.from_pandas(frame, 3)
    dsk_out = dd.from_pandas(pd.Series([0]), 1).sum()

    if out is not None:
        dsk_out = dd.from_pandas(out, 3)

    np_redfunc = getattr(np, redfunc)
    pd_redfunc = getattr(frame.__class__, redfunc)
    dsk_redfunc = getattr(dsk_in.__class__, redfunc)

    if redfunc in ['var', 'std']:
        # numpy has default ddof value 0 while
        # dask and pandas have 1, so ddof should be passed
        # explicitly when calling np.var(dask)
        np_redfunc(dsk_in, axis=axis, ddof=1, out=dsk_out)
    else:
        np_redfunc(dsk_in, axis=axis, out=dsk_out)

    assert_eq(dsk_out, pd_redfunc(frame, axis=axis))

    dsk_redfunc(dsk_in, axis=axis, split_every=False, out=dsk_out)
    assert_eq(dsk_out, pd_redfunc(frame, axis=axis))

    dsk_redfunc(dsk_in, axis=axis, split_every=2, out=dsk_out)
    assert_eq(dsk_out, pd_redfunc(frame, axis=axis))


@pytest.mark.parametrize('split_every', [False, 2])
def test_allany(split_every):
    df = pd.DataFrame(np.random.choice([True, False], size=(100, 4)),
                      columns=['A', 'B', 'C', 'D'])
    df['E'] = list('abcde') * 20
    ddf = dd.from_pandas(df, 10)

    assert_eq(ddf.all(split_every=split_every), df.all())
    assert_eq(ddf.all(axis=1, split_every=split_every), df.all(axis=1))
    assert_eq(ddf.all(axis=0, split_every=split_every), df.all(axis=0))

    assert_eq(ddf.any(split_every=split_every), df.any())
    assert_eq(ddf.any(axis=1, split_every=split_every), df.any(axis=1))
    assert_eq(ddf.any(axis=0, split_every=split_every), df.any(axis=0))

    assert_eq(ddf.A.all(split_every=split_every), df.A.all())
    assert_eq(ddf.A.any(split_every=split_every), df.A.any())

    # testing numpy functions with out param
    ddf_out_axis_default = dd.from_pandas(pd.Series([False, False, False, False, False],
                                                    index=['A', 'B', 'C', 'D', 'E']), 10)
    ddf_out_axis1 = dd.from_pandas(pd.Series(np.random.choice([True, False], size=(100,))), 10)

    # all
    ddf.all(split_every=split_every, out=ddf_out_axis_default)
    assert_eq(ddf_out_axis_default, df.all())

    ddf.all(axis=1, split_every=split_every, out=ddf_out_axis1)
    assert_eq(ddf_out_axis1, df.all(axis=1))

    ddf.all(split_every=split_every, axis=0, out=ddf_out_axis_default)
    assert_eq(ddf_out_axis_default, df.all(axis=0))

    # any
    ddf.any(split_every=split_every, out=ddf_out_axis_default)
    assert_eq(ddf_out_axis_default, df.any())

    ddf.any(axis=1, split_every=split_every, out=ddf_out_axis1)
    assert_eq(ddf_out_axis1, df.any(axis=1))

    ddf.any(split_every=split_every, axis=0, out=ddf_out_axis_default)
    assert_eq(ddf_out_axis_default, df.any(axis=0))


@pytest.mark.parametrize('split_every', [False, 2])
def test_deterministic_reduction_names(split_every):
    df = pd.DataFrame({'x': [1, 2, 3, 4], 'y': [5, 6, 7, 8]})
    ddf = dd.from_pandas(df, npartitions=2)

    for x in [ddf, ddf.x]:
        assert (x.sum(split_every=split_every)._name ==
                x.sum(split_every=split_every)._name)
        assert (x.prod(split_every=split_every)._name ==
                x.prod(split_every=split_every)._name)
        assert (x.min(split_every=split_every)._name ==
                x.min(split_every=split_every)._name)
        assert (x.max(split_every=split_every)._name ==
                x.max(split_every=split_every)._name)
        assert (x.count(split_every=split_every)._name ==
                x.count(split_every=split_every)._name)
        assert (x.std(split_every=split_every)._name ==
                x.std(split_every=split_every)._name)
        assert (x.var(split_every=split_every)._name ==
                x.var(split_every=split_every)._name)
        assert (x.sem(split_every=split_every)._name ==
                x.sem(split_every=split_every)._name)
        assert (x.mean(split_every=split_every)._name ==
                x.mean(split_every=split_every)._name)

    assert (ddf.x.nunique(split_every=split_every)._name ==
            ddf.x.nunique(split_every=split_every)._name)


def test_reduction_series_invalid_axis():
    dsk = {('x', 0): pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]},
                                  index=[0, 1, 3]),
           ('x', 1): pd.DataFrame({'a': [4, 5, 6], 'b': [3, 2, 1]},
                                  index=[5, 6, 8]),
           ('x', 2): pd.DataFrame({'a': [7, 8, 9], 'b': [0, 0, 0]},
                                  index=[9, 9, 9])}
    meta = make_meta({'a': 'i8', 'b': 'i8'}, index=pd.Index([], 'i8'))
    ddf1 = dd.DataFrame(dsk, 'x', meta, [0, 4, 9, 9])
    pdf1 = ddf1.compute()

    for axis in [1, 'columns']:
        for s in [ddf1.a, pdf1.a]:    # both must behave the same
            pytest.raises(ValueError, lambda: s.sum(axis=axis))
            pytest.raises(ValueError, lambda: s.prod(axis=axis))
            pytest.raises(ValueError, lambda: s.min(axis=axis))
            pytest.raises(ValueError, lambda: s.max(axis=axis))
            # only count doesn't have axis keyword
            pytest.raises(TypeError, lambda: s.count(axis=axis))
            pytest.raises(ValueError, lambda: s.std(axis=axis))
            pytest.raises(ValueError, lambda: s.var(axis=axis))
            pytest.raises(ValueError, lambda: s.sem(axis=axis))
            pytest.raises(ValueError, lambda: s.mean(axis=axis))


def test_reductions_non_numeric_dtypes():
    # test non-numric blocks

    def check_raises(d, p, func):
        pytest.raises((TypeError, ValueError),
                      lambda: getattr(d, func)().compute())
        pytest.raises((TypeError, ValueError),
                      lambda: getattr(p, func)())

    pds = pd.Series(['a', 'b', 'c', 'd', 'e'])
    dds = dd.from_pandas(pds, 2)
    assert_eq(dds.sum(), pds.sum())
    check_raises(dds, pds, 'prod')
    assert_eq(dds.min(), pds.min())
    assert_eq(dds.max(), pds.max())
    assert_eq(dds.count(), pds.count())
    check_raises(dds, pds, 'std')
    check_raises(dds, pds, 'var')
    check_raises(dds, pds, 'sem')
    check_raises(dds, pds, 'mean')
    assert_eq(dds.nunique(), pds.nunique())

    for pds in [pd.Series(pd.Categorical([1, 2, 3, 4, 5], ordered=True)),
                pd.Series(pd.Categorical(list('abcde'), ordered=True)),
                pd.Series(pd.date_range('2011-01-01', freq='D', periods=5))]:
        dds = dd.from_pandas(pds, 2)

        check_raises(dds, pds, 'sum')
        check_raises(dds, pds, 'prod')
        assert_eq(dds.min(), pds.min())
        assert_eq(dds.max(), pds.max())
        assert_eq(dds.count(), pds.count())
        check_raises(dds, pds, 'std')
        check_raises(dds, pds, 'var')
        check_raises(dds, pds, 'sem')
        check_raises(dds, pds, 'mean')
        assert_eq(dds.nunique(), pds.nunique())

    pds = pd.Series(pd.timedelta_range('1 days', freq='D', periods=5))
    dds = dd.from_pandas(pds, 2)
    assert_eq(dds.sum(), pds.sum())
    assert_eq(dds.min(), pds.min())
    assert_eq(dds.max(), pds.max())
    assert_eq(dds.count(), pds.count())

    # ToDo: pandas supports timedelta std, otherwise dask raises:
    # incompatible type for a datetime/timedelta operation [__pow__]
    # assert_eq(dds.std(), pds.std())
    # assert_eq(dds.var(), pds.var())

    # ToDo: pandas supports timedelta std, otherwise dask raises:
    # TypeError: unsupported operand type(s) for *: 'float' and 'Timedelta'
    # assert_eq(dds.mean(), pds.mean())

    assert_eq(dds.nunique(), pds.nunique())


@pytest.mark.parametrize('split_every', [False, 2])
def test_reductions_frame(split_every):
    dsk = {('x', 0): pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]},
                                  index=[0, 1, 3]),
           ('x', 1): pd.DataFrame({'a': [4, 5, 6], 'b': [3, 2, 1]},
                                  index=[5, 6, 8]),
           ('x', 2): pd.DataFrame({'a': [7, 8, 9], 'b': [0, 0, 0]},
                                  index=[9, 9, 9])}
    meta = make_meta({'a': 'i8', 'b': 'i8'}, index=pd.Index([], 'i8'))
    ddf1 = dd.DataFrame(dsk, 'x', meta, [0, 4, 9, 9])
    pdf1 = ddf1.compute()

    assert_eq(ddf1.sum(split_every=split_every), pdf1.sum())
    assert_eq(ddf1.prod(split_every=split_every), pdf1.prod())
    assert_eq(ddf1.min(split_every=split_every), pdf1.min())
    assert_eq(ddf1.max(split_every=split_every), pdf1.max())
    assert_eq(ddf1.count(split_every=split_every), pdf1.count())
    assert_eq(ddf1.std(split_every=split_every), pdf1.std())
    assert_eq(ddf1.var(split_every=split_every), pdf1.var())
    assert_eq(ddf1.sem(split_every=split_every), pdf1.sem())
    assert_eq(ddf1.std(ddof=0, split_every=split_every), pdf1.std(ddof=0))
    assert_eq(ddf1.var(ddof=0, split_every=split_every), pdf1.var(ddof=0))
    assert_eq(ddf1.sem(ddof=0, split_every=split_every), pdf1.sem(ddof=0))
    assert_eq(ddf1.mean(split_every=split_every), pdf1.mean())

    for axis in [0, 1, 'index', 'columns']:
        assert_eq(ddf1.sum(axis=axis, split_every=split_every),
                  pdf1.sum(axis=axis))
        assert_eq(ddf1.prod(axis=axis, split_every=split_every),
                  pdf1.prod(axis=axis))
        assert_eq(ddf1.min(axis=axis, split_every=split_every),
                  pdf1.min(axis=axis))
        assert_eq(ddf1.max(axis=axis, split_every=split_every),
                  pdf1.max(axis=axis))
        assert_eq(ddf1.count(axis=axis, split_every=split_every),
                  pdf1.count(axis=axis))
        assert_eq(ddf1.std(axis=axis, split_every=split_every),
                  pdf1.std(axis=axis))
        assert_eq(ddf1.var(axis=axis, split_every=split_every),
                  pdf1.var(axis=axis))
        assert_eq(ddf1.sem(axis=axis, split_every=split_every),
                  pdf1.sem(axis=axis))
        assert_eq(ddf1.std(axis=axis, ddof=0, split_every=split_every),
                  pdf1.std(axis=axis, ddof=0))
        assert_eq(ddf1.var(axis=axis, ddof=0, split_every=split_every),
                  pdf1.var(axis=axis, ddof=0))
        assert_eq(ddf1.sem(axis=axis, ddof=0, split_every=split_every),
                  pdf1.sem(axis=axis, ddof=0))
        assert_eq(ddf1.mean(axis=axis, split_every=split_every),
                  pdf1.mean(axis=axis))

    pytest.raises(ValueError, lambda: ddf1.sum(axis='incorrect').compute())

    # axis=0
    assert_dask_graph(ddf1.sum(split_every=split_every), 'dataframe-sum')
    assert_dask_graph(ddf1.prod(split_every=split_every), 'dataframe-prod')
    assert_dask_graph(ddf1.min(split_every=split_every), 'dataframe-min')
    assert_dask_graph(ddf1.max(split_every=split_every), 'dataframe-max')
    assert_dask_graph(ddf1.count(split_every=split_every), 'dataframe-count')
    # std, var, sem, and mean consist of sum and count operations
    assert_dask_graph(ddf1.std(split_every=split_every), 'dataframe-sum')
    assert_dask_graph(ddf1.std(split_every=split_every), 'dataframe-count')
    assert_dask_graph(ddf1.var(split_every=split_every), 'dataframe-sum')
    assert_dask_graph(ddf1.var(split_every=split_every), 'dataframe-count')
    assert_dask_graph(ddf1.sem(split_every=split_every), 'dataframe-sum')
    assert_dask_graph(ddf1.sem(split_every=split_every), 'dataframe-count')
    assert_dask_graph(ddf1.mean(split_every=split_every), 'dataframe-sum')
    assert_dask_graph(ddf1.mean(split_every=split_every), 'dataframe-count')

    # axis=1
    assert_dask_graph(ddf1.sum(axis=1, split_every=split_every),
                      'dataframe-sum')
    assert_dask_graph(ddf1.prod(axis=1, split_every=split_every),
                      'dataframe-prod')
    assert_dask_graph(ddf1.min(axis=1, split_every=split_every),
                      'dataframe-min')
    assert_dask_graph(ddf1.max(axis=1, split_every=split_every),
                      'dataframe-max')
    assert_dask_graph(ddf1.count(axis=1, split_every=split_every),
                      'dataframe-count')
    assert_dask_graph(ddf1.std(axis=1, split_every=split_every),
                      'dataframe-std')
    assert_dask_graph(ddf1.var(axis=1, split_every=split_every),
                      'dataframe-var')
    assert_dask_graph(ddf1.sem(axis=1, split_every=split_every),
                      'dataframe-sem')
    assert_dask_graph(ddf1.mean(axis=1, split_every=split_every),
                      'dataframe-mean')


def test_reductions_frame_dtypes():
    df = pd.DataFrame({'int': [1, 2, 3, 4, 5, 6, 7, 8],
                       'float': [1., 2., 3., 4., np.nan, 6., 7., 8.],
                       'dt': [pd.NaT] + [datetime(2011, i, 1) for i in range(1, 8)],
                       'str': list('abcdefgh')})
    ddf = dd.from_pandas(df, 3)
    assert_eq(df.sum(), ddf.sum())
    assert_eq(df.prod(), ddf.prod())
    assert_eq(df.min(), ddf.min())
    assert_eq(df.max(), ddf.max())
    assert_eq(df.count(), ddf.count())
    assert_eq(df.std(), ddf.std())
    assert_eq(df.var(), ddf.var())
    assert_eq(df.sem(), ddf.sem())
    assert_eq(df.std(ddof=0), ddf.std(ddof=0))
    assert_eq(df.var(ddof=0), ddf.var(ddof=0))
    assert_eq(df.sem(ddof=0), ddf.sem(ddof=0))
    assert_eq(df.mean(), ddf.mean())

    assert_eq(df._get_numeric_data(), ddf._get_numeric_data())

    numerics = ddf[['int', 'float']]
    assert numerics._get_numeric_data().dask == numerics.dask


@pytest.mark.parametrize('split_every', [False, 2])
def test_reductions_frame_nan(split_every):
    df = pd.DataFrame({'a': [1, 2, np.nan, 4, 5, 6, 7, 8],
                       'b': [1, 2, np.nan, np.nan, np.nan, 5, np.nan, np.nan],
                       'c': [np.nan] * 8})
    ddf = dd.from_pandas(df, 3)
    assert_eq(df.sum(), ddf.sum(split_every=split_every))
    assert_eq(df.prod(), ddf.prod(split_every=split_every))
    assert_eq(df.min(), ddf.min(split_every=split_every))
    assert_eq(df.max(), ddf.max(split_every=split_every))
    assert_eq(df.count(), ddf.count(split_every=split_every))
    assert_eq(df.std(), ddf.std(split_every=split_every))
    assert_eq(df.var(), ddf.var(split_every=split_every))
    assert_eq(df.sem(), ddf.sem(split_every=split_every))
    assert_eq(df.std(ddof=0), ddf.std(ddof=0, split_every=split_every))
    assert_eq(df.var(ddof=0), ddf.var(ddof=0, split_every=split_every))
    assert_eq(df.sem(ddof=0), ddf.sem(ddof=0, split_every=split_every))
    assert_eq(df.mean(), ddf.mean(split_every=split_every))

    assert_eq(df.sum(skipna=False),
              ddf.sum(skipna=False, split_every=split_every))
    assert_eq(df.prod(skipna=False),
              ddf.prod(skipna=False, split_every=split_every))
    assert_eq(df.min(skipna=False),
              ddf.min(skipna=False, split_every=split_every))
    assert_eq(df.max(skipna=False),
              ddf.max(skipna=False, split_every=split_every))
    assert_eq(df.std(skipna=False),
              ddf.std(skipna=False, split_every=split_every))
    assert_eq(df.var(skipna=False),
              ddf.var(skipna=False, split_every=split_every))
    assert_eq(df.sem(skipna=False),
              ddf.sem(skipna=False, split_every=split_every))
    assert_eq(df.std(skipna=False, ddof=0),
              ddf.std(skipna=False, ddof=0, split_every=split_every))
    assert_eq(df.var(skipna=False, ddof=0),
              ddf.var(skipna=False, ddof=0, split_every=split_every))
    assert_eq(df.sem(skipna=False, ddof=0),
              ddf.sem(skipna=False, ddof=0, split_every=split_every))
    assert_eq(df.mean(skipna=False),
              ddf.mean(skipna=False, split_every=split_every))

    assert_eq(df.sum(axis=1, skipna=False),
              ddf.sum(axis=1, skipna=False, split_every=split_every))
    assert_eq(df.prod(axis=1, skipna=False),
              ddf.prod(axis=1, skipna=False, split_every=split_every))
    assert_eq(df.min(axis=1, skipna=False),
              ddf.min(axis=1, skipna=False, split_every=split_every))
    assert_eq(df.max(axis=1, skipna=False),
              ddf.max(axis=1, skipna=False, split_every=split_every))
    assert_eq(df.std(axis=1, skipna=False),
              ddf.std(axis=1, skipna=False, split_every=split_every))
    assert_eq(df.var(axis=1, skipna=False),
              ddf.var(axis=1, skipna=False, split_every=split_every))
    assert_eq(df.sem(axis=1, skipna=False),
              ddf.sem(axis=1, skipna=False, split_every=split_every))
    assert_eq(df.std(axis=1, skipna=False, ddof=0),
              ddf.std(axis=1, skipna=False, ddof=0, split_every=split_every))
    assert_eq(df.var(axis=1, skipna=False, ddof=0),
              ddf.var(axis=1, skipna=False, ddof=0, split_every=split_every))
    assert_eq(df.sem(axis=1, skipna=False, ddof=0),
              ddf.sem(axis=1, skipna=False, ddof=0, split_every=split_every))
    assert_eq(df.mean(axis=1, skipna=False),
              ddf.mean(axis=1, skipna=False, split_every=split_every))
