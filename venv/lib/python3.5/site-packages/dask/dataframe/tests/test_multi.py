import dask.dataframe as dd
import numpy as np
import pandas as pd
import pandas.util.testing as tm

from dask.local import get_sync
from dask.base import compute_as_if_collection
from dask.dataframe.core import _Frame
from dask.dataframe.methods import concat
from dask.dataframe.multi import (align_partitions, merge_indexed_dataframes,
                                  hash_join, concat_indexed_dataframes,
                                  _maybe_align_partitions)
from dask.dataframe.utils import (assert_eq, assert_divisions, make_meta,
                                  has_known_categories, clear_known_categories)

import pytest


def test_align_partitions():
    A = pd.DataFrame({'x': [1, 2, 3, 4, 5, 6], 'y': list('abdabd')},
                     index=[10, 20, 30, 40, 50, 60])
    a = dd.repartition(A, [10, 40, 60])

    B = pd.DataFrame({'x': [1, 2, 3, 4], 'y': list('abda')},
                     index=[30, 70, 80, 100])
    b = dd.repartition(B, [30, 80, 100])

    s = dd.core.Scalar({('s', 0): 10}, 's', 'i8')

    (aa, bb), divisions, L = align_partitions(a, b)

    def _check(a, b, aa, bb):
        assert isinstance(a, dd.DataFrame)
        assert isinstance(b, dd.DataFrame)
        assert isinstance(aa, dd.DataFrame)
        assert isinstance(bb, dd.DataFrame)
        assert_eq(a, aa)
        assert_eq(b, bb)
        assert divisions == (10, 30, 40, 60, 80, 100)
        assert isinstance(L, list)
        assert len(divisions) == 1 + len(L)

    _check(a, b, aa, bb)
    assert L == [[(aa._name, 0), (bb._name, 0)],
                 [(aa._name, 1), (bb._name, 1)],
                 [(aa._name, 2), (bb._name, 2)],
                 [(aa._name, 3), (bb._name, 3)],
                 [(aa._name, 4), (bb._name, 4)]]

    (aa, ss, bb), divisions, L = align_partitions(a, s, b)
    _check(a, b, aa, bb)
    assert L == [[(aa._name, 0), None, (bb._name, 0)],
                 [(aa._name, 1), None, (bb._name, 1)],
                 [(aa._name, 2), None, (bb._name, 2)],
                 [(aa._name, 3), None, (bb._name, 3)],
                 [(aa._name, 4), None, (bb._name, 4)]]
    assert_eq(ss, 10)

    ldf = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6, 7],
                        'b': [7, 6, 5, 4, 3, 2, 1]})
    rdf = pd.DataFrame({'c': [1, 2, 3, 4, 5, 6, 7],
                        'd': [7, 6, 5, 4, 3, 2, 1]})

    for lhs, rhs in [(dd.from_pandas(ldf, 1), dd.from_pandas(rdf, 1)),
                     (dd.from_pandas(ldf, 2), dd.from_pandas(rdf, 2)),
                     (dd.from_pandas(ldf, 2), dd.from_pandas(rdf, 3)),
                     (dd.from_pandas(ldf, 3), dd.from_pandas(rdf, 2))]:
        (lresult, rresult), div, parts = align_partitions(lhs, rhs)
        assert_eq(lresult, ldf)
        assert_eq(rresult, rdf)

    # different index
    ldf = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6, 7],
                        'b': [7, 6, 5, 4, 3, 2, 1]},
                       index=list('abcdefg'))
    rdf = pd.DataFrame({'c': [1, 2, 3, 4, 5, 6, 7],
                        'd': [7, 6, 5, 4, 3, 2, 1]},
                       index=list('fghijkl'))

    for lhs, rhs in [(dd.from_pandas(ldf, 1), dd.from_pandas(rdf, 1)),
                     (dd.from_pandas(ldf, 2), dd.from_pandas(rdf, 2)),
                     (dd.from_pandas(ldf, 2), dd.from_pandas(rdf, 3)),
                     (dd.from_pandas(ldf, 3), dd.from_pandas(rdf, 2))]:
        (lresult, rresult), div, parts = align_partitions(lhs, rhs)
        assert_eq(lresult, ldf)
        assert_eq(rresult, rdf)


def test_align_partitions_unknown_divisions():
    df = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6, 7],
                       'b': [7, 6, 5, 4, 3, 2, 1]})
    # One known, one unknown
    ddf = dd.from_pandas(df, npartitions=2)
    ddf2 = dd.from_pandas(df, npartitions=2, sort=False)
    assert not ddf2.known_divisions

    with pytest.raises(ValueError):
        align_partitions(ddf, ddf2)

    # Both unknown
    ddf = dd.from_pandas(df + 1, npartitions=2, sort=False)
    ddf2 = dd.from_pandas(df, npartitions=2, sort=False)
    assert not ddf.known_divisions
    assert not ddf2.known_divisions

    with pytest.raises(ValueError):
        align_partitions(ddf, ddf2)


def test__maybe_align_partitions():
    df = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6, 7],
                       'b': [7, 6, 5, 4, 3, 2, 1]})
    # Both known, same divisions
    ddf = dd.from_pandas(df + 1, npartitions=2)
    ddf2 = dd.from_pandas(df, npartitions=2)

    a, b = _maybe_align_partitions([ddf, ddf2])
    assert a is ddf
    assert b is ddf2

    # Both unknown, same divisions
    ddf = dd.from_pandas(df + 1, npartitions=2, sort=False)
    ddf2 = dd.from_pandas(df, npartitions=2, sort=False)
    assert not ddf.known_divisions
    assert not ddf2.known_divisions

    a, b = _maybe_align_partitions([ddf, ddf2])
    assert a is ddf
    assert b is ddf2

    # Both known, different divisions
    ddf = dd.from_pandas(df + 1, npartitions=2)
    ddf2 = dd.from_pandas(df, npartitions=3)

    a, b = _maybe_align_partitions([ddf, ddf2])
    assert a.divisions == b.divisions

    # Both unknown, different divisions
    ddf = dd.from_pandas(df + 1, npartitions=2, sort=False)
    ddf2 = dd.from_pandas(df, npartitions=3, sort=False)
    assert not ddf.known_divisions
    assert not ddf2.known_divisions

    with pytest.raises(ValueError):
        _maybe_align_partitions([ddf, ddf2])

    # One known, one unknown
    ddf = dd.from_pandas(df, npartitions=2)
    ddf2 = dd.from_pandas(df, npartitions=2, sort=False)
    assert not ddf2.known_divisions

    with pytest.raises(ValueError):
        _maybe_align_partitions([ddf, ddf2])


def test_merge_indexed_dataframe_to_indexed_dataframe():
    A = pd.DataFrame({'x': [1, 2, 3, 4, 5, 6]},
                     index=[1, 2, 3, 4, 6, 7])
    a = dd.repartition(A, [1, 4, 7])

    B = pd.DataFrame({'y': list('abcdef')},
                     index=[1, 2, 4, 5, 6, 8])
    b = dd.repartition(B, [1, 2, 5, 8])

    c = merge_indexed_dataframes(a, b, how='left')
    assert c.divisions[0] == a.divisions[0]
    assert c.divisions[-1] == max(a.divisions + b.divisions)
    assert_eq(c, A.join(B))

    c = merge_indexed_dataframes(a, b, how='right')
    assert c.divisions[0] == b.divisions[0]
    assert c.divisions[-1] == b.divisions[-1]
    assert_eq(c, A.join(B, how='right'))

    c = merge_indexed_dataframes(a, b, how='inner')
    assert c.divisions[0] == 1
    assert c.divisions[-1] == max(a.divisions + b.divisions)
    assert_eq(c.compute(), A.join(B, how='inner'))

    c = merge_indexed_dataframes(a, b, how='outer')
    assert c.divisions[0] == 1
    assert c.divisions[-1] == 8
    assert_eq(c.compute(), A.join(B, how='outer'))

    assert (sorted(merge_indexed_dataframes(a, b, how='inner').dask) ==
            sorted(merge_indexed_dataframes(a, b, how='inner').dask))
    assert (sorted(merge_indexed_dataframes(a, b, how='inner').dask) !=
            sorted(merge_indexed_dataframes(a, b, how='outer').dask))


def list_eq(aa, bb):
    if isinstance(aa, dd.DataFrame):
        a = aa.compute(get=get_sync)
    else:
        a = aa
    if isinstance(bb, dd.DataFrame):
        b = bb.compute(get=get_sync)
    else:
        b = bb
    tm.assert_index_equal(a.columns, b.columns)

    if isinstance(a, pd.DataFrame):
        av = a.sort_values(list(a.columns)).values
        bv = b.sort_values(list(b.columns)).values
    else:
        av = a.sort_values().values
        bv = b.sort_values().values
    tm.assert_numpy_array_equal(av, bv)


@pytest.mark.parametrize('how', ['inner', 'left', 'right', 'outer'])
@pytest.mark.parametrize('shuffle', ['disk', 'tasks'])
def test_hash_join(how, shuffle):
    A = pd.DataFrame({'x': [1, 2, 3, 4, 5, 6], 'y': [1, 1, 2, 2, 3, 4]})
    a = dd.repartition(A, [0, 4, 5])

    B = pd.DataFrame({'y': [1, 3, 4, 4, 5, 6], 'z': [6, 5, 4, 3, 2, 1]})
    b = dd.repartition(B, [0, 2, 5])

    c = hash_join(a, 'y', b, 'y', how)

    result = c.compute()
    expected = pd.merge(A, B, how, 'y')
    list_eq(result, expected)

    # Different columns and npartitions
    c = hash_join(a, 'x', b, 'z', 'outer', npartitions=3, shuffle=shuffle)
    assert c.npartitions == 3

    result = c.compute()
    expected = pd.merge(A, B, 'outer', None, 'x', 'z')

    list_eq(result, expected)

    assert (hash_join(a, 'y', b, 'y', 'inner', shuffle=shuffle)._name ==
            hash_join(a, 'y', b, 'y', 'inner', shuffle=shuffle)._name)
    assert (hash_join(a, 'y', b, 'y', 'inner', shuffle=shuffle)._name !=
            hash_join(a, 'y', b, 'y', 'outer', shuffle=shuffle)._name)


@pytest.mark.parametrize('join', ['inner', 'outer'])
def test_indexed_concat(join):
    A = pd.DataFrame({'x': [1, 2, 3, 4, 6, 7], 'y': list('abcdef')},
                     index=[1, 2, 3, 4, 6, 7])
    a = dd.repartition(A, [1, 4, 7])

    B = pd.DataFrame({'x': [10, 20, 40, 50, 60, 80]},
                     index=[1, 2, 4, 5, 6, 8])
    b = dd.repartition(B, [1, 2, 5, 8])

    result = concat_indexed_dataframes([a, b], join=join)
    expected = pd.concat([A, B], axis=0, join=join)
    assert_eq(result, expected)

    assert (sorted(concat_indexed_dataframes([a, b], join=join).dask) ==
            sorted(concat_indexed_dataframes([a, b], join=join).dask))
    assert (sorted(concat_indexed_dataframes([a, b], join='inner').dask) !=
            sorted(concat_indexed_dataframes([a, b], join='outer').dask))


@pytest.mark.parametrize('join', ['inner', 'outer'])
def test_concat(join):
    pdf1 = pd.DataFrame({'x': [1, 2, 3, 4, 6, 7],
                         'y': list('abcdef')},
                        index=[1, 2, 3, 4, 6, 7])
    ddf1 = dd.from_pandas(pdf1, 2)
    pdf2 = pd.DataFrame({'x': [1, 2, 3, 4, 6, 7],
                         'y': list('abcdef')},
                        index=[8, 9, 10, 11, 12, 13])
    ddf2 = dd.from_pandas(pdf2, 2)

    # different columns
    pdf3 = pd.DataFrame({'x': [1, 2, 3, 4, 6, 7],
                         'z': list('abcdef')},
                        index=[8, 9, 10, 11, 12, 13])
    ddf3 = dd.from_pandas(pdf3, 2)

    for (dd1, dd2, pd1, pd2) in [(ddf1, ddf2, pdf1, pdf2),
                                 (ddf1, ddf3, pdf1, pdf3)]:
        result = dd.concat([dd1, dd2], join=join)
        expected = pd.concat([pd1, pd2], join=join)
        assert_eq(result, expected)

    # test outer only, inner has a problem on pandas side
    for (dd1, dd2, pd1, pd2) in [(ddf1, ddf2, pdf1, pdf2),
                                 (ddf1, ddf3, pdf1, pdf3),
                                 (ddf1.x, ddf2.x, pdf1.x, pdf2.x),
                                 (ddf1.x, ddf3.z, pdf1.x, pdf3.z),
                                 (ddf1.x, ddf2.x, pdf1.x, pdf2.x),
                                 (ddf1.x, ddf3.z, pdf1.x, pdf3.z)]:
        result = dd.concat([dd1, dd2])
        expected = pd.concat([pd1, pd2])
        assert_eq(result, expected)


@pytest.mark.parametrize('how', ['inner', 'outer', 'left', 'right'])
@pytest.mark.parametrize('shuffle', ['disk', 'tasks'])
def test_merge(how, shuffle):
    A = pd.DataFrame({'x': [1, 2, 3, 4, 5, 6], 'y': [1, 1, 2, 2, 3, 4]})
    a = dd.repartition(A, [0, 4, 5])

    B = pd.DataFrame({'y': [1, 3, 4, 4, 5, 6], 'z': [6, 5, 4, 3, 2, 1]})
    b = dd.repartition(B, [0, 2, 5])

    assert_eq(dd.merge(a, b, left_index=True, right_index=True,
                       how=how, shuffle=shuffle),
              pd.merge(A, B, left_index=True, right_index=True, how=how))

    result = dd.merge(a, b, on='y', how=how)
    list_eq(result, pd.merge(A, B, on='y', how=how))
    assert all(d is None for d in result.divisions)

    list_eq(dd.merge(a, b, left_on='x', right_on='z', how=how, shuffle=shuffle),
            pd.merge(A, B, left_on='x', right_on='z', how=how))
    list_eq(dd.merge(a, b, left_on='x', right_on='z', how=how,
                     suffixes=('1', '2'), shuffle=shuffle),
            pd.merge(A, B, left_on='x', right_on='z', how=how,
                     suffixes=('1', '2')))

    list_eq(dd.merge(a, b, how=how, shuffle=shuffle), pd.merge(A, B, how=how))
    list_eq(dd.merge(a, B, how=how, shuffle=shuffle), pd.merge(A, B, how=how))
    list_eq(dd.merge(A, b, how=how, shuffle=shuffle), pd.merge(A, B, how=how))
    list_eq(dd.merge(A, B, how=how, shuffle=shuffle), pd.merge(A, B, how=how))

    list_eq(dd.merge(a, b, left_index=True, right_index=True, how=how,
                     shuffle=shuffle),
            pd.merge(A, B, left_index=True, right_index=True, how=how))
    list_eq(dd.merge(a, b, left_index=True, right_index=True, how=how,
                     suffixes=('1', '2'), shuffle=shuffle),
            pd.merge(A, B, left_index=True, right_index=True, how=how,
                     suffixes=('1', '2')))

    list_eq(dd.merge(a, b, left_on='x', right_index=True, how=how,
                     shuffle=shuffle),
            pd.merge(A, B, left_on='x', right_index=True, how=how))
    list_eq(dd.merge(a, b, left_on='x', right_index=True, how=how,
                     suffixes=('1', '2'), shuffle=shuffle),
            pd.merge(A, B, left_on='x', right_index=True, how=how,
                     suffixes=('1', '2')))

    # pandas result looks buggy
    # list_eq(dd.merge(a, B, left_index=True, right_on='y'),
    #         pd.merge(A, B, left_index=True, right_on='y'))


def test_merge_tasks_passes_through():
    a = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6, 7],
                      'b': [7, 6, 5, 4, 3, 2, 1]})
    b = pd.DataFrame({'c': [1, 2, 3, 4, 5, 6, 7],
                      'd': [7, 6, 5, 4, 3, 2, 1]})

    aa = dd.from_pandas(a, npartitions=3)
    bb = dd.from_pandas(b, npartitions=2)

    cc = aa.merge(bb, left_on='a', right_on='d', shuffle='tasks')

    assert not any('partd' in k[0] for k in cc.dask)


@pytest.mark.parametrize('shuffle', ['disk', 'tasks'])
@pytest.mark.parametrize('how', ['inner', 'outer', 'left', 'right'])
def test_merge_by_index_patterns(how, shuffle):

    pdf1l = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6, 7],
                          'b': [7, 6, 5, 4, 3, 2, 1]})
    pdf1r = pd.DataFrame({'c': [1, 2, 3, 4, 5, 6, 7],
                          'd': [7, 6, 5, 4, 3, 2, 1]})

    pdf2l = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6, 7],
                          'b': [7, 6, 5, 4, 3, 2, 1]},
                         index=list('abcdefg'))
    pdf2r = pd.DataFrame({'c': [7, 6, 5, 4, 3, 2, 1],
                          'd': [7, 6, 5, 4, 3, 2, 1]},
                         index=list('abcdefg'))

    pdf3l = pdf2l
    pdf3r = pd.DataFrame({'c': [6, 7, 8, 9],
                          'd': [5, 4, 3, 2]},
                         index=list('abdg'))

    pdf4l = pdf2l
    pdf4r = pd.DataFrame({'c': [9, 10, 11, 12],
                          'd': [5, 4, 3, 2]},
                         index=list('abdg'))

    # completely different index
    pdf5l = pd.DataFrame({'a': [1, 1, 2, 2, 3, 3, 4],
                          'b': [7, 6, 5, 4, 3, 2, 1]},
                         index=list('lmnopqr'))
    pdf5r = pd.DataFrame({'c': [1, 1, 1, 1],
                          'd': [5, 4, 3, 2]},
                         index=list('abcd'))

    pdf6l = pd.DataFrame({'a': [1, 1, 2, 2, 3, 3, 4],
                          'b': [7, 6, 5, 4, 3, 2, 1]},
                         index=list('cdefghi'))
    pdf6r = pd.DataFrame({'c': [1, 2, 1, 2],
                          'd': [5, 4, 3, 2]},
                         index=list('abcd'))

    pdf7l = pd.DataFrame({'a': [1, 1, 2, 2, 3, 3, 4],
                          'b': [7, 6, 5, 4, 3, 2, 1]},
                         index=list('abcdefg'))
    pdf7r = pd.DataFrame({'c': [5, 6, 7, 8],
                          'd': [5, 4, 3, 2]},
                         index=list('fghi'))

    for pdl, pdr in [(pdf1l, pdf1r), (pdf2l, pdf2r), (pdf3l, pdf3r),
                     (pdf4l, pdf4r), (pdf5l, pdf5r), (pdf6l, pdf6r),
                     (pdf7l, pdf7r)]:

        for lpart, rpart in [(2, 2),    # same partition
                             (3, 2),    # left npartition > right npartition
                             (2, 3)]:   # left npartition < right npartition

            ddl = dd.from_pandas(pdl, lpart)
            ddr = dd.from_pandas(pdr, rpart)

            assert_eq(dd.merge(ddl, ddr, how=how, left_index=True,
                               right_index=True, shuffle=shuffle),
                      pd.merge(pdl, pdr, how=how, left_index=True,
                               right_index=True))
            assert_eq(dd.merge(ddr, ddl, how=how, left_index=True,
                               right_index=True, shuffle=shuffle),
                      pd.merge(pdr, pdl, how=how, left_index=True,
                               right_index=True))

            assert_eq(dd.merge(ddl, ddr, how=how, left_index=True,
                               right_index=True, shuffle=shuffle,
                               indicator=True),
                      pd.merge(pdl, pdr, how=how, left_index=True,
                               right_index=True, indicator=True))
            assert_eq(dd.merge(ddr, ddl, how=how, left_index=True,
                               right_index=True, shuffle=shuffle,
                               indicator=True),
                      pd.merge(pdr, pdl, how=how, left_index=True,
                               right_index=True, indicator=True))

            assert_eq(ddr.merge(ddl, how=how, left_index=True,
                                right_index=True, shuffle=shuffle),
                      pdr.merge(pdl, how=how, left_index=True,
                                right_index=True))
            assert_eq(ddl.merge(ddr, how=how, left_index=True,
                                right_index=True, shuffle=shuffle),
                      pdl.merge(pdr, how=how, left_index=True,
                                right_index=True))

            # hash join
            list_eq(dd.merge(ddl, ddr, how=how, left_on='a', right_on='c',
                             shuffle=shuffle),
                    pd.merge(pdl, pdr, how=how, left_on='a', right_on='c'))
            list_eq(dd.merge(ddl, ddr, how=how, left_on='b', right_on='d',
                             shuffle=shuffle),
                    pd.merge(pdl, pdr, how=how, left_on='b', right_on='d'))

            list_eq(dd.merge(ddr, ddl, how=how, left_on='c', right_on='a',
                             shuffle=shuffle, indicator=True),
                    pd.merge(pdr, pdl, how=how, left_on='c', right_on='a',
                             indicator=True))
            list_eq(dd.merge(ddr, ddl, how=how, left_on='d', right_on='b',
                             shuffle=shuffle, indicator=True),
                    pd.merge(pdr, pdl, how=how, left_on='d', right_on='b',
                             indicator=True))

            list_eq(dd.merge(ddr, ddl, how=how, left_on='c', right_on='a',
                             shuffle=shuffle),
                    pd.merge(pdr, pdl, how=how, left_on='c', right_on='a'))
            list_eq(dd.merge(ddr, ddl, how=how, left_on='d', right_on='b',
                             shuffle=shuffle),
                    pd.merge(pdr, pdl, how=how, left_on='d', right_on='b'))

            list_eq(ddl.merge(ddr, how=how, left_on='a', right_on='c',
                              shuffle=shuffle),
                    pdl.merge(pdr, how=how, left_on='a', right_on='c'))
            list_eq(ddl.merge(ddr, how=how, left_on='b', right_on='d',
                              shuffle=shuffle),
                    pdl.merge(pdr, how=how, left_on='b', right_on='d'))

            list_eq(ddr.merge(ddl, how=how, left_on='c', right_on='a',
                              shuffle=shuffle),
                    pdr.merge(pdl, how=how, left_on='c', right_on='a'))
            list_eq(ddr.merge(ddl, how=how, left_on='d', right_on='b',
                              shuffle=shuffle),
                    pdr.merge(pdl, how=how, left_on='d', right_on='b'))


@pytest.mark.parametrize('how', ['inner', 'outer', 'left', 'right'])
@pytest.mark.parametrize('shuffle', ['disk', 'tasks'])
def test_join_by_index_patterns(how, shuffle):

    # Similar test cases as test_merge_by_index_patterns,
    # but columns / index for join have same dtype

    pdf1l = pd.DataFrame({'a': list('abcdefg'),
                          'b': [7, 6, 5, 4, 3, 2, 1]},
                         index=list('abcdefg'))
    pdf1r = pd.DataFrame({'c': list('abcdefg'),
                          'd': [7, 6, 5, 4, 3, 2, 1]},
                         index=list('abcdefg'))

    pdf2l = pdf1l
    pdf2r = pd.DataFrame({'c': list('gfedcba'),
                          'd': [7, 6, 5, 4, 3, 2, 1]},
                         index=list('abcdefg'))

    pdf3l = pdf1l
    pdf3r = pd.DataFrame({'c': list('abdg'),
                          'd': [5, 4, 3, 2]},
                         index=list('abdg'))

    pdf4l = pd.DataFrame({'a': list('abcabce'),
                          'b': [7, 6, 5, 4, 3, 2, 1]},
                         index=list('abcdefg'))
    pdf4r = pd.DataFrame({'c': list('abda'),
                          'd': [5, 4, 3, 2]},
                         index=list('abdg'))

    # completely different index
    pdf5l = pd.DataFrame({'a': list('lmnopqr'),
                          'b': [7, 6, 5, 4, 3, 2, 1]},
                         index=list('lmnopqr'))
    pdf5r = pd.DataFrame({'c': list('abcd'),
                          'd': [5, 4, 3, 2]},
                         index=list('abcd'))

    pdf6l = pd.DataFrame({'a': list('cdefghi'),
                          'b': [7, 6, 5, 4, 3, 2, 1]},
                         index=list('cdefghi'))
    pdf6r = pd.DataFrame({'c': list('abab'),
                          'd': [5, 4, 3, 2]},
                         index=list('abcd'))

    pdf7l = pd.DataFrame({'a': list('aabbccd'),
                          'b': [7, 6, 5, 4, 3, 2, 1]},
                         index=list('abcdefg'))
    pdf7r = pd.DataFrame({'c': list('aabb'),
                          'd': [5, 4, 3, 2]},
                         index=list('fghi'))

    for pdl, pdr in [(pdf1l, pdf1r), (pdf2l, pdf2r), (pdf3l, pdf3r),
                     (pdf4l, pdf4r), (pdf5l, pdf5r), (pdf6l, pdf6r),
                     (pdf7l, pdf7r)]:

        for lpart, rpart in [(2, 2), (3, 2), (2, 3)]:

            ddl = dd.from_pandas(pdl, lpart)
            ddr = dd.from_pandas(pdr, rpart)

            assert_eq(ddl.join(ddr, how=how, shuffle=shuffle),
                      pdl.join(pdr, how=how))
            assert_eq(ddr.join(ddl, how=how, shuffle=shuffle),
                      pdr.join(pdl, how=how))

            assert_eq(ddl.join(ddr, how=how, lsuffix='l', rsuffix='r',
                               shuffle=shuffle),
                      pdl.join(pdr, how=how, lsuffix='l', rsuffix='r'))
            assert_eq(ddr.join(ddl, how=how, lsuffix='l', rsuffix='r',
                               shuffle=shuffle),
                      pdr.join(pdl, how=how, lsuffix='l', rsuffix='r'))

            """
            # temporary disabled bacause pandas may incorrectly raise
            # IndexError for empty DataFrame
            # https://github.com/pydata/pandas/pull/10826

            list_assert_eq(ddl.join(ddr, how=how, on='a', lsuffix='l', rsuffix='r'),
                    pdl.join(pdr, how=how, on='a', lsuffix='l', rsuffix='r'))

            list_eq(ddr.join(ddl, how=how, on='c', lsuffix='l', rsuffix='r'),
                    pdr.join(pdl, how=how, on='c', lsuffix='l', rsuffix='r'))

            # merge with index and columns
            list_eq(ddl.merge(ddr, how=how, left_on='a', right_index=True),
                    pdl.merge(pdr, how=how, left_on='a', right_index=True))
            list_eq(ddr.merge(ddl, how=how, left_on='c', right_index=True),
                    pdr.merge(pdl, how=how, left_on='c', right_index=True))
            list_eq(ddl.merge(ddr, how=how, left_index=True, right_on='c'),
                    pdl.merge(pdr, how=how, left_index=True, right_on='c'))
            list_eq(ddr.merge(ddl, how=how, left_index=True, right_on='a'),
                    pdr.merge(pdl, how=how, left_index=True, right_on='a'))
            """


@pytest.mark.parametrize('how', ['inner', 'outer', 'left', 'right'])
@pytest.mark.parametrize('shuffle', ['disk', 'tasks'])
def test_merge_by_multiple_columns(how, shuffle):
    # warnings here from pandas
    pdf1l = pd.DataFrame({'a': list('abcdefghij'),
                          'b': list('abcdefghij'),
                          'c': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]},
                         index=list('abcdefghij'))
    pdf1r = pd.DataFrame({'d': list('abcdefghij'),
                          'e': list('abcdefghij'),
                          'f': [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]},
                         index=list('abcdefghij'))

    pdf2l = pd.DataFrame({'a': list('abcdeabcde'),
                          'b': list('abcabcabca'),
                          'c': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]},
                         index=list('abcdefghij'))
    pdf2r = pd.DataFrame({'d': list('edcbaedcba'),
                          'e': list('aaabbbcccd'),
                          'f': [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]},
                         index=list('fghijklmno'))

    pdf3l = pd.DataFrame({'a': list('aaaaaaaaaa'),
                          'b': list('aaaaaaaaaa'),
                          'c': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]},
                         index=list('abcdefghij'))
    pdf3r = pd.DataFrame({'d': list('aaabbbccaa'),
                          'e': list('abbbbbbbbb'),
                          'f': [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]},
                         index=list('ABCDEFGHIJ'))

    for pdl, pdr in [(pdf1l, pdf1r), (pdf2l, pdf2r), (pdf3l, pdf3r)]:

        for lpart, rpart in [(2, 2), (3, 2), (2, 3)]:

            ddl = dd.from_pandas(pdl, lpart)
            ddr = dd.from_pandas(pdr, rpart)

            assert_eq(ddl.join(ddr, how=how, shuffle=shuffle),
                      pdl.join(pdr, how=how))
            assert_eq(ddr.join(ddl, how=how, shuffle=shuffle),
                      pdr.join(pdl, how=how))

            assert_eq(dd.merge(ddl, ddr, how=how, left_index=True,
                               right_index=True, shuffle=shuffle),
                      pd.merge(pdl, pdr, how=how, left_index=True,
                               right_index=True))
            assert_eq(dd.merge(ddr, ddl, how=how, left_index=True,
                               right_index=True, shuffle=shuffle),
                      pd.merge(pdr, pdl, how=how, left_index=True,
                               right_index=True))

            # hash join
            list_eq(dd.merge(ddl, ddr, how=how, left_on='a', right_on='d',
                             shuffle=shuffle),
                    pd.merge(pdl, pdr, how=how, left_on='a', right_on='d'))
            list_eq(dd.merge(ddl, ddr, how=how, left_on='b', right_on='e',
                             shuffle=shuffle),
                    pd.merge(pdl, pdr, how=how, left_on='b', right_on='e'))

            list_eq(dd.merge(ddr, ddl, how=how, left_on='d', right_on='a',
                             shuffle=shuffle),
                    pd.merge(pdr, pdl, how=how, left_on='d', right_on='a'))
            list_eq(dd.merge(ddr, ddl, how=how, left_on='e', right_on='b',
                             shuffle=shuffle),
                    pd.merge(pdr, pdl, how=how, left_on='e', right_on='b'))

            list_eq(dd.merge(ddl, ddr, how=how, left_on=['a', 'b'],
                             right_on=['d', 'e'], shuffle=shuffle),
                    pd.merge(pdl, pdr, how=how, left_on=['a', 'b'], right_on=['d', 'e']))


def test_melt():
    pdf = pd.DataFrame({'A': list('abcd') * 5,
                        'B': list('XY') * 10,
                        'C': np.random.randn(20)})
    ddf = dd.from_pandas(pdf, 4)

    list_eq(dd.melt(ddf),
            pd.melt(pdf))

    list_eq(dd.melt(ddf, id_vars='C'),
            pd.melt(pdf, id_vars='C'))
    list_eq(dd.melt(ddf, value_vars='C'),
            pd.melt(pdf, value_vars='C'))
    list_eq(dd.melt(ddf, value_vars=['A', 'C'], var_name='myvar'),
            pd.melt(pdf, value_vars=['A', 'C'], var_name='myvar'))
    list_eq(dd.melt(ddf, id_vars='B', value_vars=['A', 'C'], value_name='myval'),
            pd.melt(pdf, id_vars='B', value_vars=['A', 'C'], value_name='myval'))


def test_cheap_inner_merge_with_pandas_object():
    a = pd.DataFrame({'x': [1, 2, 3, 4, 5, 6], 'y': list('abdabd')},
                     index=[10, 20, 30, 40, 50, 60])
    da = dd.from_pandas(a, npartitions=3)

    b = pd.DataFrame({'x': [1, 2, 3, 4], 'z': list('abda')})

    dc = da.merge(b, on='x', how='inner')
    assert all('shuffle' not in k[0] for k in dc.dask)

    list_eq(da.merge(b, on='x', how='inner'),
            a.merge(b, on='x', how='inner'))


def test_cheap_single_partition_merge():
    a = pd.DataFrame({'x': [1, 2, 3, 4, 5, 6], 'y': list('abdabd')},
                     index=[10, 20, 30, 40, 50, 60])
    aa = dd.from_pandas(a, npartitions=3)

    b = pd.DataFrame({'x': [1, 2, 3, 4], 'z': list('abda')})
    bb = dd.from_pandas(b, npartitions=1, sort=False)

    cc = aa.merge(bb, on='x', how='inner')
    assert all('shuffle' not in k[0] for k in cc.dask)
    assert len(cc.dask) == len(aa.dask) * 2 + len(bb.dask)

    list_eq(aa.merge(bb, on='x', how='inner'),
            a.merge(b, on='x', how='inner'))


def test_cheap_single_partition_merge_divisions():
    a = pd.DataFrame({'x': [1, 2, 3, 4, 5, 6], 'y': list('abdabd')},
                     index=[10, 20, 30, 40, 50, 60])
    aa = dd.from_pandas(a, npartitions=3)

    b = pd.DataFrame({'x': [1, 2, 3, 4], 'z': list('abda')})
    bb = dd.from_pandas(b, npartitions=1, sort=False)

    actual = aa.merge(bb, on='x', how='inner')
    assert not actual.known_divisions
    assert_divisions(actual)

    actual = bb.merge(aa, on='x', how='inner')
    assert not actual.known_divisions
    assert_divisions(actual)


def test_cheap_single_partition_merge_on_index():
    a = pd.DataFrame({'x': [1, 2, 3, 4, 5, 6], 'y': list('abdabd')},
                     index=[10, 20, 30, 40, 50, 60])
    aa = dd.from_pandas(a, npartitions=3)

    b = pd.DataFrame({'x': [1, 2, 3, 4], 'z': list('abda')})
    bb = dd.from_pandas(b, npartitions=1, sort=False)

    actual = aa.merge(bb, left_index=True, right_on='x', how='inner')
    expected = a.merge(b, left_index=True, right_on='x', how='inner')

    assert actual.known_divisions
    assert_eq(actual, expected)

    actual = bb.merge(aa, right_index=True, left_on='x', how='inner')
    expected = b.merge(a, right_index=True, left_on='x', how='inner')

    assert actual.known_divisions
    assert_eq(actual, expected)


def test_merge_maintains_columns():
    lhs = pd.DataFrame({'A': [1, 2, 3],
                        'B': list('abc'),
                        'C': 'foo',
                        'D': 1.0},
                       columns=list('DCBA'))
    rhs = pd.DataFrame({'G': [4, 5],
                        'H': 6.0,
                        'I': 'bar',
                        'B': list('ab')},
                       columns=list('GHIB'))
    ddf = dd.from_pandas(lhs, npartitions=1)
    merged = dd.merge(ddf, rhs, on='B').compute()
    assert tuple(merged.columns) == ('D', 'C', 'B', 'A', 'G', 'H', 'I')


@pytest.mark.parametrize('shuffle', ['disk', 'tasks'])
def test_merge_index_without_divisions(shuffle):
    a = pd.DataFrame({'x': [1, 2, 3, 4, 5]}, index=[1, 2, 3, 4, 5])
    b = pd.DataFrame({'y': [1, 2, 3, 4, 5]}, index=[5, 4, 3, 2, 1])

    aa = dd.from_pandas(a, npartitions=3, sort=False)
    bb = dd.from_pandas(b, npartitions=2)

    result = aa.join(bb, how='inner', shuffle=shuffle)
    expected = a.join(b, how='inner')
    assert_eq(result, expected)


def test_half_indexed_dataframe_avoids_shuffle():
    a = pd.DataFrame({'x': np.random.randint(100, size=1000)})
    b = pd.DataFrame({'y': np.random.randint(100, size=100)},
                     index=np.random.randint(100, size=100))

    aa = dd.from_pandas(a, npartitions=100)
    bb = dd.from_pandas(b, npartitions=2)

    c = pd.merge(a, b, left_index=True, right_on='y')
    cc = dd.merge(aa, bb, left_index=True, right_on='y', shuffle='tasks')

    list_eq(c, cc)

    assert len(cc.dask) < 500


def test_errors_for_merge_on_frame_columns():
    a = pd.DataFrame({'x': [1, 2, 3, 4, 5]}, index=[1, 2, 3, 4, 5])
    b = pd.DataFrame({'y': [1, 2, 3, 4, 5]}, index=[5, 4, 3, 2, 1])

    aa = dd.from_pandas(a, npartitions=3, sort=False)
    bb = dd.from_pandas(b, npartitions=2)

    with pytest.raises(NotImplementedError):
        dd.merge(aa, bb, left_on='x', right_on=bb.y)

    with pytest.raises(NotImplementedError):
        dd.merge(aa, bb, left_on=aa.x, right_on=bb.y)


def test_concat_one_series():
    a = pd.Series([1, 2, 3, 4])
    aa = dd.from_pandas(a, npartitions=2, sort=False)

    c = dd.concat([aa], axis=0)
    assert isinstance(c, dd.Series)

    c = dd.concat([aa], axis=1)
    assert isinstance(c, dd.DataFrame)


def test_concat_unknown_divisions():
    a = pd.Series([1, 2, 3, 4])
    b = pd.Series([4, 3, 2, 1])
    aa = dd.from_pandas(a, npartitions=2, sort=False)
    bb = dd.from_pandas(b, npartitions=2, sort=False)

    assert not aa.known_divisions

    with pytest.warns(UserWarning):
        assert_eq(pd.concat([a, b], axis=1),
                  dd.concat([aa, bb], axis=1))

    cc = dd.from_pandas(b, npartitions=1, sort=False)
    with pytest.raises(ValueError):
        dd.concat([aa, cc], axis=1)


def test_concat_unknown_divisions_errors():
    a = pd.Series([1, 2, 3, 4, 5, 6])
    b = pd.Series([4, 3, 2, 1])
    aa = dd.from_pandas(a, npartitions=2, sort=False)
    bb = dd.from_pandas(b, npartitions=2, sort=False)

    with pytest.raises(ValueError):
        with pytest.warns(UserWarning):  # Concat with unknown divisions
            dd.concat([aa, bb], axis=1).compute()


def test_concat2():
    dsk = {('x', 0): pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]}),
           ('x', 1): pd.DataFrame({'a': [4, 5, 6], 'b': [3, 2, 1]}),
           ('x', 2): pd.DataFrame({'a': [7, 8, 9], 'b': [0, 0, 0]})}
    meta = make_meta({'a': 'i8', 'b': 'i8'})
    a = dd.DataFrame(dsk, 'x', meta, [None, None])
    dsk = {('y', 0): pd.DataFrame({'a': [10, 20, 30], 'b': [40, 50, 60]}),
           ('y', 1): pd.DataFrame({'a': [40, 50, 60], 'b': [30, 20, 10]}),
           ('y', 2): pd.DataFrame({'a': [70, 80, 90], 'b': [0, 0, 0]})}
    b = dd.DataFrame(dsk, 'y', meta, [None, None])

    dsk = {('y', 0): pd.DataFrame({'b': [10, 20, 30], 'c': [40, 50, 60]}),
           ('y', 1): pd.DataFrame({'b': [40, 50, 60], 'c': [30, 20, 10]})}
    meta = make_meta({'b': 'i8', 'c': 'i8'})
    c = dd.DataFrame(dsk, 'y', meta, [None, None])

    dsk = {('y', 0): pd.DataFrame({'b': [10, 20, 30], 'c': [40, 50, 60],
                                   'd': [70, 80, 90]}),
           ('y', 1): pd.DataFrame({'b': [40, 50, 60], 'c': [30, 20, 10],
                                   'd': [90, 80, 70]},
                                  index=[3, 4, 5])}
    meta = make_meta({'b': 'i8', 'c': 'i8', 'd': 'i8'},
                     index=pd.Index([], 'i8'))
    d = dd.DataFrame(dsk, 'y', meta, [0, 3, 5])

    cases = [[a, b], [a, c], [a, d]]
    assert dd.concat([a]) is a
    for case in cases:
        result = dd.concat(case)
        pdcase = [_c.compute() for _c in case]

        assert result.npartitions == case[0].npartitions + case[1].npartitions
        assert result.divisions == (None, ) * (result.npartitions + 1)
        assert_eq(pd.concat(pdcase), result)
        assert set(result.dask) == set(dd.concat(case).dask)

        result = dd.concat(case, join='inner')
        assert result.npartitions == case[0].npartitions + case[1].npartitions
        assert result.divisions == (None, ) * (result.npartitions + 1)
        assert_eq(pd.concat(pdcase, join='inner'), result)
        assert set(result.dask) == set(dd.concat(case, join='inner').dask)


def test_concat3():
    pdf1 = pd.DataFrame(np.random.randn(6, 5),
                        columns=list('ABCDE'), index=list('abcdef'))
    pdf2 = pd.DataFrame(np.random.randn(6, 5),
                        columns=list('ABCFG'), index=list('ghijkl'))
    pdf3 = pd.DataFrame(np.random.randn(6, 5),
                        columns=list('ABCHI'), index=list('mnopqr'))
    ddf1 = dd.from_pandas(pdf1, 2)
    ddf2 = dd.from_pandas(pdf2, 3)
    ddf3 = dd.from_pandas(pdf3, 2)

    result = dd.concat([ddf1, ddf2])
    assert result.divisions == ddf1.divisions[:-1] + ddf2.divisions
    assert result.npartitions == ddf1.npartitions + ddf2.npartitions
    assert_eq(result, pd.concat([pdf1, pdf2]))

    assert_eq(dd.concat([ddf1, ddf2], interleave_partitions=True),
              pd.concat([pdf1, pdf2]))

    result = dd.concat([ddf1, ddf2, ddf3])
    assert result.divisions == (ddf1.divisions[:-1] + ddf2.divisions[:-1] +
                                ddf3.divisions)
    assert result.npartitions == (ddf1.npartitions + ddf2.npartitions +
                                  ddf3.npartitions)
    assert_eq(result, pd.concat([pdf1, pdf2, pdf3]))

    assert_eq(dd.concat([ddf1, ddf2, ddf3], interleave_partitions=True),
              pd.concat([pdf1, pdf2, pdf3]))


def test_concat4_interleave_partitions():
    pdf1 = pd.DataFrame(np.random.randn(10, 5),
                        columns=list('ABCDE'), index=list('abcdefghij'))
    pdf2 = pd.DataFrame(np.random.randn(13, 5),
                        columns=list('ABCDE'), index=list('fghijklmnopqr'))
    pdf3 = pd.DataFrame(np.random.randn(13, 6),
                        columns=list('CDEXYZ'), index=list('fghijklmnopqr'))

    ddf1 = dd.from_pandas(pdf1, 2)
    ddf2 = dd.from_pandas(pdf2, 3)
    ddf3 = dd.from_pandas(pdf3, 2)

    msg = ('All inputs have known divisions which cannot be '
           'concatenated in order. Specify '
           'interleave_partitions=True to ignore order')

    cases = [[ddf1, ddf1], [ddf1, ddf2], [ddf1, ddf3], [ddf2, ddf1],
             [ddf2, ddf3], [ddf3, ddf1], [ddf3, ddf2]]
    for case in cases:
        pdcase = [c.compute() for c in case]

        with pytest.raises(ValueError) as err:
            dd.concat(case)
        assert msg in str(err.value)

        assert_eq(dd.concat(case, interleave_partitions=True),
                  pd.concat(pdcase))
        assert_eq(dd.concat(case, join='inner', interleave_partitions=True),
                  pd.concat(pdcase, join='inner'))

    msg = "'join' must be 'inner' or 'outer'"
    with pytest.raises(ValueError) as err:
        dd.concat([ddf1, ddf1], join='invalid', interleave_partitions=True)
    assert msg in str(err.value)


def test_concat5():
    pdf1 = pd.DataFrame(np.random.randn(7, 5),
                        columns=list('ABCDE'), index=list('abcdefg'))
    pdf2 = pd.DataFrame(np.random.randn(7, 6),
                        columns=list('FGHIJK'), index=list('abcdefg'))
    pdf3 = pd.DataFrame(np.random.randn(7, 6),
                        columns=list('FGHIJK'), index=list('cdefghi'))
    pdf4 = pd.DataFrame(np.random.randn(7, 5),
                        columns=list('FGHAB'), index=list('cdefghi'))
    pdf5 = pd.DataFrame(np.random.randn(7, 5),
                        columns=list('FGHAB'), index=list('fklmnop'))

    ddf1 = dd.from_pandas(pdf1, 2)
    ddf2 = dd.from_pandas(pdf2, 3)
    ddf3 = dd.from_pandas(pdf3, 2)
    ddf4 = dd.from_pandas(pdf4, 2)
    ddf5 = dd.from_pandas(pdf5, 3)

    cases = [[ddf1, ddf2], [ddf1, ddf3], [ddf1, ddf4], [ddf1, ddf5],
             [ddf3, ddf4], [ddf3, ddf5], [ddf5, ddf1, ddf4], [ddf5, ddf3],
             [ddf1.A, ddf4.A], [ddf2.F, ddf3.F], [ddf4.A, ddf5.A],
             [ddf1.A, ddf4.F], [ddf2.F, ddf3.H], [ddf4.A, ddf5.B],
             [ddf1, ddf4.A], [ddf3.F, ddf2], [ddf5, ddf1.A, ddf2]]

    for case in cases:
        pdcase = [c.compute() for c in case]

        with pytest.warns(None):
            # some cases will raise warning directly from pandas
            assert_eq(dd.concat(case, interleave_partitions=True),
                      pd.concat(pdcase))

        assert_eq(dd.concat(case, join='inner', interleave_partitions=True),
                  pd.concat(pdcase, join='inner'))

        assert_eq(dd.concat(case, axis=1), pd.concat(pdcase, axis=1))

        assert_eq(dd.concat(case, axis=1, join='inner'),
                  pd.concat(pdcase, axis=1, join='inner'))

    # Dask + pandas
    cases = [[ddf1, pdf2], [ddf1, pdf3], [pdf1, ddf4],
             [pdf1.A, ddf4.A], [ddf2.F, pdf3.F],
             [ddf1, pdf4.A], [ddf3.F, pdf2], [ddf2, pdf1, ddf3.F]]

    for case in cases:
        pdcase = [c.compute() if isinstance(c, _Frame) else c for c in case]

        assert_eq(dd.concat(case, interleave_partitions=True),
                  pd.concat(pdcase))

        assert_eq(dd.concat(case, join='inner', interleave_partitions=True),
                  pd.concat(pdcase, join='inner'))

        assert_eq(dd.concat(case, axis=1), pd.concat(pdcase, axis=1))

        assert_eq(dd.concat(case, axis=1, join='inner'),
                  pd.concat(pdcase, axis=1, join='inner'))


@pytest.mark.parametrize('known, cat_index, divisions',
                         [(True, True, False), (True, False, True),
                          (True, False, False), (False, True, False),
                          (False, False, True), (False, False, False)])
def test_concat_categorical(known, cat_index, divisions):
    frames = [pd.DataFrame({'w': list('xxxxx'),
                            'x': np.arange(5),
                            'y': list('abcbc'),
                            'z': np.arange(5, dtype='f8')}),
              pd.DataFrame({'w': list('yyyyy'),
                            'x': np.arange(5, 10),
                            'y': list('abbba'),
                            'z': np.arange(5, 10, dtype='f8')}),
              pd.DataFrame({'w': list('zzzzz'),
                            'x': np.arange(10, 15),
                            'y': list('bcbcc'),
                            'z': np.arange(10, 15, dtype='f8')})]
    for df in frames:
        df.w = df.w.astype('category')
        df.y = df.y.astype('category')

    if cat_index:
        frames = [df.set_index(df.y) for df in frames]

    dframes = [dd.from_pandas(p, npartitions=2, sort=divisions) for p in frames]

    if not known:
        dframes[0]._meta = clear_known_categories(dframes[0]._meta, ['y'],
                                                  index=True)

    def check_and_return(ddfs, dfs, join):
        sol = concat(dfs, join=join)
        res = dd.concat(ddfs, join=join, interleave_partitions=divisions)
        assert_eq(res, sol)
        if known:
            parts = compute_as_if_collection(dd.DataFrame, res.dask,
                                             res.__dask_keys__())
            for p in [i.iloc[:0] for i in parts]:
                res._meta == p  # will error if schemas don't align
        assert not cat_index or has_known_categories(res.index) == known
        return res

    for join in ['inner', 'outer']:
        # Frame
        res = check_and_return(dframes, frames, join)
        assert has_known_categories(res.w)
        assert has_known_categories(res.y) == known

        # Series
        res = check_and_return([i.y for i in dframes],
                               [i.y for i in frames], join)
        assert has_known_categories(res) == known

        # Non-cat series with cat index
        if cat_index:
            res = check_and_return([i.x for i in dframes],
                                   [i.x for i in frames], join)

        # Partition missing columns
        res = check_and_return([dframes[0][['x', 'y']]] + dframes[1:],
                               [frames[0][['x', 'y']]] + frames[1:], join)
        assert not hasattr(res, 'w') or has_known_categories(res.w)
        assert has_known_categories(res.y) == known


def test_concat_datetimeindex():
    # https://github.com/dask/dask/issues/2932
    b2 = pd.DataFrame({'x': ['a']},
                      index=pd.DatetimeIndex(['2015-03-24 00:00:16'],
                                             dtype='datetime64[ns]'))
    b3 = pd.DataFrame({'x': ['c']},
                      index=pd.DatetimeIndex(['2015-03-29 00:00:44'],
                                             dtype='datetime64[ns]'))

    b2['x'] = b2.x.astype('category').cat.set_categories(['a', 'c'])
    b3['x'] = b3.x.astype('category').cat.set_categories(['a', 'c'])

    db2 = dd.from_pandas(b2, 1)
    db3 = dd.from_pandas(b3, 1)

    result = concat([b2.iloc[:0], b3.iloc[:0]])
    assert result.index.dtype == '<M8[ns]'

    result = dd.concat([db2, db3])
    expected = pd.concat([b2, b3])
    assert_eq(result, expected)


def test_append():
    df = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6],
                       'b': [1, 2, 3, 4, 5, 6]})
    df2 = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6],
                        'b': [1, 2, 3, 4, 5, 6]},
                       index=[6, 7, 8, 9, 10, 11])
    df3 = pd.DataFrame({'b': [1, 2, 3, 4, 5, 6],
                        'c': [1, 2, 3, 4, 5, 6]},
                       index=[6, 7, 8, 9, 10, 11])

    ddf = dd.from_pandas(df, 2)
    ddf2 = dd.from_pandas(df2, 2)
    ddf3 = dd.from_pandas(df3, 2)

    s = pd.Series([7, 8], name=6, index=['a', 'b'])
    assert_eq(ddf.append(s), df.append(s))

    assert_eq(ddf.append(ddf2), df.append(df2))
    assert_eq(ddf.a.append(ddf2.a), df.a.append(df2.a))
    # different columns
    assert_eq(ddf.append(ddf3), df.append(df3))
    assert_eq(ddf.a.append(ddf3.b), df.a.append(df3.b))

    # dask + pandas
    assert_eq(ddf.append(df2), df.append(df2))
    assert_eq(ddf.a.append(df2.a), df.a.append(df2.a))

    assert_eq(ddf.append(df3), df.append(df3))
    assert_eq(ddf.a.append(df3.b), df.a.append(df3.b))

    df4 = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6],
                        'b': [1, 2, 3, 4, 5, 6]},
                       index=[4, 5, 6, 7, 8, 9])
    ddf4 = dd.from_pandas(df4, 2)
    with pytest.raises(ValueError):
        ddf.append(ddf4)


def test_append2():
    dsk = {('x', 0): pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]}),
           ('x', 1): pd.DataFrame({'a': [4, 5, 6], 'b': [3, 2, 1]}),
           ('x', 2): pd.DataFrame({'a': [7, 8, 9], 'b': [0, 0, 0]})}
    meta = make_meta({'a': 'i8', 'b': 'i8'})
    ddf1 = dd.DataFrame(dsk, 'x', meta, [None, None])

    dsk = {('y', 0): pd.DataFrame({'a': [10, 20, 30], 'b': [40, 50, 60]}),
           ('y', 1): pd.DataFrame({'a': [40, 50, 60], 'b': [30, 20, 10]}),
           ('y', 2): pd.DataFrame({'a': [70, 80, 90], 'b': [0, 0, 0]})}
    ddf2 = dd.DataFrame(dsk, 'y', meta, [None, None])

    dsk = {('y', 0): pd.DataFrame({'b': [10, 20, 30], 'c': [40, 50, 60]}),
           ('y', 1): pd.DataFrame({'b': [40, 50, 60], 'c': [30, 20, 10]})}
    meta = make_meta({'b': 'i8', 'c': 'i8'})
    ddf3 = dd.DataFrame(dsk, 'y', meta, [None, None])

    assert_eq(ddf1.append(ddf2), ddf1.compute().append(ddf2.compute()))
    assert_eq(ddf2.append(ddf1), ddf2.compute().append(ddf1.compute()))
    # Series + DataFrame
    with pytest.warns(None):
        # RuntimeWarning from pandas on comparing int and str
        assert_eq(ddf1.a.append(ddf2), ddf1.a.compute().append(ddf2.compute()))
        assert_eq(ddf2.a.append(ddf1), ddf2.a.compute().append(ddf1.compute()))

    # different columns
    assert_eq(ddf1.append(ddf3), ddf1.compute().append(ddf3.compute()))
    assert_eq(ddf3.append(ddf1), ddf3.compute().append(ddf1.compute()))
    # Series + DataFrame
    with pytest.warns(None):
        # RuntimeWarning from pandas on comparing int and str
        assert_eq(ddf1.a.append(ddf3), ddf1.a.compute().append(ddf3.compute()))
        assert_eq(ddf3.b.append(ddf1), ddf3.b.compute().append(ddf1.compute()))

    # Dask + pandas
    assert_eq(ddf1.append(ddf2.compute()), ddf1.compute().append(ddf2.compute()))
    assert_eq(ddf2.append(ddf1.compute()), ddf2.compute().append(ddf1.compute()))
    # Series + DataFrame
    with pytest.warns(None):
        # RuntimeWarning from pandas on comparing int and str
        assert_eq(ddf1.a.append(ddf2.compute()), ddf1.a.compute().append(ddf2.compute()))
        assert_eq(ddf2.a.append(ddf1.compute()), ddf2.a.compute().append(ddf1.compute()))

    # different columns
    assert_eq(ddf1.append(ddf3.compute()), ddf1.compute().append(ddf3.compute()))
    assert_eq(ddf3.append(ddf1.compute()), ddf3.compute().append(ddf1.compute()))
    # Series + DataFrame
    with pytest.warns(None):
        # RuntimeWarning from pandas on comparing int and str
        assert_eq(ddf1.a.append(ddf3.compute()), ddf1.a.compute().append(ddf3.compute()))
        assert_eq(ddf3.b.append(ddf1.compute()), ddf3.b.compute().append(ddf1.compute()))


def test_append_categorical():
    frames = [pd.DataFrame({'x': np.arange(5, 10),
                            'y': list('abbba'),
                            'z': np.arange(5, 10, dtype='f8')}),
              pd.DataFrame({'x': np.arange(10, 15),
                            'y': list('bcbcc'),
                            'z': np.arange(10, 15, dtype='f8')})]
    frames2 = []
    for df in frames:
        df.y = df.y.astype('category')
        df2 = df.copy()
        df2.y = df2.y.cat.set_categories(list('abc'))
        df.index = df.y
        frames2.append(df2.set_index(df2.y))

    df1, df2 = frames2

    for known in [True, False]:
        dframes = [dd.from_pandas(p, npartitions=2, sort=False) for p in frames]
        if not known:
            dframes[0]._meta = clear_known_categories(dframes[0]._meta,
                                                      ['y'], index=True)
        ddf1, ddf2 = dframes

        res = ddf1.append(ddf2)
        assert_eq(res, df1.append(df2))
        assert has_known_categories(res.index) == known
        assert has_known_categories(res.y) == known

        res = ddf1.y.append(ddf2.y)
        assert_eq(res, df1.y.append(df2.y))
        assert has_known_categories(res.index) == known
        assert has_known_categories(res) == known

        res = ddf1.index.append(ddf2.index)
        assert_eq(res, df1.index.append(df2.index))
        assert has_known_categories(res) == known


def test_singleton_divisions():
    df = pd.DataFrame({'x': [1, 1, 1]}, index=[1, 2, 3])
    ddf = dd.from_pandas(df, npartitions=2)
    ddf2 = ddf.set_index('x')

    joined = ddf2.join(ddf2, rsuffix='r')
    assert joined.divisions == (1, 1)
    joined.compute()
