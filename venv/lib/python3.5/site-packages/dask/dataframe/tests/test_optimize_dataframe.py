import pytest
from operator import getitem
from toolz import merge

import dask
from dask.dataframe.io import dataframe_from_ctable
import dask.dataframe as dd
import pandas as pd

dsk = {('x', 0): pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]},
                              index=[0, 1, 3]),
       ('x', 1): pd.DataFrame({'a': [4, 5, 6], 'b': [3, 2, 1]},
                              index=[5, 6, 8]),
       ('x', 2): pd.DataFrame({'a': [7, 8, 9], 'b': [0, 0, 0]},
                              index=[9, 9, 9])}
dfs = list(dsk.values())


def test_column_optimizations_with_bcolz_and_rewrite():
    bcolz = pytest.importorskip('bcolz')

    bc = bcolz.ctable([[1, 2, 3], [10, 20, 30]], names=['a', 'b'])
    for cols in [None, 'abc', ['abc']]:
        dsk2 = merge(dict((('x', i),
                          (dataframe_from_ctable, bc, slice(0, 2), cols, {}))
                          for i in [1, 2, 3]),
                     dict((('y', i),
                          (getitem, ('x', i), ['a', 'b']))
                          for i in [1, 2, 3]))

        expected = dict((('y', i), (dataframe_from_ctable,
                                    bc, slice(0, 2), ['a', 'b'], {}))
                        for i in [1, 2, 3])

        result = dd.optimize(dsk2, [('y', i) for i in [1, 2, 3]])
        assert result == expected


def test_fuse_ave_width():
    df = pd.DataFrame({'x': range(10)})
    df = dd.from_pandas(df, npartitions=5)

    s = ((df.x + 1) + (df.x + 2))

    with dask.set_options(fuse_ave_width=4):
        a = s.__dask_optimize__(s.dask, s.__dask_keys__())

    b = s.__dask_optimize__(s.dask, s.__dask_keys__())

    assert len(a) < len(b)
    assert len(a) <= 15
