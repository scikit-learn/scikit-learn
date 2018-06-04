from __future__ import print_function, division, absolute_import

import sys

from .utils import Dispatch

try:  # PyPy does not support sys.getsizeof
    sys.getsizeof(1)
    getsizeof = sys.getsizeof
except (AttributeError, TypeError):  # Monkey patch
    getsizeof = lambda x: 100


sizeof = Dispatch(name='sizeof')


@sizeof.register(object)
def sizeof_default(o):
    return getsizeof(o)


@sizeof.register(list)
@sizeof.register(tuple)
@sizeof.register(set)
@sizeof.register(frozenset)
def sizeof_python_collection(seq):
    return getsizeof(seq) + sum(map(sizeof, seq))


@sizeof.register_lazy("numpy")
def register_numpy():
    import numpy as np

    @sizeof.register(np.ndarray)
    def sizeof_numpy_ndarray(x):
        return int(x.nbytes)


@sizeof.register_lazy("pandas")
def register_pandas():
    import pandas as pd

    @sizeof.register(pd.DataFrame)
    def sizeof_pandas_dataframe(df):
        p = int(df.memory_usage(index=True).sum())
        obj = int((df.dtypes == object).sum() * len(df) * 100)
        if df.index.dtype == object:
            obj += len(df) * 100
        return int(p + obj) + 1000

    @sizeof.register(pd.Series)
    def sizeof_pandas_series(s):
        p = int(s.memory_usage(index=True))
        if s.dtype == object:
            p += len(s) * 100
        if s.index.dtype == object:
            p += len(s) * 100
        return int(p) + 1000

    @sizeof.register(pd.Index)
    def sizeof_pandas_index(i):
        p = int(i.memory_usage())
        obj = len(i) * 100 if i.dtype == object else 0
        return int(p + obj) + 1000


@sizeof.register_lazy("scipy")
def register_spmatrix():
    from scipy import sparse

    @sizeof.register(sparse.dok_matrix)
    def sizeof_spmatrix_dok(s):
        return s.__sizeof__()

    @sizeof.register(sparse.spmatrix)
    def sizeof_spmatrix(s):
        return sum(
            sizeof(v) for v in s.__dict__.values()
        )
