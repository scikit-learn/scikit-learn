import pytest
import numpy as np
import pandas as pd

import pandas.util.testing as tm
from pandas.core.indexes.api import Index, MultiIndex
from pandas.compat import lzip, long


@pytest.fixture(params=[tm.makeUnicodeIndex(100),
                        tm.makeStringIndex(100),
                        tm.makeDateIndex(100),
                        tm.makePeriodIndex(100),
                        tm.makeTimedeltaIndex(100),
                        tm.makeIntIndex(100),
                        tm.makeUIntIndex(100),
                        tm.makeFloatIndex(100),
                        Index([True, False]),
                        tm.makeCategoricalIndex(100),
                        Index([]),
                        MultiIndex.from_tuples(lzip(
                            ['foo', 'bar', 'baz'], [1, 2, 3])),
                        Index([0, 0, 1, 1, 2, 2])],
                ids=lambda x: type(x).__name__)
def indices(request):
    return request.param


@pytest.fixture(params=[1, np.array(1, dtype=np.int64)])
def one(request):
    # zero-dim integer array behaves like an integer
    return request.param


zeros = [box([0] * 5, dtype=dtype)
         for box in [pd.Index, np.array]
         for dtype in [np.int64, np.uint64, np.float64]]
zeros.extend([np.array(0, dtype=dtype)
              for dtype in [np.int64, np.uint64, np.float64]])
zeros.extend([0, 0.0, long(0)])


@pytest.fixture(params=zeros)
def zero(request):
    # For testing division by (or of) zero for Index with length 5, this
    # gives several scalar-zeros and length-5 vector-zeros
    return request.param
