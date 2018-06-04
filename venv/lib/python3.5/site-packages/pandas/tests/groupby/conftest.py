import pytest
import numpy as np
from pandas import MultiIndex, DataFrame
from pandas.util import testing as tm


@pytest.fixture
def mframe():
    index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'], ['one', 'two',
                                                              'three']],
                       labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                               [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                       names=['first', 'second'])
    return DataFrame(np.random.randn(10, 3), index=index,
                     columns=['A', 'B', 'C'])


@pytest.fixture
def df():
    return DataFrame(
        {'A': ['foo', 'bar', 'foo', 'bar', 'foo', 'bar', 'foo', 'foo'],
         'B': ['one', 'one', 'two', 'three', 'two', 'two', 'one', 'three'],
         'C': np.random.randn(8),
         'D': np.random.randn(8)})


@pytest.fixture
def ts():
    return tm.makeTimeSeries()


@pytest.fixture
def seriesd():
    return tm.getSeriesData()


@pytest.fixture
def tsd():
    return tm.getTimeSeriesData()


@pytest.fixture
def frame(seriesd):
    return DataFrame(seriesd)


@pytest.fixture
def tsframe(tsd):
    return DataFrame(tsd)


@pytest.fixture
def df_mixed_floats():
    return DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                            'foo', 'bar', 'foo', 'foo'],
                      'B': ['one', 'one', 'two', 'three',
                            'two', 'two', 'one', 'three'],
                      'C': np.random.randn(8),
                      'D': np.array(
                          np.random.randn(8), dtype='float32')})


@pytest.fixture
def three_group():
    return DataFrame({'A': ['foo', 'foo', 'foo',
                            'foo', 'bar', 'bar',
                            'bar', 'bar',
                            'foo', 'foo', 'foo'],
                      'B': ['one', 'one', 'one',
                            'two', 'one', 'one', 'one', 'two',
                            'two', 'two', 'one'],
                      'C': ['dull', 'dull', 'shiny',
                            'dull', 'dull', 'shiny', 'shiny',
                            'dull', 'shiny', 'shiny', 'shiny'],
                      'D': np.random.randn(11),
                      'E': np.random.randn(11),
                      'F': np.random.randn(11)})
