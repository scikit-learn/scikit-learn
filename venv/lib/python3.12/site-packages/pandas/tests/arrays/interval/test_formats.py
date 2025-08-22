from pandas.core.arrays import IntervalArray


def test_repr():
    # GH#25022
    arr = IntervalArray.from_tuples([(0, 1), (1, 2)])
    result = repr(arr)
    expected = (
        "<IntervalArray>\n[(0, 1], (1, 2]]\nLength: 2, dtype: interval[int64, right]"
    )
    assert result == expected
