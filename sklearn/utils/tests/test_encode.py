import numpy as np
import pytest
from numpy.testing import assert_array_equal

from sklearn.utils._encode import _unique
from sklearn.utils._encode import _encode
from sklearn.utils._encode import _encode_check_unknown


@pytest.mark.parametrize(
        "values, expected",
        [(np.array([2, 1, 3, 1, 3], dtype='int64'),
          np.array([1, 2, 3], dtype='int64')),
         (np.array(['b', 'a', 'c', 'a', 'c'], dtype=object),
          np.array(['a', 'b', 'c'], dtype=object)),
         (np.array(['b', 'a', 'c', 'a', 'c']),
          np.array(['a', 'b', 'c']))],
        ids=['int64', 'object', 'str'])
def test_encode_util(values, expected):
    uniques = _unique(values)
    assert_array_equal(uniques, expected)
    encoded = _encode(values, uniques=uniques)
    assert_array_equal(encoded, np.array([1, 0, 2, 0, 2]))


def test_encode_with_check_unknown():
    # test for the check_unknown parameter of _encode()
    uniques = np.array([1, 2, 3])
    values = np.array([1, 2, 3, 4])

    # Default is True, raise error
    with pytest.raises(ValueError,
                       match='y contains previously unseen labels'):
        _encode(values, uniques=uniques, check_unknown=True)

    # dont raise error if False
    _encode(values, uniques=uniques, check_unknown=False)

    # parameter is ignored for object dtype
    uniques = np.array(['a', 'b', 'c'], dtype=object)
    values = np.array(['a', 'b', 'c', 'd'], dtype=object)
    with pytest.raises(ValueError,
                       match='y contains previously unseen labels'):
        _encode(values, uniques=uniques, check_unknown=False)


@pytest.mark.parametrize("values, uniques, expected_diff", [
  (np.array([1, 2, 3, 4]), np.array([1, 2, 3]), [4]),
  (np.array(['a', 'b', 'c', 'd'], dtype=object),
   np.array(['a', 'b', 'c'], dtype=object),
   np.array(['d']))
])
def test_encode_check_unknown(values, uniques, expected_diff):
    diff = _encode_check_unknown(values, uniques)

    assert_array_equal(diff, expected_diff)

    diff, valid_mask = _encode_check_unknown(values, uniques, return_mask=True)

    assert_array_equal(diff, expected_diff)
    assert_array_equal(valid_mask, [True, True, True, False])
