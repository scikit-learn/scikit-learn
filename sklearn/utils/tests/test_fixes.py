# Authors: Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Justin Vincent
#          Lars Buitinck
# License: BSD 3 clause

import pickle

import numpy as np
import pytest

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_allclose

from sklearn.utils.fixes import divide
from sklearn.utils.fixes import MaskedArray
from sklearn.utils.fixes import nanmedian
from sklearn.utils.fixes import nanpercentile
from sklearn.utils.fixes import _joblib_parallel_args
from sklearn.utils.fixes import _object_dtype_isnan


def test_divide():
    assert_equal(divide(.6, 1), .600000000000)


def test_masked_array_obj_dtype_pickleable():
    marr = MaskedArray([1, None, 'a'], dtype=object)

    for mask in (True, False, [0, 1, 0]):
        marr.mask = mask
        marr_pickled = pickle.loads(pickle.dumps(marr))
        assert_array_equal(marr.data, marr_pickled.data)
        assert_array_equal(marr.mask, marr_pickled.mask)


@pytest.mark.parametrize(
    "axis, expected_median",
    [(None, 4.0),
     (0, np.array([1., 3.5, 3.5, 4., 7., np.nan])),
     (1, np.array([1., 6.]))]
)
def test_nanmedian(axis, expected_median):
    X = np.array([[1, 1, 1, 2, np.nan, np.nan],
                  [np.nan, 6, 6, 6, 7, np.nan]])
    median = nanmedian(X, axis=axis)
    if axis is None:
        assert median == pytest.approx(expected_median)
    else:
        assert_allclose(median, expected_median)


@pytest.mark.parametrize(
    "a, q, expected_percentile",
    [(np.array([1, 2, 3, np.nan]), [0, 50, 100], np.array([1., 2., 3.])),
     (np.array([1, 2, 3, np.nan]), 50, 2.),
     (np.array([np.nan, np.nan]), [0, 50], np.array([np.nan, np.nan]))]
)
def test_nanpercentile(a, q, expected_percentile):
    percentile = nanpercentile(a, q)
    assert_allclose(percentile, expected_percentile)


@pytest.mark.parametrize('joblib_version', ('0.11', '0.12.0'))
def test_joblib_parallel_args(monkeypatch, joblib_version):
    import sklearn.utils._joblib
    monkeypatch.setattr(sklearn.utils._joblib, '__version__', joblib_version)

    if joblib_version == '0.12.0':
        # arguments are simply passed through
        assert _joblib_parallel_args(prefer='threads') == {'prefer': 'threads'}
        assert _joblib_parallel_args(prefer='processes', require=None) == {
                    'prefer': 'processes', 'require': None}
        assert _joblib_parallel_args(non_existing=1) == {'non_existing': 1}
    elif joblib_version == '0.11':
        # arguments are mapped to the corresponding backend
        assert _joblib_parallel_args(prefer='threads') == {
                    'backend': 'threading'}
        assert _joblib_parallel_args(prefer='processes') == {
                    'backend': 'multiprocessing'}
        with pytest.raises(ValueError):
            _joblib_parallel_args(prefer='invalid')
        assert _joblib_parallel_args(
                prefer='processes', require='sharedmem') == {
                    'backend': 'threading'}
        with pytest.raises(ValueError):
            _joblib_parallel_args(require='invalid')
        with pytest.raises(NotImplementedError):
            _joblib_parallel_args(verbose=True)
    else:
        raise ValueError


@pytest.mark.parametrize("dtype, val", ([object, 1],
                                        [object, "a"],
                                        [float, 1]))
def test_object_dtype_isnan(dtype, val):
    X = np.array([[val, np.nan],
                  [np.nan, val]], dtype=dtype)

    expected_mask = np.array([[False, True],
                              [True, False]])

    mask = _object_dtype_isnan(X)

    assert_array_equal(mask, expected_mask)
