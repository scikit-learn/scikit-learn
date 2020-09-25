"""Test the covtype loader, if the data is available,
or if specifically requested via environment variable
(e.g. for travis cron job)."""
from functools import partial
import pytest
from sklearn.datasets.tests.test_common import check_return_X_y
from sklearn.datasets._covtype import fetch_covtype
from unittest.mock import Mock, patch
import numpy as np


def test_fetch(fetch_covtype_fxt):
    data1 = fetch_covtype_fxt(shuffle=True, random_state=42)
    data2 = fetch_covtype_fxt(shuffle=True, random_state=37)

    X1, X2 = data1['data'], data2['data']
    assert (581012, 54) == X1.shape
    assert X1.shape == X2.shape

    assert X1.sum() == X2.sum()

    y1, y2 = data1['target'], data2['target']
    assert (X1.shape[0],) == y1.shape
    assert (X1.shape[0],) == y2.shape

    # test return_X_y option
    fetch_func = partial(fetch_covtype_fxt)
    check_return_X_y(data1, fetch_func)


def test_fetch_asframe(fetch_covtype_fxt):
    pytest.importorskip("pandas")

    bunch = fetch_covtype_fxt(as_frame=True)
    assert hasattr(bunch, 'frame')
    frame = bunch.frame
    assert frame.shape == (581012, 55)
    assert bunch.data.shape == (581012, 54)
    assert bunch.target.shape == (581012,)

    column_names = set(frame.columns)

    # enumerated names are added correctly
    assert set(f"Wilderness_Area_{i}" for i in range(4)) < column_names
    assert set(f"Soil_Type_{i}" for i in range(40)) < column_names


def test_pandas_dependency_message(fetch_covtype_fxt,
                                   hide_available_pandas):
    expected_msg = ('fetch_covtype with as_frame=True'
                    ' requires pandas')
    with pytest.raises(ImportError, match=expected_msg):
        fetch_covtype_fxt(as_frame=True)

@patch('sklearn.datasets._covtype.np.genfromtxt')
@patch('sklearn.datasets._covtype.exists', return_value=False)
@patch('sklearn.datasets._covtype.makedirs', Mock())
@patch('sklearn.datasets._covtype._fetch_remote', Mock())
@patch('sklearn.datasets._covtype.GzipFile', Mock())
@patch('sklearn.datasets._covtype.remove', Mock())
def test_fetch_with_download(exists, genfromtext):
    """
        This test covers the fetch of the dataset from the network
        and the operations that are mandatory when you dont have the dataset locally
        - create the folder
        - gzip it
        - cast the numpy result
        - remove the file
        - dump in joblib
    """
    genfromtext.return_value = np.array([
        [2.596e+03, 5.100e+01, 3.000e+00, 2.580e+02, 0.000e+00, 5.100e+02,
         2.210e+02, 2.320e+02, 1.480e+02, 6.279e+03, 1.000e+00, 0.000e+00,
         0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
         0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
         0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
         0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
         0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
         1.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
         0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
         5.000e+00],

        [2.590e+03, 5.600e+01, 2.000e+00, 2.120e+02, -6.000e+00,
         3.900e+02, 2.200e+02, 2.350e+02, 1.510e+02, 6.225e+03,
         1.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
         0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
         0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
         0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
         0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
         0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
         0.000e+00, 0.000e+00, 1.000e+00, 0.000e+00, 0.000e+00,
         0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
         0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 5.000e+00]
        ])
    fetch_covtype(download_if_missing=True)
