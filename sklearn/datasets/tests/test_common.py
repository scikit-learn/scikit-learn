"""Test loaders for common functionality."""
import inspect
import os

import pytest
import numpy as np

import sklearn.datasets


KNOWN_FAILURE = {
    "return_X_y": {
        "fetch_20newsgroups": pytest.mark.xfail(
            reason="X is a list and does not have a shape argument"
        ),
        "fetch_openml": pytest.mark.xfail(
            reason="fetch_opeml requires a dataset name or id"
        ),
    },
    "as_frame": {
        "fetch_openml": pytest.mark.xfail(
            reason="fetch_opeml requires a dataset name or id"
        ),
    }
}


def check_pandas_dependency_message(fetch_func):
    try:
        import pandas  # noqa
        pytest.skip("This test requires pandas to not be installed")
    except ImportError:
        # Check that pandas is imported lazily and that an informative error
        # message is raised when pandas is missing:
        name = fetch_func.__name__
        expected_msg = f'{name} with as_frame=True requires pandas'
        with pytest.raises(ImportError, match=expected_msg):
            fetch_func(as_frame=True)


def check_return_X_y(bunch, dataset_func):
    X_y_tuple = dataset_func(return_X_y=True)
    assert isinstance(X_y_tuple, tuple)
    assert X_y_tuple[0].shape == bunch.data.shape
    assert X_y_tuple[1].shape == bunch.target.shape


def check_as_frame(bunch, dataset_func,
                   expected_data_dtype=None, expected_target_dtype=None):
    pd = pytest.importorskip('pandas')
    frame_bunch = dataset_func(as_frame=True)
    assert hasattr(frame_bunch, 'frame')
    assert isinstance(frame_bunch.frame, pd.DataFrame)
    assert isinstance(frame_bunch.data, pd.DataFrame)
    assert frame_bunch.data.shape == bunch.data.shape
    if frame_bunch.target.ndim > 1:
        assert isinstance(frame_bunch.target, pd.DataFrame)
    else:
        assert isinstance(frame_bunch.target, pd.Series)
    assert frame_bunch.target.shape[0] == bunch.target.shape[0]
    if expected_data_dtype is not None:
        assert np.all(frame_bunch.data.dtypes == expected_data_dtype)
    if expected_target_dtype is not None:
        assert np.all(frame_bunch.target.dtypes == expected_target_dtype)


def _has_network():
    return bool(os.environ.get("SKLEARN_SKIP_NETWORK_TESTS", False))


def _generate_func_supporting_param(param, dataset_type=("load", "fetch")):
    markers_known_failure = KNOWN_FAILURE.get(param, {})
    for name, obj in inspect.getmembers(sklearn.datasets):
        if not inspect.isfunction(obj):
            continue

        is_dataset_type = any([name.startswith(t) for t in dataset_type])
        is_support_param = param in inspect.signature(obj).parameters
        if is_dataset_type and is_support_param:
            # check if we should skip if we don't have network support
            marks = [pytest.mark.skipif(
                condition=name.startswith("fetch") and _has_network(),
                reason="Skip because fetcher requires internet network",
            )]
            marks.append(markers_known_failure.get(name, pytest.mark.basic))

            yield pytest.param(name, obj, marks=marks)


@pytest.mark.parametrize(
    "name, dataset_func", _generate_func_supporting_param("return_X_y")
)
def test_common_check_return_X_y(name, dataset_func):
    bunch = dataset_func()
    check_return_X_y(bunch, dataset_func)


@pytest.mark.parametrize(
    "name, dataset_func", _generate_func_supporting_param("as_frame")
)
def test_common_check_as_frame(name, dataset_func):
    bunch = dataset_func()
    check_as_frame(bunch, dataset_func)


@pytest.mark.parametrize(
    "name, dataset_func", _generate_func_supporting_param("as_frame")
)
def test_common_check_pandas_dependency(name, dataset_func):
    check_pandas_dependency_message(dataset_func)
