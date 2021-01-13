"""Test loaders for common functionality."""
import inspect

import pytest
import numpy as np

import sklearn.datasets


def is_pillow_installed():
    try:
        import PIL  # noqa
        return True
    except ImportError:
        return False


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

    # Test for return_X_y and as_frame=True
    frame_X, frame_y = dataset_func(as_frame=True, return_X_y=True)
    assert isinstance(frame_X, pd.DataFrame)
    if frame_y.ndim > 1:
        assert isinstance(frame_X, pd.DataFrame)
    else:
        assert isinstance(frame_y, pd.Series)


@pytest.fixture(params=range(7),
                ids=["20newsgroups", "20newsgroups_vectorized",
                     "california_housing", "covtype",
                     "kddcup99", "olivetti_faces", "rcv1"])
def fetch_dataset(request,
                  fetch_20newsgroups_fxt, fetch_20newsgroups_vectorized_fxt,
                  fetch_california_housing_fxt, fetch_covtype_fxt,
                  fetch_kddcup99_fxt, fetch_olivetti_faces_fxt,
                  fetch_rcv1_fxt):
    _fetch_datasets = [
        fetch_20newsgroups_fxt, fetch_20newsgroups_vectorized_fxt,
        fetch_california_housing_fxt, fetch_covtype_fxt, fetch_kddcup99_fxt,
        fetch_olivetti_faces_fxt, fetch_rcv1_fxt
    ]
    return _fetch_datasets[request.param]


def test_common_check_return_X_y_fetch(fetch_dataset):
    """Test return_X_y for fetch_* functions."""
    name = fetch_dataset.__name__
    if name == "fetch_20newsgroups":
        pytest.skip("X is a list and does not have a shape argument")
    elif name == "fetch_lfw_people" and not is_pillow_installed():
        pytest.skip("pillow is not installed")

    bunch = fetch_dataset()
    check_return_X_y(bunch, fetch_dataset)


def test_common_check_as_frame_fetch(fetch_dataset):
    """Test as_frame for fetch_* functions."""
    name = fetch_dataset.__name__
    if "as_frame" not in inspect.signature(fetch_dataset).parameters:
        pytest.skip(f"The as_frame keyword is not defined for {name}")

    bunch = fetch_dataset()
    check_as_frame(bunch, fetch_dataset)


def test_common_check_pandas_dependency_fetch(fetch_dataset):
    """Test pandas dependency messag for fetch_* functions."""
    check_pandas_dependency_message(fetch_dataset)


def _generate_load_func_supporting_param(param):
    """Generate datasets for load_* functions."""
    return inspect.getmembers(
        sklearn.datasets,
        predicate=lambda f:
            inspect.isfunction(f) and
            f.__name__.startswith("load_") and
            param in inspect.signature(f).parameters)


@pytest.mark.parametrize(
    "name, dataset_func", _generate_load_func_supporting_param("return_X_y")
)
def test_common_check_return_X_y(name, dataset_func):
    bunch = dataset_func()
    check_return_X_y(bunch, dataset_func)


@pytest.mark.parametrize(
    "name, dataset_func", _generate_load_func_supporting_param("as_frame")
)
def test_common_check_as_frame(name, dataset_func):
    bunch = dataset_func()
    check_as_frame(bunch, dataset_func)


@pytest.mark.parametrize(
    "name, dataset_func", _generate_load_func_supporting_param("as_frame")
)
def test_common_check_pandas_dependency(name, dataset_func):
    check_pandas_dependency_message(dataset_func)
