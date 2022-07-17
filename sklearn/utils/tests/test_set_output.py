import pytest

import numpy as np
from scipy.sparse import csr_matrix
from numpy.testing import assert_array_equal

from sklearn._config import config_context
from sklearn.utils.set_output import _wrap_in_pandas_container
from sklearn.utils.set_output import safe_set_output
from sklearn.utils.set_output import SetOutputMixin
from sklearn.utils.set_output import get_output_config


def test__wrap_in_pandas_container_dense():
    """Check _wrap_in_pandas_container for dense data."""
    pd = pytest.importorskip("pandas")
    X = np.asarray([[1, 0, 3], [0, 0, 1]])
    columns = np.asarray(["f0", "f1", "f2"], dtype=object)
    index = np.asarray([0, 1])

    dense_named = _wrap_in_pandas_container(X, columns=lambda: columns, index=index)
    assert isinstance(dense_named, pd.DataFrame)
    assert_array_equal(dense_named.columns, columns)
    assert_array_equal(dense_named.index, index)


def test__wrap_in_pandas_container_error_validation():
    """Check errors in _wrap_in_pandas_container."""
    X = np.asarray([[1, 0, 3], [0, 0, 1]])
    X_csr = csr_matrix(X)
    match = "Pandas output does not support sparse data"
    with pytest.raises(ValueError, match=match):
        _wrap_in_pandas_container(X_csr)


class EstimatorWithoutSetOutputAndWithoutTransform:
    pass


class EstimatorNoSetOutputWithTransform:
    def transform(self, X, y=None):
        return X  # pragma: no cover


class EstimatorWithSetOutput(SetOutputMixin):
    def transform(self, X, y=None):
        return X


def test_safe_set_output():
    """Check safe_set_output works as expected."""

    # Estimator without transform will not raise when setting set_output for
    # transform.
    est = EstimatorWithoutSetOutputAndWithoutTransform()
    safe_set_output(est, transform="pandas")

    # Estimator with transform without set_output would raise
    est = EstimatorNoSetOutputWithTransform()
    with pytest.raises(ValueError):
        safe_set_output(est, transform="pandas")

    est = EstimatorWithSetOutput()
    safe_set_output(est, transform="pandas")
    config = get_output_config(est, "transform")
    assert config["dense"] == "pandas"

    safe_set_output(est, transform="default")
    config = get_output_config(est, "transform")
    assert config["dense"] == "default"

    # transform is None is a noop, so the config remains "default"
    safe_set_output(est, transform=None)
    config = get_output_config(est, "transform")
    assert config["dense"] == "default"


def test_safe_set_output_error():
    """Check transform with invalid config."""
    X = np.asarray([[1, 0, 3], [0, 0, 1]])

    est = EstimatorWithSetOutput()
    safe_set_output(est, transform="bad")

    msg = "output config must be 'default'"
    with pytest.raises(ValueError, match=msg):
        est.transform(X)


def test_set_output_method():
    """Check that the output is pandas."""
    pd = pytest.importorskip("pandas")

    est = EstimatorWithSetOutput()
    est.set_output(transform="pandas")

    X = np.asarray([[1, 0, 3], [0, 0, 1]])
    X_trans = est.transform(X)
    assert isinstance(X_trans, pd.DataFrame)


def test_set_output_method_error():
    """Check transform fails with invalid transform."""

    est = EstimatorWithSetOutput()
    est.set_output(transform="bad")
    X = np.asarray([[1, 0, 3], [0, 0, 1]])

    msg = "output config must be 'default'"
    with pytest.raises(ValueError, match=msg):
        est.transform(X)


def test_get_output_config():
    """Check get_output_config works as expected."""

    with config_context(transform_output="pandas"):
        est = EstimatorNoSetOutputWithTransform()
        config = get_output_config(est, "transform")
        assert config["dense"] == "pandas"

        est = EstimatorWithSetOutput()
        # If estimator has not config, use global config
        config = get_output_config(est, "transform")
        assert config["dense"] == "pandas"

        # If estimator has a config, use local config
        est.set_output(transform="default")
        config = get_output_config(est, "transform")
        assert config["dense"] == "default"

    est.set_output(transform="pandas")
    config = get_output_config(est, "transform")
    assert config["dense"] == "pandas"
