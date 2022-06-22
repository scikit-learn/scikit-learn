import pytest

import numpy as np
from scipy.sparse import csr_matrix
from numpy.testing import assert_array_equal

from sklearn._config import config_context
from sklearn.utils.set_output import make_named_container
from sklearn.utils.set_output import safe_set_output
from sklearn.utils.set_output import SetOutputMixin
from sklearn.utils.set_output import get_output_config


def test_make_named_container_dense():
    """Check make_named_container for dense data."""
    pd = pytest.importorskip("pandas")
    X = np.asarray([[1, 0, 3], [0, 0, 1]])
    columns = np.asarray(["f0", "f1", "f2"], dtype=object)
    index = np.asarray([0, 1])

    X_default = make_named_container(X, dense_container="default")
    assert X_default is X

    dense_named = make_named_container(
        X, columns=lambda: columns, dense_container="pandas", index=index
    )
    assert isinstance(dense_named, pd.DataFrame)
    assert_array_equal(dense_named.columns, columns)
    assert_array_equal(dense_named.index, index)

    # Updates dataframe columns if the input is a dataframe
    new_columns = np.asarray(["g0", "g1", "g2"], dtype=object)
    dense_named_rt = make_named_container(
        dense_named, dense_container="pandas", columns=new_columns
    )
    assert dense_named_rt is dense_named
    assert_array_equal(dense_named_rt.columns, new_columns)


def test_make_named_container_error_validation():
    """Check errors in make_named_container."""
    X = np.asarray([[1, 0, 3], [0, 0, 1]])

    match = "dense_container must be 'default' or 'pandas'"
    with pytest.raises(ValueError, match=match):
        make_named_container(X, dense_container="invalid")

    X_csr = csr_matrix(X)
    match = "Pandas output does not support sparse data"
    with pytest.raises(ValueError, match=match):
        make_named_container(X_csr, dense_container="pandas")


class EstimatorWithoutSetOutputAndWithoutTransform:
    pass


class EstimatorNoSetOutputWithTransform:
    def transform(self, X, y=None):
        return X


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
    assert config.dense == "pandas"

    safe_set_output(est, transform="default")
    config = get_output_config(est, "transform")
    assert config.dense == "default"

    # transform is None is a noop, so the config remains "default"
    safe_set_output(est, transform=None)
    config = get_output_config(est, "transform")
    assert config.dense == "default"


def test_get_output_config():
    """Check get_output_config works as expected."""

    with config_context(output_transform="pandas"):
        est = EstimatorNoSetOutputWithTransform()
        config = get_output_config(est, "transform")
        assert config.dense == "pandas"

        est = EstimatorWithSetOutput()
        # If estimator has not config, use global config
        config = get_output_config(est, "transform")
        assert config.dense == "pandas"

        # If estimator has a config, use local config
        est.set_output(transform="default")
        config = get_output_config(est, "transform")
        assert config.dense == "default"

    est.set_output(transform="pandas")
    config = get_output_config(est, "transform")
    assert config.dense == "pandas"
