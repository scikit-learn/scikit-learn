import pytest

import numpy as np
from scipy.sparse import csr_matrix
from numpy.testing import assert_array_equal

from sklearn._config import config_context, get_config
from sklearn.utils.set_output import _wrap_in_pandas_container
from sklearn.utils.set_output import safe_set_output
from sklearn.utils.set_output import _SetOutputMixin
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


def test__wrap_in_pandas_container_dense_update_columns():
    pd = pytest.importorskip("pandas")
    X_df = pd.DataFrame([[1, 0, 3], [0, 0, 1]], columns=["a", "b", "c"])
    new_columns = np.asarray(["f0", "f1", "f2"], dtype=object)

    new_df = _wrap_in_pandas_container(X_df, columns=new_columns)
    assert_array_equal(new_df.columns, new_columns)


def test__wrap_in_pandas_container_error_validation():
    """Check errors in _wrap_in_pandas_container."""
    X = np.asarray([[1, 0, 3], [0, 0, 1]])
    X_csr = csr_matrix(X)
    match = "Pandas output does not support sparse data"
    with pytest.raises(ValueError, match=match):
        _wrap_in_pandas_container(X_csr, columns=["a", "b", "c"])


class EstimatorWithoutSetOutputAndWithoutTransform:
    pass


class EstimatorNoSetOutputWithTransform:
    def transform(self, X, y=None):
        return X  # pragma: no cover


class EstimatorWithSetOutput(_SetOutputMixin):
    def fit(self, X, y=None):
        self.n_features_in_ = X.shape[1]
        return self

    def transform(self, X, y=None):
        return X

    def get_feature_names_out(self, input_features=None):
        return np.asarray([f"X{i}" for i in range(self.n_features_in_)], dtype=object)


def test_safe_set_output():
    """Check safe_set_output works as expected."""

    # Estimator without transform will not raise when setting set_output for transform.
    est = EstimatorWithoutSetOutputAndWithoutTransform()
    safe_set_output(est, transform="pandas")

    # Estimator with transform but without set_output will raise
    est = EstimatorNoSetOutputWithTransform()
    with pytest.raises(ValueError, match="Unable to configure output"):
        safe_set_output(est, transform="pandas")

    est = EstimatorWithSetOutput().fit(np.asarray([[1, 2, 3]]))
    safe_set_output(est, transform="pandas")
    config = get_output_config("transform", est)
    assert config["dense"] == "pandas"

    safe_set_output(est, transform="default")
    config = get_output_config("transform", est)
    assert config["dense"] == "default"

    # transform is None is a noop, so the config remains "default"
    safe_set_output(est, transform=None)
    config = get_output_config("transform", est)
    assert config["dense"] == "default"


class EstimatorNoSetOutputWithTransformNoFeatureNamesOut(_SetOutputMixin):
    def transform(self, X, y=None):
        return X  # pragma: no cover


def test_set_output_mixin():
    """Estimator without get_feature_names_out does not define `set_output`."""
    est = EstimatorNoSetOutputWithTransformNoFeatureNamesOut()
    assert not hasattr(est, "set_output")


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

    X = np.asarray([[1, 0, 3], [0, 0, 1]])
    est = EstimatorWithSetOutput().fit(X)
    est.set_output(transform="pandas")

    X_trans = est.transform(X)
    assert isinstance(X_trans, pd.DataFrame)


def test_set_output_method_error():
    """Check transform fails with invalid transform."""

    X = np.asarray([[1, 0, 3], [0, 0, 1]])
    est = EstimatorWithSetOutput().fit(X)
    est.set_output(transform="bad")

    msg = "output config must be 'default'"
    with pytest.raises(ValueError, match=msg):
        est.transform(X)


def test_get_output_config():
    """Check get_output_config works as expected."""

    # Without a configuration set, the global config is used
    global_config = get_config()["transform_output"]
    config = get_output_config("transform")
    assert config["dense"] == global_config

    with config_context(transform_output="pandas"):
        # with estimator=None, the global config is used
        config = get_output_config("transform")
        assert config["dense"] == "pandas"

        est = EstimatorNoSetOutputWithTransform()
        config = get_output_config("transform", est)
        assert config["dense"] == "pandas"

        est = EstimatorWithSetOutput()
        # If estimator has not config, use global config
        config = get_output_config("transform", est)
        assert config["dense"] == "pandas"

        # If estimator has a config, use local config
        est.set_output(transform="default")
        config = get_output_config("transform", est)
        assert config["dense"] == "default"

    est.set_output(transform="pandas")
    config = get_output_config("transform", est)
    assert config["dense"] == "pandas"


class EstimatorWithSetOutputNoAutoWrap(_SetOutputMixin, auto_wrap_output=False):
    def transform(self, X, y=None):
        return X


def test_get_output_auto_wrap_false():
    """Check that auto_wrap_output=False does not wrap."""
    est = EstimatorWithSetOutputNoAutoWrap()
    assert not hasattr(est, "set_output")

    X = np.asarray([[1, 0, 3], [0, 0, 1]])
    assert X is est.transform(X)


def test_auto_wrap_output_errors_with_incorrect_input():
    msg = "auto_wrap_output should be a bool"
    with pytest.raises(ValueError, match=msg):

        class BadEstimator(_SetOutputMixin, auto_wrap_output="bad_parameter"):
            pass
