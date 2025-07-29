# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import re

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from sklearn import config_context
from sklearn.base import (
    BaseEstimator,
    clone,
    is_classifier,
    is_clusterer,
    is_outlier_detector,
    is_regressor,
)
from sklearn.cluster import KMeans
from sklearn.compose import make_column_transformer
from sklearn.datasets import make_classification, make_regression
from sklearn.exceptions import NotFittedError, UnsetMetadataPassedError
from sklearn.frozen import FrozenEstimator
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.neighbors import LocalOutlierFactor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import RobustScaler, StandardScaler
from sklearn.utils._testing import set_random_state
from sklearn.utils.validation import check_is_fitted


@pytest.fixture
def regression_dataset():
    return make_regression()


@pytest.fixture
def classification_dataset():
    return make_classification()


@pytest.mark.parametrize(
    "estimator, dataset",
    [
        (LinearRegression(), "regression_dataset"),
        (LogisticRegression(), "classification_dataset"),
        (make_pipeline(StandardScaler(), LinearRegression()), "regression_dataset"),
        (
            make_pipeline(StandardScaler(), LogisticRegression()),
            "classification_dataset",
        ),
        (StandardScaler(), "regression_dataset"),
        (KMeans(), "regression_dataset"),
        (LocalOutlierFactor(), "regression_dataset"),
        (
            make_column_transformer(
                (StandardScaler(), [0]),
                (RobustScaler(), [1]),
            ),
            "regression_dataset",
        ),
    ],
)
@pytest.mark.parametrize(
    "method",
    ["predict", "predict_proba", "predict_log_proba", "decision_function", "transform"],
)
def test_frozen_methods(estimator, dataset, request, method):
    """Test that frozen.fit doesn't do anything, and that all other methods are
    exposed by the frozen estimator and return the same values as the estimator.
    """
    X, y = request.getfixturevalue(dataset)
    set_random_state(estimator)
    estimator.fit(X, y)
    frozen = FrozenEstimator(estimator)
    # this should be no-op
    frozen.fit([[1]], [1])

    if hasattr(estimator, method):
        assert_array_equal(getattr(estimator, method)(X), getattr(frozen, method)(X))

    assert is_classifier(estimator) == is_classifier(frozen)
    assert is_regressor(estimator) == is_regressor(frozen)
    assert is_clusterer(estimator) == is_clusterer(frozen)
    assert is_outlier_detector(estimator) == is_outlier_detector(frozen)


@config_context(enable_metadata_routing=True)
def test_frozen_metadata_routing(regression_dataset):
    """Test that metadata routing works with frozen estimators."""

    class ConsumesMetadata(BaseEstimator):
        def __init__(self, on_fit=None, on_predict=None):
            self.on_fit = on_fit
            self.on_predict = on_predict

        def fit(self, X, y, metadata=None):
            if self.on_fit:
                assert metadata is not None
            self.fitted_ = True
            return self

        def predict(self, X, metadata=None):
            if self.on_predict:
                assert metadata is not None
            return np.ones(len(X))

    X, y = regression_dataset
    pipeline = make_pipeline(
        ConsumesMetadata(on_fit=True, on_predict=True)
        .set_fit_request(metadata=True)
        .set_predict_request(metadata=True)
    )

    pipeline.fit(X, y, metadata="test")
    frozen = FrozenEstimator(pipeline)
    pipeline.predict(X, metadata="test")
    frozen.predict(X, metadata="test")

    frozen["consumesmetadata"].set_predict_request(metadata=False)
    with pytest.raises(
        TypeError,
        match=re.escape(
            "Pipeline.predict got unexpected argument(s) {'metadata'}, which are not "
            "routed to any object."
        ),
    ):
        frozen.predict(X, metadata="test")

    frozen["consumesmetadata"].set_predict_request(metadata=None)
    with pytest.raises(UnsetMetadataPassedError):
        frozen.predict(X, metadata="test")


def test_composite_fit(classification_dataset):
    """Test that calling fit_transform and fit_predict doesn't call fit."""

    class Estimator(BaseEstimator):
        def fit(self, X, y):
            try:
                self._fit_counter += 1
            except AttributeError:
                self._fit_counter = 1
            return self

        def fit_transform(self, X, y=None):
            # only here to test that it doesn't get called
            ...  # pragma: no cover

        def fit_predict(self, X, y=None):
            # only here to test that it doesn't get called
            ...  # pragma: no cover

    X, y = classification_dataset
    est = Estimator().fit(X, y)
    frozen = FrozenEstimator(est)

    with pytest.raises(AttributeError):
        frozen.fit_predict(X, y)
    with pytest.raises(AttributeError):
        frozen.fit_transform(X, y)

    assert frozen._fit_counter == 1


def test_clone_frozen(regression_dataset):
    """Test that cloning a frozen estimator keeps the frozen state."""
    X, y = regression_dataset
    estimator = LinearRegression().fit(X, y)
    frozen = FrozenEstimator(estimator)
    cloned = clone(frozen)
    assert cloned.estimator is estimator


def test_check_is_fitted(regression_dataset):
    """Test that check_is_fitted works on frozen estimators."""
    X, y = regression_dataset

    estimator = LinearRegression()
    frozen = FrozenEstimator(estimator)
    with pytest.raises(NotFittedError):
        check_is_fitted(frozen)

    estimator = LinearRegression().fit(X, y)
    frozen = FrozenEstimator(estimator)
    check_is_fitted(frozen)


def test_frozen_tags():
    """Test that frozen estimators have the same tags as the original estimator
    except for the skip_test tag."""

    class Estimator(BaseEstimator):
        def __sklearn_tags__(self):
            tags = super().__sklearn_tags__()
            tags.input_tags.categorical = True
            return tags

    estimator = Estimator()
    frozen = FrozenEstimator(estimator)
    frozen_tags = frozen.__sklearn_tags__()
    estimator_tags = estimator.__sklearn_tags__()

    assert frozen_tags._skip_test is True
    assert estimator_tags._skip_test is False

    assert estimator_tags.input_tags.categorical is True
    assert frozen_tags.input_tags.categorical is True


def test_frozen_params():
    """Test that FrozenEstimator only exposes the estimator parameter."""
    est = LogisticRegression()
    frozen = FrozenEstimator(est)

    with pytest.raises(ValueError, match="You cannot set parameters of the inner"):
        frozen.set_params(estimator__C=1)

    assert frozen.get_params() == {"estimator": est}

    other_est = LocalOutlierFactor()
    frozen.set_params(estimator=other_est)
    assert frozen.get_params() == {"estimator": other_est}
