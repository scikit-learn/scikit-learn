# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import re

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from sklearn.base import (
    BaseEstimator,
    is_classifier,
    is_clusterer,
    is_outlier_detector,
    is_regressor,
)
from sklearn.cluster import KMeans
from sklearn.compose import make_column_transformer
from sklearn.datasets import make_classification, make_regression
from sklearn.exceptions import UnsetMetadataPassedError
from sklearn.frozen import FrozenEstimator
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.neighbors import LocalOutlierFactor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import RobustScaler, StandardScaler

REGRESSION_DATASET = make_regression()
CLASSIFICATION_DATASET = make_classification()


@pytest.mark.parametrize(
    "estimator, dataset",
    [
        (LinearRegression(), REGRESSION_DATASET),
        (LogisticRegression(), CLASSIFICATION_DATASET),
        (make_pipeline(StandardScaler(), LinearRegression()), REGRESSION_DATASET),
        (make_pipeline(StandardScaler(), LogisticRegression()), CLASSIFICATION_DATASET),
        (StandardScaler(), REGRESSION_DATASET),
        (KMeans(), REGRESSION_DATASET),
        (LocalOutlierFactor(), REGRESSION_DATASET),
        (
            make_column_transformer(
                (StandardScaler(), [0]),
                (RobustScaler(), [1]),
            ),
            REGRESSION_DATASET,
        ),
    ],
)
def test_frozen_methods(estimator, dataset):
    """Test that frozen.fit doesn't do anything, and that all other methods are
    exposed by the frozen estimator and return the same values as the estimator.
    """
    X, y = dataset
    estimator.fit(X, y)
    frozen = FrozenEstimator(estimator)
    # this should be no-op
    frozen.fit([[1]], [1])

    methods = [
        "predict",
        "predict_proba",
        "predict_log_proba",
        "decision_function",
        "transform",
    ]
    for method in methods:
        if hasattr(estimator, method):
            assert_array_equal(
                getattr(estimator, method)(X), getattr(frozen, method)(X)
            )

    assert is_classifier(estimator) == is_classifier(frozen)
    assert is_regressor(estimator) == is_regressor(frozen)
    assert is_clusterer(estimator) == is_clusterer(frozen)
    assert is_outlier_detector(estimator) == is_outlier_detector(frozen)


@pytest.mark.usefixtures("enable_slep006")
def test_frozen_metadata_routing():
    """Test that metadata routing works with frozen estimators."""

    class ConsumesMetadata(BaseEstimator):
        def __init__(self, on_fit=None, on_predict=None):
            self.on_fit = on_fit
            self.on_predict = on_predict

        def fit(self, X, y, metadata=None):
            if self.on_fit:
                assert metadata is not None
            return self

        def predict(self, X, metadata=None):
            if self.on_predict:
                assert metadata is not None
            return np.ones(len(X))

    X, y = REGRESSION_DATASET
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


def test_composite_fit():
    """Test that calling fit_transform and fit_predict doesn't call fit."""

    class Estimator(BaseEstimator):
        def fit(self, X, y):
            try:
                self._fit_counter += 1
            except AttributeError:
                self._fit_counter = 1
            return self

        def predict(self, X):
            return np.ones(len(X))

        def transform(self, X):
            return X

        def fit_transform(self, X, y=None):
            return X

        def fit_predict(self, X, y=None):
            return np.ones(len(X))

    X, y = CLASSIFICATION_DATASET
    est = Estimator().fit(X, y)
    frozen = FrozenEstimator(est)

    frozen.fit_predict(X, y)
    frozen.fit_transform(X, y)

    assert frozen._fit_counter == 1
