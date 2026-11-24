"""Check the ``sample_weight`` arrays real meta-estimators pass to consumers.

The meta-estimators are scikit-learn objects. A wrapper around each leaf
routes that leaf's metadata and records the ``sample_weight`` it forwards
to ``fit``. The wrapper requests nothing itself.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
from joblib import parallel_config
from numpy.testing import assert_array_equal

from sklearn import config_context
from sklearn.base import BaseEstimator, clone
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.utils.metadata_routing import (
    MetadataRouter,
    MethodMapping,
    process_routing,
)


class RecordFitSampleWeight(BaseEstimator):
    """Route metadata to a leaf estimator and record the forwarded ``sample_weight``.

    This wrapper requests no metadata. Parents pass ``sample_weight`` because
    the wrapped estimator requests it, and ``fit`` records the array it then
    forwards.

    Parameters
    ----------
    estimator : estimator
        Estimator whose metadata requests are forwarded.

    records : list
        Appended with a copy of each ``sample_weight`` forwarded to ``fit``.
        Clones share this list.
    """

    def __init__(self, estimator, records):
        self.estimator = estimator
        self.records = records

    def __sklearn_clone__(self):
        """Clone the wrapper and keep ``records`` as the same list.

        ``clone`` copies constructor parameters, so each cross-validation copy
        would otherwise append to its own list. The test reads the list it
        passed in.
        """
        cloned = super().__sklearn_clone__()
        cloned.records = self.records
        return cloned

    def get_metadata_routing(self):
        """Route each method's metadata to the same method of the estimator."""
        return MetadataRouter(owner=self).add(
            estimator=self.estimator,
            method_mapping=MethodMapping()
            .add(caller="fit", callee="fit")
            .add(caller="transform", callee="transform")
            .add(caller="predict", callee="predict")
            .add(caller="predict_proba", callee="predict_proba")
            .add(caller="decision_function", callee="decision_function"),
        )

    def fit(self, X, y=None, **params):
        routed_params = process_routing(self, "fit", **params)
        fit_params = routed_params["estimator"]["fit"]
        sample_weight = fit_params.get("sample_weight")
        if sample_weight is None:
            self.records.append(None)
        else:
            self.records.append(np.asarray(sample_weight).copy())

        self.estimator_ = clone(self.estimator)
        self.estimator_.fit(X, y, **fit_params)
        if hasattr(self.estimator_, "classes_"):
            self.classes_ = self.estimator_.classes_
        return self

    def __sklearn_tags__(self):
        return self.estimator.__sklearn_tags__()

    def transform(self, X, **params):
        routed_params = process_routing(self, "transform", **params)
        return self.estimator_.transform(X, **routed_params["estimator"]["transform"])

    def predict(self, X, **params):
        routed_params = process_routing(self, "predict", **params)
        return self.estimator_.predict(X, **routed_params["estimator"]["predict"])

    def predict_proba(self, X, **params):
        routed_params = process_routing(self, "predict_proba", **params)
        return self.estimator_.predict_proba(
            X, **routed_params["estimator"]["predict_proba"]
        )

    def decision_function(self, X, **params):
        routed_params = process_routing(self, "decision_function", **params)
        return self.estimator_.decision_function(
            X, **routed_params["estimator"]["decision_function"]
        )


def _assert_weights(records, expected):
    assert len(records) == len(expected)
    for recorded, target in zip(records, expected):
        assert_array_equal(recorded, target)


@config_context(enable_metadata_routing=True, enable_metadata_auto_requests=True)
def test_pipeline_gridsearch_passes_split_sample_weight():
    """Search passes each split's train weights to both pipeline steps.

    ``StandardScaler`` and ``Ridge`` are real estimators. Each is wrapped only
    so the test can see the array routed to that estimator. No ``set_*_request``
    is used. One parameter is searched, so each split is fitted once, and
    ``refit`` fits both steps again on the full ``sample_weight``.
    """
    sample_weight = np.array([1, 0, 2, 1, 3, 0, 4, 2])
    rng = np.random.RandomState(0)
    X = rng.randn(sample_weight.shape[0], 3)
    y = rng.randn(sample_weight.shape[0])
    cv = KFold(n_splits=2)
    splits = list(cv.split(X, y))

    scaler_records = []
    model_records = []
    search = GridSearchCV(
        Pipeline(
            [
                ("scale", RecordFitSampleWeight(StandardScaler(), scaler_records)),
                ("model", RecordFitSampleWeight(Ridge(), model_records)),
            ]
        ),
        param_grid={"model__estimator__alpha": [1.0]},
        scoring="neg_mean_squared_error",
        cv=cv,
        n_jobs=1,
    )

    # Threading keeps the shared record lists in this process. The loky
    # backend would append inside a worker and the test would not see them.
    with parallel_config(backend="threading", n_jobs=1):
        search.fit(X, y, sample_weight=sample_weight)

    fit_weights = [sample_weight[train] for train, _ in splits]
    fit_weights.append(sample_weight)
    _assert_weights(scaler_records, fit_weights)
    _assert_weights(model_records, fit_weights)


@config_context(enable_metadata_routing=True, enable_metadata_auto_requests=True)
def test_calibrated_classifier_passes_train_sample_weight():
    """``CalibratedClassifierCV`` fits its base estimator on each train fold's weights.

    The base estimator is a real ``LogisticRegression``. The wrapper only
    records the array it routes into ``fit``. With ``ensemble=True`` there is one
    such fit per split and no refit on the full sample.
    """
    sample_weight = np.array([1, 0, 2, 1, 3, 0, 4, 2])
    rng = np.random.RandomState(0)
    X = rng.randn(sample_weight.shape[0], 3)
    y = np.array([0, 1, 0, 1, 0, 1, 0, 1])
    cv = KFold(n_splits=2)
    splits = list(cv.split(X, y))

    records = []
    calibrated = CalibratedClassifierCV(
        RecordFitSampleWeight(LogisticRegression(), records),
        method="sigmoid",
        cv=cv,
        ensemble=True,
        n_jobs=1,
    )

    with parallel_config(backend="threading", n_jobs=1):
        calibrated.fit(X, y, sample_weight=sample_weight)

    _assert_weights(records, [sample_weight[train] for train, _ in splits])
