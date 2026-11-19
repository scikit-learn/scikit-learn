"""Check the ``sample_weight`` arrays real meta-estimators pass to consumers.

The meta-estimators are scikit-learn objects. Only the leaves record the
array they receive: a wrapper around a real estimator, and a scorer. The
wrapper stores that array and then calls the wrapped estimator's ``fit``.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from functools import partial

import numpy as np
from joblib import parallel_config
from numpy.testing import assert_array_equal

from sklearn import config_context
from sklearn.base import BaseEstimator, clone
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import mean_squared_error
from sklearn.metrics._scorer import _Scorer
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


class _SharedList(list):
    """List that ``clone`` keeps as the same object.

    ``clone`` deep-copies constructor parameters. Cross-validation fits those
    copies, so records have to land on the list that the test still holds.
    """

    def __deepcopy__(self, memo):
        memo[id(self)] = self
        return self

    def __copy__(self):
        return self


class RecordFitSampleWeight(BaseEstimator):
    """Record ``sample_weight`` and fit the wrapped estimator with it.

    Parameters
    ----------
    estimator : estimator
        Estimator fitted with the ``sample_weight`` this wrapper receives.

    records : list
        Appended with a copy of each ``sample_weight`` passed to ``fit``.
    """

    def __init__(self, estimator, records):
        self.estimator = estimator
        self.records = records

    def fit(self, X, y=None, sample_weight=None, **fit_params):
        if sample_weight is None:
            self.records.append(None)
        else:
            self.records.append(np.asarray(sample_weight).copy())

        self.estimator_ = clone(self.estimator)
        if sample_weight is not None:
            fit_params = {**fit_params, "sample_weight": sample_weight}
        self.estimator_.fit(X, y, **fit_params)
        if hasattr(self.estimator_, "classes_"):
            self.classes_ = self.estimator_.classes_
        return self

    def __sklearn_tags__(self):
        return self.estimator.__sklearn_tags__()

    def transform(self, X):
        return self.estimator_.transform(X)

    def predict(self, X):
        return self.estimator_.predict(X)

    def predict_proba(self, X):
        return self.estimator_.predict_proba(X)

    def decision_function(self, X):
        return self.estimator_.decision_function(X)


def _recording_score(y_true, y_pred, sample_weight=None, records=None):
    """Scorer metric that stores ``sample_weight`` then computes squared error."""
    records.append(np.asarray(sample_weight).copy())
    return mean_squared_error(y_true, y_pred, sample_weight=sample_weight)


class RecordingScorer(_Scorer):
    """Scorer that auto-requests ``sample_weight`` and records it.

    Parameters
    ----------
    records : list
        Appended with a copy of each ``sample_weight`` passed to the metric.
    """

    def __init__(self, records):
        super().__init__(
            score_func=partial(_recording_score, records=records),
            sign=-1,
            kwargs={},
            response_method="predict",
        )


def _assert_weights(records, expected):
    assert len(records) == len(expected)
    for recorded, target in zip(records, expected):
        assert_array_equal(recorded, target)


@config_context(enable_metadata_routing=True, enable_metadata_auto_requests=True)
def test_pipeline_gridsearch_passes_split_sample_weight():
    """Search passes train weights to each pipeline step and test weights to the scorer.

    ``StandardScaler`` and ``Ridge`` are real estimators. Each is wrapped only
    so the test can see the array ``Pipeline`` forwarded. No ``set_*_request``
    is used. One parameter is searched, so each split is fitted once, and
    ``refit`` fits both steps again on the full ``sample_weight``.
    """
    sample_weight = np.array([1, 0, 2, 1, 3, 0, 4, 2])
    rng = np.random.RandomState(0)
    X = rng.randn(sample_weight.shape[0], 3)
    y = rng.randn(sample_weight.shape[0])
    cv = KFold(n_splits=2)
    splits = list(cv.split(X, y))

    scaler_records = _SharedList()
    model_records = _SharedList()
    scorer_records = _SharedList()
    search = GridSearchCV(
        Pipeline(
            [
                ("scale", RecordFitSampleWeight(StandardScaler(), scaler_records)),
                ("model", RecordFitSampleWeight(Ridge(), model_records)),
            ]
        ),
        param_grid={"model__estimator__alpha": [1.0]},
        scoring=RecordingScorer(scorer_records),
        cv=cv,
        n_jobs=1,
    )

    # Threading keeps the shared record lists in this process. The loky
    # backend would append inside a worker and the test would not see them.
    with parallel_config(backend="threading", n_jobs=1):
        search.fit(X, y, sample_weight=sample_weight)

    fit_weights = [sample_weight[train] for train, _ in splits]
    fit_weights.append(sample_weight)
    score_weights = [sample_weight[test] for _, test in splits]
    _assert_weights(scaler_records, fit_weights)
    _assert_weights(model_records, fit_weights)
    _assert_weights(scorer_records, score_weights)


@config_context(enable_metadata_routing=True, enable_metadata_auto_requests=True)
def test_calibrated_classifier_passes_train_sample_weight():
    """``CalibratedClassifierCV`` fits its base estimator on each train fold's weights.

    The base estimator is a real ``LogisticRegression``. The wrapper only
    records the array passed into ``fit``. With ``ensemble=True`` there is one
    such fit per split and no refit on the full sample.
    """
    sample_weight = np.array([1, 0, 2, 1, 3, 0, 4, 2])
    rng = np.random.RandomState(0)
    X = rng.randn(sample_weight.shape[0], 3)
    y = np.array([0, 1, 0, 1, 0, 1, 0, 1])
    cv = KFold(n_splits=2)
    splits = list(cv.split(X, y))

    records = _SharedList()
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
