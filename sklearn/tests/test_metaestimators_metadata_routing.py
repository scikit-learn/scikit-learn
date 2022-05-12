import numpy as np
import pytest
from sklearn.base import RegressorMixin, ClassifierMixin, BaseEstimator
from sklearn.multioutput import MultiOutputRegressor, MultiOutputClassifier
from sklearn.utils.metadata_routing import MetadataRouter
from sklearn.tests.test_metadata_routing import (
    record_metadata,
    check_recorded_metadata,
    assert_request_is_empty,
)

N, M = 100, 4
X = np.random.rand(N, M)
y = np.random.randint(0, 2, size=N)
y_multi = np.random.randint(0, 2, size=(N, 3))
metadata = np.random.randint(0, 10, size=N)
sample_weight = np.random.rand(N)


class ConsumingRegressor(RegressorMixin, BaseEstimator):
    """A regressor consuming metadata."""

    def __init__(self, **kwargs):
        for param, value in kwargs.items():
            setattr(self, param, value)

    def partial_fit(self, X, y, sample_weight=None, metadata=None):
        record_metadata(
            self, "partial_fit", sample_weight=sample_weight, metadata=metadata
        )
        return self

    def fit(self, X, y, sample_weight=None, metadata=None):
        record_metadata(self, "fit", sample_weight=sample_weight, metadata=metadata)
        return self

    def predict(self, X, y, sample_weight=None, metadata=None):
        record_metadata(self, "predict", sample_weight=sample_weight, metadata=metadata)
        return np.zeros(shape=(len(X)))


class ConsumingClassifier(ClassifierMixin, BaseEstimator):
    """A regressor consuming metadata."""

    def __init__(self, **kwargs):
        for param, value in kwargs.items():
            setattr(self, param, value)

    def partial_fit(self, X, y, sample_weight=None, metadata=None):
        record_metadata(
            self, "partial_fit", sample_weight=sample_weight, metadata=metadata
        )
        self.classes_ = [1]
        return self

    def fit(self, X, y, sample_weight=None, metadata=None):
        record_metadata(self, "fit", sample_weight=sample_weight, metadata=metadata)
        self.classes_ = [1]
        return self

    def predict(self, X, y, sample_weight=None, metadata=None):
        record_metadata(self, "predict", sample_weight=sample_weight, metadata=metadata)
        return np.zeros(shape=(len(X)))

    def predict_proba(self, X, y, sample_weight=None, metadata=None):
        record_metadata(
            self, "predict_proba", sample_weight=sample_weight, metadata=metadata
        )
        return np.zeros(shape=(len(X)))

    def predict_log_proba(self, X, y, sample_weight=None, metadata=None):
        record_metadata(
            self, "predict_log_proba", sample_weight=sample_weight, metadata=metadata
        )
        return np.zeros(shape=(len(X)))


def get_empty_metaestimators():
    yield MultiOutputRegressor(estimator=ConsumingRegressor())
    yield MultiOutputClassifier(estimator=ConsumingClassifier())


@pytest.mark.parametrize(
    "metaestimator",
    get_empty_metaestimators(),
    ids=[str(x) for x in get_empty_metaestimators()],
)
def test_default_request(metaestimator):
    # Check that by default request is empty and the right type
    assert_request_is_empty(metaestimator.get_metadata_routing())
    assert isinstance(metaestimator.get_metadata_routing(), MetadataRouter)


@pytest.mark.parametrize(
    "MultiOutput, Estimator",
    [
        (MultiOutputClassifier, ConsumingClassifier),
        (MultiOutputRegressor, ConsumingRegressor),
    ],
    ids=["Classifier", "Regressor"],
)
def test_multioutput_metadata_routing(MultiOutput, Estimator):
    # Check routing of metadata
    metaest = MultiOutput(Estimator())
    with pytest.warns(
        FutureWarning,
        match=(
            "You are passing metadata for which the request values are not explicitly"
            " set. From version 1.4 this results in the following error"
        ),
    ):
        metaest.fit(X, y_multi, sample_weight=sample_weight, metadata=metadata)
        check_recorded_metadata(
            metaest.estimators_[0],
            "fit",
            sample_weight=sample_weight,
            metadata=metadata,
        )

    metaest = MultiOutput(
        Estimator()
        .set_fit_request(sample_weight=True, metadata="alias")
        .set_partial_fit_request(sample_weight=True, metadata=True)
    ).fit(X, y_multi, alias=metadata)
    check_recorded_metadata(
        metaest.estimators_[0], "fit", sample_weight=None, metadata=metadata
    )

    metaest.partial_fit(X, y_multi, metadata=metadata)
    check_recorded_metadata(
        metaest.estimators_[0], "partial_fit", sample_weight=None, metadata=metadata
    )
