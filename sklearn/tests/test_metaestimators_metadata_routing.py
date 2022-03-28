import numpy as np
from sklearn.base import RegressorMixin, BaseEstimator
from sklearn.tests.test_metadata_routing import record_metadata, check_recorded_metadata


class RegressorSubEstimator(RegressorMixin, BaseEstimator):
    """A regressor consuming metadata."""

    def __init__(self, **kwargs):
        for param, value in kwargs.items():
            setattr(self, param, value)

    def fit(self, X, y, sample_weight=None, metadata=None):
        record_metadata(self, "fit", sample_weight=sample_weight, metadata=metadata)
        return self

    def predict(self, X, y, sample_weight=None, metadata=None):
        record_metadata(self, "fit", sample_weight=sample_weight, metadata=metadata)
        return np.zeros(shape=(len(X)))
