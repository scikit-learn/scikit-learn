import re
import warnings

import numpy as np
import pytest
from sklearn.base import RegressorMixin, ClassifierMixin, BaseEstimator
from sklearn.calibration import CalibratedClassifierCV
from sklearn.exceptions import UnsetMetadataPassedError
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


class _Registry(list):
    # This list is used to get a reference to the sub-estimators, which are not
    # necessarily stored on the metaestimator. We need to override __deepcopy__
    # because the sub-estimators are probably cloned, which would result in a
    # new copy of the list.
    def __deepcopy__(self, memo):
        return self


class ConsumingRegressor(RegressorMixin, BaseEstimator):
    """A regressor consuming metadata.

    Parameters
    ----------
    registry : list or None, default=None
        If a list, that estimator will append itself to the list in order to
        have a reference to the estimator later on. Since that reference is not
        required in all tests, registration can be skipped by leaving this value
        as None.

    """
    def __init__(self, registry=None):
        self.registry = registry

    def partial_fit(self, X, y, sample_weight=None, metadata=None):
        if self.registry is not None:
            self.registry.append(self)

        record_metadata(
            self, "partial_fit", sample_weight=sample_weight, metadata=metadata
        )
        return self

    def fit(self, X, y, sample_weight=None, metadata=None):
        if self.registry is not None:
            self.registry.append(self)

        record_metadata(self, "fit", sample_weight=sample_weight, metadata=metadata)
        return self

    def predict(self, X, sample_weight=None, metadata=None):
        if self.registry is not None:
            self.registry.append(self)

        record_metadata(self, "predict", sample_weight=sample_weight, metadata=metadata)
        return np.zeros(shape=(len(X),))


class ConsumingClassifier(ClassifierMixin, BaseEstimator):
    """A classifier consuming metadata.

    Parameters
    ----------
    registry : list or None, default=None
        If a list, that estimator will append itself to the list in order to
        have a reference to the estimator later on. Since that reference is not
        required in all tests, registration can be skipped by leaving this value
        as None.

    """
    def __init__(self, registry=None):
        self.registry = registry

    def partial_fit(self, X, y, sample_weight=None, metadata=None):
        if self.registry is not None:
            self.registry.append(self)

        record_metadata(
            self, "partial_fit", sample_weight=sample_weight, metadata=metadata
        )
        self.classes_ = [0, 1]
        return self

    def fit(self, X, y, sample_weight=None, metadata=None):
        if self.registry is not None:
            self.registry.append(self)

        record_metadata(self, "fit", sample_weight=sample_weight, metadata=metadata)
        self.classes_ = [0, 1]
        return self

    def predict(self, X, sample_weight=None, metadata=None):
        if self.registry is not None:
            self.registry.append(self)

        record_metadata(self, "predict", sample_weight=sample_weight, metadata=metadata)
        return np.zeros(shape=(len(X),))

    def predict_proba(self, X, sample_weight=None, metadata=None):
        if self.registry is not None:
            self.registry.append(self)

        record_metadata(
            self, "predict_proba", sample_weight=sample_weight, metadata=metadata
        )
        return np.asarray([[0.0, 1.0]] * len(X))

    def predict_log_proba(self, X, sample_weight=None, metadata=None):
        if self.registry is not None:
            self.registry.append(self)

        record_metadata(
            self, "predict_log_proba", sample_weight=sample_weight, metadata=metadata
        )
        return np.zeros(shape=(len(X), 2))


METAESTIMATORS = [
    {
        "metaestimator": MultiOutputRegressor,
        "estimator_name": "estimator",
        "estimator": ConsumingRegressor,
        "X": X,
        "y": y_multi,
        "routing_methods": ["fit", "partial_fit"],
        "fit": {"warns_on": ["sample_weight"]},
    },
    {
        "metaestimator": MultiOutputClassifier,
        "estimator_name": "estimator",
        "estimator": ConsumingClassifier,
        "X": X,
        "y": y_multi,
        "routing_methods": ["fit", "partial_fit"],
        "fit": {"warns_on": ["sample_weight"]},
    },
    {
        "metaestimator": CalibratedClassifierCV,
        "estimator_name": "estimator",
        "estimator": ConsumingClassifier,
        "X": X,
        "y": y,
        "routing_methods": ["fit"],
        "fit": {"warns_on": "all"},
        "preserves_metadata": False,
    },
]
"""List containing all metaestimators to be tested and their settings

The keys are as follows:

- metaestimator: The metaestmator to be tested
- estimator_name: The name of the argument for the sub-estimator
- estimator: The sub-estimator
- X: X-data to fit and predict
- y: y-data to fit
- routing_methods: list of all methods to check for routing
- <method_name>: How said method behaves, e.g. "fit": {"warns_on":
  ["sample_weight"]} means fit warns when sample weights are passed but not
  requested.
- preserves_metadata: Whether the metaestimator passes the metadata to the
  sub-estimator without modification or not. If it does, we check that the
  values are identical. If it doesn', no check is performed. TODO Maybe
  something smarter could be done if the data is modified.

"""


@pytest.mark.parametrize(
    "metaestimator",
    METAESTIMATORS,
    ids=[str(row["metaestimator"].__class__.__name__) for row in METAESTIMATORS],
)
def test_default_request(metaestimator):
    # Check that by default request is empty and the right type
    cls = metaestimator["metaestimator"]
    estimator = metaestimator["estimator"]()
    estimator_name = metaestimator["estimator_name"]
    instance = cls(**{estimator_name: estimator})
    assert_request_is_empty(instance.get_metadata_routing())
    assert isinstance(instance.get_metadata_routing(), MetadataRouter)


@pytest.mark.parametrize(
    "metaestimator",
    METAESTIMATORS,
    ids=[str(row["metaestimator"].__class__.__name__) for row in METAESTIMATORS],
)
def test_warning_for_indicated_methods(metaestimator):
    # Check that the indicated methods give a warning
    # TODO: After deprecation period, always error
    cls = metaestimator["metaestimator"]
    registry = _Registry()
    estimator = metaestimator["estimator"](registry=registry)
    estimator_name = metaestimator["estimator_name"]
    instance = cls(**{estimator_name: estimator})
    X = metaestimator["X"]
    y = metaestimator["y"]
    routing_methods = metaestimator["routing_methods"]

    for method_name in routing_methods:
        if not metaestimator.get(method_name, {}).get("warns_on"):
            # this method is not expected to warn
            continue

        warn_msg = (
            "You are passing metadata for which the request values are not explicitly "
            "set: sample_weight, metadata. From version 1.4 this results in the "
            "following error: [sample_weight, metadata] are passed but are not "
            "explicitly set as requested or not for "
            f"{estimator.__class__.__name__}.{method_name}"
        )
        with pytest.warns(FutureWarning, match=re.escape(warn_msg)):
            method = getattr(instance, method_name)
            method(X, y, sample_weight=sample_weight, metadata=metadata)

        if metaestimator.get("preserves_metadata", True):
            check_recorded_metadata(
                registry[-1],
                method_name,
                sample_weight=sample_weight,
                metadata=metadata,
            )


@pytest.mark.parametrize(
    "metaestimator",
    METAESTIMATORS,
    ids=[str(row["metaestimator"].__class__.__name__) for row in METAESTIMATORS],
)
def test_error_for_other_methods(metaestimator):
    # This test complements test_warning_for_indicated_methods but checks for
    # UnsetMetadataPassedError instead of FutureWarning
    cls = metaestimator["metaestimator"]
    estimator = metaestimator["estimator"]()
    estimator_name = metaestimator["estimator_name"]
    instance = cls(**{estimator_name: estimator})
    X = metaestimator["X"]
    y = metaestimator["y"]
    routing_methods = metaestimator["routing_methods"]

    for method_name in routing_methods:
        if metaestimator.get(method_name, {}).get("warns_on"):
            # this method is not expected to warn
            continue

        msg = (
            "[sample_weight, metadata] are passed but are not explicitly set as "
            "requested or not for "
            f"{estimator.__class__.__name__}.{method_name}"
        )
        with pytest.raises(UnsetMetadataPassedError, match=re.escape(msg)):
            method = getattr(instance, method_name)
            method(X, y, sample_weight=sample_weight, metadata=metadata)


@pytest.mark.parametrize(
    "metaestimator",
    METAESTIMATORS,
    ids=[str(row["metaestimator"].__class__.__name__) for row in METAESTIMATORS],
)
def test_setting_request_removes_warning_or_error(metaestimator):
    # When the metadata is explicitly requested, there should be no warning and
    # no error.
    def set_request(estimator, method_name):
        # e.g. call set_fit_request on estimator
        method = getattr(estimator, f"set_{method_name}_request")
        method(sample_weight=True, metadata=True)

    cls = metaestimator["metaestimator"]
    estimator_name = metaestimator["estimator_name"]
    X = metaestimator["X"]
    y = metaestimator["y"]
    routing_methods = metaestimator["routing_methods"]

    for method_name in routing_methods:
        estimator = metaestimator["estimator"]()
        set_request(estimator, method_name)
        instance = cls(**{estimator_name: estimator})
        # lines below to ensure that there are no warnings
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            method = getattr(instance, method_name)
            method(X, y, sample_weight=sample_weight, metadata=metadata)
