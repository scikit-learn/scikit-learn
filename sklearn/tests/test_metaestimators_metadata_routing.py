import re
import warnings
from functools import partial

import numpy as np
import pytest
from sklearn.base import RegressorMixin, ClassifierMixin, BaseEstimator
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import BaggingClassifier, BaggingRegressor
from sklearn.exceptions import UnsetMetadataPassedError
from sklearn.multioutput import (
    MultiOutputRegressor,
    MultiOutputClassifier,
    ClassifierChain,
    RegressorChain,
)
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


record_metadata_not_default = partial(record_metadata, record_default=False)


class _Registry(list):
    # This list is used to get a reference to the sub-estimators, which are not
    # necessarily stored on the metaestimator. We need to override __deepcopy__
    # because the sub-estimators are probably cloned, which would result in a
    # new copy of the list.
    def __deepcopy__(self, memo):
        return self

    def __copy__(self):
        return self


class ConsumingRegressor(RegressorMixin, BaseEstimator):
    """A regressor consuming metadata.

    Parameters
    ----------
    registry : list, default=None
        If a list, the estimator will append itself to the list in order to have
        a reference to the estimator later on. Since that reference is not
        required in all tests, registration can be skipped by leaving this value
        as None.

    """

    def __init__(self, registry=None):
        self.registry = registry

    def partial_fit(self, X, y, sample_weight="default", metadata="default"):
        if self.registry is not None:
            self.registry.append(self)

        record_metadata_not_default(
            self, "partial_fit", sample_weight=sample_weight, metadata=metadata
        )
        return self

    def fit(self, X, y, sample_weight="default", metadata="default"):
        if self.registry is not None:
            self.registry.append(self)

        record_metadata_not_default(
            self, "fit", sample_weight=sample_weight, metadata=metadata
        )
        return self

    def predict(self, X, sample_weight="default", metadata="default"):
        if self.registry is not None:
            self.registry.append(self)

        record_metadata_not_default(
            self, "predict", sample_weight=sample_weight, metadata=metadata
        )
        return np.zeros(shape=(len(X),))


class ConsumingClassifier(ClassifierMixin, BaseEstimator):
    """A classifier consuming metadata.

    Parameters
    ----------
    registry : list, default=None
        If a list, the estimator will append itself to the list in order to have
        a reference to the estimator later on. Since that reference is not
        required in all tests, registration can be skipped by leaving this value
        as None.

    """

    def __init__(self, registry=None):
        self.registry = registry

    def partial_fit(self, X, y, sample_weight="default", metadata="default"):
        if self.registry is not None:
            self.registry.append(self)

        record_metadata_not_default(
            self, "partial_fit", sample_weight=sample_weight, metadata=metadata
        )
        self.classes_ = [0, 1]
        return self

    def fit(self, X, y, sample_weight="default", metadata="default"):
        if self.registry is not None:
            self.registry.append(self)

        record_metadata_not_default(
            self, "fit", sample_weight=sample_weight, metadata=metadata
        )
        self.classes_ = [0, 1]
        return self

    def predict(self, X, sample_weight="default", metadata="default"):
        if self.registry is not None:
            self.registry.append(self)

        record_metadata_not_default(
            self, "predict", sample_weight=sample_weight, metadata=metadata
        )
        return np.zeros(shape=(len(X),))

    def predict_proba(self, X, sample_weight="default", metadata="default"):
        if self.registry is not None:
            self.registry.append(self)

        record_metadata_not_default(
            self, "predict_proba", sample_weight=sample_weight, metadata=metadata
        )
        return np.asarray([[0.0, 1.0]] * len(X))

    def predict_log_proba(self, X, sample_weight="default", metadata="default"):
        if self.registry is not None:
            self.registry.append(self)

        record_metadata_not_default(
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
        "warns_on": {
            "fit": ["sample_weight", "metadata"],
            "partial_fit": ["sample_weight"],
        },
    },
    {
        "metaestimator": MultiOutputClassifier,
        "estimator_name": "estimator",
        "estimator": ConsumingClassifier,
        "X": X,
        "y": y_multi,
        "routing_methods": ["fit", "partial_fit"],
        "warns_on": {
            "fit": ["sample_weight", "metadata"],
            "partial_fit": ["sample_weight"],
        },
    },
    {
        "metaestimator": CalibratedClassifierCV,
        "estimator_name": "estimator",
        "estimator": ConsumingClassifier,
        "X": X,
        "y": y,
        "routing_methods": ["fit"],
        "warns_on": {"fit": ["sample_weight", "metadata"]},
        "preserves_metadata": False,  # applies CV splits
    },
    {
        "metaestimator": BaggingRegressor,
        "estimator_name": "base_estimator",
        "estimator": ConsumingRegressor,
        "X": X,
        "y": y,
        "routing_methods": ["fit"],
        "warns_on": {
            "fit": ["sample_weight"],
        },
        "preserves_metadata": False,  # applies sub-sampling
    },
    {
        "metaestimator": BaggingClassifier,
        "estimator_name": "base_estimator",
        "estimator": ConsumingClassifier,
        "X": X,
        "y": y,
        "routing_methods": ["fit"],
        "warns_on": {
            "fit": ["sample_weight"],
        },
        "preserves_metadata": False,  # applies sub-sampling
    },
    {
        "metaestimator": ClassifierChain,
        "estimator_name": "base_estimator",
        "estimator": ConsumingClassifier,
        "X": X,
        "y": y_multi,
        "routing_methods": ["fit"],
        "warns_on": {},
    },
    {
        "metaestimator": RegressorChain,
        "estimator_name": "base_estimator",
        "estimator": ConsumingRegressor,
        "X": X,
        "y": y_multi,
        "routing_methods": ["fit"],
        "warns_on": {"fit": ["sample_weight", "metadata"]},
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
- warns_on: A dict containing all methods as keys, and arguments as values,
  whose combination is supposed to result in a warning if routing is not
  requested. It is implied that all routing methods and arguments not listed
  here should result in an error.
- preserves_metadata: Whether the metaestimator passes the metadata to the
  sub-estimator without modification or not. If it does, we check that the
  values are identical. If it doesn', no check is performed. TODO Maybe
  something smarter could be done if the data is modified.

"""

# ids used for pytest fixture
METAESTIMATOR_IDS = [str(row["metaestimator"].__name__) for row in METAESTIMATORS]


@pytest.mark.parametrize(
    "metaestimator",
    METAESTIMATORS,
    ids=METAESTIMATOR_IDS,
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
    ids=METAESTIMATOR_IDS,
)
def test_warning_for_indicated_methods(metaestimator):
    # Check that the indicated methods give a warning
    # TODO: Always error for 1.4
    cls = metaestimator["metaestimator"]
    registry = _Registry()
    estimator = metaestimator["estimator"](registry=registry)
    estimator_name = metaestimator["estimator_name"]
    X = metaestimator["X"]
    y = metaestimator["y"]
    routing_methods = metaestimator["routing_methods"]
    warns_on = metaestimator["warns_on"]

    for method_name in routing_methods:
        if method_name not in warns_on:
            # this method is not expected to warn
            continue

        for key in warns_on[method_name]:
            val = {"sample_weight": sample_weight, "metadata": metadata}[key]
            kwargs = {key: val}
            warn_msg = (
                "You are passing metadata for which the request values are not"
                f" explicitly set: {key}. From version 1.4 this results in the"
                f" following error: [{key}] are passed but are not explicitly set as"
                f" requested or not for {estimator.__class__.__name__}.{method_name}"
            )

            instance = cls(**{estimator_name: estimator})
            if "fit" not in method_name:  # instance needs to be fitted first
                instance.fit(X, y)
            with pytest.warns(FutureWarning, match=re.escape(warn_msg)):
                method = getattr(instance, method_name)
                method(X, y, **kwargs)

            if metaestimator.get("preserves_metadata", True):
                # sanity check that registry is not empty, or else the test
                # passes trivially
                assert registry
                for estimator in registry:
                    check_recorded_metadata(estimator, method_name, **kwargs)
            # clear the registry since the check could be different for the next
            # method being tested
            registry.clear()


@pytest.mark.parametrize(
    "metaestimator",
    METAESTIMATORS,
    ids=METAESTIMATOR_IDS,
)
def test_error_for_other_methods(metaestimator):
    # This test complements test_warning_for_indicated_methods but checks for
    # UnsetMetadataPassedError instead of FutureWarning
    cls = metaestimator["metaestimator"]
    estimator = metaestimator["estimator"]()
    estimator_name = metaestimator["estimator_name"]
    X = metaestimator["X"]
    y = metaestimator["y"]
    routing_methods = metaestimator["routing_methods"]
    warns_on = metaestimator["warns_on"]

    for method_name in routing_methods:
        warn_args = warns_on.get(method_name, [])
        for key in ["sample_weight", "metadata"]:
            if key in warn_args:
                # this method is expected to warn for this argument, not raise
                continue

            val = {"sample_weight": sample_weight, "metadata": metadata}[key]
            kwargs = {key: val}
            msg = (
                f"[{key}] are passed but are not explicitly set as requested or not for"
                f" {estimator.__class__.__name__}.{method_name}"
            )

            instance = cls(**{estimator_name: estimator})
            if "fit" not in method_name:  # instance needs to be fitted first
                instance.fit(X, y)
            with pytest.raises(UnsetMetadataPassedError, match=re.escape(msg)):
                method = getattr(instance, method_name)
                method(X, y, **kwargs)


@pytest.mark.parametrize(
    "metaestimator",
    METAESTIMATORS,
    ids=METAESTIMATOR_IDS,
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
        with warnings.catch_warnings(record=True) as rec:
            method = getattr(instance, method_name)
            method(X, y, sample_weight=sample_weight, metadata=metadata)
            # Check that there was no FutureWarning about metadata. The exact
            # error message is not checked on purpose, because if the message is
            # changed without amending this test, the test would pass trivially.
            future_warnings = [w for w in rec if isinstance(w, FutureWarning)]
            assert not any("metadata" in w.message for w in future_warnings)
