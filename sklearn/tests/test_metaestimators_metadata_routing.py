import copy
import re

import numpy as np
import pytest

from sklearn import config_context
from sklearn.calibration import CalibratedClassifierCV
from sklearn.exceptions import UnsetMetadataPassedError
from sklearn.linear_model import LogisticRegressionCV
from sklearn.multioutput import (
    ClassifierChain,
    MultiOutputClassifier,
    MultiOutputRegressor,
    RegressorChain,
)
from sklearn.tests.metadata_routing_common import (
    ConsumingClassifier,
    ConsumingRegressor,
    ConsumingScorer,
    ConsumingSplitter,
    _Registry,
    assert_request_is_empty,
    check_recorded_metadata,
)
from sklearn.utils.metadata_routing import MetadataRouter

rng = np.random.RandomState(42)
N, M = 100, 4
X = rng.rand(N, M)
y = rng.randint(0, 2, size=N)
y_multi = rng.randint(0, 2, size=(N, 3))
metadata = rng.randint(0, 10, size=N)
sample_weight = rng.rand(N)
groups = np.array([0, 1] * (len(y) // 2))


@pytest.fixture(autouse=True)
def enable_slep006():
    """Enable SLEP006 for all tests."""
    with config_context(enable_metadata_routing=True):
        yield


METAESTIMATORS: list = [
    {
        "metaestimator": MultiOutputRegressor,
        "estimator_name": "estimator",
        "estimator": ConsumingRegressor,
        "X": X,
        "y": y_multi,
        "routing_methods": ["fit", "partial_fit"],
    },
    {
        "metaestimator": MultiOutputClassifier,
        "estimator_name": "estimator",
        "estimator": ConsumingClassifier,
        "X": X,
        "y": y_multi,
        "routing_methods": ["fit", "partial_fit"],
    },
    {
        "metaestimator": CalibratedClassifierCV,
        "estimator_name": "estimator",
        "estimator": ConsumingClassifier,
        "X": X,
        "y": y,
        "routing_methods": ["fit"],
        "preserves_metadata": False,
    },
    {
        "metaestimator": ClassifierChain,
        "estimator_name": "base_estimator",
        "estimator": ConsumingClassifier,
        "X": X,
        "y": y_multi,
        "routing_methods": ["fit"],
    },
    {
        "metaestimator": RegressorChain,
        "estimator_name": "base_estimator",
        "estimator": ConsumingRegressor,
        "X": X,
        "y": y_multi,
        "routing_methods": ["fit"],
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
- preserves_metadata: Whether the metaestimator passes the metadata to the
  sub-estimator without modification or not. If it does, we check that the
  values are identical. If it doesn't, no check is performed. TODO Maybe
  something smarter could be done if the data is modified.

"""

# ids used for pytest fixture
METAESTIMATOR_IDS = [str(row["metaestimator"].__name__) for row in METAESTIMATORS]

CV_SCORERS: list = [
    {
        "cv_estimator": LogisticRegressionCV,
        "scorer_name": "scoring",
        "routing_methods": ["fit", "score"],
    },
]

CV_SPLITTERS: list = [
    {
        "cv_estimator": LogisticRegressionCV,
        "splitter_name": "cv",
        "routing_methods": ["fit"],
    }
]

# IDs used by pytest to get meaningful verbose messages when running the tests
CV_SCORER_IDS = [x["cv_estimator"].__name__ for x in CV_SCORERS]
CV_SPLITTER_IDS = [x["cv_estimator"].__name__ for x in CV_SPLITTERS]


def test_registry_copy():
    # test that _Registry is not copied into a new instance.
    a = _Registry()
    b = _Registry()
    assert a is not b
    assert a is copy.copy(a)
    assert a is copy.deepcopy(a)


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
def test_error_on_missing_requests(metaestimator):
    # Test that a UnsetMetadataPassedError is raised when it should.
    cls = metaestimator["metaestimator"]
    estimator = metaestimator["estimator"]()
    estimator_name = metaestimator["estimator_name"]
    X = metaestimator["X"]
    y = metaestimator["y"]
    routing_methods = metaestimator["routing_methods"]

    for method_name in routing_methods:
        for key in ["sample_weight", "metadata"]:
            val = {"sample_weight": sample_weight, "metadata": metadata}[key]
            kwargs = {key: val}
            msg = (
                f"[{key}] are passed but are not explicitly set as requested or not"
                f" for {estimator.__class__.__name__}.{method_name}"
            )

            instance = cls(**{estimator_name: estimator})
            if "fit" not in method_name:  # instance needs to be fitted first
                instance.fit(X, y)  # pragma: no cover
            with pytest.raises(UnsetMetadataPassedError, match=re.escape(msg)):
                method = getattr(instance, method_name)
                method(X, y, **kwargs)


@pytest.mark.parametrize(
    "metaestimator",
    METAESTIMATORS,
    ids=METAESTIMATOR_IDS,
)
def test_setting_request_removes_error(metaestimator):
    # When the metadata is explicitly requested, there should be no errors.
    def set_request(estimator, method_name):
        # e.g. call set_fit_request on estimator
        set_request_for_method = getattr(estimator, f"set_{method_name}_request")
        set_request_for_method(sample_weight=True, metadata=True)

    cls = metaestimator["metaestimator"]
    estimator_name = metaestimator["estimator_name"]
    X = metaestimator["X"]
    y = metaestimator["y"]
    routing_methods = metaestimator["routing_methods"]
    preserves_metadata = metaestimator.get("preserves_metadata", True)

    for method_name in routing_methods:
        for key in ["sample_weight", "metadata"]:
            val = {"sample_weight": sample_weight, "metadata": metadata}[key]
            kwargs = {key: val}

            registry = _Registry()
            estimator = metaestimator["estimator"](registry=registry)
            set_request(estimator, method_name)
            instance = cls(**{estimator_name: estimator})
            method = getattr(instance, method_name)
            method(X, y, **kwargs)

            if preserves_metadata:
                # sanity check that registry is not empty, or else the test
                # passes trivially
                assert registry
                for estimator in registry:
                    check_recorded_metadata(estimator, method_name, **kwargs)


@pytest.mark.parametrize("cv_scorer", CV_SCORERS, ids=CV_SCORER_IDS)
def test_metadata_is_routed_correctly_to_scorer(cv_scorer):
    """Test that any requested metadata is correctly routed to the underlying
    scorers in CV estimators.
    """
    registry = _Registry()
    cls = cv_scorer["cv_estimator"]
    scorer_name = cv_scorer["scorer_name"]
    scorer = ConsumingScorer(registry=registry)
    scorer.set_score_request(sample_weight=True)
    routing_methods = cv_scorer["routing_methods"]

    for method_name in routing_methods:
        instance = cls(**{scorer_name: scorer})
        method = getattr(instance, method_name)
        kwargs = {"sample_weight": sample_weight}
        if "fit" not in method_name:  # instance needs to be fitted first
            instance.fit(X, y)
        method(X, y, **kwargs)
        for _scorer in registry:
            check_recorded_metadata(
                obj=_scorer,
                method="score",
                split_params=("sample_weight",),
                **kwargs,
            )


@pytest.mark.parametrize("cv_splitter", CV_SPLITTERS, ids=CV_SPLITTER_IDS)
def test_metadata_is_routed_correctly_to_splitter(cv_splitter):
    """Test that any requested metadata is correctly routed to the underlying
    splitters in CV estimators.
    """
    registry = _Registry()
    cls = cv_splitter["cv_estimator"]
    splitter_name = cv_splitter["splitter_name"]
    splitter = ConsumingSplitter(registry=registry)
    routing_methods = cv_splitter["routing_methods"]

    for method_name in routing_methods:
        instance = cls(**{splitter_name: splitter})
        method = getattr(instance, method_name)
        kwargs = {"groups": groups}
        method(X, y, **kwargs)
        for _splitter in registry:
            check_recorded_metadata(obj=_splitter, method="split", **kwargs)
