import copy
import re

import numpy as np
import pytest

from sklearn import config_context
from sklearn.base import BaseEstimator, is_classifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.compose import TransformedTargetRegressor
from sklearn.covariance import GraphicalLassoCV
from sklearn.ensemble import (
    AdaBoostClassifier,
    AdaBoostRegressor,
    BaggingClassifier,
    BaggingRegressor,
)
from sklearn.exceptions import UnsetMetadataPassedError
from sklearn.experimental import (
    enable_halving_search_cv,  # noqa: F401
    enable_iterative_imputer,  # noqa: F401
)
from sklearn.feature_selection import (
    RFE,
    RFECV,
    SelectFromModel,
    SequentialFeatureSelector,
)
from sklearn.impute import IterativeImputer
from sklearn.linear_model import (
    ElasticNetCV,
    LarsCV,
    LassoCV,
    LassoLarsCV,
    LogisticRegressionCV,
    MultiTaskElasticNetCV,
    MultiTaskLassoCV,
    OrthogonalMatchingPursuitCV,
    RANSACRegressor,
    RidgeClassifierCV,
    RidgeCV,
)
from sklearn.metrics._regression import mean_squared_error
from sklearn.metrics._scorer import make_scorer
from sklearn.model_selection import (
    FixedThresholdClassifier,
    GridSearchCV,
    GroupKFold,
    HalvingGridSearchCV,
    HalvingRandomSearchCV,
    RandomizedSearchCV,
    TunedThresholdClassifierCV,
    cross_validate,
)
from sklearn.multiclass import (
    OneVsOneClassifier,
    OneVsRestClassifier,
    OutputCodeClassifier,
)
from sklearn.multioutput import (
    ClassifierChain,
    MultiOutputClassifier,
    MultiOutputRegressor,
    RegressorChain,
)
from sklearn.semi_supervised import SelfTrainingClassifier
from sklearn.tests.metadata_routing_common import (
    ConsumingClassifier,
    ConsumingRegressor,
    ConsumingScorer,
    ConsumingSplitter,
    NonConsumingClassifier,
    NonConsumingRegressor,
    _Registry,
    assert_request_is_empty,
    check_recorded_metadata,
)
from sklearn.utils.metadata_routing import MetadataRouter

rng = np.random.RandomState(42)
N, M = 100, 4
X = rng.rand(N, M)
y = rng.randint(0, 3, size=N)
y_binary = (y >= 1).astype(int)
classes = np.unique(y)
y_multi = rng.randint(0, 3, size=(N, 3))
classes_multi = [np.unique(y_multi[:, i]) for i in range(y_multi.shape[1])]
metadata = rng.randint(0, 10, size=N)
sample_weight = rng.rand(N)
groups = rng.randint(0, 10, size=len(y))


METAESTIMATORS: list = [
    {
        "metaestimator": MultiOutputRegressor,
        "estimator_name": "estimator",
        "estimator": "regressor",
        "X": X,
        "y": y_multi,
        "estimator_routing_methods": ["fit", "partial_fit"],
    },
    {
        "metaestimator": MultiOutputClassifier,
        "estimator_name": "estimator",
        "estimator": "classifier",
        "X": X,
        "y": y_multi,
        "estimator_routing_methods": ["fit", "partial_fit"],
        "method_args": {"partial_fit": {"classes": classes_multi}},
    },
    {
        "metaestimator": CalibratedClassifierCV,
        "estimator_name": "estimator",
        "estimator": "classifier",
        "X": X,
        "y": y,
        "estimator_routing_methods": ["fit"],
        "preserves_metadata": "subset",
    },
    {
        "metaestimator": ClassifierChain,
        "estimator_name": "estimator",
        "estimator": "classifier",
        "X": X,
        "y": y_multi,
        "estimator_routing_methods": ["fit"],
    },
    {
        "metaestimator": RegressorChain,
        "estimator_name": "estimator",
        "estimator": "regressor",
        "X": X,
        "y": y_multi,
        "estimator_routing_methods": ["fit"],
    },
    {
        "metaestimator": LogisticRegressionCV,
        "X": X,
        "y": y,
        "scorer_name": "scoring",
        "scorer_routing_methods": ["fit", "score"],
        "cv_name": "cv",
        "cv_routing_methods": ["fit"],
    },
    {
        "metaestimator": GridSearchCV,
        "estimator_name": "estimator",
        "estimator": "classifier",
        "init_args": {"param_grid": {"alpha": [0.1, 0.2]}},
        "X": X,
        "y": y,
        "estimator_routing_methods": ["fit"],
        "preserves_metadata": "subset",
        "scorer_name": "scoring",
        "scorer_routing_methods": ["fit", "score"],
        "cv_name": "cv",
        "cv_routing_methods": ["fit"],
    },
    {
        "metaestimator": RandomizedSearchCV,
        "estimator_name": "estimator",
        "estimator": "classifier",
        "init_args": {"param_distributions": {"alpha": [0.1, 0.2]}},
        "X": X,
        "y": y,
        "estimator_routing_methods": ["fit"],
        "preserves_metadata": "subset",
        "scorer_name": "scoring",
        "scorer_routing_methods": ["fit", "score"],
        "cv_name": "cv",
        "cv_routing_methods": ["fit"],
    },
    {
        "metaestimator": HalvingGridSearchCV,
        "estimator_name": "estimator",
        "estimator": "classifier",
        "init_args": {"param_grid": {"alpha": [0.1, 0.2]}},
        "X": X,
        "y": y,
        "estimator_routing_methods": ["fit"],
        "preserves_metadata": "subset",
        "scorer_name": "scoring",
        "scorer_routing_methods": ["fit", "score"],
        "cv_name": "cv",
        "cv_routing_methods": ["fit"],
    },
    {
        "metaestimator": HalvingRandomSearchCV,
        "estimator_name": "estimator",
        "estimator": "classifier",
        "init_args": {"param_distributions": {"alpha": [0.1, 0.2]}},
        "X": X,
        "y": y,
        "estimator_routing_methods": ["fit"],
        "preserves_metadata": "subset",
        "scorer_name": "scoring",
        "scorer_routing_methods": ["fit", "score"],
        "cv_name": "cv",
        "cv_routing_methods": ["fit"],
    },
    {
        "metaestimator": FixedThresholdClassifier,
        "estimator_name": "estimator",
        "estimator": "classifier",
        "X": X,
        "y": y_binary,
        "estimator_routing_methods": ["fit"],
        "preserves_metadata": "subset",
    },
    {
        "metaestimator": TunedThresholdClassifierCV,
        "estimator_name": "estimator",
        "estimator": "classifier",
        "X": X,
        "y": y_binary,
        "estimator_routing_methods": ["fit"],
        "preserves_metadata": "subset",
    },
    {
        "metaestimator": OneVsRestClassifier,
        "estimator_name": "estimator",
        "estimator": "classifier",
        "X": X,
        "y": y,
        "estimator_routing_methods": ["fit", "partial_fit"],
        "method_args": {"partial_fit": {"classes": classes}},
    },
    {
        "metaestimator": OneVsOneClassifier,
        "estimator_name": "estimator",
        "estimator": "classifier",
        "X": X,
        "y": y,
        "estimator_routing_methods": ["fit", "partial_fit"],
        "preserves_metadata": "subset",
        "method_args": {"partial_fit": {"classes": classes}},
    },
    {
        "metaestimator": OutputCodeClassifier,
        "estimator_name": "estimator",
        "estimator": "classifier",
        "init_args": {"random_state": 42},
        "X": X,
        "y": y,
        "estimator_routing_methods": ["fit"],
    },
    {
        "metaestimator": SelectFromModel,
        "estimator_name": "estimator",
        "estimator": "classifier",
        "X": X,
        "y": y,
        "estimator_routing_methods": ["fit", "partial_fit"],
        "method_args": {"partial_fit": {"classes": classes}},
    },
    {
        "metaestimator": OrthogonalMatchingPursuitCV,
        "X": X,
        "y": y,
        "cv_name": "cv",
        "cv_routing_methods": ["fit"],
    },
    {
        "metaestimator": ElasticNetCV,
        "X": X,
        "y": y,
        "cv_name": "cv",
        "cv_routing_methods": ["fit"],
    },
    {
        "metaestimator": LassoCV,
        "X": X,
        "y": y,
        "cv_name": "cv",
        "cv_routing_methods": ["fit"],
    },
    {
        "metaestimator": MultiTaskElasticNetCV,
        "X": X,
        "y": y_multi,
        "cv_name": "cv",
        "cv_routing_methods": ["fit"],
    },
    {
        "metaestimator": MultiTaskLassoCV,
        "X": X,
        "y": y_multi,
        "cv_name": "cv",
        "cv_routing_methods": ["fit"],
    },
    {
        "metaestimator": LarsCV,
        "X": X,
        "y": y,
        "cv_name": "cv",
        "cv_routing_methods": ["fit"],
    },
    {
        "metaestimator": LassoLarsCV,
        "X": X,
        "y": y,
        "cv_name": "cv",
        "cv_routing_methods": ["fit"],
    },
    {
        "metaestimator": RANSACRegressor,
        "estimator_name": "estimator",
        "estimator": "regressor",
        "init_args": {"min_samples": 0.5},
        "X": X,
        "y": y,
        "preserves_metadata": "subset",
        "estimator_routing_methods": ["fit", "predict", "score"],
        "method_mapping": {"fit": ["fit", "score"]},
    },
    {
        "metaestimator": IterativeImputer,
        "estimator_name": "estimator",
        "estimator": "regressor",
        "init_args": {"skip_complete": False},
        "X": X,
        "y": y,
        "estimator_routing_methods": ["fit"],
    },
    {
        "metaestimator": BaggingClassifier,
        "estimator_name": "estimator",
        "estimator": "classifier",
        "X": X,
        "y": y,
        "preserves_metadata": False,
        "estimator_routing_methods": [
            "fit",
            "predict",
            "predict_proba",
            "predict_log_proba",
            "decision_function",
        ],
        "method_mapping": {
            "predict": ["predict", "predict_proba"],
            "predict_proba": ["predict", "predict_proba"],
            "predict_log_proba": ["predict", "predict_proba", "predict_log_proba"],
        },
    },
    {
        "metaestimator": BaggingRegressor,
        "estimator_name": "estimator",
        "estimator": "regressor",
        "X": X,
        "y": y,
        "preserves_metadata": False,
        "estimator_routing_methods": ["fit", "predict"],
    },
    {
        "metaestimator": RidgeCV,
        "X": X,
        "y": y,
        "scorer_name": "scoring",
        "scorer_routing_methods": ["fit"],
    },
    {
        "metaestimator": RidgeClassifierCV,
        "X": X,
        "y": y,
        "scorer_name": "scoring",
        "scorer_routing_methods": ["fit"],
    },
    {
        "metaestimator": RidgeCV,
        "X": X,
        "y": y,
        "scorer_name": "scoring",
        "scorer_routing_methods": ["fit"],
        "cv_name": "cv",
        "cv_routing_methods": ["fit"],
    },
    {
        "metaestimator": RidgeClassifierCV,
        "X": X,
        "y": y,
        "scorer_name": "scoring",
        "scorer_routing_methods": ["fit"],
        "cv_name": "cv",
        "cv_routing_methods": ["fit"],
    },
    {
        "metaestimator": GraphicalLassoCV,
        "X": X,
        "y": y,
        "cv_name": "cv",
        "cv_routing_methods": ["fit"],
    },
    {
        "metaestimator": TransformedTargetRegressor,
        "estimator": "regressor",
        "estimator_name": "regressor",
        "X": X,
        "y": y,
        "estimator_routing_methods": ["fit", "predict"],
    },
    {
        "metaestimator": SelfTrainingClassifier,
        "estimator_name": "estimator",
        "estimator": "classifier",
        "X": X,
        "y": y,
        "preserves_metadata": True,
        "estimator_routing_methods": [
            "fit",
            "predict",
            "predict_proba",
            "predict_log_proba",
            "decision_function",
            "score",
        ],
        "method_mapping": {"fit": ["fit", "score"]},
    },
    {
        "metaestimator": SequentialFeatureSelector,
        "estimator_name": "estimator",
        "estimator": "classifier",
        "X": X,
        "y": y,
        "estimator_routing_methods": ["fit"],
        "scorer_name": "scoring",
        "scorer_routing_methods": ["fit"],
        "cv_name": "cv",
        "cv_routing_methods": ["fit"],
    },
    {
        "metaestimator": RFE,
        "estimator": "classifier",
        "estimator_name": "estimator",
        "X": X,
        "y": y,
        "estimator_routing_methods": ["fit", "predict", "score"],
    },
    {
        "metaestimator": RFECV,
        "estimator": "classifier",
        "estimator_name": "estimator",
        "estimator_routing_methods": ["fit"],
        "cv_name": "cv",
        "cv_routing_methods": ["fit"],
        "scorer_name": "scoring",
        "scorer_routing_methods": ["fit", "score"],
        "X": X,
        "y": y,
    },
]
"""List containing all metaestimators to be tested and their settings

The keys are as follows:

- metaestimator: The metaestimator to be tested
- estimator_name: The name of the argument for the sub-estimator
- estimator: The sub-estimator type, either "regressor" or "classifier"
- init_args: The arguments to be passed to the metaestimator's constructor
- X: X-data to fit and predict
- y: y-data to fit
- estimator_routing_methods: list of all methods to check for routing metadata
  to the sub-estimator
- preserves_metadata:
    - True (default): the metaestimator passes the metadata to the
      sub-estimator without modification. We check that the values recorded by
      the sub-estimator are identical to what we've passed to the
      metaestimator.
    - False: no check is performed regarding values, we only check that a
      metadata with the expected names/keys are passed.
    - "subset": we check that the recorded metadata by the sub-estimator is a
      subset of what is passed to the metaestimator.
- scorer_name: The name of the argument for the scorer
- scorer_routing_methods: list of all methods to check for routing metadata
  to the scorer
- cv_name: The name of the argument for the CV splitter
- cv_routing_methods: list of all methods to check for routing metadata
  to the splitter
- method_args: a dict of dicts, defining extra arguments needed to be passed to
  methods, such as passing `classes` to `partial_fit`.
- method_mapping: a dict of the form `{caller: [callee1, ...]}` which signals
  which `.set_{method}_request` methods should be called to set request values.
  If not present, a one-to-one mapping is assumed.
"""

# IDs used by pytest to get meaningful verbose messages when running the tests
METAESTIMATOR_IDS = [str(row["metaestimator"].__name__) for row in METAESTIMATORS]

UNSUPPORTED_ESTIMATORS = [
    AdaBoostClassifier(),
    AdaBoostRegressor(),
]


def get_init_args(metaestimator_info, sub_estimator_consumes):
    """Get the init args for a metaestimator

    This is a helper function to get the init args for a metaestimator from
    the METAESTIMATORS list. It returns an empty dict if no init args are
    required.

    Parameters
    ----------
    metaestimator_info : dict
        The metaestimator info from METAESTIMATORS

    sub_estimator_consumes : bool
        Whether the sub-estimator consumes metadata or not.

    Returns
    -------
    kwargs : dict
        The init args for the metaestimator.

    (estimator, estimator_registry) : (estimator, registry)
        The sub-estimator and the corresponding registry.

    (scorer, scorer_registry) : (scorer, registry)
        The scorer and the corresponding registry.

    (cv, cv_registry) : (CV splitter, registry)
        The CV splitter and the corresponding registry.
    """
    kwargs = metaestimator_info.get("init_args", {})
    estimator, estimator_registry = None, None
    scorer, scorer_registry = None, None
    cv, cv_registry = None, None
    if "estimator" in metaestimator_info:
        estimator_name = metaestimator_info["estimator_name"]
        estimator_registry = _Registry()
        sub_estimator_type = metaestimator_info["estimator"]
        if sub_estimator_consumes:
            if sub_estimator_type == "regressor":
                estimator = ConsumingRegressor(estimator_registry)
            elif sub_estimator_type == "classifier":
                estimator = ConsumingClassifier(estimator_registry)
            else:
                raise ValueError("Unpermitted `sub_estimator_type`.")  # pragma: nocover
        else:
            if sub_estimator_type == "regressor":
                estimator = NonConsumingRegressor()
            elif sub_estimator_type == "classifier":
                estimator = NonConsumingClassifier()
            else:
                raise ValueError("Unpermitted `sub_estimator_type`.")  # pragma: nocover
        kwargs[estimator_name] = estimator
    if "scorer_name" in metaestimator_info:
        scorer_name = metaestimator_info["scorer_name"]
        scorer_registry = _Registry()
        scorer = ConsumingScorer(registry=scorer_registry)
        kwargs[scorer_name] = scorer
    if "cv_name" in metaestimator_info:
        cv_name = metaestimator_info["cv_name"]
        cv_registry = _Registry()
        cv = ConsumingSplitter(registry=cv_registry)
        kwargs[cv_name] = cv

    return (
        kwargs,
        (estimator, estimator_registry),
        (scorer, scorer_registry),
        (cv, cv_registry),
    )


def set_requests(obj, *, method_mapping, methods, metadata_name, value=True):
    """Call `set_{method}_request` on a list of methods from the sub-estimator.

    Parameters
    ----------
    obj : BaseEstimator
        The object for which `set_{method}_request` methods are called.

    method_mapping : dict
        The method mapping in the form of `{caller: [callee, ...]}`.
        If a "caller" is not present in the method mapping, a one-to-one mapping is
        assumed.

    methods : list of str
        The list of methods as "caller"s for which the request for the child should
        be set.

    metadata_name : str
        The name of the metadata to be routed, usually either `"metadata"` or
        `"sample_weight"` in our tests.

    value : None, bool, or str
        The request value to be set, by default it's `True`
    """
    for caller in methods:
        for callee in method_mapping.get(caller, [caller]):
            set_request_for_method = getattr(obj, f"set_{callee}_request")
            set_request_for_method(**{metadata_name: value})
            if (
                isinstance(obj, BaseEstimator)
                and is_classifier(obj)
                and callee == "partial_fit"
            ):
                set_request_for_method(classes=True)


@pytest.mark.parametrize("estimator", UNSUPPORTED_ESTIMATORS)
@config_context(enable_metadata_routing=True)
def test_unsupported_estimators_get_metadata_routing(estimator):
    """Test that get_metadata_routing is not implemented on meta-estimators for
    which we haven't implemented routing yet."""
    with pytest.raises(NotImplementedError):
        estimator.get_metadata_routing()


@pytest.mark.parametrize("estimator", UNSUPPORTED_ESTIMATORS)
@config_context(enable_metadata_routing=True)
def test_unsupported_estimators_fit_with_metadata(estimator):
    """Test that fit raises NotImplementedError when metadata routing is
    enabled and a metadata is passed on meta-estimators for which we haven't
    implemented routing yet."""
    with pytest.raises(NotImplementedError):
        try:
            estimator.fit([[1]], [1], sample_weight=[1])
        except TypeError:
            # not all meta-estimators in the list support sample_weight,
            # and for those we skip this test.
            raise NotImplementedError


@config_context(enable_metadata_routing=True)
def test_registry_copy():
    # test that _Registry is not copied into a new instance.
    a = _Registry()
    b = _Registry()
    assert a is not b
    assert a is copy.copy(a)
    assert a is copy.deepcopy(a)


@pytest.mark.parametrize("metaestimator", METAESTIMATORS, ids=METAESTIMATOR_IDS)
@config_context(enable_metadata_routing=True)
def test_default_request(metaestimator):
    # Check that by default request is empty and the right type
    metaestimator_class = metaestimator["metaestimator"]
    kwargs, *_ = get_init_args(metaestimator, sub_estimator_consumes=True)
    instance = metaestimator_class(**kwargs)
    if "cv_name" in metaestimator:
        # Our GroupCV splitters request groups by default, which we should
        # ignore in this test.
        exclude = {"splitter": ["split"]}
    else:
        exclude = None
    assert_request_is_empty(instance.get_metadata_routing(), exclude=exclude)
    assert isinstance(instance.get_metadata_routing(), MetadataRouter)


@pytest.mark.parametrize("metaestimator", METAESTIMATORS, ids=METAESTIMATOR_IDS)
@config_context(enable_metadata_routing=True)
def test_error_on_missing_requests_for_sub_estimator(metaestimator):
    # Test that a UnsetMetadataPassedError is raised when the sub-estimator's
    # requests are not set
    if "estimator" not in metaestimator:
        # This test only makes sense for metaestimators which have a
        # sub-estimator, e.g. MyMetaEstimator(estimator=MySubEstimator())
        return

    metaestimator_class = metaestimator["metaestimator"]
    X = metaestimator["X"]
    y = metaestimator["y"]
    routing_methods = metaestimator["estimator_routing_methods"]

    for method_name in routing_methods:
        for key in ["sample_weight", "metadata"]:
            kwargs, (estimator, _), (scorer, _), *_ = get_init_args(
                metaestimator, sub_estimator_consumes=True
            )
            if scorer:
                scorer.set_score_request(**{key: True})
            val = {"sample_weight": sample_weight, "metadata": metadata}[key]
            method_kwargs = {key: val}
            instance = metaestimator_class(**kwargs)
            msg = (
                f"[{key}] are passed but are not explicitly set as requested or not"
                f" requested for {estimator.__class__.__name__}.{method_name}"
            )
            with pytest.raises(UnsetMetadataPassedError, match=re.escape(msg)):
                method = getattr(instance, method_name)
                if "fit" not in method_name:
                    # set request on fit
                    set_requests(
                        estimator,
                        method_mapping=metaestimator.get("method_mapping", {}),
                        methods=["fit"],
                        metadata_name=key,
                    )
                    instance.fit(X, y, **method_kwargs)
                # making sure the requests are unset, in case they were set as a
                # side effect of setting them for fit. For instance, if method
                # mapping for fit is: `"fit": ["fit", "score"]`, that would mean
                # calling `.score` here would not raise, because we have already
                # set request value for child estimator's `score`.
                set_requests(
                    estimator,
                    method_mapping=metaestimator.get("method_mapping", {}),
                    methods=["fit"],
                    metadata_name=key,
                    value=None,
                )
                try:
                    # `fit`, `partial_fit`, 'score' accept y, others don't.
                    method(X, y, **method_kwargs)
                except TypeError:
                    method(X, **method_kwargs)


@pytest.mark.parametrize("metaestimator", METAESTIMATORS, ids=METAESTIMATOR_IDS)
@config_context(enable_metadata_routing=True)
def test_setting_request_on_sub_estimator_removes_error(metaestimator):
    # When the metadata is explicitly requested on the sub-estimator, there
    # should be no errors.
    if "estimator" not in metaestimator:
        # This test only makes sense for metaestimators which have a
        # sub-estimator, e.g. MyMetaEstimator(estimator=MySubEstimator())
        return

    metaestimator_class = metaestimator["metaestimator"]
    X = metaestimator["X"]
    y = metaestimator["y"]
    routing_methods = metaestimator["estimator_routing_methods"]
    method_mapping = metaestimator.get("method_mapping", {})
    preserves_metadata = metaestimator.get("preserves_metadata", True)

    for method_name in routing_methods:
        for key in ["sample_weight", "metadata"]:
            val = {"sample_weight": sample_weight, "metadata": metadata}[key]
            method_kwargs = {key: val}

            kwargs, (estimator, registry), (scorer, _), (cv, _) = get_init_args(
                metaestimator, sub_estimator_consumes=True
            )
            if scorer:
                set_requests(
                    scorer, method_mapping={}, methods=["score"], metadata_name=key
                )
            if cv:
                cv.set_split_request(groups=True, metadata=True)

            # `set_{method}_request({metadata}==True)` on the underlying objects
            set_requests(
                estimator,
                method_mapping=method_mapping,
                methods=[method_name],
                metadata_name=key,
            )

            instance = metaestimator_class(**kwargs)
            method = getattr(instance, method_name)
            extra_method_args = metaestimator.get("method_args", {}).get(
                method_name, {}
            )
            if "fit" not in method_name:
                # fit before calling method
                instance.fit(X, y)
            try:
                # `fit` and `partial_fit` accept y, others don't.
                method(X, y, **method_kwargs, **extra_method_args)
            except TypeError:
                method(X, **method_kwargs, **extra_method_args)

            # sanity check that registry is not empty, or else the test passes
            # trivially
            assert registry
            split_params = (
                method_kwargs.keys() if preserves_metadata == "subset" else ()
            )
            for estimator in registry:
                check_recorded_metadata(
                    estimator,
                    method=method_name,
                    parent=method_name,
                    split_params=split_params,
                    **method_kwargs,
                )


@pytest.mark.parametrize("metaestimator", METAESTIMATORS, ids=METAESTIMATOR_IDS)
@config_context(enable_metadata_routing=True)
def test_non_consuming_estimator_works(metaestimator):
    # Test that when a non-consuming estimator is given, the meta-estimator
    # works w/o setting any requests.
    # Regression test for https://github.com/scikit-learn/scikit-learn/issues/28239
    if "estimator" not in metaestimator:
        # This test only makes sense for metaestimators which have a
        # sub-estimator, e.g. MyMetaEstimator(estimator=MySubEstimator())
        return

    def set_request(estimator, method_name):
        # e.g. call set_fit_request on estimator
        if is_classifier(estimator) and method_name == "partial_fit":
            estimator.set_partial_fit_request(classes=True)

    metaestimator_class = metaestimator["metaestimator"]
    X = metaestimator["X"]
    y = metaestimator["y"]
    routing_methods = metaestimator["estimator_routing_methods"]

    for method_name in routing_methods:
        kwargs, (estimator, _), (_, _), (_, _) = get_init_args(
            metaestimator, sub_estimator_consumes=False
        )
        instance = metaestimator_class(**kwargs)
        set_request(estimator, method_name)
        method = getattr(instance, method_name)
        extra_method_args = metaestimator.get("method_args", {}).get(method_name, {})
        if "fit" not in method_name:
            instance.fit(X, y, **extra_method_args)
        # The following should pass w/o raising a routing error.
        try:
            # `fit` and `partial_fit` accept y, others don't.
            method(X, y, **extra_method_args)
        except TypeError:
            method(X, **extra_method_args)


@pytest.mark.parametrize("metaestimator", METAESTIMATORS, ids=METAESTIMATOR_IDS)
@config_context(enable_metadata_routing=True)
def test_metadata_is_routed_correctly_to_scorer(metaestimator):
    """Test that any requested metadata is correctly routed to the underlying
    scorers in CV estimators.
    """
    if "scorer_name" not in metaestimator:
        # This test only makes sense for CV estimators
        return

    metaestimator_class = metaestimator["metaestimator"]
    routing_methods = metaestimator["scorer_routing_methods"]
    method_mapping = metaestimator.get("method_mapping", {})

    for method_name in routing_methods:
        kwargs, (estimator, _), (scorer, registry), (cv, _) = get_init_args(
            metaestimator, sub_estimator_consumes=True
        )
        scorer.set_score_request(sample_weight=True)
        if cv:
            cv.set_split_request(groups=True, metadata=True)
        if estimator is not None:
            set_requests(
                estimator,
                method_mapping=method_mapping,
                methods=[method_name],
                metadata_name="sample_weight",
            )
        instance = metaestimator_class(**kwargs)
        method = getattr(instance, method_name)
        method_kwargs = {"sample_weight": sample_weight}
        if "fit" not in method_name:
            instance.fit(X, y)
        method(X, y, **method_kwargs)

        assert registry
        for _scorer in registry:
            check_recorded_metadata(
                obj=_scorer,
                method="score",
                parent=method_name,
                split_params=("sample_weight",),
                **method_kwargs,
            )


@pytest.mark.parametrize("metaestimator", METAESTIMATORS, ids=METAESTIMATOR_IDS)
@config_context(enable_metadata_routing=True)
def test_metadata_is_routed_correctly_to_splitter(metaestimator):
    """Test that any requested metadata is correctly routed to the underlying
    splitters in CV estimators.
    """
    if "cv_routing_methods" not in metaestimator:
        # This test is only for metaestimators accepting a CV splitter
        return

    metaestimator_class = metaestimator["metaestimator"]
    routing_methods = metaestimator["cv_routing_methods"]
    X_ = metaestimator["X"]
    y_ = metaestimator["y"]

    for method_name in routing_methods:
        kwargs, (estimator, _), (scorer, _), (cv, registry) = get_init_args(
            metaestimator, sub_estimator_consumes=True
        )
        if estimator:
            estimator.set_fit_request(sample_weight=False, metadata=False)
        if scorer:
            scorer.set_score_request(sample_weight=False, metadata=False)
        cv.set_split_request(groups=True, metadata=True)
        instance = metaestimator_class(**kwargs)
        method_kwargs = {"groups": groups, "metadata": metadata}
        method = getattr(instance, method_name)
        method(X_, y_, **method_kwargs)
        assert registry
        for _splitter in registry:
            check_recorded_metadata(
                obj=_splitter, method="split", parent=method_name, **method_kwargs
            )


@pytest.mark.parametrize("metaestimator", METAESTIMATORS, ids=METAESTIMATOR_IDS)
@config_context(enable_metadata_routing=True)
def test_metadata_routed_to_group_splitter(metaestimator):
    """Test that groups are routed correctly if group splitter of CV estimator is used
    within cross_validate. Regression test for issue described in PR #29634 to test that
    `ValueError: The 'groups' parameter should not be None.` is not raised."""

    if "cv_routing_methods" not in metaestimator:
        # This test is only for metaestimators accepting a CV splitter
        return

    metaestimator_class = metaestimator["metaestimator"]
    X_ = metaestimator["X"]
    y_ = metaestimator["y"]

    kwargs, *_ = get_init_args(metaestimator, sub_estimator_consumes=True)
    # remove `ConsumingSplitter` from kwargs, so 'cv' param isn't passed twice:
    kwargs.pop("cv", None)
    instance = metaestimator_class(cv=GroupKFold(n_splits=2), **kwargs)
    cross_validate(
        instance,
        X_,
        y_,
        params={"groups": groups},
        cv=GroupKFold(n_splits=2),
        scoring=make_scorer(mean_squared_error, response_method="predict"),
    )
