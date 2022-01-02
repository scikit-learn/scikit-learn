"""
Metadata Routing Utility Tests
"""

# Author: Adrin Jalali <adrin.jalali@gmail.com>
# License: BSD 3 clause

import re
import numpy as np
import pytest

from sklearn.base import BaseEstimator
from sklearn.base import ClassifierMixin
from sklearn.base import RegressorMixin
from sklearn.base import TransformerMixin
from sklearn.base import MetaEstimatorMixin
from sklearn.base import clone
from sklearn.linear_model import LinearRegression
from sklearn.utils.validation import check_is_fitted
from sklearn.utils.metadata_routing import RequestType
from sklearn.utils.metadata_routing import MetadataRequest
from sklearn.utils.metadata_routing import metadata_router_factory
from sklearn.utils.metadata_routing import MetadataRouter
from sklearn.utils.metadata_routing import MethodMapping
from sklearn.utils.metadata_routing import process_routing
from sklearn.utils._metadata_requests import MethodMetadataRequest
from sklearn.utils._metadata_requests import _MetadataRequester

N, M = 100, 4
X = np.random.rand(N, M)
y = np.random.randint(0, 2, size=N)
my_groups = np.random.randint(0, 10, size=N)
my_weights = np.random.rand(N)
my_other_weights = np.random.rand(N)


def assert_request_is_empty(metadata_request, exclude=None):
    """Check if a metadata request dict is empty.

    One can exclude a method or a list of methods from the check using the
    ``exclude`` perameter.
    """
    if isinstance(metadata_request, MetadataRouter):
        for _, route_mapping in metadata_request:
            assert_request_is_empty(route_mapping.routing)
        return

    if isinstance(metadata_request, MetadataRequest):
        metadata_request = metadata_request.serialize()
    if exclude is None:
        exclude = []
    for method, request in metadata_request.items():
        if method in exclude or method == "^type":
            continue
        props = [
            prop
            for prop, alias in request.items()
            if isinstance(alias, str)
            or RequestType(alias) != RequestType.ERROR_IF_PASSED
        ]
        assert not len(props)


def record_metadata(obj, method, **kwargs):
    """Utility function to store passed metadata to a method."""
    if not hasattr(obj, "_records"):
        setattr(obj, "_records", dict())
    obj._records[method] = kwargs


def check_recorded_metadata(obj, method, **kwargs):
    """Check whether the expected metadata is passed to the object's method."""
    records = getattr(obj, "_records", dict()).get(method, dict())
    assert set(kwargs.keys()) == set(records.keys())
    for key, value in kwargs.items():
        assert records[key] is value


class MetaRegressor(MetaEstimatorMixin, RegressorMixin, BaseEstimator):
    """A meta-regressor which is only a router."""

    def __init__(self, estimator):
        self.estimator = estimator

    @process_routing
    def fit(self, X, y, **fit_params):
        params = fit_params["params"]
        self.estimator_ = clone(self.estimator).fit(X, y, **params.estimator.fit)

    def get_metadata_routing(self):
        router = MetadataRouter().add(
            estimator=self.estimator, method_mapping="one-to-one"
        )
        return router.serialize()


class RegressorMetadata(RegressorMixin, BaseEstimator):
    """A regressor consuming a metadata."""

    def fit(self, X, y, sample_weight=None):
        record_metadata(self, "fit", sample_weight=sample_weight)
        return self

    def predict(self, X):
        return np.zeros(shape=(len(X)))


class WeightedMetaRegressor(MetaEstimatorMixin, RegressorMixin, BaseEstimator):
    """A meta-regressor which is also a consumer."""

    def __init__(self, estimator):
        self.estimator = estimator

    @process_routing
    def fit(self, X, y, sample_weight=None, **fit_params):
        record_metadata(self, "fit", sample_weight=sample_weight)
        params = fit_params["params"]
        self.estimator_ = clone(self.estimator).fit(X, y, **params.estimator.fit)
        return self

    @process_routing
    def predict(self, X, **predict_params):
        params = predict_params["params"]
        return self.estimator_.predict(X, **params.estimator.predict)

    def get_metadata_routing(self):
        router = (
            MetadataRouter()
            .add_self(self)
            .add(estimator=self.estimator, method_mapping="one-to-one")
        )
        return router.serialize()


class ClassifierNoMetadata(ClassifierMixin, BaseEstimator):
    """An estimator which accepts no metadata on any method."""

    def fit(self, X, y):
        return self

    def predict(self, X):
        return np.ones(len(X))


class ClassifierFitMetadata(ClassifierMixin, BaseEstimator):
    """An estimator accepting two metadata in its ``fit`` method."""

    def fit(self, X, y, sample_weight=None, brand=None):
        record_metadata(self, "fit", sample_weight=sample_weight, brand=brand)
        return self

    def predict(self, X):
        return np.ones(len(X))


class SimpleMetaClassifier(MetaEstimatorMixin, ClassifierMixin, BaseEstimator):
    """A meta-estimator which also consumes sample_weight itself in ``fit``."""

    def __init__(self, estimator):
        self.estimator = estimator

    @process_routing
    def fit(self, X, y, sample_weight=None, **kwargs):
        record_metadata(self, "fit", sample_weight=sample_weight)
        params = kwargs["params"]
        self.estimator_ = clone(self.estimator).fit(X, y, **params.estimator.fit)
        return self

    def get_metadata_routing(self):
        router = (
            MetadataRouter()
            .add_self(self)
            .add(estimator=self.estimator, method_mapping="fit")
        )
        return router.serialize()


class TransformerMetadata(TransformerMixin, BaseEstimator):
    """A transformer which accepts metadata on fit and transform."""

    def fit(self, X, y=None, brand=None, sample_weight=None):
        record_metadata(self, "fit", brand=brand, sample_weight=sample_weight)
        return self

    def transform(self, X, sample_weight=None):
        record_metadata(self, "transform", sample_weight=sample_weight)
        return X


class MetaTransformer(MetaEstimatorMixin, TransformerMixin, BaseEstimator):
    """A simple meta-transformer."""

    def __init__(self, transformer):
        self.transformer = transformer

    @process_routing
    def fit(self, X, y=None, **fit_params):
        params = fit_params["params"]
        self.transformer_ = clone(self.transformer).fit(X, y, **params.transformer.fit)
        return self

    @process_routing
    def transform(self, X, y=None, **transform_params):
        params = transform_params["params"]
        return self.transformer_.transform(X, **params.transformer.transform)

    def get_metadata_routing(self):
        return (
            MetadataRouter()
            .add(transformer=self.transformer, method_mapping="one-to-one")
            .serialize()
        )


class SimplePipeline(BaseEstimator):
    """A very simple pipeline, assuming the last step is always a predictor."""

    def __init__(self, steps):
        self.steps = steps

    @process_routing
    def fit(self, X, y, **fit_params):
        self.steps_ = []
        params = fit_params["params"]
        X_transformed = X
        for i, step in enumerate(self.steps[:-1]):
            transformer = clone(step).fit(
                X_transformed, y, **params.get(f"step_{i}").fit
            )
            self.steps_.append(transformer)
            X_transformed = transformer.transform(
                X_transformed, **params.get(f"step_{i}").transform
            )

        self.steps_.append(
            clone(self.steps[-1]).fit(X_transformed, y, **params.predictor.fit)
        )
        return self

    @process_routing
    def predict(self, X, **predict_params):
        check_is_fitted(self)
        X_transformed = X
        params = predict_params["params"]
        for i, step in enumerate(self.steps_[:-1]):
            X_transformed = step.transform(X, **params.get(f"step_{i}").transform)

        return self.steps_[-1].predict(X_transformed, **params.predictor.predict)

    def get_metadata_routing(self):
        router = MetadataRouter()
        for i, step in enumerate(self.steps[:-1]):
            router.add(
                **{f"step_{i}": step},
                method_mapping=MethodMapping()
                .add(method="fit", used_in="fit")
                .add(method="transform", used_in="fit")
                .add(method="transform", used_in="predict"),
            )
        router.add(predictor=self.steps[-1], method_mapping="one-to-one")
        return router.serialize()


def test_assert_request_is_empty():
    requests = MetadataRequest()
    assert_request_is_empty(requests)
    assert_request_is_empty(requests.serialize())

    requests.fit.add_request(prop="foo", alias=RequestType.ERROR_IF_PASSED)
    # this should still work, since ERROR_IF_PASSED is the default value
    assert_request_is_empty(requests)

    requests.fit.add_request(prop="bar", alias="value")
    with pytest.raises(AssertionError):
        # now requests is no more empty
        assert_request_is_empty(requests)

    # but one can exclude a method
    assert_request_is_empty(requests, exclude="fit")

    requests.score.add_request(prop="carrot", alias=RequestType.REQUESTED)
    with pytest.raises(AssertionError):
        # excluding `fit` is not enough
        assert_request_is_empty(requests, exclude="fit")

    # and excluding both fit and score would avoid an exception
    assert_request_is_empty(requests, exclude=["fit", "score"])

    # test if a router is empty
    assert_request_is_empty(
        MetadataRouter()
        .add_self(WeightedMetaRegressor(estimator=None))
        .add(method_mapping="fit", estimator=RegressorMetadata())
    )


def test_default_requests():
    class OddEstimator(BaseEstimator):
        __metadata_request__fit = {
            # set a different default request
            "sample_weight": RequestType.REQUESTED
        }  # type: ignore

    odd_request = metadata_router_factory(OddEstimator())
    assert odd_request.fit._requests == {"sample_weight": RequestType.REQUESTED}

    # check other test estimators
    assert not len(metadata_router_factory(ClassifierNoMetadata()).fit._requests)
    assert_request_is_empty(ClassifierNoMetadata().get_metadata_routing())

    trs_request = metadata_router_factory(TransformerMetadata())
    assert trs_request.fit._requests == {
        "sample_weight": RequestType(None),
        "brand": RequestType(None),
    }
    assert trs_request.transform._requests == {
        "sample_weight": RequestType(None),
    }
    assert_request_is_empty(trs_request)

    est_request = metadata_router_factory(ClassifierFitMetadata())
    assert est_request.fit._requests == {
        "sample_weight": RequestType(None),
        "brand": RequestType(None),
    }
    assert_request_is_empty(est_request)


def test_simple_metadata_routing():
    # Tests that metadata is properly routed

    # The underlying estimator doesn't accept or request metadata
    clf = SimpleMetaClassifier(estimator=ClassifierNoMetadata())
    clf.fit(X, y)

    # Meta-estimator consumes sample_weight, but doesn't forward it to the underlying
    # estimator
    clf = SimpleMetaClassifier(estimator=ClassifierNoMetadata())
    clf.fit(X, y, sample_weight=my_weights)

    # If the estimator accepts the metadata but doesn't explicitly say it doesn't
    # need it, there's an error
    clf = SimpleMetaClassifier(estimator=ClassifierFitMetadata())
    with pytest.raises(
        ValueError,
        match=(
            "sample_weight is passed but is not explicitly set as requested or not. In"
            " method: fit"
        ),
    ):
        clf.fit(X, y, sample_weight=my_weights)

    # Explicitly saying the estimator doesn't need it, makes the error go away,
    # but if a metadata is passed which is not requested by any object/estimator,
    # there will be still an error
    clf = SimpleMetaClassifier(
        estimator=ClassifierFitMetadata().fit_requests(sample_weight=False)
    )
    # this doesn't raise since SimpleMetaClassifier itself is a consumer,
    # and passing metadata to the consumer directly is fine regardless of its
    # metadata_request values.
    clf.fit(X, y, sample_weight=my_weights)
    check_recorded_metadata(clf.estimator_, "fit", sample_weight=None, brand=None)

    # Requesting a metadata will make the meta-estimator forward it correctly
    clf = SimpleMetaClassifier(
        estimator=ClassifierFitMetadata().fit_requests(sample_weight=True)
    )
    clf.fit(X, y, sample_weight=my_weights)
    check_recorded_metadata(clf.estimator_, "fit", sample_weight=my_weights, brand=None)

    # And requesting it with an alias
    clf = SimpleMetaClassifier(
        estimator=ClassifierFitMetadata().fit_requests(
            sample_weight="alternative_weight"
        )
    )
    clf.fit(X, y, alternative_weight=my_weights)
    check_recorded_metadata(clf.estimator_, "fit", sample_weight=my_weights, brand=None)


def test_nested_routing():
    # check if metadata is routed in a nested routing situation.
    pipeline = SimplePipeline(
        [
            MetaTransformer(
                transformer=TransformerMetadata()
                .fit_requests(brand=True, sample_weight=False)
                .transform_requests(sample_weight=True)
            ),
            WeightedMetaRegressor(
                estimator=RegressorMetadata().fit_requests(
                    sample_weight="inner_weights"
                )
            ).fit_requests(sample_weight="outer_weights"),
        ]
    )
    w1, w2, w3 = [1], [2], [3]
    pipeline.fit(
        X, y, brand=my_groups, sample_weight=w1, outer_weights=w2, inner_weights=w3
    )
    check_recorded_metadata(
        pipeline.steps_[0].transformer_, "fit", brand=my_groups, sample_weight=None
    )
    check_recorded_metadata(
        pipeline.steps_[0].transformer_, "transform", sample_weight=w1
    )
    check_recorded_metadata(pipeline.steps_[1], "fit", sample_weight=w2)
    check_recorded_metadata(pipeline.steps_[1].estimator_, "fit", sample_weight=w3)

    pipeline.predict(X, sample_weight=w3)
    check_recorded_metadata(
        pipeline.steps_[0].transformer_, "transform", sample_weight=w3
    )


def test_invalid_metadata():
    # check that passing wrong metadata raises an error
    trs = MetaTransformer(
        transformer=TransformerMetadata().transform_requests(sample_weight=True)
    )
    with pytest.raises(
        ValueError,
        match=(
            re.escape(
                "These passed parameters are not understood or requested by any object:"
                " {'other_param'}"
            )
        ),
    ):
        trs.fit(X, y).transform(X, other_param=my_weights)

    # passing a metadata which is not requested by any estimator should also raise
    trs = MetaTransformer(
        transformer=TransformerMetadata().transform_requests(sample_weight=False)
    )
    with pytest.raises(
        ValueError,
        match=(
            re.escape(
                "These passed parameters are not understood or requested by any object:"
                " {'sample_weight'}"
            )
        ),
    ):
        trs.fit(X, y).transform(X, sample_weight=my_weights)


def test_get_metadata_routing():
    class TestDefaultsBadMethodName(_MetadataRequester):
        __metadata_request__fit = {
            "sample_weight": RequestType.ERROR_IF_PASSED,
            "my_param": RequestType.ERROR_IF_PASSED,
        }
        __metadata_request__score = {
            "sample_weight": RequestType.ERROR_IF_PASSED,
            "my_param": True,
            "my_other_param": RequestType.ERROR_IF_PASSED,
        }
        # this will raise an error since we don't understand "other_method" as a method
        __metadata_request__other_method = {"my_param": True}

    class TestDefaults(_MetadataRequester):
        __metadata_request__fit = {
            "sample_weight": RequestType.ERROR_IF_PASSED,
            "my_other_param": RequestType.ERROR_IF_PASSED,
        }
        __metadata_request__score = {
            "sample_weight": RequestType.ERROR_IF_PASSED,
            "my_param": True,
            "my_other_param": RequestType.ERROR_IF_PASSED,
        }
        __metadata_request__predict = {"my_param": True}

    with pytest.raises(
        AttributeError, match="'MetadataRequest' object has no attribute 'other_method'"
    ):
        TestDefaultsBadMethodName().get_metadata_routing()

    expected = {
        "^type": "request",
        "score": {
            "my_param": True,
            "my_other_param": None,
            "sample_weight": None,
        },
        "fit": {
            "my_other_param": None,
            "sample_weight": None,
        },
        "partial_fit": {},
        "predict": {"my_param": True},
        "transform": {},
        "inverse_transform": {},
        "split": {},
    }
    assert TestDefaults().get_metadata_routing() == expected

    est = TestDefaults().score_requests(my_param="other_param")
    expected = {
        "^type": "request",
        "score": {
            "my_param": "other_param",
            "my_other_param": None,
            "sample_weight": None,
        },
        "fit": {
            "my_other_param": None,
            "sample_weight": None,
        },
        "partial_fit": {},
        "predict": {"my_param": True},
        "transform": {},
        "inverse_transform": {},
        "split": {},
    }
    assert est.get_metadata_routing() == expected

    est = TestDefaults().fit_requests(sample_weight=True)
    expected = {
        "^type": "request",
        "score": {
            "my_param": True,
            "my_other_param": None,
            "sample_weight": None,
        },
        "fit": {
            "my_other_param": None,
            "sample_weight": True,
        },
        "partial_fit": {},
        "predict": {"my_param": True},
        "transform": {},
        "inverse_transform": {},
        "split": {},
    }
    assert est.get_metadata_routing() == expected


def test__get_default_requests():
    # Test _get_default_requests method
    class ExplicitRequest(BaseEstimator):
        __metadata_request__fit = {"prop": RequestType.ERROR_IF_PASSED}

        def fit(self, X, y):
            return self

    assert metadata_router_factory(ExplicitRequest()).fit._requests == {
        "prop": RequestType.ERROR_IF_PASSED
    }
    assert_request_is_empty(ExplicitRequest().get_metadata_routing(), exclude="fit")

    class ExplicitRequestOverwrite(BaseEstimator):
        __metadata_request__fit = {"prop": RequestType.REQUESTED}

        def fit(self, X, y, prop=None, **kwargs):
            return self

    assert metadata_router_factory(ExplicitRequestOverwrite()).fit._requests == {
        "prop": RequestType.REQUESTED
    }
    assert_request_is_empty(
        ExplicitRequestOverwrite().get_metadata_routing(), exclude="fit"
    )

    class ImplicitRequest(BaseEstimator):
        def fit(self, X, y, prop=None, **kwargs):
            return self

    assert metadata_router_factory(ImplicitRequest()).fit._requests == {
        "prop": RequestType.ERROR_IF_PASSED
    }
    assert_request_is_empty(ImplicitRequest().get_metadata_routing(), exclude="fit")

    class ImplicitRequestRemoval(BaseEstimator):
        __metadata_request__fit = {"prop": RequestType.UNUSED}

        def fit(self, X, y, prop=None, **kwargs):
            return self

    assert metadata_router_factory(ImplicitRequestRemoval()).fit._requests == {}
    assert_request_is_empty(ImplicitRequestRemoval().get_metadata_routing())


def test_method_metadata_request():
    mmr = MethodMetadataRequest(name="fit")

    with pytest.raises(ValueError, match="alias should be either a string or"):
        mmr.add_request(prop="foo", alias=1.4)

    mmr.add_request(prop="foo", alias=None)
    assert mmr._requests == {"foo": RequestType.ERROR_IF_PASSED}
    mmr.add_request(prop="foo", alias=False)
    assert mmr._requests == {"foo": RequestType.UNREQUESTED}
    mmr.add_request(prop="foo", alias=True)
    assert mmr._requests == {"foo": RequestType.REQUESTED}
    mmr.add_request(prop="foo", alias="foo")
    assert mmr._requests == {"foo": RequestType.REQUESTED}
    mmr.add_request(prop="foo", alias="bar")
    assert mmr._requests == {"foo": "bar"}
    assert mmr._get_param_names(original_names=True) == {"foo"}
    assert mmr._get_param_names(original_names=False) == {"bar"}

    assert MethodMetadataRequest.deserialize(None, name="fit")._requests == {}
    assert MethodMetadataRequest.deserialize({"foo": None}, name="fit")._requests == {
        "foo": RequestType.ERROR_IF_PASSED
    }
    assert MethodMetadataRequest.deserialize(
        {"foo": None, "bar": None}, name="fit"
    )._requests == {
        "foo": RequestType.ERROR_IF_PASSED,
        "bar": RequestType.ERROR_IF_PASSED,
    }


def test_metadata_router_factory():
    class Consumer(BaseEstimator):
        __metadata_request__fit = {"prop": RequestType.ERROR_IF_PASSED}

    assert_request_is_empty(metadata_router_factory(None))
    assert_request_is_empty(metadata_router_factory({"^type": "request"}))
    assert_request_is_empty(metadata_router_factory({"^type": "router"}))
    with pytest.raises(ValueError, match="Cannot understand object type"):
        metadata_router_factory({"^type": "invalid"})
    assert_request_is_empty(metadata_router_factory(object()))

    mr = MetadataRequest.deserialize({"^type": "request", "fit": {"foo": "bar"}})
    mr_factory = metadata_router_factory(mr)
    assert_request_is_empty(mr_factory, exclude="fit")
    assert mr_factory.fit._requests == {"foo": "bar"}

    mr = metadata_router_factory(Consumer())
    assert_request_is_empty(mr, exclude="fit")
    assert mr.fit._requests == {"prop": RequestType.ERROR_IF_PASSED}


def test_metaestimator_warnings():
    class WeightedMetaRegressorWarn(WeightedMetaRegressor):
        __metadata_request__fit = {"sample_weight": RequestType.WARN}

    with pytest.warns(
        UserWarning, match="Support for .* has recently been added to this class"
    ):
        WeightedMetaRegressorWarn(
            estimator=LinearRegression().fit_requests(sample_weight=False)
        ).fit(X, y, sample_weight=my_weights)


def test_estimator_warnings():
    class RegressorMetadataWarn(RegressorMetadata):
        __metadata_request__fit = {"sample_weight": RequestType.WARN}

    with pytest.warns(
        UserWarning, match="Support for .* has recently been added to this class"
    ):
        MetaRegressor(estimator=RegressorMetadataWarn()).fit(
            X, y, sample_weight=my_weights
        )


@pytest.mark.parametrize(
    "obj, string",
    [
        (
            MethodMetadataRequest("fit").add_request(prop="foo", alias="bar"),
            "{'foo': 'bar'}",
        ),
        (
            MetadataRequest(),
            "{'^type': 'request', 'fit': {}, 'partial_fit': {}, 'predict': {}, 'score':"
            " {}, 'split': {}, 'transform': {}, 'inverse_transform': {}}",
        ),
        (MethodMapping.from_str("score"), "[{'method': 'score', 'used_in': 'score'}]"),
        (
            MetadataRouter().add(
                method_mapping="predict", estimator=RegressorMetadata()
            ),
            "{'^type': 'router', 'estimator': {'mapping': [{'method': 'predict',"
            " 'used_in': 'predict'}], 'routing': {'^type': 'request', 'fit':"
            " {'sample_weight': None}, 'partial_fit': {}, 'predict': {}, 'score':"
            " {'sample_weight': None}, 'split': {}, 'transform': {},"
            " 'inverse_transform': {}}}}",
        ),
    ],
)
def test_string_representations(obj, string):
    assert str(obj) == string


@pytest.mark.parametrize(
    "obj, method, inputs, err_cls, err_msg",
    [
        (
            MetadataRequest,
            "deserialize",
            {"obj": {}},
            ValueError,
            "Can only create a metadata request of type",
        ),
        (
            MetadataRequest,
            "deserialize",
            {"obj": {"^type": "request", "invalid_method": {}}},
            ValueError,
            "is not supported as a method",
        ),
        (
            MetadataRouter,
            "deserialize",
            {"obj": {}},
            ValueError,
            "Can only create a router of type",
        ),
        (
            MethodMapping(),
            "add",
            {"method": "invalid", "used_in": "fit"},
            ValueError,
            "Given method",
        ),
        (
            MethodMapping(),
            "add",
            {"method": "fit", "used_in": "invalid"},
            ValueError,
            "Given used_in",
        ),
        (
            MethodMapping,
            "from_str",
            {"route": "invalid"},
            ValueError,
            "route should be 'one-to-one' or a single method!",
        ),
        (
            MetadataRouter(),
            "add_self",
            {"obj": MetadataRouter()},
            ValueError,
            "Given object is neither a MetadataRequest nor does it implement",
        ),
    ],
)
def test_deserialize_invalid_type(obj, method, inputs, err_cls, err_msg):
    with pytest.raises(err_cls, match=err_msg):
        getattr(obj, method)(**inputs)


def test_requestmethod():
    with pytest.raises(TypeError, match="Unexpected args"):
        ClassifierFitMetadata().fit_requests(invalid=True)
