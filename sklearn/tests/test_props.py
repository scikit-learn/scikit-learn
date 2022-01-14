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
from sklearn.utils.metadata_routing import get_router_for_object
from sklearn.utils.metadata_routing import MetadataRouter
from sklearn.utils.metadata_routing import MethodMapping
from sklearn.utils.metadata_routing import process_routing
from sklearn.utils._metadata_requests import MethodMetadataRequest
from sklearn.utils._metadata_requests import _MetadataRequester
from sklearn.utils._metadata_requests import METHODS

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
            assert_request_is_empty(route_mapping.router)
        return

    exclude = [] if exclude is None else exclude
    for method in METHODS:
        if method in exclude:
            continue
        mmr = getattr(metadata_request, method)
        props = [
            prop
            for prop, alias in mmr.requests.items()
            if isinstance(alias, str)
            or RequestType(alias) != RequestType.ERROR_IF_PASSED
        ]
        assert not len(props)


def assert_request_equal(request, dictionary):
    for method, requests in dictionary.items():
        mmr = getattr(request, method)
        assert mmr.requests == requests

    empty_methods = [method for method in METHODS if method not in dictionary]
    for method in empty_methods:
        assert not len(getattr(request, method).requests)


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

    def fit(self, X, y, **fit_params):
        params = process_routing(self, "fit", fit_params)
        self.estimator_ = clone(self.estimator).fit(X, y, **params.estimator.fit)

    def get_metadata_routing(self):
        router = MetadataRouter().add(
            estimator=self.estimator, method_mapping="one-to-one"
        )
        return router


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

    def fit(self, X, y, sample_weight=None, **fit_params):
        record_metadata(self, "fit", sample_weight=sample_weight)
        params = process_routing(self, "fit", fit_params, sample_weight=sample_weight)
        self.estimator_ = clone(self.estimator).fit(X, y, **params.estimator.fit)
        return self

    def predict(self, X, **predict_params):
        params = process_routing(self, "predict", predict_params)
        return self.estimator_.predict(X, **params.estimator.predict)

    def get_metadata_routing(self):
        router = (
            MetadataRouter()
            .add_self(self)
            .add(estimator=self.estimator, method_mapping="one-to-one")
        )
        return router


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

    def fit(self, X, y, sample_weight=None, **kwargs):
        record_metadata(self, "fit", sample_weight=sample_weight)
        params = process_routing(self, "fit", kwargs, sample_weight=sample_weight)
        self.estimator_ = clone(self.estimator).fit(X, y, **params.estimator.fit)
        return self

    def get_metadata_routing(self):
        router = (
            MetadataRouter()
            .add_self(self)
            .add(estimator=self.estimator, method_mapping="fit")
        )
        return router


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

    def fit(self, X, y=None, **fit_params):
        params = process_routing(self, "fit", fit_params)
        self.transformer_ = clone(self.transformer).fit(X, y, **params.transformer.fit)
        return self

    def transform(self, X, y=None, **transform_params):
        params = process_routing(self, "transform", transform_params)
        return self.transformer_.transform(X, **params.transformer.transform)

    def get_metadata_routing(self):
        return MetadataRouter().add(
            transformer=self.transformer, method_mapping="one-to-one"
        )


class SimplePipeline(BaseEstimator):
    """A very simple pipeline, assuming the last step is always a predictor."""

    def __init__(self, steps):
        self.steps = steps

    def fit(self, X, y, **fit_params):
        self.steps_ = []
        params = process_routing(self, "fit", fit_params)
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

    def predict(self, X, **predict_params):
        check_is_fitted(self)
        X_transformed = X
        params = process_routing(self, "predict", predict_params)
        for i, step in enumerate(self.steps_[:-1]):
            X_transformed = step.transform(X, **params.get(f"step_{i}").transform)

        return self.steps_[-1].predict(X_transformed, **params.predictor.predict)

    def get_metadata_routing(self):
        router = MetadataRouter()
        for i, step in enumerate(self.steps[:-1]):
            router.add(
                **{f"step_{i}": step},
                method_mapping=MethodMapping()
                .add(callee="fit", caller="fit")
                .add(callee="transform", caller="fit")
                .add(callee="transform", caller="predict"),
            )
        router.add(predictor=self.steps[-1], method_mapping="one-to-one")
        return router


def test_assert_request_is_empty():
    requests = MetadataRequest()
    assert_request_is_empty(requests)

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

    odd_request = get_router_for_object(OddEstimator())
    assert odd_request.fit.requests == {"sample_weight": RequestType.REQUESTED}

    # check other test estimators
    assert not len(get_router_for_object(ClassifierNoMetadata()).fit.requests)
    assert_request_is_empty(ClassifierNoMetadata().get_metadata_routing())

    trs_request = get_router_for_object(TransformerMetadata())
    assert trs_request.fit.requests == {
        "sample_weight": RequestType(None),
        "brand": RequestType(None),
    }
    assert trs_request.transform.requests == {
        "sample_weight": RequestType(None),
    }
    assert_request_is_empty(trs_request)

    est_request = get_router_for_object(ClassifierFitMetadata())
    assert est_request.fit.requests == {
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
        "score": {
            "my_param": RequestType.REQUESTED,
            "my_other_param": RequestType.ERROR_IF_PASSED,
            "sample_weight": RequestType.ERROR_IF_PASSED,
        },
        "fit": {
            "my_other_param": RequestType.ERROR_IF_PASSED,
            "sample_weight": RequestType.ERROR_IF_PASSED,
        },
        "predict": {"my_param": RequestType.REQUESTED},
    }
    assert_request_equal(TestDefaults().get_metadata_routing(), expected)

    est = TestDefaults().score_requests(my_param="other_param")
    expected = {
        "score": {
            "my_param": "other_param",
            "my_other_param": RequestType.ERROR_IF_PASSED,
            "sample_weight": RequestType.ERROR_IF_PASSED,
        },
        "fit": {
            "my_other_param": RequestType.ERROR_IF_PASSED,
            "sample_weight": RequestType.ERROR_IF_PASSED,
        },
        "predict": {"my_param": RequestType.REQUESTED},
    }
    assert_request_equal(est.get_metadata_routing(), expected)

    est = TestDefaults().fit_requests(sample_weight=True)
    expected = {
        "score": {
            "my_param": RequestType.REQUESTED,
            "my_other_param": RequestType.ERROR_IF_PASSED,
            "sample_weight": RequestType.ERROR_IF_PASSED,
        },
        "fit": {
            "my_other_param": RequestType.ERROR_IF_PASSED,
            "sample_weight": RequestType.REQUESTED,
        },
        "predict": {"my_param": RequestType.REQUESTED},
    }
    assert_request_equal(est.get_metadata_routing(), expected)


def test__get_default_requests():
    # Test _get_default_requests method
    class ExplicitRequest(BaseEstimator):
        __metadata_request__fit = {"prop": RequestType.ERROR_IF_PASSED}

        def fit(self, X, y):
            return self

    assert get_router_for_object(ExplicitRequest()).fit.requests == {
        "prop": RequestType.ERROR_IF_PASSED
    }
    assert_request_is_empty(ExplicitRequest().get_metadata_routing(), exclude="fit")

    class ExplicitRequestOverwrite(BaseEstimator):
        __metadata_request__fit = {"prop": RequestType.REQUESTED}

        def fit(self, X, y, prop=None, **kwargs):
            return self

    assert get_router_for_object(ExplicitRequestOverwrite()).fit.requests == {
        "prop": RequestType.REQUESTED
    }
    assert_request_is_empty(
        ExplicitRequestOverwrite().get_metadata_routing(), exclude="fit"
    )

    class ImplicitRequest(BaseEstimator):
        def fit(self, X, y, prop=None, **kwargs):
            return self

    assert get_router_for_object(ImplicitRequest()).fit.requests == {
        "prop": RequestType.ERROR_IF_PASSED
    }
    assert_request_is_empty(ImplicitRequest().get_metadata_routing(), exclude="fit")

    class ImplicitRequestRemoval(BaseEstimator):
        __metadata_request__fit = {"prop": RequestType.UNUSED}

        def fit(self, X, y, prop=None, **kwargs):
            return self

    assert get_router_for_object(ImplicitRequestRemoval()).fit.requests == {}
    assert_request_is_empty(ImplicitRequestRemoval().get_metadata_routing())


def test_method_metadata_request():
    mmr = MethodMetadataRequest(name="fit")

    with pytest.raises(
        ValueError, match="alias should be either a valid identifier or"
    ):
        mmr.add_request(prop="foo", alias=1.4)

    mmr.add_request(prop="foo", alias=None)
    assert mmr.requests == {"foo": RequestType.ERROR_IF_PASSED}
    mmr.add_request(prop="foo", alias=False)
    assert mmr.requests == {"foo": RequestType.UNREQUESTED}
    mmr.add_request(prop="foo", alias=True)
    assert mmr.requests == {"foo": RequestType.REQUESTED}
    mmr.add_request(prop="foo", alias="foo")
    assert mmr.requests == {"foo": RequestType.REQUESTED}
    mmr.add_request(prop="foo", alias="bar")
    assert mmr.requests == {"foo": "bar"}
    assert mmr._get_param_names(original_names=True) == {"foo"}
    assert mmr._get_param_names(original_names=False) == {"bar"}


def test_get_router_for_object():
    class Consumer(BaseEstimator):
        __metadata_request__fit = {"prop": RequestType.ERROR_IF_PASSED}

    assert_request_is_empty(get_router_for_object(None))
    assert_request_is_empty(get_router_for_object(object()))

    mr = MetadataRequest()
    mr.fit.add_request(prop="foo", alias="bar")
    mr_factory = get_router_for_object(mr)
    assert_request_is_empty(mr_factory, exclude="fit")
    assert mr_factory.fit.requests == {"foo": "bar"}

    mr = get_router_for_object(Consumer())
    assert_request_is_empty(mr, exclude="fit")
    assert mr.fit.requests == {"prop": RequestType.ERROR_IF_PASSED}


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
            "{}",
        ),
        (MethodMapping.from_str("score"), "[{'callee': 'score', 'caller': 'score'}]"),
        (
            MetadataRouter().add(
                method_mapping="predict", estimator=RegressorMetadata()
            ),
            "{'estimator': {'mapping': [{'callee': 'predict', 'caller': 'predict'}],"
            " 'router': {'fit': {'sample_weight': <RequestType.ERROR_IF_PASSED: None>},"
            " 'score': {'sample_weight': <RequestType.ERROR_IF_PASSED: None>}}}}",
        ),
    ],
)
def test_string_representations(obj, string):
    assert str(obj) == string


@pytest.mark.parametrize(
    "obj, method, inputs, err_cls, err_msg",
    [
        (
            MethodMapping(),
            "add",
            {"callee": "invalid", "caller": "fit"},
            ValueError,
            "Given callee",
        ),
        (
            MethodMapping(),
            "add",
            {"callee": "fit", "caller": "invalid"},
            ValueError,
            "Given caller",
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
def test_validations(obj, method, inputs, err_cls, err_msg):
    with pytest.raises(err_cls, match=err_msg):
        getattr(obj, method)(**inputs)


def test_requestmethod():
    with pytest.raises(TypeError, match="Unexpected args"):
        ClassifierFitMetadata().fit_requests(invalid=True)
