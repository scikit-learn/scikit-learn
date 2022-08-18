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
from sklearn.utils.metadata_routing import get_routing_for_object
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
    ``exclude`` parameter.
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


def record_metadata(obj, method, record_default=True, **kwargs):
    """Utility function to store passed metadata to a method.

    If record_default is False, kwargs whose values are "default" are skipped.
    This is so that checks on keyword arguments whose default was not changed
    are skipped.

    """
    if not hasattr(obj, "_records"):
        obj._records = {}
    if not record_default:
        kwargs = {
            key: val
            for key, val in kwargs.items()
            if not isinstance(val, str) or (val != "default")
        }
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
        router = MetadataRouter(owner=self.__class__.__name__).add(
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
            MetadataRouter(owner=self.__class__.__name__)
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
            MetadataRouter(owner=self.__class__.__name__)
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
        return MetadataRouter(owner=self.__class__.__name__).add(
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
        router = MetadataRouter(owner=self.__class__.__name__)
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
    requests = MetadataRequest(owner="test")
    assert_request_is_empty(requests)

    requests.fit.add_request(param="foo", alias=RequestType.ERROR_IF_PASSED)
    # this should still work, since ERROR_IF_PASSED is the default value
    assert_request_is_empty(requests)

    requests.fit.add_request(param="bar", alias="value")
    with pytest.raises(AssertionError):
        # now requests is no more empty
        assert_request_is_empty(requests)

    # but one can exclude a method
    assert_request_is_empty(requests, exclude="fit")

    requests.score.add_request(param="carrot", alias=RequestType.REQUESTED)
    with pytest.raises(AssertionError):
        # excluding `fit` is not enough
        assert_request_is_empty(requests, exclude="fit")

    # and excluding both fit and score would avoid an exception
    assert_request_is_empty(requests, exclude=["fit", "score"])

    # test if a router is empty
    assert_request_is_empty(
        MetadataRouter(owner="test")
        .add_self(WeightedMetaRegressor(estimator=None))
        .add(method_mapping="fit", estimator=RegressorMetadata())
    )


def test_default_requests():
    class OddEstimator(BaseEstimator):
        __metadata_request__fit = {
            # set a different default request
            "sample_weight": RequestType.REQUESTED
        }  # type: ignore

    odd_request = get_routing_for_object(OddEstimator())
    assert odd_request.fit.requests == {"sample_weight": RequestType.REQUESTED}

    # check other test estimators
    assert not len(get_routing_for_object(ClassifierNoMetadata()).fit.requests)
    assert_request_is_empty(ClassifierNoMetadata().get_metadata_routing())

    trs_request = get_routing_for_object(TransformerMetadata())
    assert trs_request.fit.requests == {
        "sample_weight": RequestType(None),
        "brand": RequestType(None),
    }
    assert trs_request.transform.requests == {
        "sample_weight": RequestType(None),
    }
    assert_request_is_empty(trs_request)

    est_request = get_routing_for_object(ClassifierFitMetadata())
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
    err_message = (
        "[sample_weight] are passed but are not explicitly set as requested or"
        " not for ClassifierFitMetadata.fit"
    )
    with pytest.raises(ValueError, match=re.escape(err_message)):
        clf.fit(X, y, sample_weight=my_weights)

    # Explicitly saying the estimator doesn't need it, makes the error go away,
    # because in this case `SimpleMetaClassifier` consumes `sample_weight`. If
    # there was no consumer of sample_weight, passing it would result in an
    # error.
    clf = SimpleMetaClassifier(
        estimator=ClassifierFitMetadata().set_fit_request(sample_weight=False)
    )
    # this doesn't raise since SimpleMetaClassifier itself is a consumer,
    # and passing metadata to the consumer directly is fine regardless of its
    # metadata_request values.
    clf.fit(X, y, sample_weight=my_weights)
    check_recorded_metadata(clf.estimator_, "fit", sample_weight=None, brand=None)

    # Requesting a metadata will make the meta-estimator forward it correctly
    clf = SimpleMetaClassifier(
        estimator=ClassifierFitMetadata().set_fit_request(sample_weight=True)
    )
    clf.fit(X, y, sample_weight=my_weights)
    check_recorded_metadata(clf.estimator_, "fit", sample_weight=my_weights, brand=None)

    # And requesting it with an alias
    clf = SimpleMetaClassifier(
        estimator=ClassifierFitMetadata().set_fit_request(
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
                .set_fit_request(brand=True, sample_weight=False)
                .set_transform_request(sample_weight=True)
            ),
            WeightedMetaRegressor(
                estimator=RegressorMetadata().set_fit_request(
                    sample_weight="inner_weights"
                )
            ).set_fit_request(sample_weight="outer_weights"),
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


def test_nested_routing_conflict():
    # check if an error is raised if there's a conflict between keys
    pipeline = SimplePipeline(
        [
            MetaTransformer(
                transformer=TransformerMetadata()
                .set_fit_request(brand=True, sample_weight=False)
                .set_transform_request(sample_weight=True)
            ),
            WeightedMetaRegressor(
                estimator=RegressorMetadata().set_fit_request(sample_weight=True)
            ).set_fit_request(sample_weight="outer_weights"),
        ]
    )
    w1, w2 = [1], [2]
    with pytest.raises(
        ValueError,
        match=(
            re.escape(
                "In WeightedMetaRegressor, there is a conflict on sample_weight between"
                " what is requested for this estimator and what is requested by its"
                " children. You can resolve this conflict by using an alias for the"
                " child estimator(s) requested metadata."
            )
        ),
    ):
        pipeline.fit(X, y, brand=my_groups, sample_weight=w1, outer_weights=w2)


def test_invalid_metadata():
    # check that passing wrong metadata raises an error
    trs = MetaTransformer(
        transformer=TransformerMetadata().set_transform_request(sample_weight=True)
    )
    with pytest.raises(
        TypeError,
        match=(re.escape("transform got unexpected argument(s) {'other_param'}")),
    ):
        trs.fit(X, y).transform(X, other_param=my_weights)

    # passing a metadata which is not requested by any estimator should also raise
    trs = MetaTransformer(
        transformer=TransformerMetadata().set_transform_request(sample_weight=False)
    )
    with pytest.raises(
        TypeError,
        match=(re.escape("transform got unexpected argument(s) {'sample_weight'}")),
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

    est = TestDefaults().set_score_request(my_param="other_param")
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

    est = TestDefaults().set_fit_request(sample_weight=True)
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


def test_setting_default_requests():
    # Test _get_default_requests method
    test_cases = dict()

    class ExplicitRequest(BaseEstimator):
        # `fit` doesn't accept `props` explicitly, but we want to request it
        __metadata_request__fit = {"prop": RequestType.ERROR_IF_PASSED}

        def fit(self, X, y, **kwargs):
            return self

    test_cases[ExplicitRequest] = {"prop": RequestType.ERROR_IF_PASSED}

    class ExplicitRequestOverwrite(BaseEstimator):
        # `fit` explicitly accepts `props`, but we want to change the default
        # request value from ERROR_IF_PASSEd to REQUESTED
        __metadata_request__fit = {"prop": RequestType.REQUESTED}

        def fit(self, X, y, prop=None, **kwargs):
            return self

    test_cases[ExplicitRequestOverwrite] = {"prop": RequestType.REQUESTED}

    class ImplicitRequest(BaseEstimator):
        # `fit` requests `prop` and the default ERROR_IF_PASSED should be used
        def fit(self, X, y, prop=None, **kwargs):
            return self

    test_cases[ImplicitRequest] = {"prop": RequestType.ERROR_IF_PASSED}

    class ImplicitRequestRemoval(BaseEstimator):
        # `fit` (in this class or a parent) requests `prop`, but we don't want
        # it requested at all.
        __metadata_request__fit = {"prop": RequestType.UNUSED}

        def fit(self, X, y, prop=None, **kwargs):
            return self

    test_cases[ImplicitRequestRemoval] = {}

    for Klass, requests in test_cases.items():
        assert get_routing_for_object(Klass()).fit.requests == requests
        assert_request_is_empty(Klass().get_metadata_routing(), exclude="fit")
        Klass().fit(None, None)  # for coverage


def test_method_metadata_request():
    mmr = MethodMetadataRequest(
        router=MetadataRequest(owner="test"), owner="test", method="fit"
    )

    with pytest.raises(
        ValueError, match="alias should be either a valid identifier or"
    ):
        mmr.add_request(param="foo", alias=1.4)

    mmr.add_request(param="foo", alias=None)
    assert mmr.requests == {"foo": RequestType.ERROR_IF_PASSED}
    mmr.add_request(param="foo", alias=False)
    assert mmr.requests == {"foo": RequestType.UNREQUESTED}
    mmr.add_request(param="foo", alias=True)
    assert mmr.requests == {"foo": RequestType.REQUESTED}
    mmr.add_request(param="foo", alias="foo")
    assert mmr.requests == {"foo": RequestType.REQUESTED}
    mmr.add_request(param="foo", alias="bar")
    assert mmr.requests == {"foo": "bar"}
    assert mmr._get_param_names(return_alias=False) == {"foo"}
    assert mmr._get_param_names(return_alias=True) == {"bar"}


def test_get_routing_for_object():
    class Consumer(BaseEstimator):
        __metadata_request__fit = {"prop": RequestType.ERROR_IF_PASSED}

    assert_request_is_empty(get_routing_for_object(None))
    assert_request_is_empty(get_routing_for_object(object()))

    mr = MetadataRequest(owner="test")
    mr.fit.add_request(param="foo", alias="bar")
    mr_factory = get_routing_for_object(mr)
    assert_request_is_empty(mr_factory, exclude="fit")
    assert mr_factory.fit.requests == {"foo": "bar"}

    mr = get_routing_for_object(Consumer())
    assert_request_is_empty(mr, exclude="fit")
    assert mr.fit.requests == {"prop": RequestType.ERROR_IF_PASSED}


def test_metaestimator_warnings():
    class WeightedMetaRegressorWarn(WeightedMetaRegressor):
        __metadata_request__fit = {"sample_weight": RequestType.WARN}

    with pytest.warns(
        UserWarning, match="Support for .* has recently been added to this class"
    ):
        WeightedMetaRegressorWarn(
            estimator=LinearRegression().set_fit_request(sample_weight=False)
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
            MethodMetadataRequest(
                router=MetadataRequest(owner="test"), owner="test", method="fit"
            ).add_request(param="foo", alias="bar"),
            "{'foo': 'bar'}",
        ),
        (
            MetadataRequest(owner="test"),
            "{}",
        ),
        (MethodMapping.from_str("score"), "[{'callee': 'score', 'caller': 'score'}]"),
        (
            MetadataRouter(owner="test").add(
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
            MetadataRouter(owner="test"),
            "add_self",
            {"obj": MetadataRouter(owner="test")},
            ValueError,
            "Given `obj` is neither a `MetadataRequest` nor does it implement",
        ),
        (
            ClassifierFitMetadata(),
            "set_fit_request",
            {"invalid": True},
            TypeError,
            "Unexpected args",
        ),
    ],
)
def test_validations(obj, method, inputs, err_cls, err_msg):
    with pytest.raises(err_cls, match=err_msg):
        getattr(obj, method)(**inputs)


def test_methodmapping():
    mm = (
        MethodMapping()
        .add(caller="fit", callee="transform")
        .add(caller="fit", callee="fit")
    )

    mm_list = list(mm)
    assert mm_list[0] == ("transform", "fit")
    assert mm_list[1] == ("fit", "fit")

    mm = MethodMapping.from_str("one-to-one")
    assert (
        str(mm)
        == "[{'callee': 'fit', 'caller': 'fit'}, {'callee': 'partial_fit', 'caller':"
        " 'partial_fit'}, {'callee': 'predict', 'caller': 'predict'}, {'callee':"
        " 'predict_proba', 'caller': 'predict_proba'}, {'callee':"
        " 'predict_log_proba', 'caller': 'predict_log_proba'}, {'callee':"
        " 'decision_function', 'caller': 'decision_function'}, {'callee': 'score',"
        " 'caller': 'score'}, {'callee': 'split', 'caller': 'split'}, {'callee':"
        " 'transform', 'caller': 'transform'}, {'callee': 'inverse_transform',"
        " 'caller': 'inverse_transform'}]"
    )

    mm = MethodMapping.from_str("score")
    assert repr(mm) == "[{'callee': 'score', 'caller': 'score'}]"


def test_metadatarouter_add_self():
    # adding a MetadataRequest as `self` adds a copy
    request = MetadataRequest(owner="nested")
    request.fit.add_request(param="param", alias=True)
    router = MetadataRouter(owner="test").add_self(request)
    assert str(router._self) == str(request)
    # should be a copy, not the same object
    assert router._self is not request

    # one can add an estimator as self
    est = RegressorMetadata().set_fit_request(sample_weight="my_weights")
    router = MetadataRouter(owner="test").add_self(obj=est)
    assert str(router._self) == str(est.get_metadata_routing())
    assert router._self is not est.get_metadata_routing()

    # adding a consumer+router as self should only add the consumer part
    est = WeightedMetaRegressor(
        estimator=RegressorMetadata().set_fit_request(sample_weight="nested_weights")
    )
    router = MetadataRouter(owner="test").add_self(obj=est)
    # _get_metadata_request() returns the consumer part of the requests
    assert str(router._self) == str(est._get_metadata_request())
    # get_metadata_routing() returns the complete request set, consumer and
    # router included.
    assert str(router._self) != str(est.get_metadata_routing())
    # it should be a copy, not the same object
    assert router._self is not est._get_metadata_request()


def test_metadata_routing_add():
    # adding one with a string `method_mapping`
    router = MetadataRouter(owner="test").add(
        method_mapping="fit",
        est=RegressorMetadata().set_fit_request(sample_weight="weights"),
    )
    assert (
        str(router)
        == "{'est': {'mapping': [{'callee': 'fit', 'caller': 'fit'}], 'router': {'fit':"
        " {'sample_weight': 'weights'}, 'score': {'sample_weight':"
        " <RequestType.ERROR_IF_PASSED: None>}}}}"
    )

    # adding one with an instance of MethodMapping
    router = MetadataRouter(owner="test").add(
        method_mapping=MethodMapping().add(callee="score", caller="fit"),
        est=RegressorMetadata().set_score_request(sample_weight=True),
    )
    assert (
        str(router)
        == "{'est': {'mapping': [{'callee': 'score', 'caller': 'fit'}], 'router':"
        " {'fit': {'sample_weight': <RequestType.ERROR_IF_PASSED: None>}, 'score':"
        " {'sample_weight': <RequestType.REQUESTED: True>}}}}"
    )


def test_metadata_routing_get_param_names():
    router = (
        MetadataRouter(owner="test")
        .add_self(
            WeightedMetaRegressor(estimator=RegressorMetadata()).set_fit_request(
                sample_weight="self_weights"
            )
        )
        .add(
            method_mapping="fit",
            trs=TransformerMetadata().set_fit_request(
                sample_weight="transform_weights"
            ),
        )
    )

    assert (
        str(router)
        == "{'$self': {'fit': {'sample_weight': 'self_weights'}, 'score':"
        " {'sample_weight': <RequestType.ERROR_IF_PASSED: None>}}, 'trs':"
        " {'mapping': [{'callee': 'fit', 'caller': 'fit'}], 'router': {'fit':"
        " {'brand': <RequestType.ERROR_IF_PASSED: None>, 'sample_weight':"
        " 'transform_weights'}, 'transform': {'sample_weight':"
        " <RequestType.ERROR_IF_PASSED: None>}}}}"
    )

    assert router._get_param_names(
        method="fit", return_alias=True, ignore_self=False
    ) == {"transform_weights", "brand", "self_weights"}
    # return_alias=False will return original names for "self"
    assert router._get_param_names(
        method="fit", return_alias=False, ignore_self=False
    ) == {"sample_weight", "brand", "transform_weights"}
    # ignoring self would remove "sample_weight"
    assert router._get_param_names(
        method="fit", return_alias=False, ignore_self=True
    ) == {"brand", "transform_weights"}
    # return_alias is ignored when ignore_self=True
    assert router._get_param_names(
        method="fit", return_alias=True, ignore_self=True
    ) == router._get_param_names(method="fit", return_alias=False, ignore_self=True)


def test_warn_on_invalid_child():
    """Test that we error if the child is not known."""
    with pytest.raises(ValueError, match="Unknown child"):
        MetadataRouter(owner="test").add(
            estimator=LinearRegression(), method_mapping="one-to-one"
        ).warn_on(
            child="invalid",
            method="fit",
            params=None,
            raise_on="1.4",
        )


def test_router_deprecation_warning():
    """This test checks the warning mechanism related to `warn_on`.

    `warn_on` is there to handle backward compatibility in cases where the
    meta-estimator is already doing some routing, and SLEP006 would break
    existing user code. `warn_on` helps converting some of those errors to
    warnings.

    In different scenarios with a meta-estimator and a child estimator we test
    if the warning is raised when it should, an error raised when it should,
    and the combinations of the above cases.
    """

    class MetaEstimator(BaseEstimator, MetaEstimatorMixin):
        def __init__(self, estimator):
            self.estimator = estimator

        def fit(self, X, y, **fit_params):
            routed_params = process_routing(self, "fit", fit_params)
            self.estimator_ = clone(self.estimator).fit(
                X, y, **routed_params.estimator.fit
            )

        def predict(self, X, **predict_params):
            routed_params = process_routing(self, "predict", predict_params)
            return self.estimator_.predict(X, **routed_params.estimator.predict)

        def get_metadata_routing(self):
            return (
                MetadataRouter(owner=self.__class__.__name__)
                .add(estimator=self.estimator, method_mapping="one-to-one")
                .warn_on(
                    child="estimator",
                    method="fit",
                    params=None,
                    raise_on="1.4",
                )
            )

    class Estimator(BaseEstimator):
        def fit(self, X, y, sample_weight=None, groups=None):
            return self

        def predict(self, X, sample_weight=None):
            return np.ones(shape=len(X))

    est = MetaEstimator(estimator=Estimator())
    # the meta-estimator has set (using `warn_on`) to have a warning on `fit`.
    with pytest.warns(
        FutureWarning, match="From version 1.4 this results in the following error"
    ):
        est.fit(X, y, sample_weight=my_weights)

    err_msg = (
        "{params} are passed but are not explicitly set as requested or not for {owner}"
    )
    warn_msg = "From version 1.4 this results in the following error"
    # but predict should raise since there is no warn_on set for it.
    with pytest.raises(
        ValueError,
        match=re.escape(
            err_msg.format(params="[sample_weight]", owner="Estimator.predict")
        ),
    ):
        est.predict(X, sample_weight=my_weights)

    # In this case both a warning and an error are raised. The warning comes
    # from the MetaEstimator, and the error from WeightedMetaRegressor since it
    # doesn't have any warn_on set but sample_weight is passed.
    est = MetaEstimator(estimator=WeightedMetaRegressor(estimator=RegressorMetadata()))
    with pytest.raises(
        ValueError,
        match=re.escape(
            err_msg.format(params="[sample_weight]", owner="RegressorMetadata.fit")
        ),
    ):
        with pytest.warns(FutureWarning, match=warn_msg):
            est.fit(X, y, sample_weight=my_weights)

    class WarningWeightedMetaRegressor(WeightedMetaRegressor):
        """A WeightedMetaRegressor which warns instead."""

        def get_metadata_routing(self):
            router = (
                MetadataRouter(owner=self.__class__.__name__)
                .add_self(self)
                .add(estimator=self.estimator, method_mapping="one-to-one")
                .warn_on(
                    child="estimator",
                    method="fit",
                    params=["sample_weight"],
                    raise_on="1.4",
                )
                .warn_on(
                    child="estimator",
                    method="score",
                    params=["sample_weight"],
                    raise_on="1.4",
                )
            )
            return router

    # Now there's only a warning since both meta-estimators warn.
    est = MetaEstimator(
        estimator=WarningWeightedMetaRegressor(estimator=RegressorMetadata())
    )
    with pytest.warns(FutureWarning, match=warn_msg):
        est.fit(X, y, sample_weight=my_weights)

    # here we should raise because there is no warn_on for groups
    with pytest.raises(
        ValueError,
        match=re.escape(
            err_msg.format(params="[sample_weight, groups]", owner="Estimator.fit")
        ),
    ):
        # the sample_weight should still warn
        with pytest.warns(FutureWarning, match=warn_msg):
            WarningWeightedMetaRegressor(estimator=Estimator()).fit(
                X, y, sample_weight=my_weights, groups=1
            )

    # but if the inner estimator has a non-default request, we fall back to
    # raising an error
    est = MetaEstimator(
        estimator=WarningWeightedMetaRegressor(
            estimator=RegressorMetadata().set_fit_request(sample_weight=True)
        )
    )
    with pytest.raises(
        ValueError,
        match=re.escape(
            err_msg.format(
                params="[sample_weight]", owner="WarningWeightedMetaRegressor.fit"
            )
        ),
    ):
        est.fit(X, y, sample_weight=my_weights)


@pytest.mark.parametrize(
    "estimator, is_default_request",
    [
        (LinearRegression(), True),
        (LinearRegression().set_fit_request(sample_weight=True), False),  # type: ignore
        (WeightedMetaRegressor(estimator=LinearRegression()), True),
        (
            WeightedMetaRegressor(
                estimator=LinearRegression().set_fit_request(  # type: ignore
                    sample_weight=True
                )
            ),
            False,
        ),
        (
            WeightedMetaRegressor(
                estimator=LinearRegression()
            ).set_fit_request(  # type: ignore
                sample_weight=True
            ),
            False,
        ),
    ],
)
def test_is_default_request(estimator, is_default_request):
    """Test the `_is_default_request` machinery.

    It should be `True` only if the user hasn't changed any default values.

    Applies to both `MetadataRouter` and `MetadataRequest`.
    """
    assert estimator.get_metadata_routing()._is_default_request == is_default_request


def test_method_generation():
    # Test if all required request methods are generated.

    # TODO: these test classes can be moved to sklearn.utils._testing once we
    # have a better idea of what the commonly used classes are.
    class SimpleEstimator(BaseEstimator):
        # This class should have no set_{method}_request
        def fit(self, X, y):
            pass

        def partial_fit(self, X, y):
            pass

        def predict(self, X):
            pass

        def predict_proba(self, X):
            pass

        def predict_log_proba(self, X):
            pass

        def decision_function(self, X):
            pass

        def score(self, X, y):
            pass

        def split(self, X, y=None):
            pass

        def transform(self, X):
            pass

        def inverse_transform(self, X):
            pass

    for method in METHODS:
        assert not hasattr(SimpleEstimator(), f"set_{method}_request")

    class SimpleEstimator(BaseEstimator):
        # This class should have every set_{method}_request
        def fit(self, X, y, sample_weight=None):
            pass

        def partial_fit(self, X, y, sample_weight=None):
            pass

        def predict(self, X, sample_weight=None):
            pass

        def predict_proba(self, X, sample_weight=None):
            pass

        def predict_log_proba(self, X, sample_weight=None):
            pass

        def decision_function(self, X, sample_weight=None):
            pass

        def score(self, X, y, sample_weight=None):
            pass

        def split(self, X, y=None, sample_weight=None):
            pass

        def transform(self, X, sample_weight=None):
            pass

        def inverse_transform(self, X, sample_weight=None):
            pass

    for method in METHODS:
        assert hasattr(SimpleEstimator(), f"set_{method}_request")
