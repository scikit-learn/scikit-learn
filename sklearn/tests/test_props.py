import re
import numpy as np
import pytest

from sklearn.base import BaseEstimator
from sklearn.base import ClassifierMixin
from sklearn.base import TransformerMixin
from sklearn.base import MetaEstimatorMixin
from sklearn.base import clone
from sklearn.utils import MetadataRequest
from sklearn.utils.validation import check_is_fitted
from sklearn.utils.metadata_requests import RequestType
from sklearn.utils.metadata_requests import metadata_request_factory
from sklearn.utils.metadata_requests import MetadataRouter
from sklearn.utils.metadata_requests import MethodMetadataRequest
from sklearn.utils.metadata_requests import MethodMapping

from sklearn.base import _MetadataRequester

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

    if isinstance(metadata_request, MetadataRequest):
        metadata_request = metadata_request.to_dict()
    if exclude is None:
        exclude = []
    for method, request in metadata_request.items():
        if method in exclude or method == "_type":
            continue
        props = [
            prop
            for prop, alias in request.items()
            if isinstance(alias, str)
            or RequestType(alias) != RequestType.ERROR_IF_PASSED
        ]
        assert not len(props)


class TestEstimatorNoMetadata(ClassifierMixin, BaseEstimator):
    """An estimator which accepts no metadata on any method."""

    def fit(self, X, y):
        return self

    def predict(self, X):
        return np.ones(len(X))


class TestEstimatorFitMetadata(ClassifierMixin, BaseEstimator):
    """An estimator accepting two metadata in its ``fit`` method."""

    def __init__(self, sample_weight_none=True, brand_none=True):
        self.sample_weight_none = sample_weight_none
        self.brand_none = brand_none

    def fit(self, X, y, sample_weight=None, brand=None):
        assert (
            sample_weight is None
        ) == self.sample_weight_none, "sample_weight and sample_weight_none don't agree"
        assert (brand is None) == self.brand_none, "brand and brand_none don't agree"
        return self

    def predict(self, X):
        return np.ones(len(X))


class TestSimpleMetaEstimator(MetaEstimatorMixin, ClassifierMixin, BaseEstimator):
    """A meta-estimator which also consumes sample_weight itself in ``fit``."""

    def __init__(self, estimator, sample_weight_none):
        self.sample_weight_none = sample_weight_none
        self.estimator = estimator

    def fit(self, X, y, sample_weight=None, **kwargs):
        assert (
            sample_weight is None
        ) == self.sample_weight_none, "sample_weight and sample_weight_none don't agree"

        if sample_weight is not None:
            kwargs["sample_weight"] = sample_weight
        request_routing = metadata_request_factory(self)
        request_routing.validate_metadata(params=kwargs, method="fit")
        params = request_routing.get_params(params=kwargs, method="fit")
        self.estimator_ = clone(self.estimator).fit(X, y, **params.estimator.fit)
        return self

    def get_metadata_routing(self):
        router = (
            MetadataRouter()
            .add_self(self)
            .add(estimator=self.estimator, method_mapping="fit")
        )
        return router.to_dict()


class TestTransformer(TransformerMixin, BaseEstimator):
    """A transformer which accepts metadata on fit and transform."""

    def __init__(
        self,
        brand_none=True,
        new_param_none=True,
        fit_sample_weight_none=True,
        transform_sample_weight_none=True,
    ):
        self.brand_none = brand_none
        self.new_param_none = new_param_none
        self.fit_sample_weight_none = fit_sample_weight_none
        self.transform_sample_weight_none = transform_sample_weight_none

    def fit(self, X, y=None, brand=None, new_param=None, sample_weight=None):
        assert (
            sample_weight is None
        ) == self.fit_sample_weight_none, (
            "sample_weight and fit_sample_weight_none don't agree"
        )
        assert (
            new_param is None
        ) == self.new_param_none, "new_param and new_param_none don't agree"
        assert (brand is None) == self.brand_none, "brand and brand_none don't agree"
        return self

    def transform(self, X, y=None, sample_weight=None):
        assert (
            sample_weight is None
        ) == self.transform_sample_weight_none, (
            "sample_weight and transform_sample_weight_none don't agree"
        )
        return X


class TestMetaTransformer(MetaEstimatorMixin, TransformerMixin, BaseEstimator):
    """A simple meta-transformer."""

    def __init__(self, transformer):
        self.transformer = transformer

    def fit(self, X, y=None, **fit_params):
        request_routing = metadata_request_factory(self)
        request_routing.validate_metadata(
            method="fit",
            params=fit_params,
        )
        params = request_routing.get_params(method="fit", params=fit_params)
        self.transformer_ = clone(self.transformer).fit(X, y, **params.transformer.fit)
        return self

    def transform(self, X, y=None, **transform_params):
        request_routing = metadata_request_factory(self)
        request_routing.validate_metadata(
            method="transform",
            params=transform_params,
        )
        params = request_routing.get_params(method="transform", params=transform_params)
        return self.transformer_.transform(X, **params.transformer.transform)

    def get_metadata_routing(self):
        return (
            MetadataRouter()
            .add(transformer=self.transformer, method_mapping="one-to-one")
            .to_dict()
        )


class SimplePipeline(BaseEstimator):
    """A very simple pipeline, assuming the last step is always a predictor."""

    def __init__(self, steps):
        self.steps = steps

    def fit(self, X, y, **fit_params):
        self.steps_ = []
        request_routing = metadata_request_factory(self)
        request_routing.validate_params(params=fit_params, method="fit")
        params = request_routing.get_params(params=fit_params, method="fit")
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
        request_routing = metadata_request_factory(self)
        request_routing.validate_params(params=predict_params, method="predict")
        params = request_routing.get_params(params=predict_params, method="predict")
        for i, step in enumerate(self.steps_[:-1]):
            X_transformed = step.transform(X, **params.get(f"step_{i}").transform)

        return self.steps_[-1].predict(X_transformed, **params.predictor.transform)

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
        router.add(predictor=self.steps[-1], mapping="one-to-one")
        return router.to_dict()


def test_assert_request_is_empty():
    requests = MetadataRequest()
    assert_request_is_empty(requests)
    assert_request_is_empty(requests.to_dict())

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


def test_default_requests():
    class OddEstimator(BaseEstimator):
        __metadata_request__fit = {
            # set a different default request
            "sample_weight": RequestType.REQUESTED
        }  # type: ignore

    odd_request = metadata_request_factory(OddEstimator())
    assert odd_request.fit.requests == {"sample_weight": RequestType.REQUESTED}

    # check other test estimators
    assert not len(metadata_request_factory(TestEstimatorNoMetadata()).fit.requests)
    assert_request_is_empty(TestEstimatorNoMetadata().get_metadata_routing())

    trs_request = metadata_request_factory(TestTransformer())
    assert trs_request.fit.requests == {
        "sample_weight": RequestType(None),
        "brand": RequestType(None),
        "new_param": RequestType(None),
    }
    assert trs_request.transform.requests == {
        "sample_weight": RequestType(None),
    }
    assert_request_is_empty(trs_request)

    est_request = metadata_request_factory(TestEstimatorFitMetadata())
    assert est_request.fit.requests == {
        "sample_weight": RequestType(None),
        "brand": RequestType(None),
    }
    assert_request_is_empty(est_request)


def test_simple_metadata_routing():
    # Tests that metadata is properly routed
    # The underlying estimator doesn't accept or request metadata
    cls = TestSimpleMetaEstimator(
        estimator=TestEstimatorNoMetadata(), sample_weight_none=True
    )
    cls.fit(X, y)

    # Meta-estimator consumes sample_weight, but doesn't forward it to the underlying
    # estimator
    cls = TestSimpleMetaEstimator(
        estimator=TestEstimatorNoMetadata(), sample_weight_none=False
    )
    cls.fit(X, y, sample_weight=my_weights)

    # If the estimator accepts the metadata but doesn't explicitly say it doesn't
    # need it, there's an error
    cls = TestSimpleMetaEstimator(
        estimator=TestEstimatorFitMetadata(),
        sample_weight_none=False,
    )
    with pytest.raises(
        ValueError,
        match=(
            "sample_weight is passed but is not explicitly set as requested or not. In"
            " method: fit"
        ),
    ):
        cls.fit(X, y, sample_weight=my_weights)

    # Explicitly saying the estimator doesn't need it, makes the error go away,
    # but if a metadata is passed which is not requested by any object/estimator,
    # there will be still an error
    cls = TestSimpleMetaEstimator(
        estimator=TestEstimatorFitMetadata().fit_requests(sample_weight=False),
        sample_weight_none=False,
    )
    # this doesn't raise since TestSimpleMetaEstimator itself is a consumer,
    # and passing metadata to the consumer directly is fine regardless of its
    # metadata_request values.
    cls.fit(X, y, sample_weight=my_weights)

    # Requesting a metadata will make the meta-estimator forward it correctly
    cls = TestSimpleMetaEstimator(
        estimator=TestEstimatorFitMetadata(sample_weight_none=False).fit_requests(
            sample_weight=True
        ),
        sample_weight_none=False,
    )
    cls.fit(X, y, sample_weight=my_weights)


def test_invalid_metadata():
    # check that passing wrong metadata raises an error
    trs = TestMetaTransformer(
        transformer=TestTransformer().transform_requests(sample_weight=True)
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
    trs = TestMetaTransformer(
        transformer=TestTransformer().transform_requests(sample_weight=False)
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
        "_type": "request",
        "score": {
            "my_param": RequestType(True),
            "my_other_param": RequestType(None),
            "sample_weight": RequestType(None),
        },
        "fit": {
            "my_other_param": RequestType(None),
            "sample_weight": RequestType(None),
        },
        "partial_fit": {},
        "predict": {"my_param": RequestType(True)},
        "transform": {},
        "inverse_transform": {},
        "split": {},
    }
    assert TestDefaults().get_metadata_routing() == expected

    est = TestDefaults().score_requests(my_param="other_param")
    expected = {
        "_type": "request",
        "score": {
            "my_param": "other_param",
            "my_other_param": RequestType(None),
            "sample_weight": RequestType(None),
        },
        "fit": {
            "my_other_param": RequestType(None),
            "sample_weight": RequestType(None),
        },
        "partial_fit": {},
        "predict": {"my_param": RequestType(True)},
        "transform": {},
        "inverse_transform": {},
        "split": {},
    }
    assert est.get_metadata_routing() == expected

    est = TestDefaults().fit_requests(sample_weight=True)
    expected = {
        "_type": "request",
        "score": {
            "my_param": RequestType(True),
            "my_other_param": RequestType(None),
            "sample_weight": RequestType(None),
        },
        "fit": {
            "my_other_param": RequestType(None),
            "sample_weight": RequestType(True),
        },
        "partial_fit": {},
        "predict": {"my_param": RequestType(True)},
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

    assert metadata_request_factory(ExplicitRequest()).fit.requests == {
        "prop": RequestType.ERROR_IF_PASSED
    }
    assert_request_is_empty(ExplicitRequest().get_metadata_routing(), exclude="fit")

    class ExplicitRequestOverwrite(BaseEstimator):
        __metadata_request__fit = {"prop": RequestType.REQUESTED}

        def fit(self, X, y, prop=None, **kwargs):
            return self

    assert metadata_request_factory(ExplicitRequestOverwrite()).fit.requests == {
        "prop": RequestType.REQUESTED
    }
    assert_request_is_empty(
        ExplicitRequestOverwrite().get_metadata_routing(), exclude="fit"
    )

    class ImplicitRequest(BaseEstimator):
        def fit(self, X, y, prop=None, **kwargs):
            return self

    assert metadata_request_factory(ImplicitRequest()).fit.requests == {
        "prop": RequestType.ERROR_IF_PASSED
    }
    assert_request_is_empty(ImplicitRequest().get_metadata_routing(), exclude="fit")

    class ImplicitRequestRemoval(BaseEstimator):
        __metadata_request__fit = {"prop": RequestType.UNUSED}

        def fit(self, X, y, prop=None, **kwargs):
            return self

    assert metadata_request_factory(ImplicitRequestRemoval()).fit.requests == {}
    assert_request_is_empty(ImplicitRequestRemoval().get_metadata_routing())


def test_method_metadata_request():
    mmr = MethodMetadataRequest(name="fit")

    with pytest.raises(ValueError, match="alias should be either a string or"):
        mmr.add_request(prop="foo", alias=1.4)

    mmr.add_request(prop="foo", alias=None)
    assert mmr.requests == {"foo": RequestType.ERROR_IF_PASSED}
    mmr.add_request(prop="foo", alias=False)
    mmr.add_request(prop="foo", alias=True)
    assert mmr.requests == {"foo": RequestType.REQUESTED}

    assert MethodMetadataRequest.from_dict(None, name="fit").requests == {}
    assert MethodMetadataRequest.from_dict("foo", name="fit").requests == {
        "foo": RequestType.ERROR_IF_PASSED
    }
    assert MethodMetadataRequest.from_dict(["foo", "bar"], name="fit").requests == {
        "foo": RequestType.ERROR_IF_PASSED,
        "bar": RequestType.ERROR_IF_PASSED,
    }


def test_metadata_request_factory():
    class Consumer(BaseEstimator):
        __metadata_request__fit = {"prop": RequestType.ERROR_IF_PASSED}

    assert_request_is_empty(metadata_request_factory(None))
    assert_request_is_empty(metadata_request_factory({"_type": "request"}))
    assert_request_is_empty(metadata_request_factory({"_type": "router"}))
    assert_request_is_empty(metadata_request_factory(object()))

    mr = MetadataRequest({"fit": "foo"}, default="bar")
    mr_factory = metadata_request_factory(mr)
    assert_request_is_empty(mr_factory, exclude="fit")
    assert mr_factory.fit.requests == {"foo": "bar"}

    mr = metadata_request_factory(Consumer())
    assert_request_is_empty(mr, exclude="fit")
    assert mr.fit.requests == {"prop": RequestType.ERROR_IF_PASSED}
