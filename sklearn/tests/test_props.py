import re
import numpy as np
import pytest

from sklearn.base import BaseEstimator, ClassifierMixin, TransformerMixin
from sklearn.datasets import make_classification
from sklearn.feature_selection import SelectKBest
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import make_scorer
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_validate
from sklearn.model_selection._split import check_cv
from sklearn.pipeline import make_pipeline
from sklearn.svm import SVC
from sklearn.utils import MetadataRequest
from sklearn.utils.metadata_requests import RequestType
from sklearn.utils.metadata_requests import metadata_request_factory
from sklearn.utils.metadata_requests import MetadataRouter
from sklearn.utils.metadata_requests import MethodMetadataRequest
from sklearn.model_selection import KFold
from sklearn.model_selection import GroupKFold
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import TimeSeriesSplit
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import LeavePGroupsOut
from sklearn.model_selection import LeavePOut
from sklearn.model_selection import RepeatedKFold
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.model_selection import ShuffleSplit
from sklearn.model_selection import GroupShuffleSplit
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import PredefinedSplit

from sklearn.base import _MetadataRequester

NonGroupCVs = [
    KFold,
    LeaveOneOut,
    LeavePOut,
    RepeatedStratifiedKFold,
    RepeatedKFold,
    ShuffleSplit,
    StratifiedKFold,
    StratifiedShuffleSplit,
    PredefinedSplit,
    TimeSeriesSplit,
]

GroupCVs = [
    GroupKFold,
    LeaveOneGroupOut,
    LeavePGroupsOut,
    GroupShuffleSplit,
]


N, M = 100, 4
X = np.random.rand(N, M)
y = np.random.randint(0, 2, size=N)
my_groups = np.random.randint(0, 10, size=N)
my_weights = np.random.rand(N)
my_other_weights = np.random.rand(N)


def assert_request_is_empty(metadata_request, exclude=None):
    if isinstance(metadata_request, MetadataRequest):
        metadata_request = metadata_request.to_dict()
    if exclude is None:
        exclude = []
    for method, request in metadata_request.items():
        if method in exclude:
            continue
        props = [
            prop
            for prop, alias in request.items()
            if isinstance(alias, str)
            or RequestType(alias) != RequestType.ERROR_IF_PASSED
        ]
        assert not len(props)


class MyEst(ClassifierMixin, BaseEstimator):
    _metadata_request__sample_weight = {
        "fit": {"sample_weight": RequestType.REQUESTED}  # type: ignore
    }
    _metadata_request__brand = {"fit": {"brand": RequestType.REQUESTED}}

    def __init__(self, C=1.0):
        self.C = C

    def fit(self, X, y, **fit_params):
        metadata_request_factory(self).fit.validate_metadata(
            ignore_extras=False, kwargs=fit_params
        )
        self.svc_ = SVC(C=self.C).fit(X, y)
        return self

    def predict(self, X):
        return self.svc_.predict(X)


class MyTrs(TransformerMixin, BaseEstimator):
    def fit(self, X, y=None, brand=None, new_param=None, sample_weight=None):
        metadata_request_factory(self).fit.validate_metadata(
            ignore_extras=False,
            kwargs={
                "brand": brand,
                "new_param": new_param,
                "sample_weight": sample_weight,
            },
        )
        self._estimator = SelectKBest().fit(X, y)
        return self

    def transform(self, X, y=None):
        return self._estimator.transform(X)


def my_metric(y, y_pred, new_param):
    return balanced_accuracy_score(y, y_pred)


def test_defaults():
    assert_request_is_empty(LogisticRegression().get_metadata_request())
    # check default requests for dummy estimators
    trs_request = metadata_request_factory(MyTrs())
    assert trs_request.fit.requests == {
        "sample_weight": RequestType(None),
        "brand": RequestType(None),
        "new_param": RequestType(None),
    }
    assert_request_is_empty(trs_request)

    est_request = metadata_request_factory(MyEst())
    assert est_request.fit.requests == {
        "sample_weight": RequestType(True),
        "brand": RequestType(True),
    }
    assert_request_is_empty(est_request, exclude={"fit"})


def test_pipeline():
    X, y = make_classification()
    sw = np.random.rand(len(X))
    my_data = [5, 6]
    brand = ["my brand"]

    # MyEst is requesting "brand" but MyTrs has it as ERROR_IF_PASSED
    with pytest.raises(
        ValueError,
        match="brand is passed but is not explicitly set as requested or not.",
    ):
        clf = make_pipeline(MyTrs(), MyEst())
        clf.fit(X, y, sample_weight=sw, brand=brand)

    clf = make_pipeline(MyTrs().fit_requests(brand=True, sample_weight=True), MyEst())
    clf.fit(X, y, sample_weight=sw, brand=brand)

    with pytest.raises(ValueError, match="Metadata passed which is not understood"):
        clf.fit(X, y, sample_weight=sw, brand=brand, other_param=sw)

    trs = MyTrs().fit_requests(new_param="my_sw", sample_weight=False)

    trs_request = metadata_request_factory(trs)
    assert trs_request.fit.requests == {
        "new_param": "my_sw",
        "brand": RequestType.ERROR_IF_PASSED,
        "sample_weight": RequestType.UNREQUESTED,
    }
    assert_request_is_empty(trs_request, exclude={"fit"})

    clf = make_pipeline(trs, MyEst())
    clf.get_metadata_request()

    clf = make_pipeline(MyTrs().fit_requests(new_param="my_sw"), MyEst())
    pipe_request = metadata_request_factory(clf)
    assert pipe_request.fit.requests == {
        "my_sw": RequestType.REQUESTED,
        "sample_weight": RequestType.REQUESTED,
        "brand": RequestType.REQUESTED,
    }
    assert_request_is_empty(pipe_request, exclude={"fit"})
    with pytest.raises(
        ValueError,
        match=(
            "Error while validating fit parameters for mytrs. The underlying "
            "error message: brand is passed but is not explicitly set as "
            "requested or not."
        ),
    ):
        clf.fit(X, y, sample_weight=sw, brand=brand, my_sw=my_data)

    clf.named_steps["mytrs"].fit_requests(brand=True, sample_weight=False)
    clf.fit(X, y, sample_weight=sw, brand=brand, my_sw=my_data)

    # TODO: assert that trs did *not* receive sample_weight, but did receive
    # my_sw

    # If requested metadata is not given, no warning or error is raised
    with pytest.warns(None) as record:
        clf.fit(X, y, brand=brand)
    assert not record.list

    scorer = make_scorer(my_metric, score_params="new_param")

    param_grid = {"myest__C": [0.1, 1]}

    gs = GridSearchCV(clf, param_grid=param_grid, scoring=scorer)
    gs.get_metadata_request()
    gs.fit(X, y, new_param=brand, sample_weight=sw, my_sw=sw, brand=brand)


def test_slep_caseA():
    # Case A: weighted scoring and fitting

    # Here we presume that GroupKFold requests `groups` by default.
    # We need to explicitly request weights in make_scorer and for
    # LogisticRegressionCV. Both of these consumers understand the meaning
    # of the key "sample_weight".

    weighted_acc = make_scorer(accuracy_score, score_params=["sample_weight"])
    lr = LogisticRegressionCV(
        cv=GroupKFold(),
        scoring=weighted_acc,
    ).fit_requests(sample_weight=True)
    lr.fit(X, y, sample_weight=my_weights, groups=my_groups)
    cross_validate(
        lr,
        X,
        y,
        cv=GroupKFold(),
        props={"sample_weight": my_weights, "groups": my_groups},
        scoring=weighted_acc,
    )

    # Error handling: if props={'sample_eight': my_weights, ...} was passed,
    # cross_validate would raise an error, since 'sample_eight' was not
    # requested by any of its children.


def test_slep_caseB():
    # Case B: weighted scoring and unweighted fitting

    # Since LogisticRegressionCV requires that weights explicitly be requested,
    # removing that request means the fitting is unweighted.

    weighted_acc = make_scorer(accuracy_score, score_params=["sample_weight"])
    lr = LogisticRegressionCV(
        cv=GroupKFold(),
        scoring=weighted_acc,
    )
    with pytest.raises(
        ValueError,
        match="sample_weight is passed but is not explicitly set as requested or not",
    ):
        cross_validate(
            lr,
            X,
            y,
            cv=GroupKFold(),
            props={"sample_weight": my_weights, "groups": my_groups},
            scoring=weighted_acc,
        )

    lr.fit_requests(sample_weight=False)
    cross_validate(
        lr,
        X,
        y,
        cv=GroupKFold(),
        props={"sample_weight": my_weights, "groups": my_groups},
        scoring=weighted_acc,
    )


def test_slep_caseC():
    # Case C: unweighted feature selection

    # Like LogisticRegressionCV, SelectKBest needs to request weights
    # explicitly. Here it does not request them.

    weighted_acc = make_scorer(accuracy_score, score_params=["sample_weight"])
    lr = LogisticRegressionCV(
        cv=GroupKFold(),
        scoring=weighted_acc,
    ).fit_requests(sample_weight=True)
    sel = SelectKBest(k=2)
    pipe = make_pipeline(sel, lr)
    cross_validate(
        pipe,
        X,
        y,
        cv=GroupKFold(),
        props={"sample_weight": my_weights, "groups": my_groups},
        scoring=weighted_acc,
    )


def test_slep_caseD():
    # Case D: different scoring and fitting weights

    # Despite make_scorer and LogisticRegressionCV both expecting a key
    # sample_weight, we can use aliases to pass different weights to different
    # consumers.

    weighted_acc = make_scorer(
        accuracy_score, score_params={"sample_weight": "scoring_weight"}
    )
    lr = LogisticRegressionCV(
        cv=GroupKFold(),
        scoring=weighted_acc,
    ).fit_requests(sample_weight="fitting_weight")
    cross_validate(
        lr,
        X,
        y,
        cv=GroupKFold(),
        props={
            "scoring_weight": my_weights,
            "fitting_weight": my_other_weights,
            "groups": my_groups,
        },
        scoring=weighted_acc,
    )


@pytest.mark.parametrize("Klass", GroupCVs)
def test_group_splitter_metadata_requests(Klass):
    if Klass is LeavePGroupsOut:
        cv = Klass(n_groups=2)
    else:
        cv = Klass()
    # check the default metadata_request
    assert metadata_request_factory(cv).split.requests == {
        "groups": RequestType.REQUESTED
    }

    # test that setting split to False empties the metadata_request
    cv.split_requests(groups=None)
    assert_request_is_empty(cv.get_metadata_request())

    # set a different input name and test
    cv.split_requests(groups="my_groups")
    assert metadata_request_factory(cv).split.requests == {"groups": "my_groups"}


@pytest.mark.parametrize("Klass", NonGroupCVs)
def test_nongroup_splitter_metadata_requests(Klass):
    if Klass is LeavePOut:
        cv = Klass(p=2)
    elif Klass is PredefinedSplit:
        cv = Klass(test_fold=[1, 1, 0])
    else:
        cv = Klass()

    # check the default metadata_request
    assert_request_is_empty(cv.get_metadata_request())

    # test that setting split to False empties the metadata_request
    assert not hasattr(cv, "request_groups")


def test_invalid_arg_given():
    # tests that passing an invalid argument would raise an error
    weighted_acc = make_scorer(accuracy_score, score_params=["sample_weight"])
    model = LogisticRegression().fit_requests(sample_weight=True)
    param_grid = {"C": [0.1, 1]}
    gs = GridSearchCV(
        estimator=model,
        cv=GroupKFold(),
        scoring=weighted_acc,
        param_grid=param_grid,
    )
    gs.get_metadata_request()
    gs.fit(X, y, sample_weight=my_weights, groups=my_groups)
    with pytest.raises(ValueError, match="Metadata passed which is not understood"):
        gs.fit(
            X,
            y,
            sample_weigh=my_weights,
            groups=my_groups,
            sample_weight=my_weights,
        )

    with pytest.raises(ValueError, match="Metadata passed which is not understood"):
        gs.fit(
            X,
            y,
            sample_weigh=my_weights,
            groups=my_groups,
        )


def test_get_metadata_request():
    class TestDefaultsBadMetadataName(_MetadataRequester):
        _metadata_request__sample_weight = {
            "fit": "sample_weight",
            "score": "sample_weight",
        }

        _metadata_request__my_param = {
            "score": {"my_param": True},
            # the following method raise an error
            "other_method": {"my_param": True},
        }

        _metadata_request__my_other_param = {
            "score": "my_other_param",
            # this should raise since the name is different than the metadata
            "fit": "my_param",
        }

    class TestDefaultsBadMethodName(_MetadataRequester):
        _metadata_request__sample_weight = {
            "fit": "sample_weight",
            "score": "sample_weight",
        }

        _metadata_request__my_param = {
            "score": {"my_param": True},
            # the following method raise an error
            "other_method": {"my_param": True},
        }

        _metadata_request__my_other_param = {
            "score": "my_other_param",
            "fit": "my_other_param",
        }

    class TestDefaults(_MetadataRequester):
        _metadata_request__sample_weight = {
            "fit": "sample_weight",
            "score": "sample_weight",
        }

        _metadata_request__my_param = {
            "score": {"my_param": True},
            "predict": {"my_param": True},
        }

        _metadata_request__my_other_param = {
            "score": "my_other_param",
            "fit": "my_other_param",
        }

    with pytest.raises(ValueError, match="Expected all metadata to be called"):
        TestDefaultsBadMetadataName().get_metadata_request()

    with pytest.raises(ValueError, match="other_method is not supported as a method"):
        TestDefaultsBadMethodName().get_metadata_request()

    expected = {
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
    assert TestDefaults().get_metadata_request() == expected

    est = TestDefaults().score_requests(my_param="other_param")
    expected = {
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
    assert est.get_metadata_request() == expected

    est = TestDefaults().fit_requests(sample_weight=True)
    expected = {
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
    assert est.get_metadata_request() == expected


def test__get_default_requests():
    class ExplicitRequest(BaseEstimator):
        _metadata_request__prop = {"fit": "prop"}

        def fit(self, X, y):
            return self

    assert metadata_request_factory(ExplicitRequest()).fit.requests == {
        "prop": RequestType.ERROR_IF_PASSED
    }
    assert_request_is_empty(ExplicitRequest().get_metadata_request(), exclude="fit")

    class ExplicitRequestOverwrite(BaseEstimator):
        _metadata_request__prop = {"fit": {"prop": RequestType.REQUESTED}}

        def fit(self, X, y, prop=None, **kwargs):
            return self

    assert metadata_request_factory(ExplicitRequestOverwrite()).fit.requests == {
        "prop": RequestType.REQUESTED
    }
    assert_request_is_empty(
        ExplicitRequestOverwrite().get_metadata_request(), exclude="fit"
    )

    class ImplicitRequest(BaseEstimator):
        def fit(self, X, y, prop=None, **kwargs):
            return self

    assert metadata_request_factory(ImplicitRequest()).fit.requests == {
        "prop": RequestType.ERROR_IF_PASSED
    }
    assert_request_is_empty(ImplicitRequest().get_metadata_request(), exclude="fit")

    class ImplicitRequestRemoval(BaseEstimator):
        _metadata_request__prop = {"fit": {"prop": RequestType.UNUSED}}

        def fit(self, X, y, prop=None, **kwargs):
            return self

    assert metadata_request_factory(ImplicitRequestRemoval()).fit.requests == {}
    assert_request_is_empty(ImplicitRequestRemoval().get_metadata_request())


def test_validate():
    class ConsumerRouter(BaseEstimator):
        def __init__(self, cv=None):
            self.cv = cv

        def fit(self, X, y, sample_weight=None, **kwargs):
            kwargs["sample_weight"] = sample_weight
            metadata_request_factory(self).fit.validate_metadata(
                ignore_extras=False,
                self_metadata=super(),
                kwargs=kwargs,
            )
            return self

        def get_metadata_request(self):
            router = (
                MetadataRouter()
                .add(super(), mapping="one-to-one", overwrite=False, mask=False)
                .add(
                    check_cv(self.cv),
                    mapping={"fit": "split"},
                    overwrite=False,
                    mask=True,
                )
            )
            return router.get_metadata_request()

    err_message = "Metadata passed which is not understood: {param}. In method: fit"

    est = ConsumerRouter()
    est.fit(X=None, y=None)
    est.fit(X=None, y=None, sample_weight="test")

    with pytest.raises(
        ValueError, match=re.escape(err_message.format(param=["my_weight"]))
    ):
        est.fit(X=None, y=None, my_weight="test")

    with pytest.raises(
        ValueError, match=re.escape(err_message.format(param=["my_weight"]))
    ):
        est.fit(X=None, y=None, sample_weight="test", my_weight="test")

    est = ConsumerRouter(cv=GroupKFold())
    est.fit(X=None, y=None, groups="test")
    est.fit(X=None, y=None, sample_weight="test", groups="test")

    with pytest.raises(
        ValueError, match=re.escape(err_message.format(param=["my_weight"]))
    ):
        est.fit(X=None, y=None, my_weight="test")

    with pytest.raises(
        ValueError, match=re.escape(err_message.format(param=["my_weight"]))
    ):
        est.fit(X=None, y=None, sample_weight="test", my_weight="test")

    est = ConsumerRouter(cv=GroupKFold().split_requests(groups="my_groups"))
    est.fit(X=None, y=None, sample_weight="test", my_groups="test")

    with pytest.raises(
        ValueError, match=re.escape(err_message.format(param=["groups"]))
    ):
        est.fit(X=None, y=None, groups="test")

    with pytest.raises(
        ValueError, match=re.escape(err_message.format(param=["groups"]))
    ):
        est.fit(X=None, y=None, sample_weight="test", my_groups="test", groups="test")

    with pytest.raises(
        ValueError, match=re.escape(err_message.format(param=["groups"]))
    ):
        est.fit(X=None, y=None, sample_weight="test", groups="test")

    est = ConsumerRouter(
        cv=GroupKFold().split_requests(groups="my_groups")
    ).fit_requests(sample_weight="my_weight")
    est.fit(X=None, y=None, sample_weight="test", my_groups="test")
    est.fit(X=None, y=None, my_groups="test")

    with pytest.raises(
        ValueError, match=re.escape(err_message.format(param=["groups"]))
    ):
        est.fit(X=None, y=None, groups="test")

    with pytest.raises(
        ValueError, match=re.escape(err_message.format(param=["groups"]))
    ):
        est.fit(X=None, y=None, sample_weight="test", my_groups="test", groups="test")

    with pytest.raises(
        ValueError, match=re.escape(err_message.format(param=["groups"]))
    ):
        est.fit(X=None, y=None, sample_weight="test", groups="test")

    with pytest.raises(
        ValueError, match=re.escape(err_message.format(param=["my_weights"]))
    ):
        est.fit(X=None, y=None, my_weights="test", my_groups="test")

    with pytest.raises(
        ValueError, match=re.escape(err_message.format(param=["groups", "my_weights"]))
    ):
        est.fit(X=None, y=None, my_weights="test", groups="test")


def test_method_metadata_request():
    mmr = MethodMetadataRequest(name="fit")
    with pytest.raises(
        ValueError,
        match="overwrite can only be one of {True, False, 'smart', 'ignore'}.",
    ):
        mmr.add_request(prop="test", alias=None, overwrite="test")

    with pytest.raises(ValueError, match="Expected all metadata to be called test"):
        mmr.add_request(prop="foo", alias="bar", expected_metadata="test")

    with pytest.raises(ValueError, match="Aliasing is not allowed"):
        mmr.add_request(prop="foo", alias="bar", allow_aliasing=False)

    with pytest.raises(ValueError, match="alias should be either a string or"):
        mmr.add_request(prop="foo", alias=1.4)

    mmr.add_request(prop="foo", alias=None)
    assert mmr.requests == {"foo": RequestType.ERROR_IF_PASSED}
    with pytest.raises(ValueError, match="foo is already requested"):
        mmr.add_request(prop="foo", alias=True)
    with pytest.raises(ValueError, match="foo is already requested"):
        mmr.add_request(prop="foo", alias=True)
    mmr.add_request(prop="foo", alias=True, overwrite="smart")
    assert mmr.requests == {"foo": RequestType.REQUESTED}

    with pytest.raises(ValueError, match="Can only add another MethodMetadataRequest"):
        mmr.merge_method_request({})

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
        _metadata_request__prop = {"fit": "prop"}

    assert_request_is_empty(metadata_request_factory(None))
    assert_request_is_empty(metadata_request_factory({}))
    assert_request_is_empty(metadata_request_factory(object()))

    mr = MetadataRequest({"fit": "foo"}, default="bar")
    mr_factory = metadata_request_factory(mr)
    assert_request_is_empty(mr_factory, exclude="fit")
    assert mr_factory.fit.requests == {"foo": "bar"}

    mr = metadata_request_factory(Consumer())
    assert_request_is_empty(mr, exclude="fit")
    assert mr.fit.requests == {"prop": RequestType.ERROR_IF_PASSED}
