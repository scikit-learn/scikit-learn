import numpy as np
import pytest

from sklearn.base import BaseEstimator, ClassifierMixin, TransformerMixin
from sklearn.base import MetadataConsumer, SampleWeightConsumer
from sklearn.datasets import make_classification
from sklearn.feature_selection import SelectKBest
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import make_scorer
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_validate
from sklearn.pipeline import make_pipeline
from sklearn.svm import SVC
from sklearn.utils import _standardize_metadata_request
from sklearn.utils import MetadataRequest
from sklearn.utils.metadata_requests import RequestType
from sklearn.utils import build_method_metadata_params
from sklearn.utils import build_router_metadata_request

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
        self.get_metadata_request(
            output="MetadataRequest"
        ).fit.validate_metadata(ignore_extras=False, **fit_params)
        self.svc_ = SVC(C=self.C).fit(X, y)
        return self

    def predict(self, X):
        return self.svc_.predict(X)

    def predict_proba(self, X):
        return self.svc_predict_proba(X)


class StuffConsumer(MetadataConsumer):
    _metadata_request__new_param = {"fit": "new_param"}
    _metadata_request__brand = {"fit": "brand"}


class MyTrs(
    SampleWeightConsumer, StuffConsumer, TransformerMixin, BaseEstimator
):
    def fit(self, X, y=None, **fit_params):
        self.get_metadata_request(
            output="MetadataRequest"
        ).fit.validate_metadata(ignore_extras=False, **fit_params)
        self._estimator = SelectKBest().fit(X, y)
        return self

    def transform(self, X, y=None):
        return self._estimator.transform(X)


def my_metric(y, y_pred, new_param):
    return balanced_accuracy_score(y, y_pred)


def test_defaults():
    assert_request_is_empty(LogisticRegression().get_metadata_request())
    # check default requests for dummy estimators
    trs_request = MyTrs().get_metadata_request(output="MetadataRequest")
    assert trs_request.fit.requests == {
        "sample_weight": RequestType(None),
        "brand": RequestType(None),
        "new_param": RequestType(None),
    }
    assert_request_is_empty(trs_request)

    est_request = MyEst().get_metadata_request(output="MetadataRequest")
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
        match=(
            "brand is passed but is not explicitly set as requested or not."
        ),
    ):
        clf = make_pipeline(MyTrs(), MyEst())
        clf.fit(X, y, sample_weight=sw, brand=brand)

    clf = make_pipeline(
        MyTrs().fit_requests(brand=True, sample_weight=True), MyEst()
    )
    clf.fit(X, y, sample_weight=sw, brand=brand)

    with pytest.raises(
        ValueError, match="Metadata passed which is not understood"
    ):
        clf.fit(X, y, sample_weight=sw, brand=brand, other_param=sw)

    trs = MyTrs().fit_requests(new_param="my_sw", sample_weight=False)

    trs_request = trs.get_metadata_request(output="MetadataRequest")
    assert trs_request.fit.requests == {
        "new_param": "my_sw",
        "brand": RequestType.ERROR_IF_PASSED,
        "sample_weight": RequestType.UNREQUESTED,
    }
    assert_request_is_empty(trs_request, exclude={"fit"})

    clf = make_pipeline(trs, MyEst())
    pipe_request = clf.get_metadata_request(output="MetadataRequest")
    assert pipe_request.fit.to_dict() == {
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

    clf.named_steps["mytrs"].fit_requests(brand=True)
    clf.fit(X, y, sample_weight=sw, brand=brand, my_sw=my_data)

    # TODO: assert that trs did *not* receive sample_weight, but did receive
    # my_sw

    # If requested metadata is not given, no warning or error is raised
    with pytest.warns(None) as record:
        clf.fit(X, y, brand=brand)
    assert not record.list

    scorer = make_scorer(my_metric, request_props="new_param")

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

    weighted_acc = make_scorer(accuracy_score, request_props=["sample_weight"])
    lr = LogisticRegressionCV(
        cv=GroupKFold(),
        scoring=weighted_acc,
    ).request_sample_weight(fit="sample_weight")
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

    weighted_acc = make_scorer(accuracy_score, request_props=["sample_weight"])
    lr = LogisticRegressionCV(
        cv=GroupKFold(),
        scoring=weighted_acc,
    )
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

    weighted_acc = make_scorer(accuracy_score, request_props=["sample_weight"])
    lr = LogisticRegressionCV(
        cv=GroupKFold(),
        scoring=weighted_acc,
    ).request_sample_weight(fit=True)
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
        accuracy_score, request_props={"scoring_weight": "sample_weight"}
    )
    lr = LogisticRegressionCV(
        cv=GroupKFold(),
        scoring=weighted_acc,
    ).request_sample_weight(fit="fitting_weight")
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
    assert cv.get_metadata_request() == _standardize_metadata_request(
        {"split": ["groups"]}
    )

    # test that setting split to False empties the metadata_request
    cv.request_groups(split=False)
    assert_request_is_empty(cv.get_metadata_request())

    # set a different input name and test
    cv.request_groups(split="my_groups")
    assert cv.get_metadata_request() == _standardize_metadata_request(
        {"split": {"my_groups": "groups"}}
    )


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
    weighted_acc = make_scorer(accuracy_score, request_props=["sample_weight"])
    model = LogisticRegression().request_sample_weight(fit=True)
    param_grid = {"C": [0.1, 1]}
    gs = GridSearchCV(
        estimator=model,
        cv=GroupKFold(),
        scoring=weighted_acc,
        param_grid=param_grid,
    )
    gs.fit(X, y, sample_weight=my_weights, groups=my_groups)
    with pytest.raises(ValueError, match="Requested properties are"):
        gs.fit(
            X,
            y,
            sample_weigh=my_weights,
            groups=my_groups,
            sample_weight=my_weights,
        )

    with pytest.raises(ValueError, match="Requested properties are"):
        gs.fit(
            X,
            y,
            sample_weigh=my_weights,
            groups=my_groups,
        )


def test_get_metadata_request():
    class TestDefaultsBadMetadataName(SampleWeightConsumer, _MetadataRequester):
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

    class TestDefaultsBadMethodName(SampleWeightConsumer, _MetadataRequester):
        _metadata_request__my_param = {
            "score": {"my_param": True},
            # the following method raise an error
            "other_method": {"my_param": True},
        }

        _metadata_request__my_other_param = {
            "score": "my_other_param",
            "fit": "my_other_param",
        }

    class TestDefaults(_MetadataRequester, SampleWeightConsumer):
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

    with pytest.raises(
        ValueError, match="other_method is not supported as a method"
    ):
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


def test_build_router_metadata_request():
    cv = GroupKFold().request_groups(split="my_groups")
    est = LogisticRegression().request_sample_weight(fit=True, score=True)
    scorer = make_scorer(
        accuracy_score, request_props={"my_weights": "sample_weight"}
    )
    got = build_router_metadata_request(
        children={"base": est, "cv": cv, "scorer": scorer},
        routing=[
            ("base", "all", "all"),
            ("cv", "split", "split"),
            ("scorer", "fit", "score"),
            ("scorer", "score", "score"),
        ],
    )
    expected = {
        "score": {
            "my_weights": {"my_weights"},
            "sample_weight": {"sample_weight"},
        },
        "fit": {
            "sample_weight": {"sample_weight"},
            "my_weights": {"my_weights"},
        },
        "partial_fit": {},
        "predict": {},
        "transform": {},
        "inverse_transform": {},
        "split": {"my_groups": {"my_groups"}},
    }
    assert expected == got


def test_build_method_metadata_params():
    cv = GroupKFold().request_groups(split="my_groups")
    est = LogisticRegression().request_sample_weight(fit=True, score=True)
    scorer = make_scorer(
        accuracy_score, request_props={"my_weights": "sample_weight"}
    )

    with pytest.raises(ValueError, match="Conflicting parameters"):
        build_method_metadata_params(
            children={"base": est, "cv": cv, "scorer": scorer},
            routing=[
                ("base", "all", "all"),
                ("cv", "split", "split"),
                ("scorer", "fit", "score"),
                ("scorer", "score", "score"),
            ],
            metadata={
                "my_weights": [1],
                "sample_weight": [2],
                "my_groups": [3],
            },
        )

    # in fit
    got = build_method_metadata_params(
        children={"base": est, "scorer": scorer},
        routing=[
            ("base", "fit", "fit"),
            ("scorer", "score", "score"),
        ],
        metadata={"my_weights": [1], "sample_weight": [2], "my_groups": [3]},
    )
    expected = {
        "score": {
            "sample_weight": [1],
        },
        "fit": {
            "sample_weight": [2],
        },
        "partial_fit": {},
        "predict": {},
        "transform": {},
        "inverse_transform": {},
        "split": {},
    }
    assert expected == got
