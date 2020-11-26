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
from sklearn.utils import _validate_required_props

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
    if exclude is None:
        exclude = []
    assert not any(
        method_request
        for method_name, method_request in metadata_request.items()
        if method_name not in exclude
    )


class MyEst(ClassifierMixin, BaseEstimator):
    def __init__(self, C=1.0):
        self._metadata_request = {"fit": ["sample_weight", "brand"]}
        self.C = C

    def fit(self, X, y, **fit_params):
        _validate_required_props(self.get_metadata_request().fit, fit_params)
        assert set(fit_params.keys()) <= set(
            [list(x)[0] for x in self.get_metadata_request().fit.values()]
        )
        self.svc_ = SVC(C=self.C).fit(X, y)
        return self

    def predict(self, X):
        return self.svc_.predict(X)

    def predict_proba(self, X):
        return self.svc_predict_proba(X)


class StuffConsumer(MetadataConsumer):
    def request_new_param(self, *, fit=True):
        self._request_key_for_method(
            method="fit", param="new_param", user_provides=fit
        )
        return self

    def request_brand(self, *, fit=True):
        self._request_key_for_method(
            method="fit", param="brand", user_provides=fit
        )
        return self


class MyTrs(
    SampleWeightConsumer, StuffConsumer, TransformerMixin, BaseEstimator
):
    def __init__(self):
        self._metadata_request = {"fit": ["sample_weight"]}

    def fit(self, X, y=None, **fit_params):
        # extract the values from the metadata requests. Since this is not a
        # meta-estimator, the values will always have exactly 1 member in the
        # corresponding set.
        req_props = [list(x)[0] for x in
                     self.get_metadata_request().fit.values()]
        assert set(fit_params.keys()) <= set(req_props)
        self._estimator = SelectKBest().fit(X, y)
        return self

    def transform(self, X, y=None):
        return self._estimator.transform(X)


def my_metric(y, y_pred, new_param):
    return balanced_accuracy_score(y, y_pred)


def test_defaults():
    assert_request_is_empty(LogisticRegression().get_metadata_request())
    # check default requests for dummy estimators
    trs = MyTrs()
    trs_request = _standardize_metadata_request(trs.get_metadata_request())
    assert trs_request.fit == {"sample_weight": {"sample_weight"}}
    assert_request_is_empty(trs_request, exclude={"fit"})

    est = MyEst()
    est_request = _standardize_metadata_request(est.get_metadata_request())
    assert est_request.fit == {
        "sample_weight": {"sample_weight"},
        "brand": {"brand"},
    }
    assert_request_is_empty(est_request, exclude={"fit"})


def test_pipeline():
    X, y = make_classification()
    sw = np.random.rand(len(X))
    my_data = [5, 6]
    brand = ["my brand"]

    clf = make_pipeline(MyTrs(), MyEst())
    clf.fit(X, y, sample_weight=sw, brand=brand)
    with pytest.raises(
        ValueError, match="Requested properties are.*other_param"
    ):
        clf.fit(X, y, sample_weight=sw, brand=brand, other_param=sw)

    trs = MyTrs().request_new_param(fit="my_sw")
    trs.request_sample_weight(fit=False)

    trs_request = _standardize_metadata_request(trs.get_metadata_request())
    assert trs_request.fit == {"my_sw": {"new_param"}}
    assert_request_is_empty(trs_request, exclude={"fit"})

    clf = make_pipeline(trs, MyEst())
    pipe_request = _standardize_metadata_request(clf.get_metadata_request())
    print(pipe_request.fit)
    assert pipe_request.fit == {
        "my_sw": {"my_sw"},
        "sample_weight": {"sample_weight"},
        "brand": {"brand"},
    }
    assert_request_is_empty(pipe_request, exclude={"fit"})
    clf.fit(X, y, sample_weight=sw, brand=brand, my_sw=my_data)

    # TODO: assert that trs did *not* receive sample_weight, but did receive
    # my_sw

    # If requested metadata is not given, no warning or error is raised
    with pytest.warns(None) as record:
        clf.fit(X, y, brand=brand)
    assert not record.list

    scorer = make_scorer(my_metric, request_props="new_param")

    param_grid = {"myest__C": [0.1, 1]}

    print("@" * 150 + " GS")
    gs = GridSearchCV(clf, param_grid=param_grid, scoring=scorer)
    print("GS props request: ", gs.get_metadata_request())
    gs.fit(X, y, new_param=brand, sample_weight=sw, my_sw=sw, brand=brand)


def test_slep_caseA():
    # Case A: weighted scoring and fitting

    # Here we presume that GroupKFold requests `groups` by default.
    # We need to explicitly request weights in make_scorer and for
    # LogisticRegressionCV. Both of these consumers understand the meaning
    # of the key "sample_weight".

    weighted_acc = make_scorer(accuracy_score, request_props=["sample_weight"])
    lr = LogisticRegressionCV(
        cv=GroupKFold(), scoring=weighted_acc,
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
    lr = LogisticRegressionCV(cv=GroupKFold(), scoring=weighted_acc,)
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
        cv=GroupKFold(), scoring=weighted_acc,
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
        cv=GroupKFold(), scoring=weighted_acc,
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
