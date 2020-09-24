import numpy as np
import pytest

from sklearn.feature_selection import SelectKBest
from sklearn.pipeline import make_pipeline
from sklearn.base import BaseEstimator, ClassifierMixin, TransformerMixin
from sklearn.utils.validation import _validate_required_props
from sklearn.datasets import make_classification
from sklearn.metrics import make_scorer, balanced_accuracy_score
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.utils import _standardize_metadata_request


def assert_request_is_empty(metadata_request, exclude=None):
    if exclude is None:
        exclude = []
    assert not any(method_request
                   for method_name, method_request in metadata_request.items()
                   if method_name not in exclude)


class MyEst(ClassifierMixin, BaseEstimator):
    def __init__(self, C=1.0):
        self._metadata_request = {'fit': ['sample_weight', 'brand']}
        self.C = C

    def fit(self, X, y, **fit_params):
        _validate_required_props(self.get_metadata_request().fit, fit_params)
        assert set(fit_params.keys()) <= \
            set(self.get_metadata_request().fit.values())
        self.svc_ = SVC(C=self.C).fit(X, y)
        return self

    def predict(self, X):
        return self.svc_.predict(X)

    def predict_proba(self, X):
        return self.svc_predict_proba(X)


class MyTrs(TransformerMixin, BaseEstimator):
    def __init__(self):
        self._metadata_request = {'fit': ['sample_weight']}

    def fit(self, X, y=None, **fit_params):
        req_props = self.get_metadata_request().fit
        _validate_required_props(req_props, fit_params)
        self._estimator = SelectKBest().fit(X, y)
        assert set(fit_params.keys()) <= set(req_props.values())
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
    assert est_request.fit == {"sample_weight": {"sample_weight"},
                               "brand": {"brand"}}
    assert_request_is_empty(est_request, exclude={"fit"})


def test_pipeline():
    X, y = make_classification()
    sw = np.random.rand(len(X))
    my_data = [5, 6]
    brand = ['my brand']

    clf = make_pipeline(MyTrs(), MyEst())
    clf.fit(X, y, sample_weight=sw, brand=brand)
    with pytest.raises(ValueError,
                       match="Requested properties are.*other_param"):
        clf.fit(X, y, sample_weight=sw, brand=brand, other_param=sw)

    trs = MyTrs().set_metadata_request(
        {'fit': {'my_sw': 'new_param'}})

    trs_request = _standardize_metadata_request(trs.get_metadata_request())
    assert trs_request.fit == {"my_sw": {"new_param"}}
    assert_request_is_empty(trs_request, exclude={"fit"})

    clf = make_pipeline(trs, MyEst())
    pipe_request = _standardize_metadata_request(clf.get_metadata_request())
    print(pipe_request.fit)
    assert pipe_request.fit == {'my_sw': {'my_sw'},
                                'sample_weight': {'sample_weight'},
                                'brand': {'brand'}}
    assert_request_is_empty(pipe_request, exclude={"fit"})
    clf.fit(X, y, sample_weight=sw, brand=brand, my_sw=my_data)

    # TODO: assert that trs did *not* receive sample_weight, but did receive
    # my_sw

    # If requested metadata is not given, no warning or error is raised
    with pytest.warns(None) as record:
        clf.fit(X, y, brand=brand)
    assert not record.list

    scorer = make_scorer(my_scorer, request_props="new_param")

    param_grid = {'myest__C': [0.1, 1]}

    print("@" * 150 + " GS")
    gs = GridSearchCV(clf, param_grid=param_grid, scoring=scorer)
    print("GS props request: ", gs.get_metadata_request())
    gs.fit(X, y, new_param=brand, sample_weight=sw, my_sw=sw, brand=brand)
