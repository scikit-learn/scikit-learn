import numpy as np
import pytest

from sklearn.feature_selection import SelectKBest
from sklearn.pipeline import make_pipeline
from sklearn.base import BaseEstimator, ClassifierMixin, TransformerMixin
from sklearn.utils.validation import _validate_required_props
from sklearn.datasets import make_classification
from sklearn.metrics import get_scorer, make_scorer
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC


class MyEst(ClassifierMixin, BaseEstimator):
    def __init__(self, C=1.0):
        self._props_request = {'fit': ['sample_weight', 'brand']}
        self.C = C

    def fit(self, X, y, **fit_params):
        _validate_required_props(self.get_props_request().fit, fit_params)
        assert set(fit_params.keys()) == \
            set(self.get_props_request().fit.values())
        self.svc_ = SVC(C=self.C).fit(X, y)
        return self

    def predict(self, X):
        return self.svc_.predict(X)

    def predict_proba(self, X):
        return self.svc_predict_proba(X)


class MyTrs(TransformerMixin, BaseEstimator):
    def __init__(self):
        self._props_request = {'fit': ['sample_weight']}

    def fit(self, X, y=None, **fit_params):
        req_props = self.get_props_request().fit
        _validate_required_props(req_props, fit_params)
        self._estimator = SelectKBest().fit(X, y)
        assert set(fit_params.keys()) == set(req_props.values())
        return self

    def transform(self, X, y=None):
        return self._estimator.transform(X)


def my_scorer(y, y_pred, **kwargs):
    assert kwargs.keys() == ["new_param"]
    return get_scorer("balanced_accuracy")(y, y_pred)


def test_pipeline():
    X, y = make_classification()
    sw = np.random.rand(len(X))
    my_data = [5, 6]
    brand = ['my brand']

    #clf = make_pipeline(MyTrs(), MyEst())
    #print("=======")
    #print(clf.get_props_request())
    #print("=======")
    #clf.fit(X, y, sample_weight=sw, brand=brand)

    trs = MyTrs().set_props_request(
        None).set_props_request(
        {'fit': {'my_sw': 'new_param'}})
    clf = make_pipeline(trs, MyEst())
    clf.fit(X, y, sample_weight=sw, brand=brand, my_sw=my_data)

    with pytest.raises(ValueError, match="Requested properties are"):
        clf.fit(X, y, brand=brand)

    scorer = make_scorer(my_scorer, request_props="new_param")

    param_grid = {'myest__C': [0.1, 1]}

    print("@" * 150 + " GS")
    gs = GridSearchCV(clf, param_grid=param_grid, scoring=scorer)
    print("GS props request: ", gs.get_props_request())
    gs.fit(X, y, new_param=brand, sample_weight=sw, my_sw=sw, brand=brand)
