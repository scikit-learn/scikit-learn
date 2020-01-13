import numpy as np
import pytest

from sklearn.feature_selection import SelectKBest
from sklearn.pipeline import make_pipeline
from sklearn.base import BaseEstimator, ClassifierMixin, TransformerMixin
from sklearn.utils.validation import _validate_required_props
from sklearn.datasets import make_classification


class MyEst(ClassifierMixin, BaseEstimator):
    def __init__(self):
        self._props_request = {'fit': ['sample_weight', 'brand']}

    def fit(self, X, y, **fit_params):
        _validate_required_props(self, fit_params, 'fit')
        assert fit_params.keys() == \
            self._get_props_request_mapping('fit').keys()
        return self


class MyTrs(TransformerMixin, BaseEstimator):
    def __init__(self):
        self._props_request = {'fit': ['sample_weight']}

    def fit(self, X, y=None, **fit_params):
        _validate_required_props(self, fit_params, 'fit')
        self._estimator = SelectKBest().fit(X, y)
        assert fit_params.keys() == \
            self._get_props_request_mapping('fit').keys()
        return self

    def transform(self, X, y=None):
        return self._estimator.transform(X)


def test_pipeline():
    X, y = make_classification()
    sw = np.random.rand(len(X))
    my_data = [5, 6]
    brand = ['my brand']

    clf = make_pipeline(MyTrs(), MyEst())
    clf.fit(X, y, sample_weight=sw, brand=brand)

    trs = MyTrs().set_props_request(
        None).set_props_request(
        {'fit': {'new_param': 'my_sw'}})
    clf = make_pipeline(trs, MyEst())
    clf.fit(X, y, sample_weight=sw, brand=brand, my_sw=my_data)

    with pytest.raises(ValueError, match="Requested properties are"):
        clf.fit(X, y, brand=brand)
