"""Common tests for metaestimators"""

import functools

import numpy as np

from sklearn.base import BaseEstimator
from sklearn.externals.six import iterkeys
from sklearn.datasets import make_classification
from sklearn.utils.testing import assert_true, assert_false, assert_raises
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV
from sklearn.feature_selection import RFE, RFECV
from sklearn.ensemble import BaggingClassifier


class DelegatorData(object):
    def __init__(self, name, construct, skip_methods=(),
                 fit_args=make_classification()):
        self.name = name
        self.construct = construct
        self.fit_args = fit_args
        self.skip_methods = skip_methods


DELEGATING_METAESTIMATORS = [
    DelegatorData('Pipeline', lambda est: Pipeline([('est', est)])),
    DelegatorData('GridSearchCV',
                  lambda est: GridSearchCV(
                      est, param_grid={'param': [5]}, cv=2),
                  skip_methods=['score']),
    DelegatorData('RandomizedSearchCV',
                  lambda est: RandomizedSearchCV(
                      est, param_distributions={'param': [5]}, cv=2, n_iter=1),
                  skip_methods=['score']),
    DelegatorData('RFE', RFE,
                  skip_methods=['transform', 'inverse_transform', 'score']),
    DelegatorData('RFECV', RFECV,
                  skip_methods=['transform', 'inverse_transform', 'score']),
    DelegatorData('BaggingClassifier', BaggingClassifier,
                  skip_methods=['transform', 'inverse_transform', 'score',
                                'predict_proba', 'predict_log_proba', 'predict'])
]


def test_metaestimator_delegation():
    # Ensures specified metaestimators have methods iff subestimator does
    def hides(method):
        @property
        def wrapper(obj):
            if obj.hidden_method == method.__name__:
                raise AttributeError('%r is hidden' % obj.hidden_method)
            return functools.partial(method, obj)
        return wrapper

    class SubEstimator(BaseEstimator):
        def __init__(self, param=1, hidden_method=None):
            self.param = param
            self.hidden_method = hidden_method

        def fit(self, X, y=None, *args, **kwargs):
            self.coef_ = np.arange(X.shape[1])
            return True

        def _check_fit(self):
            if not hasattr(self, 'coef_'):
                raise RuntimeError('Estimator is not fit')

        @hides
        def inverse_transform(self, X, *args, **kwargs):
            self._check_fit()
            return X

        @hides
        def transform(self, X, *args, **kwargs):
            self._check_fit()
            return X

        @hides
        def predict(self, X, *args, **kwargs):
            self._check_fit()
            return np.ones(X.shape[0])

        @hides
        def predict_proba(self, X, *args, **kwargs):
            self._check_fit()
            return np.ones(X.shape[0])

        @hides
        def predict_log_proba(self, X, *args, **kwargs):
            self._check_fit()
            return np.ones(X.shape[0])

        @hides
        def decision_function(self, X, *args, **kwargs):
            self._check_fit()
            return np.ones(X.shape[0])

        @hides
        def score(self, X, *args, **kwargs):
            self._check_fit()
            return 1.0

    methods = [k for k in iterkeys(SubEstimator.__dict__)
               if not k.startswith('_') and not k.startswith('fit')]
    methods.sort()

    for delegator_data in DELEGATING_METAESTIMATORS:
        delegate = SubEstimator()
        delegator = delegator_data.construct(delegate)
        for method in methods:
            if method in delegator_data.skip_methods:
                continue
            assert_true(hasattr(delegate, method))
            assert_true(hasattr(delegator, method),
                        msg="%s does not have method %r when its delegate does"
                            % (delegator_data.name, method))
            # delegation before fit raises an exception
            assert_raises(Exception, getattr(delegator, method),
                          delegator_data.fit_args[0])

        delegator.fit(*delegator_data.fit_args)
        for method in methods:
            if method in delegator_data.skip_methods:
                continue
            # smoke test delegation
            getattr(delegator, method)(delegator_data.fit_args[0])

        for method in methods:
            if method in delegator_data.skip_methods:
                continue
            delegate = SubEstimator(hidden_method=method)
            delegator = delegator_data.construct(delegate)
            assert_false(hasattr(delegate, method))
            assert_false(hasattr(delegator, method),
                         msg="%s has method %r when its delegate does not"
                             % (delegator_data.name, method))
