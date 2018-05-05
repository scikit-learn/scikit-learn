"""Test the stacking classifier and regressor."""

# Authors: Guillaume Lemaitre <g.lemaitre58@gmail.com>
# License: BSD 3 clause

import pytest
import numpy as np

from sklearn.base import BaseEstimator
from sklearn.base import ClassifierMixin
from sklearn.base import RegressorMixin

from sklearn.datasets import load_iris
from sklearn.datasets import load_diabetes

from sklearn.dummy import DummyClassifier
from sklearn.dummy import DummyRegressor
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LinearRegression
from sklearn.svm import LinearSVC
from sklearn.svm import LinearSVR
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor

from sklearn.ensemble import StackingClassifier
from sklearn.ensemble import StackingRegressor

from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import LeaveOneOut

X_diabetes, y_diabetes = load_diabetes(return_X_y=True)
X_iris, y_iris = load_iris(return_X_y=True)


@pytest.mark.parametrize(
    "cv",
    [3, StratifiedKFold(shuffle=True, random_state=42), LeaveOneOut()]
)
@pytest.mark.parametrize(
    "final_estimator",
    [None, RandomForestClassifier(random_state=42)]
)
def test_stacking_classifier_iris(cv, final_estimator):
    X, y = load_iris(return_X_y=True)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)
    estimators = [('lr', LogisticRegression()), ('svc', LinearSVC())]
    clf = StackingClassifier(estimators=estimators,
                             final_estimator=final_estimator,
                             cv=cv, random_state=42)
    clf.fit(X_train, y_train)
    clf.predict(X_test)
    clf.predict_proba(X_test)
    assert clf.score(X_test, y_test) > 0.8

    clf.set_params(lr=None)
    clf.fit(X_train, y_train)
    clf.predict(X_test)
    clf.predict_proba(X_test)


@pytest.mark.parametrize(
    "cv",
    [3, StratifiedKFold(shuffle=True, random_state=42), LeaveOneOut()]
)
@pytest.mark.parametrize(
    "final_estimator",
    [None, RandomForestRegressor(random_state=42)]
)
def test_stacking_regressor_diabetes(cv, final_estimator):
    X, y = load_diabetes(return_X_y=True)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)
    estimators = [('lr', LinearRegression()), ('svr', LinearSVR())]
    reg = StackingRegressor(estimators=estimators,
                            final_estimator=final_estimator,
                            cv=cv, random_state=42)
    reg.fit(X_train, y_train)
    reg.predict(X_test)
    assert reg.score(X_test, y_test) < 0.5

    reg.set_params(lr=None)
    reg.fit(X_train, y_train)
    reg.predict(X_test)


class NoWeightRegressor(BaseEstimator, RegressorMixin):
    def __init__(self):
        self.reg = DummyRegressor()

    def fit(self, X, y):
        return self.reg.fit(X, y)


class NoWeightClassifier(BaseEstimator, ClassifierMixin):
    def __init__(self):
        self.clf = DummyClassifier()

    def fit(self, X, y):
        return self.clf.fit(X, y)


@pytest.mark.parametrize(
    "X, y, estimators, methods, final_estimator, type_err, msg_err",
    [(X_iris, y_iris, None, 'auto', RandomForestClassifier(),
      AttributeError, 'Invalid `estimators`'),
     (X_iris, y_iris, [('lr', LogisticRegression()), ('svm', LinearSVC())],
      'random', RandomForestClassifier(),
      AttributeError, 'When "method_estimators" is a string'),
     (X_iris, y_iris, [('lr', LogisticRegression()), ('svm', LinearSVC())],
      ['predict'], RandomForestClassifier(),
      AttributeError, 'When "method_estimators" is a list'),
     (X_iris, y_iris, [('lr', LinearRegression()),
                       ('svm', LinearSVR())],
      ['predict', 'predict_proba'], None,
      ValueError, 'does not implement the method'),
     (X_iris, y_iris, [('lr', LogisticRegression()),
                       ('cor', NoWeightClassifier())],
      'auto', None, ValueError, 'does not support sample weight'),
     (X_iris, y_iris, [('lr', None), ('svm', None)],
      'auto', None, ValueError, 'All estimators are None'),
     (X_iris, y_iris, [('lr', LogisticRegression()), ('svm', LinearSVC())],
      'auto', RandomForestRegressor(),
      AttributeError, 'attribute should be a classifier.')]
)
def test_stacking_classifier_error(X, y, estimators, methods, final_estimator,
                                   type_err, msg_err):
    with pytest.raises(type_err, match=msg_err):
        clf = StackingClassifier(estimators=estimators,
                                 method_estimators=methods,
                                 final_estimator=final_estimator)
        clf.fit(X, y, sample_weight=np.ones(X.shape[0]))


@pytest.mark.parametrize(
    "X, y, estimators, methods, final_estimator, type_err, msg_err",
    [(X_diabetes, y_diabetes, None, 'auto', RandomForestRegressor(),
      AttributeError, 'Invalid `estimators`'),
     (X_diabetes, y_diabetes, [('lr', LinearRegression()),
                               ('svm', LinearSVR())],
      'random', RandomForestRegressor(),
      AttributeError, 'When "method_estimators" is a string'),
     (X_diabetes, y_diabetes, [('lr', LinearRegression()),
                               ('svm', LinearSVR())],
      ['predict'], RandomForestRegressor(),
      AttributeError, 'When "method_estimators" is a list'),
     (X_diabetes, y_diabetes, [('lr', LinearRegression()),
                               ('svm', LinearSVR())],
      ['predict', 'predict_proba'], None,
      ValueError, 'does not implement the method'),
     (X_diabetes, y_diabetes, [('lr', LinearRegression()),
                               ('cor', NoWeightRegressor())],
      'auto', None, ValueError, 'does not support sample weight'),
     (X_diabetes, y_diabetes, [('lr', None), ('svm', None)],
      'auto', None, ValueError, 'All estimators are None'),
     (X_diabetes, y_diabetes, [('lr', LinearRegression()),
                               ('svm', LinearSVR())],
      'auto', RandomForestClassifier(),
      AttributeError, 'attribute should be a regressor.')]
)
def test_stacking_regressor_error(X, y, estimators, methods, final_estimator,
                                  type_err, msg_err):
    with pytest.raises(type_err, match=msg_err):
        reg = StackingRegressor(estimators=estimators,
                                method_estimators=methods,
                                final_estimator=final_estimator)
        reg.fit(X, y, sample_weight=np.ones(X.shape[0]))


@pytest.mark.parametrize(
    "stacking_estimator",
    [StackingClassifier(estimators=[('lr', LogisticRegression()),
                                    ('svm', LinearSVC())]),
     StackingRegressor(estimators=[('lr', LinearRegression()),
                                   ('svm', LinearSVR())])]
)
def test_stacking_named_estimators(stacking_estimator):
    estimators = stacking_estimator.named_estimators
    assert len(estimators) == 2
    assert list(estimators.keys()) == ['lr', 'svm']


@pytest.mark.parametrize(
    "stacking_estimator",
    [StackingClassifier(estimators=[('lr', LogisticRegression()),
                                    ('svm', LinearSVC())]),
     StackingRegressor(estimators=[('lr', LinearRegression()),
                                   ('svm', LinearSVR())])]
)
def test_stacking_set_get_params(stacking_estimator):
    params = stacking_estimator.get_params()
    assert 'lr' in list(params.keys())
    assert 'svm' in list(params.keys())

    stacking_estimator.set_params(lr=None)
    params = stacking_estimator.get_params()
    assert params['lr'] is None
