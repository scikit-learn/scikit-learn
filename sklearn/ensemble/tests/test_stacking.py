"""Test the stacking classifier and regressor."""

# Authors: Guillaume Lemaitre <g.lemaitre58@gmail.com>
# License: BSD 3 clause

import pytest
import numpy as np

from sklearn.base import BaseEstimator
from sklearn.base import ClassifierMixin
from sklearn.base import RegressorMixin
from sklearn.base import clone

from sklearn.datasets import load_iris
from sklearn.datasets import load_diabetes

from sklearn.dummy import DummyClassifier
from sklearn.dummy import DummyRegressor
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LinearRegression
from sklearn.svm import LinearSVC
from sklearn.svm import LinearSVR
from sklearn.tree import DecisionTreeClassifier
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import scale

from sklearn.ensemble import StackingClassifier
from sklearn.ensemble import StackingRegressor

from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import KFold

from sklearn.utils.testing import assert_allclose
from sklearn.utils.estimator_checks import check_estimator
from sklearn.utils.estimator_checks import check_no_attributes_set_in_init

X_diabetes, y_diabetes = load_diabetes(return_X_y=True)
X_iris, y_iris = load_iris(return_X_y=True)


@pytest.mark.parametrize(
    "cv", [3, StratifiedKFold(n_splits=3, shuffle=True, random_state=42)]
)
@pytest.mark.parametrize(
    "final_estimator", [None, RandomForestClassifier(random_state=42)]
)
def test_stacking_classifier_iris(cv, final_estimator):

    # prescale the data to avoid convergence warning without using a pipeline
    # for later assert
    X_train, X_test, y_train, y_test = train_test_split(
        scale(X_iris), y_iris, stratify=y_iris, random_state=42
    )
    estimators = [('lr', LogisticRegression()),
                  ('svc', LinearSVC())]
    clf = StackingClassifier(
        estimators=estimators, final_estimator=final_estimator, cv=cv
    )
    clf.fit(X_train, y_train)
    clf.predict(X_test)
    clf.predict_proba(X_test)
    assert clf.score(X_test, y_test) > 0.8

    X_trans = clf.transform(X_test)
    assert X_trans.shape[1] == 6

    clf.set_params(lr=None)
    clf.fit(X_train, y_train)
    clf.predict(X_test)
    clf.predict_proba(X_test)

    X_trans = clf.transform(X_test)
    assert X_trans.shape[1] == 3


@pytest.mark.parametrize(
    "estimators",
    [[('lr', None), ('svc', LinearSVC(random_state=0))],
     [('lr', 'drop'), ('svc', LinearSVC(random_state=0))]]
)
def test_stacking_classifier_drop_estimator(estimators):
    # prescale the data to avoid convergence warning without using a pipeline
    # for later assert
    X_train, X_test, y_train, _ = train_test_split(
        scale(X_iris), y_iris, stratify=y_iris, random_state=42
    )
    rf = RandomForestClassifier(n_estimators=10, random_state=42)
    clf = StackingClassifier(
        estimators=[('svc', LinearSVC(random_state=0))],
        final_estimator=rf, cv=5
    )
    clf_drop = StackingClassifier(
        estimators=estimators, final_estimator=rf, cv=5
    )

    clf.fit(X_train, y_train)
    clf_drop.fit(X_train, y_train)
    assert_allclose(clf.predict(X_test), clf_drop.predict(X_test))
    assert_allclose(clf.predict_proba(X_test), clf_drop.predict_proba(X_test))
    assert_allclose(clf.transform(X_test), clf_drop.transform(X_test))


@pytest.mark.parametrize(
    "estimators",
    [[('lr', None), ('svr', LinearSVR(random_state=0))],
     [('lr', 'drop'), ('svr', LinearSVR(random_state=0))]]
)
def test_stacking_regressor_drop_estimator(estimators):
    # prescale the data to avoid convergence warning without using a pipeline
    # for later assert
    X_train, X_test, y_train, _ = train_test_split(
        scale(X_diabetes), y_diabetes, random_state=42
    )
    rf = RandomForestRegressor(n_estimators=10, random_state=42)
    reg = StackingRegressor(
        estimators=[('svr', LinearSVR(random_state=0))],
        final_estimator=rf, cv=5
    )
    reg_drop = StackingRegressor(
        estimators=estimators, final_estimator=rf, cv=5
    )

    reg.fit(X_train, y_train)
    reg_drop.fit(X_train, y_train)
    assert_allclose(reg.predict(X_test), reg_drop.predict(X_test))
    assert_allclose(reg.transform(X_test), reg_drop.transform(X_test))


@pytest.mark.parametrize(
    "cv", [3, KFold(n_splits=3, shuffle=True, random_state=42)]
)
@pytest.mark.parametrize(
    "final_estimator, predict_params",
    [(None, {}),
     (RandomForestRegressor(random_state=42), {}),
     (DummyRegressor(), {'return_std': True})]
)
def test_stacking_regressor_diabetes(cv, final_estimator, predict_params):
    # prescale the data to avoid convergence warning without using a pipeline
    # for later assert
    X_train, X_test, y_train, _ = train_test_split(
        scale(X_diabetes), y_diabetes, random_state=42
    )
    estimators = [('lr', LinearRegression()), ('svr', LinearSVR())]
    reg = StackingRegressor(
        estimators=estimators, final_estimator=final_estimator, cv=cv
    )
    reg.fit(X_train, y_train)
    result = reg.predict(X_test, **predict_params)
    expected_result_length = 2 if predict_params else 1
    if predict_params:
        assert len(result) == expected_result_length

    X_trans = reg.transform(X_test)
    assert X_trans.shape[1] == 2

    reg.set_params(lr=None)
    reg.fit(X_train, y_train)
    reg.predict(X_test)

    X_trans = reg.transform(X_test)
    assert X_trans.shape[1] == 1


def test_stacking_classifier_drop_binary_prob():
    # check that classifier will drop one of the probability column for
    # binary classification problem

    # Select only the 2 first classes
    X_, y_ = scale(X_iris[:100]), y_iris[:100]

    estimators = [
        ('lr', LogisticRegression()), ('rf', RandomForestClassifier())
    ]
    clf = StackingClassifier(estimators=estimators)
    clf.fit(X_, y_)
    X_meta = clf.transform(X_)
    assert X_meta.shape[1] == 2


class NoWeightRegressor(BaseEstimator, RegressorMixin):
    def fit(self, X, y):
        self.reg = DummyRegressor()
        return self.reg.fit(X, y)

    def predict(self, X):
        return np.ones(X.shape[0])


class NoWeightClassifier(BaseEstimator, ClassifierMixin):
    def fit(self, X, y):
        self.clf = DummyClassifier()
        return self.clf.fit(X, y)


@pytest.mark.parametrize(
    "y, params, type_err, msg_err",
    [(y_iris,
      {'estimators': None},
      ValueError, "Invalid 'estimators' attribute,"),
     (y_iris,
      {'estimators': []},
      ValueError, "Invalid 'estimators' attribute,"),
     (y_iris,
      {'estimators': [('lr', LinearRegression()),
                      ('svm', LinearSVC(max_iter=5e4))],
       'stack_method': 'predict_proba'},
      ValueError, 'does not implement the method'),
     (y_iris,
      {'estimators': [('lr', LogisticRegression()),
                      ('cor', NoWeightClassifier())]},
      TypeError, 'does not support sample weight'),
     (y_iris,
      {'estimators': [('lr', None), ('svm', None)]},
      ValueError, 'All estimators are None'),
     (y_iris,
      {'estimators': [('lr', 'drop'), ('svm', None)]},
      ValueError, 'All estimators are None'),
     (y_iris,
      {'estimators': [('lr', 'drop'), ('svm', 'drop')]},
      ValueError, 'All estimators are None'),
     (y_iris,
      {'estimators': [('lr', LogisticRegression()), ('svm', LinearSVC())],
       'final_estimator': RandomForestRegressor()},
      ValueError, 'parameter should be a classifier.')]
)
def test_stacking_classifier_error(y, params, type_err, msg_err):
    with pytest.raises(type_err, match=msg_err):
        clf = StackingClassifier(**params, cv=3)
        clf.fit(
            scale(X_iris), y, sample_weight=np.ones(X_iris.shape[0])
        )


@pytest.mark.parametrize(
    "y, params, type_err, msg_err",
    [(y_diabetes,
      {'estimators': None},
      ValueError, "Invalid 'estimators' attribute,"),
     (y_diabetes,
      {'estimators': []},
      ValueError, "Invalid 'estimators' attribute,"),
     (y_diabetes,
      {'estimators': [('lr', LogisticRegression()), ('svm', LinearSVR())],
       'stack_method': 'predict_proba'},
      ValueError, 'does not implement the method'),
     (y_diabetes,
      {'estimators': [('lr', LinearRegression()),
                      ('cor', NoWeightRegressor())]},
      TypeError, 'does not support sample weight'),
     (y_diabetes,
      {'estimators': [('lr', None), ('svm', None)]},
      ValueError, 'All estimators are None'),
     (y_diabetes,
      {'estimators': [('lr', 'drop'), ('svm', None)]},
      ValueError, 'All estimators are None'),
     (y_diabetes,
      {'estimators': [('lr', 'drop'), ('svm', 'drop')]},
      ValueError, 'All estimators are None'),
     (y_diabetes,
      {'estimators': [('lr', LinearRegression()), ('svm', LinearSVR())],
       'final_estimator': RandomForestClassifier()},
      ValueError, 'parameter should be a regressor.')]
)
def test_stacking_regressor_error(y, params, type_err, msg_err):
    with pytest.raises(type_err, match=msg_err):
        reg = StackingRegressor(**params, cv=3)
        reg.fit(
            scale(X_diabetes), y, sample_weight=np.ones(X_diabetes.shape[0])
        )


@pytest.mark.parametrize(
    "stacking_estimator",
    [StackingClassifier(estimators=[('lr', LogisticRegression()),
                                    ('svm', LinearSVC())]),
     StackingRegressor(estimators=[('lr', LinearRegression()),
                                   ('svm', LinearSVR(max_iter=1e4))])]
)
def test_stacking_named_estimators(stacking_estimator):
    stacking_estimator.fit(scale(X_iris), y_iris)
    estimators = stacking_estimator.named_estimators_
    assert len(estimators) == 2
    assert sorted(list(estimators.keys())) == sorted(['lr', 'svm'])


@pytest.mark.parametrize(
    "stacking_estimator",
    [StackingClassifier(estimators=[('lr', LogisticRegression()),
                                    ('rf', RandomForestClassifier()),
                                    ('svm', LinearSVC())]),
     StackingRegressor(estimators=[('lr', LinearRegression()),
                                   ('rf', RandomForestRegressor()),
                                   ('svm', LinearSVR(max_iter=1e4))])]
)
def test_stacking_named_estimators_dropped(stacking_estimator):
    stacking_estimator.set_params(rf='drop')
    stacking_estimator.fit(scale(X_iris), y_iris)
    estimators = stacking_estimator.named_estimators_
    assert 'rf' not in estimators.keys()
    assert len(estimators) == 2
    assert sorted(list(estimators.keys())) == sorted(['lr', 'svm'])


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


@pytest.mark.parametrize(
    "estimator, X, y",
    [(StackingClassifier(
        estimators=[('lr', LogisticRegression(random_state=0)),
                    ('svm', LinearSVC(random_state=0))]),
      X_iris[:100], y_iris[:100]),  # keep only classes 0 and 1
     (StackingRegressor(
         estimators=[('lr', LinearRegression()),
                     ('svm', LinearSVR(random_state=0))]),
      X_diabetes, y_diabetes)],
    ids=['StackingClassifier', 'StackingRegressor']
)
def test_stacking_randomness(estimator, X, y):
    estimator_full = clone(estimator)
    estimator_full.set_params(
        cv=KFold(shuffle=True, random_state=np.random.RandomState(0))
    )

    estimator_drop = clone(estimator)
    estimator_drop.set_params(lr='drop')
    estimator_drop.set_params(
        cv=KFold(shuffle=True, random_state=np.random.RandomState(0))
    )

    assert_allclose(
        estimator_full.fit(X, y).transform(X)[:, 1:],
        estimator_drop.fit(X, y).transform(X)
    )


@pytest.mark.parametrize(
    "estimator",
    [StackingClassifier(
        estimators=[('lr', LogisticRegression(random_state=0)),
                    ('tree', DecisionTreeClassifier(random_state=0))]),
     StackingRegressor(
         estimators=[('lr', LinearRegression()),
                     ('tree', DecisionTreeRegressor(random_state=0))])],
    ids=['StackingClassifier', 'StackingRegressor']
)
def test_check_estimators_stacking_estimator(estimator):
    # FIXME: to be removed when meta-estimators can specified themselves
    # their testing parameters (for required parameters).
    check_estimator(estimator)
    check_no_attributes_set_in_init(estimator.__class__.__name__, estimator)
