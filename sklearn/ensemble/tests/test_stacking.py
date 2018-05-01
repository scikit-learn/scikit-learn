import pytest

from sklearn.datasets import load_iris
from sklearn.datasets import load_diabetes

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

X_diabetes, y_diabetes = load_diabetes(return_X_y=True)
X_iris, y_iris = load_iris(return_X_y=True)


@pytest.mark.parametrize(
    "cv",
    [3, StratifiedKFold(shuffle=True, random_state=42)]
)
@pytest.mark.parametrize(
    "meta_estimator",
    [None, RandomForestClassifier(random_state=42)]
)
def test_stacking_classifier_iris(cv, meta_estimator):
    X, y = load_iris(return_X_y=True)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)
    base_estimators = [('lr', LogisticRegression()),
                       ('svc', LinearSVC())]
    clf = StackingClassifier(base_estimators=base_estimators,
                             stacking_estimator=meta_estimator,
                             cv=cv, random_state=42)
    clf.fit(X_train, y_train)
    clf.predict(X_test)
    clf.predict_proba(X_test)
    assert clf.score(X_test, y_test) > 0.8


@pytest.mark.parametrize(
    "cv",
    [3, StratifiedKFold(shuffle=True, random_state=42)]
)
@pytest.mark.parametrize(
    "meta_estimator",
    [None, RandomForestRegressor(random_state=42)]
)
def test_stacking_regressor_diabetes(cv, meta_estimator):
    X, y = load_diabetes(return_X_y=True)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)
    base_estimators = [('lr', LinearRegression()),
                       ('svr', LinearSVR())]
    clf = StackingRegressor(base_estimators=base_estimators,
                            stacking_estimator=meta_estimator,
                            cv=cv, random_state=42)
    clf.fit(X_train, y_train)
    clf.predict(X_test)
    assert clf.score(X_test, y_test) < 0.5


@pytest.mark.parametrize(
    "X, y, base_estimators, methods, stacking_estimator, msg_err",
    [(X_iris, y_iris, None, 'auto', RandomForestClassifier(),
      'Invalid `base_estimators`'),
     (X_iris, y_iris, [('lr', LogisticRegression()), ('svm', LinearSVC())],
      'random', RandomForestClassifier(),
      'When "method_base_estimators" is a string'),
     (X_iris, y_iris, [('lr', LogisticRegression()), ('svm', LinearSVC())],
      ['predict'], RandomForestClassifier(),
      'When "method_base_estimators" is a list')]
)
def test_stacking_classifier_error(X, y, base_estimators, methods,
                                   stacking_estimator, msg_err):
    with pytest.raises(AttributeError, match=msg_err):
        clf = StackingClassifier(base_estimators=base_estimators,
                                 method_base_estimators=methods,
                                 stacking_estimator=stacking_estimator)
        clf.fit(X, y)
