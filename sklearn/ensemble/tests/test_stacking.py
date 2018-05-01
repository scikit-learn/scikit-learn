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
    stacked_estimators = [('lr', LogisticRegression()),
                          ('svc', LinearSVC())]
    clf = StackingClassifier(estimators=stacked_estimators,
                             meta_estimator=meta_estimator,
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
    stacked_estimators = [('lr', LinearRegression()),
                          ('svr', LinearSVR())]
    clf = StackingRegressor(estimators=stacked_estimators,
                            meta_estimator=meta_estimator,
                            cv=cv, random_state=42)
    clf.fit(X_train, y_train)
    clf.predict(X_test)
    assert clf.score(X_test, y_test) < 0.5
