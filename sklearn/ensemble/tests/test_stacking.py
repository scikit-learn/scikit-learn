import pytest

from sklearn.datasets import load_iris

from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import StackingClassifier

from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold


@pytest.mark.parametrize(
    "cv",
    [3, StratifiedKFold(shuffle=True, random_state=42)]
)
def test_iris(cv):
    X, y = load_iris(return_X_y=True)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)
    stacked_estimators = [('lr', LogisticRegression()),
                          ('svc', LinearSVC())]
    clf = StackingClassifier(estimators=stacked_estimators,
                             classifier=RandomForestClassifier(),
                             cv=cv, random_state=42)
    clf.fit(X_train, y_train)
    clf.predict(X_test)
    clf.predict_proba(X_test)
    assert clf.score(X_test, y_test) > 0.8
