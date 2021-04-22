import pytest

from sklearn.datasets import make_classification
from sklearn.datasets import make_regression
from sklearn.ensemble import GAMBoostingClassifier
from sklearn.ensemble import GAMBoostingRegressor
from sklearn.model_selection import train_test_split


def test_smoke_test_regression():
    """Smoke test for GAM regression"""
    X, y = make_regression(random_state=0, n_samples=1_000, n_features=6)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

    gam = GAMBoostingRegressor(random_state=0, max_depth=2)
    gam.fit(X_train, y_train)
    assert gam.score(X_test, y_test) >= 0.8


@pytest.mark.parametrize("n_classes", [2, 3])
def test_smoke_test_classification(n_classes):
    """Smoke test for GAM regression"""
    X, y = make_classification(random_state=0, n_samples=1_000, n_features=6,
                               n_informative=3, n_classes=n_classes)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0,
                                                        stratify=y)

    gam = GAMBoostingClassifier(random_state=0, max_depth=2)
    gam.fit(X_train, y_train)
    assert gam.score(X_test, y_test) >= 0.8
