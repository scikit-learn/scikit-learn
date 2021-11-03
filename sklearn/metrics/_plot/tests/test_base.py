import pytest

from sklearn.datasets import load_iris
from sklearn.exceptions import NotFittedError
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor

from sklearn.metrics._plot.base import _check_estimator_target

X, y = load_iris(return_X_y=True)
X_binary, y_binary = X[:100], y[:100]


@pytest.mark.parametrize(
    "estimator, target, err_type, err_msg",
    [
        (
            DecisionTreeClassifier(),
            y_binary,
            NotFittedError,
            "This DecisionTreeClassifier instance is not fitted yet",
        ),
        (
            DecisionTreeRegressor().fit(X_binary, y_binary),
            y_binary,
            ValueError,
            "This plotting functionalities only support a binary classifier",
        ),
        (
            DecisionTreeClassifier().fit(X, y),
            y,
            ValueError,
            "This DecisionTreeClassifier instance is not a binary classifier",
        ),
        (
            DecisionTreeClassifier().fit(X_binary, y_binary),
            y,
            ValueError,
            "The target y is not binary",
        ),
    ],
)
def test_check_estimator_target(estimator, target, err_type, err_msg):
    """Check that we raise the expected error when checking the estimator and target."""
    with pytest.raises(err_type, match=err_msg):
        _check_estimator_target(estimator, target)
