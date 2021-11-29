import pytest

from sklearn.datasets import load_iris
from sklearn.exceptions import NotFittedError
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor

from sklearn.metrics._plot.base import _check_estimator_and_target_is_binary

X, y = load_iris(return_X_y=True)
X_binary, y_binary = X[:100], y[:100]


@pytest.mark.parametrize(
    "estimator, target, target_type, err_type, err_msg",
    [
        (
            DecisionTreeClassifier(),
            y_binary,
            None,
            NotFittedError,
            "This DecisionTreeClassifier instance is not fitted yet",
        ),
        (
            DecisionTreeClassifier(),
            y_binary,
            "binary",
            NotFittedError,
            "This DecisionTreeClassifier instance is not fitted yet",
        ),
        (
            DecisionTreeRegressor().fit(X_binary, y_binary),
            y_binary,
            None,
            ValueError,
            "This plotting functionalities only support a binary classifier",
        ),
        (
            DecisionTreeRegressor().fit(X_binary, y_binary),
            y_binary,
            "binary",
            ValueError,
            "This plotting functionalities only support a binary classifier",
        ),
        (
            DecisionTreeClassifier().fit(X, y),
            y,
            None,
            ValueError,
            "This DecisionTreeClassifier instance is not a binary classifier",
        ),
        (
            DecisionTreeClassifier().fit(X, y),
            y,
            "multiclass",
            ValueError,
            "This DecisionTreeClassifier instance is not a binary classifier",
        ),
        (
            DecisionTreeClassifier().fit(X, y),
            y_binary,
            "multiclass",
            ValueError,
            "This DecisionTreeClassifier instance is not a binary classifier",
        ),
        (
            DecisionTreeClassifier().fit(X_binary, y_binary),
            y,
            None,
            ValueError,
            "The target y is not binary",
        ),
        (
            DecisionTreeClassifier().fit(X_binary, y_binary),
            y,
            "multiclass",
            ValueError,
            "The target y is not binary",
        ),
    ],
)
def test_check_estimator_and_target_is_binary(
    estimator, target, target_type, err_type, err_msg
):
    """Check that we raise the expected error when checking the estimator and target."""
    with pytest.raises(err_type, match=err_msg):
        _check_estimator_and_target_is_binary(estimator, target, target_type)
