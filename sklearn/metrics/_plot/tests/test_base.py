import numpy as np
import pytest

from sklearn.datasets import load_iris
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor

from sklearn.metrics._plot.base import _get_response_binary

X, y = load_iris(return_X_y=True)
X_binary, y_binary = X[:100], y[:100]


@pytest.mark.parametrize(
    "estimator, X, y, err_msg, params",
    [
        (
            DecisionTreeRegressor(),
            X_binary,
            y_binary,
            "Expected 'estimator' to be a binary classifier",
            {"response_method": "auto"},
        ),
        (
            DecisionTreeClassifier(),
            X_binary,
            y_binary,
            r"pos_label=unknown is not a valid label: It should be one of \[0 1\]",
            {"response_method": "auto", "pos_label": "unknown"},
        ),
        (
            DecisionTreeClassifier(),
            X,
            y,
            "Expected 'estimator' to be a binary classifier, but got"
            " DecisionTreeClassifier",
            {"response_method": "predict_proba"},
        ),
    ],
)
def test_get_response_error(estimator, X, y, err_msg, params):
    """Check that we raise the proper error messages in `_get_response_binary`."""

    estimator.fit(X, y)
    with pytest.raises(ValueError, match=err_msg):
        _get_response_binary(estimator, X, **params)


def test_get_response_predict_proba():
    """Check the behaviour of `_get_response_binary` using `predict_proba`."""
    classifier = DecisionTreeClassifier().fit(X_binary, y_binary)
    y_proba, pos_label = _get_response_binary(
        classifier, X_binary, response_method="predict_proba"
    )
    np.testing.assert_allclose(y_proba, classifier.predict_proba(X_binary)[:, 1])
    assert pos_label == 1

    y_proba, pos_label = _get_response_binary(
        classifier, X_binary, response_method="predict_proba", pos_label=0
    )
    np.testing.assert_allclose(y_proba, classifier.predict_proba(X_binary)[:, 0])
    assert pos_label == 0


def test_get_response_decision_function():
    """Check the behaviour of `_get_response_binary` using `decision_function`."""
    classifier = LogisticRegression().fit(X_binary, y_binary)
    y_score, pos_label = _get_response_binary(
        classifier, X_binary, response_method="decision_function"
    )
    np.testing.assert_allclose(y_score, classifier.decision_function(X_binary))
    assert pos_label == 1

    y_score, pos_label = _get_response_binary(
        classifier, X_binary, response_method="decision_function", pos_label=0
    )
    np.testing.assert_allclose(y_score, classifier.decision_function(X_binary) * -1)
    assert pos_label == 0
