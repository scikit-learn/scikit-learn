"""
Test reproducing issue #32394 (AttributeError in is_regressor)
and validating the fix for missing __sklearn_tags__.
"""

import pytest
from sklearn.base import is_regressor
from sklearn.utils._tags import get_tags


def test_get_tags_handles_missing_dunder_sklearn_tags():
    """Ensure get_tags works even when estimator lacks __sklearn_tags__."""

    class DummyEstimator:
        def fit(self, X, y=None): return self
        def predict(self, X): return X

    # Should not raise
    tags = get_tags(DummyEstimator())
    assert isinstance(tags, dict), "get_tags should return a dictionary"
    assert "non_deterministic" in tags, "Expected default tags to include 'non_deterministic'"


def test_is_regressor_does_not_raise_when_missing_tags():
    """Ensure is_regressor does not raise AttributeError for basic estimators."""

    class MyEstimator:
        def fit(self, X, y=None): return self
        def predict(self, X): return X

    # Should not raise any error
    try:
        result = is_regressor(MyEstimator())
        assert isinstance(result, bool)
    except AttributeError as e:
        pytest.fail(f"is_regressor raised unexpected AttributeError: {e}")
