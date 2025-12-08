import numpy as np

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.compose import ColumnTransformer
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import Normalizer, StandardScaler
from sklearn.utils._repr_html.estimator import estimator_html_repr
from sklearn.utils._repr_html.features import _features_html
from sklearn.utils._testing import MinimalTransformer

ct = ColumnTransformer([("norm1", Normalizer(), [0, 1])])
rng = np.random.RandomState(42)


def test_n_features_not_fitted():
    assert "2 features" not in estimator_html_repr(ct)
    assert "x0" not in estimator_html_repr(ct)
    assert "x1" not in estimator_html_repr(ct)


def test_n_features_fitted():
    X = np.array([[0, 2], [1, 1]])
    ct.fit(X)
    assert "2 features" in estimator_html_repr(ct)
    assert "x0" in estimator_html_repr(ct)
    assert "x1" in estimator_html_repr(ct)


def test_n_features_with_Transformer():
    class Trans(TransformerMixin, BaseEstimator):
        def fit(self, X, y=None):
            return self

        def transform(self, X, y=None):
            # 1D Series -> 2D DataFrame
            if hasattr(X, "to_frame"):
                return X.to_frame()
            # 1D array -> 2D array
            if getattr(X, "ndim", 2) == 1:
                return np.atleast_2d(X).T
            return X

    X = np.array([[0, 2], [1, 1]])
    ct = ColumnTransformer([("trans", Trans(), [0, 1])])
    ct.fit(X)
    estimator_html_repr(ct)


def test_n_features_with_MinimalTransformer():
    X, y = np.array([[0, 1], [1, 1]]), np.array([[0, 1]])
    ct = ColumnTransformer([("minimal", MinimalTransformer(), [0, 1])])
    model = Pipeline([("transformer", ct)])
    model.fit(X, y)

    estimator_html_repr(model)


def test_features_html_no_fitted_class():
    """Test that features HTML is generated without fitted CSS class."""
    features = ["feature1", "feature2", "feature3"]
    html = _features_html(features)

    assert "3 features" in html
    assert "feature1" in html
    assert "feature2" in html
    assert "feature3" in html
    assert 'class="features "' in html or 'class="features"' in html


def test_features_html_with_fitted_class():
    """Test that features HTML includes fitted CSS class when specified."""
    features = ["feature1", "feature2"]
    html = _features_html(features, is_fitted_css_class="fitted")

    assert "2 features" in html
    assert "fitted" in html
    assert 'class="features fitted"' in html


def test_features_html_correct_feature_names():
    """Test that feature names are correctly escaped and displayed in table."""
    features = ["age", "income", "<script>alert('xss')</script>"]
    html = _features_html(features)

    assert "age" in html
    assert "income" in html
    # Check that HTML is escaped
    assert "&lt;script&gt;" in html
    assert "<script>alert" not in html
    assert "3 features" in html


def test_features_html_column_transformer_feature_count():
    """Test feature count with ColumnTransformer in Pipeline."""

    # Create a ColumnTransformer with multiple transformers
    ct = ColumnTransformer(
        [
            ("scaler", StandardScaler(), [0, 1]),
            ("pca", PCA(n_components=2, random_state=42), [2, 3, 4]),
        ]
    )

    # Fit the transformer
    X = rng.randn(10, 5)
    ct.fit(X)

    # Get feature names
    feature_names = ct.get_feature_names_out()
    html = _features_html(feature_names)

    # Should have 4 features (2 from scaler + 2 from PCA)
    assert f"{len(feature_names)} features" in html
    assert len(feature_names) == 4


def test_features_html_with_minimal_transformer():
    """Test that _features_html works with MinimalTransformer (non-full API)."""

    # Create a pipeline with MinimalTransformer
    pipe = Pipeline([("minimal", MinimalTransformer()), ("scaler", StandardScaler())])

    # Fit the pipeline
    X = rng.randn(10, 3)
    pipe.fit(X)

    # MinimalTransformer passes through features, StandardScaler preserves them
    # This should not crash
    try:
        feature_names = pipe.get_feature_names_out()
        html = _features_html(feature_names)
        assert f"{len(feature_names)} features" in html
    except AttributeError:
        # If get_feature_names_out is not available, test with manual features
        features = [f"x{i}" for i in range(3)]
        html = _features_html(features)
        assert "3 features" in html


def test_features_html_empty_features():
    """Test that _features_html handles empty feature list."""
    features = []
    html = _features_html(features)

    assert "0 features" in html
    assert "<tbody>" in html


def test_features_html_special_characters():
    """Test that special characters in feature names are properly escaped."""
    features = ["feature&1", 'feature"2', "feature'3", "feature>4", "feature<5"]
    html = _features_html(features)

    assert "&amp;" in html
    assert "&quot;" in html or "&#x27;" in html
    assert "&gt;" in html
    assert "&lt;" in html
    assert "5 features" in html


def test_features_html_structure():
    """Test that HTML structure contains expected elements."""
    features = ["feat1", "feat2"]
    html = _features_html(features)

    # Check for key structural elements
    assert "<details>" in html
    assert "<summary>" in html
    assert "</summary>" in html
    assert "</details>" in html
    assert '<table class="features-table">' in html
    assert "<tbody>" in html
    assert "</tbody>" in html
    assert '<i class="copy-paste-icon"' in html
    assert "copyFeatureNamesToClipboard" in html
