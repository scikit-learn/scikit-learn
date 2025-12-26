import numpy as np
import pytest

from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import Normalizer, StandardScaler
from sklearn.utils._repr_html.estimator import estimator_html_repr
from sklearn.utils._repr_html.features import _features_html
from sklearn.utils._testing import MinimalTransformer

ct = ColumnTransformer([("norm1", Normalizer(), [0, 1])])
rng = np.random.RandomState(42)


def test_n_features_not_fitted():
    out = estimator_html_repr(ct)
    assert "2 features" not in out
    assert "x0" not in out
    assert "x1" not in out
    assert "<div class='features fitted'>" not in out


def test_with_MinimalTransformer():
    X, y = np.array([[0, 1], [1, 1]]), np.array([[0, 1]])
    ct = ColumnTransformer(
        [("minimal", MinimalTransformer(), [0])], remainder="passthrough"
    )
    model = Pipeline([("transformer", ct)])
    model.fit(X, y)

    out = estimator_html_repr(model)
    assert "1 features" in out
    assert "x0" in out
    assert "passthrough" in out


@pytest.mark.parametrize(
    "pandas, feature_cols",
    [
        (True, ["A", "B"]),
        (False, ["x0", "x1"]),
    ],
)
def test_estimator_html_repr_col_names(pandas, feature_cols):
    """Test features names are kept with pandas col names and generic."""
    if pandas:
        pd = pytest.importorskip("pandas")
        X = pd.DataFrame({"A": [0, 2], "B": [1, 1]})
    else:
        X = np.array([[0, 2], [1, 1]])

    ct.fit(X)
    out = estimator_html_repr(ct)
    assert feature_cols[0] in out
    assert feature_cols[1] in out


def test_features_html_with_fitted_class():
    """Test that features HTML includes fitted CSS class when specified."""
    features = ["feature1", "feature2"]
    html = _features_html(features, is_fitted_css_class="fitted")

    assert "2 features" in html
    assert 'class="features fitted"' in html


def test_features_html_correct_feature_names():
    """Test that feature names are correctly escaped and displayed in table."""
    features = ["age", "income", "<script>alert('xss')</script>"]
    html = _features_html(features)

    assert "age" in html
    assert "income" in html

    assert "&lt;script&gt;" in html
    assert "<script>alert" not in html
    assert "3 features" in html


def test_features_html_with_pipeline():
    """Test works with MinimalTransformer in a pipeline."""

    pipe = Pipeline([("minimal", MinimalTransformer()), ("scaler", StandardScaler())])

    X = rng.randn(10, 3)
    pipe.fit(X)
    html = estimator_html_repr(pipe)
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

    assert "<details>" in html
    assert "<summary>" in html
    assert "</summary>" in html
    assert "</details>" in html
    assert '<table class="features-table">' in html
    assert "<tbody>" in html
    assert "</tbody>" in html
    assert '<i class="copy-paste-icon"' in html
    assert "copyFeatureNamesToClipboard" in html
