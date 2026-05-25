import numpy as np
import pytest

from sklearn.base import BaseEstimator, TransformerMixin, clone
from sklearn.compose import ColumnTransformer, make_column_transformer
from sklearn.datasets import load_iris
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import TunedThresholdClassifierCV
from sklearn.pipeline import FeatureUnion, Pipeline, make_pipeline
from sklearn.preprocessing import (
    FunctionTransformer,
    MinMaxScaler,
    Normalizer,
    StandardScaler,
)
from sklearn.utils._repr_html.estimator import estimator_html_repr
from sklearn.utils._repr_html.features import _features_html
from sklearn.utils._testing import MinimalTransformer

ct = ColumnTransformer([("norm1", Normalizer(), [0, 1])], remainder="passthrough")
ct2 = FeatureUnion(
    [("pca", PCA(n_components=1)), ("svd", TruncatedSVD(n_components=2))]
)
rng = np.random.RandomState(42)


def test_n_features_not_fitted():
    out = estimator_html_repr(ct)
    assert "2 features" not in out
    assert "x0" not in out
    assert "x1" not in out
    assert "<div class='features fitted'>" not in out


def test_with_MinimalTransformer():
    """Test works with MinimalTransformer in a pipeline
    (doesn't inherit from BaseEstimator)"""
    X, y = np.array([[0, 1], [1, 1]]), np.array([[0, 1]])

    model = Pipeline([("transformer", MinimalTransformer())])
    model.fit(X, y)
    out = estimator_html_repr(model)
    assert "MinimalTransformer" in out


@pytest.mark.parametrize(
    "pandas, feature_cols",
    [
        (True, ["Feature A", "Feature B"]),
        (False, ["x0", "x1"]),
    ],
)
def test_estimator_html_repr_col_names(pandas, feature_cols):
    """Test features names are kept with pandas col names and generic."""
    if pandas:
        pd = pytest.importorskip("pandas")
        X = pd.DataFrame({"Feature A": [0, 2], "Feature B": [1, 1]})
    else:
        X = np.array([[0, 2], [1, 1]])

    ct_cloned = clone(ct).fit(X)
    out = estimator_html_repr(ct_cloned)
    assert feature_cols[0] in out
    assert feature_cols[1] in out


@pytest.mark.parametrize(
    "pandas, total_output_features",
    [
        (True, ["norm1__A", "norm1__B", "remainder__C"]),
        (False, ["norm1__x0", "norm1__x1", "remainder__x2"]),
    ],
)
def test_estimator_html_repr_total_feature_names(pandas, total_output_features):
    """Test features names are kept with pandas col names and generic."""
    if pandas:
        pd = pytest.importorskip("pandas")
        X = pd.DataFrame({"A": [0, 2, 3], "B": [1, 1, 3], "C": [3, 5, 4]})
    else:
        X = np.array([[0, 2, 3], [1, 1, 3], [3, 5, 4]])

    ct_cloned = clone(ct).fit(X)
    out = estimator_html_repr(ct_cloned)

    assert "<div class='total_features'>" in out
    assert "3 features</div>" in out
    for feature_name in total_output_features:
        assert feature_name in out


def test_estimator_html_col_names_featureunion():
    X = [[0.0, 1.0, 3], [2.0, 2.0, 5]]
    ct2_cloned = clone(ct2)
    ct2_cloned.fit_transform(X)
    out = estimator_html_repr(ct2_cloned)

    assert "pca__pca0" in out
    assert "svd__truncatedsvd0" in out
    assert "2 features" in out


def test_features_html_with_pipeline():
    """Test works with MinimalTransformer in a pipeline and scaler
    to test number of features."""

    pipe = Pipeline([("minimal", MinimalTransformer()), ("scaler", StandardScaler())])

    X = rng.randn(10, 3)
    pipe.fit(X)
    html = estimator_html_repr(pipe)
    assert "3 features" in html


def test_countvectorizer_output_features():
    """Non-regression test for
    https://github.com/scikit-learn/scikit-learn/issues/33772"""
    corpus = [
        "cat",
        "dog",
        "mouse",
        "bird",
    ]
    vectorizer = CountVectorizer()
    vectorizer.fit_transform(corpus)
    html = estimator_html_repr(vectorizer)
    assert "4 features" in html


def test_meta_estimator_output_features():
    """Non-regression test for
    https://github.com/scikit-learn/scikit-learn/issues/33887
    """
    pytest.importorskip("pandas")
    X, y = load_iris(return_X_y=True, as_frame=True)
    X, y = X.iloc[:100], y.iloc[:100]

    preprocessor = make_column_transformer(
        (StandardScaler(), [0, 1]),
        (MinMaxScaler(), [2, 3]),
    )
    estimator = make_pipeline(preprocessor, LogisticRegression())
    meta_estimator = TunedThresholdClassifierCV(
        estimator, store_cv_results=True, random_state=0
    ).fit(X, y)
    html = estimator_html_repr(meta_estimator)
    assert "4 features" in html


def test_get_feature_names_out_exception():
    """Non-regression test for
    https://github.com/scikit-learn/scikit-learn/issues/33887
    Testing that error in _get_feature_names_out doesn't break
    and we still get an HTML display with no number of features.
    """

    X = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])

    union = FeatureUnion(
        [
            ("pca", PCA(n_components=1)),
            ("identity", FunctionTransformer()),
        ]
    )
    Xt = union.fit_transform(X)
    html = estimator_html_repr(union)
    assert "1 feature" in html
    assert "<div> 0 features</div>" not in html
    assert "identity" in html


def test_single_estimator_get_feature_names_out_exception():
    """Non-regression test for
    https://github.com/scikit-learn/scikit-learn/issues/33887
    Testing that error in _get_feature_names_out doesn't break
    hitting single block except branch"""

    class BrokenTransformer(TransformerMixin, BaseEstimator):
        def fit(self, X, y=None):
            self.n_features_in_ = X.shape[1]
            return self

        def get_feature_names_out(self, input_features=None):
            raise RuntimeError("Simulated failure")

    X = np.array([[1, 2], [3, 4]])
    t = BrokenTransformer()
    t.fit(X)
    html = estimator_html_repr(t)
    assert "BrokenTransformer" in html


def test_features_html_empty_features():
    """Test that _features_html handles empty feature list."""
    features = []
    html = _features_html(features)

    assert "<tbody>" in html


def test_features_html_special_characters():
    """Test that special characters in feature names are properly escaped."""
    features = [
        "feature&1",
        'feature"2',
        "feature'3",
        "feature>4",
        "feature<5",
        "<script>alert('xss')</script>",
    ]
    html = _features_html(features)

    assert "&amp;" in html
    assert "&lt;script&gt;alert(&#x27;xss&#x27;)&lt;/script&gt;" in html
    assert "&gt;" in html
    assert "&lt;" in html
    assert "6 features" in html


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
