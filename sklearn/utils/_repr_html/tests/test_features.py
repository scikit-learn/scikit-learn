import numpy as np

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import Normalizer
from sklearn.utils._repr_html.estimator import estimator_html_repr
from sklearn.utils._testing import MinimalTransformer

ct = ColumnTransformer([("norm1", Normalizer(), [0, 1])])


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
