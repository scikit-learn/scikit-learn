"""Tests for permutation_importance with Polars DataFrames."""

import numpy as np
import pytest

from sklearn import set_config
from sklearn.compose import ColumnTransformer
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

pl = pytest.importorskip("polars")


def test_permutation_importance_polars_dataframe():
    """Test that permutation_importance works with Polars DataFrames."""
    set_config(transform_output="polars")

    # Create sample data
    X = np.random.randn(100, 5)
    y = (X[:, 0] + X[:, 1] > 0).astype(int)
    feature_names = [f"feature_{i}" for i in range(5)]

    X_pl = pl.DataFrame({name: X[:, i] for i, name in enumerate(feature_names)})

    # Create pipeline with ColumnTransformer using string column names
    preprocessor = ColumnTransformer(
        transformers=[("scaler", StandardScaler(), feature_names)]
    )

    pipeline = Pipeline(
        [
            ("preprocessor", preprocessor),
            ("classifier", LogisticRegression(random_state=42)),
        ]
    )

    # Fit model
    pipeline.fit(X_pl, y)

    # Calculate permutation importance - should not raise
    result = permutation_importance(pipeline, X_pl, y, n_repeats=5, random_state=42)

    # Check result structure
    assert hasattr(result, "importances_mean")
    assert hasattr(result, "importances_std")
    assert hasattr(result, "importances")
    assert result.importances_mean.shape == (5,)
    assert result.importances.shape == (5, 5)
