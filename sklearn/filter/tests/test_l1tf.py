"""L1 trend filtering test suite"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import pytest
import numpy as np
from sklearn.pipeline import Pipeline
from sklearn.filter import L1TrendFiltering


@pytest.fixture
def sample_data():
    """
    Provides sample data for testing.
    """
    np.random.seed(42)
    n = 100
    x = np.linspace(0, 10, n)
    y = np.sin(x) + 0.5 * np.random.randn(n)  # Sine wave with noise
    return y


def test_l1_trend_filtering_basic(sample_data):
    """
    Test the basic functionality of the L1TrendFiltering transformer.
    """
    transformer = L1TrendFiltering(lambda_=10)
    result = transformer.transform(sample_data)

    assert result.shape == sample_data.shape, "Output shape should match input shape"
    assert not np.allclose(result, sample_data), "Filtered data should not exactly match noisy input"
    assert np.all(np.isfinite(result)), "All output values should be finite"


def test_l1_trend_filtering_pipeline(sample_data):
    """
    Test L1TrendFiltering within a scikit-learn pipeline.
    """
    pipeline = Pipeline([
        ('l1_filter', L1TrendFiltering(lambda_=10)),
    ])
    result = pipeline.fit_transform(sample_data)

    assert result.shape == sample_data.shape, "Pipeline output shape should match input shape"
    assert not np.allclose(result, sample_data), "Filtered data should not exactly match noisy input"


def test_l1_trend_filtering_with_zero_lambda(sample_data):
    """
    Test L1TrendFiltering with lambda set to 0, which should ideally return the input unchanged.
    """
    transformer = L1TrendFiltering(lambda_=0)
    result = transformer.transform(sample_data)

    assert np.allclose(result, sample_data), "With lambda=0, the output should match the input"


def test_l1_trend_filtering_extreme_lambda(sample_data):
    """
    Test L1TrendFiltering with a very large lambda value, which should produce a flat line.
    """
    transformer = L1TrendFiltering(lambda_=1e6)
    result = transformer.transform(sample_data)

    assert np.allclose(result, np.mean(sample_data) * np.ones_like(sample_data), atol=1e-2), \
        "With large lambda, the output should approach a flat line"


def test_l1_trend_filtering_invalid_lambda():
    """
    Test L1TrendFiltering with invalid lambda values.
    """
    with pytest.raises(ValueError, match=".*lambda.*positive.*"):
        L1TrendFiltering(lambda_=-1)


def test_l1_trend_filtering_tolerance(sample_data):
    """
    Test L1TrendFiltering with a custom tolerance.
    """
    transformer = L1TrendFiltering(lambda_=10, tol=1e-6)
    result = transformer.transform(sample_data)

    assert result.shape == sample_data.shape, "Output shape should match input shape"
    assert np.all(np.isfinite(result)), "All output values should be finite"


def test_l1_trend_filtering_multidimensional_input():
    """
    Test L1TrendFiltering with multidimensional input.
    """
    data = np.random.randn(100, 3)  # 3 features
    transformer = L1TrendFiltering(lambda_=10)
    result = transformer.transform(data)

    assert result.shape == data.shape, "Output shape should match input shape for multidimensional input"
    assert np.all(np.isfinite(result)), "All output values should be finite"


def test_l1_trend_filtering_zero_input():
    """
    Test L1TrendFiltering with an input of all zeros.
    """
    data = np.zeros(100)
    transformer = L1TrendFiltering(lambda_=10)
    result = transformer.transform(data)

    assert np.allclose(result, 0), "Output should remain zero for zero input"


def test_l1_trend_filtering_short_signal():
    """
    Test L1TrendFiltering with a very short signal (edge case).
    """
    data = np.array([1.0, 2.0, 3.0])
    transformer = L1TrendFiltering(lambda_=10)
    result = transformer.transform(data)

    assert result.shape == data.shape, "Output shape should match input shape for short signals"


def test_l1_trend_filtering_large_input():
    """
    Test L1TrendFiltering with a very large input signal.
    """
    data = np.random.randn(10000)  # Large input
    transformer = L1TrendFiltering(lambda_=10)
    result = transformer.transform(data)

    assert result.shape == data.shape, "Output shape should match input shape for large signals"
    assert np.all(np.isfinite(result)), "All output values should be finite"
