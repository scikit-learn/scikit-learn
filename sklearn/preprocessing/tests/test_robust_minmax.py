import numpy as np

from sklearn.preprocessing import robust_minmax_scale

def test_basic_scaling():
    X = np.array([[1], [5], [100]])
    X_scaled = robust_minmax_scale(X)
    assert X_scaled.min() >= 0
    assert X_scaled.max() <= 1


def test_resistant_to_outliers():
    X = np.array([[1], [2], [3], [4], [1000]])
    X_scaled = robust_minmax_scale(X)
    # Check that the outlier is clipped at 1
    assert X_scaled[-1] == 1.0
