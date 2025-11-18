import numpy as np

from sklearn.metrics import accuracy_score, decision_threshold_curve


def test_decision_threshold_curve():
    """Dummy test, just to check function works."""
    y_true = np.array([0, 0, 1, 1])
    y_score = np.array([0.1, 0.4, 0.35, 0.8])
    scores, thresholds = decision_threshold_curve(accuracy_score, y_true, y_score)
