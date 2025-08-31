import numpy as np
import pytest

from sklearn.datasets import make_classification
from sklearn.inspection import DecisionBoundaryDisplay
from sklearn.linear_model import LogisticRegression


@pytest.mark.parametrize(
    "response_method",
    ["predict", "predict_proba", "decision_function"],
)
def test_decision_boundary_display_plotting_with_sample_weight(response_method):
    """Check that DecisionBoundaryDisplay.plot works with sample_weight (gh-27462)."""
    X, y = make_classification(
        n_samples=50, n_features=2, n_classes=2, random_state=42
    )
    sample_weight = np.random.RandomState(0).rand(y.shape[0])

    clf = LogisticRegression().fit(X, y, sample_weight=sample_weight)

    disp = DecisionBoundaryDisplay.from_estimator(
        clf,
        X,
        response_method=response_method,
    )
    disp.plot(X, y, sample_weight=sample_weight)
