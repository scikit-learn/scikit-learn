import numpy as np
import pytest

from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split

from ..cap_curve import CumulativeAccuracyDisplay


@pytest.fixture
def binary_classification_dataset():
    X, y = make_classification(
        n_samples=100,
        n_features=2,
        n_informative=2,
        n_redundant=0,
        n_repeated=0,
        n_classes=2,
        random_state=42,
    )
    return train_test_split(X, y, test_size=0.2, random_state=42)


@pytest.fixture
def logistic_regression_model(binary_classification_dataset):
    X_train, _, y_train, _ = binary_classification_dataset
    clf = LogisticRegression(max_iter=1000)
    clf.fit(X_train, y_train)
    return clf


@pytest.mark.parametrize("normalize_scale", [True, False])
@pytest.mark.parametrize("plot_chance_level", [True, False])
def test_cumulative_accuracy_display_from_predictions(
    binary_classification_dataset, normalize_scale, plot_chance_level
):
    _, X_test, _, y_test = binary_classification_dataset
    y_scores = np.random.rand(len(y_test))

    cap_display = CumulativeAccuracyDisplay.from_predictions(
        y_test,
        y_scores,
        normalize_scale=normalize_scale,
        plot_chance_level=plot_chance_level,
        name="Test Classifier",
    )

    assert cap_display is not None
    assert hasattr(cap_display, "line_"), "The display must have a line attribute"
    assert hasattr(cap_display, "ax_"), "The display must have an ax attribute"
    assert hasattr(cap_display, "figure_"), "The display must have a figure attribute"
    if plot_chance_level:
        assert (
            cap_display.chance_level_ is not None
        ), "Chance level line should be present"


@pytest.mark.parametrize("normalize_scale", [True, False])
@pytest.mark.parametrize("plot_chance_level", [True, False])
def test_cumulative_accuracy_display_from_estimator(
    logistic_regression_model,
    binary_classification_dataset,
    normalize_scale,
    plot_chance_level,
):
    _, X_test, _, y_test = binary_classification_dataset

    cap_display = CumulativeAccuracyDisplay.from_estimator(
        logistic_regression_model,
        X_test,
        y_test,
        normalize_scale=normalize_scale,
        plot_chance_level=plot_chance_level,
        name="Logistic Regression",
    )

    assert cap_display is not None
    assert hasattr(cap_display, "line_"), "The display must have a line attribute"
    assert hasattr(cap_display, "ax_"), "The display must have an ax attribute"
    assert hasattr(cap_display, "figure_"), "The display must have a figure attribute"
    if plot_chance_level:
        assert (
            cap_display.chance_level_ is not None
        ), "Chance level line should be present"
