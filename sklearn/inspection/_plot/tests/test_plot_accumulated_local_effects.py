import pytest
import numpy as np

from sklearn.datasets import make_classification, make_regression
from sklearn.inspection import AccumulatedLocalEffectsDisplay
from sklearn.linear_model import LinearRegression, LogisticRegression


def test_ale_display_from_estimator():
    """Verify AccumulatedLocalEffectsDisplay.from_estimator creates valid plots."""
    pytest.importorskip("matplotlib")
    import matplotlib.pyplot as plt

    rng = np.random.RandomState(42)
    X, y = make_regression(n_samples=150, n_features=3, random_state=rng)
    model = LinearRegression().fit(X, y)

    disp = AccumulatedLocalEffectsDisplay.from_estimator(model, X, features=[0, 1])

    assert len(disp.lines_) == 2
    assert disp.axes_.shape == (1, 2)
    assert disp.figure_ is not None

    plt.close("all")


def test_ale_display_custom_ax():
    """Verify plotting into a pre-existing matplotlib axes."""
    pytest.importorskip("matplotlib")
    import matplotlib.pyplot as plt

    rng = np.random.RandomState(42)
    X, y = make_regression(n_samples=100, n_features=2, random_state=rng)
    model = LinearRegression().fit(X, y)

    fig, ax = plt.subplots(figsize=(6, 4))
    disp = AccumulatedLocalEffectsDisplay.from_estimator(model, X, features=[0], ax=ax)

    assert disp.lines_[0] in ax.get_lines()
    assert ax.get_xlabel() == "x0"
    assert ax.get_ylabel() == "Accumulated Local Effect"

    plt.close("all")


def test_ale_display_multiclass():
    """Verify target parameter selects correct class in multiclass display."""
    pytest.importorskip("matplotlib")
    import matplotlib.pyplot as plt

    rng = np.random.RandomState(42)
    X, y = make_classification(
        n_samples=150, n_features=5, n_informative=3, n_redundant=1, n_classes=3, n_clusters_per_class=1, random_state=rng
    )
    clf = LogisticRegression().fit(X, y)

    disp = AccumulatedLocalEffectsDisplay.from_estimator(clf, X, features=[0], target=1)
    assert disp.target_idx == 1
    assert len(disp.lines_) == 1

    plt.close("all")
