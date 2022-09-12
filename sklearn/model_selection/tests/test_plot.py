import pytest

from sklearn.datasets import load_iris
from sklearn.tree import DecisionTreeClassifier
from sklearn.utils import shuffle
from sklearn.utils._testing import assert_allclose, assert_array_equal

from sklearn.model_selection import learning_curve
from sklearn.model_selection import LearningCurveDisplay


@pytest.fixture
def data():
    return shuffle(*load_iris(return_X_y=True), random_state=0)


def test_learnig_curve_display_default_usage(pyplot, data):
    """Check the default usage of the LearningCurveDisplay class."""
    X, y = data
    estimator = DecisionTreeClassifier(random_state=0)

    train_sizes = [0.3, 0.6, 0.9]
    display = LearningCurveDisplay.from_estimator(
        estimator, X, y, train_sizes=train_sizes
    )

    import matplotlib as mpl

    assert display.errorbar_ is None
    assert isinstance(display.lines_, list)
    for line in display.lines_:
        assert isinstance(line, mpl.lines.Line2D)
    assert isinstance(display.fill_between_, list)
    for fill in display.fill_between_:
        assert isinstance(fill, mpl.collections.PolyCollection)
    assert display.score_name == "Score"
    assert display.ax_.get_xlabel() == "Number of samples in the training set"
    assert display.ax_.get_ylabel() == "Score"

    train_sizes_abs, train_scores, test_scores = learning_curve(
        estimator, X, y, train_sizes=train_sizes
    )

    assert_array_equal(display.train_sizes, train_sizes_abs)
    assert_allclose(display.train_scores, train_scores)
    assert_allclose(display.test_scores, test_scores)
