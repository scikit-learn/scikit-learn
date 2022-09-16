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


def test_learning_curve_display_default_usage(pyplot, data):
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
        assert fill.get_alpha() == 0.5

    assert display.score_name == "Score"
    assert display.ax_.get_xlabel() == "Number of samples in the training set"
    assert display.ax_.get_ylabel() == "Score"

    _, legend_labels = display.ax_.get_legend_handles_labels()
    # one label for the line and the fill-between
    assert legend_labels == ["Test metric", "Test metric"]

    train_sizes_abs, train_scores, test_scores = learning_curve(
        estimator, X, y, train_sizes=train_sizes
    )

    assert_array_equal(display.train_sizes, train_sizes_abs)
    assert_allclose(display.train_scores, train_scores)
    assert_allclose(display.test_scores, test_scores)


def test_learning_curve_display_negate_score(pyplot, data):
    """Check the behaviour of the `negate_score` parameter calling `from_estimator` and
    `plot`.
    """
    X, y = data
    estimator = DecisionTreeClassifier(max_depth=1, random_state=0)

    train_sizes = [0.3, 0.6, 0.9]
    negate_score = False
    display = LearningCurveDisplay.from_estimator(
        estimator,
        X,
        y,
        train_sizes=train_sizes,
        negate_score=negate_score,
    )

    positive_scores = display.lines_[0].get_data()[1]
    assert (positive_scores >= 0).all()
    assert display.ax_.get_ylabel() == "Score"

    negate_score = True
    display = LearningCurveDisplay.from_estimator(
        estimator, X, y, train_sizes=train_sizes, negate_score=negate_score
    )

    negative_scores = display.lines_[0].get_data()[1]
    assert (negative_scores <= 0).all()
    assert_allclose(negative_scores, -positive_scores)
    assert display.ax_.get_ylabel() == "Score"

    negate_score = False
    display = LearningCurveDisplay.from_estimator(
        estimator,
        X,
        y,
        train_sizes=train_sizes,
        negate_score=negate_score,
    )
    assert display.ax_.get_ylabel() == "Score"
    display.plot(negate_score=not negate_score)
    assert display.ax_.get_ylabel() == "Score"
    assert (display.lines_[0].get_data()[1] < 0).all()
