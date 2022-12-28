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


@pytest.mark.parametrize(
    "params, err_type, err_msg",
    [
        ({"std_display_style": "invalid"}, ValueError, "Unknown std_display_style:"),
        ({"score_type": "invalid"}, ValueError, "Unknown score_type:"),
    ],
)
def test_learning_curve_display_parameters_validation(
    pyplot, data, params, err_type, err_msg
):
    """Check that we raise a proper error when passing invalid parameters."""
    X, y = data
    estimator = DecisionTreeClassifier(random_state=0)

    train_sizes = [0.3, 0.6, 0.9]
    with pytest.raises(err_type, match=err_msg):
        LearningCurveDisplay.from_estimator(
            estimator, X, y, train_sizes=train_sizes, **params
        )


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
    assert legend_labels == ["Testing metric"]

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


@pytest.mark.parametrize(
    "score_name, ylabel", [(None, "Score"), ("Accuracy", "Accuracy")]
)
def test_learning_curve_display_score_name(pyplot, data, score_name, ylabel):
    """Check that we can overwrite the default score name shown on the y-axis."""
    X, y = data
    estimator = DecisionTreeClassifier(random_state=0)

    train_sizes = [0.3, 0.6, 0.9]
    display = LearningCurveDisplay.from_estimator(
        estimator, X, y, train_sizes=train_sizes, score_name=score_name
    )

    assert display.ax_.get_ylabel() == ylabel
    X, y = data
    estimator = DecisionTreeClassifier(max_depth=1, random_state=0)

    train_sizes = [0.3, 0.6, 0.9]
    display = LearningCurveDisplay.from_estimator(
        estimator, X, y, train_sizes=train_sizes, score_name=score_name
    )

    assert display.score_name == ylabel


@pytest.mark.parametrize("std_display_style", (None, "errorbar"))
def test_learning_curve_display_score_type(pyplot, data, std_display_style):
    """Check the behaviour of setting the `score_type` parameter."""
    X, y = data
    estimator = DecisionTreeClassifier(random_state=0)

    train_sizes = [0.3, 0.6, 0.9]
    train_sizes_abs, train_scores, test_scores = learning_curve(
        estimator, X, y, train_sizes=train_sizes
    )

    score_type = "train"
    display = LearningCurveDisplay.from_estimator(
        estimator,
        X,
        y,
        train_sizes=train_sizes,
        score_type=score_type,
        std_display_style=std_display_style,
    )

    _, legend_label = display.ax_.get_legend_handles_labels()
    assert legend_label == ["Training metric"]

    if std_display_style is None:
        assert len(display.lines_) == 1
        assert display.errorbar_ is None
        x_data, y_data = display.lines_[0].get_data()
    else:
        assert display.lines_ is None
        assert len(display.errorbar_) == 1
        x_data, y_data = display.errorbar_[0].lines[0].get_data()

    assert_array_equal(x_data, train_sizes_abs)
    assert_allclose(y_data, train_scores.mean(axis=1))

    score_type = "test"
    display = LearningCurveDisplay.from_estimator(
        estimator,
        X,
        y,
        train_sizes=train_sizes,
        score_type=score_type,
        std_display_style=std_display_style,
    )

    _, legend_label = display.ax_.get_legend_handles_labels()
    assert legend_label == ["Testing metric"]

    if std_display_style is None:
        assert len(display.lines_) == 1
        assert display.errorbar_ is None
        x_data, y_data = display.lines_[0].get_data()
    else:
        assert display.lines_ is None
        assert len(display.errorbar_) == 1
        x_data, y_data = display.errorbar_[0].lines[0].get_data()

    assert_array_equal(x_data, train_sizes_abs)
    assert_allclose(y_data, test_scores.mean(axis=1))

    score_type = "both"
    display = LearningCurveDisplay.from_estimator(
        estimator,
        X,
        y,
        train_sizes=train_sizes,
        score_type=score_type,
        std_display_style=std_display_style,
    )

    _, legend_label = display.ax_.get_legend_handles_labels()
    assert legend_label == ["Training metric", "Testing metric"]

    if std_display_style is None:
        assert len(display.lines_) == 2
        assert display.errorbar_ is None
        x_data_train, y_data_train = display.lines_[0].get_data()
        x_data_test, y_data_test = display.lines_[1].get_data()
    else:
        assert display.lines_ is None
        assert len(display.errorbar_) == 2
        x_data_train, y_data_train = display.errorbar_[0].lines[0].get_data()
        x_data_test, y_data_test = display.errorbar_[1].lines[0].get_data()

    assert_array_equal(x_data_train, train_sizes_abs)
    assert_allclose(y_data_train, train_scores.mean(axis=1))
    assert_array_equal(x_data_test, train_sizes_abs)
    assert_allclose(y_data_test, test_scores.mean(axis=1))


def test_learning_curve_display_log_scale(pyplot, data):
    """Check the behaviour of the parameter `log_scale`."""
    X, y = data
    estimator = DecisionTreeClassifier(random_state=0)

    train_sizes = [0.3, 0.6, 0.9]
    display = LearningCurveDisplay.from_estimator(
        estimator, X, y, train_sizes=train_sizes, log_scale=True
    )

    assert display.ax_.get_xscale() == "log"
    assert display.ax_.get_yscale() == "linear"

    display = LearningCurveDisplay.from_estimator(
        estimator, X, y, train_sizes=train_sizes, log_scale=False
    )

    assert display.ax_.get_xscale() == "linear"
    assert display.ax_.get_yscale() == "linear"


def test_learning_curve_display_std_display_style(pyplot, data):
    """Check the behaviour of the parameter `std_display_style`."""
    X, y = data
    estimator = DecisionTreeClassifier(random_state=0)

    import matplotlib as mpl

    train_sizes = [0.3, 0.6, 0.9]
    std_display_style = None
    display = LearningCurveDisplay.from_estimator(
        estimator,
        X,
        y,
        train_sizes=train_sizes,
        std_display_style=std_display_style,
    )

    assert len(display.lines_) == 1
    assert isinstance(display.lines_[0], mpl.lines.Line2D)
    assert display.errorbar_ is None
    assert display.fill_between_ is None
    _, legend_label = display.ax_.get_legend_handles_labels()
    assert len(legend_label) == 1

    std_display_style = "fill_between"
    display = LearningCurveDisplay.from_estimator(
        estimator,
        X,
        y,
        train_sizes=train_sizes,
        std_display_style=std_display_style,
    )

    assert len(display.lines_) == 1
    assert isinstance(display.lines_[0], mpl.lines.Line2D)
    assert display.errorbar_ is None
    assert len(display.fill_between_) == 1
    assert isinstance(display.fill_between_[0], mpl.collections.PolyCollection)
    _, legend_label = display.ax_.get_legend_handles_labels()
    assert len(legend_label) == 1

    std_display_style = "errorbar"
    display = LearningCurveDisplay.from_estimator(
        estimator,
        X,
        y,
        train_sizes=train_sizes,
        std_display_style=std_display_style,
    )

    assert display.lines_ is None
    assert len(display.errorbar_) == 1
    assert isinstance(display.errorbar_[0], mpl.container.ErrorbarContainer)
    assert display.fill_between_ is None
    _, legend_label = display.ax_.get_legend_handles_labels()
    assert len(legend_label) == 1


def test_learning_curve_display_plot_kwargs(pyplot, data):
    """Check the behaviour of the different plotting keyword arguments: `line_kw`,
    `fill_between_kw`, and `errorbar_kw`."""
    X, y = data
    estimator = DecisionTreeClassifier(random_state=0)

    train_sizes = [0.3, 0.6, 0.9]
    std_display_style = "fill_between"
    line_kw = {"color": "red"}
    fill_between_kw = {"color": "red", "alpha": 1.0}
    display = LearningCurveDisplay.from_estimator(
        estimator,
        X,
        y,
        train_sizes=train_sizes,
        std_display_style=std_display_style,
        line_kw=line_kw,
        fill_between_kw=fill_between_kw,
    )

    assert display.lines_[0].get_color() == "red"
    assert_allclose(
        display.fill_between_[0].get_facecolor(),
        [[1.0, 0.0, 0.0, 1.0]],  # trust me, it's red
    )

    std_display_style = "errorbar"
    errorbar_kw = {"color": "red"}
    display = LearningCurveDisplay.from_estimator(
        estimator,
        X,
        y,
        train_sizes=train_sizes,
        std_display_style=std_display_style,
        errorbar_kw=errorbar_kw,
    )

    assert display.errorbar_[0].lines[0].get_color() == "red"
