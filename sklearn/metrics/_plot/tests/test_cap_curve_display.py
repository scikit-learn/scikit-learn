import matplotlib.pyplot as plt
import numpy as np
import pytest

from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split

from ..cap_curve import CAPCurveDisplay


@pytest.fixture
def data_binary():
    X, y = make_classification(
        n_samples=100,
        n_features=2,
        n_informative=2,
        n_redundant=0,
        n_repeated=0,
        n_classes=2,
        random_state=42,
    )
    return X, y


@pytest.fixture
def logistic_regression_model(data_binary):
    X, y = data_binary
    X_train, _, y_train, _ = train_test_split(X, y, test_size=0.2, random_state=42)
    clf = LogisticRegression(max_iter=1000)
    clf.fit(X_train, y_train)
    return clf


@pytest.mark.parametrize("normalize_scale", [True, False])
@pytest.mark.parametrize("plot_chance_level", [True, False])
@pytest.mark.parametrize("name", ["Logistic Regression", None])
@pytest.mark.parametrize("chance_level_kw", [{"color": "red", "lw": 3}, None])
def test_cumulative_accuracy_display_from_predictions(
    data_binary,
    normalize_scale,
    plot_chance_level,
    name,
    chance_level_kw,
):
    X, y = data_binary
    y_scores = np.random.rand(len(y))

    cap_display = CAPCurveDisplay.from_predictions(
        y,
        y_scores,
        normalize_scale=normalize_scale,
        plot_chance_level=plot_chance_level,
        name=name,
        chance_level_kw=chance_level_kw,
    )

    assert cap_display is not None
    assert hasattr(cap_display, "line_"), "The display must have a line attribute"
    assert hasattr(cap_display, "ax_"), "The display must have an ax attribute"
    assert hasattr(cap_display, "figure_"), "The display must have a figure attribute"
    assert hasattr(
        cap_display, "pos_label"
    ), "The display must have a pos_label attribute"
    assert hasattr(
        cap_display, "y_true_cumulative"
    ), "The display must have a y_true_cumulative attribute"
    if plot_chance_level:
        assert (
            cap_display.chance_level_ is not None
        ), "Chance level line should be present"

    # Matplotlib raises a warning when opening a large number
    # of figures without closing them.
    plt.close()


@pytest.mark.parametrize(
    "response_method", ["auto", "predict_proba", "decision_function"]
)
@pytest.mark.parametrize("normalize_scale", [True, False])
@pytest.mark.parametrize("plot_chance_level", [True, False])
@pytest.mark.parametrize("name", ["Logistic Regression", None])
def test_cumulative_accuracy_display_from_estimator(
    logistic_regression_model,
    data_binary,
    normalize_scale,
    plot_chance_level,
    name,
    response_method,
):
    X, y = data_binary

    cap_display = CAPCurveDisplay.from_estimator(
        logistic_regression_model,
        X,
        y,
        normalize_scale=normalize_scale,
        plot_chance_level=plot_chance_level,
        name=name,
        response_method=response_method,
    )

    assert cap_display is not None
    assert hasattr(cap_display, "line_"), "The display must have a line attribute"
    assert hasattr(cap_display, "ax_"), "The display must have an ax attribute"
    assert hasattr(cap_display, "figure_"), "The display must have a figure attribute"
    assert hasattr(
        cap_display, "y_true_cumulative"
    ), "The display must have a y_true_cumulative attribute"
    assert hasattr(
        cap_display.y_true_cumulative, "shape"
    ), "y_true_cumulative must have a shape attribute"
    assert hasattr(
        cap_display.y_true_cumulative, "dtype"
    ), "y_true_cumulative must have a dtype attribute"
    if plot_chance_level:
        assert (
            cap_display.chance_level_ is not None
        ), "Chance level line should be present"

    # Matplotlib raises a warning when opening a large number
    # of figures without closing them.
    plt.close()


@pytest.mark.parametrize(
    "response_method", ["auto", "predict_proba", "decision_function"]
)
@pytest.mark.parametrize("normalize_scale", [True, False])
@pytest.mark.parametrize("name", ["Logistic Regression", None])
def test_display_from_estimator_and_from_prediction(
    logistic_regression_model,
    data_binary,
    normalize_scale,
    name,
    response_method,
):
    X, y = data_binary
    y_pred = logistic_regression_model.decision_function(X)

    cap_display_from_estimator = CAPCurveDisplay.from_estimator(
        logistic_regression_model,
        X,
        y,
        normalize_scale=normalize_scale,
        name=name,
        response_method=response_method,
    )

    cap_display_from_prediction = CAPCurveDisplay.from_predictions(
        y,
        y_pred,
        normalize_scale=normalize_scale,
        name=name,
    )

    np.testing.assert_allclose(
        cap_display_from_estimator.y_true_cumulative,
        cap_display_from_prediction.y_true_cumulative,
        rtol=1e-6,
        atol=1e-8,
    )


def test_invalid_response_method(
    logistic_regression_model,
    data_binary,
):
    with pytest.raises(ValueError):
        X, y = data_binary

        _ = CAPCurveDisplay.from_estimator(
            logistic_regression_model, X, y, response_method="invalid input"
        )


@pytest.mark.parametrize("plot_chance_level", [True, False])
@pytest.mark.parametrize(
    "chance_level_kw",
    [
        None,
        {"linewidth": 1, "color": "red", "linestyle": "-", "label": "DummyEstimator"},
        {"lw": 1, "c": "red", "ls": "-", "label": "DummyEstimator"},
    ],
)
@pytest.mark.parametrize(
    "constructor_name",
    ["from_estimator", "from_predictions"],
)
def test_cap_curve_chance_level_line(
    logistic_regression_model,
    data_binary,
    constructor_name,
    plot_chance_level,
    chance_level_kw,
):
    X, y = data_binary
    y_pred = getattr(logistic_regression_model, "predict_proba")(X)
    y_pred = y_pred if y_pred.ndim == 1 else y_pred[:, 1]

    if constructor_name == "from_estimator":
        display = CAPCurveDisplay.from_estimator(
            logistic_regression_model,
            X,
            y,
            alpha=0.8,
            plot_chance_level=plot_chance_level,
            chance_level_kw=chance_level_kw,
        )
    else:
        display = CAPCurveDisplay.from_predictions(
            y,
            y_pred,
            alpha=0.8,
            plot_chance_level=plot_chance_level,
            chance_level_kw=chance_level_kw,
        )

    import matplotlib as mpl  # noqa

    assert isinstance(display.line_, mpl.lines.Line2D)
    assert display.line_.get_alpha() == 0.8
    assert isinstance(display.ax_, mpl.axes.Axes)
    assert isinstance(display.figure_, mpl.figure.Figure)

    if plot_chance_level:
        assert isinstance(display.chance_level_, mpl.lines.Line2D)
        assert tuple(display.chance_level_.get_xdata()) == (0, 1)
        assert tuple(display.chance_level_.get_ydata()) == (0, 1)
    else:
        assert display.chance_level_ is None

    # Checking for chance level line styles
    if plot_chance_level and chance_level_kw is None:
        assert display.chance_level_.get_color() == "k"
        assert display.chance_level_.get_linestyle() == "--"
        assert display.chance_level_.get_label() == "Chance level"
    elif plot_chance_level:
        assert display.chance_level_.get_label() == chance_level_kw["label"]
        if "c" in chance_level_kw:
            assert display.chance_level_.get_color() == chance_level_kw["c"]
        else:
            assert display.chance_level_.get_color() == chance_level_kw["color"]
        if "lw" in chance_level_kw:
            assert display.chance_level_.get_linewidth() == chance_level_kw["lw"]
        else:
            assert display.chance_level_.get_linewidth() == chance_level_kw["linewidth"]
        if "ls" in chance_level_kw:
            assert display.chance_level_.get_linestyle() == chance_level_kw["ls"]
        else:
            assert display.chance_level_.get_linestyle() == chance_level_kw["linestyle"]


@pytest.mark.parametrize("response_method", ["predict_proba", "decision_function"])
@pytest.mark.parametrize("with_sample_weight", [True, False])
@pytest.mark.parametrize("with_strings", [True, False])
@pytest.mark.parametrize(
    "constructor_name, default_name",
    [
        ("from_estimator", "LogisticRegression"),
        ("from_predictions", "Classifier"),
    ],
)
def test_roc_curve_display_plotting(
    logistic_regression_model,
    data_binary,
    response_method,
    constructor_name,
    default_name,
    with_sample_weight,
    with_strings,
):
    X, y = data_binary

    pos_label = None
    if with_strings:
        y = np.array(["c", "b"])[y]
        pos_label = "c"

    if with_sample_weight:
        rng = np.random.RandomState(42)
        sample_weight = rng.randint(1, 4, size=(X.shape[0]))
    else:
        sample_weight = None

    lr = LogisticRegression()
    lr.fit(X, y)

    y_pred = getattr(lr, response_method)(X)
    y_pred = y_pred if y_pred.ndim == 1 else y_pred[:, 1]

    if constructor_name == "from_estimator":
        display = CAPCurveDisplay.from_estimator(
            lr,
            X,
            y,
            pos_label=pos_label,
            sample_weight=sample_weight,
            alpha=0.8,
        )
    else:
        display = CAPCurveDisplay.from_predictions(
            y,
            y_pred,
            pos_label=pos_label,
            sample_weight=sample_weight,
            alpha=0.8,
        )

    assert (
        display.estimator_name == default_name
    ), f"{display.estimator_name} != {default_name}"
