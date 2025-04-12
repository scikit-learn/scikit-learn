import re

import numpy as np
import pytest

from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.datasets import make_classification
from sklearn.exceptions import NotFittedError
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.svm import LinearSVC

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
@pytest.mark.parametrize("y_pred_dtype", [np.float32, np.float64])
def test_cumulative_accuracy_display_from_predictions(
    pyplot,
    data_binary,
    normalize_scale,
    plot_chance_level,
    name,
    chance_level_kw,
    y_pred_dtype,
):
    _, y = data_binary
    y_scores = np.random.rand(len(y)).astype(y_pred_dtype)

    cap_display = CAPCurveDisplay.from_predictions(
        y,
        y_scores,
        normalize_scale=normalize_scale,
        plot_chance_level=plot_chance_level,
        name=name,
        chance_level_kw=chance_level_kw,
    )

    assert cap_display is not None
    assert hasattr(
        cap_display, "y_true_cumulative"
    ), "The display must have a y_true_cumulative attribute"
    assert hasattr(
        cap_display, "cumulative_total"
    ), "The display must have a cumulative_total attribute"
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

    y_true_cumulative = cap_display.y_true_cumulative
    cumulative_total = cap_display.cumulative_total
    assert y_true_cumulative.dtype == cumulative_total.dtype, (
        f"Dtype mismatch: y_true_cumulative.dtype={y_true_cumulative.dtype}, "
        f"cumulative_total.dtype={cumulative_total.dtype}"
    )


@pytest.mark.parametrize(
    "response_method", ["auto", "predict_proba", "decision_function"]
)
@pytest.mark.parametrize("normalize_scale", [True, False])
@pytest.mark.parametrize("plot_chance_level", [True, False])
@pytest.mark.parametrize("name", ["Logistic Regression", None])
def test_cumulative_accuracy_display_from_estimator(
    pyplot,
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


@pytest.mark.parametrize(
    "response_method", ["auto", "predict_proba", "decision_function"]
)
@pytest.mark.parametrize("normalize_scale", [True, False])
@pytest.mark.parametrize("name", ["Logistic Regression", None])
def test_display_from_estimator_and_from_prediction(
    pyplot,
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


def test_invalid_response_method(pyplot, logistic_regression_model, data_binary):
    X, y = data_binary
    pattern = r"response_method must be in:.*"

    with pytest.raises(AttributeError, match=pattern):
        _ = CAPCurveDisplay.from_estimator(
            logistic_regression_model, X, y, response_method="invalid input"
        )


def test_unfitted_estimator(pyplot, data_binary):
    pattern = (
        r"This .* instance is not fitted yet\. Call 'fit' with "
        r"appropriate arguments before using this estimator\."
    )
    with pytest.raises(NotFittedError, match=pattern):
        X, y = data_binary

        lr_model_nofit = LogisticRegression(max_iter=1000)
        _ = CAPCurveDisplay.from_estimator(lr_model_nofit, X, y)


def test_estimator_has_too_many_classes(pyplot):
    X, y = make_classification(
        n_samples=100,
        n_features=5,
        n_informative=3,
        n_classes=3,
        random_state=42,
    )

    lr_model_n_classes = LogisticRegression(max_iter=1000)
    lr_model_n_classes.fit(X, y)

    pattern = re.escape(
        "Expected 'estimator' to be a binary classifier. Got 3 classes instead."
    )

    with pytest.raises(ValueError, match=pattern):
        _ = CAPCurveDisplay.from_estimator(lr_model_n_classes, X, y)


@pytest.mark.parametrize("plot_chance_level", [True, False])
@pytest.mark.parametrize("plot_perfect", [True, False])
@pytest.mark.parametrize(
    "chance_level_kw",
    [
        None,
        {"linewidth": 1, "color": "red", "linestyle": "-", "label": "DummyEstimator"},
        {"lw": 1, "c": "red", "ls": "-", "label": "DummyEstimator"},
    ],
)
@pytest.mark.parametrize(
    "perfect_level_kw",
    [
        None,
        {
            "linewidth": 2,
            "color": "green",
            "linestyle": "-.",
            "label": "PerfectEstimator",
        },
        {"lw": 1.5, "c": "blue", "ls": "-", "label": "PerfectEstimator"},
    ],
)
@pytest.mark.parametrize(
    "constructor_name",
    ["from_estimator", "from_predictions"],
)
@pytest.mark.parametrize("pos_label", [1, 0, None])
def test_cap_curve_chance_level_line(
    pyplot,
    logistic_regression_model,
    data_binary,
    constructor_name,
    plot_chance_level,
    chance_level_kw,
    plot_perfect,
    perfect_level_kw,
    pos_label,
):
    X, y = data_binary
    y_pred = getattr(logistic_regression_model, "predict_proba")(X)
    if y_pred.ndim != 1:
        if pos_label is not None:
            class_index = np.where(logistic_regression_model.classes_ == pos_label)[0][
                0
            ]
        else:
            class_index = -1
        y_pred = y_pred[:, class_index]

    if constructor_name == "from_estimator":
        display = CAPCurveDisplay.from_estimator(
            logistic_regression_model,
            X,
            y,
            alpha=0.8,
            plot_chance_level=plot_chance_level,
            chance_level_kw=chance_level_kw,
            plot_perfect=plot_perfect,
            perfect_level_kw=perfect_level_kw,
        )
    else:
        display = CAPCurveDisplay.from_predictions(
            y,
            y_pred,
            alpha=0.8,
            plot_chance_level=plot_chance_level,
            chance_level_kw=chance_level_kw,
            plot_perfect=plot_perfect,
            perfect_level_kw=perfect_level_kw,
        )

    import matplotlib as mpl

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

    if plot_perfect:
        assert isinstance(display.perfect_level_, mpl.lines.Line2D)
        xdata = tuple(display.perfect_level_.get_xdata())
        assert xdata[0] == 0
        assert xdata[-1] == 1
        ydata = tuple(display.perfect_level_.get_ydata())
        assert ydata[0] == 0
        assert ydata[-1] == 1
    else:
        assert display.perfect_level_ is None

    if plot_perfect and perfect_level_kw is None:
        assert display.perfect_level_.get_color() == "black"
        assert display.perfect_level_.get_linestyle() == ":"
        assert display.perfect_level_.get_label() == "Perfect predictions"
    elif plot_perfect:
        assert display.perfect_level_.get_label() == perfect_level_kw["label"]
        if "c" in perfect_level_kw:
            assert display.perfect_level_.get_color() == perfect_level_kw["c"]
        else:
            assert display.perfect_level_.get_color() == perfect_level_kw["color"]
        if "lw" in perfect_level_kw:
            assert display.perfect_level_.get_linewidth() == perfect_level_kw["lw"]
        else:
            assert (
                display.perfect_level_.get_linewidth() == perfect_level_kw["linewidth"]
            )
        if "ls" in perfect_level_kw:
            assert display.perfect_level_.get_linestyle() == perfect_level_kw["ls"]
        else:
            assert (
                display.perfect_level_.get_linestyle() == perfect_level_kw["linestyle"]
            )


@pytest.mark.parametrize(
    "response_method", ["predict_proba", "decision_function", "auto"]
)
@pytest.mark.parametrize("with_sample_weight", [True, False])
@pytest.mark.parametrize(
    ("pos_label", "with_strings"),
    [
        ("a", True),
        ("b", True),
        (None, False),
    ],
)
@pytest.mark.parametrize(
    "constructor_name, default_name",
    [
        ("from_estimator", "LogisticRegression"),
        ("from_predictions", "Classifier"),
    ],
)
def test_cap_curve_between_chance_and_perfect(
    pyplot,
    data_binary,
    response_method,
    constructor_name,
    default_name,
    with_sample_weight,
    with_strings,
    pos_label,
):
    X, y = data_binary

    if with_strings:
        y = np.array(["a", "b"])[y]

    if with_sample_weight:
        rng = np.random.RandomState(42)
        sample_weight = rng.randint(1, 4, size=(X.shape[0]))
    else:
        sample_weight = None

    lr = LogisticRegression()
    lr.fit(X, y)

    if response_method == "auto":
        if hasattr(lr, "predict_proba"):
            response_method = "predict_proba"
        elif hasattr(lr, "decision_function"):
            response_method = "decision_function"

    y_pred = getattr(lr, response_method)(X)
    if y_pred.ndim != 1:
        if pos_label is not None:
            class_index = np.where(lr.classes_ == pos_label)[0][0]
        else:
            class_index = -1
        y_pred = y_pred[:, class_index]
    else:
        if pos_label is not None:
            if pos_label != lr.classes_[1]:
                y_pred = -y_pred

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

    # Check that CAP curve is on or above the chance level line
    cumulative_total = display.cumulative_total
    y_true_cumulative = display.y_true_cumulative
    assert np.all(y_true_cumulative >= cumulative_total - 1e-8), (
        f"CAP curve dips below the chance line.\n"
        f"x: {cumulative_total}\n"
        f"y: {y_true_cumulative}"
    )

    # Check that CAP curve is on or below the perfect prediction line
    perfect_y = np.interp(
        display.cumulative_total, display._perfect_x, display._perfect_y
    )
    assert np.all(y_true_cumulative <= perfect_y + 1e-8), (
        "CAP curve rises above the perfect prediction line."
        f"x: {perfect_y}\n"
        f"y: {y_true_cumulative}"
    )


@pytest.mark.parametrize(
    "response_method", ["predict_proba", "decision_function", "auto"]
)
def test_with_linear_svc(pyplot, data_binary, response_method):
    X, y = data_binary
    svc = LinearSVC().fit(X, y)

    if response_method == "predict_proba":
        match = "LinearSVC has none of the following attributes: predict_proba."
        with pytest.raises(
            AttributeError,
            match=match,
        ):
            CAPCurveDisplay.from_estimator(svc, X, y, response_method=response_method)
    else:
        CAPCurveDisplay.from_estimator(svc, X, y, response_method=response_method)


def test_classifier_without_expected_methods(
    pyplot,
    data_binary,
    logistic_regression_model,
):
    class InvalidClassifier(ClassifierMixin, BaseEstimator):
        def __init__(self, *, param=1):
            self.param = param

        def fit(self, X, y=None):
            self.is_fitted_ = True
            return self

        def predict(self, X):
            return np.full(shape=X.shape[0], fill_value=self.param)

    X, y = data_binary
    invalid_classif = InvalidClassifier().fit(X, y)

    match = "Estimator does not have a predict_proba or decision_function method."
    with pytest.raises(AttributeError, match=match):
        CAPCurveDisplay.from_estimator(invalid_classif, X, y)
