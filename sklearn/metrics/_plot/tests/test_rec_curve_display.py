import numpy as np
import pytest
from numpy.testing import assert_allclose

from sklearn import clone
from sklearn.compose import make_column_transformer
from sklearn.datasets import make_regression
from sklearn.exceptions import NotFittedError
from sklearn.linear_model import LinearRegression
from sklearn.metrics import RecCurveDisplay, rec_curve
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


@pytest.fixture(scope="module")
def regression_data():
    X, y = make_regression(
        n_samples=200,
        n_features=20,
        n_informative=5,
        noise=1,
        random_state=42,
    )
    return X, y


@pytest.mark.parametrize("plot_const_predictor", [True, False])
@pytest.mark.parametrize("clip_max_const_error", [True, False])
@pytest.mark.parametrize("constant_predictor", [None, "mean", "median"])
@pytest.mark.parametrize("loss", [None, "squared", "absolute"])
@pytest.mark.parametrize(
    "constructor_name, default_name",
    [
        ("from_estimator", "LinearRegression"),
        ("from_predictions", "Model"),
    ],
)
def test_roc_curve_display_plotting(
    pyplot,
    regression_data,
    plot_const_predictor,
    clip_max_const_error,
    constant_predictor,
    loss,
    constructor_name,
    default_name,
):
    X, y = regression_data

    lr = LinearRegression()
    lr.fit(X, y)
    y_score = lr.predict(X)

    ctor_kwargs = dict(
        plot_const_predictor=plot_const_predictor,
        clip_max_const_error=clip_max_const_error,
    )
    if loss is not None:
        ctor_kwargs |= dict(loss=loss)
    if constant_predictor is not None:
        ctor_kwargs |= dict(constant_predictor=constant_predictor)

    if constructor_name == "from_estimator":
        display = RecCurveDisplay.from_estimator(lr, X, y, **ctor_kwargs)
    else:
        display = RecCurveDisplay.from_predictions(y, y_score, **ctor_kwargs)

    if loss is None:
        deviations, accuracy = rec_curve(y, y_score)
    else:
        deviations, accuracy = rec_curve(y, y_score, loss=loss)

    assert_allclose(display.deviations, deviations)
    assert_allclose(display.accuracy, accuracy)

    if plot_const_predictor or clip_max_const_error:
        const_predictor_labels = {
            "mean": "Mean Predictor",
            "median": "Median Predictor",
        }
        const_predictor_key = constant_predictor or (
            "mean" if loss == "squared" else "median"
        )
        assert (
            display.constant_predictor_name
            == const_predictor_labels[const_predictor_key]
        )

        const_value = np.mean(y) if const_predictor_key == "mean" else np.median(y)
        if loss is None:
            exp_const_deviations, exp_const_accuracy = rec_curve(y, const_value)
        else:
            exp_const_deviations, exp_const_accuracy = rec_curve(
                y, const_value, loss=loss
            )
        assert_allclose(display.constant_predictor_deviations, exp_const_deviations)
        assert_allclose(display.constant_predictor_accuracy, exp_const_accuracy)

    assert display.estimator_name == default_name
    assert display.loss == (loss or "absolute")

    import matplotlib as mpl

    assert isinstance(display.line_, mpl.lines.Line2D)
    assert isinstance(display.ax_, mpl.axes.Axes)
    assert isinstance(display.figure_, mpl.figure.Figure)

    if plot_const_predictor:
        assert isinstance(display.constant_predictor_line_, mpl.lines.Line2D)

    expected_label = f"{default_name}"
    assert display.line_.get_label() == expected_label

    expected_ylabel = "Accuracy (Fraction of samples)"
    expected_xlabel = f"Error Tolerance (Deviation - {loss or 'absolute'} loss)"

    assert display.ax_.get_ylabel() == expected_ylabel
    assert display.ax_.get_xlabel() == expected_xlabel


@pytest.mark.parametrize(
    "clf",
    [
        LinearRegression(),
        make_pipeline(StandardScaler(), LinearRegression()),
        make_pipeline(
            make_column_transformer((StandardScaler(), [0, 1])), LinearRegression()
        ),
    ],
)
@pytest.mark.parametrize("constructor_name", ["from_estimator", "from_predictions"])
def test_rec_curve_display_complex_pipeline(
    pyplot, regression_data, clf, constructor_name
):
    """Check the behaviour with complex pipeline."""
    X, y = regression_data

    clf = clone(clf)

    if constructor_name == "from_estimator":
        with pytest.raises(NotFittedError):
            RecCurveDisplay.from_estimator(clf, X, y)

    clf.fit(X, y)

    if constructor_name == "from_estimator":
        display = RecCurveDisplay.from_estimator(clf, X, y)
        name = clf.__class__.__name__
    else:
        display = RecCurveDisplay.from_predictions(y, y)
        name = "Model"

    assert name in display.line_.get_label()
    assert display.estimator_name == name
