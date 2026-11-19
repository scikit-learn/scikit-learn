import numpy as np
import pytest
from numpy.testing import assert_allclose

from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import DetCurveDisplay, det_curve
from sklearn.model_selection import cross_validate


@pytest.fixture(scope="module")
def data_binary():
    X, y = make_classification(
        n_samples=200,
        n_features=20,
        n_informative=5,
        n_redundant=2,
        flip_y=0.1,
        class_sep=0.8,
        random_state=42,
    )
    return X, y


def _check_figure_axes_and_labels(display, pos_label):
    """Check mpl axes and figure labels are correct."""
    import matplotlib as mpl

    assert isinstance(display.ax_, mpl.axes.Axes)
    assert isinstance(display.figure_, mpl.figure.Figure)

    expected_pos_label = 1 if pos_label is None else pos_label
    expected_ylabel = f"False Negative Rate (Positive label: {expected_pos_label})"
    expected_xlabel = f"False Positive Rate (Positive label: {expected_pos_label})"
    assert display.ax_.get_ylabel() == expected_ylabel
    assert display.ax_.get_xlabel() == expected_xlabel


@pytest.mark.parametrize("constructor_name", ["from_estimator", "from_predictions"])
@pytest.mark.parametrize("response_method", ["predict_proba", "decision_function"])
@pytest.mark.parametrize("with_sample_weight", [True, False])
@pytest.mark.parametrize("drop_intermediate", [True, False])
@pytest.mark.parametrize("with_strings", [True, False])
def test_det_curve_display(
    pyplot,
    data_binary,
    constructor_name,
    response_method,
    with_sample_weight,
    drop_intermediate,
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
    y_score = getattr(lr, response_method)(X)
    if y_score.ndim == 2:
        y_score = y_score[:, 1]
    # safe guard for the binary if/else construction
    assert constructor_name in ("from_estimator", "from_predictions")

    common_kwargs = {
        "name": lr.__class__.__name__,
        "curve_kwargs": {"alpha": 0.8},
        "sample_weight": sample_weight,
        "drop_intermediate": drop_intermediate,
        "pos_label": pos_label,
    }
    if constructor_name == "from_estimator":
        disp = DetCurveDisplay.from_estimator(lr, X, y, **common_kwargs)
    else:
        disp = DetCurveDisplay.from_predictions(y, y_score, **common_kwargs)

    fpr, fnr, _ = det_curve(
        y,
        y_score,
        sample_weight=sample_weight,
        drop_intermediate=drop_intermediate,
        pos_label=pos_label,
    )

    assert_allclose(disp.fpr, fpr, atol=1e-7)
    assert_allclose(disp.fnr, fnr, atol=1e-7)

    assert disp.name == "LogisticRegression"

    import matplotlib as mpl

    _check_figure_axes_and_labels(disp, pos_label)
    assert isinstance(disp.line_, mpl.lines.Line2D)
    assert disp.line_.get_alpha() == 0.8
    assert disp.line_.get_label() == "LogisticRegression"


def test_det_curve_display_from_cv_results(pyplot, data_binary):
    """Check one `from_cv_results` call forwards its arguments to `det_curve`."""
    X, y = data_binary
    y = np.array(["c", "b"])[y]
    sample_weight = np.ones(X.shape[0])
    sample_weight[::2] = 2

    cv_results = cross_validate(
        LogisticRegression(), X, y, cv=3, return_estimator=True, return_indices=True
    )
    display = DetCurveDisplay.from_cv_results(
        cv_results,
        X,
        y,
        sample_weight=sample_weight,
        drop_intermediate=False,
        response_method="decision_function",
        pos_label="c",
        name="test",
    )

    for fpr, fnr, estimator, test_indices in zip(
        display.fpr, display.fnr, cv_results["estimator"], cv_results["indices"]["test"]
    ):
        y_score = estimator.decision_function(X[test_indices])
        expected_fpr, expected_fnr, _ = det_curve(
            y[test_indices],
            y_score,
            sample_weight=sample_weight[test_indices],
            drop_intermediate=False,
            pos_label="c",
        )
        assert_allclose(fpr, expected_fpr)
        assert_allclose(fnr, expected_fnr)

    assert display.pos_label == "c"
    assert [line.get_label() for line in display.line_] == [
        "test",
        "_child1",
        "_child2",
    ]
    assert display.line_[0].get_linestyle() == "--"


def test_det_curve_plot_parameter_length_validation(pyplot):
    """Check `plot` rejects `fpr`, `fnr`, and `name` lists of different lengths."""
    fpr = [np.array([0.5, 0.5, 0.0])] * 2
    with pytest.raises(ValueError, match="self.fpr and self.fnr"):
        DetCurveDisplay(fpr=fpr, fnr=[np.array([0.0, 1.0, 1.0])]).plot()

    fnr = [np.array([0.0, 1.0, 1.0])] * 2
    with pytest.raises(ValueError, match="'name'"):
        DetCurveDisplay(fpr=fpr, fnr=fnr, name=["a", "b", "c"]).plot()
