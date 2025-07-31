import numpy as np
import pytest
from numpy.testing import assert_allclose

from sklearn.datasets import load_iris
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import DetCurveDisplay, det_curve


@pytest.mark.parametrize("constructor_name", ["from_estimator", "from_predictions"])
@pytest.mark.parametrize("response_method", ["predict_proba", "decision_function"])
@pytest.mark.parametrize("with_sample_weight", [True, False])
@pytest.mark.parametrize("drop_intermediate", [True, False])
@pytest.mark.parametrize("with_strings", [True, False])
def test_det_curve_display(
    pyplot,
    constructor_name,
    response_method,
    with_sample_weight,
    drop_intermediate,
    with_strings,
):
    X, y = load_iris(return_X_y=True)
    # Binarize the data with only the two first classes
    X, y = X[y < 2], y[y < 2]

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
        "alpha": 0.8,
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

    assert disp.estimator_name == "LogisticRegression"

    # cannot fail thanks to pyplot fixture
    import matplotlib as mpl

    assert isinstance(disp.line_, mpl.lines.Line2D)
    assert disp.line_.get_alpha() == 0.8
    assert isinstance(disp.ax_, mpl.axes.Axes)
    assert isinstance(disp.figure_, mpl.figure.Figure)
    assert disp.line_.get_label() == "LogisticRegression"

    expected_pos_label = 1 if pos_label is None else pos_label
    expected_ylabel = f"False Negative Rate (Positive label: {expected_pos_label})"
    expected_xlabel = f"False Positive Rate (Positive label: {expected_pos_label})"
    assert disp.ax_.get_ylabel() == expected_ylabel
    assert disp.ax_.get_xlabel() == expected_xlabel


@pytest.mark.parametrize(
    "constructor_name, expected_clf_name",
    [
        ("from_estimator", "LogisticRegression"),
        ("from_predictions", "Classifier"),
    ],
)
def test_det_curve_display_default_name(
    pyplot,
    constructor_name,
    expected_clf_name,
):
    # Check the default name display in the figure when `name` is not provided
    X, y = load_iris(return_X_y=True)
    # Binarize the data with only the two first classes
    X, y = X[y < 2], y[y < 2]

    lr = LogisticRegression().fit(X, y)
    y_score = lr.predict_proba(X)[:, 1]

    if constructor_name == "from_estimator":
        disp = DetCurveDisplay.from_estimator(lr, X, y)
    else:
        disp = DetCurveDisplay.from_predictions(y, y_score)

    assert disp.estimator_name == expected_clf_name
    assert disp.line_.get_label() == expected_clf_name


# TODO(1.10): remove
def test_y_score_and_y_pred_specified_error(pyplot):
    """1. Check that an error is raised when both y_score and y_pred are specified.
    2. Check that a warning is raised when y_pred is specified.
    """
    y_true = np.array([0, 0, 1, 1])
    y_score = np.array([0.1, 0.4, 0.35, 0.8])
    y_pred = np.array([0.2, 0.3, 0.5, 0.1])

    with pytest.raises(
        ValueError, match="`y_pred` and `y_score` cannot be both specified"
    ):
        DetCurveDisplay.from_predictions(y_true, y_score=y_score, y_pred=y_pred)

    with pytest.warns(FutureWarning, match="y_pred was deprecated in 1.8"):
        DetCurveDisplay.from_predictions(y_true, y_pred=y_score)
