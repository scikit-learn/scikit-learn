import numpy as np
import pytest
from numpy.testing import assert_allclose

from sklearn.datasets import load_iris
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import DetCurveDisplay, det_curve
from sklearn.model_selection import cross_validate
from sklearn.utils import _safe_indexing
from sklearn.utils._response import _get_response_values_binary


@pytest.fixture(scope="module")
def data_binary():
    X, y = load_iris(return_X_y=True)
    # Binarize the data with only the two first classes
    return X[y < 2], y[y < 2]


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


@pytest.mark.parametrize(
    "params, err_msg",
    [
        (
            {
                "fpr": [np.array([0.5, 0.5, 0]), np.array([0.6, 0.5, 0])],
                "fnr": [np.array([0, 1, 1])],
                "name": None,
            },
            "self.fpr and self.fnr from `DetCurveDisplay` initialization,",
        ),
        (
            {
                "fpr": [np.array([0.5, 0.5, 0])],
                "fnr": [np.array([0, 1, 1]), np.array([0, 1, 1])],
                "name": None,
            },
            "self.fpr and self.fnr from `DetCurveDisplay`",
        ),
        (
            {
                "fpr": [np.array([0.5, 0.5, 0]), np.array([0.5, 0.5, 0])],
                "fnr": [np.array([0, 1, 1]), np.array([0, 1, 1])],
                "name": ["curve1", "curve2", "curve3"],
            },
            (
                r"self\.fpr, self\.fnr and 'name' \(or self\.name\).*Got: self.fpr: 2, "
                r"self.fnr: 2, 'name' \(or self.name\): 3"
            ),
        ),
        (
            {
                "fpr": [np.array([0.5, 0.5, 0]), np.array([0.5, 0.5, 0])],
                "fnr": [np.array([0, 1, 1]), np.array([0, 1, 1])],
                # List of length 1 is always allowed
                "name": ["curve1"],
            },
            None,
        ),
    ],
)
def test_det_curve_plot_parameter_length_validation(pyplot, params, err_msg):
    """Check `plot` parameter length validation performed correctly."""
    display = DetCurveDisplay(**params)
    if err_msg:
        with pytest.raises(ValueError, match=err_msg):
            display.plot()
    else:
        display.plot()


@pytest.mark.parametrize(
    "curve_kwargs",
    [
        None,
        {"color": "blue"},
        [{"color": "blue"}, {"color": "green"}, {"color": "red"}],
    ],
)
@pytest.mark.parametrize("drop_intermediate", [True, False])
@pytest.mark.parametrize("response_method", ["predict_proba", "decision_function"])
@pytest.mark.parametrize("with_sample_weight", [True, False])
@pytest.mark.parametrize("with_strings", [True, False])
def test_det_curve_display_plotting_from_cv_results(
    pyplot,
    data_binary,
    with_strings,
    with_sample_weight,
    response_method,
    drop_intermediate,
    curve_kwargs,
):
    """Check overall plotting of `from_cv_results`."""
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

    cv_results = cross_validate(
        LogisticRegression(), X, y, cv=3, return_estimator=True, return_indices=True
    )
    display = DetCurveDisplay.from_cv_results(
        cv_results,
        X,
        y,
        sample_weight=sample_weight,
        drop_intermediate=drop_intermediate,
        response_method=response_method,
        pos_label=pos_label,
        name="test",
        curve_kwargs=curve_kwargs,
    )

    for idx, (estimator, test_indices) in enumerate(
        zip(cv_results["estimator"], cv_results["indices"]["test"])
    ):
        y_true = _safe_indexing(y, test_indices)
        y_score = _get_response_values_binary(
            estimator,
            _safe_indexing(X, test_indices),
            response_method=response_method,
            pos_label=pos_label,
        )[0]
        sample_weight_fold = (
            None
            if sample_weight is None
            else _safe_indexing(sample_weight, test_indices)
        )
        fpr, fnr, _ = det_curve(
            y_true,
            y_score,
            sample_weight=sample_weight_fold,
            drop_intermediate=drop_intermediate,
            pos_label=pos_label,
        )
        assert_allclose(display.fpr[idx], fpr)
        assert_allclose(display.fnr[idx], fnr)

    assert display.name == "test"

    import matplotlib as mpl

    _check_figure_axes_and_labels(display, pos_label)
    aggregate_expected_labels = ["test", "_child1", "_child2"]
    for idx, line in enumerate(display.line_):
        assert isinstance(line, mpl.lines.Line2D)
        assert line.get_alpha() == 0.5
        assert line.get_linestyle() == "--"
        if isinstance(curve_kwargs, list):
            assert line.get_label() == "test"
        else:
            assert line.get_label() == aggregate_expected_labels[idx]


@pytest.mark.parametrize(
    "curve_kwargs",
    [None, {"color": "red"}, [{"c": "red"}, {"c": "green"}, {"c": "yellow"}]],
)
@pytest.mark.parametrize("name", [None, "single", ["one", "two", "three"]])
def test_det_curve_plot_legend_label(pyplot, name, curve_kwargs):
    """Check legend labels for every `curve_kwargs` and `name` combination.

    DET curves have no summary metric, so labels contain only `name`.
    """
    fpr = [np.array([0.5, 0.5, 0.0])] * 3
    fnr = [np.array([0.0, 1.0, 1.0])] * 3
    if not isinstance(curve_kwargs, list) and isinstance(name, list):
        with pytest.raises(ValueError, match="To avoid labeling individual curves"):
            DetCurveDisplay(fpr=fpr, fnr=fnr).plot(name=name, curve_kwargs=curve_kwargs)
        return

    display = DetCurveDisplay(fpr=fpr, fnr=fnr).plot(
        name=name, curve_kwargs=curve_kwargs
    )
    legend = display.ax_.get_legend()
    if legend is None:
        assert name is None
        return

    legend_labels = [text.get_text() for text in legend.get_texts()]
    if isinstance(curve_kwargs, list):
        assert len(legend_labels) == 3
        for idx, label in enumerate(legend_labels):
            if isinstance(name, str):
                assert label == "single"
            else:
                assert label == name[idx]
    else:
        assert len(legend_labels) == 1
        assert legend_labels[0] == "single"
