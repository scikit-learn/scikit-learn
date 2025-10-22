import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.integrate import trapezoid

from sklearn import clone
from sklearn.compose import make_column_transformer
from sklearn.datasets import make_classification
from sklearn.exceptions import NotFittedError
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import RocCurveDisplay, auc, roc_curve
from sklearn.metrics._plot.tests.test_common_curve_display import (
    _check_pos_label_statistics,
)
from sklearn.model_selection import cross_validate
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.utils import _safe_indexing
from sklearn.utils._response import _get_response_values_binary


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
    """Check mpl axes and figure defaults are correct."""
    import matplotlib as mpl

    assert isinstance(display.ax_, mpl.axes.Axes)
    assert isinstance(display.figure_, mpl.figure.Figure)
    assert display.ax_.get_adjustable() == "box"
    assert display.ax_.get_aspect() in ("equal", 1.0)
    assert display.ax_.get_xlim() == display.ax_.get_ylim() == (-0.01, 1.01)

    expected_pos_label = 1 if pos_label is None else pos_label
    expected_ylabel = f"True Positive Rate (Positive label: {expected_pos_label})"
    expected_xlabel = f"False Positive Rate (Positive label: {expected_pos_label})"

    assert display.ax_.get_ylabel() == expected_ylabel
    assert display.ax_.get_xlabel() == expected_xlabel


@pytest.mark.parametrize("response_method", ["predict_proba", "decision_function"])
@pytest.mark.parametrize("with_sample_weight", [True, False])
@pytest.mark.parametrize("drop_intermediate", [True, False])
@pytest.mark.parametrize("with_strings", [True, False])
@pytest.mark.parametrize(
    "constructor_name, default_name",
    [
        ("from_estimator", "LogisticRegression"),
        ("from_predictions", "Classifier"),
    ],
)
def test_roc_curve_display_plotting(
    pyplot,
    response_method,
    data_binary,
    with_sample_weight,
    drop_intermediate,
    with_strings,
    constructor_name,
    default_name,
):
    """Check the overall plotting behaviour for single curve."""
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
    y_score = y_score if y_score.ndim == 1 else y_score[:, 1]

    if constructor_name == "from_estimator":
        display = RocCurveDisplay.from_estimator(
            lr,
            X,
            y,
            sample_weight=sample_weight,
            drop_intermediate=drop_intermediate,
            pos_label=pos_label,
            curve_kwargs={"alpha": 0.8},
        )
    else:
        display = RocCurveDisplay.from_predictions(
            y,
            y_score,
            sample_weight=sample_weight,
            drop_intermediate=drop_intermediate,
            pos_label=pos_label,
            curve_kwargs={"alpha": 0.8},
        )

    fpr, tpr, _ = roc_curve(
        y,
        y_score,
        sample_weight=sample_weight,
        drop_intermediate=drop_intermediate,
        pos_label=pos_label,
    )

    assert_allclose(display.roc_auc, auc(fpr, tpr))
    assert_allclose(display.fpr, fpr)
    assert_allclose(display.tpr, tpr)

    assert display.name == default_name

    import matplotlib as mpl

    _check_figure_axes_and_labels(display, pos_label)
    assert isinstance(display.line_, mpl.lines.Line2D)
    assert display.line_.get_alpha() == 0.8

    expected_label = f"{default_name} (AUC = {display.roc_auc:.2f})"
    assert display.line_.get_label() == expected_label


@pytest.mark.parametrize(
    "params, err_msg",
    [
        (
            {
                "fpr": [np.array([0, 0.5, 1]), np.array([0, 0.5, 1])],
                "tpr": [np.array([0, 0.5, 1])],
                "roc_auc": None,
                "name": None,
            },
            "self.fpr and self.tpr from `RocCurveDisplay` initialization,",
        ),
        (
            {
                "fpr": [np.array([0, 0.5, 1])],
                "tpr": [np.array([0, 0.5, 1]), np.array([0, 0.5, 1])],
                "roc_auc": [0.8, 0.9],
                "name": None,
            },
            "self.fpr, self.tpr and self.roc_auc from `RocCurveDisplay`",
        ),
        (
            {
                "fpr": [np.array([0, 0.5, 1]), np.array([0, 0.5, 1])],
                "tpr": [np.array([0, 0.5, 1]), np.array([0, 0.5, 1])],
                "roc_auc": [0.8],
                "name": None,
            },
            "Got: self.fpr: 2, self.tpr: 2, self.roc_auc: 1",
        ),
        (
            {
                "fpr": [np.array([0, 0.5, 1]), np.array([0, 0.5, 1])],
                "tpr": [np.array([0, 0.5, 1]), np.array([0, 0.5, 1])],
                "roc_auc": [0.8, 0.9],
                "name": ["curve1", "curve2", "curve3"],
            },
            r"self.fpr, self.tpr, self.roc_auc and 'name' \(or self.name\)",
        ),
        (
            {
                "fpr": [np.array([0, 0.5, 1]), np.array([0, 0.5, 1])],
                "tpr": [np.array([0, 0.5, 1]), np.array([0, 0.5, 1])],
                "roc_auc": [0.8, 0.9],
                # List of length 1 is always allowed
                "name": ["curve1"],
            },
            None,
        ),
    ],
)
def test_roc_curve_plot_parameter_length_validation(pyplot, params, err_msg):
    """Check `plot` parameter length validation performed correctly."""
    display = RocCurveDisplay(**params)
    if err_msg:
        with pytest.raises(ValueError, match=err_msg):
            display.plot()
    else:
        # No error should be raised
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
def test_roc_curve_display_plotting_from_cv_results(
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
    display = RocCurveDisplay.from_cv_results(
        cv_results,
        X,
        y,
        sample_weight=sample_weight,
        drop_intermediate=drop_intermediate,
        response_method=response_method,
        pos_label=pos_label,
        curve_kwargs=curve_kwargs,
    )

    for idx, (estimator, test_indices) in enumerate(
        zip(cv_results["estimator"], cv_results["indices"]["test"])
    ):
        y_true = _safe_indexing(y, test_indices)
        y_pred = _get_response_values_binary(
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
        fpr, tpr, _ = roc_curve(
            y_true,
            y_pred,
            sample_weight=sample_weight_fold,
            drop_intermediate=drop_intermediate,
            pos_label=pos_label,
        )
        assert_allclose(display.roc_auc[idx], auc(fpr, tpr))
        assert_allclose(display.fpr[idx], fpr)
        assert_allclose(display.tpr[idx], tpr)

    assert display.name is None

    import matplotlib as mpl

    _check_figure_axes_and_labels(display, pos_label)
    if with_sample_weight:
        aggregate_expected_labels = ["AUC = 0.64 +/- 0.04", "_child1", "_child2"]
    else:
        aggregate_expected_labels = ["AUC = 0.61 +/- 0.05", "_child1", "_child2"]
    for idx, line in enumerate(display.line_):
        assert isinstance(line, mpl.lines.Line2D)
        # Default alpha for `from_cv_results`
        line.get_alpha() == 0.5
        if isinstance(curve_kwargs, list):
            # Each individual curve labelled
            assert line.get_label() == f"AUC = {display.roc_auc[idx]:.2f}"
        else:
            # Single aggregate label
            assert line.get_label() == aggregate_expected_labels[idx]


def _check_chance_level(plot_chance_level, chance_level_kw, display):
    """Check chance level line and line styles correct."""
    import matplotlib as mpl

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
        assert display.chance_level_.get_label() == "Chance level (AUC = 0.5)"
    elif plot_chance_level:
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


@pytest.mark.parametrize("plot_chance_level", [True, False])
@pytest.mark.parametrize("label", [None, "Test Label"])
@pytest.mark.parametrize(
    "chance_level_kw",
    [
        None,
        {"linewidth": 1, "color": "red", "linestyle": "-", "label": "DummyEstimator"},
        {"lw": 1, "c": "red", "ls": "-", "label": "DummyEstimator"},
        {"lw": 1, "color": "blue", "ls": "-", "label": None},
    ],
)
@pytest.mark.parametrize("constructor_name", ["from_estimator", "from_predictions"])
def test_roc_curve_chance_level_line(
    pyplot,
    data_binary,
    plot_chance_level,
    chance_level_kw,
    label,
    constructor_name,
):
    """Check chance level plotting behavior of `from_predictions`, `from_estimator`."""
    X, y = data_binary

    lr = LogisticRegression()
    lr.fit(X, y)

    y_score = getattr(lr, "predict_proba")(X)
    y_score = y_score if y_score.ndim == 1 else y_score[:, 1]

    if constructor_name == "from_estimator":
        display = RocCurveDisplay.from_estimator(
            lr,
            X,
            y,
            curve_kwargs={"alpha": 0.8, "label": label},
            plot_chance_level=plot_chance_level,
            chance_level_kw=chance_level_kw,
        )
    else:
        display = RocCurveDisplay.from_predictions(
            y,
            y_score,
            curve_kwargs={"alpha": 0.8, "label": label},
            plot_chance_level=plot_chance_level,
            chance_level_kw=chance_level_kw,
        )

    import matplotlib as mpl

    assert isinstance(display.line_, mpl.lines.Line2D)
    assert display.line_.get_alpha() == 0.8
    assert isinstance(display.ax_, mpl.axes.Axes)
    assert isinstance(display.figure_, mpl.figure.Figure)

    _check_chance_level(plot_chance_level, chance_level_kw, display)

    # Checking for legend behaviour
    if plot_chance_level and chance_level_kw is not None:
        if label is not None or chance_level_kw.get("label") is not None:
            legend = display.ax_.get_legend()
            assert legend is not None  #  Legend should be present if any label is set
            legend_labels = [text.get_text() for text in legend.get_texts()]
            if label is not None:
                assert label in legend_labels
            if chance_level_kw.get("label") is not None:
                assert chance_level_kw["label"] in legend_labels
        else:
            assert display.ax_.get_legend() is None


@pytest.mark.parametrize("plot_chance_level", [True, False])
@pytest.mark.parametrize(
    "chance_level_kw",
    [
        None,
        {"linewidth": 1, "color": "red", "linestyle": "-", "label": "DummyEstimator"},
        {"lw": 1, "c": "red", "ls": "-", "label": "DummyEstimator"},
        {"lw": 1, "color": "blue", "ls": "-", "label": None},
    ],
)
@pytest.mark.parametrize("curve_kwargs", [None, {"alpha": 0.8}])
def test_roc_curve_chance_level_line_from_cv_results(
    pyplot,
    data_binary,
    plot_chance_level,
    chance_level_kw,
    curve_kwargs,
):
    """Check chance level plotting behavior with `from_cv_results`."""
    X, y = data_binary
    n_cv = 3
    cv_results = cross_validate(
        LogisticRegression(), X, y, cv=n_cv, return_estimator=True, return_indices=True
    )

    display = RocCurveDisplay.from_cv_results(
        cv_results,
        X,
        y,
        plot_chance_level=plot_chance_level,
        chance_level_kwargs=chance_level_kw,
        curve_kwargs=curve_kwargs,
    )

    import matplotlib as mpl

    assert all(isinstance(line, mpl.lines.Line2D) for line in display.line_)
    # Ensure both curve line kwargs passed correctly as well
    if curve_kwargs:
        assert all(line.get_alpha() == 0.8 for line in display.line_)
    assert isinstance(display.ax_, mpl.axes.Axes)
    assert isinstance(display.figure_, mpl.figure.Figure)

    _check_chance_level(plot_chance_level, chance_level_kw, display)

    legend = display.ax_.get_legend()
    # There is always a legend, to indicate each 'Fold' curve
    assert legend is not None
    legend_labels = [text.get_text() for text in legend.get_texts()]
    if plot_chance_level and chance_level_kw is not None:
        if chance_level_kw.get("label") is not None:
            assert chance_level_kw["label"] in legend_labels
        else:
            assert len(legend_labels) == 1


@pytest.mark.parametrize(
    "clf",
    [
        LogisticRegression(),
        make_pipeline(StandardScaler(), LogisticRegression()),
        make_pipeline(
            make_column_transformer((StandardScaler(), [0, 1])), LogisticRegression()
        ),
    ],
)
@pytest.mark.parametrize("constructor_name", ["from_estimator", "from_predictions"])
def test_roc_curve_display_complex_pipeline(pyplot, data_binary, clf, constructor_name):
    """Check the behaviour with complex pipeline."""
    X, y = data_binary

    clf = clone(clf)

    if constructor_name == "from_estimator":
        with pytest.raises(NotFittedError):
            RocCurveDisplay.from_estimator(clf, X, y)

    clf.fit(X, y)

    if constructor_name == "from_estimator":
        display = RocCurveDisplay.from_estimator(clf, X, y)
        name = clf.__class__.__name__
    else:
        display = RocCurveDisplay.from_predictions(y, y)
        name = "Classifier"

    assert name in display.line_.get_label()
    assert display.name == name


@pytest.mark.parametrize(
    "roc_auc, name, curve_kwargs, expected_labels",
    [
        ([0.9, 0.8], None, None, ["AUC = 0.85 +/- 0.05", "_child1"]),
        ([0.9, 0.8], "Est name", None, ["Est name (AUC = 0.85 +/- 0.05)", "_child1"]),
        (
            [0.8, 0.7],
            ["fold1", "fold2"],
            [{"c": "blue"}, {"c": "red"}],
            ["fold1 (AUC = 0.80)", "fold2 (AUC = 0.70)"],
        ),
        (None, ["fold1", "fold2"], [{"c": "blue"}, {"c": "red"}], ["fold1", "fold2"]),
    ],
)
def test_roc_curve_display_default_labels(
    pyplot, roc_auc, name, curve_kwargs, expected_labels
):
    """Check the default labels used in the display."""
    fpr = [np.array([0, 0.5, 1]), np.array([0, 0.3, 1])]
    tpr = [np.array([0, 0.5, 1]), np.array([0, 0.3, 1])]
    disp = RocCurveDisplay(fpr=fpr, tpr=tpr, roc_auc=roc_auc, name=name).plot(
        curve_kwargs=curve_kwargs
    )
    for idx, expected_label in enumerate(expected_labels):
        assert disp.line_[idx].get_label() == expected_label


@pytest.mark.parametrize("response_method", ["predict_proba", "decision_function"])
@pytest.mark.parametrize(
    "constructor_name", ["from_estimator", "from_predictions", "from_cv_results"]
)
def test_plot_roc_curve_pos_label(pyplot, response_method, constructor_name):
    """Check switching `pos_label` give correct statistics, using imbalanced data."""

    def _check_auc(display, constructor_name, pos_label):
        roc_auc_limit = 0.95679
        roc_auc_limit_multi = [0.97007, 0.985915, 0.980952]

        if constructor_name == "from_cv_results":
            for idx, roc_auc in enumerate(display.roc_auc):
                assert roc_auc == pytest.approx(roc_auc_limit_multi[idx])
        else:
            assert display.roc_auc == pytest.approx(roc_auc_limit)
            assert trapezoid(display.tpr, display.fpr) == pytest.approx(roc_auc_limit)

    _check_pos_label_statistics(
        RocCurveDisplay, response_method, constructor_name, _check_auc
    )


# TODO(1.9): remove
def test_y_score_and_y_pred_specified_error(pyplot):
    """1. Check that an error is raised when both y_score and y_pred are specified.
    2. Check that a warning is raised when y_pred is specified.
    """
    y_true = np.array([0, 1, 1, 0])
    y_score = np.array([0.1, 0.4, 0.35, 0.8])
    y_pred = np.array([0.2, 0.3, 0.5, 0.1])

    with pytest.raises(
        ValueError, match="`y_pred` and `y_score` cannot be both specified"
    ):
        RocCurveDisplay.from_predictions(y_true, y_score=y_score, y_pred=y_pred)

    with pytest.warns(FutureWarning, match="y_pred was deprecated in 1.7"):
        display_y_pred = RocCurveDisplay.from_predictions(y_true, y_pred=y_score)
    desired_fpr, desired_fnr, _ = roc_curve(y_true, y_score)
    assert_allclose(display_y_pred.fpr, desired_fpr)
    assert_allclose(display_y_pred.tpr, desired_fnr)

    display_y_score = RocCurveDisplay.from_predictions(y_true, y_score)
    assert_allclose(display_y_score.fpr, desired_fpr)
    assert_allclose(display_y_score.tpr, desired_fnr)


@pytest.mark.parametrize("despine", [True, False])
@pytest.mark.parametrize(
    "constructor_name", ["from_estimator", "from_predictions", "from_cv_results"]
)
def test_plot_roc_curve_despine(pyplot, data_binary, despine, constructor_name):
    # Check that the despine keyword is working correctly
    X, y = data_binary

    lr = LogisticRegression().fit(X, y)
    lr.fit(X, y)
    cv_results = cross_validate(
        LogisticRegression(), X, y, cv=3, return_estimator=True, return_indices=True
    )

    y_pred = lr.decision_function(X)

    # safe guard for the if/else construction
    assert constructor_name in ("from_estimator", "from_predictions", "from_cv_results")

    if constructor_name == "from_estimator":
        display = RocCurveDisplay.from_estimator(lr, X, y, despine=despine)
    elif constructor_name == "from_predictions":
        display = RocCurveDisplay.from_predictions(y, y_pred, despine=despine)
    else:
        display = RocCurveDisplay.from_cv_results(cv_results, X, y, despine=despine)

    for s in ["top", "right"]:
        assert display.ax_.spines[s].get_visible() is not despine

    if despine:
        for s in ["bottom", "left"]:
            assert display.ax_.spines[s].get_bounds() == (0, 1)
