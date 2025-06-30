from collections.abc import Mapping

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.integrate import trapezoid

from sklearn import clone
from sklearn.compose import make_column_transformer
from sklearn.datasets import load_breast_cancer, make_classification
from sklearn.exceptions import NotFittedError
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import RocCurveDisplay, auc, roc_curve
from sklearn.model_selection import cross_validate, train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.utils import _safe_indexing, shuffle
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


def test_validate_plot_params(pyplot):
    """Check `_validate_plot_params` returns the correct variables."""
    fpr = np.array([0, 0.5, 1])
    tpr = [np.array([0, 0.5, 1])]
    roc_auc = None
    name = "test_curve"

    # Initialize display with test inputs
    display = RocCurveDisplay(
        fpr=fpr,
        tpr=tpr,
        roc_auc=roc_auc,
        name=name,
        pos_label=None,
    )
    fpr_out, tpr_out, roc_auc_out, name_out = display._validate_plot_params(
        ax=None, name=None
    )

    assert isinstance(fpr_out, list)
    assert isinstance(tpr_out, list)
    assert len(fpr_out) == 1
    assert len(tpr_out) == 1
    assert roc_auc_out is None
    assert name_out == ["test_curve"]


def test_roc_curve_from_cv_results_param_validation(pyplot, data_binary):
    """Check parameter validation is correct."""
    X, y = data_binary

    # `cv_results` missing key
    cv_results_no_est = cross_validate(
        LogisticRegression(), X, y, cv=3, return_estimator=True, return_indices=False
    )
    cv_results_no_indices = cross_validate(
        LogisticRegression(), X, y, cv=3, return_estimator=True, return_indices=False
    )
    for cv_results in (cv_results_no_est, cv_results_no_indices):
        with pytest.raises(
            ValueError,
            match="`cv_results` does not contain one of the following required",
        ):
            RocCurveDisplay.from_cv_results(cv_results, X, y)

    cv_results = cross_validate(
        LogisticRegression(), X, y, cv=3, return_estimator=True, return_indices=True
    )

    # `X` wrong length
    with pytest.raises(ValueError, match="`X` does not contain the correct"):
        RocCurveDisplay.from_cv_results(cv_results, X[:10, :], y)

    # `y` not binary
    y_multi = y.copy()
    y_multi[0] = 2
    with pytest.raises(ValueError, match="The target `y` is not binary."):
        RocCurveDisplay.from_cv_results(cv_results, X, y_multi)

    # input inconsistent length
    with pytest.raises(ValueError, match="Found input variables with inconsistent"):
        RocCurveDisplay.from_cv_results(cv_results, X, y[:10])
    with pytest.raises(ValueError, match="Found input variables with inconsistent"):
        RocCurveDisplay.from_cv_results(cv_results, X, y, sample_weight=[1, 2])

    # `pos_label` inconsistency
    y_multi[y_multi == 1] = 2
    with pytest.raises(ValueError, match=r"y takes value in \{0, 2\}"):
        RocCurveDisplay.from_cv_results(cv_results, X, y_multi)

    # `name` is list while `curve_kwargs` is None or dict
    for curve_kwargs in (None, {"alpha": 0.2}):
        with pytest.raises(ValueError, match="To avoid labeling individual curves"):
            RocCurveDisplay.from_cv_results(
                cv_results,
                X,
                y,
                name=["one", "two", "three"],
                curve_kwargs=curve_kwargs,
            )

    # `curve_kwargs` incorrect length
    with pytest.raises(ValueError, match="`curve_kwargs` must be None, a dictionary"):
        RocCurveDisplay.from_cv_results(cv_results, X, y, curve_kwargs=[{"alpha": 1}])

    # `curve_kwargs` both alias provided
    with pytest.raises(TypeError, match="Got both c and"):
        RocCurveDisplay.from_cv_results(
            cv_results, X, y, curve_kwargs={"c": "blue", "color": "red"}
        )


@pytest.mark.parametrize(
    "curve_kwargs",
    [None, {"alpha": 0.2}, [{"alpha": 0.2}, {"alpha": 0.3}, {"alpha": 0.4}]],
)
def test_roc_curve_display_from_cv_results_curve_kwargs(
    pyplot, data_binary, curve_kwargs
):
    """Check `curve_kwargs` correctly passed."""
    X, y = data_binary
    n_cv = 3
    cv_results = cross_validate(
        LogisticRegression(), X, y, cv=n_cv, return_estimator=True, return_indices=True
    )
    display = RocCurveDisplay.from_cv_results(
        cv_results,
        X,
        y,
        curve_kwargs=curve_kwargs,
    )
    if curve_kwargs is None:
        # Default `alpha` used
        assert all(line.get_alpha() == 0.5 for line in display.line_)
    elif isinstance(curve_kwargs, Mapping):
        # `alpha` from dict used for all curves
        assert all(line.get_alpha() == 0.2 for line in display.line_)
    else:
        # Different `alpha` used for each curve
        assert all(
            line.get_alpha() == curve_kwargs[i]["alpha"]
            for i, line in enumerate(display.line_)
        )


# TODO(1.9): Remove in 1.9
def test_roc_curve_display_estimator_name_deprecation(pyplot):
    """Check deprecation of `estimator_name`."""
    fpr = np.array([0, 0.5, 1])
    tpr = np.array([0, 0.5, 1])
    with pytest.warns(FutureWarning, match="`estimator_name` is deprecated in"):
        RocCurveDisplay(fpr=fpr, tpr=tpr, estimator_name="test")


# TODO(1.9): Remove in 1.9
@pytest.mark.parametrize(
    "constructor_name", ["from_estimator", "from_predictions", "plot"]
)
def test_roc_curve_display_kwargs_deprecation(pyplot, data_binary, constructor_name):
    """Check **kwargs deprecated correctly in favour of `curve_kwargs`."""
    X, y = data_binary
    lr = LogisticRegression()
    lr.fit(X, y)
    fpr = np.array([0, 0.5, 1])
    tpr = np.array([0, 0.5, 1])

    # Error when both `curve_kwargs` and `**kwargs` provided
    with pytest.raises(ValueError, match="Cannot provide both `curve_kwargs`"):
        if constructor_name == "from_estimator":
            RocCurveDisplay.from_estimator(
                lr, X, y, curve_kwargs={"alpha": 1}, label="test"
            )
        elif constructor_name == "from_predictions":
            RocCurveDisplay.from_predictions(
                y, y, curve_kwargs={"alpha": 1}, label="test"
            )
        else:
            RocCurveDisplay(fpr=fpr, tpr=tpr).plot(
                curve_kwargs={"alpha": 1}, label="test"
            )

    # Warning when `**kwargs`` provided
    with pytest.warns(FutureWarning, match=r"`\*\*kwargs` is deprecated and will be"):
        if constructor_name == "from_estimator":
            RocCurveDisplay.from_estimator(lr, X, y, label="test")
        elif constructor_name == "from_predictions":
            RocCurveDisplay.from_predictions(y, y, label="test")
        else:
            RocCurveDisplay(fpr=fpr, tpr=tpr).plot(label="test")


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


@pytest.mark.parametrize("roc_auc", [[1.0, 1.0, 1.0], None])
@pytest.mark.parametrize(
    "curve_kwargs",
    [None, {"color": "red"}, [{"c": "red"}, {"c": "green"}, {"c": "yellow"}]],
)
@pytest.mark.parametrize("name", [None, "single", ["one", "two", "three"]])
def test_roc_curve_plot_legend_label(pyplot, data_binary, name, curve_kwargs, roc_auc):
    """Check legend label correct with all `curve_kwargs`, `name` combinations."""
    fpr = [np.array([0, 0.5, 1]), np.array([0, 0.5, 1]), np.array([0, 0.5, 1])]
    tpr = [np.array([0, 0.5, 1]), np.array([0, 0.5, 1]), np.array([0, 0.5, 1])]
    if not isinstance(curve_kwargs, list) and isinstance(name, list):
        with pytest.raises(ValueError, match="To avoid labeling individual curves"):
            RocCurveDisplay(fpr=fpr, tpr=tpr, roc_auc=roc_auc).plot(
                name=name, curve_kwargs=curve_kwargs
            )

    else:
        display = RocCurveDisplay(fpr=fpr, tpr=tpr, roc_auc=roc_auc).plot(
            name=name, curve_kwargs=curve_kwargs
        )
        legend = display.ax_.get_legend()
        if legend is None:
            # No legend is created, exit test early
            assert name is None
            assert roc_auc is None
            return
        else:
            legend_labels = [text.get_text() for text in legend.get_texts()]

        if isinstance(curve_kwargs, list):
            # Multiple labels in legend
            assert len(legend_labels) == 3
            for idx, label in enumerate(legend_labels):
                if name is None:
                    expected_label = "AUC = 1.00" if roc_auc else None
                    assert label == expected_label
                elif isinstance(name, str):
                    expected_label = "single (AUC = 1.00)" if roc_auc else "single"
                    assert label == expected_label
                else:
                    # `name` is a list of different strings
                    expected_label = (
                        f"{name[idx]} (AUC = 1.00)" if roc_auc else f"{name[idx]}"
                    )
                    assert label == expected_label
        else:
            # Single label in legend
            assert len(legend_labels) == 1
            if name is None:
                expected_label = "AUC = 1.00 +/- 0.00" if roc_auc else None
                assert legend_labels[0] == expected_label
            else:
                # name is single string
                expected_label = "single (AUC = 1.00 +/- 0.00)" if roc_auc else "single"
                assert legend_labels[0] == expected_label


@pytest.mark.parametrize(
    "curve_kwargs",
    [None, {"color": "red"}, [{"c": "red"}, {"c": "green"}, {"c": "yellow"}]],
)
@pytest.mark.parametrize("name", [None, "single", ["one", "two", "three"]])
def test_roc_curve_from_cv_results_legend_label(
    pyplot, data_binary, name, curve_kwargs
):
    """Check legend label correct with all `curve_kwargs`, `name` combinations."""
    X, y = data_binary
    n_cv = 3
    cv_results = cross_validate(
        LogisticRegression(), X, y, cv=n_cv, return_estimator=True, return_indices=True
    )

    if not isinstance(curve_kwargs, list) and isinstance(name, list):
        with pytest.raises(ValueError, match="To avoid labeling individual curves"):
            RocCurveDisplay.from_cv_results(
                cv_results, X, y, name=name, curve_kwargs=curve_kwargs
            )
    else:
        display = RocCurveDisplay.from_cv_results(
            cv_results, X, y, name=name, curve_kwargs=curve_kwargs
        )

        legend = display.ax_.get_legend()
        legend_labels = [text.get_text() for text in legend.get_texts()]
        if isinstance(curve_kwargs, list):
            # Multiple labels in legend
            assert len(legend_labels) == 3
            auc = ["0.62", "0.66", "0.55"]
            for idx, label in enumerate(legend_labels):
                if name is None:
                    assert label == f"AUC = {auc[idx]}"
                elif isinstance(name, str):
                    assert label == f"single (AUC = {auc[idx]})"
                else:
                    # `name` is a list of different strings
                    assert label == f"{name[idx]} (AUC = {auc[idx]})"
        else:
            # Single label in legend
            assert len(legend_labels) == 1
            if name is None:
                assert legend_labels[0] == "AUC = 0.61 +/- 0.05"
            else:
                # name is single string
                assert legend_labels[0] == "single (AUC = 0.61 +/- 0.05)"


@pytest.mark.parametrize(
    "curve_kwargs",
    [None, {"color": "red"}, [{"c": "red"}, {"c": "green"}, {"c": "yellow"}]],
)
def test_roc_curve_from_cv_results_curve_kwargs(pyplot, data_binary, curve_kwargs):
    """Check line kwargs passed correctly in `from_cv_results`."""

    X, y = data_binary
    cv_results = cross_validate(
        LogisticRegression(), X, y, cv=3, return_estimator=True, return_indices=True
    )
    display = RocCurveDisplay.from_cv_results(
        cv_results, X, y, curve_kwargs=curve_kwargs
    )

    for idx, line in enumerate(display.line_):
        color = line.get_color()
        if curve_kwargs is None:
            # Default color
            assert color == "blue"
        elif isinstance(curve_kwargs, Mapping):
            # All curves "red"
            assert color == "red"
        else:
            assert color == curve_kwargs[idx]["c"]


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


def _check_auc(display, constructor_name):
    roc_auc_limit = 0.95679
    roc_auc_limit_multi = [0.97007, 0.985915, 0.980952]

    if constructor_name == "from_cv_results":
        for idx, roc_auc in enumerate(display.roc_auc):
            assert roc_auc == pytest.approx(roc_auc_limit_multi[idx])
    else:
        assert display.roc_auc == pytest.approx(roc_auc_limit)
        assert trapezoid(display.tpr, display.fpr) == pytest.approx(roc_auc_limit)


@pytest.mark.parametrize("response_method", ["predict_proba", "decision_function"])
@pytest.mark.parametrize(
    "constructor_name", ["from_estimator", "from_predictions", "from_cv_results"]
)
def test_plot_roc_curve_pos_label(pyplot, response_method, constructor_name):
    # check that we can provide the positive label and display the proper
    # statistics
    X, y = load_breast_cancer(return_X_y=True)
    # create an highly imbalanced
    idx_positive = np.flatnonzero(y == 1)
    idx_negative = np.flatnonzero(y == 0)
    idx_selected = np.hstack([idx_negative, idx_positive[:25]])
    X, y = X[idx_selected], y[idx_selected]
    X, y = shuffle(X, y, random_state=42)
    # only use 2 features to make the problem even harder
    X = X[:, :2]
    y = np.array(["cancer" if c == 1 else "not cancer" for c in y], dtype=object)
    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        stratify=y,
        random_state=0,
    )

    classifier = LogisticRegression()
    classifier.fit(X_train, y_train)
    cv_results = cross_validate(
        LogisticRegression(), X, y, cv=3, return_estimator=True, return_indices=True
    )

    # Sanity check to be sure the positive class is `classes_[0]`
    # Class imbalance ensures a large difference in prediction values between classes,
    # allowing us to catch errors when we switch `pos_label`
    assert classifier.classes_.tolist() == ["cancer", "not cancer"]

    y_score = getattr(classifier, response_method)(X_test)
    # we select the corresponding probability columns or reverse the decision
    # function otherwise
    y_score_cancer = -1 * y_score if y_score.ndim == 1 else y_score[:, 0]
    y_score_not_cancer = y_score if y_score.ndim == 1 else y_score[:, 1]

    pos_label = "cancer"
    y_score = y_score_cancer
    if constructor_name == "from_estimator":
        display = RocCurveDisplay.from_estimator(
            classifier,
            X_test,
            y_test,
            pos_label=pos_label,
            response_method=response_method,
        )
    elif constructor_name == "from_predictions":
        display = RocCurveDisplay.from_predictions(
            y_test,
            y_score,
            pos_label=pos_label,
        )
    else:
        display = RocCurveDisplay.from_cv_results(
            cv_results,
            X,
            y,
            response_method=response_method,
            pos_label=pos_label,
        )

    _check_auc(display, constructor_name)

    pos_label = "not cancer"
    y_score = y_score_not_cancer
    if constructor_name == "from_estimator":
        display = RocCurveDisplay.from_estimator(
            classifier,
            X_test,
            y_test,
            response_method=response_method,
            pos_label=pos_label,
        )
    elif constructor_name == "from_predictions":
        display = RocCurveDisplay.from_predictions(
            y_test,
            y_score,
            pos_label=pos_label,
        )
    else:
        display = RocCurveDisplay.from_cv_results(
            cv_results,
            X,
            y,
            response_method=response_method,
            pos_label=pos_label,
        )

    _check_auc(display, constructor_name)


# TODO(1.9): remove
def test_y_score_and_y_pred_specified_error():
    """Check that an error is raised when both y_score and y_pred are specified."""
    y_true = np.array([0, 1, 1, 0])
    y_score = np.array([0.1, 0.4, 0.35, 0.8])
    y_pred = np.array([0.2, 0.3, 0.5, 0.1])

    with pytest.raises(
        ValueError, match="`y_pred` and `y_score` cannot be both specified"
    ):
        RocCurveDisplay.from_predictions(y_true, y_score=y_score, y_pred=y_pred)


# TODO(1.9): remove
def test_y_pred_deprecation_warning(pyplot):
    """Check that a warning is raised when y_pred is specified."""
    y_true = np.array([0, 1, 1, 0])
    y_score = np.array([0.1, 0.4, 0.35, 0.8])

    with pytest.warns(FutureWarning, match="y_pred is deprecated in 1.7"):
        display_y_pred = RocCurveDisplay.from_predictions(y_true, y_pred=y_score)

    assert_allclose(display_y_pred.fpr, [0, 0.5, 0.5, 1])
    assert_allclose(display_y_pred.tpr, [0, 0, 1, 1])

    display_y_score = RocCurveDisplay.from_predictions(y_true, y_score)
    assert_allclose(display_y_score.fpr, [0, 0.5, 0.5, 1])
    assert_allclose(display_y_score.tpr, [0, 0, 1, 1])


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
