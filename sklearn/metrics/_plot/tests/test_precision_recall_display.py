from collections import Counter

import numpy as np
import pytest
from scipy.integrate import trapezoid

from sklearn.compose import make_column_transformer
from sklearn.datasets import load_breast_cancer, make_classification
from sklearn.exceptions import NotFittedError
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    PrecisionRecallDisplay,
    average_precision_score,
    precision_recall_curve,
)
from sklearn.metrics._plot.tests.test_common_curve_display import (
    _check_from_cv_results_legend_label,
    _check_from_cv_results_param_validation,
    _check_plot_legend_label,
    _check_validate_plot_params,
)
from sklearn.model_selection import cross_validate, train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.utils import _safe_indexing, shuffle


def _check_figure_axes_and_labels(display, pos_label):
    """Check mpl figure and axes are correct."""
    import matplotlib as mpl

    assert isinstance(display.ax_, mpl.axes.Axes)
    assert isinstance(display.figure_, mpl.figure.Figure)

    assert display.ax_.get_xlabel() == f"Recall (Positive label: {pos_label})"
    assert display.ax_.get_ylabel() == f"Precision (Positive label: {pos_label})"
    assert display.ax_.get_adjustable() == "box"
    assert display.ax_.get_aspect() in ("equal", 1.0)
    assert display.ax_.get_xlim() == display.ax_.get_ylim() == (-0.01, 1.01)


@pytest.mark.parametrize("constructor_name", ["from_estimator", "from_predictions"])
@pytest.mark.parametrize("response_method", ["predict_proba", "decision_function"])
@pytest.mark.parametrize("drop_intermediate", [True, False])
@pytest.mark.parametrize("with_sample_weight", [True, False])
def test_precision_recall_display_plotting(
    pyplot,
    constructor_name,
    response_method,
    drop_intermediate,
    with_sample_weight,
):
    """Check the overall plotting rendering."""
    X, y = make_classification(n_classes=2, n_samples=50, random_state=0)
    pos_label = 1

    classifier = LogisticRegression().fit(X, y)
    classifier.fit(X, y)

    if with_sample_weight:
        rng = np.random.RandomState(42)
        sample_weight = rng.randint(1, 4, size=(X.shape[0]))
    else:
        sample_weight = None

    y_pred = getattr(classifier, response_method)(X)
    y_pred = y_pred if y_pred.ndim == 1 else y_pred[:, 1]

    # safe guard for the binary if/else construction
    assert constructor_name in ("from_estimator", "from_predictions")

    if constructor_name == "from_estimator":
        display = PrecisionRecallDisplay.from_estimator(
            classifier,
            X,
            y,
            sample_weight=sample_weight,
            response_method=response_method,
            drop_intermediate=drop_intermediate,
        )
    else:
        display = PrecisionRecallDisplay.from_predictions(
            y,
            y_pred,
            sample_weight=sample_weight,
            pos_label=pos_label,
            drop_intermediate=drop_intermediate,
        )

    precision, recall, _ = precision_recall_curve(
        y,
        y_pred,
        pos_label=pos_label,
        sample_weight=sample_weight,
        drop_intermediate=drop_intermediate,
    )
    average_precision = average_precision_score(
        y, y_pred, pos_label=pos_label, sample_weight=sample_weight
    )

    np.testing.assert_allclose(display.precision, precision)
    np.testing.assert_allclose(display.recall, recall)
    assert display.average_precision == pytest.approx(average_precision)

    import matplotlib as mpl

    _check_figure_axes_and_labels(display, pos_label)
    assert isinstance(display.line_, mpl.lines.Line2D)

    # plotting passing some new parameters
    display.plot(name="MySpecialEstimator", curve_kwargs={"alpha": 0.8})
    expected_label = f"MySpecialEstimator (AP = {average_precision:0.2f})"
    assert display.line_.get_label() == expected_label
    assert display.line_.get_alpha() == pytest.approx(0.8)

    # Check that the chance level line is not plotted by default
    assert display.chance_level_ is None


@pytest.mark.parametrize("response_method", ["predict_proba", "decision_function"])
@pytest.mark.parametrize("drop_intermediate", [True, False])
@pytest.mark.parametrize("with_sample_weight", [True, False])
def test_precision_recall_display_from_cv_results_plotting(
    pyplot, response_method, drop_intermediate, with_sample_weight
):
    """Check the overall plotting of `from_cv_results`."""
    X, y = make_classification(n_classes=2, n_samples=50, random_state=0)
    pos_label = 1

    cv_results = cross_validate(
        LogisticRegression(), X, y, cv=3, return_estimator=True, return_indices=True
    )

    if with_sample_weight:
        rng = np.random.RandomState(42)
        sample_weight = rng.randint(1, 4, size=(X.shape[0]))
    else:
        sample_weight = None

    display = PrecisionRecallDisplay.from_cv_results(
        cv_results,
        X,
        y,
        sample_weight=sample_weight,
        response_method=response_method,
        drop_intermediate=drop_intermediate,
        pos_label=pos_label,
    )

    for idx, (estimator, test_indices) in enumerate(
        zip(cv_results["estimator"], cv_results["indices"]["test"])
    ):
        y_true = _safe_indexing(y, test_indices)
        y_pred = getattr(estimator, response_method)(_safe_indexing(X, test_indices))
        y_pred = y_pred if y_pred.ndim == 1 else y_pred[:, 1]
        sample_weight_test = (
            _safe_indexing(sample_weight, test_indices)
            if sample_weight is not None
            else None
        )
        precision, recall, _ = precision_recall_curve(
            y_true,
            y_pred,
            pos_label=pos_label,
            drop_intermediate=drop_intermediate,
            sample_weight=sample_weight_test,
        )
        average_precision = average_precision_score(
            y_true, y_pred, pos_label=pos_label, sample_weight=sample_weight_test
        )

        np.testing.assert_allclose(display.precision[idx], precision)
        np.testing.assert_allclose(display.recall[idx], recall)
        assert display.average_precision[idx] == pytest.approx(average_precision)

        import matplotlib as mpl

        assert isinstance(display.line_[idx], mpl.lines.Line2D)

    _check_figure_axes_and_labels(display, pos_label)
    # Check that the chance level line is not plotted by default
    assert display.chance_level_ is None


@pytest.mark.parametrize(
    "params, err_msg",
    [
        (
            {
                "precision": [np.array([1, 0.5, 0]), np.array([1, 0.5, 0])],
                "recall": [np.array([0, 0.5, 1])],
                "average_precision": None,
                "prevalence_pos_label": None,
                "name": None,
            },
            "self.precision and self.recall from `PrecisionRecallDisplay`",
        ),
        (
            {
                "precision": [np.array([1, 0.5, 0])],
                "recall": [np.array([0, 0.5, 1]), np.array([0, 0.5, 1])],
                "average_precision": [0.8, 0.9],
                "prevalence_pos_label": None,
                "name": None,
            },
            "self.precision, self.recall and self.average_precision",
        ),
        (
            {
                "precision": [np.array([1, 0.5, 0])],
                "recall": [np.array([0, 0.5, 1]), np.array([0, 0.5, 1])],
                "average_precision": [0.8, 0.9],
                "prevalence_pos_label": [0.5, 0.5, 0.5],
                "name": None,
            },
            (
                "self.precision, self.recall, self.average_precision and "
                "self.prevalence_pos_label"
            ),
        ),
        (
            {
                "precision": [np.array([1, 0.5, 0]), np.array([1, 0.5, 0])],
                "recall": [np.array([0, 0.5, 1]), np.array([0, 0.5, 1])],
                "average_precision": [0.8],
                "prevalence_pos_label": [0.8, 0.6, 0.5],
                "name": None,
            },
            (
                "Got: self.precision: 2, self.recall: 2, self.average_precision: 1, "
                "self.prevalence_pos_label: 3"
            ),
        ),
        (
            {
                "precision": [np.array([1, 0.5, 0]), np.array([1, 0.5, 0])],
                "recall": [np.array([0, 0.5, 1]), np.array([0, 0.5, 1])],
                "average_precision": [0.8, 0.9],
                "prevalence_pos_label": None,
                "name": ["curve1", "curve2", "curve3"],
            },
            (
                "self.precision, self.recall, self.average_precision and 'name' "
                r"\(or self.name\)"
            ),
        ),
        (
            {
                "precision": [np.array([1, 0.5, 0]), np.array([1, 0.5, 0])],
                "recall": [np.array([0, 0.5, 1]), np.array([0, 0.5, 1])],
                "average_precision": [0.8, 0.9],
                "prevalence_pos_label": [0.5, 0.4],
                # List of length 1 is always allowed
                "name": ["curve1"],
            },
            None,
        ),
    ],
)
def test_precison_recall_plot_parameter_length_validation(pyplot, params, err_msg):
    """Check `plot` parameter length validation performed correctly."""
    display = PrecisionRecallDisplay(**params)
    if err_msg:
        with pytest.raises(ValueError, match=err_msg):
            display.plot()
    else:
        # No error should be raised
        display.plot()


def test_validate_plot_params(pyplot):
    """Check `_validate_plot_params` returns the correct variables."""
    display_args = {
        "precision": np.array([1, 0.5, 0]),
        "recall": [np.array([0, 0.5, 1])],
        "average_precision": None,
        "name": "test_curve",
        "prevalence_pos_label": 0.5,
    }

    _check_validate_plot_params(PrecisionRecallDisplay, display_args)


def test_precision_recall_from_cv_results_param_validation(pyplot):
    """Check parameter validation is correct."""
    data = make_classification(n_classes=2, n_samples=50, random_state=0)
    _check_from_cv_results_param_validation(data, PrecisionRecallDisplay)


@pytest.mark.parametrize("chance_level_kw", [None, {"color": "r"}, {"c": "r"}])
@pytest.mark.parametrize("constructor_name", ["from_estimator", "from_predictions"])
def test_precision_recall_chance_level_line(pyplot, chance_level_kw, constructor_name):
    """Check chance level plotting behavior, for `from_estimator`/`from_predictions`."""
    X, y = make_classification(n_classes=2, n_samples=50, random_state=0)
    pos_prevalence = Counter(y)[1] / len(y)

    lr = LogisticRegression()
    y_pred = lr.fit(X, y).predict_proba(X)[:, 1]

    if constructor_name == "from_estimator":
        display = PrecisionRecallDisplay.from_estimator(
            lr,
            X,
            y,
            plot_chance_level=True,
            chance_level_kw=chance_level_kw,
        )
    else:
        display = PrecisionRecallDisplay.from_predictions(
            y,
            y_pred,
            plot_chance_level=True,
            chance_level_kw=chance_level_kw,
        )

    import matplotlib as mpl

    assert isinstance(display.chance_level_, mpl.lines.Line2D)
    assert tuple(display.chance_level_.get_xdata()) == (0, 1)
    assert tuple(display.chance_level_.get_ydata()) == (pos_prevalence, pos_prevalence)

    # Checking for chance level line styles
    if chance_level_kw is None:
        assert display.chance_level_.get_color() == "k"
    else:
        assert display.chance_level_.get_color() == "r"

    assert display.chance_level_.get_label() == f"Chance level (AP = {pos_prevalence})"


@pytest.mark.parametrize("chance_level_kw", [None, {"color": "r"}, {"c": "r"}])
def test_precision_recall_chance_level_line_from_cv_results(pyplot, chance_level_kw):
    """Check chance level plotting behavior for `from_cv_results`."""
    # Note a separate chance line is plotted for each cv split
    X, y = make_classification(n_classes=2, n_samples=50, random_state=0)
    n_cv = 3
    cv_results = cross_validate(
        LogisticRegression(), X, y, cv=n_cv, return_estimator=True, return_indices=True
    )

    display = PrecisionRecallDisplay.from_cv_results(
        cv_results,
        X,
        y,
        plot_chance_level=True,
        chance_level_kwargs=chance_level_kw,
    )

    import matplotlib as mpl

    pos_prevalence_folds = []
    for idx in range(n_cv):
        assert isinstance(display.chance_level_[idx], mpl.lines.Line2D)
        assert tuple(display.chance_level_[idx].get_xdata()) == (0, 1)
        test_indices = cv_results["indices"]["test"][idx]
        pos_prevalence = Counter(_safe_indexing(y, test_indices))[1] / len(test_indices)
        pos_prevalence_folds.append(pos_prevalence)
        assert tuple(display.chance_level_[idx].get_ydata()) == (
            pos_prevalence,
            pos_prevalence,
        )

        # Checking for chance level line styles
        if chance_level_kw is None:
            assert display.chance_level_[idx].get_color() == "k"
        else:
            assert display.chance_level_[idx].get_color() == "r"

    for idx in range(n_cv):
        # Only the first chance line should have a label
        if idx == 0:
            assert display.chance_level_[idx].get_label() == (
                f"Chance level (AP = {np.mean(pos_prevalence_folds):0.2f} +/- "
                f"{np.std(pos_prevalence_folds):0.2f})"
            )
        else:
            assert display.chance_level_[idx].get_label() == f"_child{3 + idx}"


@pytest.mark.parametrize(
    "constructor_name, default_label",
    [
        ("from_estimator", "LogisticRegression (AP = {:.2f})"),
        ("from_predictions", "Classifier (AP = {:.2f})"),
    ],
)
def test_precision_recall_display_name(pyplot, constructor_name, default_label):
    """Check the behaviour of the name parameters"""
    X, y = make_classification(n_classes=2, n_samples=100, random_state=0)
    pos_label = 1

    classifier = LogisticRegression().fit(X, y)
    classifier.fit(X, y)

    y_pred = classifier.predict_proba(X)[:, pos_label]

    # safe guard for the binary if/else construction
    assert constructor_name in ("from_estimator", "from_predictions")

    if constructor_name == "from_estimator":
        display = PrecisionRecallDisplay.from_estimator(classifier, X, y)
    else:
        display = PrecisionRecallDisplay.from_predictions(
            y, y_pred, pos_label=pos_label
        )

    average_precision = average_precision_score(y, y_pred, pos_label=pos_label)

    # check that the default name is used
    assert display.line_.get_label() == default_label.format(average_precision)

    # check that the name can be set
    display.plot(name="MySpecialEstimator")
    assert (
        display.line_.get_label()
        == f"MySpecialEstimator (AP = {average_precision:.2f})"
    )


@pytest.mark.parametrize("average_precision", [[1.0, 1.0, 1.0], None])
@pytest.mark.parametrize(
    "curve_kwargs",
    [None, {"color": "red"}, [{"c": "red"}, {"c": "green"}, {"c": "yellow"}]],
)
@pytest.mark.parametrize("name", [None, "single", ["one", "two", "three"]])
def test_precision_recall_plot_legend_label(
    pyplot, name, curve_kwargs, average_precision
):
    """Check legend label correct with all `curve_kwargs`, `name` combinations."""
    precision = [np.array([1, 0.5, 0]), np.array([1, 0.5, 0]), np.array([1, 0.5, 0])]
    recall = [np.array([0, 0.5, 1]), np.array([0, 0.5, 1]), np.array([0, 0.5, 1])]

    _check_plot_legend_label(
        PrecisionRecallDisplay,
        {
            "precision": precision,
            "recall": recall,
            "average_precision": average_precision,
        },
        name,
        curve_kwargs,
        average_precision,
        "AP",
    )


@pytest.mark.parametrize(
    "curve_kwargs",
    [None, {"color": "red"}, [{"c": "red"}, {"c": "green"}, {"c": "yellow"}]],
)
@pytest.mark.parametrize("name", [None, "single", ["one", "two", "three"]])
def test_precision_recall_from_cv_results_legend_label(pyplot, name, curve_kwargs):
    """Check legend label correct with all `curve_kwargs`, `name` combinations."""
    X, y = make_classification(n_classes=2, n_samples=50, random_state=0)
    cv_results = cross_validate(
        LogisticRegression(), X, y, cv=3, return_estimator=True, return_indices=True
    )

    _check_from_cv_results_legend_label(
        PrecisionRecallDisplay,
        cv_results,
        X,
        y,
        name,
        curve_kwargs,
        [0.97, 1.00, 1.00],
        "AP",
    )


@pytest.mark.parametrize(
    "clf",
    [
        make_pipeline(StandardScaler(), LogisticRegression()),
        make_pipeline(
            make_column_transformer((StandardScaler(), [0, 1])), LogisticRegression()
        ),
    ],
)
def test_precision_recall_display_pipeline(pyplot, clf):
    X, y = make_classification(n_classes=2, n_samples=50, random_state=0)
    with pytest.raises(NotFittedError):
        PrecisionRecallDisplay.from_estimator(clf, X, y)
    clf.fit(X, y)
    display = PrecisionRecallDisplay.from_estimator(clf, X, y)
    assert display.name == clf.__class__.__name__


def test_precision_recall_display_string_labels(pyplot):
    # regression test #15738
    cancer = load_breast_cancer()
    X, y = cancer.data, cancer.target_names[cancer.target]

    lr = make_pipeline(StandardScaler(), LogisticRegression())
    lr.fit(X, y)
    for klass in cancer.target_names:
        assert klass in lr.classes_
    display = PrecisionRecallDisplay.from_estimator(lr, X, y)

    y_pred = lr.predict_proba(X)[:, 1]
    avg_prec = average_precision_score(y, y_pred, pos_label=lr.classes_[1])

    assert display.average_precision == pytest.approx(avg_prec)
    assert display.name == lr.__class__.__name__

    err_msg = r"y_true takes value in {'benign', 'malignant'}"
    with pytest.raises(ValueError, match=err_msg):
        PrecisionRecallDisplay.from_predictions(y, y_pred)

    display = PrecisionRecallDisplay.from_predictions(
        y, y_pred, pos_label=lr.classes_[1]
    )
    assert display.average_precision == pytest.approx(avg_prec)


@pytest.mark.parametrize(
    "average_precision, name, expected_label",
    [
        (0.9, None, "AP = 0.90"),
        (None, "my_est", "my_est"),
        (0.8, "my_est2", "my_est2 (AP = 0.80)"),
    ],
)
def test_default_labels(pyplot, average_precision, name, expected_label):
    """Check the default labels used in the display."""
    precision = np.array([1, 0.5, 0])
    recall = np.array([0, 0.5, 1])
    display = PrecisionRecallDisplay(
        precision,
        recall,
        average_precision=average_precision,
        name=name,
    )
    display.plot()
    assert display.line_.get_label() == expected_label


@pytest.mark.parametrize("constructor_name", ["from_estimator", "from_predictions"])
@pytest.mark.parametrize("response_method", ["predict_proba", "decision_function"])
def test_plot_precision_recall_pos_label(pyplot, constructor_name, response_method):
    # check that we can provide the positive label and display the proper
    # statistics
    X, y = load_breast_cancer(return_X_y=True)
    # create an highly imbalanced version of the breast cancer dataset
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

    # sanity check to be sure the positive class is classes_[0] and that we
    # are betrayed by the class imbalance
    assert classifier.classes_.tolist() == ["cancer", "not cancer"]

    y_pred = getattr(classifier, response_method)(X_test)
    # we select the corresponding probability columns or reverse the decision
    #  function otherwise
    y_pred_cancer = -1 * y_pred if y_pred.ndim == 1 else y_pred[:, 0]
    y_pred_not_cancer = y_pred if y_pred.ndim == 1 else y_pred[:, 1]

    if constructor_name == "from_estimator":
        display = PrecisionRecallDisplay.from_estimator(
            classifier,
            X_test,
            y_test,
            pos_label="cancer",
            response_method=response_method,
        )
    else:
        display = PrecisionRecallDisplay.from_predictions(
            y_test,
            y_pred_cancer,
            pos_label="cancer",
        )
    # we should obtain the statistics of the "cancer" class
    avg_prec_limit = 0.65
    assert display.average_precision < avg_prec_limit
    assert -trapezoid(display.precision, display.recall) < avg_prec_limit

    # otherwise we should obtain the statistics of the "not cancer" class
    if constructor_name == "from_estimator":
        display = PrecisionRecallDisplay.from_estimator(
            classifier,
            X_test,
            y_test,
            response_method=response_method,
            pos_label="not cancer",
        )
    else:
        display = PrecisionRecallDisplay.from_predictions(
            y_test,
            y_pred_not_cancer,
            pos_label="not cancer",
        )
    avg_prec_limit = 0.95
    assert display.average_precision > avg_prec_limit
    assert -trapezoid(display.precision, display.recall) > avg_prec_limit


@pytest.mark.parametrize("constructor_name", ["from_estimator", "from_predictions"])
def test_precision_recall_prevalence_pos_label_reusable(pyplot, constructor_name):
    # Check that even if one passes plot_chance_level=False the first time
    # one can still call disp.plot with plot_chance_level=True and get the
    # chance level line
    X, y = make_classification(n_classes=2, n_samples=50, random_state=0)

    lr = LogisticRegression()
    y_pred = lr.fit(X, y).predict_proba(X)[:, 1]

    if constructor_name == "from_estimator":
        display = PrecisionRecallDisplay.from_estimator(
            lr, X, y, plot_chance_level=False
        )
    else:
        display = PrecisionRecallDisplay.from_predictions(
            y, y_pred, plot_chance_level=False
        )
    assert display.chance_level_ is None

    import matplotlib as mpl

    # When calling from_estimator or from_predictions,
    # prevalence_pos_label should have been set, so that directly
    # calling plot_chance_level=True should plot the chance level line
    display.plot(plot_chance_level=True)
    assert isinstance(display.chance_level_, mpl.lines.Line2D)


def test_precision_recall_raise_no_prevalence(pyplot):
    # Check that raises correctly when plotting chance level with
    # no prvelance_pos_label is provided
    precision = np.array([1, 0.5, 0])
    recall = np.array([0, 0.5, 1])
    display = PrecisionRecallDisplay(precision, recall)

    msg = (
        "You must provide prevalence_pos_label when constructing the "
        "PrecisionRecallDisplay object in order to plot the chance "
        "level line. Alternatively, you may use "
        "PrecisionRecallDisplay.from_estimator or "
        "PrecisionRecallDisplay.from_predictions "
        "to automatically set prevalence_pos_label"
    )

    with pytest.raises(ValueError, match=msg):
        display.plot(plot_chance_level=True)


@pytest.mark.parametrize("despine", [True, False])
@pytest.mark.parametrize(
    "constructor_name", ["from_estimator", "from_predictions", "from_cv_results"]
)
def test_plot_precision_recall_despine(pyplot, despine, constructor_name):
    # Check that the despine keyword is working correctly
    X, y = make_classification(n_classes=2, n_samples=50, random_state=0)

    clf = LogisticRegression().fit(X, y)
    clf.fit(X, y)
    cv_results = cross_validate(
        LogisticRegression(), X, y, cv=3, return_estimator=True, return_indices=True
    )

    y_pred = clf.decision_function(X)

    if constructor_name == "from_estimator":
        display = PrecisionRecallDisplay.from_estimator(clf, X, y, despine=despine)
    elif constructor_name == "from_predictions":
        display = PrecisionRecallDisplay.from_predictions(y, y_pred, despine=despine)
    else:
        display = PrecisionRecallDisplay.from_cv_results(
            cv_results, X, y, despine=despine
        )

    for s in ["top", "right"]:
        assert display.ax_.spines[s].get_visible() is not despine

    if despine:
        for s in ["bottom", "left"]:
            assert display.ax_.spines[s].get_bounds() == (0, 1)
