from collections import Counter

import numpy as np
import pytest
from numpy.testing import assert_allclose

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
    _check_pos_label_statistics,
)
from sklearn.model_selection import cross_validate
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.utils import _safe_indexing
from sklearn.utils._response import _get_response_values_binary


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

    y_score = getattr(classifier, response_method)(X)
    y_score = y_score if y_score.ndim == 1 else y_score[:, pos_label]

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
            y_score,
            sample_weight=sample_weight,
            pos_label=pos_label,
            drop_intermediate=drop_intermediate,
        )

    precision, recall, _ = precision_recall_curve(
        y,
        y_score,
        pos_label=pos_label,
        sample_weight=sample_weight,
        drop_intermediate=drop_intermediate,
    )
    average_precision = average_precision_score(
        y, y_score, pos_label=pos_label, sample_weight=sample_weight
    )

    assert_allclose(display.precision, precision)
    assert_allclose(display.recall, recall)
    assert display.average_precision == pytest.approx(average_precision)

    import matplotlib as mpl

    _check_figure_axes_and_labels(display, pos_label)
    assert isinstance(display.line_, mpl.lines.Line2D)
    # Check default curve kwarg
    assert display.line_.get_drawstyle() == "steps-post"

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
        y_score = getattr(estimator, response_method)(_safe_indexing(X, test_indices))
        y_score = y_score if y_score.ndim == 1 else y_score[:, 1]
        sample_weight_test = (
            _safe_indexing(sample_weight, test_indices)
            if sample_weight is not None
            else None
        )
        precision, recall, _ = precision_recall_curve(
            y_true,
            y_score,
            pos_label=pos_label,
            drop_intermediate=drop_intermediate,
            sample_weight=sample_weight_test,
        )
        average_precision = average_precision_score(
            y_true, y_score, pos_label=pos_label, sample_weight=sample_weight_test
        )

        assert_allclose(display.precision[idx], precision)
        assert_allclose(display.recall[idx], recall)
        assert display.average_precision[idx] == pytest.approx(average_precision)

        import matplotlib as mpl

        assert isinstance(display.line_[idx], mpl.lines.Line2D)
        # Check default curve kwarg
        assert display.line_[idx].get_drawstyle() == "steps-post"

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


@pytest.mark.parametrize("plot_chance_level", [True, False])
@pytest.mark.parametrize("chance_level_kw", [None, {"color": "r"}, {"c": "r"}])
@pytest.mark.parametrize("constructor_name", ["from_estimator", "from_predictions"])
def test_precision_recall_chance_level_line(
    pyplot, plot_chance_level, chance_level_kw, constructor_name
):
    """Check chance level plotting behavior, for `from_estimator`/`from_predictions`."""
    X, y = make_classification(n_classes=2, n_samples=50, random_state=0)
    pos_prevalence = Counter(y)[1] / len(y)

    lr = LogisticRegression()
    y_score = lr.fit(X, y).predict_proba(X)[:, 1]

    if constructor_name == "from_estimator":
        display = PrecisionRecallDisplay.from_estimator(
            lr,
            X,
            y,
            plot_chance_level=plot_chance_level,
            chance_level_kw=chance_level_kw,
        )
    else:
        display = PrecisionRecallDisplay.from_predictions(
            y,
            y_score,
            plot_chance_level=plot_chance_level,
            chance_level_kw=chance_level_kw,
        )

    import matplotlib as mpl

    if not plot_chance_level:
        assert display.chance_level_ is None
        # Early return if chance level not plotted
        return

    assert isinstance(display.chance_level_, mpl.lines.Line2D)
    assert tuple(display.chance_level_.get_xdata()) == (0, 1)
    assert tuple(display.chance_level_.get_ydata()) == (pos_prevalence, pos_prevalence)

    # Checking for chance level line styles
    if chance_level_kw is None:
        assert display.chance_level_.get_color() == "k"
    else:
        assert display.chance_level_.get_color() == "r"

    assert display.chance_level_.get_label() == f"Chance level (AP = {pos_prevalence})"


@pytest.mark.parametrize("plot_chance_level", [True, False])
@pytest.mark.parametrize("chance_level_kw", [None, {"color": "r"}, {"c": "r"}])
def test_precision_recall_chance_level_line_from_cv_results(
    pyplot, plot_chance_level, chance_level_kw
):
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
        plot_chance_level=plot_chance_level,
        chance_level_kwargs=chance_level_kw,
    )

    import matplotlib as mpl

    if not plot_chance_level:
        assert display.chance_level_ is None
        # Early return if chance level not plotted
        return

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
        ("from_cv_results", "AP = {:.2f} +/- {:.2f}"),
    ],
)
def test_precision_recall_display_name(pyplot, constructor_name, default_label):
    """Check the behaviour of the name parameters"""
    X, y = make_classification(n_classes=2, n_samples=100, random_state=0)
    pos_label = 1

    classifier = LogisticRegression()
    n_cv = 3
    cv_results = cross_validate(
        classifier, X, y, cv=n_cv, return_estimator=True, return_indices=True
    )
    classifier.fit(X, y)
    y_score = classifier.predict_proba(X)[:, pos_label]

    if constructor_name == "from_estimator":
        display = PrecisionRecallDisplay.from_estimator(classifier, X, y)
    elif constructor_name == "from_predictions":
        display = PrecisionRecallDisplay.from_predictions(
            y, y_score, pos_label=pos_label
        )
    else:
        display = PrecisionRecallDisplay.from_cv_results(cv_results, X, y)

    if constructor_name == "from_cv_results":
        average_precision = []
        for idx in range(n_cv):
            test_indices = cv_results["indices"]["test"][idx]
            y_score, _ = _get_response_values_binary(
                cv_results["estimator"][idx],
                _safe_indexing(X, test_indices),
                response_method="auto",
            )
            average_precision.append(
                average_precision_score(
                    _safe_indexing(y, test_indices), y_score, pos_label=pos_label
                )
            )
        # By default, only the first curve is labelled
        assert display.line_[0].get_label() == default_label.format(
            np.mean(average_precision), np.std(average_precision)
        )

        # check that the name can be set
        display.plot(name="MySpecialEstimator")
        # Sets only first labelled curve
        assert display.line_[0].get_label() == (
            f"MySpecialEstimator (AP = {np.mean(average_precision):.2f} +/- "
            f"{np.std(average_precision):.2f})"
        )
    else:
        average_precision = average_precision_score(y, y_score, pos_label=pos_label)

        # check that the default name is used
        assert display.line_.get_label() == default_label.format(average_precision)

        # check that the name can be set
        display.plot(name="MySpecialEstimator")
        assert (
            display.line_.get_label()
            == f"MySpecialEstimator (AP = {average_precision:.2f})"
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
    n_cv = 3
    cv_results = cross_validate(
        lr, X, y, cv=n_cv, return_estimator=True, return_indices=True
    )
    lr.fit(X, y)
    for klass in cancer.target_names:
        assert klass in lr.classes_

    # `from_estimator`
    display = PrecisionRecallDisplay.from_estimator(lr, X, y)

    y_score = lr.predict_proba(X)[:, 1]
    avg_prec = average_precision_score(y, y_score, pos_label=lr.classes_[1])

    assert display.average_precision == pytest.approx(avg_prec)
    assert display.name == lr.__class__.__name__

    # `from_predictions`
    err_msg = r"y_true takes value in {'benign', 'malignant'}"
    with pytest.raises(ValueError, match=err_msg):
        PrecisionRecallDisplay.from_predictions(y, y_score)

    display = PrecisionRecallDisplay.from_predictions(
        y, y_score, pos_label=lr.classes_[1]
    )
    assert display.average_precision == pytest.approx(avg_prec)

    # `from_cv_results`
    display = PrecisionRecallDisplay.from_cv_results(cv_results, X, y)
    average_precision = []
    for idx in range(n_cv):
        test_indices = cv_results["indices"]["test"][idx]
        y_pred, _ = _get_response_values_binary(
            cv_results["estimator"][idx],
            _safe_indexing(X, test_indices),
            response_method="auto",
        )
        # Note `pos_label` cannot be `None` (default=1), unlike other metrics
        average_precision.append(
            average_precision_score(
                _safe_indexing(y, test_indices),
                y_pred,
                pos_label=cv_results["estimator"][idx].classes_[1],
            )
        )
    assert_allclose(display.average_precision, average_precision)


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
    """Check switching `pos_label` give correct statistics, using imbalanced data."""

    def _check_average_precision(display, constructor_name, pos_label):
        if pos_label == "cancer":
            avg_prec_limit = 0.6338
            avg_prec_limit_multi = [0.8189, 0.8802, 0.8795]
        else:
            avg_prec_limit = 0.9953
            avg_prec_limit_multi = [0.9966, 0.9984, 0.9976]

        def average_precision_uninterpolated(precision, recall):
            return -np.sum(np.diff(recall) * np.array(precision)[:-1])

        if constructor_name == "from_cv_results":
            for idx, average_precision in enumerate(display.average_precision):
                assert average_precision == pytest.approx(
                    avg_prec_limit_multi[idx], rel=1e-3
                )
                assert average_precision_uninterpolated(
                    display.precision[idx], display.recall[idx]
                ) == pytest.approx(avg_prec_limit_multi[idx], rel=1e-3)
        else:
            assert display.average_precision == pytest.approx(avg_prec_limit, rel=1e-3)
            assert average_precision_uninterpolated(
                display.precision, display.recall
            ) == pytest.approx(avg_prec_limit, rel=1e-3)

    _check_pos_label_statistics(
        PrecisionRecallDisplay,
        response_method,
        constructor_name,
        _check_average_precision,
    )


@pytest.mark.parametrize(
    "constructor_name", ["from_estimator", "from_predictions", "from_cv_results"]
)
def test_precision_recall_prevalence_pos_label_reusable(pyplot, constructor_name):
    # Check that even if one passes plot_chance_level=False the first time
    # one can still call disp.plot with plot_chance_level=True and get the
    # chance level line
    X, y = make_classification(n_classes=2, n_samples=50, random_state=0)

    lr = LogisticRegression()
    n_cv = 3
    cv_results = cross_validate(
        lr, X, y, cv=n_cv, return_estimator=True, return_indices=True
    )
    y_score = lr.fit(X, y).predict_proba(X)[:, 1]

    if constructor_name == "from_estimator":
        display = PrecisionRecallDisplay.from_estimator(
            lr, X, y, plot_chance_level=False
        )
    elif constructor_name == "from_predictions":
        display = PrecisionRecallDisplay.from_predictions(
            y, y_score, plot_chance_level=False
        )
    else:
        display = PrecisionRecallDisplay.from_cv_results(
            cv_results, X, y, plot_chance_level=False
        )
    assert display.chance_level_ is None

    import matplotlib as mpl

    # When calling from_estimator or from_predictions,
    # prevalence_pos_label should have been set, so that directly
    # calling plot_chance_level=True should plot the chance level line
    display.plot(plot_chance_level=True)
    if constructor_name == "from_cv_results":
        for idx in range(n_cv):
            assert isinstance(display.chance_level_[idx], mpl.lines.Line2D)
    else:
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

    y_score = clf.decision_function(X)

    if constructor_name == "from_estimator":
        display = PrecisionRecallDisplay.from_estimator(clf, X, y, despine=despine)
    elif constructor_name == "from_predictions":
        display = PrecisionRecallDisplay.from_predictions(y, y_score, despine=despine)
    else:
        display = PrecisionRecallDisplay.from_cv_results(
            cv_results, X, y, despine=despine
        )

    for s in ["top", "right"]:
        assert display.ax_.spines[s].get_visible() is not despine

    if despine:
        for s in ["bottom", "left"]:
            assert display.ax_.spines[s].get_bounds() == (0, 1)


# TODO(1.10): remove
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
        PrecisionRecallDisplay.from_predictions(y_true, y_score=y_score, y_pred=y_pred)

    with pytest.warns(FutureWarning, match="y_pred was deprecated in 1.8"):
        PrecisionRecallDisplay.from_predictions(y_true, y_pred=y_score)
