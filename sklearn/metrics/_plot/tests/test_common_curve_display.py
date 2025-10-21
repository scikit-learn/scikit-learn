from collections.abc import Mapping

import numpy as np
import pytest

from sklearn.base import BaseEstimator, ClassifierMixin, clone
from sklearn.calibration import CalibrationDisplay
from sklearn.compose import make_column_transformer
from sklearn.datasets import load_breast_cancer, load_iris
from sklearn.exceptions import NotFittedError
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    ConfusionMatrixDisplay,
    DetCurveDisplay,
    PrecisionRecallDisplay,
    PredictionErrorDisplay,
    RocCurveDisplay,
)
from sklearn.model_selection import cross_validate, train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor
from sklearn.utils import shuffle


@pytest.fixture(scope="module")
def data():
    return load_iris(return_X_y=True)


@pytest.fixture(scope="module")
def data_binary(data):
    X, y = data
    return X[y < 2], y[y < 2]


@pytest.mark.parametrize(
    "Display",
    [CalibrationDisplay, DetCurveDisplay, PrecisionRecallDisplay, RocCurveDisplay],
)
def test_display_curve_error_classifier(pyplot, data, data_binary, Display):
    """Check that a proper error is raised when only binary classification is
    supported."""
    X, y = data
    X_binary, y_binary = data_binary
    clf = DecisionTreeClassifier().fit(X, y)

    # Case 1: multiclass classifier with multiclass target
    msg = "Expected 'estimator' to be a binary classifier. Got 3 classes instead."
    with pytest.raises(ValueError, match=msg):
        Display.from_estimator(clf, X, y)

    # Case 2: multiclass classifier with binary target
    with pytest.raises(ValueError, match=msg):
        Display.from_estimator(clf, X_binary, y_binary)

    # Case 3: binary classifier with multiclass target
    clf = DecisionTreeClassifier().fit(X_binary, y_binary)
    msg = "The target y is not binary. Got multiclass type of target."
    with pytest.raises(ValueError, match=msg):
        Display.from_estimator(clf, X, y)


@pytest.mark.parametrize(
    "Display",
    [CalibrationDisplay, DetCurveDisplay, PrecisionRecallDisplay, RocCurveDisplay],
)
def test_display_curve_error_regression(pyplot, data_binary, Display):
    """Check that we raise an error with regressor."""

    # Case 1: regressor
    X, y = data_binary
    regressor = DecisionTreeRegressor().fit(X, y)

    msg = "Expected 'estimator' to be a binary classifier. Got DecisionTreeRegressor"
    with pytest.raises(ValueError, match=msg):
        Display.from_estimator(regressor, X, y)

    # Case 2: regression target
    classifier = DecisionTreeClassifier().fit(X, y)
    # Force `y_true` to be seen as a regression problem
    y = y + 0.5
    msg = "The target y is not binary. Got continuous type of target."
    with pytest.raises(ValueError, match=msg):
        Display.from_estimator(classifier, X, y)
    with pytest.raises(ValueError, match=msg):
        Display.from_predictions(y, regressor.fit(X, y).predict(X))


@pytest.mark.parametrize(
    "response_method, msg",
    [
        (
            "predict_proba",
            "MyClassifier has none of the following attributes: predict_proba.",
        ),
        (
            "decision_function",
            "MyClassifier has none of the following attributes: decision_function.",
        ),
        (
            "auto",
            (
                "MyClassifier has none of the following attributes: predict_proba,"
                " decision_function."
            ),
        ),
        (
            "bad_method",
            "MyClassifier has none of the following attributes: bad_method.",
        ),
    ],
)
@pytest.mark.parametrize(
    "Display", [DetCurveDisplay, PrecisionRecallDisplay, RocCurveDisplay]
)
def test_display_curve_error_no_response(
    pyplot,
    data_binary,
    response_method,
    msg,
    Display,
):
    """Check that a proper error is raised when the response method requested
    is not defined for the given trained classifier."""
    X, y = data_binary

    class MyClassifier(ClassifierMixin, BaseEstimator):
        def fit(self, X, y):
            self.classes_ = [0, 1]
            return self

    clf = MyClassifier().fit(X, y)

    with pytest.raises(AttributeError, match=msg):
        Display.from_estimator(clf, X, y, response_method=response_method)


@pytest.mark.parametrize(
    "Display", [DetCurveDisplay, PrecisionRecallDisplay, RocCurveDisplay]
)
@pytest.mark.parametrize("constructor_name", ["from_estimator", "from_predictions"])
def test_display_curve_estimator_name_multiple_calls(
    pyplot,
    data_binary,
    Display,
    constructor_name,
):
    """Check that passing `name` when calling `plot` will overwrite the original name
    in the legend."""
    X, y = data_binary
    clf_name = "my hand-crafted name"
    clf = LogisticRegression().fit(X, y)
    y_pred = clf.predict_proba(X)[:, 1]

    # safe guard for the binary if/else construction
    assert constructor_name in ("from_estimator", "from_predictions")

    if constructor_name == "from_estimator":
        disp = Display.from_estimator(clf, X, y, name=clf_name)
    else:
        disp = Display.from_predictions(y, y_pred, name=clf_name)
    # TODO: Clean-up once `estimator_name` deprecated in all displays
    if Display in (PrecisionRecallDisplay, RocCurveDisplay):
        assert disp.name == clf_name
    else:
        assert disp.estimator_name == clf_name
    pyplot.close("all")
    disp.plot()
    assert clf_name in disp.line_.get_label()
    pyplot.close("all")
    clf_name = "another_name"
    disp.plot(name=clf_name)
    assert clf_name in disp.line_.get_label()


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
@pytest.mark.parametrize(
    "Display", [DetCurveDisplay, PrecisionRecallDisplay, RocCurveDisplay]
)
def test_display_curve_not_fitted_errors_old_name(pyplot, data_binary, clf, Display):
    """Check that a proper error is raised when the classifier is not
    fitted."""
    X, y = data_binary
    # clone since we parametrize the test and the classifier will be fitted
    # when testing the second and subsequent plotting function
    model = clone(clf)
    with pytest.raises(NotFittedError):
        Display.from_estimator(model, X, y)
    model.fit(X, y)
    disp = Display.from_estimator(model, X, y)
    assert model.__class__.__name__ in disp.line_.get_label()
    # TODO: Clean-up once `estimator_name` deprecated in all displays
    if Display in (PrecisionRecallDisplay, RocCurveDisplay):
        assert disp.name == model.__class__.__name__
    else:
        assert disp.estimator_name == model.__class__.__name__


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
@pytest.mark.parametrize("Display", [RocCurveDisplay])
def test_display_curve_not_fitted_errors(pyplot, data_binary, clf, Display):
    """Check that a proper error is raised when the classifier is not fitted."""
    X, y = data_binary
    # clone since we parametrize the test and the classifier will be fitted
    # when testing the second and subsequent plotting function
    model = clone(clf)
    with pytest.raises(NotFittedError):
        Display.from_estimator(model, X, y)
    model.fit(X, y)
    disp = Display.from_estimator(model, X, y)
    assert model.__class__.__name__ in disp.line_.get_label()
    assert disp.name == model.__class__.__name__


@pytest.mark.parametrize(
    "Display", [DetCurveDisplay, PrecisionRecallDisplay, RocCurveDisplay]
)
def test_display_curve_n_samples_consistency(pyplot, data_binary, Display):
    """Check the error raised when `y_pred` or `sample_weight` have inconsistent
    length."""
    X, y = data_binary
    classifier = DecisionTreeClassifier().fit(X, y)

    msg = "Found input variables with inconsistent numbers of samples"
    with pytest.raises(ValueError, match=msg):
        Display.from_estimator(classifier, X[:-2], y)
    with pytest.raises(ValueError, match=msg):
        Display.from_estimator(classifier, X, y[:-2])
    with pytest.raises(ValueError, match=msg):
        Display.from_estimator(classifier, X, y, sample_weight=np.ones(X.shape[0] - 2))


@pytest.mark.parametrize(
    "Display", [DetCurveDisplay, PrecisionRecallDisplay, RocCurveDisplay]
)
def test_display_curve_error_pos_label(pyplot, data_binary, Display):
    """Check consistence of error message when `pos_label` should be specified."""
    X, y = data_binary
    y = y + 10

    classifier = DecisionTreeClassifier().fit(X, y)
    y_pred = classifier.predict_proba(X)[:, -1]
    msg = r"y_true takes value in {10, 11} and pos_label is not specified"
    with pytest.raises(ValueError, match=msg):
        Display.from_predictions(y, y_pred)


@pytest.mark.parametrize(
    "Display",
    [
        CalibrationDisplay,
        DetCurveDisplay,
        PrecisionRecallDisplay,
        RocCurveDisplay,
        PredictionErrorDisplay,
        ConfusionMatrixDisplay,
    ],
)
@pytest.mark.parametrize(
    "constructor",
    ["from_predictions", "from_estimator"],
)
def test_classifier_display_curve_named_constructor_return_type(
    pyplot, data_binary, Display, constructor
):
    """Check that named constructors return the correct type when subclassed.

    Non-regression test for:
    https://github.com/scikit-learn/scikit-learn/pull/27675
    """
    X, y = data_binary

    # This can be anything - we just need to check the named constructor return
    # type so the only requirement here is instantiating the class without error
    y_pred = y

    classifier = LogisticRegression().fit(X, y)

    class SubclassOfDisplay(Display):
        pass

    if constructor == "from_predictions":
        curve = SubclassOfDisplay.from_predictions(y, y_pred)
    else:  # constructor == "from_estimator"
        curve = SubclassOfDisplay.from_estimator(classifier, X, y)

    assert isinstance(curve, SubclassOfDisplay)


def _check_validate_plot_params(display_class, display_args):
    """Helper to check `_validate_plot_params` returns the correct variables.

    `display_args` should be given in the same order as output by
    `_validate_plot_params`. All `display_args` should be for a single curve.
    """
    display = display_class(**display_args)
    results = display._validate_plot_params(ax=None, name=None)

    # Check if the number of parameters match
    assert len(results) == len(display_args)

    for idx, (param, value) in enumerate(display_args.items()):
        if param == "name":
            assert results[idx] == [value] if isinstance(value, str) else value
        elif value is None:
            assert results[idx] is None
        else:
            assert isinstance(results[idx], list)
            assert len(results[idx]) == 1


def _check_from_cv_results_param_validation(data, display_class):
    """Check parameter validation is correct."""
    X, y = data

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
            display_class.from_cv_results(cv_results, X, y)

    cv_results = cross_validate(
        LogisticRegression(), X, y, cv=3, return_estimator=True, return_indices=True
    )

    # `X` wrong length
    with pytest.raises(ValueError, match="`X` does not contain the correct"):
        display_class.from_cv_results(cv_results, X[:10, :], y)

    # `y` not binary
    y_multi = y.copy()
    y_multi[0] = 2
    with pytest.raises(ValueError, match="The target `y` is not binary."):
        display_class.from_cv_results(cv_results, X, y_multi)

    # input inconsistent length
    with pytest.raises(ValueError, match="Found input variables with inconsistent"):
        display_class.from_cv_results(cv_results, X, y[:10])
    with pytest.raises(ValueError, match="Found input variables with inconsistent"):
        display_class.from_cv_results(cv_results, X, y, sample_weight=[1, 2])

    # `name` is list while `curve_kwargs` is None or dict
    for curve_kwargs in (None, {"alpha": 0.2}):
        with pytest.raises(ValueError, match="To avoid labeling individual curves"):
            display_class.from_cv_results(
                cv_results,
                X,
                y,
                name=["one", "two", "three"],
                curve_kwargs=curve_kwargs,
            )

    # `curve_kwargs` incorrect length
    with pytest.raises(ValueError, match="`curve_kwargs` must be None, a dictionary"):
        display_class.from_cv_results(cv_results, X, y, curve_kwargs=[{"alpha": 1}])

    # `curve_kwargs` both alias provided
    with pytest.raises(TypeError, match="Got both c and"):
        display_class.from_cv_results(
            cv_results, X, y, curve_kwargs={"c": "blue", "color": "red"}
        )


def _check_plot_legend_label(
    display_class, display_args, name, curve_kwargs, auc_metric, auc_metric_name
):
    """Check legend label correct with all `curve_kwargs`, `name` combinations.

    Checks `from_estimator` and `from_predictions` methods, when plotting multiple
    curves.

    Parameters
    ----------
    display_class : class
        The display class to test (e.g., RocCurveDisplay, PrecisionRecallDisplay).

    display_args : dict
        Dictionary of arguments to instantiate the display class.

    name : str, list of str, or None
        The name parameter to pass to the plot method.

    curve_kwargs : dict, list of dict, or None
        The curve_kwargs parameter to pass to the plot method.

    auc_metric : list of float or None
        The area under curve metric value (e.g., ROC AUC, average precision) to be
        displayed in the legend. When a list, should be all 1.0.

    auc_metric_name : str
        The name of the metric to display in the legend (e.g., "AUC", "AP").
    """
    if not isinstance(curve_kwargs, list) and isinstance(name, list):
        with pytest.raises(ValueError, match="To avoid labeling individual curves"):
            display_class(**display_args).plot(name=name, curve_kwargs=curve_kwargs)
        return

    display = display_class(**display_args).plot(name=name, curve_kwargs=curve_kwargs)
    legend = display.ax_.get_legend()
    if legend is None:
        # No legend is created, exit test early
        assert name is None
        assert auc_metric is None
        return
    else:
        legend_labels = [text.get_text() for text in legend.get_texts()]

    if isinstance(curve_kwargs, list):
        # Multiple labels in legend
        assert len(legend_labels) == 3
        for idx, label in enumerate(legend_labels):
            if name is None:
                expected_label = f"{auc_metric_name} = 1.00" if auc_metric else None
                assert label == expected_label
            elif isinstance(name, str):
                expected_label = (
                    f"single ({auc_metric_name} = 1.00)" if auc_metric else "single"
                )
                assert label == expected_label
            else:
                # `name` is a list of different strings
                expected_label = (
                    f"{name[idx]} ({auc_metric_name} = 1.00)"
                    if auc_metric
                    else f"{name[idx]}"
                )
                assert label == expected_label
    else:
        # Single label in legend
        assert len(legend_labels) == 1
        if name is None:
            expected_label = (
                f"{auc_metric_name} = 1.00 +/- 0.00" if auc_metric else None
            )
            assert legend_labels[0] == expected_label
        else:
            # name is single string
            expected_label = (
                f"single ({auc_metric_name} = 1.00 +/- 0.00)"
                if auc_metric
                else "single"
            )
            assert legend_labels[0] == expected_label


def _check_from_cv_results_legend_label(
    display_class, cv_results, X, y, name, curve_kwargs, auc_metrics, auc_metric_name
):
    """Check legend label correct with all `curve_kwargs`, `name` combinations.

    This function verifies that the legend labels in a Display object created from
    cross-validation results are correctly formatted based on the provided parameters.

    Parameters
    ----------
    display_class : class
        The display class to test (e.g., RocCurveDisplay, PrecisionRecallDisplay).

    cv_results : dict
        The cross-validation results dictionary containing the estimators and indices.
        Should have 3 cv folds.

    X : array-like of shape (n_samples, n_features)
        The input samples.

    y : array-like of shape (n_samples,)
        The target values (binary).

    name : str, list of str, or None
        The `name` parameter to pass to `from_cv_results`.

    curve_kwargs : dict, list of dict, or None
        The `curve_kwargs` parameter to pass to `from_cv_results`.

    auc_metrics : list of float
        The area under curve metric values (e.g., ROC AUC, average precision) to be
        displayed in the legend.

    auc_metric_name : str
        The name of the metric to display in the legend (e.g., "AUC", "AP").
    """
    if not isinstance(curve_kwargs, list) and isinstance(name, list):
        with pytest.raises(ValueError, match="To avoid labeling individual curves"):
            display_class.from_cv_results(
                cv_results, X, y, name=name, curve_kwargs=curve_kwargs
            )
    else:
        display = display_class.from_cv_results(
            cv_results, X, y, name=name, curve_kwargs=curve_kwargs
        )

        legend = display.ax_.get_legend()
        legend_labels = [text.get_text() for text in legend.get_texts()]
        if isinstance(curve_kwargs, list):
            # Multiple labels in legend
            assert len(legend_labels) == 3
            for idx, label in enumerate(legend_labels):
                if name is None:
                    assert label == f"{auc_metric_name} = {auc_metrics[idx]:.2f}"
                elif isinstance(name, str):
                    assert (
                        label == f"single ({auc_metric_name} = {auc_metrics[idx]:.2f})"
                    )
                else:
                    # `name` is a list of different strings
                    assert (
                        label
                        == f"{name[idx]} ({auc_metric_name} = {auc_metrics[idx]:.2f})"
                    )
        else:
            # Single label in legend
            assert len(legend_labels) == 1
            if name is None:
                assert legend_labels[0] == (
                    f"{auc_metric_name} = {np.mean(auc_metrics):.2f} +/- "
                    f"{np.std(auc_metrics):.2f}"
                )
            else:
                # name is single string
                assert legend_labels[0] == (
                    f"single ({auc_metric_name} = {np.mean(auc_metrics):.2f} +/- "
                    f"{np.std(auc_metrics):.2f})"
                )


def _check_from_cv_results_curve_kwargs(display_class, cv_results, X, y, curve_kwargs):
    """Check `curve_kwargs` correctly passed in `from_cv_results`.

    `curve_kwargs` should only set "alpha" kwarg.
    """
    display = display_class.from_cv_results(
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
        assert all(line.get_alpha() == curve_kwargs["alpha"] for line in display.line_)
    else:
        # Different `alpha` used for each curve
        assert all(
            line.get_alpha() == curve_kwargs[i]["alpha"]
            for i, line in enumerate(display.line_)
        )


def _check_display_kwargs_deprecation(
    display_class, constructor_name, clf, X, y, display_args
):
    """Check **kwargs deprecated correctly in favour of `curve_kwargs`."""
    # Error when both `curve_kwargs` and `**kwargs` provided
    with pytest.raises(ValueError, match="Cannot provide both `curve_kwargs`"):
        if constructor_name == "from_estimator":
            display_class.from_estimator(
                clf, X, y, curve_kwargs={"alpha": 1}, label="test"
            )
        elif constructor_name == "from_predictions":
            display_class.from_predictions(
                y, y, curve_kwargs={"alpha": 1}, label="test"
            )
        else:
            display_class(**display_args).plot(curve_kwargs={"alpha": 1}, label="test")

    # Warning when `**kwargs`` provided
    with pytest.warns(FutureWarning, match=r"`\*\*kwargs` is deprecated and will be"):
        if constructor_name == "from_estimator":
            display_class.from_estimator(clf, X, y, label="test")
        elif constructor_name == "from_predictions":
            display_class.from_predictions(y, y, label="test")
        else:
            display_class(**display_args).plot(label="test")


def _check_pos_label_statistics(
    display_class, response_method, constructor_name, check_metric
):
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
        display = display_class.from_estimator(
            classifier,
            X_test,
            y_test,
            pos_label=pos_label,
            response_method=response_method,
        )
    elif constructor_name == "from_predictions":
        display = display_class.from_predictions(
            y_test,
            y_score,
            pos_label=pos_label,
        )
    else:
        display = display_class.from_cv_results(
            cv_results,
            X,
            y,
            response_method=response_method,
            pos_label=pos_label,
        )

    check_metric(display, constructor_name, pos_label)

    pos_label = "not cancer"
    y_score = y_score_not_cancer
    if constructor_name == "from_estimator":
        display = display_class.from_estimator(
            classifier,
            X_test,
            y_test,
            response_method=response_method,
            pos_label=pos_label,
        )
    elif constructor_name == "from_predictions":
        display = display_class.from_predictions(
            y_test,
            y_score,
            pos_label=pos_label,
        )
    else:
        display = display_class.from_cv_results(
            cv_results,
            X,
            y,
            response_method=response_method,
            pos_label=pos_label,
        )

    check_metric(display, constructor_name, pos_label)


@pytest.mark.parametrize(
    "Display, display_kwargs",
    [
        # TODO(1.10): Remove
        (
            PrecisionRecallDisplay,
            {"precision": np.array([1, 0.5, 0]), "recall": np.array([0, 0.5, 1])},
        ),
        # TODO(1.9): Remove
        (RocCurveDisplay, {"fpr": np.array([0, 0.5, 1]), "tpr": np.array([0, 0.5, 1])}),
    ],
)
def test_display_estimator_name_deprecation(pyplot, Display, display_kwargs):
    """Check deprecation of `estimator_name`."""
    with pytest.warns(FutureWarning, match="`estimator_name` is deprecated in"):
        Display(**display_kwargs, estimator_name="test")
