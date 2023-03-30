import pytest
import numpy as np
from numpy.testing import assert_allclose

from sklearn.compose import make_column_transformer
from sklearn.datasets import load_iris

from sklearn.datasets import load_breast_cancer
from sklearn.exceptions import NotFittedError
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.metrics import roc_curve
from sklearn.metrics import auc

from sklearn.model_selection import cross_validate, train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.utils import shuffle

from sklearn.metrics import RocCurveDisplay, MultiRocCurveDisplay


@pytest.fixture(scope="module")
def data():
    return load_iris(return_X_y=True)


@pytest.fixture(scope="module")
def data_binary(data):
    X, y = data
    return X[y < 2], y[y < 2]


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
    """Check the overall plotting behaviour."""
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

    y_pred = getattr(lr, response_method)(X)
    y_pred = y_pred if y_pred.ndim == 1 else y_pred[:, 1]

    if constructor_name == "from_estimator":
        display = RocCurveDisplay.from_estimator(
            lr,
            X,
            y,
            sample_weight=sample_weight,
            drop_intermediate=drop_intermediate,
            pos_label=pos_label,
            alpha=0.8,
        )
    else:
        display = RocCurveDisplay.from_predictions(
            y,
            y_pred,
            sample_weight=sample_weight,
            drop_intermediate=drop_intermediate,
            pos_label=pos_label,
            alpha=0.8,
        )

    fpr, tpr, _ = roc_curve(
        y,
        y_pred,
        sample_weight=sample_weight,
        drop_intermediate=drop_intermediate,
        pos_label=pos_label,
    )

    assert_allclose(display.roc_auc, auc(fpr, tpr))
    assert_allclose(display.fpr, fpr)
    assert_allclose(display.tpr, tpr)

    assert display.estimator_name == default_name

    import matplotlib as mpl  # noqal

    assert isinstance(display.line_, mpl.lines.Line2D)
    assert display.line_.get_alpha() == 0.8
    assert isinstance(display.ax_, mpl.axes.Axes)
    assert isinstance(display.figure_, mpl.figure.Figure)

    expected_label = f"{default_name} (AUC = {display.roc_auc:.2f})"
    assert display.line_.get_label() == expected_label

    expected_pos_label = 1 if pos_label is None else pos_label
    expected_ylabel = f"True Positive Rate (Positive label: {expected_pos_label})"
    expected_xlabel = f"False Positive Rate (Positive label: {expected_pos_label})"

    assert display.ax_.get_ylabel() == expected_ylabel
    assert display.ax_.get_xlabel() == expected_xlabel


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
    assert display.estimator_name == name


@pytest.mark.parametrize(
    "roc_auc, estimator_name, expected_label",
    [
        (0.9, None, "AUC = 0.90"),
        (None, "my_est", "my_est"),
        (0.8, "my_est2", "my_est2 (AUC = 0.80)"),
    ],
)
def test_roc_curve_display_default_labels(
    pyplot, roc_auc, estimator_name, expected_label
):
    """Check the default labels used in the display."""
    fpr = np.array([0, 0.5, 1])
    tpr = np.array([0, 0.5, 1])
    disp = RocCurveDisplay(
        fpr=fpr,
        tpr=tpr,
        roc_auc=roc_auc,
        estimator_name=estimator_name,
    ).plot()
    assert disp.line_.get_label() == expected_label


@pytest.mark.parametrize("response_method", ["predict_proba", "decision_function"])
@pytest.mark.parametrize("constructor_name", ["from_estimator", "from_predictions"])
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

    # sanity check to be sure the positive class is classes_[0] and that we
    # are betrayed by the class imbalance
    assert classifier.classes_.tolist() == ["cancer", "not cancer"]

    y_pred = getattr(classifier, response_method)(X_test)
    # we select the corresponding probability columns or reverse the decision
    # function otherwise
    y_pred_cancer = -1 * y_pred if y_pred.ndim == 1 else y_pred[:, 0]
    y_pred_not_cancer = y_pred if y_pred.ndim == 1 else y_pred[:, 1]

    if constructor_name == "from_estimator":
        display = RocCurveDisplay.from_estimator(
            classifier,
            X_test,
            y_test,
            pos_label="cancer",
            response_method=response_method,
        )
    else:
        display = RocCurveDisplay.from_predictions(
            y_test,
            y_pred_cancer,
            pos_label="cancer",
        )

    roc_auc_limit = 0.95679

    assert display.roc_auc == pytest.approx(roc_auc_limit)
    assert np.trapz(display.tpr, display.fpr) == pytest.approx(roc_auc_limit)

    if constructor_name == "from_estimator":
        display = RocCurveDisplay.from_estimator(
            classifier,
            X_test,
            y_test,
            response_method=response_method,
            pos_label="not cancer",
        )
    else:
        display = RocCurveDisplay.from_predictions(
            y_test,
            y_pred_not_cancer,
            pos_label="not cancer",
        )

    assert display.roc_auc == pytest.approx(roc_auc_limit)
    assert np.trapz(display.tpr, display.fpr) == pytest.approx(roc_auc_limit)


@pytest.mark.parametrize(
    "param",
    [
        {"return_estimator": False, "return_indices": True},
        {"return_estimator": True, "return_indices": False},
    ],
)
def test_from_cv_results_missing_dict_key(pyplot, data_binary, param):
    """Check that we raise an error if `estimator` or `indices` are missing in
    the dictionary returned by `cross_validate`."""
    X, y = data_binary
    model = make_pipeline(StandardScaler(), LogisticRegression())
    cv_results = cross_validate(model, X, y, cv=3, **param)

    err_msg = "cv_results does not contain one of the following required keys"
    with pytest.raises(ValueError, match=err_msg):
        RocCurveDisplay.from_cv_results(cv_results, X, y)


def test_from_cv_results_inconsistent_n_samples(pyplot, data_binary):
    """Check that we raise an error if the data given in `cross_validate` and
    `from_cv_results` are inconsistent regarding the number of samples."""
    X, y = data_binary
    model = make_pipeline(StandardScaler(), LogisticRegression())
    cv_results = cross_validate(
        model, X, y, cv=3, return_estimator=True, return_indices=True
    )

    err_msg = "X does not contain the correct number of samples"
    with pytest.raises(ValueError, match=err_msg):
        RocCurveDisplay.from_cv_results(cv_results, X[:-1], y[:-1])


def test_from_cv_results_wrong_estimator_type(pyplot, data_binary):
    """Check that we raise an error if the type of estimator in `cv_results` is not
    a binary classifier."""
    X, y = data_binary
    model = make_pipeline(StandardScaler(), LinearRegression())
    cv_results = cross_validate(
        model, X, y, cv=3, return_estimator=True, return_indices=True
    )

    err_msg = "Expected 'estimator' to be a binary classifier."
    with pytest.raises(ValueError, match=err_msg):
        RocCurveDisplay.from_cv_results(cv_results, X, y)


def test_from_cv_results_wrong_length_fold_name(pyplot, data_binary):
    """Check that we raise an error if the length of the fold name is not
    consistent with the number of estimators in `cv_results`."""
    X, y = data_binary
    model = make_pipeline(StandardScaler(), LogisticRegression())

    n_fold = 3
    cv_results = cross_validate(
        model, X, y, cv=3, return_estimator=True, return_indices=True
    )

    fold_name = [f"fold_{i}" for i in range(n_fold + 1)]
    err_msg = "When `fold_name` is provided, it must have the same length as "
    with pytest.raises(ValueError, match=err_msg):
        RocCurveDisplay.from_cv_results(cv_results, X, y, fold_name=fold_name)

    # this error could be also raised when invoking `plot` using the display
    display = RocCurveDisplay.from_cv_results(cv_results, X, y, kind="folds")
    with pytest.raises(ValueError, match=err_msg):
        display.plot(fold_name=fold_name)


def test_from_cv_results_wrong_kind(pyplot, data_binary):
    """Check that we raise an error if `kind` is unknown."""
    X, y = data_binary
    model = make_pipeline(StandardScaler(), LogisticRegression())

    cv_results = cross_validate(
        model, X, y, cv=3, return_estimator=True, return_indices=True
    )

    err_msg = "Parameter `kind` must be one of"
    with pytest.raises(ValueError, match=err_msg):
        RocCurveDisplay.from_cv_results(cv_results, X, y, kind="unknown")


def test_from_cv_results_wrong_length_fold_line_kw(pyplot, data_binary):
    """Check that we raise an error if the length of `fold_line_kw` is not
    consistent with the number of estimators in `cv_results`."""
    X, y = data_binary
    model = make_pipeline(StandardScaler(), LogisticRegression())

    n_fold = 3
    cv_results = cross_validate(
        model, X, y, cv=3, return_estimator=True, return_indices=True
    )

    fold_line_kw = [{}] * (n_fold + 1)
    err_msg = "When `fold_line_kw` is a list, it must have the same length as "
    with pytest.raises(ValueError, match=err_msg):
        RocCurveDisplay.from_cv_results(
            cv_results, X, y, kind="folds", fold_line_kw=fold_line_kw
        )


def test_from_cv_results_default(pyplot, data_binary):
    """Check default behaviour of `RocCurveDisplay.from_cv_results`."""
    X, y = data_binary
    model = make_pipeline(StandardScaler(), LogisticRegression())
    n_folds = 3
    cv_results = cross_validate(
        model, X, y, cv=n_folds, return_estimator=True, return_indices=True
    )
    display = RocCurveDisplay.from_cv_results(cv_results, X, y, kind="folds")
    assert isinstance(display, MultiRocCurveDisplay)
    assert all(isinstance(d, RocCurveDisplay) for d in display.displays)

    legend = display.ax_.get_legend()
    assert not legend.get_title().get_text()
    assert len(legend.get_texts()) == n_folds
    for i, text in enumerate(legend.get_texts()):
        assert text.get_text().startswith(f"ROC fold #{i}")

    display = RocCurveDisplay.from_cv_results(cv_results, X, y, kind="both")

    legend = display.ax_.get_legend()
    assert legend.get_title().get_text() == "Uncertainties via cross-validation"
    assert len(legend.get_texts()) == n_folds + 2  # adding mean and std dev
    for i in range(n_folds):
        assert legend.get_texts()[i].get_text().startswith(f"ROC fold #{i}")
    assert legend.get_texts()[n_folds].get_text().startswith("Mean ROC")
    assert "1 std. dev." in legend.get_texts()[n_folds + 1].get_text()


def test_from_cv_results_attributes(pyplot, data_binary):
    """Check that `MultiRocCurveDisplay` contains the correct attributes with
    the expected behaviour."""
    X, y = data_binary
    model = make_pipeline(StandardScaler(), LogisticRegression())
    n_folds = 3

    cv_results = cross_validate(
        model, X, y, cv=n_folds, return_estimator=True, return_indices=True
    )
    display = RocCurveDisplay.from_cv_results(cv_results, X, y, kind="folds")

    import matplotlib as mpl  # noqal

    assert isinstance(display.ax_, mpl.axes.Axes)
    assert isinstance(display.figure_, mpl.figure.Figure)

    assert isinstance(display.fold_lines_, np.ndarray)
    assert len(display.fold_lines_) == n_folds
    assert all(isinstance(line, mpl.lines.Line2D) for line in display.fold_lines_)
    assert display.mean_line_ is None
    assert display.std_area_ is None
    assert display.chance_level_ is None

    display = RocCurveDisplay.from_cv_results(cv_results, X, y, kind="aggregate")

    assert isinstance(display.ax_, mpl.axes.Axes)
    assert isinstance(display.figure_, mpl.figure.Figure)

    assert display.fold_lines_ is None
    assert isinstance(display.mean_line_, mpl.lines.Line2D)
    assert isinstance(display.std_area_, mpl.collections.PolyCollection)
    assert display.chance_level_ is None

    display = RocCurveDisplay.from_cv_results(
        cv_results, X, y, kind="both", plot_chance_level=True
    )

    assert isinstance(display.ax_, mpl.axes.Axes)
    assert isinstance(display.figure_, mpl.figure.Figure)

    assert isinstance(display.fold_lines_, np.ndarray)
    assert len(display.fold_lines_) == n_folds
    assert all(isinstance(line, mpl.lines.Line2D) for line in display.fold_lines_)
    assert isinstance(display.mean_line_, mpl.lines.Line2D)
    assert isinstance(display.std_area_, mpl.collections.PolyCollection)
    assert isinstance(display.chance_level_, mpl.lines.Line2D)


def test_from_cv_results_names(pyplot, data_binary):
    """Check that passing names for each ROC curves is behaving as expected."""
    X, y = data_binary
    model = make_pipeline(StandardScaler(), LogisticRegression())
    n_folds = 3

    cv_results = cross_validate(
        model, X, y, cv=n_folds, return_estimator=True, return_indices=True
    )
    fold_name = [f"From CV #{i}" for i in range(n_folds)]
    aggregate_name = "Mean from CV"
    display = RocCurveDisplay.from_cv_results(
        cv_results,
        X,
        y,
        kind="both",
        fold_name=fold_name,
        aggregate_name=aggregate_name,
    )

    legend = display.ax_.get_legend()
    for i in range(n_folds):
        assert fold_name[i] in legend.get_texts()[i].get_text()
    assert aggregate_name in legend.get_texts()[n_folds].get_text()

    # passing name in `plot` should override the name passed in `from_cv_results`
    fold_name = [f"From plot #{i}" for i in range(n_folds)]
    aggregate_name = "Mean from plot"
    display.plot(kind="both", fold_name=fold_name, aggregate_name=aggregate_name)

    legend = display.ax_.get_legend()
    for i in range(n_folds):
        assert fold_name[i] in legend.get_texts()[i].get_text()
    assert aggregate_name in legend.get_texts()[n_folds].get_text()


def test_from_cv_results_effect_kwargs(pyplot, data_binary):
    """Check that the different keywords arguments are behaving as expected."""
    X, y = data_binary
    model = make_pipeline(StandardScaler(), LogisticRegression())

    cv_results = cross_validate(
        model, X, y, cv=3, return_estimator=True, return_indices=True
    )

    fold_line_kw = {"color": "red", "linestyle": "--", "linewidth": 2}
    display = RocCurveDisplay.from_cv_results(
        cv_results, X, y, fold_line_kw=fold_line_kw
    )

    for line in display.fold_lines_:
        assert line.get_color() == fold_line_kw["color"]
        assert line.get_linestyle() == fold_line_kw["linestyle"]
        assert line.get_linewidth() == fold_line_kw["linewidth"]

    fold_line_kw = [{"color": "red"}, {"color": "green"}, {"color": "blue"}]
    display = RocCurveDisplay.from_cv_results(
        cv_results, X, y, fold_line_kw=fold_line_kw
    )

    for line, kwargs in zip(display.fold_lines_, fold_line_kw):
        assert line.get_color() == kwargs["color"]

    aggregate_line_kw = {"color": "red", "linestyle": "--", "linewidth": 2}
    aggregate_uncertainty_kw = {"alpha": 0.5, "color": "red"}
    display = RocCurveDisplay.from_cv_results(
        cv_results,
        X,
        y,
        aggregate_line_kw=aggregate_line_kw,
        aggregate_uncertainty_kw=aggregate_uncertainty_kw,
        kind="both",
    )

    assert display.mean_line_.get_color() == aggregate_line_kw["color"]
    assert display.mean_line_.get_linestyle() == aggregate_line_kw["linestyle"]
    assert display.mean_line_.get_linewidth() == aggregate_line_kw["linewidth"]

    # `facecolor` is encoded in RGBA
    assert_allclose(
        display.std_area_.get_facecolor(),
        [[1, 0, 0, aggregate_uncertainty_kw["alpha"]]],
    )

    chance_level_kw = {"color": "red", "label": "random"}
    display = RocCurveDisplay.from_cv_results(
        cv_results, X, y, plot_chance_level=True, chance_level_kw=chance_level_kw
    )

    assert display.chance_level_.get_color() == chance_level_kw["color"]
    assert display.chance_level_.get_label() == chance_level_kw["label"]
