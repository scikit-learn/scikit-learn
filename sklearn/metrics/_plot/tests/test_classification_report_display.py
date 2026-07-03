# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
import pytest
from numpy.testing import assert_allclose

from sklearn.datasets import make_classification
from sklearn.exceptions import NotFittedError
from sklearn.metrics import ClassificationReportDisplay, classification_report
from sklearn.svm import SVC, SVR

# Columns that carry a 0-1 score and share the colormap.
SCORE_COLUMNS = ["precision", "recall", "f1-score"]


@pytest.fixture
def data():
    X, y = make_classification(
        n_samples=200, n_informative=5, n_classes=3, random_state=0
    )
    return X, y


@pytest.fixture
def fitted_clf(data):
    X, y = data
    return SVC().fit(X, y)


def test_classification_report_display_validation(pyplot, data):
    """from_estimator rejects unfitted estimators and regressors."""
    X, y = data

    with pytest.raises(NotFittedError):
        ClassificationReportDisplay.from_estimator(SVC(), X, y)

    regressor = SVR().fit(X, y)
    err_msg = "ClassificationReportDisplay.from_estimator only supports classifiers"
    with pytest.raises(ValueError, match=err_msg):
        ClassificationReportDisplay.from_estimator(regressor, X, y)


def test_from_dict_matches_report(pyplot, data):
    """Constructing from an output_dict renders exactly those score values."""
    X, y = data
    y_pred = SVC().fit(X, y).predict(X)
    report = classification_report(y, y_pred, output_dict=True)

    disp = ClassificationReportDisplay(report).plot()

    # rows: each class + macro avg + weighted avg (accuracy scalar excluded)
    expected_rows = [k for k in report if k != "accuracy"]
    assert list(disp.display_labels) == expected_rows

    # score cells match report; im_ covers only the 3 score columns
    assert disp.im_.get_array().shape == (len(expected_rows), len(SCORE_COLUMNS))
    for i, row in enumerate(expected_rows):
        for j, col in enumerate(SCORE_COLUMNS):
            assert_allclose(disp.im_.get_array()[i, j], report[row][col])


def test_from_predictions_matches_from_dict(pyplot, data):
    """from_predictions produces the same score matrix as the dict constructor."""
    X, y = data
    y_pred = SVC().fit(X, y).predict(X)
    report = classification_report(y, y_pred, output_dict=True)

    disp_pred = ClassificationReportDisplay.from_predictions(y, y_pred)
    disp_dict = ClassificationReportDisplay(report).plot()

    assert_allclose(disp_pred.im_.get_array(), disp_dict.im_.get_array())


def test_from_estimator_matches_from_predictions(pyplot, data, fitted_clf):
    """from_estimator agrees with from_predictions on the classifier's output."""
    X, y = data
    y_pred = fitted_clf.predict(X)

    disp_est = ClassificationReportDisplay.from_estimator(fitted_clf, X, y)
    disp_pred = ClassificationReportDisplay.from_predictions(y, y_pred)

    assert_allclose(disp_est.im_.get_array(), disp_pred.im_.get_array())


def test_color_scale_fixed_to_unit_interval(pyplot, data):
    """Colors reflect scores (0-1), never support magnitude (issue #16880)."""
    X, y = make_classification(
        n_samples=500, weights=[0.9, 0.1], n_classes=2, random_state=0
    )
    y_pred = SVC().fit(X, y).predict(X)

    disp = ClassificationReportDisplay.from_predictions(y, y_pred)

    assert disp.im_.get_clim() == (0.0, 1.0)


def test_target_names_used_as_row_labels(pyplot, data):
    """target_names propagate to the y tick labels."""
    X, y = data
    y_pred = SVC().fit(X, y).predict(X)
    names = ["setosa", "versicolor", "virginica"]

    disp = ClassificationReportDisplay.from_predictions(y, y_pred, target_names=names)

    ytick_labels = [t.get_text() for t in disp.ax_.get_yticklabels()]
    for name in names:
        assert name in ytick_labels


def test_include_values_toggles_text_and_support(pyplot, data):
    """include_values controls annotations; support is a text-only 4th column."""
    X, y = data
    y_pred = SVC().fit(X, y).predict(X)
    report = classification_report(y, y_pred, output_dict=True)
    n_rows = len([k for k in report if k != "accuracy"])

    disp = ClassificationReportDisplay.from_predictions(y, y_pred)
    # 3 score columns + 1 support column
    assert disp.text_.shape == (n_rows, len(SCORE_COLUMNS) + 1)
    # support column holds integer counts, not color-mapped values
    xtick_labels = [t.get_text() for t in disp.ax_.get_xticklabels()]
    assert "support" in xtick_labels

    disp_no_values = ClassificationReportDisplay.from_predictions(
        y, y_pred, include_values=False
    )
    assert disp_no_values.text_ is None


def test_colorbar_toggle(pyplot, data):
    """colorbar=False adds no colorbar."""
    X, y = data
    y_pred = SVC().fit(X, y).predict(X)

    disp = ClassificationReportDisplay.from_predictions(y, y_pred, colorbar=False)
    assert disp.figure_.axes == [disp.ax_]

    disp_cbar = ClassificationReportDisplay.from_predictions(y, y_pred, colorbar=True)
    assert len(disp_cbar.figure_.axes) == 2


def test_plot_returns_display_with_attributes(pyplot, data):
    """plot() stores the standard Display attributes."""
    X, y = data
    y_pred = SVC().fit(X, y).predict(X)
    report = classification_report(y, y_pred, output_dict=True)

    disp = ClassificationReportDisplay(report).plot()

    assert disp.figure_ is not None
    assert disp.ax_ is not None
    assert disp.im_ is not None
    assert isinstance(disp.display_labels, (list, np.ndarray))
