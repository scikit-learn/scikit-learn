from numpy.testing import (
    assert_allclose,
    assert_array_equal,
)
import pytest

from sklearn.datasets import make_classification
from sklearn.svm import SVC
from sklearn.svm import SVR

from sklearn.metrics import ConfusionMatrixDisplay
from sklearn.metrics import confusion_matrix


# TODO: Remove when https://github.com/numpy/numpy/issues/14397 is resolved
pytestmark = pytest.mark.filterwarnings(
    "ignore:In future, it will be an error for 'np.bool_':DeprecationWarning:"
    "matplotlib.*"
)


@pytest.fixture(scope="module")
def n_classes():
    return 5


@pytest.fixture(scope="module")
def data(n_classes):
    X, y = make_classification(
        n_samples=100, n_informative=5, n_classes=n_classes, random_state=0
    )
    return X, y


@pytest.fixture(scope="module")
def fitted_clf(data):
    return SVC(kernel="linear", C=0.01).fit(*data)


@pytest.fixture(scope="module")
def y_pred(data, fitted_clf):
    X, _ = data
    return fitted_clf.predict(X)


@pytest.fixture(scope="module")
def estimator_api_params(fitted_clf, data):
    X, y = data
    return ConfusionMatrixDisplay.from_estimator, fitted_clf, X, y


@pytest.fixture(scope="module")
def predictions_api_params(data, y_pred):
    _, y = data
    return ConfusionMatrixDisplay.from_predictions, y, y_pred


@pytest.fixture(
    params=["estimator_api_params", "predictions_api_params"],
    scope="module",
)
def confusion_matrix_display_fxt(request):
    return request.getfixturevalue(request.param)


def test_confusion_matrix_display_error_on_regressor(pyplot, data):
    X, y = data
    regressor = SVR().fit(X, y)

    err_msg = "ConfusionMatrixDisplay.from_estimator only supports classifiers"
    with pytest.raises(ValueError, match=err_msg):
        ConfusionMatrixDisplay.from_estimator(regressor, X, y)

    y_pred = regressor.predict(X)
    err_msg = "Got `y_pred` of type 'continuous'"
    with pytest.raises(ValueError, match=err_msg):
        ConfusionMatrixDisplay.from_predictions(y, y_pred)


def test_confusion_matrix_display_invalid_option(
    pyplot, confusion_matrix_display_fxt
):
    constructor, *params = confusion_matrix_display_fxt
    extra_params = {"normalize": "invalid"}

    err_msg = r"normalize must be one of \{'true', 'pred', 'all', None\}"
    with pytest.raises(ValueError, match=err_msg):
        constructor(*params, **extra_params)


@pytest.mark.parametrize("with_labels", [True, False])
@pytest.mark.parametrize("with_display_labels", [True, False])
def test_confusion_matrix_display_custom_labels(
    pyplot,
    confusion_matrix_display_fxt,
    data,
    y_pred,
    n_classes,
    with_labels,
    with_display_labels,
):

    X, y = data
    constructor, *params = confusion_matrix_display_fxt

    ax = pyplot.gca()
    labels = [2, 1, 0, 3, 4] if with_labels else None
    display_labels = ["b", "d", "a", "e", "f"] if with_display_labels else None

    cm = confusion_matrix(y, y_pred, labels=labels)
    disp = constructor(
        *params,
        ax=ax,
        display_labels=display_labels,
        labels=labels
    )
    assert_allclose(disp.confusion_matrix, cm)

    if with_display_labels:
        expected_display_labels = display_labels
    elif with_labels:
        expected_display_labels = labels
    else:
        expected_display_labels = list(range(n_classes))

    expected_display_labels_str = [str(name)
                                   for name in expected_display_labels]

    x_ticks = [tick.get_text() for tick in disp.ax_.get_xticklabels()]
    y_ticks = [tick.get_text() for tick in disp.ax_.get_yticklabels()]

    assert_array_equal(disp.display_labels, expected_display_labels)
    assert_array_equal(x_ticks, expected_display_labels_str)
    assert_array_equal(y_ticks, expected_display_labels_str)
