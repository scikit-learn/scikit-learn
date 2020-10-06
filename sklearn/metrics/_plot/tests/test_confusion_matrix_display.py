from numpy.testing import (
    assert_allclose,
    assert_array_equal,
)
import numpy as np
import pytest

from sklearn.datasets import make_classification
from sklearn.compose import make_column_transformer
from sklearn.exceptions import NotFittedError
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
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


def test_confusion_matrix_display_error_on_regressor(pyplot, data, y_pred):
    X, y = data
    regressor = SVR().fit(X, y)

    err_msg = "ConfusionMatrixDisplay.from_estimator only supports classifiers"
    with pytest.raises(ValueError, match=err_msg):
        ConfusionMatrixDisplay.from_estimator(regressor, X, y)

    y_pred_regressor = regressor.predict(X)

    err_msg = (
        "Only 'binary' and 'multiclass' classification problems are supported"
    )
    with pytest.raises(ValueError, match=err_msg):
        # Force `y_true` to be seen as a regression problem
        ConfusionMatrixDisplay.from_predictions(y + 0.5, y_pred)

    err_msg = "Got `y_pred` of type 'continuous'"
    with pytest.raises(ValueError, match=err_msg):
        ConfusionMatrixDisplay.from_predictions(y, y_pred_regressor)

    err_msg = "Found input variables with inconsistent numbers of samples"
    with pytest.raises(ValueError, match=err_msg):
        ConfusionMatrixDisplay.from_predictions(y, y_pred[::2])


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


@pytest.mark.parametrize("normalize", ["true", "pred", "all", None])
@pytest.mark.parametrize("include_values", [True, False])
def test_confusion_matrix_display_plotting(
    pyplot,
    confusion_matrix_display_fxt,
    data,
    y_pred,
    n_classes,
    normalize,
    include_values,
):
    X, y = data
    constructor, *params = confusion_matrix_display_fxt
    ax = pyplot.gca()
    cmap = "plasma"

    cm = confusion_matrix(y, y_pred)
    disp = constructor(
        *params,
        normalize=normalize,
        cmap=cmap,
        ax=ax,
        include_values=include_values
    )

    assert disp.ax_ == ax

    if normalize == "true":
        cm = cm / cm.sum(axis=1, keepdims=True)
    elif normalize == "pred":
        cm = cm / cm.sum(axis=0, keepdims=True)
    elif normalize == "all":
        cm = cm / cm.sum()

    assert_allclose(disp.confusion_matrix, cm)
    import matplotlib as mpl

    assert isinstance(disp.im_, mpl.image.AxesImage)
    assert disp.im_.get_cmap().name == cmap
    assert isinstance(disp.ax_, pyplot.Axes)
    assert isinstance(disp.figure_, pyplot.Figure)

    assert disp.ax_.get_ylabel() == "True label"
    assert disp.ax_.get_xlabel() == "Predicted label"

    x_ticks = [tick.get_text() for tick in disp.ax_.get_xticklabels()]
    y_ticks = [tick.get_text() for tick in disp.ax_.get_yticklabels()]

    expected_display_labels = list(range(n_classes))

    expected_display_labels_str = [
        str(name) for name in expected_display_labels
    ]

    assert_array_equal(disp.display_labels, expected_display_labels)
    assert_array_equal(x_ticks, expected_display_labels_str)
    assert_array_equal(y_ticks, expected_display_labels_str)

    image_data = disp.im_.get_array().data
    assert_allclose(image_data, cm)

    if include_values:
        assert disp.text_.shape == (n_classes, n_classes)
        fmt = ".2g"
        expected_text = np.array([format(v, fmt) for v in cm.ravel(order="C")])
        text_text = np.array(
            [t.get_text() for t in disp.text_.ravel(order="C")]
        )
        assert_array_equal(expected_text, text_text)
    else:
        assert disp.text_ is None


def test_confusion_matrix_display(
    pyplot, confusion_matrix_display_fxt, data, y_pred, n_classes
):
    X, y = data
    constructor, *params = confusion_matrix_display_fxt

    cm = confusion_matrix(y, y_pred)
    disp = constructor(
        *params,
        normalize=None,
        include_values=True,
        cmap="viridis",
        xticks_rotation=45.0,
    )

    assert_allclose(disp.confusion_matrix, cm)
    assert disp.text_.shape == (n_classes, n_classes)

    rotations = [tick.get_rotation() for tick in disp.ax_.get_xticklabels()]
    assert_allclose(rotations, 45.0)

    image_data = disp.im_.get_array().data
    assert_allclose(image_data, cm)

    disp.plot(cmap="plasma")
    assert disp.im_.get_cmap().name == "plasma"

    disp.plot(include_values=False)
    assert disp.text_ is None

    disp.plot(xticks_rotation=90.0)
    rotations = [tick.get_rotation() for tick in disp.ax_.get_xticklabels()]
    assert_allclose(rotations, 90.0)

    disp.plot(values_format="e")
    expected_text = np.array([format(v, "e") for v in cm.ravel(order="C")])
    text_text = np.array([t.get_text() for t in disp.text_.ravel(order="C")])
    assert_array_equal(expected_text, text_text)


def test_confusion_matrix_contrast(pyplot):
    # make sure text color is appropriate depending on background

    cm = np.eye(2) / 2
    disp = ConfusionMatrixDisplay(cm, display_labels=[0, 1])

    disp.plot(cmap=pyplot.cm.gray)
    # diagonal text is black
    assert_allclose(disp.text_[0, 0].get_color(), [0.0, 0.0, 0.0, 1.0])
    assert_allclose(disp.text_[1, 1].get_color(), [0.0, 0.0, 0.0, 1.0])

    # off-diagonal text is white
    assert_allclose(disp.text_[0, 1].get_color(), [1.0, 1.0, 1.0, 1.0])
    assert_allclose(disp.text_[1, 0].get_color(), [1.0, 1.0, 1.0, 1.0])

    disp.plot(cmap=pyplot.cm.gray_r)
    # diagonal text is white
    assert_allclose(disp.text_[0, 1].get_color(), [0.0, 0.0, 0.0, 1.0])
    assert_allclose(disp.text_[1, 0].get_color(), [0.0, 0.0, 0.0, 1.0])

    # off-diagonal text is black
    assert_allclose(disp.text_[0, 0].get_color(), [1.0, 1.0, 1.0, 1.0])
    assert_allclose(disp.text_[1, 1].get_color(), [1.0, 1.0, 1.0, 1.0])

    # Regression test for #15920
    cm = np.array([[19, 34], [32, 58]])
    disp = ConfusionMatrixDisplay(cm, display_labels=[0, 1])

    disp.plot(cmap=pyplot.cm.Blues)
    min_color = pyplot.cm.Blues(0)
    max_color = pyplot.cm.Blues(255)
    assert_allclose(disp.text_[0, 0].get_color(), max_color)
    assert_allclose(disp.text_[0, 1].get_color(), max_color)
    assert_allclose(disp.text_[1, 0].get_color(), max_color)
    assert_allclose(disp.text_[1, 1].get_color(), min_color)


@pytest.mark.parametrize(
    "clf",
    [
        LogisticRegression(),
        make_pipeline(StandardScaler(), LogisticRegression()),
        make_pipeline(
            make_column_transformer((StandardScaler(), [0, 1])),
            LogisticRegression(),
        ),
    ],
    ids=["clf", "pipeline-clf", "pipeline-column_transformer-clf"]
)
def test_confusion_matrix_pipeline(pyplot, clf, data, n_classes):
    X, y = data
    with pytest.raises(NotFittedError):
        ConfusionMatrixDisplay.from_estimator(clf, X, y)
    clf.fit(X, y)
    y_pred = clf.predict(X)

    disp = ConfusionMatrixDisplay.from_estimator(clf, X, y)
    cm = confusion_matrix(y, y_pred)

    assert_allclose(disp.confusion_matrix, cm)
    assert disp.text_.shape == (n_classes, n_classes)
