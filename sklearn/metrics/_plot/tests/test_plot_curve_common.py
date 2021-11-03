import numpy as np
import pytest

from sklearn.base import ClassifierMixin
from sklearn.base import clone
from sklearn.compose import make_column_transformer
from sklearn.datasets import load_iris
from sklearn.exceptions import NotFittedError
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor

from sklearn.metrics import plot_det_curve, plot_roc_curve, plot_precision_recall_curve

pytestmark = pytest.mark.filterwarnings(
    "ignore:Function plot_roc_curve is deprecated",
    "ignore:Function plot_det_curve is deprecated",
    "ignore:Function plot_precision_recall_curve is deprecated",
)


@pytest.fixture(scope="module")
def data():
    return load_iris(return_X_y=True)


@pytest.fixture(scope="module")
def data_binary(data):
    X, y = data
    return X[y < 2], y[y < 2]


@pytest.mark.parametrize(
    "plot_func", [plot_det_curve, plot_roc_curve, plot_precision_recall_curve]
)
def test_plot_curve_error_classifier(pyplot, data, data_binary, plot_func):
    """Check that a proper error is raised when only binary classification is
    supported."""
    X, y = data
    X_binary, y_binary = data_binary

    # Case 1: multiclass classifier with multiclass target
    clf = DecisionTreeClassifier().fit(X, y)

    msg = (
        "This DecisionTreeClassifier instance is not a binary classifier. It was "
        f"fitted on multiclass problem with {len(np.unique(y))} classes."
    )
    with pytest.raises(ValueError, match=msg):
        plot_func(clf, X, y)

    # Case 2: multiclass classifier with binary target
    with pytest.raises(ValueError, match=msg):
        plot_func(clf, X_binary, y_binary)

    # Case 3: binary classifier with multiclass target
    clf = DecisionTreeClassifier().fit(X_binary, y_binary)
    msg = "The target y is not binary. Got multiclass type of target."
    with pytest.raises(ValueError, match=msg):
        plot_func(clf, X, y)


@pytest.mark.parametrize(
    "plot_func", [plot_det_curve, plot_roc_curve, plot_precision_recall_curve]
)
def test_plot_curve_error_regression(pyplot, data_binary, plot_func):
    """Check that we raise an error with regressor."""

    # Case 1: regressor
    X, y = data_binary
    regressor = DecisionTreeRegressor().fit(X, y)

    msg = (
        "This plotting functionalities only support a binary classifier. Got a "
        "DecisionTreeRegressor instead."
    )
    with pytest.raises(ValueError, match=msg):
        plot_func(regressor, X, y)

    # Case 2: regression target
    classifier = DecisionTreeClassifier().fit(X, y)
    # Force `y_true` to be seen as a regression problem
    y = y + 0.5
    msg = "The target y is not binary. Got continuous type of target."
    with pytest.raises(ValueError, match=msg):
        plot_func(classifier, X, y)


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
            "MyClassifier has none of the following attributes: predict_proba, "
            "decision_function.",
        ),
        (
            "bad_method",
            "MyClassifier has none of the following attributes: bad_method",
        ),
    ],
)
@pytest.mark.parametrize(
    "plot_func", [plot_det_curve, plot_roc_curve, plot_precision_recall_curve]
)
def test_plot_curve_error_no_response(
    pyplot,
    data_binary,
    response_method,
    msg,
    plot_func,
):
    X, y = data_binary

    class MyClassifier(ClassifierMixin):
        def fit(self, X, y):
            self.classes_ = [0, 1]
            return self

    clf = MyClassifier().fit(X, y)

    with pytest.raises(AttributeError, match=msg):
        plot_func(clf, X, y, response_method=response_method)


@pytest.mark.parametrize(
    "plot_func", [plot_det_curve, plot_roc_curve, plot_precision_recall_curve]
)
def test_plot_curve_estimator_name_multiple_calls(pyplot, data_binary, plot_func):
    # non-regression test checking that the `name` used when calling
    # `plot_func` is used as well when calling `disp.plot()`
    X, y = data_binary
    clf_name = "my hand-crafted name"
    clf = LogisticRegression().fit(X, y)
    disp = plot_func(clf, X, y, name=clf_name)
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
    "plot_func", [plot_det_curve, plot_roc_curve, plot_precision_recall_curve]
)
def test_plot_det_curve_not_fitted_errors(pyplot, data_binary, clf, plot_func):
    X, y = data_binary
    # clone since we parametrize the test and the classifier will be fitted
    # when testing the second and subsequent plotting function
    model = clone(clf)
    with pytest.raises(NotFittedError):
        plot_func(model, X, y)
    model.fit(X, y)
    disp = plot_func(model, X, y)
    assert model.__class__.__name__ in disp.line_.get_label()
    assert disp.estimator_name == model.__class__.__name__
