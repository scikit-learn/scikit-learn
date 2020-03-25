import pytest
from numpy.testing import assert_allclose
import numpy as np

from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import plot_roc_curve
from sklearn.datasets import load_iris
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc
from sklearn.base import ClassifierMixin
from sklearn.exceptions import NotFittedError
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.compose import make_column_transformer


# TODO: Remove when https://github.com/numpy/numpy/issues/14397 is resolved
pytestmark = pytest.mark.filterwarnings(
    "ignore:In future, it will be an error for 'np.bool_':DeprecationWarning:"
    "matplotlib.*")


@pytest.fixture(scope="module")
def data():
    return load_iris(return_X_y=True)


@pytest.fixture(scope="module")
def data_binary(data):
    X, y = data
    return X[y < 2], y[y < 2]


def test_plot_roc_curve_error_non_binary(pyplot, data):
    X, y = data
    clf = DecisionTreeClassifier()
    clf.fit(X, y)

    msg = "DecisionTreeClassifier should be a binary classifier"
    with pytest.raises(ValueError, match=msg):
        plot_roc_curve(clf, X, y)


@pytest.mark.parametrize(
    "response_method, msg",
    [("predict_proba", "response method predict_proba is not defined in "
                       "MyClassifier"),
     ("decision_function", "response method decision_function is not defined "
                           "in MyClassifier"),
     ("auto", "response method decision_function or predict_proba is not "
              "defined in MyClassifier"),
     ("bad_method", "response_method must be 'predict_proba', "
                    "'decision_function' or 'auto'")])
def test_plot_roc_curve_error_no_response(pyplot, data_binary, response_method,
                                          msg):
    X, y = data_binary

    class MyClassifier(ClassifierMixin):
        def fit(self, X, y):
            self.classes_ = [0, 1]
            return self

    clf = MyClassifier().fit(X, y)

    with pytest.raises(ValueError, match=msg):
        plot_roc_curve(clf, X, y, response_method=response_method)


@pytest.mark.parametrize("response_method",
                         ["predict_proba", "decision_function"])
@pytest.mark.parametrize("with_sample_weight", [True, False])
@pytest.mark.parametrize("drop_intermediate", [True, False])
@pytest.mark.parametrize("with_strings", [True, False])
def test_plot_roc_curve(pyplot, response_method, data_binary,
                        with_sample_weight, drop_intermediate,
                        with_strings):
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

    viz = plot_roc_curve(lr, X, y, alpha=0.8, sample_weight=sample_weight,
                         drop_intermediate=drop_intermediate)

    y_pred = getattr(lr, response_method)(X)
    if y_pred.ndim == 2:
        y_pred = y_pred[:, 1]

    fpr, tpr, _ = roc_curve(y, y_pred, sample_weight=sample_weight,
                            drop_intermediate=drop_intermediate,
                            pos_label=pos_label)

    assert_allclose(viz.roc_auc, auc(fpr, tpr))
    assert_allclose(viz.fpr, fpr)
    assert_allclose(viz.tpr, tpr)

    assert viz.estimator_name == "LogisticRegression"

    # cannot fail thanks to pyplot fixture
    import matplotlib as mpl  # noqal
    assert isinstance(viz.line_, mpl.lines.Line2D)
    assert viz.line_.get_alpha() == 0.8
    assert isinstance(viz.ax_, mpl.axes.Axes)
    assert isinstance(viz.figure_, mpl.figure.Figure)

    expected_label = "LogisticRegression (AUC = {:0.2f})".format(viz.roc_auc)
    assert viz.line_.get_label() == expected_label
    assert viz.ax_.get_ylabel() == "True Positive Rate"
    assert viz.ax_.get_xlabel() == "False Positive Rate"


@pytest.mark.parametrize(
    "clf", [LogisticRegression(),
            make_pipeline(StandardScaler(), LogisticRegression()),
            make_pipeline(make_column_transformer((StandardScaler(), [0, 1])),
                          LogisticRegression())])
def test_roc_curve_not_fitted_errors(pyplot, data_binary, clf):
    X, y = data_binary
    with pytest.raises(NotFittedError):
        plot_roc_curve(clf, X, y)
    clf.fit(X, y)
    disp = plot_roc_curve(clf, X, y)
    assert clf.__class__.__name__ in disp.line_.get_label()
    assert disp.estimator_name == clf.__class__.__name__


def test_plot_roc_curve_estimator_name_multiple_calls(pyplot, data_binary):
    # non-regression test checking that the `name` used when calling
    # `plot_roc_curve` is used as well when calling `disp.plot()`
    X, y = data_binary
    clf_name = "my hand-crafted name"
    clf = LogisticRegression().fit(X, y)
    disp = plot_roc_curve(clf, X, y, name=clf_name)
    assert disp.estimator_name == clf_name
    pyplot.close("all")
    disp.plot()
    assert clf_name in disp.line_.get_label()
    pyplot.close("all")
    clf_name = "another_name"
    disp.plot(name=clf_name)
    assert clf_name in disp.line_.get_label()
