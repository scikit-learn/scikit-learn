import pytest
from numpy.testing import assert_allclose

from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import plot_roc_curve
from sklearn.datasets import load_iris
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc


@pytest.fixture(scope="module")
def data():
    X, y = load_iris(return_X_y=True)
    return X, y


@pytest.fixture(scope="module")
def data_binary(data):
    X, y = data
    return X[y < 2], y[y < 2]


def test_plot_roc_curve_error_non_binary(data):
    X, y = data
    clf = DecisionTreeClassifier()
    clf.fit(X, y)

    msg = "Estimator must be a binary classifier"
    with pytest.raises(ValueError, match=msg):
        plot_roc_curve(clf, X, y)


@pytest.mark.parametrize("response_method",
                         ["predict_proba", "decision_function", "auto"])
def test_plot_roc_curve_error_no_response(data_binary, response_method):
    X, y = data_binary

    class MyClassifier:
        pass

    clf = MyClassifier()

    with pytest.raises(ValueError):
        plot_roc_curve(clf, X, y, response_method=response_method)


@pytest.mark.parametrize("response_method",
                         ["predict_proba", "decision_function"])
def test_plot_roc_curve(pyplot, response_method, data_binary):
    X, y = data_binary

    lr = LogisticRegression()
    lr.fit(X, y)

    viz = plot_roc_curve(lr, X, y, alpha=0.8)

    y_pred = getattr(lr, response_method)(X)
    if y_pred.ndim == 2:
        y_pred = y_pred[:, 1]

    fpr, tpr, _ = roc_curve(y, y_pred)

    assert_allclose(viz.auc_, auc(fpr, tpr))
    assert_allclose(viz.fpr_, fpr)
    assert_allclose(viz.tpr_, tpr)

    assert viz.estimator_name_ == "LogisticRegression"

    import matplotlib as mpl
    assert isinstance(viz.line_, mpl.lines.Line2D)
    assert viz.line_.get_alpha() == 0.8
    assert isinstance(viz.ax_, mpl.axes.Axes)
    assert isinstance(viz.figure_, mpl.figure.Figure)

    expected_label = "LogisticRegression (AUC = {:0.2f})".format(viz.auc_)
    assert viz.line_.get_label() == expected_label
