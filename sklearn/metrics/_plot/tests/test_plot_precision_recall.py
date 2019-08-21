import pytest
from numpy.testing import assert_allclose

from sklearn.metrics import plot_precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve
from sklearn.datasets import load_iris
from sklearn.tree import DecisionTreeClassifier
from sklearn.linear_model import LogisticRegression


@pytest.fixture(scope="module")
def data():
    return load_iris(return_X_y=True)


@pytest.fixture(scope="module")
def data_binary(data):
    X, y = data
    return X[y < 2], y[y < 2]


def test_error_non_binary(data):
    X, y = data
    clf = DecisionTreeClassifier()
    clf.fit(X, y)

    msg = "Estimator should solve a binary classification problem"
    with pytest.raises(ValueError, match=msg):
        plot_precision_recall_curve(clf, X, y)


@pytest.mark.parametrize(
    "response_method, msg",
    [("predict_proba", "response method predict_proba is not defined"),
     ("decision_function", "response method decision_function is not defined"),
     ("auto", "response methods not defined"),
     ("bad_method", "response_method must be 'predict_proba', "
                    "'decision_function' or 'auto'")])
def test_error_no_response(data_binary, response_method, msg):
    X, y = data_binary

    class MyClassifier:
        pass

    clf = MyClassifier()

    with pytest.raises(ValueError, match=msg):
        plot_precision_recall_curve(clf, X, y, response_method=response_method)


@pytest.mark.parametrize("response_method",
                         ["predict_proba", "decision_function"])
def test_plot_roc_curve(pyplot, response_method, data_binary):
    X, y = data_binary

    lr = LogisticRegression()
    lr.fit(X, y)

    viz = plot_precision_recall_curve(lr, X, y, alpha=0.8)

    y_pred = getattr(lr, response_method)(X)
    if y_pred.ndim == 2:
        y_pred = y_pred[:, 1]

    fpr, tpr, _ = precision_recall_curve(y, y_pred)

    assert_allclose(viz.roc_auc, average_precision_score(fpr, tpr))
    assert_allclose(viz.fpr, fpr)
    assert_allclose(viz.tpr, tpr)

    assert viz.estimator_name == "LogisticRegression"

    # cannot fail thanks to pyplot fixture
    import matplotlib as mpl  # noqal
    assert isinstance(viz.line_, mpl.lines.Line2D)
    assert viz.line_.get_alpha() == 0.8
    assert isinstance(viz.ax_, mpl.axes.Axes)
    assert isinstance(viz.figure_, mpl.figure.Figure)

    expected_label = "LogisticRegression (MAP = {:0.2f})".format(viz.roc_auc)
    assert viz.line_.get_label() == expected_label
    assert viz.ax_.get_xlabel() == "Recall"
    assert viz.ax_.get_ylabel() == "Precision"
