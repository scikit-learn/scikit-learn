import pytest
import numpy as np
from numpy.testing import assert_allclose

from sklearn.metrics import plot_precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve
from sklearn.datasets import load_breast_cancer
from sklearn.datasets import load_iris
from sklearn.tree import DecisionTreeClassifier
from sklearn.linear_model import LogisticRegression


@pytest.fixture(scope="module")
def data_binary():
    return load_breast_cancer(return_X_y=True)


def test_error_non_binary(pyplot):
    X, y = load_iris(return_X_y=True)
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
def test_error_no_response(pyplot, data_binary, response_method, msg):
    X, y = data_binary

    class MyClassifier:
        pass

    clf = MyClassifier()

    with pytest.raises(ValueError, match=msg):
        plot_precision_recall_curve(clf, X, y, response_method=response_method)


@pytest.mark.parametrize("response_method",
                         ["predict_proba", "decision_function"])
@pytest.mark.parametrize("with_sample_weight", [True, False])
def test_plot_precision_recall(pyplot, response_method, data_binary,
                               with_sample_weight):
    X, y = data_binary

    lr = LogisticRegression()
    lr.fit(X, y)

    if with_sample_weight:
        rng = np.random.RandomState(42)
        sample_weight = rng.randint(0, 4, size=X.shape[0])
    else:
        sample_weight = None

    viz = plot_precision_recall_curve(lr, X, y, alpha=0.8,
                                      sample_weight=sample_weight)

    y_score = getattr(lr, response_method)(X)
    if y_score.ndim == 2:
        y_score = y_score[:, 1]

    prec, recall, _ = precision_recall_curve(y, y_score,
                                             sample_weight=sample_weight)
    avg_prec = average_precision_score(y, y_score, sample_weight=sample_weight)

    assert_allclose(viz.precision, prec)
    assert_allclose(viz.recall, recall)
    assert_allclose(viz.average_precision, avg_prec)

    assert viz.estimator_name == "LogisticRegression"

    # cannot fail thanks to pyplot fixture
    import matplotlib as mpl  # noqa
    assert isinstance(viz.line_, mpl.lines.Line2D)
    assert viz.line_.get_alpha() == 0.8
    assert isinstance(viz.ax_, mpl.axes.Axes)
    assert isinstance(viz.figure_, mpl.figure.Figure)

    expected_label = "LogisticRegression (AP = {:0.2f})".format(avg_prec)
    assert viz.line_.get_label() == expected_label
    assert viz.ax_.get_xlabel() == "Recall"
    assert viz.ax_.get_ylabel() == "Precision"
    assert_allclose(viz.ax_.get_xlim(), [0.0, 1.0])
    assert_allclose(viz.ax_.get_ylim(), [0.0, 1.05])
