import pytest
import numpy as np
from numpy.testing import assert_allclose

from sklearn.base import BaseEstimator
from sklearn.metrics import plot_precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve
from sklearn.datasets import load_breast_cancer
from sklearn.datasets import load_iris
from sklearn.tree import DecisionTreeClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.exceptions import NotFittedError


@pytest.fixture(scope="module")
def data_binary():
    return load_breast_cancer(return_X_y=True)


def test_error_non_binary(pyplot):
    X, y = load_iris(return_X_y=True)
    clf = DecisionTreeClassifier()
    clf.fit(X, y)

    msg = "multiclass format is not supported"
    with pytest.raises(ValueError, match=msg):
        plot_precision_recall_curve(clf, X, y)

    msg = "Estimator should solve a binary classification problem"
    y_binary = y >= 1
    with pytest.raises(ValueError, match=msg):
        plot_precision_recall_curve(clf, X, y_binary)


def test_unfitted_classifier(pyplot, data_binary):
    X, y = data_binary
    clf = DecisionTreeClassifier()
    with pytest.raises(NotFittedError):
        plot_precision_recall_curve(clf, X, y)


@pytest.mark.parametrize(
    "response_method, msg",
    [("predict_proba", "response method predict_proba not defined for "
                       "estimator MyClassifier"),
     ("decision_function", "response method decision_function not defined for "
                           "estimator MyClassifier"),
     ("auto", "response method decision_function or predict_proba not defined "
              "for estimator MyClassifier"),
     ("bad_method", "response_method must be 'predict_proba', "
                    "'decision_function' or 'auto'")])
def test_error_no_response(pyplot, data_binary, response_method, msg):
    X, y = data_binary

    class MyClassifier(BaseEstimator):
        def fit(self, X, y):
            self.fitted_ = True
            return self

    clf = MyClassifier().fit(X, y)

    with pytest.raises(ValueError, match=msg):
        plot_precision_recall_curve(clf, X, y, response_method=response_method)


@pytest.mark.parametrize("response_method",
                         ["predict_proba", "decision_function"])
@pytest.mark.parametrize("with_sample_weight", [True, False])
def test_plot_precision_recall(pyplot, response_method, data_binary,
                               with_sample_weight):
    X, y = data_binary

    lr = LogisticRegression().fit(X, y)

    if with_sample_weight:
        rng = np.random.RandomState(42)
        sample_weight = rng.randint(0, 4, size=X.shape[0])
    else:
        sample_weight = None

    disp = plot_precision_recall_curve(lr, X, y, alpha=0.8,
                                       response_method=response_method,
                                       sample_weight=sample_weight)

    y_score = getattr(lr, response_method)(X)
    if y_score.ndim == 2:
        y_score = y_score[:, 1]

    prec, recall, _ = precision_recall_curve(y, y_score,
                                             sample_weight=sample_weight)
    avg_prec = average_precision_score(y, y_score, sample_weight=sample_weight)

    assert_allclose(disp.precision, prec)
    assert_allclose(disp.recall, recall)
    assert disp.average_precision == pytest.approx(avg_prec)

    assert disp.estimator_name == "LogisticRegression"

    # cannot fail thanks to pyplot fixture
    import matplotlib as mpl  # noqa
    assert isinstance(disp.line_, mpl.lines.Line2D)
    assert disp.line_.get_alpha() == 0.8
    assert isinstance(disp.ax_, mpl.axes.Axes)
    assert isinstance(disp.figure_, mpl.figure.Figure)

    expected_label = "LogisticRegression (AP = {:0.2f})".format(avg_prec)
    assert disp.line_.get_label() == expected_label
    assert disp.ax_.get_xlabel() == "Recall"
    assert disp.ax_.get_ylabel() == "Precision"
    assert_allclose(disp.ax_.get_xlim(), [0.0, 1.0])
    assert_allclose(disp.ax_.get_ylim(), [0.0, 1.05])

    # draw again with another label
    disp.plot(label_name="MySpecialEstimator")
    expected_label = "MySpecialEstimator (AP = {:0.2f})".format(avg_prec)
    assert disp.line_.get_label() == expected_label
