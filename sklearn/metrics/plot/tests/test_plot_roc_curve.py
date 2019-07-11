import pytest
from numpy.testing import assert_allclose

from sklearn.metrics import plot_roc_curve
from sklearn.datasets import load_iris
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve


@pytest.fixture(scope="module")
def data():
    X, y = load_iris(return_X_y=True)
    return X, y


def test_plot_roc_curve_error_non_binary(pyplot, data):
    X, y = data
    lr = LogisticRegression(solver="lbfgs",
                            multi_class="multinomial")
    lr.fit(X, y)

    msg = "Estimator must be a binary classifier"
    with pytest.raises(ValueError, match=msg):
        plot_roc_curve(lr, X, y)


@pytest.mark.parametrize("response_method",
                         ["predict_proba", "decision_function"])
@pytest.mark.parametrize("pos_label", [0, 1])
@pytest.mark.parametrize("drop_intermediate", [True, False])
def test_plot_roc_curve(pyplot, response_method, data,
                        pos_label, drop_intermediate):
    X, y = data
    X, y = X[y < 2], y[y < 2]

    lr = LogisticRegression(solver="lbfgs")
    lr.fit(X, y)

    viz = plot_roc_curve(lr, X, y, pos_label=pos_label,
                         drop_intermediate=drop_intermediate)

    y_pred = getattr(lr, response_method)(X)
    if y_pred.ndim == 2:
        y_pred = y_pred[:, 1]

    fpr, tpr, _ = roc_curve(y, y_pred, pos_label=pos_label,
                            drop_intermediate=drop_intermediate)

    assert viz.label_ == "LogisticRegression"

    assert_allclose(viz.fpr_, fpr)
    assert_allclose(viz.tpr_, tpr)

    import matplotlib as mpl
    assert isinstance(viz.line_, mpl.lines.Line2D)
    assert isinstance(viz.ax_, mpl.axes.Axes)
    assert isinstance(viz.figure_, mpl.figure.Figure)
