import pytest
from numpy.testing import assert_allclose

from sklearn.datasets import load_iris
from sklearn.tree import DecisionTreeClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.base import ClassifierMixin
from sklearn.calibration import calibration_curve
from sklearn.metrics import plot_calibration_curve
from sklearn.exceptions import NotFittedError
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.compose import make_column_transformer


@pytest.fixture(scope="module")
def data():
    return load_iris(return_X_y=True)


@pytest.fixture(scope="module")
def data_binary(data):
    X, y = data
    return X[y < 2], y[y < 2]


def test_plot_calibration_curve_error_non_binary(pyplot, data):
    X, y = data
    clf = DecisionTreeClassifier()
    clf.fit(X, y)

    msg = "Only binary classification is supported."
    with pytest.raises(ValueError, match=msg):
        plot_calibration_curve(clf, X, y)


def test_plot_calibration_curve_no_predict_proba(pyplot, data_binary):
    X, y = data_binary

    class MyClassifier(ClassifierMixin):
        def fit(self, X, y):
            self.classes_ = [0, 1]
            return self

    clf = MyClassifier().fit(X, y)

    msg = "Response method 'predict_proba' not defined in"
    with pytest.raises(ValueError, match=msg):
        plot_calibration_curve(clf, X, y)


def test_plot_calibration_curve_not_fitted(pyplot, data_binary):
    X, y = data_binary
    clf = LogisticRegression()

    with pytest.raises(NotFittedError):
        plot_calibration_curve(clf, X, y)


@pytest.mark.parametrize("n_bins", [5, 10])
@pytest.mark.parametrize("strategy", ["uniform", "quantile"])
def test_plot_calibration_curve(pyplot, data_binary, n_bins, strategy):
    X, y = data_binary
    lr = LogisticRegression().fit(X, y)

    viz = plot_calibration_curve(
        lr, X, y, n_bins=n_bins, strategy=strategy, alpha=0.8
    )

    y_proba = lr.predict_proba(X)[:, 1]
    prob_true, prob_pred = calibration_curve(
        y, y_proba, n_bins=n_bins, strategy=strategy
    )

    assert_allclose(viz.prob_true, prob_true)
    assert_allclose(viz.prob_pred, prob_pred)

    assert viz.estimator_name == "LogisticRegression"

    # cannot fail thanks to pyplot fixture
    import matplotlib as mpl  # noqa
    assert isinstance(viz.line_, mpl.lines.Line2D)
    assert viz.line_.get_alpha() == 0.8
    assert isinstance(viz.ax_, mpl.axes.Axes)
    assert isinstance(viz.figure_, mpl.figure.Figure)

    assert viz.line_.get_label() == "LogisticRegression"
    assert viz.ax_.get_xlabel() == "Mean predicted probability"
    assert viz.ax_.get_ylabel() == "Fraction of positives"


@pytest.mark.parametrize(
    "clf", [make_pipeline(StandardScaler(), LogisticRegression()),
            make_pipeline(make_column_transformer((StandardScaler(), [0, 1])),
                          LogisticRegression())])
def test_plot_calibration_curve_pipeline(pyplot, data_binary, clf):
    X, y = data_binary
    clf.fit(X, y)
    viz = plot_calibration_curve(clf, X, y)
    assert clf.__class__.__name__ in viz.line_.get_label()
    assert viz.estimator_name == clf.__class__.__name__


def test_plot_roc_curve_estimator_name_multiple_calls(pyplot, data_binary):
    # non-regression test checking that the `name` used when calling
    # `plot_calibration_curve` is used as well when calling `viz.plot()`
    X, y = data_binary
    clf_name = "my hand-crafted name"
    clf = LogisticRegression().fit(X, y)
    viz = plot_calibration_curve(clf, X, y, name=clf_name)
    assert viz.estimator_name == clf_name
    pyplot.close("all")
    viz.plot()
    assert clf_name in viz.line_.get_label()
    pyplot.close("all")
    clf_name = "another_name"
    viz.plot(name=clf_name)
    assert clf_name in viz.line_.get_label()


def test_plot_calibration_curve_ref_line(pyplot, data_binary):
    # Check that `ref_line` only appears once
    X, y = data_binary
    lr = LogisticRegression().fit(X, y)
    dt = DecisionTreeClassifier().fit(X, y)

    viz = plot_calibration_curve(lr, X, y)
    viz2 = plot_calibration_curve(dt, X, y, ax=viz.ax_)

    labels = viz2.ax_.get_legend_handles_labels()[1]
    assert labels.count('Perfectly calibrated') == 1
