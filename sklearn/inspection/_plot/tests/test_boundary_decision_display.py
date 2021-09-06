import pytest

from sklearn.base import BaseEstimator
from sklearn.base import ClassifierMixin
from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression

from sklearn.inspection import DecisionBoundaryDisplay
from sklearn.inspection._plot.decision_boundary import _check_boundary_response_method


# TODO: Remove when https://github.com/numpy/numpy/issues/14397 is resolved
pytestmark = pytest.mark.filterwarnings(
    "ignore:In future, it will be an error for 'np.bool_':DeprecationWarning:"
    "matplotlib.*"
)


@pytest.fixture(scope="module")
def data():
    X, y = make_classification(
        n_informative=1,
        n_redundant=1,
        n_clusters_per_class=1,
        n_features=2,
        random_state=42,
    )
    return X, y


@pytest.fixture(scope="module")
def fitted_clf(data):
    return LogisticRegression().fit(*data)


def test_check_boundary_response_method_auto():
    class A:
        def decision_function(self):
            pass

    a_inst = A()
    method = _check_boundary_response_method(a_inst, "auto")
    assert method == a_inst.decision_function

    class B:
        def predict_proba(self):
            pass

    b_inst = B()
    method = _check_boundary_response_method(b_inst, "auto")
    assert method == b_inst.predict_proba

    class C:
        def predict_proba(self):
            pass

        def decision_function(self):
            pass

    c_inst = C()
    method = _check_boundary_response_method(c_inst, "auto")
    assert method == c_inst.decision_function

    class D:
        def predict(self):
            pass

    d_inst = D()
    method = _check_boundary_response_method(d_inst, "auto")
    assert method == d_inst.predict


@pytest.mark.parametrize(
    "response_method", ["auto", "predict_proba", "decision_function"]
)
def test_multiclass_error(pyplot, response_method):
    X, y = make_classification(n_classes=3, n_informative=3, random_state=0)
    X = X[:, [0, 1]]
    lr = LogisticRegression().fit(X, y)

    msg = "Multiclass classifiers are only supported when response_method='predict'"
    with pytest.raises(ValueError, match=msg):
        DecisionBoundaryDisplay.from_estimator(lr, X, response_method=response_method)


@pytest.mark.parametrize(
    "kwargs, error_msg",
    [
        (
            {"plot_method": "hello_world"},
            r"plot_method must be one of contourf, contour, pcolormesh. Got hello_world"
            r" instead.",
        ),
        (
            {"grid_resolution": 1},
            r"grid_resolution must be greater than 1. Got 1 instead",
        ),
        (
            {"grid_resolution": -1},
            r"grid_resolution must be greater than 1. Got -1 instead",
        ),
        ({"eps": -1.1}, r"eps must be greater than or equal to 0. Got -1.1 instead"),
    ],
)
def test_input_validation_errors(pyplot, kwargs, error_msg, fitted_clf, data):
    X, _ = data
    with pytest.raises(ValueError, match=error_msg):
        DecisionBoundaryDisplay.from_estimator(fitted_clf, X, **kwargs)


def test_display_plot_input_error(pyplot, fitted_clf, data):
    X, y = data
    disp = DecisionBoundaryDisplay.from_estimator(fitted_clf, X, grid_resolution=5)

    with pytest.raises(ValueError, match="plot_method must be 'contourf'"):
        disp.plot(plot_method="hello_world")


@pytest.mark.parametrize(
    "response_method", ["auto", "predict", "predict_proba", "decision_function"]
)
@pytest.mark.parametrize("plot_method", ["contourf", "contour"])
def test_decision_boundary_display(
    pyplot, fitted_clf, data, response_method, plot_method
):
    fig, ax = pyplot.subplots()
    eps = 2.0
    X, y = data
    disp = DecisionBoundaryDisplay.from_estimator(
        fitted_clf,
        X,
        grid_resolution=5,
        response_method=response_method,
        plot_method=plot_method,
        eps=eps,
        ax=ax,
    )
    assert isinstance(disp.surface_, pyplot.matplotlib.contour.QuadContourSet)
    assert disp.ax_ == ax
    assert disp.figure_ == fig

    x0, x1 = X[:, 0], X[:, 1]

    x0_min, x0_max = x0.min() - eps, x0.max() + eps
    x1_min, x1_max = x1.min() - eps, x1.max() + eps

    assert disp.xx0.min() == pytest.approx(x0_min)
    assert disp.xx0.max() == pytest.approx(x0_max)
    assert disp.xx1.min() == pytest.approx(x1_min)
    assert disp.xx1.max() == pytest.approx(x1_max)

    fig2, ax2 = pyplot.subplots()
    # change plotting method for second plot
    disp.plot(plot_method="pcolormesh", ax=ax2, shading="auto")
    assert isinstance(disp.surface_, pyplot.matplotlib.collections.QuadMesh)
    assert disp.ax_ == ax2
    assert disp.figure_ == fig2


@pytest.mark.parametrize(
    "response_method, msg",
    [
        (
            "predict_proba",
            "response method predict_proba is not defined in MyClassifier",
        ),
        (
            "decision_function",
            "response method decision_function is not defined in MyClassifier",
        ),
        (
            "auto",
            "response method decision_function, predict_proba, or predict "
            "is not defined in MyClassifier",
        ),
        (
            "bad_method",
            "response_method must be one of predict_proba, decision_function, auto,"
            " predict",
        ),
    ],
)
def test_error_bad_response(pyplot, response_method, msg, data):
    X, y = data

    class MyClassifier(BaseEstimator, ClassifierMixin):
        def fit(self, X, y):
            self.fitted_ = True
            self.classes_ = [0, 1]
            return self

    clf = MyClassifier().fit(X, y)

    with pytest.raises(ValueError, match=msg):
        DecisionBoundaryDisplay.from_estimator(clf, X, response_method=response_method)


def test_dataframe_labels_used(pyplot, data, fitted_clf):
    pd = pytest.importorskip("pandas")
    df = pd.DataFrame(data[0], columns=["col_x", "col_y"])

    # pandas column names are used by default
    _, ax = pyplot.subplots()
    disp = DecisionBoundaryDisplay.from_estimator(fitted_clf, df, ax=ax)
    assert ax.get_xlabel() == "col_x"
    assert ax.get_ylabel() == "col_y"

    # second call to plot will have the names
    fig, ax = pyplot.subplots()
    disp.plot(ax=ax)
    assert ax.get_xlabel() == "col_x"
    assert ax.get_ylabel() == "col_y"

    # axes with a label will not get overridden
    fig, ax = pyplot.subplots()
    ax.set(xlabel="hello", ylabel="world")
    disp.plot(ax=ax)
    assert ax.get_xlabel() == "hello"
    assert ax.get_ylabel() == "world"

    # labels get overriden only if provided to the `plot` method
    disp.plot(ax=ax, xlabel="overwritten_x", ylabel="overwritten_y")
    assert ax.get_xlabel() == "overwritten_x"
    assert ax.get_ylabel() == "overwritten_y"

    # labels do not get inferred if provided to `from_estimator`
    _, ax = pyplot.subplots()
    disp = DecisionBoundaryDisplay.from_estimator(
        fitted_clf, df, ax=ax, xlabel="overwritten_x", ylabel="overwritten_y"
    )
    assert ax.get_xlabel() == "overwritten_x"
    assert ax.get_ylabel() == "overwritten_y"
