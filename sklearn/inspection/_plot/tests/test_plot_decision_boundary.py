import pytest

from sklearn.base import BaseEstimator
from sklearn.base import ClassifierMixin
from sklearn.inspection import plot_decision_boundary
from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression


# TODO: Remove when https://github.com/numpy/numpy/issues/14397 is resolved
pytestmark = pytest.mark.filterwarnings(
    "ignore:In future, it will be an error for 'np.bool_':DeprecationWarning:"
    "matplotlib.*")


@pytest.fixture(scope="module")
def data():
    X, y = make_classification(n_informative=1, n_redundant=1,
                               n_clusters_per_class=1, n_features=2)
    return X, y


@pytest.fixture(scope="module")
def fitted_clf(data):
    return LogisticRegression().fit(*data)


@pytest.mark.parametrize("response_method",
                         ['auto', 'predict_proba', 'decision_function'])
def test_multiclass_error(pyplot, response_method):
    X, y = make_classification(n_classes=3, n_informative=3, random_state=0)
    X = X[:, [0, 1]]
    lr = LogisticRegression().fit(X, y)

    msg = ("multiclass classifers are only supported when "
           "response_method='predict'")
    with pytest.raises(ValueError, match=msg):
        plot_decision_boundary(lr, X, response_method=response_method)


@pytest.mark.parametrize("kwargs, error_msg", [
    ({"plot_method": "hello_world"},
     r"plot_method must be 'contourf',"),
    ({"grid_resolution": 1},
     r"grid_resolution must be greater than 1"),
    ({"grid_resolution": -1},
     r"grid_resolution must be greater than 1"),
    ({"eps": -1.1},
     r"eps must be greater than or equal to 0")
])
def test_input_validation_errors(pyplot, kwargs, error_msg, fitted_clf, data):
    X, _ = data
    with pytest.raises(ValueError, match=error_msg):
        plot_decision_boundary(fitted_clf, X, **kwargs)


def test_display_plot_input_error(pyplot, fitted_clf, data):
    X, y = data
    disp = plot_decision_boundary(fitted_clf, X, grid_resolution=5)

    with pytest.raises(ValueError, match="plot_method must be 'contourf'"):
        disp.plot(plot_method="hello_world")


@pytest.mark.parametrize("response_method", ['auto', 'predict',
                                             'predict_proba',
                                             'decision_function'])
@pytest.mark.parametrize("plot_method", ['contourf', 'contour'])
def test_plot_decision_boundary(pyplot, fitted_clf, data,
                                response_method, plot_method):
    fig, ax = pyplot.subplots()
    eps = 2.0
    X, y = data
    disp = plot_decision_boundary(fitted_clf, X, grid_resolution=5,
                                  response_method=response_method,
                                  plot_method=plot_method,
                                  eps=eps, ax=ax)
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

    # change plotting method for second plot
    disp.plot(plot_method='pcolormesh')
    assert isinstance(disp.surface_, pyplot.matplotlib.collections.QuadMesh)


@pytest.mark.parametrize(
    "response_method, msg",
    [("predict_proba", "response method predict_proba is not defined in "
                       "MyClassifier"),
     ("decision_function", "response method decision_function is not defined "
                           "in MyClassifier"),
     ("auto", "response method decision_function or predict_proba is not "
              "defined in MyClassifier"),
     ("bad_method", "response_method must be 'predict_proba', "
                    "'decision_function', 'predict', or 'auto'")])
def test_error_bad_response(pyplot, response_method, msg):
    X, y = make_classification(n_classes=2, n_samples=50, random_state=0)

    class MyClassifier(BaseEstimator, ClassifierMixin):
        def fit(self, X, y):
            self.fitted_ = True
            self.classes_ = [0, 1]
            return self

    clf = MyClassifier().fit(X, y)

    with pytest.raises(ValueError, match=msg):
        plot_decision_boundary(clf, X, response_method=response_method)
