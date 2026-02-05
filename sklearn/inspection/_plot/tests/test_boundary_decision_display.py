import warnings

import numpy as np
import pytest

from sklearn.base import (
    BaseEstimator,
    ClassifierMixin,
)
from sklearn.cluster import KMeans
from sklearn.datasets import (
    load_diabetes,
    load_iris,
    make_blobs,
    make_classification,
    make_multilabel_classification,
)
from sklearn.ensemble import IsolationForest
from sklearn.inspection import DecisionBoundaryDisplay
from sklearn.inspection._plot.decision_boundary import _check_boundary_response_method
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, scale
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor
from sklearn.utils._testing import (
    _convert_container,
    assert_allclose,
    assert_array_equal,
)
from sklearn.utils.fixes import parse_version

X, y = make_classification(
    n_informative=1,
    n_redundant=1,
    n_clusters_per_class=1,
    n_features=2,
    random_state=42,
)


def load_iris_2d_scaled():
    X, y = load_iris(return_X_y=True)
    X = scale(X)[:, :2]
    return X, y


@pytest.fixture(scope="module")
def fitted_clf():
    return LogisticRegression().fit(X, y)


def test_input_data_dimension(pyplot):
    """Check that we raise an error when `X` does not have exactly 2 features."""
    X, y = make_classification(n_samples=10, n_features=4, random_state=0)

    clf = LogisticRegression().fit(X, y)
    msg = "n_features must be equal to 2. Got 4 instead."
    with pytest.raises(ValueError, match=msg):
        DecisionBoundaryDisplay.from_estimator(estimator=clf, X=X)


def test_check_boundary_response_method_error():
    """Check error raised for multi-output multi-class classifiers by
    `_check_boundary_response_method`.
    """

    class MultiLabelClassifier:
        classes_ = [np.array([0, 1]), np.array([0, 1])]

    err_msg = "Multi-label and multi-output multi-class classifiers are not supported"
    with pytest.raises(ValueError, match=err_msg):
        _check_boundary_response_method(MultiLabelClassifier(), "predict")


@pytest.mark.parametrize(
    "estimator, response_method, expected_prediction_method",
    [
        (DecisionTreeRegressor(), "predict", "predict"),
        (DecisionTreeRegressor(), "auto", "predict"),
        (LogisticRegression().fit(*load_iris_2d_scaled()), "predict", "predict"),
        (
            LogisticRegression().fit(*load_iris_2d_scaled()),
            "auto",
            ["decision_function", "predict_proba", "predict"],
        ),
        (
            LogisticRegression().fit(*load_iris_2d_scaled()),
            "predict_proba",
            "predict_proba",
        ),
        (
            LogisticRegression().fit(*load_iris_2d_scaled()),
            "decision_function",
            "decision_function",
        ),
        (
            LogisticRegression().fit(X, y),
            "auto",
            ["decision_function", "predict_proba", "predict"],
        ),
        (LogisticRegression().fit(X, y), "predict", "predict"),
        (
            LogisticRegression().fit(X, y),
            ["predict_proba", "decision_function"],
            ["predict_proba", "decision_function"],
        ),
    ],
)
def test_check_boundary_response_method(
    estimator, response_method, expected_prediction_method
):
    """Check the behaviour of `_check_boundary_response_method` for the supported
    cases.
    """
    prediction_method = _check_boundary_response_method(estimator, response_method)
    assert prediction_method == expected_prediction_method


def test_multiclass_predict(pyplot):
    """Check multiclass `response=predict` gives expected results."""
    grid_resolution = 10
    eps = 1.0
    X, y = make_classification(n_classes=3, n_informative=3, random_state=0)
    X = X[:, [0, 1]]
    lr = LogisticRegression(random_state=0).fit(X, y)

    disp = DecisionBoundaryDisplay.from_estimator(
        lr, X, response_method="predict", grid_resolution=grid_resolution, eps=1.0
    )

    x0_min, x0_max = X[:, 0].min() - eps, X[:, 0].max() + eps
    x1_min, x1_max = X[:, 1].min() - eps, X[:, 1].max() + eps
    xx0, xx1 = np.meshgrid(
        np.linspace(x0_min, x0_max, grid_resolution),
        np.linspace(x1_min, x1_max, grid_resolution),
    )
    response = lr.predict(np.c_[xx0.ravel(), xx1.ravel()])
    assert_allclose(disp.response, response.reshape(xx0.shape))
    assert_allclose(disp.xx0, xx0)
    assert_allclose(disp.xx1, xx1)


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
def test_input_validation_errors(pyplot, kwargs, error_msg, fitted_clf):
    """Check input validation from_estimator."""
    with pytest.raises(ValueError, match=error_msg):
        DecisionBoundaryDisplay.from_estimator(fitted_clf, X, **kwargs)


@pytest.mark.parametrize(
    "kwargs, error_msg",
    [
        (
            {"multiclass_colors": {"dict": "not_list"}},
            "'multiclass_colors' must be a list or a str.",
        ),
        ({"multiclass_colors": "not_cmap"}, "it must be a valid Matplotlib colormap"),
        ({"multiclass_colors": ["red", "green"]}, "it must be of the same length"),
        (
            {"multiclass_colors": ["red", "green", "not color"]},
            "it can only contain valid Matplotlib color names",
        ),
    ],
)
def test_input_validation_errors_multiclass_colors(pyplot, kwargs, error_msg):
    """Check input validation for `multiclass_colors` in `from_estimator`."""
    X, y = load_iris_2d_scaled()
    clf = LogisticRegression().fit(X, y)
    with pytest.raises(ValueError, match=error_msg):
        DecisionBoundaryDisplay.from_estimator(clf, X, **kwargs)


def test_display_plot_input_error(pyplot, fitted_clf):
    """Check input validation for `plot`."""
    disp = DecisionBoundaryDisplay.from_estimator(fitted_clf, X, grid_resolution=5)

    with pytest.raises(ValueError, match="plot_method must be 'contourf'"):
        disp.plot(plot_method="hello_world")


@pytest.mark.parametrize(
    "response_method", ["auto", "predict", "predict_proba", "decision_function"]
)
@pytest.mark.parametrize("plot_method", ["contourf", "contour"])
def test_decision_boundary_display_classifier(
    pyplot, fitted_clf, response_method, plot_method
):
    """Check that decision boundary is correct."""
    fig, ax = pyplot.subplots()
    eps = 2.0
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


@pytest.mark.parametrize("response_method", ["auto", "predict", "decision_function"])
@pytest.mark.parametrize("plot_method", ["contourf", "contour", "pcolormesh"])
def test_decision_boundary_display_outlier_detector(
    pyplot, response_method, plot_method
):
    """Check that decision boundary is correct for outlier detector."""
    fig, ax = pyplot.subplots()
    eps = 2.0
    outlier_detector = IsolationForest(random_state=0).fit(X, y)
    disp = DecisionBoundaryDisplay.from_estimator(
        outlier_detector,
        X,
        grid_resolution=5,
        response_method=response_method,
        plot_method=plot_method,
        eps=eps,
        ax=ax,
    )
    if plot_method == "pcolormesh":
        assert isinstance(disp.surface_, pyplot.matplotlib.collections.QuadMesh)
    else:
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


@pytest.mark.parametrize("response_method", ["auto", "predict"])
@pytest.mark.parametrize("plot_method", ["contourf", "contour"])
def test_decision_boundary_display_regressor(pyplot, response_method, plot_method):
    """Check that we can display the decision boundary for a regressor."""
    X, y = load_diabetes(return_X_y=True)
    X = X[:, :2]
    tree = DecisionTreeRegressor().fit(X, y)
    fig, ax = pyplot.subplots()
    eps = 2.0
    disp = DecisionBoundaryDisplay.from_estimator(
        tree,
        X,
        response_method=response_method,
        ax=ax,
        eps=eps,
        plot_method=plot_method,
    )
    if disp.n_classes == 2 or plot_method == "contour":
        assert isinstance(disp.surface_, pyplot.matplotlib.contour.QuadContourSet)
    else:
        assert isinstance(disp.surface_, list)
        for surface in disp.surface_:
            assert isinstance(surface, pyplot.matplotlib.contour.QuadContourSet)
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
            "MyClassifier has none of the following attributes: predict_proba",
        ),
        (
            "decision_function",
            "MyClassifier has none of the following attributes: decision_function",
        ),
        (
            "auto",
            (
                "MyClassifier has none of the following attributes: decision_function, "
                "predict_proba, predict"
            ),
        ),
        (
            "bad_method",
            "MyClassifier has none of the following attributes: bad_method",
        ),
    ],
)
def test_error_bad_response(pyplot, response_method, msg):
    """Check errors for bad response."""

    class MyClassifier(ClassifierMixin, BaseEstimator):
        def fit(self, X, y):
            self.fitted_ = True
            self.classes_ = [0, 1]
            return self

    clf = MyClassifier().fit(X, y)

    with pytest.raises(AttributeError, match=msg):
        DecisionBoundaryDisplay.from_estimator(clf, X, response_method=response_method)


@pytest.mark.parametrize("response_method", ["auto", "predict", "predict_proba"])
def test_multilabel_classifier_error(pyplot, response_method):
    """Check that multilabel classifier raises correct error."""
    X, y = make_multilabel_classification(random_state=0)
    X = X[:, :2]
    tree = DecisionTreeClassifier().fit(X, y)

    msg = "Multi-label and multi-output multi-class classifiers are not supported"
    with pytest.raises(ValueError, match=msg):
        DecisionBoundaryDisplay.from_estimator(
            tree,
            X,
            response_method=response_method,
        )


@pytest.mark.parametrize("response_method", ["auto", "predict", "predict_proba"])
def test_multi_output_multi_class_classifier_error(pyplot, response_method):
    """Check that multi-output multi-class classifier raises correct error."""
    X = np.asarray([[0, 1], [1, 2]])
    y = np.asarray([["tree", "cat"], ["cat", "tree"]])
    tree = DecisionTreeClassifier().fit(X, y)

    msg = "Multi-label and multi-output multi-class classifiers are not supported"
    with pytest.raises(ValueError, match=msg):
        DecisionBoundaryDisplay.from_estimator(
            tree,
            X,
            response_method=response_method,
        )


def test_multioutput_regressor_error(pyplot):
    """Check that multioutput regressor raises correct error."""
    X = np.asarray([[0, 1], [1, 2]])
    y = np.asarray([[0, 1], [4, 1]])
    tree = DecisionTreeRegressor().fit(X, y)
    with pytest.raises(ValueError, match="Multi-output regressors are not supported"):
        DecisionBoundaryDisplay.from_estimator(tree, X, response_method="predict")


@pytest.mark.filterwarnings(
    # We expect to raise the following warning because the classifier is fit on a
    # NumPy array
    "ignore:X has feature names, but LogisticRegression was fitted without"
)
def test_dataframe_labels_used(pyplot, fitted_clf):
    """Check that column names are used for pandas."""
    pd = pytest.importorskip("pandas")
    df = pd.DataFrame(X, columns=["col_x", "col_y"])

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

    # labels get overridden only if provided to the `plot` method
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


def test_string_target(pyplot):
    """Check that decision boundary works with classifiers trained on string labels."""
    iris = load_iris()
    X = iris.data[:, [0, 1]]

    # Use strings as target
    y = iris.target_names[iris.target]
    log_reg = LogisticRegression().fit(X, y)

    # Does not raise
    DecisionBoundaryDisplay.from_estimator(
        log_reg,
        X,
        grid_resolution=5,
        response_method="predict",
    )


@pytest.mark.parametrize("constructor_name", ["pandas", "polars"])
def test_dataframe_support(pyplot, constructor_name):
    """Check that passing a dataframe at fit and to the Display does not
    raise warnings.

    Non-regression test for:
    * https://github.com/scikit-learn/scikit-learn/issues/23311
    * https://github.com/scikit-learn/scikit-learn/issues/28717
    """
    df = _convert_container(
        X, constructor_name=constructor_name, columns_name=["col_x", "col_y"]
    )
    estimator = LogisticRegression().fit(df, y)

    with warnings.catch_warnings():
        # no warnings linked to feature names validation should be raised
        warnings.simplefilter("error", UserWarning)
        DecisionBoundaryDisplay.from_estimator(estimator, df, response_method="predict")


@pytest.mark.parametrize("response_method", ["predict_proba", "decision_function"])
def test_class_of_interest_binary(pyplot, response_method):
    """Check the behaviour of passing `class_of_interest` for plotting the output of
    `predict_proba` and `decision_function` in the binary case.
    """
    iris = load_iris()
    X = iris.data[:100, :2]
    y = iris.target[:100]
    assert_array_equal(np.unique(y), [0, 1])

    estimator = LogisticRegression().fit(X, y)
    # We will check that `class_of_interest=None` is equivalent to
    # `class_of_interest=estimator.classes_[1]`
    disp_default = DecisionBoundaryDisplay.from_estimator(
        estimator,
        X,
        response_method=response_method,
        class_of_interest=None,
    )
    disp_class_1 = DecisionBoundaryDisplay.from_estimator(
        estimator,
        X,
        response_method=response_method,
        class_of_interest=estimator.classes_[1],
    )

    assert_allclose(disp_default.response, disp_class_1.response)

    # we can check that `_get_response_values` modifies the response when targeting
    # the other class, i.e. 1 - p(y=1|x) for `predict_proba` and -decision_function
    # for `decision_function`.
    disp_class_0 = DecisionBoundaryDisplay.from_estimator(
        estimator,
        X,
        response_method=response_method,
        class_of_interest=estimator.classes_[0],
    )

    if response_method == "predict_proba":
        assert_allclose(disp_default.response, 1 - disp_class_0.response)
    else:
        assert response_method == "decision_function"
        assert_allclose(disp_default.response, -disp_class_0.response)


@pytest.mark.parametrize("response_method", ["predict_proba", "decision_function"])
def test_class_of_interest_multiclass(pyplot, response_method):
    """Check the behaviour of passing `class_of_interest` for plotting the output of
    `predict_proba` and `decision_function` in the multiclass case.
    """
    iris = load_iris()
    X = iris.data[:, :2]
    y = iris.target  # the target are numerical labels
    class_of_interest_idx = 2

    estimator = LogisticRegression().fit(X, y)
    disp = DecisionBoundaryDisplay.from_estimator(
        estimator,
        X,
        response_method=response_method,
        class_of_interest=class_of_interest_idx,
    )

    # we will check that we plot the expected values as response
    grid = np.concatenate([disp.xx0.reshape(-1, 1), disp.xx1.reshape(-1, 1)], axis=1)
    response = getattr(estimator, response_method)(grid)[:, class_of_interest_idx]
    assert_allclose(response.reshape(*disp.response.shape), disp.response)

    # make the same test but this time using target as strings
    y = iris.target_names[iris.target]
    estimator = LogisticRegression().fit(X, y)

    disp = DecisionBoundaryDisplay.from_estimator(
        estimator,
        X,
        response_method=response_method,
        class_of_interest=iris.target_names[class_of_interest_idx],
    )

    grid = np.concatenate([disp.xx0.reshape(-1, 1), disp.xx1.reshape(-1, 1)], axis=1)
    response = getattr(estimator, response_method)(grid)[:, class_of_interest_idx]
    assert_allclose(response.reshape(*disp.response.shape), disp.response)

    # check that we raise an error for unknown labels
    # this test should already be handled in `_get_response_values` but we can have this
    # test here as well
    err_msg = "class_of_interest=2 is not a valid label: It should be one of"
    with pytest.raises(ValueError, match=err_msg):
        DecisionBoundaryDisplay.from_estimator(
            estimator,
            X,
            response_method=response_method,
            class_of_interest=class_of_interest_idx,
        )


@pytest.mark.parametrize("response_method", ["predict_proba", "decision_function"])
def test_multiclass_plot_max_class(pyplot, response_method):
    """Check plot correct when plotting max multiclass class."""
    import matplotlib as mpl

    # In matplotlib < v3.5, default value of `pcolormesh(shading)` is 'flat', which
    # results in the last row and column being dropped. Thus older versions produce
    # a 99x99 grid, while newer versions produce a 100x100 grid.
    if parse_version(mpl.__version__) < parse_version("3.5"):
        pytest.skip("`pcolormesh` in Matplotlib >= 3.5 gives smaller grid size.")

    X, y = load_iris_2d_scaled()
    clf = LogisticRegression().fit(X, y)

    disp = DecisionBoundaryDisplay.from_estimator(
        clf,
        X,
        plot_method="pcolormesh",
        response_method=response_method,
    )

    grid = np.concatenate([disp.xx0.reshape(-1, 1), disp.xx1.reshape(-1, 1)], axis=1)
    response = getattr(clf, response_method)(grid).reshape(*disp.response.shape)
    assert_allclose(response, disp.response)

    assert len(disp.surface_) == len(clf.classes_)
    # Get which class has highest response and check it is plotted
    highest_class = np.argmax(response, axis=2)
    for idx, quadmesh in enumerate(disp.surface_):
        # Note quadmesh mask is True (i.e. masked) when `idx` is NOT the highest class
        assert_array_equal(
            highest_class != idx,
            quadmesh.get_array().mask.reshape(*highest_class.shape),
        )


@pytest.mark.parametrize(
    "multiclass_colors, n_classes",
    [
        (None, 3),
        (None, 11),
        ("plasma", 3),
        ("Blues", 3),
        (["red", "green", "blue"], 3),
    ],
)
@pytest.mark.parametrize(
    "response_method", ["decision_function", "predict_proba", "predict"]
)
@pytest.mark.parametrize("plot_method", ["contourf", "contour", "pcolormesh"])
def test_multiclass_colors_cmap(
    pyplot,
    n_classes,
    response_method,
    plot_method,
    multiclass_colors,
):
    """Check correct cmap used for all `multiclass_colors` inputs."""
    import matplotlib as mpl

    if parse_version(mpl.__version__) < parse_version("3.5"):
        pytest.skip(
            "Matplotlib >= 3.5 is needed for `==` to check equivalence of colormaps"
        )

    X, y = make_blobs(n_samples=150, centers=n_classes, n_features=2, random_state=42)
    clf = LogisticRegression().fit(X, y)

    disp = DecisionBoundaryDisplay.from_estimator(
        clf,
        X,
        response_method=response_method,
        plot_method=plot_method,
        multiclass_colors=multiclass_colors,
    )

    if multiclass_colors is None:
        if len(clf.classes_) <= 10:
            colors = mpl.pyplot.get_cmap("tab10", 10).colors[: len(clf.classes_)]
        else:
            cmap = mpl.pyplot.get_cmap("gist_rainbow", len(clf.classes_))
            colors = cmap(np.linspace(0, 1, len(clf.classes_)))
    elif multiclass_colors == "plasma":
        colors = mpl.pyplot.get_cmap(multiclass_colors, len(clf.classes_)).colors
    elif multiclass_colors == "Blues":
        cmap = mpl.pyplot.get_cmap(multiclass_colors, len(clf.classes_))
        colors = cmap(np.linspace(0, 1, len(clf.classes_)))
    else:
        colors = [mpl.colors.to_rgba(color) for color in multiclass_colors]

    # Make sure the colormap has enough distinct colors.
    assert disp.n_classes == len(np.unique(colors, axis=0))

    if response_method == "predict":
        cmap = mpl.colors.ListedColormap(colors)
        assert disp.surface_.cmap == cmap
    elif plot_method != "contour":
        cmaps = [
            mpl.colors.LinearSegmentedColormap.from_list(
                f"colormap_{class_idx}", [(1.0, 1.0, 1.0, 1.0), (r, g, b, 1.0)]
            )
            for class_idx, (r, g, b, _) in enumerate(colors)
        ]
        # Make sure every class has its own surface.
        assert len(disp.surface_) == disp.n_classes

        for idx, quad in enumerate(disp.surface_):
            assert quad.cmap == cmaps[idx]
    else:
        assert_allclose(disp.surface_.colors, colors)

    # non-regression test for issue #32866 (currently still fails)
    # if hasattr(disp.surface_, "levels"):
    #    assert len(disp.surface_.levels) >= disp.n_classes


# estimator classes for non-regression test cases for issue #33194
class CustomBinaryEstimator(BaseEstimator):
    def fit(self, X, y):
        self.fitted_ = True
        return self

    def predict(self, X):
        return np.arange(X.shape[0]) % 2


class CustomMulticlassEstimator(BaseEstimator):
    def fit(self, X, y):
        self.fitted_ = True
        return self

    def predict(self, X):
        return np.arange(X.shape[0]) % 7


class CustomContinuousEstimator(BaseEstimator):
    def fit(self, X, y):
        self.fitted_ = True
        return self

    def predict(self, X):
        return np.arange(X.shape[0]) * 0.5


@pytest.mark.parametrize(
    "estimator, n_blobs, expected_n_classes",
    [
        (DecisionTreeClassifier(random_state=0), 7, 7),
        (DecisionTreeClassifier(random_state=0), 2, 2),
        (KMeans(n_clusters=7, random_state=0), 7, 7),
        (KMeans(n_clusters=2, random_state=0), 2, 2),
        (DecisionTreeRegressor(random_state=0), 7, 2),
        (IsolationForest(random_state=0), 7, 2),
        (CustomBinaryEstimator(), 2, 2),
        (CustomMulticlassEstimator(), 7, 7),
        (CustomContinuousEstimator(), 7, 2),
        (
            Pipeline(
                [
                    ("scale", StandardScaler()),
                    ("dt", DecisionTreeClassifier(random_state=0)),
                ]
            ),
            7,
            7,
        ),
        # non-regression test case for issue #33194
        (
            Pipeline(
                [
                    ("scale", StandardScaler()),
                    ("kmeans", KMeans(n_clusters=7, random_state=0)),
                ]
            ),
            7,
            7,
        ),
        (
            Pipeline(
                [
                    ("scale", StandardScaler()),
                    ("reg", DecisionTreeRegressor(random_state=0)),
                ]
            ),
            7,
            2,
        ),
        (
            Pipeline(
                [
                    ("scale", StandardScaler()),
                    ("kmeans", IsolationForest(random_state=0)),
                ]
            ),
            7,
            2,
        ),
    ],
)
def test_n_classes_attribute(pyplot, estimator, n_blobs, expected_n_classes):
    """Check that `n_classes` is set correctly.

    Introduced in https://github.com/scikit-learn/scikit-learn/pull/33015.
    """

    X, y = make_blobs(n_samples=150, centers=n_blobs, n_features=2, random_state=42)
    clf = estimator.fit(X, y)
    disp = DecisionBoundaryDisplay.from_estimator(clf, X, response_method="predict")
    assert disp.n_classes == expected_n_classes

    # Test that setting class_of_interest always converts to a binary problem.
    disp_coi = DecisionBoundaryDisplay.from_estimator(
        clf, X, class_of_interest=y[0], response_method="predict"
    )
    assert disp_coi.n_classes == 2


def test_n_classes_raises_if_not_inferrable(pyplot):
    """Check behaviour if `n_classes` can't be inferred.

    Non-regression test for issue #33194.
    """

    class CustomUnknownEstimator(BaseEstimator):
        def fit(self, X, y):
            self.fitted_ = True
            return self

        def predict(self, X):
            return np.array(0)

    X, y = load_iris_2d_scaled()
    est = CustomUnknownEstimator().fit(X, y)
    msg = "Number of classes or labels cannot be inferred from CustomUnknownEstimator"
    with pytest.raises(ValueError, match=msg):
        DecisionBoundaryDisplay.from_estimator(est, X, response_method="predict")


def test_cmap_and_colors_logic(pyplot):
    """Check the handling logic for `cmap` and `colors`."""
    X, y = load_iris_2d_scaled()
    clf = LogisticRegression().fit(X, y)

    with pytest.warns(
        UserWarning,
        match="'cmap' is ignored in favor of 'multiclass_colors'",
    ):
        DecisionBoundaryDisplay.from_estimator(
            clf,
            X,
            multiclass_colors="plasma",
            cmap="Blues",
        )

    with pytest.warns(
        UserWarning,
        match="'colors' is ignored in favor of 'multiclass_colors'",
    ):
        DecisionBoundaryDisplay.from_estimator(
            clf,
            X,
            multiclass_colors="plasma",
            colors="blue",
        )


def test_subclass_named_constructors_return_type_is_subclass(pyplot):
    """Check that named constructors return the correct type when subclassed.

    Non-regression test for:
    https://github.com/scikit-learn/scikit-learn/pull/27675
    """
    clf = LogisticRegression().fit(X, y)

    class SubclassOfDisplay(DecisionBoundaryDisplay):
        pass

    curve = SubclassOfDisplay.from_estimator(estimator=clf, X=X)

    assert isinstance(curve, SubclassOfDisplay)
