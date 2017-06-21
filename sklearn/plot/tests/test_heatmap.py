from sklearn.plot import plot_heatmap
from sklearn.plot import plot_confusion_matrix
from sklearn.plot import plot_gridsearch_results
from sklearn.datasets import load_iris
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from sklearn.utils.testing import SkipTest
from sklearn.utils.testing import assert_raises
import numpy as np


def test_heatmap():
    try:
        import matplotlib
    except ImportError:
        raise SkipTest("Not testing plot_heatmap, matplotlib not installed.")

    import matplotlib.pyplot as plt

    with matplotlib.rc_context(rc={'backend': 'Agg', 'interactive': False}):
        plt.figure()
        rng = np.random.RandomState(0)
        X = rng.normal(size=(10, 5))
        # use mixture of default values and keyword args
        plot_heatmap(X, ylabel="y-axis",
                     xticklabels=["a", "b", "c", "d", "efgh"],
                     cmap="Paired", ax=plt.gca())

        plt.draw()
        plt.close()


def test_confusion_matrix():
    try:
        import matplotlib
    except ImportError:
        raise SkipTest("Not testing plot_heatmap, matplotlib not installed.")

    import matplotlib.pyplot as plt

    with matplotlib.rc_context(rc={'backend': 'Agg', 'interactive': False}):
        plt.figure()
        cnf_matrix = np.random.randomint(1, 10, size=(2, 2))

        # use mixture of default values and keyword args
        plot_confusion_matrix(cnf_matrix, classes=["dummay1", "dummy2"],
                              cmap="Paired", ax=plt.gca())

        plt.draw()
        plt.close()


def test_gridsearch_results():
    try:
        import matplotlib
    except ImportError:
        raise SkipTest("Not testing plot_heatmap, matplotlib not installed.")

    import matplotlib.pyplot as plt

    iris = load_iris()
    X = iris.data
    y = iris.target

    # We only keep the first two features in X and sub-sample the dataset to
    # keep only 2 classes and make it a binary classification problem.

    X_2d = X[:, :2]
    X_2d = X_2d[y > 0]
    y_2d = y[y > 0]
    y_2d -= 1

    # Define parameters:
    C_range = np.logspace(-2, 10, 2)
    gamma_range = np.logspace(-9, 3, 2)
    tol_range = [1e-3, 1e-4]

    # Test 1D case:
    param_grid = dict(gamma=gamma_range)
    grid = GridSearchCV(SVC(), param_grid=param_grid, cv=3)
    grid.fit(X, y)

    with matplotlib.rc_context(rc={'backend': 'Agg', 'interactive': False}):
        plt.figure()
        plot_gridsearch_results(grid.cv_results_)
        plt.draw()
        plt.close()

    # Test 2D case:
    param_grid = dict(gamma=gamma_range, C=C_range)
    grid = GridSearchCV(SVC(), param_grid=param_grid, cv=3)
    grid.fit(X, y)

    with matplotlib.rc_context(rc={'backend': 'Agg', 'interactive': False}):
        plt.figure()
        plot_gridsearch_results(grid.cv_results_)
        plt.draw()
        plt.close()

    # Test 3D case:
    param_grid = dict(gamma=gamma_range, C=C_range, tol=tol_range)
    grid = GridSearchCV(SVC(), param_grid=param_grid, cv=3)
    grid.fit(X, y)

    with matplotlib.rc_context(rc={'backend': 'Agg', 'interactive': False}):
        plt.figure()
        assert_raises(ValueError, plot_gridsearch_results, grid.cv_results_)
        plt.draw()
        plt.close()
