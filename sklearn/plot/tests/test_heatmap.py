from sklearn.plot import plot_heatmap
from sklearn.plot import plot_confusion_matrix
from sklearn.plot import plot_gridsearch_results
from sklearn.datasets import make_blobs
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from sklearn.utils.testing import SkipTest
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raise_message
import numpy as np
from numpy.random import (RandomState,
                          randint)


def test_heatmap():
    try:
        import matplotlib
    except ImportError:
        raise SkipTest("Not testing plot_heatmap, matplotlib not installed.")

    import matplotlib.pyplot as plt

    with matplotlib.rc_context(rc={'backend': 'Agg', 'interactive': False}):
        plt.figure()
        rng = RandomState(0)
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
        array1 = randint(1, 3, size=20)
        array2 = randint(1, 3, size=20)

        # plot un-normalized matrix
        plot_confusion_matrix(array1, array2, classes=["dummy1", "dummy2"],
                              cmap="Paired", ax=plt.gca())

        # plot normalized matrix
        plot_confusion_matrix(array1, array2, normalize=True,
                              classes=["dummay1", "dummy2"],
                              cmap="Paired", ax=plt.gca())

        # plot without passing classes explicitly
        plot_confusion_matrix(array1, array2,
                              cmap="Paired", ax=plt.gca())

        # y having different value than classes should raise error
        expected_msg = ("`classes=[1,2]` are not a superset of the unique",
                        "values of y_true and y_pred which are [1,2,3]")
        assert_raise_message(ValueError, expected_msg,
                             plot_confusion_matrix, array1, array2,
                             classes=[1, 2], ax=plt.gca())

        plt.draw()
        plt.close()


def test_gridsearch_results_1d():
    try:
        import matplotlib
    except ImportError:
        raise SkipTest("Not testing plot_heatmap, matplotlib not installed.")

    import matplotlib.pyplot as plt

    X, y = make_blobs(n_samples=20, centers=2, n_features=3,
                      random_state=0)

    # Define parameters:
    C_range = np.logspace(-2, 10, 2)

    # Test 1D case:
    param_grid = dict(C=C_range)
    grid = GridSearchCV(SVC(), param_grid=param_grid, cv=3)
    grid.fit(X, y)

    with matplotlib.rc_context(rc={'backend': 'Agg', 'interactive': False}):
        plt.figure()
        plot_gridsearch_results(grid.cv_results_)
        plt.draw()
        plt.close()


def test_gridsearch_results_2d():
    try:
        import matplotlib
    except ImportError:
        raise SkipTest("Not testing plot_heatmap, matplotlib not installed.")

    import matplotlib.pyplot as plt

    X, y = make_blobs(n_samples=20, centers=2, n_features=3,
                      random_state=0)

    # Define parameters:
    C_range = np.logspace(-2, 10, 2)
    gamma_range = np.logspace(-9, 3, 2)

    # Test 1D case:
    param_grid = dict(gamma=gamma_range, C=C_range)
    grid = GridSearchCV(SVC(), param_grid=param_grid, cv=3)
    grid.fit(X, y)

    with matplotlib.rc_context(rc={'backend': 'Agg', 'interactive': False}):
        plt.figure()
        plot_gridsearch_results(grid.cv_results_)
        plt.draw()
        plt.close()


def test_gridsearch_results_3d():
    try:
        import matplotlib
    except ImportError:
        raise SkipTest("Not testing plot_heatmap, matplotlib not installed.")

    import matplotlib.pyplot as plt

    X, y = make_blobs(n_samples=20, centers=2, n_features=3,
                      random_state=0)

    # Define parameters:
    C_range = np.logspace(-2, 10, 2)
    gamma_range = np.logspace(-9, 3, 2)
    tol_range = [1e-3, 1e-4]

    # Test 1D case:
    param_grid = dict(gamma=gamma_range, C=C_range, tol=tol_range)
    grid = GridSearchCV(SVC(), param_grid=param_grid, cv=3)
    grid.fit(X, y)

    with matplotlib.rc_context(rc={'backend': 'Agg', 'interactive': False}):
        plt.figure()
        assert_raises(ValueError, plot_gridsearch_results, grid.cv_results_)
        plt.draw()
        plt.close()
