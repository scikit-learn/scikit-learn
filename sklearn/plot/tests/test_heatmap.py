from sklearn.plot import plot_heatmap
from sklearn.utils.testing import SkipTest
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
