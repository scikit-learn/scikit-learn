import numpy as np
from ...utils import check_matplotlib_support


class DecisionBoundaryDisplay:
    def __init__(self, xx0, xx1, response):
        self.xx0 = xx0
        self.xx1 = xx1
        self.response = response

    def plot(self, ax=None, plot_method='contourf', **kwargs):
        check_matplotlib_support('DecisionBoundaryDisplay.plot')
        import matplotlib.pyplot as plt  # noqa

        if ax is None:
            _, ax = plt.subplots()

        plot_func = getattr(ax, plot_method)
        self.surface_ = plot_func(self.xx0, self.xx1,
                                  self.response, **kwargs)
        self.ax_ = ax
        self.figure_ = ax.figure
        return self


def plot_decision_boundary(est, X, grid_resolution=100,
                           features=(0, 1),
                           response_method='decision_function',
                           plot_method='contourf',
                           ax=None,
                           **kwargs):
    check_matplotlib_support('plot_decision_boundary')

    x0, x1 = X[:, features[0]], X[:, features[1]]

    x0_min, x0_max = x0.min() - 1, x0.max() + 1
    x1_min, x1_max = x1.min() - 1, x1.max() + 1

    xx0, xx1 = np.meshgrid(np.linspace(x0_min, x0_max, grid_resolution),
                           np.linspace(x1_min, x1_max, grid_resolution))

    response_func = getattr(est, response_method)
    response = response_func(np.c_[xx0.ravel(), xx1.ravel()])
    display = DecisionBoundaryDisplay(xx0=xx0, xx1=xx1,
                                      response=response.reshape(xx0.shape))

    return display.plot(ax=ax, plot_method=plot_method, **kwargs)
