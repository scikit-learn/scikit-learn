import numpy as np
from ...utils import check_matplotlib_support
from ...utils import _safe_indexing


def _check_boundary_response_method(estimator, response_method):
    """Return prediction method from the response_method for decision boundary

    Parameters
    ----------
    estimator: object
        Estimator to check

    response_method: {'auto', 'predict_proba', 'decision_function', 'predict'}
        Specifies whether to use :term:`predict_proba`,
        :term:`decision_function`, :term:`predict` as the target response.
        If set to 'auto', the response method is tried in the following order:
        :term:`predict_proba`, :term:`decision_function`, :term:`predict`.

    Returns
    -------
    prediction_method: callable
        prediction method of estimator
    """

    if response_method not in ("predict_proba", "decision_function",
                               "auto", "predict"):
        raise ValueError("response_method must be 'predict_proba', "
                         "'decision_function', 'predict', or 'auto'")

    error_msg = "response method {} is not defined in {}"
    if response_method != "auto":
        if not hasattr(estimator, response_method):
            raise ValueError(error_msg.format(response_method,
                                              estimator.__class__.__name__))
        return getattr(estimator, response_method)
    elif hasattr(estimator, 'decision_function'):
        return getattr(estimator, 'decision_function')
    elif hasattr(estimator, 'predict_proba'):
        return getattr(estimator, 'predict_proba')
    elif hasattr(estimator, 'predict'):
        return getattr(estimator, 'predict')

    raise ValueError(error_msg.format(
        "decision_function, predict_proba, or predict",
        estimator.__class__.__name__))


class DecisionBoundaryDisplay:
    """Decisions Boundary visualization.

    It is recommend to use :func:`~sklearn.inspection.plot_decision_boundary`
    to create a :class:`DecisionBoundaryDisplay`. All parameters are stored as
    attributes.

    Read more in the :ref:`User Guide <visualizations>`.

    Parameters
    ----------
    xx0 : ndarray of shape (grid_resolution, grid_resolution)
        First output of :func:`meshgrid <numpy.meshgrid>`.

    xx1 : ndarray of shape (grid_resolution, grid_resolution)
        Second output of :func:`meshgrid <numpy.meshgrid>`.

    response : ndarray of shape (grid_resolution, grid_resolution)
        Values of the response function.

    xlabel : str, default=""
        Default label to place on x axis.

    ylabel : str, default=""
        DEfault label to place on y axis.

    Attributes
    ----------
    surface_ : matplotlib `QuadContourSet` or `QuadMesh`
        If `plot_method` is 'contour' or 'contourf', `surface_` is a
        :class:`QuadContourSet <matplotlib.contour.QuadContourSet>`. If
        `plot_method is `pcolormesh`, `surface_` is a
        :class:`QuadMesh <matplotlib.collections.QuadMesh>`.

    ax_ : matplotlib Axes
        Axes with confusion matrix.

    figure_ : matplotlib Figure
        Figure containing the confusion matrix.
    """
    def __init__(self, xx0, xx1, response, xlabel=None, ylabel=None):
        self.xx0 = xx0
        self.xx1 = xx1
        self.response = response
        self.xlabel = xlabel
        self.ylabel = ylabel

    def plot(self, plot_method='contourf', ax=None, **kwargs):
        """Plot visualization.

        Parameters
        ----------
        plot_method : {'contourf', 'contour', 'pcolormesh'}, default='contourf'
            Plotting method to call when plotting the response. Please refer
            to the following matplotlib documentation for details:
            :func:`contourf <matplotlib.pyplot.contourf>`,
            :func:`contour <matplotlib.pyplot.contour>`,
            :func:`pcolomesh <matplotlib.pyplot.pcolomesh>`.

        ax : Matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        **kwargs : dict
            Additional keyword arguments to be pased to the `plot_method`.

        Returns
        -------
        display: :class:`~sklearn.inspection.DecisionBoundaryDisplay`
        """
        check_matplotlib_support('DecisionBoundaryDisplay.plot')
        import matplotlib.pyplot as plt  # noqa

        if plot_method not in ('contourf', 'contour', 'pcolormesh'):
            raise ValueError("plot_method must be 'contourf', 'contour', or "
                             "'pcolormesh'")

        if ax is None:
            _, ax = plt.subplots()

        plot_func = getattr(ax, plot_method)
        self.surface_ = plot_func(self.xx0, self.xx1, self.response, **kwargs)

        if not ax.get_xlabel():
            ax.set_xlabel(self.xlabel)
        if not ax.get_ylabel():
            ax.set_ylabel(self.ylabel)

        self.ax_ = ax
        self.figure_ = ax.figure
        return self


def plot_decision_boundary(estimator, X, grid_resolution=100, eps=1.0,
                           plot_method='contourf', response_method='auto',
                           ax=None, **kwargs):
    """Plot Decision Boundary.

    Please see examples below for usage.

    Parameters
    ----------
    estimator : estimator instance
        Trained estimator.

    X : ndarray or pandas dataframe of shape (n_samples, 2)
        Input values.

    grid_resolution : int, default=100
        The number of equally spaced points to evaluate the response function.

    eps : float, default=1.0
        Extends the minimum and maximum values of X for evaluating the
        response function.

    plot_method : {'contourf', 'contour', 'pcolormesh'}, default='contourf'
        Plotting method to call when plotting the response. Please refer
        to the following matplotlib documentation for details:
        :func:`contourf <matplotlib.pyplot.contourf>`,
        :func:`contour <matplotlib.pyplot.contour>`,
        :func:`pcolomesh <matplotlib.pyplot.pcolomesh>`.

    response_method : {'auto', 'predict_proba', 'decision_function', \
            'predict'}, defaul='auto'
        Specifies whether to use :term:`predict_proba`,
        :term:`decision_function`, :term:`predict` as the target response.
        If set to 'auto', the response method is tried in the following order:
        :term:`predict_proba`, :term:`decision_function`, :term:`predict`.

    ax : Matplotlib axes, default=None
        Axes object to plot on. If `None`, a new figure and axes is
        created.

    **kwargs : dict
        Additional keyword arguments to be pased to the `plot_method`.

    Returns
    -------
    display: :class:`~sklearn.inspection.DecisionBoundaryDisplay`
    """
    check_matplotlib_support('plot_decision_boundary')

    if not grid_resolution > 1:
        raise ValueError("grid_resolution must be greater than 1")

    if not eps >= 0:
        raise ValueError("eps must be greater than or equal to 0")

    if plot_method not in ('contourf', 'contour', 'pcolormesh'):
        raise ValueError("plot_method must be 'contourf', 'contour', or "
                         "'pcolormesh'")

    pred_func = _check_boundary_response_method(estimator, response_method)

    x0, x1 = _safe_indexing(X, 0, axis=1), _safe_indexing(X, 1, axis=1)

    x0_min, x0_max = x0.min() - eps, x0.max() + eps
    x1_min, x1_max = x1.min() - eps, x1.max() + eps

    xx0, xx1 = np.meshgrid(np.linspace(x0_min, x0_max, grid_resolution),
                           np.linspace(x1_min, x1_max, grid_resolution))
    response = pred_func(np.c_[xx0.ravel(), xx1.ravel()])

    if response.ndim != 1:
        if response.shape[1] != 2:
            raise ValueError("multiclass classifiers are only supported when "
                             "response_method='predict'")
        response = response[:, 1]

    if hasattr(X, "columns"):
        xlabel, ylabel = X.columns[0], X.columns[1]
    else:
        xlabel, ylabel = "", ""

    display = DecisionBoundaryDisplay(xx0=xx0, xx1=xx1,
                                      response=response.reshape(xx0.shape),
                                      xlabel=xlabel, ylabel=ylabel)
    return display.plot(ax=ax, plot_method=plot_method, **kwargs)
