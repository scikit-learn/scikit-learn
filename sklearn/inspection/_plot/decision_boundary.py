import numpy as np

from ...utils import check_matplotlib_support
from ...utils import _safe_indexing


def _check_boundary_response_method(estimator, response_method):
    """Return prediction method from the `response_method` for decision boundary.

    Parameters
    ----------
    estimator : object
        Estimator to check.

    response_method : {'auto', 'predict_proba', 'decision_function', 'predict'}
        Specifies whether to use :term:`predict_proba`,
        :term:`decision_function`, :term:`predict` as the target response.
        If set to 'auto', the response method is tried in the following order:
        :term:`predict_proba`, :term:`decision_function`, :term:`predict`.

    Returns
    -------
    prediction_method: callable
        Prediction method of estimator.
    """

    possible_response_methods = (
        "predict_proba",
        "decision_function",
        "auto",
        "predict",
    )
    if response_method not in possible_response_methods:
        raise ValueError(
            f"response_method must be one of {', '.join(possible_response_methods)}"
        )

    error_msg = "response method {} is not defined in {}"
    if response_method != "auto":
        if not hasattr(estimator, response_method):
            raise ValueError(
                error_msg.format(response_method, estimator.__class__.__name__)
            )
        return getattr(estimator, response_method)
    elif hasattr(estimator, "decision_function"):
        return getattr(estimator, "decision_function")
    elif hasattr(estimator, "predict_proba"):
        return getattr(estimator, "predict_proba")
    elif hasattr(estimator, "predict"):
        return getattr(estimator, "predict")

    raise ValueError(
        error_msg.format(
            "decision_function, predict_proba, or predict", estimator.__class__.__name__
        )
    )


class DecisionBoundaryDisplay:
    """Decisions boundary visualization.

    It is recommend to use
    :func:`~sklearn.inspection.DecisionBoundaryDisplay.from_estimator`
    to create a :class:`DecisionBoundaryDisplay`. All parameters are stored as
    attributes.

    Read more in the :ref:`User Guide <visualizations>`.

    .. versionadded:: 1.0

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
        Default label to place on y axis.

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

    def __init__(self, *, xx0, xx1, response, xlabel=None, ylabel=None):
        self.xx0 = xx0
        self.xx1 = xx1
        self.response = response
        self.xlabel = xlabel
        self.ylabel = ylabel

    def plot(self, plot_method="contourf", ax=None, xlabel=None, ylabel=None, **kwargs):
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

        xlabel : str, default=None
            Overwrite the x-axis label.

        ylabel : str, default=None
            Overwrite the y-axis label.

        **kwargs : dict
            Additional keyword arguments to be passed to the `plot_method`.

        Returns
        -------
        display: :class:`~sklearn.inspection.DecisionBoundaryDisplay`
        """
        check_matplotlib_support("DecisionBoundaryDisplay.plot")
        import matplotlib.pyplot as plt  # noqa

        if plot_method not in ("contourf", "contour", "pcolormesh"):
            raise ValueError(
                "plot_method must be 'contourf', 'contour', or 'pcolormesh'"
            )

        if ax is None:
            _, ax = plt.subplots()

        plot_func = getattr(ax, plot_method)
        self.surface_ = plot_func(self.xx0, self.xx1, self.response, **kwargs)

        if xlabel is not None or not ax.get_xlabel():
            xlabel = self.xlabel if xlabel is None else xlabel
            ax.set_xlabel(xlabel)
        if ylabel is not None or not ax.get_ylabel():
            ylabel = self.ylabel if ylabel is None else ylabel
            ax.set_ylabel(ylabel)

        self.ax_ = ax
        self.figure_ = ax.figure
        return self

    @classmethod
    def from_estimator(
        cls,
        estimator,
        X,
        *,
        grid_resolution=100,
        eps=1.0,
        plot_method="contourf",
        response_method="auto",
        xlabel=None,
        ylabel=None,
        ax=None,
        **kwargs,
    ):
        """Plot decision boundary given an estimator.

        Read more in the :ref:`User Guide <visualizations>`.

        .. versionadded:: 1.0

        Parameters
        ----------
        estimator : object
            Trained estimator used to plot the decision boundary.

        X : {array-like, sparse matrix, dataframe} of shape (n_samples, 2)
            Input data that should be only 2-dimensional.

        grid_resolution : int, default=100
            Number of grid points to use for plotting decision boundary.
            Higher values will make the plot look nicer but be slower to
            render.

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
                'predict'}, default='auto'
            Specifies whether to use :term:`predict_proba`,
            :term:`decision_function`, :term:`predict` as the target response.
            If set to 'auto', the response method is tried in the following order:
            :term:`predict_proba`, :term:`decision_function`, :term:`predict`.

        xlabel : str, default=None
            The label used for the x-axis. If `None`, an attempt is made to
            extract a label from `X` if it is a dataframe, otherwise an empty
            string is used.

        ylabel : str, default=None
            The label used for the y-axis. If `None`, an attempt is made to
            extract a label from `X` if it is a dataframe, otherwise an empty
            string is used.

        ax : Matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        **kwargs : dict
            Additional keyword arguments to be passed to the
            `plot_method`.

        Returns
        -------
        display : :class:`~sklearn.inspection.DecisionBoundaryDisplay`
            Object that stores the result.

        See Also
        --------
        DecisionBoundaryDisplay : Decision boundary visualization.
        ConfusionMatrixDisplay.from_estimator : Plot the confusion matrix
            given an estimator, the data, and the label.
        ConfusionMatrixDisplay.from_predictions : Plot the confusion matrix
            given the true and predicted labels.

        Examples
        --------
        >>> import matplotlib.pyplot as plt
        >>> from sklearn.datasets import load_iris
        >>> from sklearn.linear_model import LogisticRegression
        >>> from sklearn.inspection import DecisionBoundaryDisplay
        >>> iris = load_iris()
        >>> X = iris.data[:, :2]
        >>> classifier = LogisticRegression().fit(X, iris.target)
        >>> disp = DecisionBoundaryDisplay.from_estimator(
        ...     classifier, X, response_method="predict",
        ...     xlabel=iris.feature_names[0], ylabel=iris.feature_names[1],
        ...     alpha=0.5,
        ... )
        >>> disp.ax_.scatter(X[:, 0], X[:, 1], c=iris.target, edgecolor="k")
        <...>
        >>> plt.show()
        """
        check_matplotlib_support(f"{cls.__name__}.from_estimator")

        if not grid_resolution > 1:
            raise ValueError(
                "grid_resolution must be greater than 1. Got"
                f" {grid_resolution} instead."
            )

        if not eps >= 0:
            raise ValueError(
                f"eps must be greater than or equal to 0. Got {eps} instead."
            )

        possible_plot_methods = ("contourf", "contour", "pcolormesh")
        if plot_method not in possible_plot_methods:
            avaliable_methods = ", ".join(possible_plot_methods)
            raise ValueError(
                f"plot_method must be one of {avaliable_methods}. "
                f"Got {plot_method} instead."
            )

        x0, x1 = _safe_indexing(X, 0, axis=1), _safe_indexing(X, 1, axis=1)

        x0_min, x0_max = x0.min() - eps, x0.max() + eps
        x1_min, x1_max = x1.min() - eps, x1.max() + eps

        xx0, xx1 = np.meshgrid(
            np.linspace(x0_min, x0_max, grid_resolution),
            np.linspace(x1_min, x1_max, grid_resolution),
        )

        pred_func = _check_boundary_response_method(estimator, response_method)
        response = pred_func(np.c_[xx0.ravel(), xx1.ravel()])

        if response.ndim != 1:
            if response.shape[1] != 2:
                raise ValueError(
                    "Multiclass classifiers are only supported when "
                    "response_method='predict'"
                )
            response = response[:, 1]

        if xlabel is not None:
            xlabel = xlabel
        else:
            xlabel = X.columns[0] if hasattr(X, "columns") else ""

        if ylabel is not None:
            ylabel = ylabel
        else:
            ylabel = X.columns[1] if hasattr(X, "columns") else ""

        display = DecisionBoundaryDisplay(
            xx0=xx0,
            xx1=xx1,
            response=response.reshape(xx0.shape),
            xlabel=xlabel,
            ylabel=ylabel,
        )
        return display.plot(ax=ax, plot_method=plot_method, **kwargs)
