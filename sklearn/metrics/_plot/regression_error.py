import numbers

import numpy as np

from ...utils import check_matplotlib_support
from ...utils import check_random_state


class PredictionErrorDisplay:
    """Prediction error visualization.

    It is recommended to use :func:`~sklearn.metrics.plot_prediction_error` to
    create a visualizer. All parameters are stored as attributes.

    Read more in the :ref:`User Guide <visualizations>`.

    Parameters
    ----------

    Attributes
    ----------

    Examples
    --------
    """

    def __init__(
        self,
        *,
        y_true,
        y_pred,
        subsample=1000,
        with_residuals=False,
        random_state=None,
    ):
        self.y_true = y_true
        self.y_pred = y_pred
        self.subsample = subsample
        self.with_residuals = with_residuals
        self.random_state = random_state

    def plot(
        self,
        ax=None,
        *,
        scatter_kwargs=None,
        line_kwargs=None,
        residuals_kwargs=None,
    ):
        """Plot visualization.

        Extra keyword arguments will be passed to matplotlib's ``plot``.

        Parameters
        ----------
        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        Returns
        -------
        display : :class:`~sklearn.metrics.plot.PredictionErrorDisplay`
            Object that stores computed values.
        """
        check_matplotlib_support("PredictionErrorDisplay.plot")

        import matplotlib.pyplot as plt

        random_state = check_random_state(self.random_state)

        subsample = self.subsample
        if isinstance(subsample, numbers.Integral):
            if subsample <= 0:
                raise ValueError(
                    f"When an integer, subsample={subsample} should be "
                    f"positive."
                )
        elif isinstance(subsample, numbers.Real):
            if subsample <= 0 or subsample >= 1:
                raise ValueError(
                    f"When a floating-point, subsample={subsample} should"
                    f" be in the (0, 1) range."
                )
            subsample = int(len(self.y_true) * subsample)

        if subsample is not None and subsample < len(self.y_true):
            indices = random_state.choice(
                np.arange(len(self.y_true)), size=subsample
            )
            y_true, y_pred = self.y_true[indices], self.y_pred[indices]
        else:
            y_true, y_pred = self.y_true, self.y_pred

        if scatter_kwargs is None:
            scatter_kwargs = {}
        if line_kwargs is None:
            line_kwargs = {}
        if residuals_kwargs is None:
            residuals_kwargs = {}

        default_scatter_kwargs = {"color": "black", "alpha": 0.8}
        default_line_kwargs = {
            "color": "tab:blue",
            "alpha": 0.7,
            "linestyle": "--",
        }
        default_residuals_kwargs = {"color": "red", "alpha": 0.5}

        scatter_kwargs = {**default_scatter_kwargs, **scatter_kwargs}
        line_kwargs = {**default_line_kwargs, **line_kwargs}
        residuals_kwargs = {**default_residuals_kwargs, **residuals_kwargs}

        if ax is None:
            _, ax = plt.subplots()

        if self.with_residuals:
            for actual, predicted in zip(y_true, y_pred):
                residual_line = ax.plot(
                    [actual, actual], [actual, predicted], **residuals_kwargs,
                )

        max_value = max(np.max(y_true), np.max(y_pred))
        min_value = min(np.min(y_true), np.min(y_pred))
        perfect_fit_line = ax.plot(
            [min_value, max_value], [min_value, max_value], **line_kwargs
        )

        ax.scatter(y_true, y_pred, **scatter_kwargs)

        xlabel, ylabel = "Actual values", "Predicted values"
        ax.set(xlabel=xlabel, ylabel=ylabel)
        ax.axis("square")
        ax.set_xticks(np.linspace(min_value, max_value, num=5))
        ax.set_yticks(np.linspace(min_value, max_value, num=5))

        legend_obj = [perfect_fit_line[0]]
        legend_name = ["Perfect fit"]
        if self.with_residuals:
            legend_obj.append(residual_line[0])
            legend_name.append("Residuals")
        ax.legend(legend_obj, legend_name)

        self.ax_ = ax
        self.figure_ = ax.figure

        return self


def plot_prediction_error(
    estimator,
    X,
    y,
    *,
    subsample=None,
    with_residuals=False,
    ax=None,
    scatter_kwargs=None,
    line_kwargs=None,
    residuals_kwargs=None,
):
    """Plot the prediction error of a regressor.

    Extra keyword arguments will be passed to matplotlib's ``plot``.

    Read more in the :ref:`User Guide <visualizations>`.

    Parameters
    ----------
    estimator : estimator instance
        Fitted regressor or a fitted :class:`~sklearn.pipeline.Pipeline`
        in which the last estimator is a regressor.

    X : {array-like, sparse matrix} of shape (n_samples, n_features)
        Input values.

    y : array-like of shape (n_samples,)
        Target values.

    subsample : float, int or None, default=1000
        Sampling the samples to be shown on the scatter plot. If `float`, it
        should be between 0 and 1 and represents the proportion of the original
        dataset. If `int`, it represents the number of samples display on the
        scatter plot. If `None`, no subsampling will be applied. by default,
        a 1000 samples or less will be displayed.

    with_residuals : bool, default=False
        Whether or not to display the residuals on the plot for each sample.

    ax : matplotlib axes, default=None
        Axes object to plot on. If `None`, a new figure and axes is created.

    Returns
    -------
    display : :class:`~sklearn.metrics.PredictionErrorDisplay`
        Object that stores the computed values.

    Examples
    --------
    """
    check_matplotlib_support("plot_prediction_error")

    y_pred = estimator.predict(X)

    viz = PredictionErrorDisplay(
        y_true=y,
        y_pred=y_pred,
        subsample=subsample,
        with_residuals=with_residuals,
    )

    return viz.plot(
        ax=ax,
        scatter_kwargs=scatter_kwargs,
        line_kwargs=line_kwargs,
        residuals_kwargs=residuals_kwargs,
    )
