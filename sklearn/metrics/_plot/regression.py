import numbers

import numpy as np

from ...utils import check_matplotlib_support
from ...utils import check_random_state
from ...utils import _safe_indexing


class PredictionErrorDisplay:
    """Prediction error visualization.

    It is recommended to use
    :func:`~sklearn.metrics.PredictionErrorDisplay.from_estimator` or
    :func:`~sklearn.metrics.PredictionErrorDisplay.from_predictions` to
    create a visualizer. All parameters are stored as attributes.

    Read more in the :ref:`User Guide <visualizations>`.

    .. versionadded:: 1.0

    Parameters
    ----------
    y_true : ndarray of shape (n_samples,)
        True values.

    y_pred : ndarray of shape (n_samples,)
        Prediction values.

    scores : dict, default=None
        Dictionary where the key is the name of the metric displayed and the
        value is the metric value.

    with_residuals : bool, default=False
        Whether or not to display the residuals on the plot for each sample.

    Attributes
    ----------
    line_ : matplotlib Artist
        Diagonal curve.

    residual_lines_ : matplotlib Artist or None
        Residual lines. If `with_residuals=False`, then it is set to `None`.

    scatter : matplotlib Artist
        Scatter data points.

    ax_ : matplotlib Axes
        Axes with the different matplotlib axis.

    figure_ : matplotlib Figure
        Figure containing the scatter and lines.

    See Also
    --------
    PredictionDisplay.from_estimator : Prediction error visualization
        given an estimator and some data.
    PredictionDisplay.from_predictions : Prediction error visualization
        given the true and predicted targets.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from sklearn.datasets import load_diabetes
    >>> from sklearn.linear_model import Ridge
    >>> from sklearn.metrics import PredictionErrorDisplay
    >>> X, y = load_diabetes(return_X_y=True)
    >>> ridge = Ridge().fit(X, y)
    >>> y_pred = ridge.predict(X)
    >>> display = PredictionErrorDisplay(y_true=y, y_pred=y_pred)
    >>> display.plot()
    <...>
    >>> plt.show()
    """

    def __init__(
        self,
        *,
        y_true,
        y_pred,
        scores=None,
        with_residuals=False,
    ):
        self.y_true = y_true
        self.y_pred = y_pred
        self.scores = scores
        self.with_residuals = with_residuals

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

        scatter_kwargs : dict, default=None
            Dictionary with keywords passed to the `matplotlib.pyplot.scatter`
            call.

        line_kwargs : dict, default=None
            Dictionary with keyword passed to the `matplotlib.pyplot.plot`
            call to draw the diagonal line.

        residuals_kwargs : dict, default=None
            Dictionary with keyword passed to the `matplotlib.pyplot.plot`
            call to draw the residual lines. Only taken into account when
            `with_residuals=True`.

        Returns
        -------
        display : :class:`~sklearn.metrics.plot.PredictionErrorDisplay`
            Object that stores computed values.
        """
        check_matplotlib_support("PredictionErrorDisplay.plot")

        import matplotlib.pyplot as plt

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
            self.residual_lines_ = []
            for actual, predicted in zip(self.y_true, self.y_pred):
                residual_line = ax.plot(
                    [actual, actual],
                    [actual, predicted],
                    **residuals_kwargs,
                )
                self.residual_lines_ += residual_line
        else:
            self.residual_lines_ = None

        max_value = max(np.max(self.y_true), np.max(self.y_pred))
        min_value = min(np.min(self.y_true), np.min(self.y_pred))
        self.line_ = ax.plot(
            [min_value, max_value], [min_value, max_value], **line_kwargs
        )[0]

        self.scatter_ = ax.scatter(self.y_true, self.y_pred, **scatter_kwargs)

        xlabel, ylabel = "Actual values", "Predicted values"
        ax.set(xlabel=xlabel, ylabel=ylabel)
        ax.set_aspect("equal", adjustable="datalim")
        ax.set_xticks(np.linspace(min_value, max_value, num=5))
        ax.set_yticks(np.linspace(min_value, max_value, num=5))

        if self.scores is not None:
            extra = plt.Rectangle(
                (0, 0), 0, 0, fc="w", fill=False, edgecolor="none", linewidth=0
            )
            scoring_legend = "\n".join(
                [f"{name} = {value}" for name, value in self.scores.items()],
            )
            ax.legend([extra], [scoring_legend])

        self.ax_ = ax
        self.figure_ = ax.figure

        return self

    @classmethod
    def from_estimator(
        cls,
        estimator,
        X,
        y,
        *,
        scores=None,
        subsample=1_000,
        random_state=None,
        with_residuals=False,
        ax=None,
        scatter_kwargs=None,
        line_kwargs=None,
        residuals_kwargs=None,
    ):
        """Plot the prediction error given a regressor and some data.

        Read more in the :ref:`User Guide <visualizations>`.

        .. versionadded:: 1.0

        Parameters
        ----------
        estimator : estimator instance
            Fitted regressor or a fitted :class:`~sklearn.pipeline.Pipeline`
            in which the last estimator is a regressor.

        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Input values.

        y : array-like of shape (n_samples,)
            Target values.

        scores : dict, default=None
            Dictionary containing scores that will be shown on the legend.
            Key are the score names, values are the score values.

        subsample : float, int or None, default=1_000
            Sampling the samples to be shown on the scatter plot. If `float`,
            it should be between 0 and 1 and represents the proportion of the
            original dataset. If `int`, it represents the number of samples
            display on the scatter plot. If `None`, no subsampling will be
            applied. by default, a 1000 samples or less will be displayed.

        random_state : int or RandomState, default=None
            Controls the randomness when `subsample` is not `None`.
            See :term:`Glossary <random_state>` for details.

        with_residuals : bool, default=False
            Whether or not to display the residuals on the plot for each
            sample.

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        scatter_kwargs : dict, default=None
            Dictionary with keywords passed to the `matplotlib.pyplot.scatter`
            call.

        line_kwargs : dict, default=None
            Dictionary with keyword passed to the `matplotlib.pyplot.plot`
            call to draw the diagonal line.

        residuals_kwargs : dict, default=None
            Dictionary with keyword passed to the `matplotlib.pyplot.plot`
            call to draw the residual lines. Only taken into account when
            `with_residuals=True`.

        Returns
        -------
        display : :class:`~sklearn.metrics.PredictionErrorDisplay`
            Object that stores the computed values.

        See Also
        --------
        PredictionErrorDisplay : Prediction error visualization for regression.
        PredictionDisplay.from_predictions : Prediction error visualization
            given the true and predicted targets.

        Examples
        --------
        >>> import matplotlib.pyplot as plt
        >>> from sklearn.datasets import load_diabetes
        >>> from sklearn.linear_model import Ridge
        >>> from sklearn.metrics import PredictionErrorDisplay
        >>> X, y = load_diabetes(return_X_y=True)
        >>> ridge = Ridge().fit(X, y)
        >>> disp = PredictionErrorDisplay.from_estimator(ridge, X, y)
        >>> plt.show()
        """
        check_matplotlib_support(f"{cls.__name__}.from_estimator")

        y_pred = estimator.predict(X)

        return cls.from_predictions(
            y_true=y,
            y_pred=y_pred,
            scores=scores,
            subsample=subsample,
            random_state=random_state,
            with_residuals=with_residuals,
            ax=ax,
            scatter_kwargs=scatter_kwargs,
            line_kwargs=line_kwargs,
            residuals_kwargs=residuals_kwargs,
        )

    @classmethod
    def from_predictions(
        cls,
        y_true,
        y_pred,
        *,
        scores=None,
        subsample=1_000,
        random_state=None,
        with_residuals=False,
        ax=None,
        scatter_kwargs=None,
        line_kwargs=None,
        residuals_kwargs=None,
    ):
        """Plot the prediction error given the true and predicted targets.

        Read more in the :ref:`User Guide <visualizations>`.

        .. versionadded:: 1.0

        Parameters
        ----------
        y_true : array-like of shape (n_samples,)
            True target values.

        y_pred : array-like of shape (n_samples,)
            Predicted target values.

        scores : dict, default=None
            Dictionary containing scores that will be shown on the legend.
            Key are the score names, values are the score values.

        subsample : float, int or None, default=1_000
            Sampling the samples to be shown on the scatter plot. If `float`,
            it should be between 0 and 1 and represents the proportion of the
            original dataset. If `int`, it represents the number of samples
            display on the scatter plot. If `None`, no subsampling will be
            applied. by default, a 1000 samples or less will be displayed.

        random_state : int or RandomState, default=None
            Controls the randomness when `subsample` is not `None`.
            See :term:`Glossary <random_state>` for details.

        with_residuals : bool, default=False
            Whether or not to display the residuals on the plot for each
            sample.

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        scatter_kwargs : dict, default=None
            Dictionary with keywords passed to the `matplotlib.pyplot.scatter`
            call.

        line_kwargs : dict, default=None
            Dictionary with keyword passed to the `matplotlib.pyplot.plot`
            call to draw the diagonal line.

        residuals_kwargs : dict, default=None
            Dictionary with keyword passed to the `matplotlib.pyplot.plot`
            call to draw the residual lines. Only taken into account when
            `with_residuals=True`.

        Returns
        -------
        display : :class:`~sklearn.metrics.PredictionErrorDisplay`
            Object that stores the computed values.

        See Also
        --------
        PredictionErrorDisplay : Prediction error visualization for regression.
        PredictionDisplay.from_estimator : Prediction error visualization
            given an estimator and some data.

        Examples
        --------
        >>> import matplotlib.pyplot as plt
        >>> from sklearn.datasets import load_diabetes
        >>> from sklearn.linear_model import Ridge
        >>> from sklearn.metrics import PredictionErrorDisplay
        >>> X, y = load_diabetes(return_X_y=True)
        >>> ridge = Ridge().fit(X, y)
        >>> y_pred = ridge.predict(X)
        >>> disp = PredictionErrorDisplay.from_predictions(y_true=y, y_pred=y_pred)
        >>> plt.show()
        """
        check_matplotlib_support(f"{cls.__name__}.from_predictions")

        random_state = check_random_state(random_state)

        n_samples = len(y_true)
        if isinstance(subsample, numbers.Integral):
            if subsample <= 0:
                raise ValueError(
                    f"When an integer, subsample={subsample} should be positive."
                )
        elif isinstance(subsample, numbers.Real):
            if subsample <= 0 or subsample >= 1:
                raise ValueError(
                    f"When a floating-point, subsample={subsample} should"
                    " be in the (0, 1) range."
                )
            subsample = int(n_samples * subsample)

        if subsample is not None and subsample < n_samples:
            indices = random_state.choice(np.arange(n_samples), size=subsample)
            y_true = _safe_indexing(y_true, indices, axis=0)
            y_pred = _safe_indexing(y_pred, indices, axis=0)

        viz = PredictionErrorDisplay(
            y_true=y_true,
            y_pred=y_pred,
            scores=scores,
            with_residuals=with_residuals,
        )

        return viz.plot(
            ax=ax,
            scatter_kwargs=scatter_kwargs,
            line_kwargs=line_kwargs,
            residuals_kwargs=residuals_kwargs,
        )
