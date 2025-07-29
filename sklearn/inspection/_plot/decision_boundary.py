# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import warnings

import numpy as np

from ...base import is_regressor
from ...preprocessing import LabelEncoder
from ...utils import _safe_indexing
from ...utils._optional_dependencies import check_matplotlib_support
from ...utils._response import _get_response_values
from ...utils._set_output import _get_adapter_from_container
from ...utils.validation import (
    _is_arraylike_not_scalar,
    _is_pandas_df,
    _is_polars_df,
    _num_features,
    check_is_fitted,
)


def _check_boundary_response_method(estimator, response_method, class_of_interest):
    """Validate the response methods to be used with the fitted estimator.

    Parameters
    ----------
    estimator : object
        Fitted estimator to check.

    response_method : {'auto', 'decision_function', 'predict_proba', 'predict'}
        Specifies whether to use :term:`decision_function`, :term:`predict_proba`,
        :term:`predict` as the target response. If set to 'auto', the response method is
        tried in the before mentioned order.

    class_of_interest : int, float, bool, str or None
        The class considered when plotting the decision. Cannot be None if
        multiclass and `response_method` is 'predict_proba' or 'decision_function'.

        .. versionadded:: 1.4

    Returns
    -------
    prediction_method : list of str or str
        The name or list of names of the response methods to use.
    """
    has_classes = hasattr(estimator, "classes_")
    if has_classes and _is_arraylike_not_scalar(estimator.classes_[0]):
        msg = "Multi-label and multi-output multi-class classifiers are not supported"
        raise ValueError(msg)

    if response_method == "auto":
        if is_regressor(estimator):
            prediction_method = "predict"
        else:
            prediction_method = ["decision_function", "predict_proba", "predict"]
    else:
        prediction_method = response_method

    return prediction_method


class DecisionBoundaryDisplay:
    """Decisions boundary visualization.

    It is recommended to use
    :func:`~sklearn.inspection.DecisionBoundaryDisplay.from_estimator`
    to create a :class:`DecisionBoundaryDisplay`. All parameters are stored as
    attributes.

    Read more in the :ref:`User Guide <visualizations>`.

    For a detailed example comparing the decision boundaries of multinomial and
    one-vs-rest logistic regression, please see
    :ref:`sphx_glr_auto_examples_linear_model_plot_logistic_multinomial.py`.

    .. versionadded:: 1.1

    Parameters
    ----------
    xx0 : ndarray of shape (grid_resolution, grid_resolution)
        First output of :func:`meshgrid <numpy.meshgrid>`.

    xx1 : ndarray of shape (grid_resolution, grid_resolution)
        Second output of :func:`meshgrid <numpy.meshgrid>`.

    response : ndarray of shape (grid_resolution, grid_resolution) or \
            (grid_resolution, grid_resolution, n_classes)
        Values of the response function.

    multiclass_colors : list of str or str, default=None
        Specifies how to color each class when plotting all classes of multiclass
        problem. Ignored for binary problems and multiclass problems when plotting a
        single prediction value per point.
        Possible inputs are:

        * list: list of Matplotlib
          `color <https://matplotlib.org/stable/users/explain/colors/colors.html#colors-def>`_
          strings, of length `n_classes`
        * str: name of :class:`matplotlib.colors.Colormap`
        * None: 'viridis' colormap is used to sample colors

        Single color colormaps will be generated from the colors in the list or
        colors taken from the colormap and passed to the `cmap` parameter of
        the `plot_method`.

        .. versionadded:: 1.7

    xlabel : str, default=None
        Default label to place on x axis.

    ylabel : str, default=None
        Default label to place on y axis.

    Attributes
    ----------
    surface_ : matplotlib `QuadContourSet` or `QuadMesh` or list of such objects
        If `plot_method` is 'contour' or 'contourf', `surface_` is
        :class:`QuadContourSet <matplotlib.contour.QuadContourSet>`. If
        `plot_method` is 'pcolormesh', `surface_` is
        :class:`QuadMesh <matplotlib.collections.QuadMesh>`.

    multiclass_colors_ : array of shape (n_classes, 4)
        Colors used to plot each class in multiclass problems.
        Only defined when `color_of_interest` is None.

        .. versionadded:: 1.7

    ax_ : matplotlib Axes
        Axes with decision boundary.

    figure_ : matplotlib Figure
        Figure containing the decision boundary.

    See Also
    --------
    DecisionBoundaryDisplay.from_estimator : Plot decision boundary given an estimator.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> from sklearn.datasets import load_iris
    >>> from sklearn.inspection import DecisionBoundaryDisplay
    >>> from sklearn.tree import DecisionTreeClassifier
    >>> iris = load_iris()
    >>> feature_1, feature_2 = np.meshgrid(
    ...     np.linspace(iris.data[:, 0].min(), iris.data[:, 0].max()),
    ...     np.linspace(iris.data[:, 1].min(), iris.data[:, 1].max())
    ... )
    >>> grid = np.vstack([feature_1.ravel(), feature_2.ravel()]).T
    >>> tree = DecisionTreeClassifier().fit(iris.data[:, :2], iris.target)
    >>> y_pred = np.reshape(tree.predict(grid), feature_1.shape)
    >>> display = DecisionBoundaryDisplay(
    ...     xx0=feature_1, xx1=feature_2, response=y_pred
    ... )
    >>> display.plot()
    <...>
    >>> display.ax_.scatter(
    ...     iris.data[:, 0], iris.data[:, 1], c=iris.target, edgecolor="black"
    ... )
    <...>
    >>> plt.show()
    """

    def __init__(
        self, *, xx0, xx1, response, multiclass_colors=None, xlabel=None, ylabel=None
    ):
        self.xx0 = xx0
        self.xx1 = xx1
        self.response = response
        self.multiclass_colors = multiclass_colors
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
            :func:`pcolormesh <matplotlib.pyplot.pcolormesh>`.

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
            Object that stores computed values.
        """
        check_matplotlib_support("DecisionBoundaryDisplay.plot")
        import matplotlib as mpl
        import matplotlib.pyplot as plt

        if plot_method not in ("contourf", "contour", "pcolormesh"):
            raise ValueError(
                "plot_method must be 'contourf', 'contour', or 'pcolormesh'. "
                f"Got {plot_method} instead."
            )

        if ax is None:
            _, ax = plt.subplots()

        plot_func = getattr(ax, plot_method)
        if self.response.ndim == 2:
            self.surface_ = plot_func(self.xx0, self.xx1, self.response, **kwargs)
        else:  # self.response.ndim == 3
            n_responses = self.response.shape[-1]
            for kwarg in ("cmap", "colors"):
                if kwarg in kwargs:
                    warnings.warn(
                        f"'{kwarg}' is ignored in favor of 'multiclass_colors' "
                        "in the multiclass case when the response method is "
                        "'decision_function' or 'predict_proba'."
                    )
                    del kwargs[kwarg]

            if self.multiclass_colors is None or isinstance(
                self.multiclass_colors, str
            ):
                if self.multiclass_colors is None:
                    cmap = "tab10" if n_responses <= 10 else "gist_rainbow"
                else:
                    cmap = self.multiclass_colors

                # Special case for the tab10 and tab20 colormaps that encode a
                # discrete set of colors that are easily distinguishable
                # contrary to other colormaps that are continuous.
                if cmap == "tab10" and n_responses <= 10:
                    colors = plt.get_cmap("tab10", 10).colors[:n_responses]
                elif cmap == "tab20" and n_responses <= 20:
                    colors = plt.get_cmap("tab20", 20).colors[:n_responses]
                else:
                    cmap = plt.get_cmap(cmap, n_responses)
                    if not hasattr(cmap, "colors"):
                        # For LinearSegmentedColormap
                        colors = cmap(np.linspace(0, 1, n_responses))
                    else:
                        colors = cmap.colors
            elif isinstance(self.multiclass_colors, list):
                colors = [mpl.colors.to_rgba(color) for color in self.multiclass_colors]
            else:
                raise ValueError("'multiclass_colors' must be a list or a str.")

            self.multiclass_colors_ = colors
            if plot_method == "contour":
                # Plot only argmax map for contour
                class_map = self.response.argmax(axis=2)
                self.surface_ = plot_func(
                    self.xx0, self.xx1, class_map, colors=colors, **kwargs
                )
            else:
                multiclass_cmaps = [
                    mpl.colors.LinearSegmentedColormap.from_list(
                        f"colormap_{class_idx}", [(1.0, 1.0, 1.0, 1.0), (r, g, b, 1.0)]
                    )
                    for class_idx, (r, g, b, _) in enumerate(colors)
                ]

                self.surface_ = []
                for class_idx, cmap in enumerate(multiclass_cmaps):
                    response = np.ma.array(
                        self.response[:, :, class_idx],
                        mask=~(self.response.argmax(axis=2) == class_idx),
                    )
                    self.surface_.append(
                        plot_func(self.xx0, self.xx1, response, cmap=cmap, **kwargs)
                    )

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
        class_of_interest=None,
        multiclass_colors=None,
        xlabel=None,
        ylabel=None,
        ax=None,
        **kwargs,
    ):
        """Plot decision boundary given an estimator.

        Read more in the :ref:`User Guide <visualizations>`.

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
            :func:`pcolormesh <matplotlib.pyplot.pcolormesh>`.

        response_method : {'auto', 'decision_function', 'predict_proba', \
                'predict'}, default='auto'
            Specifies whether to use :term:`decision_function`,
            :term:`predict_proba` or :term:`predict` as the target response.
            If set to 'auto', the response method is tried in the order as
            listed above.

            .. versionchanged:: 1.6
                For multiclass problems, 'auto' no longer defaults to 'predict'.

        class_of_interest : int, float, bool or str, default=None
            The class to be plotted when `response_method` is 'predict_proba'
            or 'decision_function'. If None, `estimator.classes_[1]` is considered
            the positive class for binary classifiers. For multiclass
            classifiers, if None, all classes will be represented in the
            decision boundary plot; the class with the highest response value
            at each point is plotted. The color of each class can be set via
            `multiclass_colors`.

            .. versionadded:: 1.4

        multiclass_colors : list of str, or str, default=None
            Specifies how to color each class when plotting multiclass
            'predict_proba' or 'decision_function' and `class_of_interest` is
            None. Ignored in all other cases.

            Possible inputs are:

            * list: list of Matplotlib
              `color <https://matplotlib.org/stable/users/explain/colors/colors.html#colors-def>`_
              strings, of length `n_classes`
            * str: name of :class:`matplotlib.colors.Colormap`
            * None: 'tab10' colormap is used to sample colors if the number of
                classes is less than or equal to 10, otherwise 'gist_rainbow'
                colormap.

            Single color colormaps will be generated from the colors in the list or
            colors taken from the colormap, and passed to the `cmap` parameter of
            the `plot_method`.

            .. versionadded:: 1.7

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
        sklearn.metrics.ConfusionMatrixDisplay.from_estimator : Plot the
            confusion matrix given an estimator, the data, and the label.
        sklearn.metrics.ConfusionMatrixDisplay.from_predictions : Plot the
            confusion matrix given the true and predicted labels.

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
        check_is_fitted(estimator)
        import matplotlib as mpl

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
            available_methods = ", ".join(possible_plot_methods)
            raise ValueError(
                f"plot_method must be one of {available_methods}. "
                f"Got {plot_method} instead."
            )

        num_features = _num_features(X)
        if num_features != 2:
            raise ValueError(
                f"n_features must be equal to 2. Got {num_features} instead."
            )

        if (
            response_method in ("predict_proba", "decision_function", "auto")
            and multiclass_colors is not None
            and hasattr(estimator, "classes_")
            and (n_classes := len(estimator.classes_)) > 2
        ):
            if isinstance(multiclass_colors, list):
                if len(multiclass_colors) != n_classes:
                    raise ValueError(
                        "When 'multiclass_colors' is a list, it must be of the same "
                        f"length as 'estimator.classes_' ({n_classes}), got: "
                        f"{len(multiclass_colors)}."
                    )
                elif any(
                    not mpl.colors.is_color_like(col) for col in multiclass_colors
                ):
                    raise ValueError(
                        "When 'multiclass_colors' is a list, it can only contain valid"
                        f" Matplotlib color names. Got: {multiclass_colors}"
                    )
            if isinstance(multiclass_colors, str):
                if multiclass_colors not in mpl.pyplot.colormaps():
                    raise ValueError(
                        "When 'multiclass_colors' is a string, it must be a valid "
                        f"Matplotlib colormap. Got: {multiclass_colors}"
                    )

        x0, x1 = _safe_indexing(X, 0, axis=1), _safe_indexing(X, 1, axis=1)

        x0_min, x0_max = x0.min() - eps, x0.max() + eps
        x1_min, x1_max = x1.min() - eps, x1.max() + eps

        xx0, xx1 = np.meshgrid(
            np.linspace(x0_min, x0_max, grid_resolution),
            np.linspace(x1_min, x1_max, grid_resolution),
        )

        X_grid = np.c_[xx0.ravel(), xx1.ravel()]
        if _is_pandas_df(X) or _is_polars_df(X):
            adapter = _get_adapter_from_container(X)
            X_grid = adapter.create_container(
                X_grid,
                X_grid,
                columns=X.columns,
            )

        prediction_method = _check_boundary_response_method(
            estimator, response_method, class_of_interest
        )
        try:
            response, _, response_method_used = _get_response_values(
                estimator,
                X_grid,
                response_method=prediction_method,
                pos_label=class_of_interest,
                return_response_method_used=True,
            )
        except ValueError as exc:
            if "is not a valid label" in str(exc):
                # re-raise a more informative error message since `pos_label` is unknown
                # to our user when interacting with
                # `DecisionBoundaryDisplay.from_estimator`
                raise ValueError(
                    f"class_of_interest={class_of_interest} is not a valid label: It "
                    f"should be one of {estimator.classes_}"
                ) from exc
            raise

        # convert classes predictions into integers
        if response_method_used == "predict" and hasattr(estimator, "classes_"):
            encoder = LabelEncoder()
            encoder.classes_ = estimator.classes_
            response = encoder.transform(response)

        if response.ndim == 1:
            response = response.reshape(*xx0.shape)
        else:
            if is_regressor(estimator):
                raise ValueError("Multi-output regressors are not supported")

            if class_of_interest is not None:
                # For the multiclass case, `_get_response_values` returns the response
                # as-is. Thus, we have a column per class and we need to select the
                # column corresponding to the positive class.
                col_idx = np.flatnonzero(estimator.classes_ == class_of_interest)[0]
                response = response[:, col_idx].reshape(*xx0.shape)
            else:
                response = response.reshape(*xx0.shape, response.shape[-1])

        if xlabel is None:
            xlabel = X.columns[0] if hasattr(X, "columns") else ""

        if ylabel is None:
            ylabel = X.columns[1] if hasattr(X, "columns") else ""

        display = cls(
            xx0=xx0,
            xx1=xx1,
            response=response,
            multiclass_colors=multiclass_colors,
            xlabel=xlabel,
            ylabel=ylabel,
        )
        return display.plot(ax=ax, plot_method=plot_method, **kwargs)
