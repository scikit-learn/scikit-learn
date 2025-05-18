# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np

from ...base import is_regressor  # To check if estimator is a regressor
from ...metrics._regression_characteristic import rec_curve
from ...utils._optional_dependencies import check_matplotlib_support
from ...utils.validation import check_is_fitted


class RecCurveDisplay:
    """Regression Error Characteristic (REC) Curve visualization.

    It is recommended to use :func:`~sklearn.metrics.RecCurveDisplay.from_estimator`
    or :func:`~sklearn.metrics.RecCurveDisplay.from_predictions` to create
    a visualizer. All parameters are stored as attributes.


    Parameters
    ----------
    deviations : ndarray
        Sorted unique error tolerance values (x-coordinates).
    accuracy : ndarray
        Corresponding accuracy values (y-coordinates).
    estimator_name : str, default=None
        Name of the estimator. If `None`, then the name will be `"Model"`.
    loss : {'absolute', 'squared'}, default='absolute'
        The loss function used to compute the REC curve.
    constant_predictor_deviations : ndarray, default=None
        Deviations for the constant predictor's REC curve. `None` if not plotted
        or not computed.
    constant_predictor_accuracy : ndarray, default=None
        Accuracy for the constant predictor's REC curve. `None` if not plotted
        or not computed.
    constant_predictor_name : str, default=None
        Name of the constant predictor. `None` if not plotted or not computed.

    Attributes
    ----------
    line_ : matplotlib Artist
        REC curve.
    ax_ : matplotlib Axes
        Axes with REC curve.
    figure_ : matplotlib Figure
        Figure containing the curve.
    constant_predictor_line_ : matplotlib Artist, default=None
        Constant predictor REC curve. Only defined if a constant predictor
        was plotted.

    See Also
    --------
    rec_curve : Compute Regression Error Characteristic (REC) curve.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> from sklearn.linear_model import LinearRegression
    >>> X = np.array([[1], [2], [3], [4], [5]])
    >>> y = np.array([1, 2.5, 3, 4.5, 5])
    >>> estimator = LinearRegression().fit(X, y)
    >>> display = RecCurveDisplay.from_estimator(estimator, X, y, loss='absolute')
    <...>
    >>> y_pred = estimator.predict(X)
    >>> display_pred = RecCurveDisplay.from_predictions(
    ...     y, y_pred, loss='squared', name="My Model", plot_const_predictor=False
    ... )
    <...>
    """

    def __init__(
        self,
        *,
        deviations,
        accuracy,
        estimator_name=None,
        loss=None,
        constant_predictor_deviations=None,
        constant_predictor_accuracy=None,
        constant_predictor_name=None,
    ):
        self.deviations = deviations
        self.accuracy = accuracy
        self.estimator_name = estimator_name if estimator_name is not None else "Model"
        self.loss = loss

        self.constant_predictor_deviations = constant_predictor_deviations
        self.constant_predictor_accuracy = constant_predictor_accuracy
        self.constant_predictor_name = constant_predictor_name

    def plot(
        self,
        ax=None,
        *,
        name=None,
        plot_const_predictor=True,
        clip_max_const_error=True,
        **kwargs,
    ):
        """Plot visualization.

        Parameters
        ----------
        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is created.
        name : str, default=None
            Name of REC curve for labeling. If `None`, use the name stored in
            `estimator_name`.
        **kwargs : dict
            Keyword arguments to be passed to `matplotlib.pyplot.plot` for the
            main REC curve.

        Returns
        -------
        display : :class:`~sklearn.metrics.RecCurveDisplay`
            Object that stores computed values.
        """

        check_matplotlib_support(f"{self.__class__.__name__}.plot")
        import matplotlib.pyplot as plt

        if ax is None:
            self.figure_, self.ax_ = plt.subplots()
        else:
            self.ax_ = ax
            self.figure_ = self.ax_.figure

        plot_name = name if name is not None else self.estimator_name
        line_kwargs = {}
        if "label" not in kwargs:  # Allow user to override label
            line_kwargs["label"] = plot_name
        line_kwargs.update(kwargs)

        if self.constant_predictor_deviations is not None:
            max_const_error = max(self.constant_predictor_deviations)
        elif clip_max_const_error:
            raise ValueError(
                "clip_max_const_error is True, but no constant deviations were given."
            )

        if clip_max_const_error:
            mask = self.deviations <= max_const_error
            self.line_, *_ = self.ax_.plot(
                self.deviations[mask], self.accuracy[mask], **line_kwargs
            )
        else:
            self.line_, *_ = self.ax_.plot(
                self.deviations, self.accuracy, **line_kwargs
            )

        # Plot constant predictor if its data exists and the flag is set
        if (
            plot_const_predictor
            and self.constant_predictor_deviations is not None
            and self.constant_predictor_accuracy is not None
        ):
            cp_name = (
                self.constant_predictor_name
                if self.constant_predictor_name
                else "Constant Predictor"
            )
            # Default style for constant predictor, can be overridden if needed
            cp_kwargs = {"label": cp_name, "linestyle": "--"}
            self.constant_predictor_line_, *_ = self.ax_.plot(
                self.constant_predictor_deviations,
                self.constant_predictor_accuracy,
                **cp_kwargs,
            )

        self.ax_.set_xlabel(f"Error Tolerance (Deviation - {self.loss} loss)")
        self.ax_.set_ylabel("Accuracy (Fraction of samples)")
        self.ax_.set_title("Regression Error Characteristic (REC) Curve")
        self.ax_.legend(loc="lower right")
        self.ax_.grid(True)

        return self

    @classmethod
    def from_estimator(
        cls,
        estimator,
        X,
        y,
        *,
        loss="absolute",
        constant_predictor=None,
        plot_const_predictor=True,
        clip_max_const_error=True,
        name=None,
        ax=None,
        **kwargs,
    ):
        """Create a REC Curve display from an estimator.

        Parameters
        ----------
        estimator : object
            Fitted estimator or a fitted :class:`~sklearn.pipeline.Pipeline`
            in which the last estimator is a regressor.
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Input values.
        y : array-like of shape (n_samples,)
            Target values.
        loss : {'absolute', 'squared'}, default='absolute'
            The loss function to use for calculating deviations.
        constant_predictor : {'mean', 'median', None}, default=None
            The type of constant predictor to plot as a baseline.
            If 'mean', uses the mean of `y_true`.
            If 'median', uses the median of `y_true`.
            If `None`, chooses 'mean' for 'squared' loss and 'median' for
            'absolute' loss, as these are the optimal constant predictors
            for these losses.
        plot_const_predictor : bool, default=True
            Whether to compute and plot the REC curve for the constant predictor.
        clip_max_const_error : bool, default=True
            If `True`, the x-axis (error tolerance) will be cut off at the
            maximum error achieved by the constant predictor.
        name : str, default=None
            Name for the REC curve. If `None`, the estimator's class name will be used.
        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is created.
        **kwargs : dict
            Keyword arguments to be passed to `matplotlib.pyplot.plot` for the
            estimator's REC curve.

        Returns
        -------
        display : :class:`~sklearn.metrics.RecCurveDisplay`
            Object that stores computed values.
        """
        check_is_fitted(estimator)
        if not is_regressor(estimator):
            raise TypeError(f"{estimator.__class__.__name__} is not a regressor.")

        y_pred = estimator.predict(X)

        if name is None:
            name = estimator.__class__.__name__

        # Call from_predictions with validated parameters
        return cls.from_predictions(
            y_true=y,
            y_pred=y_pred,
            loss=loss,
            constant_predictor=constant_predictor,
            plot_const_predictor=plot_const_predictor,
            clip_max_const_error=clip_max_const_error,
            name=name,
            ax=ax,
            **kwargs,
        )

    @classmethod
    def from_predictions(
        cls,
        y_true,
        y_pred,
        *,
        loss="absolute",
        constant_predictor=None,
        plot_const_predictor=True,
        clip_max_const_error=True,
        name=None,
        ax=None,
        **kwargs,
    ):
        """Plot REC curve given true and predicted values.

        Parameters
        ----------
        y_true : array-like of shape (n_samples,)
            True target values.
        y_pred : array-like of shape (n_samples,) or scalar
            Estimated target values. The `rec_curve` function handles scalar `y_pred`.
        loss : {'absolute', 'squared'}, default='absolute'
            The loss function to use for calculating deviations.
        constant_predictor : {'mean', 'median', None}, default=None
            The type of constant predictor to plot as a baseline.
            If 'mean', uses the mean of `y_true`.
            If 'median', uses the median of `y_true`.
            If `None`, chooses 'mean' for 'squared' loss and 'median' for
            'absolute' loss.
        plot_const_predictor : bool, default=True
            Whether to compute and plot the REC curve for the constant predictor.
        clip_max_const_error : bool, default=True
            If `True`, the x-axis (error tolerance) will be cut off at the
            maximum error achieved by the constant predictor. This is only
            effective if a constant predictor is computed.
        name : str, default=None
            Name for the REC curve. If `None`, will be "Model".
        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is created.
        **kwargs : dict
            Keyword arguments to be passed to `matplotlib.pyplot.plot` for the
            main REC curve.

        Returns
        -------
        display : :class:`~sklearn.metrics.RecCurveDisplay`
            Object that stores computed values.
        """
        main_devs, main_acc = rec_curve(y_true, y_pred, loss=loss)

        # Convert to NumPy arrays for consistent handling in display class,
        # as rec_curve might accept and return array_api specific arrays.
        y_true_np = np.asarray(y_true)
        main_devs_np = np.asarray(main_devs)
        main_acc_np = np.asarray(main_acc)

        # Determine the constant predictor type if not specified
        actual_constant_predictor_type = constant_predictor
        if actual_constant_predictor_type is None:
            if loss == "squared":
                actual_constant_predictor_type = "mean"
            else:  # loss == "absolute"
                actual_constant_predictor_type = "median"

        # Compute constant predictor data if needed for plotting or cutoff
        if plot_const_predictor or clip_max_const_error:
            if actual_constant_predictor_type == "mean":
                constant_value = np.mean(y_true_np)
                cp_name_val = "Mean Predictor"
            elif actual_constant_predictor_type == "median":
                constant_value = np.median(y_true_np)
                cp_name_val = "Median Predictor"

            cp_devs, cp_accs = rec_curve(y_true, constant_value, loss=loss)
            cp_devs_np = np.asarray(cp_devs)
            cp_accs_np = np.asarray(cp_accs)
        else:
            cp_devs_np, cp_accs_np, cp_name_val = None, None, None

        display_name = name if name is not None else "Model"

        obj = RecCurveDisplay(
            deviations=main_devs_np,
            accuracy=main_acc_np,
            estimator_name=display_name,
            loss=loss,  # loss is already validated
            constant_predictor_deviations=cp_devs_np,
            constant_predictor_accuracy=cp_accs_np,
            constant_predictor_name=cp_name_val,
        )
        return obj.plot(
            ax=ax,
            plot_const_predictor=plot_const_predictor,
            clip_max_const_error=clip_max_const_error,
            **kwargs,
        )
