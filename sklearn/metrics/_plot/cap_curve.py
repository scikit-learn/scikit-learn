# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np

from sklearn.base import is_classifier
from sklearn.utils.multiclass import type_of_target

from ...utils._plotting import (
    _BinaryClassifierCurveDisplayMixin,
    _despine,
    _validate_style_kwargs,
)


class CAPCurveDisplay(_BinaryClassifierCurveDisplayMixin):
    """Cumulative Accuracy Profile (CAP) Curve visualization.

    It is recommended to use
    :func:`~sklearn.metrics.CAPCurveDisplay.from_estimator` or
    :func:`~sklearn.metrics.CAPCurveDisplay.from_predictions` to create
    a :class:`~sklearn.metrics.CAPCurveDisplay`. All parameters are
    stored as attributes.

    Read more in the :ref:`User Guide <visualizations>`.

    Parameters
    ----------
    y_true_cumulative : ndarray
        Cumulative number of true positives, ordered by predictions. If
        `sample_weight is not `None`, the positive cases are multiplied by
        their weight before being included into the cumulative sum.

    cumulative_total : ndarray
        Cumulative number or fraction of cases examined, ordered by
        predictions. If `sample_weight` is `None`, each case contributes
        equally. Otherwise this is the cumulative sum of weights,
        normalized if `normalize_scale=True`.

    estimator_name : str, default=None
        Name of estimator. If None, the estimator name is not shown.

    pos_label : int, float, bool or str, default=None
        The class considered as the positive class when computing the metrics.
        By default, `estimators.classes_[1]` is considered as the positive class.

    Attributes
    ----------
    line_ : matplotlib Artist
        CAP Curve.

    chance_level_ : matplotlib Artist
        Curve of the independent classifier.

    perfect_level_ : matplotlib Artist
        Curve of the perfect classifier.

    ax_ : matplotlib Axes
        Axes with CAP Curve.

    figure_ : matplotlib Figure
        Figure containing the curve.

    See Also
    --------
    CAPCurveDisplay.from_estimator : Plot Cumulative Accuracy Profile (CAP)
        CAP curve given the estimator and data.
    CAPCurveDisplay.from_predictions : Plot Cumulative Accuracy Profile (CAP)
        CAP curve given the observations and predicted values.
    RocCurveDisplay : Receiver Operating Characteristic (ROC) curve.
    PrecisionRecallDisplay : Precision-Recall curve.
    DetCurveDisplay : Detection Error Tradeoff (DET) curve.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> from sklearn import metrics
    >>> y_observed = np.array([0, 0, 1, 1])
    >>> y_pred = np.array([0.1, 0.4, 0.35, 0.8])
    >>> display = metrics.CAPCurveDisplay.from_predictions(y_observed, y_pred)
    >>> display.plot()
    <...>
    >>> plt.show()
    """

    def __init__(
        self,
        *,
        y_true_cumulative,
        cumulative_total,
        estimator_name=None,
        pos_label=None,
    ):
        self.estimator_name = estimator_name
        self.y_true_cumulative = y_true_cumulative
        self.cumulative_total = cumulative_total
        self.pos_label = pos_label

    def plot(
        self,
        ax=None,
        *,
        normalize_scale=True,
        name=None,
        plot_chance_level=True,
        chance_level_kwargs=None,
        plot_perfect=True,
        perfect_level_kwargs=None,
        despine=False,
        **kwargs,
    ):
        """Plot visualization.

        Extra keyword arguments will be passed to matplotlib's ``plot``.

        Parameters
        ----------
        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is created.

        normalize_scale : bool, default=True
            Whether to normalize values between 0 and 1 for the plot.

        name : str, default=None
            Name of CAP Curve for labeling. If `None`, use `self.estimator_name`.

        plot_chance_level : bool, default=True
            Whether to plot the expected curve of a classifier whose
            predictions are independent of the feature values passed to
            `decision_function` or `predict_proba`.

        chance_level_kwargs : dict, default=None
            Keyword arguments to be passed to matplotlib's `plot` for rendering
            the chance level line.

        plot_perfect : bool, default=True
            Whether to plot the perfect model line.

        perfect_level_kwargs : dict, default=None
            Keyword arguments to be passed to matplotlib's `plot` for rendering
            the perfect line.

        despine : bool, default=False
            Whether to remove the top and right spines from the plot.

        **kwargs : dict
            Keyword arguments to be passed to matplotlib's `plot`.

        Returns
        -------
        display : :class:`~sklearn.metrics.CAPCurveDisplay`
            Object that stores computed values.
        """

        self.ax_, self.figure_, name = self._validate_plot_params(ax=ax, name=name)

        x_max = self.cumulative_total[-1]
        y_max = self.y_true_cumulative[-1]
        if normalize_scale:
            if x_max != 0:
                self.cumulative_total = self.cumulative_total / x_max
            if y_max != 0:
                self.y_true_cumulative = self.y_true_cumulative / y_max
            self.ax_.set_xlim(-0.01, 1.01)
            self.ax_.set_ylim(-0.01, 1.01)
        else:
            self.ax_.set_xlim(0, x_max)
            self.ax_.set_ylim(0, y_max)

        line_kwargs = {"label": name} if name is not None else {}
        line_kwargs.update(**kwargs)

        default_perfect_level_line_kw = {
            "label": "Perfect predictions",
            "linestyle": ":",
            "color": "black",
        }

        if perfect_level_kwargs is None:
            perfect_level_kwargs = {}

        perfect_level_line_kw = _validate_style_kwargs(
            default_perfect_level_line_kw, perfect_level_kwargs
        )

        default_chance_level_line_kw = {
            "label": "Chance level",
            "color": "k",
            "linestyle": "--",
        }

        if chance_level_kwargs is None:
            chance_level_kwargs = {}

        chance_level_line_kw = _validate_style_kwargs(
            default_chance_level_line_kw, chance_level_kwargs
        )

        (self.line_,) = self.ax_.plot(
            self.cumulative_total, self.y_true_cumulative, **line_kwargs
        )

        if plot_chance_level:
            if normalize_scale:
                x_vals, y_vals = (0, 1), (0, 1)
            else:
                x_vals, y_vals = (0, x_max), (0, y_max)
            (self.chance_level_,) = self.ax_.plot(
                x_vals, y_vals, **chance_level_line_kw
            )
        else:
            self.chance_level_ = None

        if plot_perfect:
            pos_rate = y_max / x_max
            if normalize_scale:
                x_perfect = [0.0, pos_rate, 1.0]
                y_perfect = [0.0, 1.0, 1.0]
            else:
                x_perfect = [0, y_max, x_max]
                y_perfect = [0, y_max, y_max]

            (self.perfect_level_,) = self.ax_.plot(
                x_perfect, y_perfect, **perfect_level_line_kw
            )

            self._perfect_x = np.array(x_perfect)
            self._perfect_y = np.array(y_perfect)
        else:
            self.perfect_level_ = None
            self._perfect_x = None
            self._perfect_y = None

        info_pos_label = (
            f" (Positive label: {self.pos_label})" if self.pos_label is not None else ""
        )

        xlabel = "Total Cases Examined"
        ylabel = "Cumulative True Positives" + info_pos_label
        self.ax_.set(xlabel=xlabel, ylabel=ylabel, aspect="equal")

        if despine:
            _despine(self.ax_)

        if "label" in line_kwargs:
            self.ax_.legend(loc="lower right")

        return self

    @classmethod
    def from_predictions(
        cls,
        y_true,
        y_pred,
        *,
        sample_weight=None,
        pos_label=None,
        normalize_scale=True,
        plot_chance_level=True,
        chance_level_kwargs=None,
        plot_perfect=True,
        perfect_level_kwargs=None,
        despine=False,
        name=None,
        ax=None,
        **kwargs,
    ):
        """Create the Cumulative Accuracy Profile.

        This is also known as a (cumulative) gain curve. It can be used with
        either classifiers or regressors (with positive targets only). Other
        types of estimators will give rise to an error.

        Parameters
        ----------
        y_true : array-like of shape (n_samples,)
            True labels (for classifiers) or true values (for regressors).

        y_pred : array-like of shape (n_samples,)
            Target scores, can either be probability estimates of the positive
            class, confidence values, non-thresholded measure of decisions
            (as returned by :term:`decision_function` on some classifiers) or the
            :term:`predict` output in the case of a regressor.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights. If not `None`, the `y_true_cumulative` array is the
            cumulative sum of the weights of the positive cases (ordered by
            predictions) and `cumulative_total` is the cumulative sum of all
            cases (ordered by predictions).

        pos_label : int, float, bool or str, default=None
            The label of the positive class. When `pos_label=None`, if `y_true`
            is in {-1, 1} or {0, 1}, `pos_label` is set to 1, otherwise an
            error will be raised. Should not be used in the case of regression.

        normalize_scale : bool, default=True
            Whether to normalize values between 0 and 1 for the plot.

        plot_chance_level : bool, default=True
            Whether to plot the expected curve of a classifier whose
            predictions are independent of the feature values passed to
            `decision_function` or `predict_proba`.

        chance_level_kwargs : dict, default=None
            Keyword arguments to be passed to matplotlib's `plot` for rendering
            the chance level line.

        plot_perfect : bool, default=True
            Whether to plot the perfect model line.

        perfect_level_kwargs : dict, default=None
            Keyword arguments to be passed to matplotlib's `plot` for rendering
            a curve representing the perfect estimator.

        despine : bool, default=False
            Whether to remove the top and right spines from the plot.

        name : str, default=None
            Name of CAP curve for labeling in the legend. If `None`, name will be set to
            `"Classifier"` or `"Regressor"`.

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is created.

        **kwargs : dict
            Additional keywords arguments passed to matplotlib `plot` function.

        Returns
        -------
        display : :class:`~sklearn.metrics.CAPCurveDisplay`
            Object that stores computed values.
        """
        target_type = type_of_target(y_true)
        if target_type == "continuous" and pos_label is None:
            if np.any(y_true < 0):
                raise ValueError(
                    "`y_true` contains negative values, which isn't allowed for "
                    "continuous targets. If your data shouldn't be treated as "
                    "continuous, try converting the values to integers or strings "
                    "instead."
                )

            pos_label_validated = None

            if name is None:
                name = "Regressor"

            if sample_weight is None:
                sample_weight = np.ones_like(y_true, dtype=np.float64)

            sorted_indices = np.argsort(y_pred)
            y_true_sorted = y_true[sorted_indices]
            sample_weight_sorted = sample_weight[sorted_indices]

            weighted_y_true_sorted = y_true_sorted * sample_weight_sorted
            y_true_cumulative = np.cumsum(weighted_y_true_sorted)
            cumulative_total = np.cumsum(sample_weight_sorted)

        else:
            pos_label_validated, name = cls._validate_from_predictions_params(
                y_true,
                y_pred,
                sample_weight=sample_weight,
                pos_label=pos_label,
                name=name,
            )

            if sample_weight is None:
                sample_weight = np.ones_like(y_true, dtype=y_pred.dtype)

            y_true_bool = y_true == pos_label_validated

            sorted_indices = np.argsort(y_pred)[::-1]
            y_true_bool_sorted = y_true_bool[sorted_indices]
            sample_weight_sorted = sample_weight[sorted_indices]

            y_true_cumulative = np.cumsum(y_true_bool_sorted * sample_weight_sorted)
            cumulative_total = np.cumsum(sample_weight_sorted)

        viz = cls(
            y_true_cumulative=y_true_cumulative,
            cumulative_total=cumulative_total,
            estimator_name=name,
            pos_label=pos_label_validated,
        )

        return viz.plot(
            ax=ax,
            name=name,
            normalize_scale=normalize_scale,
            plot_chance_level=plot_chance_level,
            chance_level_kwargs=chance_level_kwargs,
            plot_perfect=plot_perfect,
            perfect_level_kwargs=perfect_level_kwargs,
            despine=despine,
            **kwargs,
        )

    @classmethod
    def from_estimator(
        cls,
        estimator,
        X,
        y,
        *,
        sample_weight=None,
        response_method="auto",
        pos_label=None,
        normalize_scale=True,
        plot_chance_level=True,
        chance_level_kwargs=None,
        plot_perfect=True,
        perfect_level_kwargs=None,
        despine=False,
        name=None,
        ax=None,
        **kwargs,
    ):
        """Create the Cumulative Accuracy Profile.

        This is also known as a (cumulative) gain curve. It can be used with
        either classifiers or regressors (with positive outcomes only). Other
        types of estimators will give rise to an error.

        Parameters
        ----------
        estimator : BaseEstimator
            A fitted estimator object (classifier or regressor with positive
            outcomes) implementing :term:`predict`, :term:`predict_proba`, or
            :term:`decision_function`. Multiclass classifiers are not supported.

        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Input values.

        y : array-like of shape (n_samples,)
            Target values.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights. If not `None`, the `y_true_cumulative` array is the
            cumulative sum of the weights of the positive cases (ordered by
            predictions) and `cumulative_total` is the cumulative sum of all
            cases (ordered by predictions).

        response_method : {'predict_proba', 'decision_function', 'auto'} \
                default='auto'
            Specifies whether to use :term:`predict_proba` or
            :term:`decision_function` as the target response. For regressors
            this parameter is ignored and the response is always the output of
            :term:`predict`. By default, :term:`predict_proba` is tried first
            and we revert to :term:`decision_function` if it doesn't exist.

        pos_label : int, float, bool or str, default=None
            The class considered as the positive class when computing CAP.
            By default, `estimator.classes_[1]` is considered as the positive class.

        normalize_scale : bool, default=True
            Whether to normalize values between 0 and 1 for the plot.

        plot_chance_level : bool, default=True
            Whether to plot the expected curve of a classifier whose
            predictions are independent of the feature values passed to
            `decision_function` or `predict_proba`.

        chance_level_kwargs : dict, default=None
            Keyword arguments to be passed to matplotlib's `plot` for rendering
            the chance level line.

        plot_perfect : bool, default=True
            Whether to plot the perfect model line.

        perfect_level_kwargs : dict, default=None
            Keyword arguments to be passed to matplotlib's `plot` for rendering
            the perfect line.

        despine : bool, default=False
            Whether to remove the top and right spines from the plot.

        name : str, default=None
            Name of CAP Curve for labeling. If `None`, use the name of the estimator.

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is created.

        **kwargs : dict
            Keyword arguments to be passed to matplotlib's `plot`.

        Returns
        -------
        display : :class:`~sklearn.metrics.CAPCurveDisplay`
            The CAP Curve display.
        """
        if is_classifier(estimator):
            y_pred, pos_label, name = cls._validate_and_get_response_values(
                estimator,
                X,
                y,
                response_method=response_method,
                pos_label=pos_label,
                name=name,
            )
        else:
            y_pred = estimator.predict(X)
            if name is None:
                name = type(estimator).__name__

        return cls.from_predictions(
            y_true=y,
            y_pred=y_pred,
            sample_weight=sample_weight,
            name=name,
            normalize_scale=normalize_scale,
            plot_chance_level=plot_chance_level,
            chance_level_kwargs=chance_level_kwargs,
            plot_perfect=plot_perfect,
            perfect_level_kwargs=perfect_level_kwargs,
            despine=despine,
            ax=ax,
            pos_label=pos_label,
            **kwargs,
        )
