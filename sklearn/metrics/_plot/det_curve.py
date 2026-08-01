# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
import scipy as sp

from sklearn.metrics._ranking import det_curve
from sklearn.utils._plotting import (
    _BinaryClassifierCurveDisplayMixin,
    _deprecate_estimator_name,
    _deprecate_y_pred_parameter,
)


class DetCurveDisplay(_BinaryClassifierCurveDisplayMixin):
    """Detection Error Tradeoff (DET) curve visualization.

    It is recommended to use :func:`~sklearn.metrics.DetCurveDisplay.from_estimator`
    or :func:`~sklearn.metrics.DetCurveDisplay.from_predictions` to create a
    visualizer. All parameters are stored as attributes.

    For general information regarding `scikit-learn` visualization tools, see
    the :ref:`Visualization Guide <visualizations>`.
    For guidance on interpreting these plots, refer to the
    :ref:`Model Evaluation Guide <det_curve>`.

    .. versionadded:: 0.24

    Parameters
    ----------
    fpr : ndarray
        False positive rate.

    fnr : ndarray
        False negative rate.

    name : str, default=None
        Name for labeling the legend entry. If None, the estimator name is not shown.

        .. versionchanged:: 1.10
            `estimator_name` was deprecated in favor of `name`.

    pos_label : int, float, bool or str, default=None
        The label of the positive class. If not `None`, this value is displayed in
        the x- and y-axes labels.

    estimator_name : str, default=None
        Name of estimator. If None, the estimator name is not shown.

        .. deprecated:: 1.10
            `estimator_name` is deprecated and will be removed in 1.12. Use `name`
            instead.

    Attributes
    ----------
    line_ : matplotlib Artist
        DET Curve.

    ax_ : matplotlib Axes
        Axes with DET Curve.

    figure_ : matplotlib Figure
        Figure containing the curve.

    See Also
    --------
    det_curve : Compute error rates for different probability thresholds.
    DetCurveDisplay.from_estimator : Plot DET curve given an estimator and
        some data.
    DetCurveDisplay.from_predictions : Plot DET curve given the true and
        predicted labels.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from sklearn.datasets import make_classification
    >>> from sklearn.metrics import det_curve, DetCurveDisplay
    >>> from sklearn.model_selection import train_test_split
    >>> from sklearn.svm import SVC
    >>> X, y = make_classification(n_samples=1000, random_state=0)
    >>> X_train, X_test, y_train, y_test = train_test_split(
    ...     X, y, test_size=0.4, random_state=0)
    >>> clf = SVC(random_state=0).fit(X_train, y_train)
    >>> y_score = clf.decision_function(X_test)
    >>> fpr, fnr, _ = det_curve(y_test, y_score)
    >>> display = DetCurveDisplay(
    ...     fpr=fpr, fnr=fnr, name="SVC"
    ... )
    >>> display.plot()
    <...>
    >>> plt.show()
    """

    def __init__(
        self, *, fpr, fnr, name=None, pos_label=None, estimator_name="deprecated"
    ):
        self.fpr = fpr
        self.fnr = fnr
        self.name = _deprecate_estimator_name(estimator_name, name, "1.10")
        self.pos_label = pos_label

    @classmethod
    def from_estimator(
        cls,
        estimator,
        X,
        y,
        *,
        sample_weight=None,
        drop_intermediate=True,
        response_method="auto",
        pos_label=None,
        name=None,
        ax=None,
        curve_kwargs=None,
        **kwargs,
    ):
        """Plot DET curve given an estimator and data.

        For general information regarding `scikit-learn` visualization tools, see
        the :ref:`Visualization Guide <visualizations>`.
        For guidance on interpreting these plots, refer to the
        :ref:`Model Evaluation Guide <det_curve>`.

        .. versionadded:: 1.0

        Parameters
        ----------
        estimator : estimator instance
            Fitted classifier or a fitted :class:`~sklearn.pipeline.Pipeline`
            in which the last estimator is a classifier.

        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Input values.

        y : array-like of shape (n_samples,)
            Target values.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights.

        drop_intermediate : bool, default=True
            Whether to drop thresholds where true positives (tp) do not change
            from the previous or subsequent threshold. All points with the same
            tp value have the same `fnr` and thus same y coordinate.

            .. versionadded:: 1.7

        response_method : {'predict_proba', 'decision_function', 'auto'} \
                default='auto'
            Specifies whether to use :term:`predict_proba` or
            :term:`decision_function` as the predicted target response. If set
            to 'auto', :term:`predict_proba` is tried first and if it does not
            exist :term:`decision_function` is tried next.

        pos_label : int, float, bool or str, default=None
            The label of the positive class. By default, `estimators.classes_[1]`
            is considered as the positive class.

        name : str, default=None
            Name of DET curve for labeling. If `None`, use the name of the
            estimator.

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        curve_kwargs : dict, default=None
            Keywords arguments to be passed to matplotlib's `plot` function.

            .. versionadded:: 1.10

        **kwargs : dict
            Additional keywords arguments passed to matplotlib `plot` function.

            .. deprecated:: 1.10
                `**kwargs` is deprecated and will be removed in 1.12. Pass matplotlib
                arguments to `curve_kwargs` as a dictionary instead.

        Returns
        -------
        display : :class:`~sklearn.metrics.DetCurveDisplay`
            Object that stores computed values.

        See Also
        --------
        det_curve : Compute error rates for different probability thresholds.
        DetCurveDisplay.from_predictions : Plot DET curve given the true and
            predicted labels.

        Examples
        --------
        >>> import matplotlib.pyplot as plt
        >>> from sklearn.datasets import make_classification
        >>> from sklearn.metrics import DetCurveDisplay
        >>> from sklearn.model_selection import train_test_split
        >>> from sklearn.svm import SVC
        >>> X, y = make_classification(n_samples=1000, random_state=0)
        >>> X_train, X_test, y_train, y_test = train_test_split(
        ...     X, y, test_size=0.4, random_state=0)
        >>> clf = SVC(random_state=0).fit(X_train, y_train)
        >>> DetCurveDisplay.from_estimator(
        ...    clf, X_test, y_test)
        <...>
        >>> plt.show()
        """
        y_score, pos_label, name = cls._validate_and_get_response_values(
            estimator,
            X,
            y,
            response_method=response_method,
            pos_label=pos_label,
            name=name,
        )

        return cls.from_predictions(
            y_true=y,
            y_score=y_score,
            sample_weight=sample_weight,
            drop_intermediate=drop_intermediate,
            name=name,
            ax=ax,
            curve_kwargs=curve_kwargs,
            pos_label=pos_label,
            **kwargs,
        )

    @classmethod
    def from_predictions(
        cls,
        y_true,
        y_score=None,
        *,
        sample_weight=None,
        drop_intermediate=True,
        pos_label=None,
        name=None,
        ax=None,
        curve_kwargs=None,
        y_pred="deprecated",
        **kwargs,
    ):
        """Plot the DET curve given the true and predicted labels.

        For general information regarding `scikit-learn` visualization tools, see
        the :ref:`Visualization Guide <visualizations>`.
        For guidance on interpreting these plots, refer to the
        :ref:`Model Evaluation Guide <det_curve>`.

        .. versionadded:: 1.0

        Parameters
        ----------
        y_true : array-like of shape (n_samples,)
            True labels.

        y_score : array-like of shape (n_samples,)
            Target scores, can either be probability estimates of the positive
            class or non-thresholded decision values (as returned by
            :term:`decision_function` on some classifiers).

            .. versionadded:: 1.8
                `y_pred` has been renamed to `y_score`.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights.

        drop_intermediate : bool, default=True
            Whether to drop thresholds where true positives (tp) do not change
            from the previous or subsequent threshold. All points with the same
            tp value have the same `fnr` and thus same y coordinate.

            .. versionadded:: 1.7

        pos_label : int, float, bool or str, default=None
            The label of the positive class. When `pos_label=None`, if `y_true`
            is in {-1, 1} or {0, 1}, `pos_label` is set to 1, otherwise an
            error will be raised.

        name : str, default=None
            Name of DET curve for labeling. If `None`, name will be set to
            `"Classifier"`.

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        curve_kwargs : dict, default=None
            Keywords arguments to be passed to matplotlib's `plot` function.

            .. versionadded:: 1.10

        y_pred : array-like of shape (n_samples,)
            Target scores, can either be probability estimates of the positive
            class or non-thresholded decision values (as returned by
            :term:`decision_function` on some classifiers).

            .. deprecated:: 1.8
                `y_pred` is deprecated and will be removed in 1.10. Use
                `y_score` instead.

        **kwargs : dict
            Additional keywords arguments passed to matplotlib `plot` function.

            .. deprecated:: 1.10
                `**kwargs` is deprecated and will be removed in 1.12. Pass matplotlib
                arguments to `curve_kwargs` as a dictionary instead.

        Returns
        -------
        display : :class:`~sklearn.metrics.DetCurveDisplay`
            Object that stores computed values.

        See Also
        --------
        det_curve : Compute error rates for different probability thresholds.
        DetCurveDisplay.from_estimator : Plot DET curve given an estimator and
            some data.

        Examples
        --------
        >>> import matplotlib.pyplot as plt
        >>> from sklearn.datasets import make_classification
        >>> from sklearn.metrics import DetCurveDisplay
        >>> from sklearn.model_selection import train_test_split
        >>> from sklearn.svm import SVC
        >>> X, y = make_classification(n_samples=1000, random_state=0)
        >>> X_train, X_test, y_train, y_test = train_test_split(
        ...     X, y, test_size=0.4, random_state=0)
        >>> clf = SVC(random_state=0).fit(X_train, y_train)
        >>> y_score = clf.decision_function(X_test)
        >>> DetCurveDisplay.from_predictions(
        ...    y_test, y_score)
        <...>
        >>> plt.show()
        """
        y_score = _deprecate_y_pred_parameter(y_score, y_pred, "1.8")
        pos_label_validated, name = cls._validate_from_predictions_params(
            y_true, y_score, sample_weight=sample_weight, pos_label=pos_label, name=name
        )

        fpr, fnr, _ = det_curve(
            y_true,
            y_score,
            pos_label=pos_label,
            sample_weight=sample_weight,
            drop_intermediate=drop_intermediate,
        )

        viz = cls(
            fpr=fpr,
            fnr=fnr,
            name=name,
            pos_label=pos_label_validated,
        )

        return viz.plot(ax=ax, name=name, curve_kwargs=curve_kwargs, **kwargs)

    def plot(self, ax=None, *, name=None, curve_kwargs=None, **kwargs):
        """Plot visualization.

        Parameters
        ----------
        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        name : str, default=None
            Name of DET curve for labeling. If `None`, use `name` provided at
            `DetCurveDisplay` initialization, otherwise no labeling is shown.

        curve_kwargs : dict, default=None
            Keywords arguments to be passed to matplotlib's `plot` function.

            .. versionadded:: 1.10

        **kwargs : dict
            Additional keywords arguments passed to matplotlib `plot` function.

            .. deprecated:: 1.10
                `**kwargs` is deprecated and will be removed in 1.12. Pass matplotlib
                arguments to `curve_kwargs` as a dictionary instead.

        Returns
        -------
        display : :class:`~sklearn.metrics.DetCurveDisplay`
            Object that stores computed values.
        """
        self.ax_, self.figure_, name = self._validate_plot_params(ax=ax, name=name)

        _, legend_metric = self._get_legend_metric(curve_kwargs, 1, None)
        (curve_kwargs,) = self._validate_curve_kwargs(
            1,
            name,
            legend_metric,
            "",
            curve_kwargs=curve_kwargs,
            removed_version="1.12",
            **kwargs,
        )

        # We have the following bounds:
        # sp.stats.norm.ppf(0.0) = -np.inf
        # sp.stats.norm.ppf(1.0) = np.inf
        # We therefore clip to eps and 1 - eps to not provide infinity to matplotlib.
        eps = np.finfo(self.fpr.dtype).eps
        self.fpr = self.fpr.clip(eps, 1 - eps)
        self.fnr = self.fnr.clip(eps, 1 - eps)

        (self.line_,) = self.ax_.plot(
            sp.stats.norm.ppf(self.fpr),
            sp.stats.norm.ppf(self.fnr),
            **curve_kwargs,
        )
        info_pos_label = (
            f" (Positive label: {self.pos_label})" if self.pos_label is not None else ""
        )

        xlabel = "False Positive Rate" + info_pos_label
        ylabel = "False Negative Rate" + info_pos_label
        self.ax_.set(xlabel=xlabel, ylabel=ylabel)

        if curve_kwargs.get("label") is not None:
            self.ax_.legend(loc="lower right")

        ticks = [0.001, 0.01, 0.05, 0.20, 0.5, 0.80, 0.95, 0.99, 0.999]
        tick_locations = sp.stats.norm.ppf(ticks)
        tick_labels = [
            "{:.0%}".format(s) if (100 * s).is_integer() else "{:.1%}".format(s)
            for s in ticks
        ]
        self.ax_.set_xticks(tick_locations)
        self.ax_.set_xticklabels(tick_labels)
        self.ax_.set_xlim(-3, 3)
        self.ax_.set_yticks(tick_locations)
        self.ax_.set_yticklabels(tick_labels)
        self.ax_.set_ylim(-3, 3)

        return self
