# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
import scipy as sp

from sklearn.metrics._ranking import det_curve
from sklearn.utils import _safe_indexing
from sklearn.utils._plotting import (
    _BinaryClassifierCurveDisplayMixin,
    _check_param_lengths,
    _convert_to_list_leaving_none,
    _deprecate_estimator_name,
)
from sklearn.utils._response import _get_response_values_binary


class DetCurveDisplay(_BinaryClassifierCurveDisplayMixin):
    """Detection Error Tradeoff (DET) curve visualization.

    It is recommended to use
    :func:`~sklearn.metrics.DetCurveDisplay.from_estimator`,
    :func:`~sklearn.metrics.DetCurveDisplay.from_predictions` or
    :func:`~sklearn.metrics.DetCurveDisplay.from_cv_results` to create a
    :class:`~sklearn.metrics.DetCurveDisplay`. All parameters are stored as
    attributes.

    For general information regarding `scikit-learn` visualization tools, see
    the :ref:`Visualization Guide <visualizations>`.
    For guidance on interpreting these plots, refer to the
    :ref:`Model Evaluation Guide <det_curve>`.

    Parameters
    ----------
    fpr : ndarray or list of ndarrays
        False positive rate. Each ndarray should contain values for a single curve.
        If plotting multiple curves, `fpr` should be a list of the same length as
        `fnr`.

        .. versionchanged:: 1.11
            Now accepts a list for plotting multiple curves.

    fnr : ndarray or list of ndarrays
        False negative rate. Each ndarray should contain values for a single curve.
        If plotting multiple curves, `fnr` should be a list of the same length as
        `fpr`.

        .. versionchanged:: 1.11
            Now accepts a list for plotting multiple curves.

    name : str or list of str, default=None
        Name for labeling legend entries. The number of legend entries is determined
        by the `curve_kwargs` passed to `plot`, and is not affected by `name`.

        If a string is provided, it will be used to either label the single legend
        entry or if there are multiple legend entries, label each individual curve
        with the same name.

        If a list is provided, it will be used to label each curve individually.
        Passing a list will raise an error if `curve_kwargs` is not a list to avoid
        labeling individual curves that have the same appearance.

        If `None`, no name is shown in the legend.

        .. versionchanged:: 1.10
            `estimator_name` was deprecated in favor of `name`.

        .. versionchanged:: 1.11
            Now accepts a list for plotting multiple curves.

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
    line_ : matplotlib Artist or list of matplotlib Artists
        DET Curve.

        .. versionchanged:: 1.11
            This attribute can now be a list of Artists when multiple curves
            are plotted.

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
    DetCurveDisplay.from_cv_results : Plot multi-fold DET curves given
        cross-validation results.

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

    def _validate_plot_params(self, *, ax, name):
        self.ax_, self.figure_, name = super()._validate_plot_params(ax=ax, name=name)

        fpr = _convert_to_list_leaving_none(self.fpr)
        fnr = _convert_to_list_leaving_none(self.fnr)
        name = _convert_to_list_leaving_none(name)

        optional = {}
        if isinstance(name, list) and len(name) != 1:
            optional.update({"'name' (or self.name)": name})
        _check_param_lengths(
            required={"self.fpr": fpr, "self.fnr": fnr},
            optional=optional,
            class_name="DetCurveDisplay",
        )
        return fpr, fnr, name

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
        DetCurveDisplay.from_cv_results : Plot multi-fold DET curves given
            cross-validation results.

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
        y_score,
        *,
        sample_weight=None,
        drop_intermediate=True,
        pos_label=None,
        name=None,
        ax=None,
        curve_kwargs=None,
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
        DetCurveDisplay.from_cv_results : Plot multi-fold DET curves given
            cross-validation results.

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

    @classmethod
    def from_cv_results(
        cls,
        cv_results,
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
    ):
        """Create a multi-fold DET curve display given cross-validation results.

        .. versionadded:: 1.11

        Parameters
        ----------
        cv_results : dict
            Dictionary as returned by :func:`~sklearn.model_selection.cross_validate`
            using `return_estimator=True` and `return_indices=True` (i.e., dictionary
            should contain the keys "estimator" and "indices").

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

        response_method : {'predict_proba', 'decision_function', 'auto'} \
                default='auto'
            Specifies whether to use :term:`predict_proba` or
            :term:`decision_function` as the target response. If set to 'auto',
            :term:`predict_proba` is tried first and if it does not exist
            :term:`decision_function` is tried next.

        pos_label : int, float, bool or str, default=None
            The class considered as the positive class when computing the DET
            curve. By default, `estimators.classes_[1]` is considered as the
            positive class.

        name : str or list of str, default=None
            Name for labeling legend entries. The number of legend entries
            is determined by `curve_kwargs`, and is not affected by `name`.

            If a string is provided, it will be used to either label the single legend
            entry or if there are multiple legend entries, label each individual curve
            with the same name.

            If a list is provided, it will be used to label each curve individually.
            Passing a list will raise an error if `curve_kwargs` is not a list to avoid
            labeling individual curves that have the same appearance.

            If `None`, no name is shown in the legend.

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        curve_kwargs : dict or list of dict, default=None
            Keywords arguments to be passed to matplotlib's `plot` function
            to draw individual DET curves. If a list is provided the
            parameters are applied to the DET curves of each CV fold
            sequentially and a legend entry is added for each curve.
            If a single dictionary is provided, the same parameters are applied
            to all DET curves and a single legend entry for all curves is added.

        Returns
        -------
        display : :class:`~sklearn.metrics.DetCurveDisplay`
            The multi-fold DET curve display.

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
        >>> from sklearn.metrics import DetCurveDisplay
        >>> from sklearn.model_selection import cross_validate
        >>> from sklearn.svm import SVC
        >>> X, y = make_classification(n_samples=1000, random_state=0)
        >>> clf = SVC(random_state=0)
        >>> cv_results = cross_validate(
        ...     clf, X, y, cv=3, return_estimator=True, return_indices=True)
        >>> n_splits = len(cv_results["estimator"])
        >>> colors = ["#3f90da", "#ffa90e", "#bd1f01"]
        >>> DetCurveDisplay.from_cv_results(
        ...     cv_results,
        ...     X,
        ...     y,
        ...     name=[f"SVC fold {idx}" for idx in range(n_splits)],
        ...     curve_kwargs=[{"color": color} for color in colors],
        ... )
        <...>
        >>> plt.show()
        """
        cls._validate_from_cv_results_params(
            cv_results,
            X,
            y,
            sample_weight=sample_weight,
        )

        fpr_folds, fnr_folds = [], []
        for estimator, test_indices in zip(
            cv_results["estimator"], cv_results["indices"]["test"]
        ):
            y_true = _safe_indexing(y, test_indices)
            y_score, pos_label_ = _get_response_values_binary(
                estimator,
                _safe_indexing(X, test_indices),
                response_method=response_method,
                pos_label=pos_label,
            )
            sample_weight_fold = (
                None
                if sample_weight is None
                else _safe_indexing(sample_weight, test_indices)
            )
            fpr, fnr, _ = det_curve(
                y_true,
                y_score,
                pos_label=pos_label_,
                sample_weight=sample_weight_fold,
                drop_intermediate=drop_intermediate,
            )

            fpr_folds.append(fpr)
            fnr_folds.append(fnr)

        viz = cls(
            fpr=fpr_folds,
            fnr=fnr_folds,
            name=name,
            pos_label=pos_label_,
        )
        return viz.plot(ax=ax, curve_kwargs=curve_kwargs)

    def plot(self, ax=None, *, name=None, curve_kwargs=None, **kwargs):
        """Plot visualization.

        Parameters
        ----------
        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        name : str or list of str, default=None
            Name for labeling legend entries. The number of legend entries is
            determined by the `curve_kwargs` passed to `plot`, and is not affected
            by `name`.

            If a string is provided, it will be used to either label the single legend
            entry or if there are multiple legend entries, label each individual curve
            with the same name.

            If a list is provided, it will be used to label each curve individually.
            Passing a list will raise an error if `curve_kwargs` is not a list to avoid
            labeling individual curves that have the same appearance.

            If `None`, set to `name` provided at `DetCurveDisplay` initialization. If
            still `None`, no name is shown in the legend.

            .. versionchanged:: 1.11
                Now accepts a list for plotting multiple curves.

        curve_kwargs : dict or list of dict, default=None
            Keywords arguments to be passed to matplotlib's `plot` function
            to draw individual DET curves. For single curve plotting, should be
            a dictionary. For multi-curve plotting, if a list is provided the
            parameters are applied to each DET curve sequentially and a legend
            entry is added for each curve. If a single dictionary is provided,
            the same parameters are applied to all DET curves and a single legend
            entry for all curves is added.

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
        fpr, fnr, name = self._validate_plot_params(ax=ax, name=name)
        n_curves = len(fpr)
        _, legend_metric = self._get_legend_metric(curve_kwargs, n_curves, None)
        curve_kwargs = self._validate_curve_kwargs(
            n_curves,
            name,
            legend_metric,
            "",
            curve_kwargs=curve_kwargs,
            default_multi_curve_kwargs={
                "alpha": 0.5,
                "linestyle": "--",
                "color": "blue",
            },
            removed_version="1.12",
            **kwargs,
        )

        # We have the following bounds:
        # sp.stats.norm.ppf(0.0) = -np.inf
        # sp.stats.norm.ppf(1.0) = np.inf
        # We therefore clip to eps and 1 - eps to not provide infinity to matplotlib.
        fpr = [
            np.clip(f, np.finfo(f.dtype).eps, 1 - np.finfo(f.dtype).eps) for f in fpr
        ]
        fnr = [
            np.clip(f, np.finfo(f.dtype).eps, 1 - np.finfo(f.dtype).eps) for f in fnr
        ]
        single_curve = not isinstance(self.fpr, list)
        self.fpr = fpr[0] if single_curve else fpr
        self.fnr = fnr[0] if single_curve else fnr

        self.line_ = []
        for fpr_curve, fnr_curve, line_kw in zip(fpr, fnr, curve_kwargs):
            self.line_.extend(
                self.ax_.plot(
                    sp.stats.norm.ppf(fpr_curve),
                    sp.stats.norm.ppf(fnr_curve),
                    **line_kw,
                )
            )
        # Return single artist if only one curve is plotted
        if len(self.line_) == 1:
            self.line_ = self.line_[0]
        info_pos_label = (
            f" (Positive label: {self.pos_label})" if self.pos_label is not None else ""
        )

        xlabel = "False Positive Rate" + info_pos_label
        ylabel = "False Negative Rate" + info_pos_label
        self.ax_.set(xlabel=xlabel, ylabel=ylabel)

        if curve_kwargs[0].get("label") is not None:
            self.ax_.legend(loc="upper right")

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
