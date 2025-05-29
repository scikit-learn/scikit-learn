# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause


import warnings

import numpy as np

from ...utils import _safe_indexing
from ...utils._plotting import (
    _BinaryClassifierCurveDisplayMixin,
    _check_param_lengths,
    _convert_to_list_leaving_none,
    _deprecate_estimator_name,
    _despine,
    _validate_style_kwargs,
)
from ...utils._response import _get_response_values_binary
from .._ranking import auc, roc_curve


class RocCurveDisplay(_BinaryClassifierCurveDisplayMixin):
    """ROC Curve visualization.

    It is recommended to use
    :func:`~sklearn.metrics.RocCurveDisplay.from_estimator` or
    :func:`~sklearn.metrics.RocCurveDisplay.from_predictions` or
    :func:`~sklearn.metrics.RocCurveDisplay.from_cv_results` to create
    a :class:`~sklearn.metrics.RocCurveDisplay`. All parameters are
    stored as attributes.

    For general information regarding `scikit-learn` visualization tools, see
    the :ref:`Visualization Guide <visualizations>`.
    For guidance on interpreting these plots, refer to the :ref:`Model
    Evaluation Guide <roc_metrics>`.

    Parameters
    ----------
    fpr : ndarray or list of ndarrays
        False positive rates. Each ndarray should contain values for a single curve.
        If plotting multiple curves, list should be of same length as `tpr`.

        .. versionchanged:: 1.7
            Now accepts a list for plotting multiple curves.

    tpr : ndarray or list of ndarrays
        True positive rates. Each ndarray should contain values for a single curve.
        If plotting multiple curves, list should be of same length as `fpr`.

        .. versionchanged:: 1.7
            Now accepts a list for plotting multiple curves.

    roc_auc : float or list of floats, default=None
        Area under ROC curve, used for labeling each curve in the legend.
        If plotting multiple curves, should be a list of the same length as `fpr`
        and `tpr`. If `None`, ROC AUC scores are not shown in the legend.

        .. versionchanged:: 1.7
            Now accepts a list for plotting multiple curves.

    name : str or list of str, default=None
        Name for labeling legend entries. The number of legend entries
        is determined by the `curve_kwargs` passed to `plot`.
        To label each curve, provide a list of strings. To avoid labeling
        individual curves that have the same appearance, this cannot be used in
        conjunction with `curve_kwargs` being a dictionary or None. If a
        string is provided, it will be used to either label the single legend entry
        or if there are multiple legend entries, label each individual curve with
        the same name. If `None`, set to `name` provided at `RocCurveDisplay`
        initialization. If still `None`, no name is shown in the legend.

        .. versionadded:: 1.7

    pos_label : int, float, bool or str, default=None
        The class considered as the positive class when computing the roc auc
        metrics. By default, `estimators.classes_[1]` is considered
        as the positive class.

        .. versionadded:: 0.24

    estimator_name : str, default=None
        Name of estimator. If None, the estimator name is not shown.

        .. deprecated:: 1.7
            `estimator_name` is deprecated and will be removed in 1.9. Use `name`
            instead.

    Attributes
    ----------
    line_ : matplotlib Artist or list of matplotlib Artists
        ROC Curves.

        .. versionchanged:: 1.7
            This attribute can now be a list of Artists, for when multiple curves
            are plotted.

    chance_level_ : matplotlib Artist or None
        The chance level line. It is `None` if the chance level is not plotted.

        .. versionadded:: 1.3

    ax_ : matplotlib Axes
        Axes with ROC Curve.

    figure_ : matplotlib Figure
        Figure containing the curve.

    See Also
    --------
    roc_curve : Compute Receiver operating characteristic (ROC) curve.
    RocCurveDisplay.from_estimator : Plot Receiver Operating Characteristic
        (ROC) curve given an estimator and some data.
    RocCurveDisplay.from_predictions : Plot Receiver Operating Characteristic
        (ROC) curve given the true and predicted values.
    roc_auc_score : Compute the area under the ROC curve.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> from sklearn import metrics
    >>> y_true = np.array([0, 0, 1, 1])
    >>> y_score = np.array([0.1, 0.4, 0.35, 0.8])
    >>> fpr, tpr, thresholds = metrics.roc_curve(y_true, y_score)
    >>> roc_auc = metrics.auc(fpr, tpr)
    >>> display = metrics.RocCurveDisplay(fpr=fpr, tpr=tpr, roc_auc=roc_auc,
    ...                                   name='example estimator')
    >>> display.plot()
    <...>
    >>> plt.show()
    """

    def __init__(
        self,
        *,
        fpr,
        tpr,
        roc_auc=None,
        name=None,
        pos_label=None,
        estimator_name="deprecated",
    ):
        self.fpr = fpr
        self.tpr = tpr
        self.roc_auc = roc_auc
        self.name = _deprecate_estimator_name(estimator_name, name, "1.7")
        self.pos_label = pos_label

    def _validate_plot_params(self, *, ax, name):
        self.ax_, self.figure_, name = super()._validate_plot_params(ax=ax, name=name)

        fpr = _convert_to_list_leaving_none(self.fpr)
        tpr = _convert_to_list_leaving_none(self.tpr)
        roc_auc = _convert_to_list_leaving_none(self.roc_auc)
        name = _convert_to_list_leaving_none(name)

        optional = {"self.roc_auc": roc_auc}
        if isinstance(name, list) and len(name) != 1:
            optional.update({"'name' (or self.name)": name})
        _check_param_lengths(
            required={"self.fpr": fpr, "self.tpr": tpr},
            optional=optional,
            class_name="RocCurveDisplay",
        )
        return fpr, tpr, roc_auc, name

    def plot(
        self,
        ax=None,
        *,
        name=None,
        curve_kwargs=None,
        plot_chance_level=False,
        chance_level_kw=None,
        despine=False,
        **kwargs,
    ):
        """Plot visualization.

        Parameters
        ----------
        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        name : str or list of str, default=None
            Name for labeling legend entries. The number of legend entries
            is determined by `curve_kwargs`.
            To label each curve, provide a list of strings. To avoid labeling
            individual curves that have the same appearance, this cannot be used in
            conjunction with `curve_kwargs` being a dictionary or None. If a
            string is provided, it will be used to either label the single legend entry
            or if there are multiple legend entries, label each individual curve with
            the same name. If `None`, set to `name` provided at `RocCurveDisplay`
            initialization. If still `None`, no name is shown in the legend.

            .. versionadded:: 1.7

        curve_kwargs : dict or list of dict, default=None
            Keywords arguments to be passed to matplotlib's `plot` function
            to draw individual ROC curves. For single curve plotting, should be
            a dictionary. For multi-curve plotting, if a list is provided the
            parameters are applied to the ROC curves of each CV fold
            sequentially and a legend entry is added for each curve.
            If a single dictionary is provided, the same parameters are applied
            to all ROC curves and a single legend entry for all curves is added,
            labeled with the mean ROC AUC score.

            .. versionadded:: 1.7

        plot_chance_level : bool, default=False
            Whether to plot the chance level.

            .. versionadded:: 1.3

        chance_level_kw : dict, default=None
            Keyword arguments to be passed to matplotlib's `plot` for rendering
            the chance level line.

            .. versionadded:: 1.3

        despine : bool, default=False
            Whether to remove the top and right spines from the plot.

            .. versionadded:: 1.6

        **kwargs : dict
            Keyword arguments to be passed to matplotlib's `plot`.

            .. deprecated:: 1.7
                kwargs is deprecated and will be removed in 1.9. Pass matplotlib
                arguments to `curve_kwargs` as a dictionary instead.

        Returns
        -------
        display : :class:`~sklearn.metrics.RocCurveDisplay`
            Object that stores computed values.
        """
        fpr, tpr, roc_auc, name = self._validate_plot_params(ax=ax, name=name)
        n_curves = len(fpr)
        if not isinstance(curve_kwargs, list) and n_curves > 1:
            if roc_auc:
                legend_metric = {"mean": np.mean(roc_auc), "std": np.std(roc_auc)}
            else:
                legend_metric = {"mean": None, "std": None}
        else:
            roc_auc = roc_auc if roc_auc is not None else [None] * n_curves
            legend_metric = {"metric": roc_auc}

        curve_kwargs = self._validate_curve_kwargs(
            n_curves,
            name,
            legend_metric,
            "AUC",
            curve_kwargs=curve_kwargs,
            **kwargs,
        )

        default_chance_level_line_kw = {
            "label": "Chance level (AUC = 0.5)",
            "color": "k",
            "linestyle": "--",
        }

        if chance_level_kw is None:
            chance_level_kw = {}

        chance_level_kw = _validate_style_kwargs(
            default_chance_level_line_kw, chance_level_kw
        )

        self.line_ = []
        for fpr, tpr, line_kw in zip(fpr, tpr, curve_kwargs):
            self.line_.extend(self.ax_.plot(fpr, tpr, **line_kw))
        # Return single artist if only one curve is plotted
        if len(self.line_) == 1:
            self.line_ = self.line_[0]

        info_pos_label = (
            f" (Positive label: {self.pos_label})" if self.pos_label is not None else ""
        )

        xlabel = "False Positive Rate" + info_pos_label
        ylabel = "True Positive Rate" + info_pos_label
        self.ax_.set(
            xlabel=xlabel,
            xlim=(-0.01, 1.01),
            ylabel=ylabel,
            ylim=(-0.01, 1.01),
            aspect="equal",
        )

        if plot_chance_level:
            (self.chance_level_,) = self.ax_.plot((0, 1), (0, 1), **chance_level_kw)
        else:
            self.chance_level_ = None

        if despine:
            _despine(self.ax_)

        if curve_kwargs[0].get("label") is not None or (
            plot_chance_level and chance_level_kw.get("label") is not None
        ):
            self.ax_.legend(loc="lower right")

        return self

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
        plot_chance_level=False,
        chance_level_kw=None,
        despine=False,
        **kwargs,
    ):
        """Create a ROC Curve display from an estimator.

        For general information regarding `scikit-learn` visualization tools,
        see the :ref:`Visualization Guide <visualizations>`.
        For guidance on interpreting these plots, refer to the :ref:`Model
        Evaluation Guide <roc_metrics>`.

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
            Whether to drop thresholds where the resulting point is collinear
            with its neighbors in ROC space. This has no effect on the ROC AUC
            or visual shape of the curve, but reduces the number of plotted
            points.

        response_method : {'predict_proba', 'decision_function', 'auto'} \
                default='auto'
            Specifies whether to use :term:`predict_proba` or
            :term:`decision_function` as the target response. If set to 'auto',
            :term:`predict_proba` is tried first and if it does not exist
            :term:`decision_function` is tried next.

        pos_label : int, float, bool or str, default=None
            The class considered as the positive class when computing the ROC AUC.
            By default, `estimators.classes_[1]` is considered
            as the positive class.

        name : str, default=None
            Name of ROC Curve for labeling. If `None`, use the name of the
            estimator.

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is created.

        curve_kwargs : dict, default=None
            Keywords arguments to be passed to matplotlib's `plot` function.

            .. versionadded:: 1.7

        plot_chance_level : bool, default=False
            Whether to plot the chance level.

            .. versionadded:: 1.3

        chance_level_kw : dict, default=None
            Keyword arguments to be passed to matplotlib's `plot` for rendering
            the chance level line.

            .. versionadded:: 1.3

        despine : bool, default=False
            Whether to remove the top and right spines from the plot.

            .. versionadded:: 1.6

        **kwargs : dict
            Keyword arguments to be passed to matplotlib's `plot`.

            .. deprecated:: 1.7
                kwargs is deprecated and will be removed in 1.9. Pass matplotlib
                arguments to `curve_kwargs` as a dictionary instead.

        Returns
        -------
        display : :class:`~sklearn.metrics.RocCurveDisplay`
            The ROC Curve display.

        See Also
        --------
        roc_curve : Compute Receiver operating characteristic (ROC) curve.
        RocCurveDisplay.from_predictions : ROC Curve visualization given the
            probabilities of scores of a classifier.
        roc_auc_score : Compute the area under the ROC curve.

        Examples
        --------
        >>> import matplotlib.pyplot as plt
        >>> from sklearn.datasets import make_classification
        >>> from sklearn.metrics import RocCurveDisplay
        >>> from sklearn.model_selection import train_test_split
        >>> from sklearn.svm import SVC
        >>> X, y = make_classification(random_state=0)
        >>> X_train, X_test, y_train, y_test = train_test_split(
        ...     X, y, random_state=0)
        >>> clf = SVC(random_state=0).fit(X_train, y_train)
        >>> RocCurveDisplay.from_estimator(
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
            pos_label=pos_label,
            curve_kwargs=curve_kwargs,
            plot_chance_level=plot_chance_level,
            chance_level_kw=chance_level_kw,
            despine=despine,
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
        plot_chance_level=False,
        chance_level_kw=None,
        despine=False,
        y_pred="deprecated",
        **kwargs,
    ):
        """Plot ROC curve given the true and predicted values.

        For general information regarding `scikit-learn` visualization tools,
        see the :ref:`Visualization Guide <visualizations>`.
        For guidance on interpreting these plots, refer to the :ref:`Model
        Evaluation Guide <roc_metrics>`.

        .. versionadded:: 1.0

        Parameters
        ----------
        y_true : array-like of shape (n_samples,)
            True labels.

        y_score : array-like of shape (n_samples,)
            Target scores, can either be probability estimates of the positive
            class, confidence values, or non-thresholded measure of decisions
            (as returned by “decision_function” on some classifiers).

            .. versionadded:: 1.7
                `y_pred` has been renamed to `y_score`.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights.

        drop_intermediate : bool, default=True
            Whether to drop thresholds where the resulting point is collinear
            with its neighbors in ROC space. This has no effect on the ROC AUC
            or visual shape of the curve, but reduces the number of plotted
            points.

        pos_label : int, float, bool or str, default=None
            The label of the positive class when computing the ROC AUC.
            When `pos_label=None`, if `y_true` is in {-1, 1} or {0, 1}, `pos_label`
            is set to 1, otherwise an error will be raised.

        name : str, default=None
            Name of ROC curve for legend labeling. If `None`, name will be set to
            `"Classifier"`.

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        curve_kwargs : dict, default=None
            Keywords arguments to be passed to matplotlib's `plot` function.

            .. versionadded:: 1.7

        plot_chance_level : bool, default=False
            Whether to plot the chance level.

            .. versionadded:: 1.3

        chance_level_kw : dict, default=None
            Keyword arguments to be passed to matplotlib's `plot` for rendering
            the chance level line.

            .. versionadded:: 1.3

        despine : bool, default=False
            Whether to remove the top and right spines from the plot.

            .. versionadded:: 1.6

        y_pred : array-like of shape (n_samples,)
            Target scores, can either be probability estimates of the positive
            class, confidence values, or non-thresholded measure of decisions
            (as returned by “decision_function” on some classifiers).

            .. deprecated:: 1.7
                `y_pred` is deprecated and will be removed in 1.9. Use
                `y_score` instead.

        **kwargs : dict
            Additional keywords arguments passed to matplotlib `plot` function.

            .. deprecated:: 1.7
                kwargs is deprecated and will be removed in 1.9. Pass matplotlib
                arguments to `curve_kwargs` as a dictionary instead.

        Returns
        -------
        display : :class:`~sklearn.metrics.RocCurveDisplay`
            Object that stores computed values.

        See Also
        --------
        roc_curve : Compute Receiver operating characteristic (ROC) curve.
        RocCurveDisplay.from_estimator : ROC Curve visualization given an
            estimator and some data.
        roc_auc_score : Compute the area under the ROC curve.

        Examples
        --------
        >>> import matplotlib.pyplot as plt
        >>> from sklearn.datasets import make_classification
        >>> from sklearn.metrics import RocCurveDisplay
        >>> from sklearn.model_selection import train_test_split
        >>> from sklearn.svm import SVC
        >>> X, y = make_classification(random_state=0)
        >>> X_train, X_test, y_train, y_test = train_test_split(
        ...     X, y, random_state=0)
        >>> clf = SVC(random_state=0).fit(X_train, y_train)
        >>> y_score = clf.decision_function(X_test)
        >>> RocCurveDisplay.from_predictions(y_test, y_score)
        <...>
        >>> plt.show()
        """
        # TODO(1.9): remove after the end of the deprecation period of `y_pred`
        if y_score is not None and not (
            isinstance(y_pred, str) and y_pred == "deprecated"
        ):
            raise ValueError(
                "`y_pred` and `y_score` cannot be both specified. Please use `y_score`"
                " only as `y_pred` is deprecated in 1.7 and will be removed in 1.9."
            )
        if not (isinstance(y_pred, str) and y_pred == "deprecated"):
            warnings.warn(
                (
                    "y_pred is deprecated in 1.7 and will be removed in 1.9. "
                    "Please use `y_score` instead."
                ),
                FutureWarning,
            )
            y_score = y_pred

        pos_label_validated, name = cls._validate_from_predictions_params(
            y_true, y_score, sample_weight=sample_weight, pos_label=pos_label, name=name
        )

        fpr, tpr, _ = roc_curve(
            y_true,
            y_score,
            pos_label=pos_label,
            sample_weight=sample_weight,
            drop_intermediate=drop_intermediate,
        )
        roc_auc = auc(fpr, tpr)

        viz = cls(
            fpr=fpr,
            tpr=tpr,
            roc_auc=roc_auc,
            name=name,
            pos_label=pos_label_validated,
        )

        return viz.plot(
            ax=ax,
            curve_kwargs=curve_kwargs,
            plot_chance_level=plot_chance_level,
            chance_level_kw=chance_level_kw,
            despine=despine,
            **kwargs,
        )

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
        ax=None,
        name=None,
        curve_kwargs=None,
        plot_chance_level=False,
        chance_level_kwargs=None,
        despine=False,
    ):
        """Create a multi-fold ROC curve display given cross-validation results.

        .. versionadded:: 1.7

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
            Whether to drop some suboptimal thresholds which would not appear
            on a plotted ROC curve. This is useful in order to create lighter
            ROC curves.

        response_method : {'predict_proba', 'decision_function', 'auto'} \
                default='auto'
            Specifies whether to use :term:`predict_proba` or
            :term:`decision_function` as the target response. If set to 'auto',
            :term:`predict_proba` is tried first and if it does not exist
            :term:`decision_function` is tried next.

        pos_label : int, float, bool or str, default=None
            The class considered as the positive class when computing the ROC AUC
            metrics. By default, `estimators.classes_[1]` is considered
            as the positive class.

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        name : str or list of str, default=None
            Name for labeling legend entries. The number of legend entries
            is determined by `curve_kwargs`.
            To label each curve, provide a list of strings. To avoid labeling
            individual curves that have the same appearance, this cannot be used in
            conjunction with `curve_kwargs` being a dictionary or None. If a
            string is provided, it will be used to either label the single legend entry
            or if there are multiple legend entries, label each individual curve with
            the same name. If `None`, no name is shown in the legend.

        curve_kwargs : dict or list of dict, default=None
            Keywords arguments to be passed to matplotlib's `plot` function
            to draw individual ROC curves. If a list is provided the
            parameters are applied to the ROC curves of each CV fold
            sequentially and a legend entry is added for each curve.
            If a single dictionary is provided, the same parameters are applied
            to all ROC curves and a single legend entry for all curves is added,
            labeled with the mean ROC AUC score.

        plot_chance_level : bool, default=False
            Whether to plot the chance level.

        chance_level_kwargs : dict, default=None
            Keyword arguments to be passed to matplotlib's `plot` for rendering
            the chance level line.

        despine : bool, default=False
            Whether to remove the top and right spines from the plot.

        Returns
        -------
        display : :class:`~sklearn.metrics.RocCurveDisplay`
            The multi-fold ROC curve display.

        See Also
        --------
        roc_curve : Compute Receiver operating characteristic (ROC) curve.
            RocCurveDisplay.from_estimator : ROC Curve visualization given an
            estimator and some data.
        RocCurveDisplay.from_predictions : ROC Curve visualization given the
            probabilities of scores of a classifier.
        roc_auc_score : Compute the area under the ROC curve.

        Examples
        --------
        >>> import matplotlib.pyplot as plt
        >>> from sklearn.datasets import make_classification
        >>> from sklearn.metrics import RocCurveDisplay
        >>> from sklearn.model_selection import cross_validate
        >>> from sklearn.svm import SVC
        >>> X, y = make_classification(random_state=0)
        >>> clf = SVC(random_state=0)
        >>> cv_results = cross_validate(
        ...     clf, X, y, cv=3, return_estimator=True, return_indices=True)
        >>> RocCurveDisplay.from_cv_results(cv_results, X, y)
        <...>
        >>> plt.show()
        """
        pos_label_ = cls._validate_from_cv_results_params(
            cv_results,
            X,
            y,
            sample_weight=sample_weight,
            pos_label=pos_label,
        )

        fpr_folds, tpr_folds, auc_folds = [], [], []
        for estimator, test_indices in zip(
            cv_results["estimator"], cv_results["indices"]["test"]
        ):
            y_true = _safe_indexing(y, test_indices)
            y_pred, _ = _get_response_values_binary(
                estimator,
                _safe_indexing(X, test_indices),
                response_method=response_method,
                pos_label=pos_label_,
            )
            sample_weight_fold = (
                None
                if sample_weight is None
                else _safe_indexing(sample_weight, test_indices)
            )
            fpr, tpr, _ = roc_curve(
                y_true,
                y_pred,
                pos_label=pos_label_,
                sample_weight=sample_weight_fold,
                drop_intermediate=drop_intermediate,
            )
            roc_auc = auc(fpr, tpr)

            fpr_folds.append(fpr)
            tpr_folds.append(tpr)
            auc_folds.append(roc_auc)

        viz = cls(
            fpr=fpr_folds,
            tpr=tpr_folds,
            name=name,
            roc_auc=auc_folds,
            pos_label=pos_label_,
        )
        return viz.plot(
            ax=ax,
            curve_kwargs=curve_kwargs,
            plot_chance_level=plot_chance_level,
            chance_level_kw=chance_level_kwargs,
            despine=despine,
        )
