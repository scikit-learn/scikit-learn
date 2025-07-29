# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from collections import Counter

import numpy as np

from sklearn.metrics._ranking import average_precision_score, precision_recall_curve
from sklearn.utils import _safe_indexing
from sklearn.utils._plotting import (
    _BinaryClassifierCurveDisplayMixin,
    _check_param_lengths,
    _convert_to_list_leaving_none,
    _deprecate_estimator_name,
    _deprecate_y_pred_parameter,
    _despine,
    _validate_style_kwargs,
)
from sklearn.utils._response import _get_response_values_binary


class PrecisionRecallDisplay(_BinaryClassifierCurveDisplayMixin):
    """Precision Recall visualization.

    It is recommended to use
    :func:`~sklearn.metrics.PrecisionRecallDisplay.from_estimator` or
    :func:`~sklearn.metrics.PrecisionRecallDisplay.from_predictions` to create
    a :class:`~sklearn.metrics.PrecisionRecallDisplay`. All parameters are
    stored as attributes.

    For general information regarding `scikit-learn` visualization tools, see
    the :ref:`Visualization Guide <visualizations>`.
    For guidance on interpreting these plots, refer to the :ref:`Model
    Evaluation Guide <precision_recall_f_measure_metrics>`.

    Parameters
    ----------
    precision : ndarray or list of ndarrays
        Precision values. Each ndarray should contain values for a single curve.
        If plotting multiple curves, list should be of same length as
        and `recall`.

        .. versionchanged:: 1.8
            Now accepts a list for plotting multiple curves.

    recall : ndarray or list of ndarrays
        Recall values. Each ndarray should contain values for a single curve.
        If plotting multiple curves, list should be of same length as
        and `precision`.

        .. versionchanged:: 1.8
            Now accepts a list for plotting multiple curves.

    average_precision : float or list of floats, default=None
        Average precision, used for labeling each curve in the legend.
        If plotting multiple curves, should be a list of the same length as `precision`
        and `recall`. If `None`, average precision values are not shown in the legend.

        .. versionchanged:: 1.8
            Now accepts a list for plotting multiple curves.

    name : str or list of str, default=None
        Name for labeling legend entries. The number of legend entries is determined
        by the `curve_kwargs` passed to `plot`, and is not affected by `name`.
        To label each curve, provide a list of strings. To avoid labeling
        individual curves that have the same appearance, this cannot be used in
        conjunction with `curve_kwargs` being a dictionary or None. If a
        string is provided, it will be used to either label the single legend entry
        or if there are multiple legend entries, label each individual curve with
        the same name. If still `None`, no name is shown in the legend.

        .. versionadded:: 1.8

    pos_label : int, float, bool or str, default=None
        The class considered the positive class when precision and recall metrics
        computed. If not `None`, this value is displayed in the x- and y-axes labels.

        .. versionadded:: 0.24

    prevalence_pos_label : float or list of floats, default=None
        The prevalence of the positive label. It is used for plotting the
        chance level lines. If None, no chance level line will be plotted
        even if `plot_chance_level` is set to True when plotting.

        .. versionadded:: 1.3

    estimator_name : str, default=None
        Name of estimator. If None, the estimator name is not shown.

        .. deprecated:: 1.8
            `estimator_name` is deprecated and will be removed in 1.10. Use `name`
            instead.

    Attributes
    ----------
    line_ : matplotlib Artist or list of Artists
        Precision recall curve.

        .. versionchanged:: 1.8
            This attribute can now be a list of Artists, for when multiple curves
            are plotted.

    chance_level_ : matplotlib Artist or list of Artists or None
        Chance level lines. It is `None` if the chance level is not plotted.

        .. versionadded:: 1.3

        .. versionchanged:: 1.7
            This attribute can now be a list of Artists, for when multiple curves
            are plotted.

    ax_ : matplotlib Axes
        Axes with precision recall curve.

    figure_ : matplotlib Figure
        Figure containing the curve.

    See Also
    --------
    precision_recall_curve : Compute precision-recall pairs for different
        probability thresholds.
    PrecisionRecallDisplay.from_estimator : Plot Precision Recall Curve given
        a binary classifier.
    PrecisionRecallDisplay.from_predictions : Plot Precision Recall Curve
        using predictions from a binary classifier.

    Notes
    -----
    The average precision (cf. :func:`~sklearn.metrics.average_precision_score`) in
    scikit-learn is computed without any interpolation. To be consistent with
    this metric, the precision-recall curve is plotted without any
    interpolation as well (step-wise style).

    You can change this style by passing the keyword argument
    `drawstyle="default"` in :meth:`plot`, :meth:`from_estimator`, or
    :meth:`from_predictions`. However, the curve will not be strictly
    consistent with the reported average precision.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from sklearn.datasets import make_classification
    >>> from sklearn.metrics import (precision_recall_curve,
    ...                              PrecisionRecallDisplay)
    >>> from sklearn.model_selection import train_test_split
    >>> from sklearn.svm import SVC
    >>> X, y = make_classification(random_state=0)
    >>> X_train, X_test, y_train, y_test = train_test_split(X, y,
    ...                                                     random_state=0)
    >>> clf = SVC(random_state=0)
    >>> clf.fit(X_train, y_train)
    SVC(random_state=0)
    >>> predictions = clf.predict(X_test)
    >>> precision, recall, _ = precision_recall_curve(y_test, predictions)
    >>> disp = PrecisionRecallDisplay(precision=precision, recall=recall)
    >>> disp.plot()
    <...>
    >>> plt.show()
    """

    def __init__(
        self,
        precision,
        recall,
        *,
        average_precision=None,
        name=None,
        pos_label=None,
        prevalence_pos_label=None,
        estimator_name="deprecated",
    ):
        self.precision = precision
        self.recall = recall
        self.average_precision = average_precision
        self.name = _deprecate_estimator_name(estimator_name, name, "1.8")
        self.pos_label = pos_label
        self.prevalence_pos_label = prevalence_pos_label

    def _validate_plot_params(self, *, ax, name):
        self.ax_, self.figure_, name = super()._validate_plot_params(ax=ax, name=name)

        precision = _convert_to_list_leaving_none(self.precision)
        recall = _convert_to_list_leaving_none(self.recall)
        average_precision = _convert_to_list_leaving_none(self.average_precision)
        prevalence_pos_label = _convert_to_list_leaving_none(self.prevalence_pos_label)
        name = _convert_to_list_leaving_none(name)

        optional = {
            "self.average_precision": average_precision,
            "self.prevalence_pos_label": prevalence_pos_label,
        }
        if isinstance(name, list) and len(name) != 1:
            optional.update({"'name' (or self.name)": name})
        _check_param_lengths(
            required={"self.precision": precision, "self.recall": recall},
            optional=optional,
            class_name="PrecisionRecallDisplay",
        )
        return precision, recall, average_precision, name, prevalence_pos_label

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
        ax : Matplotlib Axes, default=None
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
            the same name. If `None`, set to `name` provided at `PrecisionRecallDisplay`
            initialization. If still `None`, no name is shown in the legend.

            .. versionchanged:: 1.8
                Now accepts a list for plotting multiple curves.

        curve_kwargs : dict or list of dict, default=None
            Keywords arguments to be passed to matplotlib's `plot` function
            to draw individual precision-recall curves. For single curve plotting, this
            should be a dictionary. For multi-curve plotting, if a list is provided,
            the parameters are applied to each precision-recall curve
            sequentially and a legend entry is added for each curve.
            If a single dictionary is provided, the same parameters are applied
            to all curves and a single legend entry for all curves is added,
            labeled with the mean average precision.

            .. versionadded:: 1.8

        plot_chance_level : bool, default=False
            Whether to plot the chance level. The chance level is the prevalence
            of the positive label computed from the data passed during
            :meth:`from_estimator` or :meth:`from_predictions` call.

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

            .. deprecated:: 1.8
                kwargs is deprecated and will be removed in 1.10. Pass matplotlib
                arguments to `curve_kwargs` as a dictionary instead.

        Returns
        -------
        display : :class:`~sklearn.metrics.PrecisionRecallDisplay`
            Object that stores computed values.

        Notes
        -----
        The average precision (cf. :func:`~sklearn.metrics.average_precision_score`)
        in scikit-learn is computed without any interpolation. To be consistent
        with this metric, the precision-recall curve is plotted without any
        interpolation as well (step-wise style).

        You can change this style by passing the keyword argument
        `drawstyle="default"`. However, the curve will not be strictly
        consistent with the reported average precision.
        """
        precision, recall, average_precision, name, prevalence_pos_label = (
            self._validate_plot_params(ax=ax, name=name)
        )
        n_curves = len(precision)
        if not isinstance(curve_kwargs, list) and n_curves > 1:
            if average_precision:
                legend_metric = {
                    "mean": np.mean(average_precision),
                    "std": np.std(average_precision),
                }
            else:
                legend_metric = {"mean": None, "std": None}
        else:
            average_precision = (
                average_precision
                if average_precision is not None
                else [None] * n_curves
            )
            legend_metric = {"metric": average_precision}

        curve_kwargs = self._validate_curve_kwargs(
            n_curves,
            name,
            legend_metric,
            "AP",
            curve_kwargs=curve_kwargs,
            default_curve_kwargs={"drawstyle": "steps-post"},
            removed_version="1.10",
            **kwargs,
        )
        self.line_ = []
        for recall_val, precision_val, curve_kwarg in zip(
            recall, precision, curve_kwargs
        ):
            self.line_.extend(self.ax_.plot(recall_val, precision_val, **curve_kwarg))
        # Return single artist if only one curve is plotted
        if len(self.line_) == 1:
            self.line_ = self.line_[0]

        info_pos_label = (
            f" (Positive label: {self.pos_label})" if self.pos_label is not None else ""
        )

        xlabel = "Recall" + info_pos_label
        ylabel = "Precision" + info_pos_label
        self.ax_.set(
            xlabel=xlabel,
            xlim=(-0.01, 1.01),
            ylabel=ylabel,
            ylim=(-0.01, 1.01),
            aspect="equal",
        )

        if plot_chance_level:
            if self.prevalence_pos_label is None:
                raise ValueError(
                    "You must provide prevalence_pos_label when constructing the "
                    "PrecisionRecallDisplay object in order to plot the chance "
                    "level line. Alternatively, you may use "
                    "PrecisionRecallDisplay.from_estimator or "
                    "PrecisionRecallDisplay.from_predictions "
                    "to automatically set prevalence_pos_label"
                )

            default_chance_level_kwargs = {
                "color": "k",
                "linestyle": "--",
            }
            if n_curves > 1:
                default_chance_level_kwargs["alpha"] = 0.3

            if chance_level_kw is None:
                chance_level_kw = {}

            chance_level_line_kw = _validate_style_kwargs(
                default_chance_level_kwargs, chance_level_kw
            )
            self.chance_level_ = []
            for prevalence_pos_label_val in prevalence_pos_label:
                self.chance_level_.extend(
                    self.ax_.plot(
                        (0, 1),
                        (prevalence_pos_label_val, prevalence_pos_label_val),
                        **chance_level_line_kw,
                    )
                )

            if len(self.chance_level_) == 1:
                # Return single artist if only one curve is plotted
                self.chance_level_ = self.chance_level_[0]
                if "label" not in chance_level_line_kw:
                    self.chance_level_.set_label(
                        f"Chance level (AP = {prevalence_pos_label[0]:0.2f})"
                    )
            else:
                if "label" not in chance_level_line_kw:
                    # Only label first curve with mean AP, to get single legend entry
                    self.chance_level_[0].set_label(
                        f"Chance level (AP = {np.mean(prevalence_pos_label):0.2f} "
                        f"+/- {np.std(prevalence_pos_label):0.2f})"
                    )
        else:
            self.chance_level_ = None

        if despine:
            _despine(self.ax_)

        # Note: if 'label' present in one `line_kwargs`, it should be present in all
        if curve_kwargs[0].get("label") is not None or (
            plot_chance_level and chance_level_kw.get("label") is not None
        ):
            self.ax_.legend(loc="lower left")

        return self

    @classmethod
    def from_estimator(
        cls,
        estimator,
        X,
        y,
        *,
        sample_weight=None,
        drop_intermediate=False,
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
        """Plot precision-recall curve given an estimator and some data.

        For general information regarding `scikit-learn` visualization tools, see
        the :ref:`Visualization Guide <visualizations>`.
        For guidance on interpreting these plots, refer to the :ref:`Model
        Evaluation Guide <precision_recall_f_measure_metrics>`.

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

        drop_intermediate : bool, default=False
            Whether to drop some suboptimal thresholds which would not appear
            on a plotted precision-recall curve. This is useful in order to
            create lighter precision-recall curves.

            .. versionadded:: 1.3

        response_method : {'predict_proba', 'decision_function', 'auto'}, \
            default='auto'
            Specifies whether to use :term:`predict_proba` or
            :term:`decision_function` as the target response. If set to 'auto',
            :term:`predict_proba` is tried first and if it does not exist
            :term:`decision_function` is tried next.

        pos_label : int, float, bool or str, default=None
            The class considered as the positive class when computing the
            precision and recall metrics. By default, `estimators.classes_[1]`
            is considered as the positive class.

        name : str, default=None
            Name for labeling curve. If `None`, no name is used.

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is created.

        curve_kwargs : dict, default=None
            Keywords arguments to be passed to matplotlib's `plot` function.

            .. versionadded:: 1.8

        plot_chance_level : bool, default=False
            Whether to plot the chance level. The chance level is the prevalence
            of the positive label computed from the data passed during
            :meth:`from_estimator` or :meth:`from_predictions` call.

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

            .. deprecated:: 1.8
                kwargs is deprecated and will be removed in 1.10. Pass matplotlib
                arguments to `curve_kwargs` as a dictionary instead.

        Returns
        -------
        display : :class:`~sklearn.metrics.PrecisionRecallDisplay`

        See Also
        --------
        PrecisionRecallDisplay.from_predictions : Plot precision-recall curve
            using estimated probabilities or output of decision function.

        Notes
        -----
        The average precision (cf. :func:`~sklearn.metrics.average_precision_score`)
        in scikit-learn is computed without any interpolation. To be consistent
        with this metric, the precision-recall curve is plotted without any
        interpolation as well (step-wise style).

        You can change this style by passing the keyword argument
        `drawstyle="default"`. However, the curve will not be strictly
        consistent with the reported average precision.

        Examples
        --------
        >>> import matplotlib.pyplot as plt
        >>> from sklearn.datasets import make_classification
        >>> from sklearn.metrics import PrecisionRecallDisplay
        >>> from sklearn.model_selection import train_test_split
        >>> from sklearn.linear_model import LogisticRegression
        >>> X, y = make_classification(random_state=0)
        >>> X_train, X_test, y_train, y_test = train_test_split(
        ...         X, y, random_state=0)
        >>> clf = LogisticRegression()
        >>> clf.fit(X_train, y_train)
        LogisticRegression()
        >>> PrecisionRecallDisplay.from_estimator(
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
            y,
            y_score,
            sample_weight=sample_weight,
            drop_intermediate=drop_intermediate,
            pos_label=pos_label,
            name=name,
            ax=ax,
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
        drop_intermediate=False,
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
        """Plot precision-recall curve given binary class predictions.

        For general information regarding `scikit-learn` visualization tools, see
        the :ref:`Visualization Guide <visualizations>`.
        For guidance on interpreting these plots, refer to the :ref:`Model
        Evaluation Guide <precision_recall_f_measure_metrics>`.

        Parameters
        ----------
        y_true : array-like of shape (n_samples,)
            True binary labels.

        y_score : array-like of shape (n_samples,)
            Estimated probabilities or output of decision function.

            .. versionadded:: 1.8
                `y_pred` has been renamed to `y_score`.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights.

        drop_intermediate : bool, default=False
            Whether to drop some suboptimal thresholds which would not appear
            on a plotted precision-recall curve. This is useful in order to
            create lighter precision-recall curves.

            .. versionadded:: 1.3

        pos_label : int, float, bool or str, default=None
            The class considered as the positive class when computing the
            precision and recall metrics. When `pos_label=None`, if `y_true` is
            in {-1, 1} or {0, 1}, `pos_label` is set to 1, otherwise an error
            will be raised.

        name : str, default=None
            Name for labeling curve. If `None`, name will be set to
            `"Classifier"`.

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is created.

        curve_kwargs : dict, default=None
            Keywords arguments to be passed to matplotlib's `plot` function.

            .. versionadded:: 1.8

        plot_chance_level : bool, default=False
            Whether to plot the chance level. The chance level is the prevalence
            of the positive label computed from the data passed during
            :meth:`from_estimator` or :meth:`from_predictions` call.

            .. versionadded:: 1.3

        chance_level_kw : dict, default=None
            Keyword arguments to be passed to matplotlib's `plot` for rendering
            the chance level line.

            .. versionadded:: 1.3

        despine : bool, default=False
            Whether to remove the top and right spines from the plot.

            .. versionadded:: 1.6

        y_pred : array-like of shape (n_samples,)
            Estimated probabilities or output of decision function.

            .. deprecated:: 1.8
                `y_pred` is deprecated and will be removed in 1.10. Use
                `y_score` instead.

        **kwargs : dict
            Keyword arguments to be passed to matplotlib's `plot`.

            .. deprecated:: 1.8
                kwargs is deprecated and will be removed in 1.10. Pass matplotlib
                arguments to `curve_kwargs` as a dictionary instead.

        Returns
        -------
        display : :class:`~sklearn.metrics.PrecisionRecallDisplay`

        See Also
        --------
        PrecisionRecallDisplay.from_estimator : Plot precision-recall curve
            using an estimator.

        Notes
        -----
        The average precision (cf. :func:`~sklearn.metrics.average_precision_score`)
        in scikit-learn is computed without any interpolation. To be consistent
        with this metric, the precision-recall curve is plotted without any
        interpolation as well (step-wise style).

        You can change this style by passing the keyword argument
        `drawstyle="default"`. However, the curve will not be strictly
        consistent with the reported average precision.

        Examples
        --------
        >>> import matplotlib.pyplot as plt
        >>> from sklearn.datasets import make_classification
        >>> from sklearn.metrics import PrecisionRecallDisplay
        >>> from sklearn.model_selection import train_test_split
        >>> from sklearn.linear_model import LogisticRegression
        >>> X, y = make_classification(random_state=0)
        >>> X_train, X_test, y_train, y_test = train_test_split(
        ...         X, y, random_state=0)
        >>> clf = LogisticRegression()
        >>> clf.fit(X_train, y_train)
        LogisticRegression()
        >>> y_score = clf.predict_proba(X_test)[:, 1]
        >>> PrecisionRecallDisplay.from_predictions(
        ...    y_test, y_score)
        <...>
        >>> plt.show()
        """
        y_score = _deprecate_y_pred_parameter(y_score, y_pred, "1.8")
        pos_label, name = cls._validate_from_predictions_params(
            y_true, y_score, sample_weight=sample_weight, pos_label=pos_label, name=name
        )

        precision, recall, _ = precision_recall_curve(
            y_true,
            y_score,
            pos_label=pos_label,
            sample_weight=sample_weight,
            drop_intermediate=drop_intermediate,
        )
        average_precision = average_precision_score(
            y_true, y_score, pos_label=pos_label, sample_weight=sample_weight
        )

        class_count = Counter(y_true)
        prevalence_pos_label = class_count[pos_label] / sum(class_count.values())

        viz = cls(
            precision=precision,
            recall=recall,
            average_precision=average_precision,
            name=name,
            pos_label=pos_label,
            prevalence_pos_label=prevalence_pos_label,
        )

        return viz.plot(
            ax=ax,
            name=name,
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
        name=None,
        ax=None,
        curve_kwargs=None,
        plot_chance_level=False,
        chance_level_kwargs=None,
        despine=False,
    ):
        """Plot multi-fold precision-recall curves given cross-validation results.

        .. versionadded:: 1.8

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
            on a plotted precision-recall curve. This is useful in order to
            create lighter precision-recall curves.

        response_method : {'predict_proba', 'decision_function', 'auto'} \
                default='auto'
            Specifies whether to use :term:`predict_proba` or
            :term:`decision_function` as the target response. If set to 'auto',
            :term:`predict_proba` is tried first and if it does not exist
            :term:`decision_function` is tried next.

        pos_label : int, float, bool or str, default=None
            The class considered as the positive class when computing the precision
            and recall metrics. By default, `estimators.classes_[1]` is considered
            as the positive class.

        name : str or list of str, default=None
            Name for labeling legend entries. The number of legend entries
            is determined by `curve_kwargs`, and is not affected by `name`.
            To label each curve, provide a list of strings. To avoid labeling
            individual curves that have the same appearance, this cannot be used in
            conjunction with `curve_kwargs` being a dictionary or None. If a
            string is provided, it will be used to either label the single legend entry
            or if there are multiple legend entries, label each individual curve with
            the same name. If `None`, no name is shown in the legend.

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        curve_kwargs : dict or list of dict, default=None
            Dictionary with keywords passed to the matplotlib's `plot` function
            to draw the individual precision-recall curves. If a list is provided, the
            parameters are applied to the precision-recall curves of each CV fold
            sequentially. If a single dictionary is provided, the same
            parameters are applied to all precision-recall curves.

        plot_chance_level : bool, default=False
            Whether to plot the chance level lines.

        chance_level_kwargs : dict, default=None
            Keyword arguments to be passed to matplotlib's `plot` for rendering
            the chance level lines.

        despine : bool, default=False
            Whether to remove the top and right spines from the plot.

        Returns
        -------
        display : :class:`~sklearn.metrics.PrecisionRecallDisplay`

        See Also
        --------
        PrecisionRecallDisplay.from_predictions : Plot precision-recall curve
            using estimated probabilities or output of decision function.
        PrecisionRecallDisplay.from_estimator : Plot precision-recall curve
            using an estimator.
        precision_recall_curve : Compute precision-recall pairs for different
            probability thresholds.
        average_precision_score : Compute average precision (AP) from prediction scores.

        Notes
        -----
        The average precision (cf. :func:`~sklearn.metrics.average_precision_score`)
        in scikit-learn is computed without any interpolation. To be consistent
        with this metric, the precision-recall curve is plotted without any
        interpolation as well (step-wise style).

        You can change this style by passing the keyword argument
        `drawstyle="default"`. However, the curve will not be strictly
        consistent with the reported average precision.

        Examples
        --------
        >>> import matplotlib.pyplot as plt
        >>> from sklearn.datasets import make_classification
        >>> from sklearn.metrics import PrecisionRecallDisplay
        >>> from sklearn.model_selection import cross_validate
        >>> from sklearn.svm import SVC
        >>> X, y = make_classification(random_state=0)
        >>> clf = SVC(random_state=0)
        >>> cv_results = cross_validate(
        ...     clf, X, y, cv=3, return_estimator=True, return_indices=True)
        >>> PrecisionRecallDisplay.from_cv_results(cv_results, X, y)
        <...>
        >>> plt.show()
        """
        cls._validate_from_cv_results_params(
            cv_results, X, y, sample_weight=sample_weight
        )

        precision_folds, recall_folds, ap_folds, prevalence_pos_label_folds = (
            [],
            [],
            [],
            [],
        )
        for estimator, test_indices in zip(
            cv_results["estimator"], cv_results["indices"]["test"]
        ):
            y_true = _safe_indexing(y, test_indices)
            y_pred, pos_label_ = _get_response_values_binary(
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
            precision, recall, _ = precision_recall_curve(
                y_true,
                y_pred,
                pos_label=pos_label_,
                sample_weight=sample_weight_fold,
                drop_intermediate=drop_intermediate,
            )
            # Note `pos_label` cannot be `None` (default=1), unlike other metrics
            # such as roc_auc
            average_precision = average_precision_score(
                y_true, y_pred, pos_label=pos_label_, sample_weight=sample_weight_fold
            )
            class_count = Counter(y_true)
            # would `y_true.shape[0]` be faster?
            prevalence_pos_label = class_count[pos_label_] / sum(class_count.values())

            precision_folds.append(precision)
            recall_folds.append(recall)
            ap_folds.append(average_precision)
            prevalence_pos_label_folds.append(prevalence_pos_label)

        viz = cls(
            precision=precision_folds,
            recall=recall_folds,
            average_precision=ap_folds,
            pos_label=pos_label_,
            prevalence_pos_label=prevalence_pos_label_folds,
        )
        return viz.plot(
            ax=ax,
            name=name,
            curve_kwargs=curve_kwargs,
            plot_chance_level=plot_chance_level,
            chance_level_kw=chance_level_kwargs,
            despine=despine,
        )
