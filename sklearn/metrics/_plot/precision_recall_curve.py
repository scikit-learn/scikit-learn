# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from collections import Counter
from collections.abc import Mapping

from ...utils import _safe_indexing
from ...utils._plotting import (
    _BinaryClassifierCurveDisplayMixin,
    _despine,
    _validate_style_kwargs,
)
from .._ranking import average_precision_score, precision_recall_curve
from ...utils._response import _get_response_values_binary
from ...utils.validation import _num_samples


class PrecisionRecallDisplay(_BinaryClassifierCurveDisplayMixin):
    """Precision Recall visualization.

    It is recommend to use
    :func:`~sklearn.metrics.PrecisionRecallDisplay.from_estimator` or
    :func:`~sklearn.metrics.PrecisionRecallDisplay.from_predictions` to create
    a :class:`~sklearn.metrics.PrecisionRecallDisplay`. All parameters are
    stored as attributes.

    Read more in the :ref:`User Guide <visualizations>`.

    Parameters
    ----------
    precision : ndarray or list of ndarray
        Precision values. When plotting multiple precision-recall curves, `precision`
        and `recall` should be lists of the same length.

    recall : ndarray or list of ndarray
        Recall values. When plotting multiple precision-recall curves, `precision`
        and `recall` should be lists of the same length.

    average_precision : float or list of floats, default=None
        Average precision.  When plotting multiple precision-recall curves, can be
        a list of the same length as `precision` and `recall`.
        If None, no average precision score is shown.

    curve_name : str or list of str, default=None
        Label for the precision-recall curve. For multiple precision-recall curves,
        `curve_name` can be a list of the same length as `precision` and `recall`.
        If None, no name is shown.

    pos_label : int, float, bool or str, default=None
        The class considered as the positive class. If None, the class will not
        be shown in the legend.

        .. versionadded:: 0.24

    prevalence_pos_label : float, default=None
        The prevalence of the positive label. It is used for plotting the
        chance level line. If None, the chance level line will not be plotted
        even if `plot_chance_level` is set to True when plotting.

        .. versionadded:: 1.3

    Attributes
    ----------
    line_ : matplotlib Artist or list of Artists
        Precision recall curve.

    chance_level_ : matplotlib Artist or None
        The chance level line. It is `None` if the chance level is not plotted.

        .. versionadded:: 1.3

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
        curve_name=None,
        pos_label=None,
        prevalence_pos_label=None,
    ):
        # `curve_name` is a more generalizable term (vs `estimator_name`), especially
        # for multi cv curves
        # where we want to name each curve by the fold number.
        # It does mean we have to change `_validate_plot_params`, which means
        # we have to change the name for all other classes with
        # `_BinaryClassifierCurveDisplayMixin`
        # Opted for singular form as `precision`, `recall` etc are all singular too
        self.precision = precision
        self.recall = recall
        self.average_precision = average_precision
        self.curve_name = curve_name
        self.pos_label = pos_label
        self.prevalence_pos_label = prevalence_pos_label

    def _get_default_line_kwargs(self, average_precision, name):
        # I think we could factorize this for all display classes
        default_line_kwargs = {"drawstyle": "steps-post"}
        if average_precision is not None and name is not None:
            default_line_kwargs["label"] = f"{name} (AP = {average_precision:0.2f})"
        elif average_precision is not None:
            default_line_kwargs["label"] = f"AP = {average_precision:0.2f}"
        elif name is not None:
            default_line_kwargs["label"] = name
        return default_line_kwargs

    def plot(
        self,
        ax=None,
        *,
        name=None,
        plot_chance_level=False,
        chance_level_kw=None,
        despine=False,
        fold_line_kw=None,
        **kwargs,
    ):
        """Plot visualization.

        Extra keyword arguments will be passed to matplotlib's `plot`.

        Parameters
        ----------
        ax : Matplotlib Axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        name : str, default=None
            Name of precision recall curve for labeling. If `None`, use
            `self.curve_name` if not `None`, otherwise no labeling is shown.

        despine : bool, default=False
            Whether to remove the top and right spines from the plot.

            .. versionadded:: 1.6

        plot_chance_level : bool, default=False
            Whether to plot the chance level. The chance level is the prevalence
            of the positive label computed from the data passed during
            :meth:`from_estimator` or :meth:`from_predictions` call.

            .. versionadded:: 1.3

        chance_level_kw : dict, default=None
            Keyword arguments to be passed to matplotlib's `plot` for rendering
            the chance level line.

            .. versionadded:: 1.3

        fold_line_kw : dict or list of dict, default=None
            Dictionary with keywords passed to the matplotlib's `plot` function
            to draw the individual precision-recall curves. If a list is provided,
            the parameters are applied to the precision-recall curves of each fold
            sequentially. If a single dictionary is provided, the same parameters
            are applied to all precision-recall curves. Ignored for single curve
            plots.

        **kwargs : dict
            For a single curve plots only, keyword arguments to be passed to
            matplotlib's `plot`. Ignored for multi-curve plots.

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
        name_ = name if name else self.curve_name
        # If multi-curve, ensure all args are of the right length
        req_multi = [
            input for input in (self.precision, self.recall) if isinstance(input, list)
        ]
        if req_multi and (
            (len(req_multi) != 2) or len({len(arg) for arg in req_multi}) > 1
        ):
            raise ValueError(
                "When plotting multiple precision-recall curves, `self.precision` "
                "and `self.recall` should both be lists of the same length."
            )
        n_multi = len(self.precision) if req_multi else None
        if req_multi:
            for name, param in zip(
                ["self.average_precision", "`name` or `self.curve_name`"],
                (self.average_precision, name_),
            ):
                if not (
                    (isinstance(param, list) and len(param) != n_multi)
                    or param is not None
                ):
                    raise ValueError(
                        f"For multi precision-recall curves, {name} must either be "
                        "a list of the same length as `self.precision` and "
                        "`self.recall`, or None."
                    )

        n_multi = len(self.precision) if req_multi else None
        self.ax_, self.figure_, name = self._validate_plot_params(
            ax=ax,
            name=name,
            n_multi=n_multi,
            curve_type="PR",
        )

        if n_multi:
            if fold_line_kw is None:
                fold_line_kw = [
                    {"alpha": 0.5, "color": "tab:blue", "linestyle": "--"}
                ] * n_multi
            elif isinstance(fold_line_kw, Mapping):
                fold_line_kw = [fold_line_kw] * n_multi
            elif len(fold_line_kw) != n_multi:
                raise ValueError(
                    "When `fold_line_kw` is a list, it must have the same length as "
                    "the number of precision-recall curves to be plotted."
                )
            name_ = [name_] * n_multi if name_ is None else name_
            average_precision_ = (
                [None] * n_multi
                if self.average_precision is None
                else self.average_precision
            )
            line_kwargs = []
            for fold_idx, (curve_name, curve_ap) in enumerate(
                zip(name_, average_precision_)
            ):
                default_line_kwargs = self._get_default_line_kwargs(
                    curve_ap, curve_name
                )
                line_kwargs.append(
                    _validate_style_kwargs(default_line_kwargs, fold_line_kw[fold_idx])
                )
        else:
            default_line_kwargs = self._get_default_line_kwargs(
                self.average_precision, name_
            )
            line_kwargs = _validate_style_kwargs(default_line_kwargs, kwargs)

        if n_multi:
            self.line_ = []
            for recall, precision, line_kw in zip(
                self.recall, self.precision, line_kwargs
            ):
                self.line_.extend(self.ax_.plot(recall, precision, **line_kw))
        else:
            (self.line_,) = self.ax_.plot(self.recall, self.precision, **line_kwargs)

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

            default_chance_level_line_kw = {
                "label": f"Chance level (AP = {self.prevalence_pos_label:0.2f})",
                "color": "k",
                "linestyle": "--",
            }

            if chance_level_kw is None:
                chance_level_kw = {}

            chance_level_line_kw = _validate_style_kwargs(
                default_chance_level_line_kw, chance_level_kw
            )

            (self.chance_level_,) = self.ax_.plot(
                (0, 1),
                (self.prevalence_pos_label, self.prevalence_pos_label),
                **chance_level_line_kw,
            )
        else:
            self.chance_level_ = None

        if despine:
            _despine(self.ax_)

        if "label" in line_kwargs or plot_chance_level:
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
        pos_label=None,
        drop_intermediate=False,
        response_method="auto",
        name=None,
        ax=None,
        plot_chance_level=False,
        chance_level_kw=None,
        despine=False,
        **kwargs,
    ):
        """Plot precision-recall curve given an estimator and some data.

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

        pos_label : int, float, bool or str, default=None
            The class considered as the positive class when computing the
            precision and recall metrics. By default, `estimators.classes_[1]`
            is considered as the positive class.

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

        name : str, default=None
            Name for labeling curve. If `None`, no name is used.

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is created.

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
        y_pred, pos_label, name = cls._validate_and_get_response_values(
            estimator,
            X,
            y,
            response_method=response_method,
            pos_label=pos_label,
            name=name,
        )

        return cls.from_predictions(
            y,
            y_pred,
            sample_weight=sample_weight,
            name=name,
            pos_label=pos_label,
            drop_intermediate=drop_intermediate,
            ax=ax,
            plot_chance_level=plot_chance_level,
            chance_level_kw=chance_level_kw,
            despine=despine,
            **kwargs,
        )

    @classmethod
    def from_predictions(
        cls,
        y_true,
        y_pred,
        *,
        sample_weight=None,
        pos_label=None,
        drop_intermediate=False,
        name=None,
        ax=None,
        plot_chance_level=False,
        chance_level_kw=None,
        despine=False,
        **kwargs,
    ):
        """Plot precision-recall curve given binary class predictions.

        Parameters
        ----------
        y_true : array-like of shape (n_samples,)
            True binary labels.

        y_pred : array-like of shape (n_samples,)
            Estimated probabilities or output of decision function.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights.

        pos_label : int, float, bool or str, default=None
            The class considered as the positive class when computing the
            precision and recall metrics.

        drop_intermediate : bool, default=False
            Whether to drop some suboptimal thresholds which would not appear
            on a plotted precision-recall curve. This is useful in order to
            create lighter precision-recall curves.

            .. versionadded:: 1.3

        name : str, default=None
            Name for labeling curve. If `None`, name will be set to
            `"Classifier"`.

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is created.

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
        >>> y_pred = clf.predict_proba(X_test)[:, 1]
        >>> PrecisionRecallDisplay.from_predictions(
        ...    y_test, y_pred)
        <...>
        >>> plt.show()
        """
        pos_label, name = cls._validate_from_predictions_params(
            y_true, y_pred, sample_weight=sample_weight, pos_label=pos_label, name=name
        )

        precision, recall, _ = precision_recall_curve(
            y_true,
            y_pred,
            pos_label=pos_label,
            sample_weight=sample_weight,
            drop_intermediate=drop_intermediate,
        )
        average_precision = average_precision_score(
            y_true, y_pred, pos_label=pos_label, sample_weight=sample_weight
        )

        class_count = Counter(y_true)
        prevalence_pos_label = class_count[pos_label] / sum(class_count.values())

        viz = cls(
            precision=precision,
            recall=recall,
            average_precision=average_precision,
            curve_name=name,
            pos_label=pos_label,
            prevalence_pos_label=prevalence_pos_label,
        )

        return viz.plot(
            ax=ax,
            name=name,
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
        pos_label=None,
        drop_intermediate=True,
        response_method="auto",
        ax=None,
        despine=False,
        fold_name=None,
        fold_line_kw=None,
        plot_chance_level=False,
        chance_level_kw=None,
        prevalence_pos_label=None,
    ):
        """Plot multi-fold precision-recall curves given cross-validation results.

        .. versionadded:: 1.7

        Parameters
        ----------
        cv_results : dict
            Dictionary as returned by :func:`~sklearn.model_selection.cross_validate`
            using `return_estimator=True` and `return_indices=True`.

        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Input values.

        y : array-like of shape (n_samples,)
            Target values.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights.

        pos_label : str or int, default=None
            The class considered as the positive class when computing the precision
            and recall metrics. By default, `estimators.classes_[1]` is considered
            as the positive class.

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

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        despine : bool, default=False
            Whether to remove the top and right spines from the plot.

        fold_name : list of str, default=None
            Name used in the legend for each individual precision-recall curve. If
            `None`, the name will be set to "PR fold {N}" where N is the index of the
            CV fold.

        fold_line_kw : dict or list of dict, default=None
            Dictionary with keywords passed to the matplotlib's `plot` function
            to draw the individual precision-recall curves. If a list is provided, the
            parameters are applied to the precision-recall curves of each CV fold
            sequentially. If a single dictionary is provided, the same
            parameters are applied to all precision-recall curves.

        plot_chance_level : bool, default=False
            Whether to plot the chance level.

        chance_level_kw : dict, default=None
            Keyword arguments to be passed to matplotlib's `plot` for rendering
            the chance level line.

        prevalence_pos_label : float, default=None
            The prevalence of the positive label. It is used for plotting the
            chance level line. If None, the chance level line will not be plotted
            even if `plot_chance_level` is set to True when plotting.

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
            probability thresholds
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
        required_keys = {"estimator", "indices"}
        if not all(key in cv_results for key in required_keys):
            raise ValueError(
                "cv_results does not contain one of the following required keys: "
                f"{required_keys}. Set explicitly the parameters return_estimator=True "
                "and return_indices=True to the function cross_validate."
            )

        train_size, test_size = (
            len(cv_results["indices"]["train"][0]),
            len(cv_results["indices"]["test"][0]),
        )

        if _num_samples(X) != train_size + test_size:
            raise ValueError(
                "X does not contain the correct number of samples. "
                f"Expected {train_size + test_size}, got {_num_samples(X)}."
            )

        if fold_name is None:
            # create iterable of same length as number of precision-recall curves
            fold_name_ = [None] * len(cv_results["estimator"])
        elif fold_name is not None and len(fold_name) != len(cv_results["estimator"]):
            raise ValueError(
                "When `fold_name` is provided, it must have the same length as "
                "the number of precision-recall curves to be plotted. Got "
                f"{len(fold_name)} names instead of {len(cv_results['estimator'])}."
            )
        else:
            fold_name_ = fold_name

        if fold_line_kw is None:
            fold_line_kw = [
                {"alpha": 0.5, "color": "tab:blue", "linestyle": "--"}
            ] * len(cv_results["estimator"])
        elif isinstance(fold_line_kw, Mapping):
            fold_line_kw = [fold_line_kw] * len(cv_results["estimator"])
        elif len(fold_line_kw) != len(cv_results["estimator"]):
            raise ValueError(
                "When `fold_line_kw` is a list, it must have the same length as "
                "the number of precision-recall curves to be plotted."
            )

        precision_all = []
        recall_all = []
        ap_all = []
        for estimator, test_indices, name in zip(
            cv_results["estimator"], cv_results["indices"]["test"], fold_name_
        ):
            y_true = _safe_indexing(y, test_indices)
            y_pred, pos_label_ = _get_response_values_binary(
                estimator,
                _safe_indexing(X, test_indices),
                response_method=response_method,
                pos_label=pos_label,
            )
            # Should we use `_validate_from_predictions_params` here?
            # The check would technically only be needed once though
            precision, recall, _ = precision_recall_curve(
                y_true,
                y_pred,
                pos_label=pos_label,
                sample_weight=sample_weight,
                drop_intermediate=drop_intermediate,
            )
            # Note `pos_label` cannot be `None` (default=1), unlike other metrics
            # such as roc_auc
            average_precision = average_precision_score(
                y_true, y_pred, pos_label=pos_label_, sample_weight=sample_weight
            )
            # Append all
            precision_all.append(precision)
            recall_all.append(recall)
            ap_all.append(average_precision)

        # Used all data provided to compute prevalence here, not sure if this is
        # misleading and thus chance line should be avoided completely in
        # multi-curve plots
        class_count = Counter(y)
        prevalence_pos_label = class_count[pos_label] / sum(class_count.values())

        viz = cls(
            precision=precision_all,
            recall=recall_all,
            average_precision=ap_all,
            curve_name=name,
            pos_label=pos_label,
            prevalence_pos_label=prevalence_pos_label,
        )
        return viz.plot(
            ax=ax,
            name=fold_name_,
            despine=despine,
            plot_chance_level=plot_chance_level,
            chance_level_kw=chance_level_kw,
            fold_line_kw=fold_line_kw,
        )
