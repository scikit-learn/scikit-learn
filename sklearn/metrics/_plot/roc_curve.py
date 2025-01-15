# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ...utils import _safe_indexing, deprecated
from ...utils._plotting import (
    _BinaryClassifierCurveDisplayMixin,
    _check_param_lengths,
    _deprecate_singular,
    _despine,
    _process_fold_names_line_kwargs,
    _validate_style_kwargs,
)
from ...utils._response import _get_response_values_binary
from ...utils.validation import _num_samples
from .._ranking import auc, roc_curve


class RocCurveDisplay(_BinaryClassifierCurveDisplayMixin):
    """ROC Curve visualization.

    It is recommend to use
    :func:`~sklearn.metrics.RocCurveDisplay.from_estimator` or
    :func:`~sklearn.metrics.RocCurveDisplay.from_predictions` or
    :func:`~sklearn.metrics.RocCurveDisplay.from_cv_results` to create
    a :class:`~sklearn.metrics.RocCurveDisplay`. All parameters are
    stored as attributes.

    Read more in the :ref:`User Guide <visualizations>`.

    Parameters
    ----------
    fprs : list of ndarray
        False positive rates. Each ndarray should contain values for a single curve.
        If plotting multiple curves, list should be of same length as
        and `tprs`.

    tprs : list of ndarray
        True positive rates. Each ndarray should contain values for a single curve.
        If plotting multiple curves, list should be of same length as
        and `fprs`.

    roc_aucs : list of floats, default=None
        Area under ROC curve. Should be list of the same length as `fprs` and
        `tprs` or None, in which case no area under ROC curve score is shown.

    names : list of str, default=None
        Label for the ROC curve. Should be list of the same length as
        `fprs` and `tprs` or None, in which case no name is shown.

    pos_label : int, float, bool or str, default=None
        The class considered as the positive class when computing the roc auc
        metrics. By default, `estimators.classes_[1]` is considered
        as the positive class.

        .. versionadded:: 0.24

    fpr : ndarray or list of ndarray
        False positive rate. When plotting multiple ROC curves, `fpr` and `tpr` should
        lists of the same length.

        .. deprecated:: 1.7
            `fpr` is deprecated in 1.7 and will be removed in 1.9.
            Use `fprs` instead.

    tpr : ndarray or list of ndarray
        True positive rate. When plotting multiple ROC curves, `fpr` and `tpr` should
        lists of the same length.

        .. deprecated:: 1.7
            `tpr` is deprecated in 1.7 and will be removed in 1.9.
            Use `tprs` instead.

    roc_auc : float or list of floats, default=None
        Area under ROC curve. When plotting multiple ROC curves, can be a list
        of the same length as `fpr` and `tpr`.
        If None, no roc_auc score is shown.

        .. deprecated:: 1.7
            `roc_auc` is deprecated in 1.7 and will be removed in 1.9.
            Use `roc_aucs` instead.

    name : str, default=None
        Label for the ROC curve. For multiple ROC curves, `name` can be a list
        of the same length as `tpr` and `fpr`.
        If None, no name is shown.

        .. deprecated:: 1.7
            `name` is deprecated in 1.7 and will be removed in 1.9.
            Use `names` instead.

    Attributes
    ----------
    lines_ : list of matplotlib Artists
        ROC Curves.

    chance_level_ : matplotlib Artist or None
        The chance level line. It is `None` if the chance level is not plotted.

        .. versionadded:: 1.3

    ax_ : matplotlib Axes
        Axes with ROC Curve.

    figure_ : matplotlib Figure
        Figure containing the curve.

    line_ : matplotlib Artist
        ROC Curve.

        .. deprecated:: 1.7
            `line_` is deprecated in 1.7 and will be removed in 1.9. Use `lines_`
            instead.

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
    >>> y = np.array([0, 0, 1, 1])
    >>> pred = np.array([0.1, 0.4, 0.35, 0.8])
    >>> fpr, tpr, thresholds = metrics.roc_curve(y, pred)
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
        fprs,
        tprs,
        roc_aucs=None,
        names=None,
        pos_label=None,
        fpr="deprecated",
        tpr="deprecated",
        roc_auc="deprecated",
        name="deprecated",
    ):
        self.fprs = _deprecate_singular(fpr, fprs, "fpr")
        self.tprs = _deprecate_singular(tpr, tprs, "tpr")
        self.roc_aucs = _deprecate_singular(roc_auc, roc_aucs, "roc_auc")
        self.names = _deprecate_singular(name, names, "name")
        self.pos_label = pos_label

    def plot(
        self,
        ax=None,
        *,
        names=None,
        plot_chance_level=False,
        chance_level_kw=None,
        despine=False,
        fold_line_kwargs=None,
        name="deprecated",
        **kwargs,
    ):
        """Plot visualization.

        For single curve plots, extra keyword arguments will be passed to
        matplotlib's ``plot``.

        Parameters
        ----------
        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        names : list of str, default=None
            Names of each ROC curve for labeling each curve in the legend.
            If `None`, use `names` provided at `RocCurveDisplay` initialization. If
            also not provided at initialization, no legend is added.

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

        fold_line_kwargs : dict or list of dict, default=None
            Dictionary with keywords passed to the matplotlib's `plot` function
            to draw the individual ROC curves. If a list is provided, the
            parameters are applied to the ROC curves sequentially. If a single
            dictionary is provided, the same parameters are applied to all ROC
            curves. Ignored for single curve plots.

            .. versionadded:: 1.7

        name : str, default=None
            Name of ROC Curve for labeling. If `None`, use `estimator_name` if
            not `None`, otherwise no labeling is shown.

            .. deprecated:: 1.7
                `name` is deprecated in 1.7 and will be removed in 1.9.
                Use `names` instead.

        **kwargs : dict
            For a single curve plots only, keyword arguments to be passed to
            matplotlib's `plot`. Ignored for multi-curve plots.
            (Note req for backwards compat, maybe not ideal?)

        Returns
        -------
        display : :class:`~sklearn.metrics.RocCurveDisplay`
            Object that stores computed values.
        """
        names = _deprecate_singular(name, names, "name")
        # Not sure about this, as ideally we would check params are correct first??
        self.ax_, self.figure_, names_ = self._validate_plot_params(ax=ax, name=names)
        _check_param_lengths(
            {"self.fprs": self.fprs, "self.tprs": self.tprs},
            {"roc_aucs": self.roc_aucs, "self.names (or names from `plot`)": names_},
            "RocCurveDisplay",
        )

        n_curves = len(self.fprs)
        line_kwargs = self._get_line_kwargs(
            n_curves, names_, self.roc_aucs, fold_line_kwargs, **kwargs
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

        self.lines_ = []
        for fpr, tpr, line_kw in zip(self.fprs, self.tprs, line_kwargs):
            self.lines_.extend(self.ax_.plot(fpr, tpr, **line_kw))

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

        if (
            line_kwargs[0].get("label") is not None
            or chance_level_kw.get("label") is not None
        ):
            self.ax_.legend(loc="lower right")

        return self

    # TODO(1.9): Remove
    # Is it worth adding a global ignore for mypy error?
    # mypy error: Decorated property not supported
    @deprecated(  # type: ignore
        "Attribute `line_` is deprecated in 1.7 and will be removed in "
        "1.9. Use `lines_` instead."
    )
    @property
    def line_(self):
        return self.lines_[0]

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
        plot_chance_level=False,
        chance_level_kw=None,
        despine=False,
        **kwargs,
    ):
        """Create a ROC Curve display from an estimator.

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
            The class considered as the positive class when computing the roc auc
            metrics. By default, `estimators.classes_[1]` is considered
            as the positive class.

        name : str, default=None
            Name of ROC Curve for labeling. If `None`, use the name of the
            estimator.

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is created.

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
        y_pred, pos_label, name = cls._validate_and_get_response_values(
            estimator,
            X,
            y,
            response_method=response_method,
            pos_label=pos_label,
            name=name,
        )

        return cls.from_predictions(
            y_true=y,
            y_pred=y_pred,
            sample_weight=sample_weight,
            drop_intermediate=drop_intermediate,
            name=name,
            ax=ax,
            pos_label=pos_label,
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
        drop_intermediate=True,
        pos_label=None,
        name=None,
        ax=None,
        plot_chance_level=False,
        chance_level_kw=None,
        despine=False,
        **kwargs,
    ):
        """Plot ROC curve given the true and predicted values.

        Read more in the :ref:`User Guide <visualizations>`.

        .. versionadded:: 1.0

        Parameters
        ----------
        y_true : array-like of shape (n_samples,)
            True labels.

        y_pred : array-like of shape (n_samples,)
            Target scores, can either be probability estimates of the positive
            class, confidence values, or non-thresholded measure of decisions
            (as returned by “decision_function” on some classifiers).

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights.

        drop_intermediate : bool, default=True
            Whether to drop some suboptimal thresholds which would not appear
            on a plotted ROC curve. This is useful in order to create lighter
            ROC curves.

        pos_label : int, float, bool or str, default=None
            The label of the positive class. When `pos_label=None`, if `y_true`
            is in {-1, 1} or {0, 1}, `pos_label` is set to 1, otherwise an
            error will be raised.

        name : str, default=None
            Name of ROC curve for labeling. If `None`, name will be set to
            `"Classifier"`.

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

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
            Additional keywords arguments passed to matplotlib `plot` function.

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
        >>> y_pred = clf.decision_function(X_test)
        >>> RocCurveDisplay.from_predictions(
        ...    y_test, y_pred)
        <...>
        >>> plt.show()
        """
        pos_label_validated, name = cls._validate_from_predictions_params(
            y_true, y_pred, sample_weight=sample_weight, pos_label=pos_label, name=name
        )

        fpr, tpr, _ = roc_curve(
            y_true,
            y_pred,
            pos_label=pos_label,
            sample_weight=sample_weight,
            drop_intermediate=drop_intermediate,
        )
        roc_auc = auc(fpr, tpr)

        viz = cls(
            fprs=[fpr],
            tprs=[tpr],
            roc_aucs=[roc_auc],
            names=[name],
            pos_label=pos_label_validated,
        )

        return viz.plot(
            ax=ax,
            # Should we provide `name` to both `cls` and `plot` or just `cls`?
            names=[name],
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
        fold_names=None,
        fold_line_kwargs=None,
        plot_chance_level=False,
        chance_level_kw=None,
    ):
        """Create a multi-fold ROC curve display given cross-validation results.

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

        pos_label : str or int, default=None
            The class considered as the positive class when computing the roc auc
            metrics. By default, `estimators.classes_[1]` is considered
            as the positive class.

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        fold_names : list of str, default=None
            Names used in the legend for each individual ROC curve. If `None`,
            the name will be set to "ROC fold <N>" where N is the index of the
            CV fold.

        fold_line_kwargs : dict or list of dict, default=None
            Dictionary with keywords passed to the matplotlib's `plot` function
            to draw the individual ROC curves. If a list is provided, the
            parameters are applied to the ROC curves of each CV fold
            sequentially. If a single dictionary is provided, the same
            parameters are applied to all ROC curves.

        plot_chance_level : bool, default=False
            Whether to plot the chance level.

        chance_level_kw : dict, default=None
            Keyword arguments to be passed to matplotlib's `plot` for rendering
            the chance level line.

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

        fold_names_, fold_line_kwargs_ = _process_fold_names_line_kwargs(
            len(cv_results["estimator"]), fold_names, fold_line_kwargs
        )

        fpr_all = []
        tpr_all = []
        auc_all = []
        for estimator, test_indices in zip(
            cv_results["estimator"], cv_results["indices"]["test"]
        ):
            y_true = _safe_indexing(y, test_indices)
            y_pred = _get_response_values_binary(
                estimator,
                _safe_indexing(X, test_indices),
                response_method=response_method,
                pos_label=pos_label,
            )[0]
            # Should we use `_validate_from_predictions_params` here?
            # The check would technically only be needed once though
            fpr, tpr, _ = roc_curve(
                y_true,
                y_pred,
                pos_label=pos_label,
                sample_weight=sample_weight,
                drop_intermediate=drop_intermediate,
            )
            roc_auc = auc(fpr, tpr)
            # Append all
            fpr_all.append(fpr)
            tpr_all.append(tpr)
            auc_all.append(roc_auc)

        viz = cls(
            fprs=fpr_all,
            tprs=tpr_all,
            names=fold_names_,
            roc_aucs=auc_all,
            pos_label=pos_label,
        )
        return viz.plot(
            ax=ax,
            fold_line_kwargs=fold_line_kwargs_,
            plot_chance_level=plot_chance_level,
            chance_level_kw=chance_level_kw,
        )
