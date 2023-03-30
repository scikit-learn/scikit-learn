from collections.abc import Mapping
from functools import cached_property

import numpy as np
from .. import auc
from .. import roc_curve
from .._base import _check_pos_label_consistency

from ...utils import check_matplotlib_support, _safe_indexing
from ...utils._response import _get_response_values_binary
from ...utils.validation import _num_samples


class RocCurveDisplay:
    """ROC Curve visualization.

    It is recommended to use
    :meth:`~sklearn.metrics.RocCurveDisplay.from_cv_results`,
    :meth:`~sklearn.metrics.RocCurveDisplay.from_estimator` or
    :meth:`~sklearn.metrics.RocCurveDisplay.from_predictions` to create
    a :class:`~sklearn.metrics.RocCurveDisplay`. All parameters are
    stored as attributes.

    Read more in the :ref:`User Guide <visualizations>`.

    Parameters
    ----------
    fpr : ndarray
        False positive rate.

    tpr : ndarray
        True positive rate.

    roc_auc : float, default=None
        Area under ROC curve. If None, the roc_auc score is not shown.

    estimator_name : str, default=None
        Name of estimator. If None, the estimator name is not shown.

    pos_label : str or int, default=None
        The class considered as the positive class when computing the roc auc
        metrics. By default, `estimators.classes_[1]` is considered
        as the positive class.

        .. versionadded:: 0.24

    Attributes
    ----------
    line_ : matplotlib Artist
        ROC Curve.

    ax_ : matplotlib Axes
        Axes with ROC Curve.

    figure_ : matplotlib Figure
        Figure containing the curve.

    See Also
    --------
    roc_curve : Compute Receiver operating characteristic (ROC) curve.
    RocCurveDisplay.from_cv_results : ROC Curve visualization given the
        results of a cross-validation.
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
    >>> fpr, tpr, _ = metrics.roc_curve(y, pred)
    >>> roc_auc = metrics.auc(fpr, tpr)
    >>> display = metrics.RocCurveDisplay(fpr=fpr, tpr=tpr,
    ...                                   roc_auc=roc_auc,
    ...                                   estimator_name='example estimator')
    >>> display.plot()
    <...>
    >>> plt.show()
    """

    def __init__(self, *, fpr, tpr, roc_auc=None, estimator_name=None, pos_label=None):
        self.estimator_name = estimator_name
        self.fpr = fpr
        self.tpr = tpr
        self.roc_auc = roc_auc
        self.pos_label = pos_label

    def plot(self, ax=None, *, name=None, **kwargs):
        """Plot visualization.

        Extra keyword arguments will be passed to matplotlib's ``plot``.

        Parameters
        ----------
        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        name : str, default=None
            Name of ROC Curve for labeling. If `None`, use `estimator_name` if
            not `None`, otherwise no labeling is shown.

        **kwargs : dict
            Keyword arguments to be passed to matplotlib's `plot`.

        Returns
        -------
        display : :class:`~sklearn.metrics.plot.RocCurveDisplay`
            Object that stores computed values.
        """
        check_matplotlib_support("RocCurveDisplay.plot")

        name = self.estimator_name if name is None else name

        line_kwargs = {}
        if self.roc_auc is not None and name is not None:
            line_kwargs["label"] = f"{name} (AUC = {self.roc_auc:0.2f})"
        elif self.roc_auc is not None:
            line_kwargs["label"] = f"AUC = {self.roc_auc:0.2f}"
        elif name is not None:
            line_kwargs["label"] = name

        line_kwargs.update(**kwargs)

        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots()

        (self.line_,) = ax.plot(self.fpr, self.tpr, **line_kwargs)
        info_pos_label = (
            f" (Positive label: {self.pos_label})" if self.pos_label is not None else ""
        )

        xlabel = "False Positive Rate" + info_pos_label
        ylabel = "True Positive Rate" + info_pos_label
        ax.set(xlabel=xlabel, ylabel=ylabel)

        if "label" in line_kwargs:
            ax.legend(loc="lower right")

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
        sample_weight=None,
        drop_intermediate=True,
        response_method="auto",
        pos_label=None,
        name=None,
        ax=None,
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

        pos_label : str or int, default=None
            The class considered as the positive class when computing the roc auc
            metrics. By default, `estimators.classes_[1]` is considered
            as the positive class.

        name : str, default=None
            Name of ROC Curve for labeling. If `None`, use the name of the
            estimator.

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is created.

        **kwargs : dict
            Keyword arguments to be passed to matplotlib's `plot`.

        Returns
        -------
        display : :class:`~sklearn.metrics.plot.RocCurveDisplay`
            The ROC Curve display.

        See Also
        --------
        roc_curve : Compute Receiver operating characteristic (ROC) curve.
        RocCurveDisplay.from_cv_results : ROC Curve visualization given the
            results of a cross-validation.
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
        check_matplotlib_support(f"{cls.__name__}.from_estimator")

        name = estimator.__class__.__name__ if name is None else name

        y_pred, pos_label = _get_response_values_binary(
            estimator,
            X,
            response_method=response_method,
            pos_label=pos_label,
        )

        return cls.from_predictions(
            y_true=y,
            y_pred=y_pred,
            sample_weight=sample_weight,
            drop_intermediate=drop_intermediate,
            name=name,
            ax=ax,
            pos_label=pos_label,
            **kwargs,
        )

    @staticmethod
    def _create_display_from_predictions(
        *,
        y_true,
        y_pred,
        pos_label,
        sample_weight,
        drop_intermediate,
        name,
    ):
        """Private function that create the display given the predictions.
        No plotting is intended here.
        """
        fpr, tpr, _ = roc_curve(
            y_true,
            y_pred,
            pos_label=pos_label,
            sample_weight=sample_weight,
            drop_intermediate=drop_intermediate,
        )
        roc_auc = auc(fpr, tpr)

        name = "Classifier" if name is None else name
        pos_label = _check_pos_label_consistency(pos_label, y_true)

        viz = RocCurveDisplay(
            fpr=fpr,
            tpr=tpr,
            roc_auc=roc_auc,
            estimator_name=name,
            pos_label=pos_label,
        )
        return viz

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

        pos_label : str or int, default=None
            The label of the positive class. When `pos_label=None`, if `y_true`
            is in {-1, 1} or {0, 1}, `pos_label` is set to 1, otherwise an
            error will be raised.

        name : str, default=None
            Name of ROC curve for labeling. If `None`, name will be set to
            `"Classifier"`.

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        **kwargs : dict
            Additional keywords arguments passed to matplotlib `plot` function.

        Returns
        -------
        display : :class:`~sklearn.metrics.RocCurveDisplay`
            Object that stores computed values.

        See Also
        --------
        roc_curve : Compute Receiver operating characteristic (ROC) curve.
        RocCurveDisplay.from_cv_results : ROC Curve visualization given the
            results of a cross-validation.
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
        check_matplotlib_support(f"{cls.__name__}.from_predictions")

        return cls._create_display_from_predictions(
            y_true=y_true,
            y_pred=y_pred,
            pos_label=pos_label,
            sample_weight=sample_weight,
            drop_intermediate=drop_intermediate,
            name=name,
        ).plot(ax=ax, name=name, **kwargs)

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
        kind="folds",
        fold_name=None,
        fold_line_kw=None,
        aggregate_name=None,
        aggregate_line_kw=None,
        aggregate_uncertainty_kw=None,
        plot_chance_level=False,
        chance_level_kw=None,
    ):
        """Create a multi-fold ROC curve display given cross-validation results.

        .. versionadded:: 1.3

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

        kind : {"folds", "aggregate", "both"}, default="folds"
            Controls the type of display to show.

            - `"folds"`: show individual ROC curves for each CV fold.
            - `"aggregate"`: show the mean and standard deviation of the ROC
                curves across CV folds.
            - `"both"`: show both individual ROC curves and the mean and
                standard deviation.

        fold_name : list of str, default=None
            Name used in the legend for each individual ROC curve. If `None`,
            the name will be set to "ROC fold #N" where N is the index of the
            CV fold.

        fold_line_kw : dict or list of dict, default=None
            Dictionary with keywords passed to the matplotlib's `plot` function
            to draw the individual ROC curves. If a list is provided, the
            parameters are applied to the ROC curves of each CV fold
            sequentially. If a single dictionary is provided, the same
            parameters are applied to all ROC curves.

        aggregate_name : str, default=None
            Name used in the legend for the mean ROC curve. If `None`, the name
            will be set to "Mean ROC".

        aggregate_line_kw : dict, default=None
            Dictionary with keywords passed to the matplotlib's `plot` function
            to draw the mean ROC curve.

        aggregate_uncertainty_kw : dict, default=None
            Dictionary with keywords passed to the matplotlib's `fill_between`
            function to draw the standard deviation area.

        plot_chance_level : bool, default=False
            Whether to plot the chance level.

        chance_level_kw : dict, default=None
            Keyword arguments to be passed to matplotlib's `plot` for rendering
            the chance level line.

        Returns
        -------
        display : :class:`~sklearn.metrics.MultiRocCurveDisplay`
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
        >>> RocCurveDisplay.from_cv_results(cv_results, X, y, kind="both")
        <...>
        >>> plt.show()
        """
        check_matplotlib_support(f"{cls.__name__}.from_cv_results")

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
            # create an iterable of the same length as the number of ROC curves
            fold_name_ = [None] * len(cv_results["estimator"])
        elif fold_name is not None and len(fold_name) != len(cv_results["estimator"]):
            raise ValueError(
                "When `fold_name` is provided, it must have the same length as "
                f"the number of ROC curves to be plotted. Got {len(fold_name)} names "
                f"instead of {len(cv_results['estimator'])}."
            )
        else:
            fold_name_ = fold_name

        displays = []
        for fold_id, (estimator, test_indices, name) in enumerate(
            zip(cv_results["estimator"], cv_results["indices"]["test"], fold_name_)
        ):
            y_true = _safe_indexing(y, test_indices)
            y_pred = _get_response_values_binary(
                estimator,
                _safe_indexing(X, test_indices),
                response_method=response_method,
                pos_label=pos_label,
            )[0]

            displays.append(
                cls._create_display_from_predictions(
                    y_true=y_true,
                    y_pred=y_pred,
                    sample_weight=sample_weight,
                    drop_intermediate=drop_intermediate,
                    pos_label=pos_label,
                    name=f"ROC fold #{fold_id}" if name is None else name,
                )
            )

        viz = MultiRocCurveDisplay(displays=displays)
        return viz.plot(
            ax=ax,
            kind=kind,
            fold_name=fold_name,
            fold_line_kw=fold_line_kw,
            aggregate_name=aggregate_name,
            aggregate_line_kw=aggregate_line_kw,
            aggregate_uncertainty_kw=aggregate_uncertainty_kw,
            plot_chance_level=plot_chance_level,
            chance_level_kw=chance_level_kw,
        )


class MultiRocCurveDisplay:
    """Multiple ROC curves visualization.

    It is recommended to use
    :meth:`~sklearn.metrics.RocCurveDisplay.from_cv_results` to create a
    :class:`~sklearn.metrics.MultiRocCurveDisplay`. All parameters are stored as
    attributes.

    Read more in the :ref:`User Guide <visualizations>`.

    .. versionadded:: 1.3

    Parameters
    ----------
    displays : list of :class:`~sklearn.metrics.RocCurveDisplay`
        Display for each ROC curve of each CV fold.

    Attributes
    ----------
    ax_ : matplotlib Axes
        Axes with ROC Curve.

    fold_lines_ : ndarray of matplotlib Artist or None
        Artist representing each individual ROC curve of each CV fold. If
        `None`, no individual ROC curves are plotted.

    mean_line_ : matplotlib Artist or None
        Artist representing the mean ROC curve. If `None`, no mean ROC curve is
        plotted.

    std_area_ : matplotlib Artist or None
        Artist representing the standard deviation area. If `None`, no standard
        deviation area is plotted.

    chance_level_ : matplotlib Artist or None
        Artist representing the chance level. If `None`, no chance level is
        plotted.

    figure_ : matplotlib Figure
        Figure containing the curve.

    See Also
    --------
    roc_curve : Compute Receiver operating characteristic (ROC) curve.
    RocCurveDisplay.from_cv_results : ROC Curve visualization given the
        results of a cross-validation.
    RocCurveDisplay.from_estimator : Plot Receiver Operating Characteristic
        (ROC) curve given an estimator and some data.
    RocCurveDisplay.from_predictions : Plot Receiver Operating Characteristic
        (ROC) curve given the true and predicted values.
    roc_auc_score : Compute the area under the ROC curve.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from sklearn.datasets import make_classification
    >>> from sklearn.linear_model import LogisticRegression
    >>> from sklearn.tree import DecisionTreeClassifier
    >>> from sklearn.metrics import RocCurveDisplay, MultiRocCurveDisplay
    >>> X, y = make_classification(random_state=0)
    >>> log_reg = LogisticRegression().fit(X, y)
    >>> tree = DecisionTreeClassifier(random_state=0).fit(X, y)
    >>> displays = [
    ...     RocCurveDisplay.from_estimator(est, X, y) for est in (log_reg, tree)]
    >>> MultiRocCurveDisplay(displays=displays).plot()
    <...>
    >>> plt.show()
    """

    def __init__(self, *, displays):
        self.displays = displays

    @cached_property
    def _means_and_stds(self):
        """Cached property to compute means and standard deviations of the statistics
        provided in all displays.

        The caching mechanism avoids recomputing the statistics at each call of
        :func:`plot`.
        """
        mean_fpr = np.linspace(0, 1, 100)
        tpr_per_fold, auc_per_fold = [], []
        for display in self.displays:
            interpolated_tpr = np.interp(mean_fpr, display.fpr, display.tpr)
            interpolated_tpr[0] = 0.0
            tpr_per_fold.append(interpolated_tpr)
            auc_per_fold.append(display.roc_auc)
        mean_tpr = np.mean(tpr_per_fold, axis=0)
        mean_tpr[-1] = 1.0
        std_tpr = np.std(tpr_per_fold, axis=0)
        mean_auc = auc(mean_fpr, mean_tpr)
        std_auc = np.std(auc_per_fold)
        return mean_fpr, mean_tpr, std_tpr, mean_auc, std_auc

    def plot(
        self,
        ax=None,
        kind="folds",
        fold_name=None,
        fold_line_kw=None,
        aggregate_name=None,
        aggregate_line_kw=None,
        aggregate_uncertainty_kw=None,
        plot_chance_level=False,
        chance_level_kw=None,
    ):
        """Plot visualization.

        Parameters
        ----------
        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        kind : {"folds", "aggregate", "both"}, default="folds"
            Controls the type of display to show.

            - `"folds"`: show individual ROC curves for each CV fold.
            - `"aggregate"`: show the mean and standard deviation of the ROC
                curves across CV folds.
            - `"both"`: show both individual ROC curves and the mean and
                standard deviation.

        fold_name : list of str, default=None
            Name used in the legend for each individual ROC curve. If `None`,
            the name will be set to "ROC fold #N" where N is the index of the
            CV fold.

        fold_line_kw : dict or list of dict, default=None
            Dictionary with keywords passed to the matplotlib's `plot` function
            to draw the individual ROC curves. If a list is provided, the
            parameters are applied to the ROC curves of each CV fold
            sequentially. If a single dictionary is provided, the same
            parameters are applied to all ROC curves.

        aggregate_name : str, default=None
            Name used in the legend for the mean ROC curve. If `None`, the name
            will be set to "Mean ROC".

        aggregate_line_kw : dict, default=None
            Dictionary with keywords passed to the matplotlib's `plot` function
            to draw the mean ROC curve.

        aggregate_uncertainty_kw : dict, default=None
            Dictionary with keywords passed to the matplotlib's `fill_between`
            function to draw the standard deviation area.

        plot_chance_level : bool, default=False
            Whether to plot the chance level.

        chance_level_kw : dict, default=None
            Keyword arguments to be passed to matplotlib's `plot` for rendering
            the chance level line.

        Returns
        -------
        display : :class:`~sklearn.metrics.MultiRocCurveDisplay`
            Object that stores computed values.
        """
        check_matplotlib_support(f"{self.__class__.__name__}.plot")

        accepted_kinds = ("folds", "aggregate", "both")
        if kind not in accepted_kinds:
            raise ValueError(
                f"Parameter `kind` must be one of {accepted_kinds}. Got"
                f" {kind!r} instead."
            )

        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots()

        if kind in ("folds", "both"):
            if fold_line_kw is None:
                fold_line_kw = [
                    {"alpha": 0.5, "color": "tab:blue", "linestyle": "--"}
                ] * len(self.displays)
            elif isinstance(fold_line_kw, Mapping):
                fold_line_kw = [fold_line_kw] * len(self.displays)
            elif len(fold_line_kw) != len(self.displays):
                raise ValueError(
                    "When `fold_line_kw` is a list, it must have the same length as "
                    "the number of ROC curves to be plotted."
                )

            if fold_name is None:
                fold_name = [None] * len(self.displays)
            elif fold_name is not None and len(fold_name) != len(self.displays):
                raise ValueError(
                    "When `fold_name` is provided, it must have the same length as "
                    "the number of ROC curves to be plotted."
                )

            self.fold_lines_ = []
            for display, name, kwargs in zip(self.displays, fold_name, fold_line_kw):
                display.plot(ax=ax, name=name, **kwargs)
                self.fold_lines_.append(display.line_)
            self.fold_lines_ = np.array(self.fold_lines_)
        else:
            self.fold_lines_ = None

        if kind in ("aggregate", "both"):
            if aggregate_line_kw is None:
                aggregate_line_kw = {"color": "tab:blue", "linewidth": 2}
            if aggregate_uncertainty_kw is None:
                aggregate_uncertainty_kw = {"alpha": 0.1, "color": "tab:blue"}

            mean_fpr, mean_tpr, std_tpr, mean_auc, std_auc = self._means_and_stds
            if aggregate_name is None:
                name = f"Mean ROC (AUC = {mean_auc:.2f} +/- {std_auc:.2f})"
            else:
                name = f"{aggregate_name} (AUC = {mean_auc:.2f} +/- {std_auc:.2f})"

            (self.mean_line_,) = ax.plot(
                mean_fpr, mean_tpr, label=name, **aggregate_line_kw
            )

            tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
            tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
            self.std_area_ = ax.fill_between(
                mean_fpr,
                tprs_lower,
                tprs_upper,
                label=r"$\pm$ 1 std. dev.",
                **aggregate_uncertainty_kw,
            )
            legend_title = "Uncertainties via cross-validation"
        else:
            self.mean_line_ = None
            self.std_area_ = None
            legend_title = None

        if plot_chance_level:
            chance_level_line_kw = {
                "label": "Chance level (AUC = 0.5)",
                "color": "k",
                "linestyle": "--",
            }
            if chance_level_kw is not None:
                chance_level_line_kw.update(chance_level_kw)

            (self.chance_level_,) = ax.plot([0, 1], [0, 1], **chance_level_line_kw)
        else:
            self.chance_level_ = None

        ax.legend(loc="lower right", title=legend_title)
        ax.set_aspect("equal")
        ax.set_xlim((-0.01, 1.01))
        ax.set_ylim((-0.01, 1.01))

        self.ax_ = ax
        self.figure_ = ax.figure
        return self
