from .base import _get_response

from .. import auc
from .. import roc_curve

from ...utils import check_matplotlib_support
from ...utils.validation import _deprecate_positional_args


class RocCurveDisplay:
    """ROC Curve visualization.

    It is recommend to use :func:`~sklearn.metrics.plot_roc_curve` to create a
    visualizer. All parameters are stored as attributes.

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
    plot_roc_curve : Plot Receiver operating characteristic (ROC) curve.
    roc_auc_score : Compute the area under the ROC curve.

    Examples
    --------
    >>> import matplotlib.pyplot as plt  # doctest: +SKIP
    >>> import numpy as np
    >>> from sklearn import metrics
    >>> y = np.array([0, 0, 1, 1])
    >>> pred = np.array([0.1, 0.4, 0.35, 0.8])
    >>> fpr, tpr, thresholds = metrics.roc_curve(y, pred)
    >>> roc_auc = metrics.auc(fpr, tpr)
    >>> display = metrics.RocCurveDisplay(fpr=fpr, tpr=tpr, roc_auc=roc_auc,\
                                          estimator_name='example estimator')
    >>> display.plot()  # doctest: +SKIP
    >>> plt.show()      # doctest: +SKIP
    """
    @_deprecate_positional_args
    def __init__(self, *, fpr, tpr,
                 roc_auc=None, estimator_name=None, pos_label=None):
        self.estimator_name = estimator_name
        self.fpr = fpr
        self.tpr = tpr
        self.roc_auc = roc_auc
        self.pos_label = pos_label

    @_deprecate_positional_args
    def plot(self, ax=None, *, name=None, **kwargs):
        """Plot visualization

        Extra keyword arguments will be passed to matplotlib's ``plot``.

        Parameters
        ----------
        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        name : str, default=None
            Name of ROC Curve for labeling. If `None`, use the name of the
            estimator.

        Returns
        -------
        display : :class:`~sklearn.metrics.plot.RocCurveDisplay`
            Object that stores computed values.
        """
        check_matplotlib_support('RocCurveDisplay.plot')

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

        self.line_, = ax.plot(self.fpr, self.tpr, **line_kwargs)
        info_pos_label = (f" (Positive label: {self.pos_label})"
                          if self.pos_label is not None else "")

        xlabel = "False Positive Rate" + info_pos_label
        ylabel = "True Positive Rate" + info_pos_label
        ax.set(xlabel=xlabel, ylabel=ylabel)

        if "label" in line_kwargs:
            ax.legend(loc="lower right")

        self.ax_ = ax
        self.figure_ = ax.figure
        return self


@_deprecate_positional_args
def plot_roc_curve(estimator, X, y, *, sample_weight=None,
                   drop_intermediate=True, response_method="auto",
                   name=None, ax=None, pos_label=None, **kwargs):
    """Plot Receiver operating characteristic (ROC) curve.

    Extra keyword arguments will be passed to matplotlib's `plot`.

    Read more in the :ref:`User Guide <visualizations>`.

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

    drop_intermediate : boolean, default=True
        Whether to drop some suboptimal thresholds which would not appear
        on a plotted ROC curve. This is useful in order to create lighter
        ROC curves.

    response_method : {'predict_proba', 'decision_function', 'auto'} \
    default='auto'
        Specifies whether to use :term:`predict_proba` or
        :term:`decision_function` as the target response. If set to 'auto',
        :term:`predict_proba` is tried first and if it does not exist
        :term:`decision_function` is tried next.

    name : str, default=None
        Name of ROC Curve for labeling. If `None`, use the name of the
        estimator.

    ax : matplotlib axes, default=None
        Axes object to plot on. If `None`, a new figure and axes is created.

    pos_label : str or int, default=None
        The class considered as the positive class when computing the roc auc
        metrics. By default, `estimators.classes_[1]` is considered
        as the positive class.

        .. versionadded:: 0.24

    Returns
    -------
    display : :class:`~sklearn.metrics.RocCurveDisplay`
        Object that stores computed values.

    See Also
    --------
    roc_curve : Compute Receiver operating characteristic (ROC) curve.
    RocCurveDisplay : ROC Curve visualization.
    roc_auc_score : Compute the area under the ROC curve.

    Examples
    --------
    >>> import matplotlib.pyplot as plt  # doctest: +SKIP
    >>> from sklearn import datasets, metrics, model_selection, svm
    >>> X, y = datasets.make_classification(random_state=0)
    >>> X_train, X_test, y_train, y_test = model_selection.train_test_split(
    ...     X, y, random_state=0)
    >>> clf = svm.SVC(random_state=0)
    >>> clf.fit(X_train, y_train)
    SVC(random_state=0)
    >>> metrics.plot_roc_curve(clf, X_test, y_test)  # doctest: +SKIP
    >>> plt.show()                                   # doctest: +SKIP
    """
    check_matplotlib_support('plot_roc_curve')

    y_pred, pos_label = _get_response(
        X, estimator, response_method, pos_label=pos_label)

    fpr, tpr, _ = roc_curve(y, y_pred, pos_label=pos_label,
                            sample_weight=sample_weight,
                            drop_intermediate=drop_intermediate)
    roc_auc = auc(fpr, tpr)

    name = estimator.__class__.__name__ if name is None else name

    viz = RocCurveDisplay(
        fpr=fpr,
        tpr=tpr,
        roc_auc=roc_auc,
        estimator_name=name,
        pos_label=pos_label
    )

    return viz.plot(ax=ax, name=name, **kwargs)
