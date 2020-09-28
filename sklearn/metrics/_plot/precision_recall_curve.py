from .base import _get_response

from .. import average_precision_score
from .. import precision_recall_curve

from ...utils import check_matplotlib_support
from ...utils.validation import _deprecate_positional_args


class PrecisionRecallDisplay:
    """Precision Recall visualization.

    It is recommend to use :func:`~sklearn.metrics.plot_precision_recall_curve`
    to create a visualizer. All parameters are stored as attributes.

    Read more in the :ref:`User Guide <visualizations>`.

    Parameters
    -----------
    precision : ndarray
        Precision values.

    recall : ndarray
        Recall values.

    average_precision : float, default=None
        Average precision. If None, the average precision is not shown.

    estimator_name : str, default=None
        Name of estimator. If None, then the estimator name is not shown.

    pos_label : str or int, default=None
        The class considered as the positive class. If None, the class will not
        be shown in the legend.

        .. versionadded:: 0.24

    Attributes
    ----------
    line_ : matplotlib Artist
        Precision recall curve.

    ax_ : matplotlib Axes
        Axes with precision recall curve.

    figure_ : matplotlib Figure
        Figure containing the curve.

    See Also
    --------
    precision_recall_curve : Compute precision-recall pairs for different
        probability thresholds.
    plot_precision_recall_curve : Plot Precision Recall Curve for binary
        classifiers.

    Examples
    --------
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
    >>> disp.plot() # doctest: +SKIP
    """
    @_deprecate_positional_args
    def __init__(self, precision, recall, *,
                 average_precision=None, estimator_name=None, pos_label=None):
        self.estimator_name = estimator_name
        self.precision = precision
        self.recall = recall
        self.average_precision = average_precision
        self.pos_label = pos_label

    @_deprecate_positional_args
    def plot(self, ax=None, *, name=None, **kwargs):
        """Plot visualization.

        Extra keyword arguments will be passed to matplotlib's `plot`.

        Parameters
        ----------
        ax : Matplotlib Axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        name : str, default=None
            Name of precision recall curve for labeling. If `None`, use the
            name of the estimator.

        **kwargs : dict
            Keyword arguments to be passed to matplotlib's `plot`.

        Returns
        -------
        display : :class:`~sklearn.metrics.PrecisionRecallDisplay`
            Object that stores computed values.
        """
        check_matplotlib_support("PrecisionRecallDisplay.plot")

        name = self.estimator_name if name is None else name

        line_kwargs = {"drawstyle": "steps-post"}
        if self.average_precision is not None and name is not None:
            line_kwargs["label"] = (f"{name} (AP = "
                                    f"{self.average_precision:0.2f})")
        elif self.average_precision is not None:
            line_kwargs["label"] = (f"AP = "
                                    f"{self.average_precision:0.2f}")
        elif name is not None:
            line_kwargs["label"] = name
        line_kwargs.update(**kwargs)

        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots()

        self.line_, = ax.plot(self.recall, self.precision, **line_kwargs)
        info_pos_label = (f" (Positive label: {self.pos_label})"
                          if self.pos_label is not None else "")

        xlabel = "Recall" + info_pos_label
        ylabel = "Precision" + info_pos_label
        ax.set(xlabel=xlabel, ylabel=ylabel)

        if "label" in line_kwargs:
            ax.legend(loc="lower left")

        self.ax_ = ax
        self.figure_ = ax.figure
        return self


@_deprecate_positional_args
def plot_precision_recall_curve(estimator, X, y, *,
                                sample_weight=None, response_method="auto",
                                name=None, ax=None, pos_label=None, **kwargs):
    """Plot Precision Recall Curve for binary classifiers.

    Extra keyword arguments will be passed to matplotlib's `plot`.

    Read more in the :ref:`User Guide <precision_recall_f_measure_metrics>`.

    Parameters
    ----------
    estimator : estimator instance
        Fitted classifier or a fitted :class:`~sklearn.pipeline.Pipeline`
        in which the last estimator is a classifier.

    X : {array-like, sparse matrix} of shape (n_samples, n_features)
        Input values.

    y : array-like of shape (n_samples,)
        Binary target values.

    sample_weight : array-like of shape (n_samples,), default=None
        Sample weights.

    response_method : {'predict_proba', 'decision_function', 'auto'}, \
                      default='auto'
        Specifies whether to use :term:`predict_proba` or
        :term:`decision_function` as the target response. If set to 'auto',
        :term:`predict_proba` is tried first and if it does not exist
        :term:`decision_function` is tried next.

    name : str, default=None
        Name for labeling curve. If `None`, the name of the
        estimator is used.

    ax : matplotlib axes, default=None
        Axes object to plot on. If `None`, a new figure and axes is created.

    pos_label : str or int, default=None
        The class considered as the positive class when computing the precision
        and recall metrics. By default, `estimators.classes_[1]` is considered
        as the positive class.

        .. versionadded:: 0.24

    **kwargs : dict
        Keyword arguments to be passed to matplotlib's `plot`.

    Returns
    -------
    display : :class:`~sklearn.metrics.PrecisionRecallDisplay`
        Object that stores computed values.

    See Also
    --------
    precision_recall_curve : Compute precision-recall pairs for different
        probability thresholds.
    PrecisionRecallDisplay : Precision Recall visualization.
    """
    check_matplotlib_support("plot_precision_recall_curve")

    y_pred, pos_label = _get_response(
        X, estimator, response_method, pos_label=pos_label)

    precision, recall, _ = precision_recall_curve(y, y_pred,
                                                  pos_label=pos_label,
                                                  sample_weight=sample_weight)
    average_precision = average_precision_score(y, y_pred,
                                                pos_label=pos_label,
                                                sample_weight=sample_weight)

    name = name if name is not None else estimator.__class__.__name__

    viz = PrecisionRecallDisplay(
        precision=precision,
        recall=recall,
        average_precision=average_precision,
        estimator_name=name,
        pos_label=pos_label,
    )

    return viz.plot(ax=ax, name=name, **kwargs)
