import scipy as sp

from .._base import _get_response

from .. import det_curve

from ...base import is_classifier
from ...utils import check_matplotlib_support


class DetCurveDisplay:
    """DET curve visualization.

    It is recommend to use :func:`~sklearn.metrics.plot_det_curve` to create a
    visualizer. All parameters are stored as attributes.

    Read more in the :ref:`User Guide <visualizations>`.

    .. versionadded:: 0.24

    Parameters
    ----------
    fpr : ndarray
        False positive rate.

    tpr : ndarray
        True positive rate.

    estimator_name : str, default=None
        Name of estimator. If None, the estimator name is not shown.

    pos_label : str or int, default=None
        The label of the positive class.

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
    plot_det_curve : Plot detection error tradeoff (DET) curve.

    Examples
    --------
    >>> import matplotlib.pyplot as plt  # doctest: +SKIP
    >>> import numpy as np
    >>> from sklearn import metrics
    >>> y = np.array([0, 0, 1, 1])
    >>> pred = np.array([0.1, 0.4, 0.35, 0.8])
    >>> fpr, fnr, thresholds = metrics.det_curve(y, pred)
    >>> display = metrics.DetCurveDisplay(
    ...     fpr=fpr, fnr=fnr, estimator_name='example estimator'
    ... )
    >>> display.plot()  # doctest: +SKIP
    >>> plt.show()      # doctest: +SKIP
    """
    def __init__(self, *, fpr, fnr, estimator_name=None, pos_label=None):
        self.fpr = fpr
        self.fnr = fnr
        self.estimator_name = estimator_name
        self.pos_label = pos_label

    def plot(self, ax=None, *, name=None, **kwargs):
        """Plot visualization.

        Parameters
        ----------
        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        name : str, default=None
            Name of DET curve for labeling. If `None`, use the name of the
            estimator.

        Returns
        -------
        display : :class:`~sklearn.metrics.plot.DetCurveDisplay`
            Object that stores computed values.
        """
        check_matplotlib_support('DetCurveDisplay.plot')

        name = self.estimator_name if name is None else name
        line_kwargs = {} if name is None else {"label": name}
        line_kwargs.update(**kwargs)

        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots()

        self.line_, = ax.plot(
            sp.stats.norm.ppf(self.fpr),
            sp.stats.norm.ppf(self.fnr),
            **line_kwargs,
        )
        info_pos_label = (f" (Positive label: {self.pos_label})"
                          if self.pos_label is not None else "")

        xlabel = "False Positive Rate" + info_pos_label
        ylabel = "False Negative Rate" + info_pos_label
        ax.set(xlabel=xlabel, ylabel=ylabel)

        if "label" in line_kwargs:
            ax.legend(loc="lower right")

        ticks = [0.001, 0.01, 0.05, 0.20, 0.5, 0.80, 0.95, 0.99, 0.999]
        tick_locations = sp.stats.norm.ppf(ticks)
        tick_labels = [
            '{:.0%}'.format(s) if (100*s).is_integer() else '{:.1%}'.format(s)
            for s in ticks
        ]
        ax.set_xticks(tick_locations)
        ax.set_xticklabels(tick_labels)
        ax.set_xlim(-3, 3)
        ax.set_yticks(tick_locations)
        ax.set_yticklabels(tick_labels)
        ax.set_ylim(-3, 3)

        self.ax_ = ax
        self.figure_ = ax.figure
        return self


def plot_det_curve(
    estimator,
    X,
    y,
    *,
    sample_weight=None,
    response_method="auto",
    name=None,
    ax=None,
    pos_label=None,
    **kwargs
):
    """Plot detection error tradeoff (DET) curve.

    Extra keyword arguments will be passed to matplotlib's `plot`.

    Read more in the :ref:`User Guide <visualizations>`.

    .. versionadded:: 0.24

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

    response_method : {'predict_proba', 'decision_function', 'auto'} \
            default='auto'
        Specifies whether to use :term:`predict_proba` or
        :term:`decision_function` as the predicted target response. If set to
        'auto', :term:`predict_proba` is tried first and if it does not exist
        :term:`decision_function` is tried next.

    name : str, default=None
        Name of DET curve for labeling. If `None`, use the name of the
        estimator.

    ax : matplotlib axes, default=None
        Axes object to plot on. If `None`, a new figure and axes is created.

    pos_label : str or int, default=None
        The label of the positive class.
        When `pos_label=None`, if `y_true` is in {-1, 1} or {0, 1},
        `pos_label` is set to 1, otherwise an error will be raised.

    Returns
    -------
    display : :class:`~sklearn.metrics.DetCurveDisplay`
        Object that stores computed values.

    See Also
    --------
    det_curve : Compute error rates for different probability thresholds.
    DetCurveDisplay : DET curve visualization.
    plot_roc_curve : Plot Receiver operating characteristic (ROC) curve.

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
    >>> metrics.plot_det_curve(clf, X_test, y_test)  # doctest: +SKIP
    >>> plt.show()                                   # doctest: +SKIP
    """
    check_matplotlib_support('plot_det_curve')

    if not is_classifier(estimator):
        raise ValueError(
            f"{estimator.__class__.__name__} should be a binary classifier."
        )

    y_pred, pos_label = _get_response(
        estimator, X, y, response_method, pos_label=pos_label
    )

    fpr, fnr, _ = det_curve(
        y, y_pred, pos_label=pos_label, sample_weight=sample_weight,
    )

    name = estimator.__class__.__name__ if name is None else name

    viz = DetCurveDisplay(
        fpr=fpr,
        fnr=fnr,
        estimator_name=name,
        pos_label=pos_label
    )

    return viz.plot(ax=ax, name=name, **kwargs)
