from .. import auc
from .. import roc_curve

from ...utils import check_matplotlib_support


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

    roc_auc : float
        Area under ROC curve.

    estimator_name : str
        Name of estimator.

    Attributes
    ----------
    line_ : matplotlib Artist
        ROC Curve.

    ax_ : matplotlib Axes
        Axes with ROC Curve.

    figure_ : matplotlib Figure
        Figure containing the curve.

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

    def __init__(self, fpr, tpr, roc_auc, estimator_name):
        self.fpr = fpr
        self.tpr = tpr
        self.roc_auc = roc_auc
        self.estimator_name = estimator_name

    def plot(self, ax=None, name=None, **kwargs):
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
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots()

        name = self.estimator_name if name is None else name

        line_kwargs = {
            'label': "{} (AUC = {:0.2f})".format(name, self.roc_auc)
        }
        line_kwargs.update(**kwargs)

        self.line_ = ax.plot(self.fpr, self.tpr, **line_kwargs)[0]
        ax.set_xlabel("False Positive Rate")
        ax.set_ylabel("True Positive Rate")
        ax.legend(loc='lower right')

        self.ax_ = ax
        self.figure_ = ax.figure
        return self


def plot_roc_curve(estimator, X, y, pos_label=None, sample_weight=None,
                   drop_intermediate=True, response_method="auto",
                   name=None, ax=None, **kwargs):
    """Plot Receiver operating characteristic (ROC) curve.

    Extra keyword arguments will be passed to matplotlib's `plot`.

    Read more in the :ref:`User Guide <visualizations>`.

    Parameters
    ----------
    estimator : estimator instance
        Trained binary classifier.

    X : {array-like, sparse matrix} of shape (n_samples, n_features)
        Input values.

    y : array-like of shape (n_samples,)
        Target values.

    pos_label : int or str, default=None
        Label of the positive class.
        By default, pos_label is inferred automatically by taking the last
        class label from the estimator.classes_ attribute.

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

    Returns
    -------
    display : :class:`~sklearn.metrics.RocCurveDisplay`
        Object that stores computed values.

    Examples
    --------
    >>> import matplotlib.pyplot as plt  # doctest: +SKIP
    >>> from sklearn import datasets, metrics, model_selection, svm
    >>> X, y = datasets.make_classification(random_state=0)
    >>> X_train, X_test, y_train, y_test = model_selection.train_test_split(\
            X, y, random_state=0)
    >>> clf = svm.SVC(random_state=0)
    >>> clf.fit(X_train, y_train)
    SVC(random_state=0)
    >>> metrics.plot_roc_curve(clf, X_test, y_test)  # doctest: +SKIP
    >>> plt.show()                                   # doctest: +SKIP
    """
    check_matplotlib_support('plot_roc_curve')

    if response_method not in ("predict_proba", "decision_function", "auto"):
        raise ValueError("response_method must be 'predict_proba', "
                         "'decision_function' or 'auto'")

    if response_method != "auto":
        prediction_method = getattr(estimator, response_method, None)
        if prediction_method is None:
            raise ValueError(
                "response method {} is not defined".format(response_method))
    else:
        predict_proba = getattr(estimator, 'predict_proba', None)
        decision_function = getattr(estimator, 'decision_function', None)
        prediction_method = predict_proba or decision_function

        if prediction_method is None:
            raise ValueError('response methods not defined')

    y_pred = prediction_method(X)
    estimator_name = estimator.__class__.__name__
    estimator_classes = getattr(estimator, "classes_", [])

    if len(estimator_classes) != 2:
        raise ValueError("Estimator {} is not a binary classifier: "
                         "its classes_ attribute is set to: {}".format(
                             estimator_name, estimator_classes))

    if y_pred.ndim != 1:
        if y_pred.shape[1] > 2:
            raise ValueError("Predictions by {}.{} should have shape ({}, 2),"
                             " got {}.".format(
                                 estimator_name, prediction_method.__name__,
                                 y_pred.shape[0], y_pred.shape))
        y_pred = y_pred[:, 1]
    if pos_label is None:
        pos_label = estimator_classes[1]
    fpr, tpr, _ = roc_curve(y, y_pred, pos_label=pos_label,
                            sample_weight=sample_weight,
                            drop_intermediate=drop_intermediate)
    roc_auc = auc(fpr, tpr)
    viz = RocCurveDisplay(fpr, tpr, roc_auc, estimator_name)
    return viz.plot(ax=ax, name=name, **kwargs)
