import numpy as np

from sklearn.base import is_classifier
from sklearn.preprocessing import label_binarize
from sklearn.utils.multiclass import type_of_target


def _check_classifier_response_method(estimator, response_method):
    """Return prediction method from the response_method

    Parameters
    ----------
    estimator: object
        Classifier to check

    response_method: {'auto', 'predict_proba', 'decision_function'}
        Specifies whether to use :term:`predict_proba` or
        :term:`decision_function` as the target response. If set to 'auto',
        :term:`predict_proba` is tried first and if it does not exist
        :term:`decision_function` is tried next.

    Returns
    -------
    prediction_method: callable
        prediction method of estimator
    """

    if response_method not in ("predict_proba", "decision_function", "auto"):
        raise ValueError("response_method must be 'predict_proba', "
                         "'decision_function' or 'auto'")

    error_msg = "response method {} is not defined in {}"
    if response_method != "auto":
        prediction_method = getattr(estimator, response_method, None)
        if prediction_method is None:
            raise ValueError(error_msg.format(response_method,
                                              estimator.__class__.__name__))
    else:
        predict_proba = getattr(estimator, 'predict_proba', None)
        decision_function = getattr(estimator, 'decision_function', None)
        prediction_method = predict_proba or decision_function
        if prediction_method is None:
            raise ValueError(error_msg.format(
                "decision_function or predict_proba",
                estimator.__class__.__name__))

    return prediction_method


def _get_response(X, estimator, response_method,
                  n_classes=2, pos_label=None):
    """Return response and positive label.

    Parameters
    ----------
    X : {array-like, sparse matrix} of shape (n_samples, n_features)
        Input values.

    estimator : estimator instance
        Fitted classifier or a fitted :class:`~sklearn.pipeline.Pipeline`
        in which the last estimator is a classifier.

    response_method: {'auto', 'predict_proba', 'decision_function'}
        Specifies whether to use :term:`predict_proba` or
        :term:`decision_function` as the target response. If set to 'auto',
        :term:`predict_proba` is tried first and if it does not exist
        :term:`decision_function` is tried next.

    pos_label : str or int, default=None
        The class considered as the positive class when computing
        the metrics. By default, `estimators.classes_[1]` is
        considered as the positive class.

    Returns
    -------
    y_pred: ndarray of shape (n_samples,)
        Target scores calculated from the provided response_method
        and pos_label.

    pos_label: str or int
        The class considered as the positive class when computing
        the metrics.
    """
    classification_error = (
        "{} should be a classifier".format(estimator.__class__.__name__)
    )

    binary_classification_error = (
        "{} should be a binary classifier".format(estimator.__class__.__name__)
    )

    if not is_classifier(estimator):
        raise ValueError(classification_error)

    prediction_method = _check_classifier_response_method(
        estimator, response_method)

    y_pred = prediction_method(X)

    if n_classes > 2:
        return y_pred, None

    if pos_label is not None and pos_label not in estimator.classes_:
        raise ValueError(
            f"The class provided by 'pos_label' is unknown. Got "
            f"{pos_label} instead of one of {estimator.classes_}"
        )

    if y_pred.ndim != 1:  # `predict_proba`
        if y_pred.shape[1] != 2:
            raise ValueError(binary_classification_error)
        if pos_label is None:
            pos_label = estimator.classes_[1]
            y_pred = y_pred[:, 1]
        else:
            class_idx = np.flatnonzero(estimator.classes_ == pos_label)
            y_pred = y_pred[:, class_idx]
    else:
        if pos_label is None:
            pos_label = estimator.classes_[1]
        elif pos_label == estimator.classes_[0]:
            y_pred *= -1

    return y_pred, pos_label


def _plot_curve(plot_curve_func,
                estimator, X, y, *,
                response_method="auto",
                name=None, ax=None, pos_label=None, **kwargs):
    """Plot curve base function.

    Extra keyword arguments will be passed to matplotlib's `plot`.

    Read more in the :ref:`User Guide <visualizations>`.

    Parameters
    ----------
    plot_curve_func : callable
        Plot curve function (computes the plot metrics and returns a display).
        Can be either `_get_roc_curve_display` or
        `_get_precision_recall_display`.

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
        Name of the curve for labeling. If `None`, use the name of the
        estimator.

    ax : Matplotlib axes or array-like of Matplotlib axes, default=None
        Axes object to plot on. If `None`, a new figure and axes is created.
        For the multiclass cenario:
        - If a single axis is passed in, all plots are plotted in
          the same axis.
        - If an array-like of axes are passed in, the roc curve
          plots will be drawn directly into these axes.

    pos_label : str or int, default=None
        The class considered as the positive class when computing the roc auc
        metrics. By default, `estimators.classes_[1]` is considered
        as the positive class.

        .. versionadded:: 0.24

    Returns
    -------
    display : :class:`~sklearn.metrics.RocCurveDisplay`
              or :class:`~sklearn.metrics.PrecisionRecallDisplay`
        Object or array-like of object that stores computed values.

    """
    import matplotlib.pyplot as plt

    n_classes = len(np.unique(y)) if y.ndim == 1 else y.shape[1]

    y_type = type_of_target(y)

    y_pred, pos_label = _get_response(
        X, estimator, response_method,
        n_classes=n_classes, pos_label=pos_label)

    name = estimator.__class__.__name__ if name is None else name

    # Early exit if the axes does not have the correct number of axes
    if ax is not None and not isinstance(ax, plt.Axes):
        axes = np.asarray(ax, dtype=object)
        if axes.size != n_classes:
            raise ValueError("Expected ax to have {} axes, got {}".format(
                n_classes, axes.size))

    if n_classes == 2:

        viz = plot_curve_func(y, y_pred, pos_label=pos_label,
                              y_type=y_type, name=name)

        return viz.plot(ax=ax, name=name, **kwargs)
    else:
        # binarize if y is a vector
        if y.ndim == 1:
            y = label_binarize(y, classes=np.unique(y))

        if ax is None:
            fig, ax = plt.subplots(1, n_classes)

        vizs = []

        for i in range(n_classes):
            viz = plot_curve_func(y[:, i], y_pred[:, i],
                                  y_type=y_type, name=name)

            axes = ax if isinstance(ax, plt.Axes) else ax[i]

            viz.plot(ax=axes, name='{} for class {}'.format(name, i), **kwargs)

            vizs.append(viz)

        return vizs
