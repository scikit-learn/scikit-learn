from .. import average_precision_score
from .. import precision_recall_curve

from ...utils import check_matplotlib_support
from ...utils.multiclass import type_of_target
from ...utils.validation import check_is_fitted


class PrecisionRecallDisplay:
    """Precision Recall visualization.

    It is recommend to use `sklearn.metrics.plot_precision_recall_curve` to
    create a visualizer. All parameters are stored as attributes.

    Read more in the :ref:`User Guide <visualizations>`.

    Parameters
    -----------
    precision : ndarray of shape (n_thresholds + 1,)
        Precision values.

    recall : ndarray of shape (n_thresholds + 1,)
        Recall values.

    average_precision : float
        Average precision.

    estimator_name : str
        Name of estimator.

    Attributes
    ----------
    line_ : matplotlib Artist
        Precision recall curve.

    ax_ : matplotlib Axes
        Axes with precision recall curve.

    figure_ : matplotlib Figure
        Figure containing the curve.
    """

    def __init__(self, precision, recall, average_precision, estimator_name):
        self.precision = precision
        self.recall = recall
        self.average_precision = average_precision
        self.estimator_name = estimator_name

    def plot(self, ax=None, label_name=None, **kwargs):
        """Plot visualization.

        Extra keyword arguments will be passed to matplotlib's `plot`.

        Parameters
        ----------
        ax : Matplotlib Axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        label_name : str, default=None
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
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots()

        label_name = self.estimator_name if label_name is None else label_name

        line_kwargs = {
            "label": "{} (AP = {:0.2f})".format(label_name,
                                                self.average_precision),
            "drawstyle": "steps-post"
        }
        line_kwargs.update(**kwargs)

        self.line_, = ax.plot(self.recall, self.precision, **line_kwargs)
        ax.set(xlabel="Recall", ylabel="Precision", ylim=[0.0, 1.05],
               xlim=[0.0, 1.0])
        ax.legend(loc='lower left')

        self.ax_ = ax
        self.figure_ = ax.figure
        return self


def plot_precision_recall_curve(estimator, X, y, pos_label=None,
                                sample_weight=None, response_method="auto",
                                label_name=None, ax=None, **kwargs):
    """Plot Precision Recall Curve for binary classifers.

    Extra keyword arguments will be passed to matplotlib's `plot`.

    Read more in the :ref:`User Guide <visualizations>`.

    Parameters
    ----------
    estimator : estimator instance
        Trained classifier.

    X : {array-like, sparse matrix} of shape (n_samples, n_features)
        Input values.

    y : array-like of shape (n_samples,)
        Binary target values.

    pos_label : int or str, default=None
        The label of the positive class.
        When `pos_label=None`, if y_true is in {-1, 1} or {0, 1},
        `pos_label` is set to 1, otherwise an error will be raised.

    sample_weight : array-like of shape (n_samples,), default=None
        Sample weights.

    response_method : {'predict_proba', 'decision_function', 'auto'}, \
                      default='auto'
        Specifies whether to use :term:`predict_proba` or
        :term:`decision_function` as the target response. If set to 'auto',
        :term:`predict_proba` is tried first and if it does not exist
        :term:`decision_function` is tried next.

    label_name : str, default=None
        Name for labeling curve. If `None`, the name of the
        estimator is used.

    ax : matplotlib axes, default=None
        Axes object to plot on. If `None`, a new figure and axes is created.

    **kwargs : dict
        Keyword arguments to be passed to matplotlib's `plot`.

    Returns
    -------
    display : :class:`~sklearn.metrics.PrecisionRecallDisplay`
        Object that stores computed values.
    """
    check_matplotlib_support("plot_precision_recall_curve")
    check_is_fitted(estimator)

    if response_method not in ("predict_proba", "decision_function", "auto"):
        raise ValueError("response_method must be 'predict_proba', "
                         "'decision_function' or 'auto'")

    type_y = type_of_target(y)
    if type_y != 'binary':
        raise ValueError("{} format is not supported".format(type_y))

    error_msg = "response method {} not defined for estimator {}"
    if response_method != "auto":
        prediction_method = getattr(estimator, response_method, None)
        if prediction_method is None:
            raise ValueError(error_msg.format(response_method,
                                              estimator.__class__.__name__))
        is_predict_proba = response_method == 'predict_proba'
    else:
        predict_proba = getattr(estimator, 'predict_proba', None)
        decision_function = getattr(estimator, 'decision_function', None)
        prediction_method = predict_proba or decision_function
        if prediction_method is None:
            raise ValueError(error_msg.format(
                "decision_function or predict_proba",
                estimator.__class__.__name__))
        is_predict_proba = prediction_method == predict_proba

    y_pred = prediction_method(X)

    if is_predict_proba and y_pred.ndim != 1:
        if y_pred.shape[1] > 2:
            raise ValueError("Estimator should solve a "
                             "binary classification problem")
        y_pred = y_pred[:, 1]

    precision, recall, _ = precision_recall_curve(y, y_pred,
                                                  pos_label=pos_label,
                                                  sample_weight=sample_weight)
    average_precision = average_precision_score(y, y_pred,
                                                sample_weight=sample_weight)
    viz = PrecisionRecallDisplay(precision, recall, average_precision,
                                 estimator.__class__.__name__)
    return viz.plot(ax=ax, label_name=label_name, **kwargs)
