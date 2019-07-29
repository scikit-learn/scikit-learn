from .. import auc
from .. import roc_curve

from ...utils import check_matplotlib_support


class RocCurveVisualizer:
    """ROC Curve visualization.

    Parameters
    ----------
    fpr : ndarray
        False positive rate.
    tpr : ndarray
        True positive rate.
    auc_roc : float
        Area under ROC curve.
    estimator_name : str
        Name of estimator.

    Attributes
    ----------
    line_ : matplotlib Artist
        ROC Curve.
    ax_ : matplotlib Axes
        Axes with ROC Curve
    figure_ : matplotlib Figure
        Figure containing the curve
    """

    def __init__(self, fpr, tpr, auc_roc, estimator_name):
        self.fpr = fpr
        self.tpr = tpr
        self.auc_roc = auc_roc
        self.estimator_name = estimator_name

    def plot(self, ax=None, name=None, **kwargs):
        """Plot visualization

        Extra keyword arguments will be passed to matplotlib's ``plot``.

        Parameters
        ----------
        ax : Matplotlib Axes or None, default=None
            Axes object to plot on.

        name : str or None, default=None
            Name of ROC Curve for labeling. If `None`, use the name of the
            estimator.
        """
        check_matplotlib_support('plot_roc_curve')
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots()

        name = self.estimator_name if name is None else name

        if 'label' not in kwargs:
            label = "{} (AUC = {:0.2f})".format(name, self.auc_roc)
            kwargs['label'] = label
        self.line_ = ax.plot(self.fpr, self.tpr, **kwargs)[0]
        ax.set_xlabel("False Positive Rate")
        ax.set_ylabel("True Positive Rate")
        ax.legend()

        self.ax_ = ax
        self.figure_ = ax.figure
        return self


def plot_roc_curve(estimator, X, y, pos_label=None, sample_weight=None,
                   drop_intermediate=True, response_method="auto",
                   name=None, ax=None, **kwargs):
    """Plot Receiver operating characteristic (ROC) curve

    Extra keyword arguments will be passed to matplotlib's `plot`.

    Parameters
    ----------
    estimator : estimator instance
        Trained classifier.

    X : {array-like, sparse matrix}, shape (n_samples, n_features)
        Input values.

    y : array-like, shape (n_samples, )
        Target values.

    pos_label : int or str, default=None
        The label of the positive class.
        When `pos_label=None`, if y_true is in {-1, 1} or {0, 1},
        `pos_label` is set to 1, otherwise an error will be raised.

    sample_weight : array-like, shape (n_samples, ) or None, default=None
        Sample weights.

    drop_intermediate : boolean, default=True
        Whether to drop some suboptimal thresholds which would not appear
        on a plotted ROC curve. This is useful in order to create lighter
        ROC curves.

    response_method : 'predict_proba', 'decision_function', or 'auto' \
    default='auto'
        Specifies whether to use `predict_proba` or `decision_function` as the
        target response. If set to 'auto', `predict_proba` is tried first
        and if it does not exist `decision_function` is tried next.

    name : str or None, default=None
        Name of ROC Curve for labeling. If `None`, use the name of the
        estimator.

    ax : matplotlib axes, default=None
        axes object to plot on

    Returns
    -------
    viz : :class:`sklearn.metrics.plot.RocCurveVisualizer`
        object that stores computed values
    """
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

    if y_pred.ndim != 1:
        if y_pred.shape[1] > 2:
            raise ValueError("Estimator must be a binary classifier")
        y_pred = y_pred[:, 1]
    fpr, tpr, _ = roc_curve(y, y_pred, pos_label=pos_label,
                            drop_intermediate=drop_intermediate)
    auc_roc = auc(fpr, tpr)
    viz = RocCurveVisualizer(fpr, tpr, auc_roc, estimator.__class__.__name__)
    return viz.plot(ax=ax, name=name, **kwargs)
