from .. import roc_curve


class RocCurveViz:
    """ROC Curve visualization

    Parameters
    ----------

    estimator : estimator instance
        Trained classifier.

    X : {array-like, sparse matrix}, shape (n_samples, n_features)
        Input values.

    y : array-like, shape (n_samples,)
        Target values.

    pos_label : int or str, default=None
        The label of the positive class.
        When ``pos_label=None``, if y_true is in {-1, 1} or {0, 1},
        ``pos_label`` is set to 1, otherwise an error will be raised.

    sample_weight : array-like of shape = [n_samples], optional
        Sample weights.

    drop_intermediate : boolean, optional (default=True)
        Whether to drop some suboptimal thresholds which would not appear
        on a plotted ROC curve. This is useful in order to create lighter
        ROC curves.

    Attributes
    ----------
    hello_ : int
    """

    def __init__(self, estimator, X, y, *,
                 pos_label=None,
                 sample_weight=None,
                 drop_intermediate=True,
                 response_method="predict_proba"):
        """Computes and stores values needed for visualization"""

        if y_pred.ndim == 2:
            y_pred = y_pred[:,1]

    def plot(self, ax=None):
        """Plot visualization

        Parameters
        ----------
        ax : Matplotlib Axes, optional (default=None)
            axes object to plot on
        """
        return


def plot_roc_curve(estimator,
                   X,
                   y,
                   pos_label=None,
                   sample_weight=None,
                   drop_intermediate=True,
                   response_method="predict_proba",
                   ax=None):
    """Plot Receiver operating characteristic (ROC) curve

    Note: this implementation is restricted to the binary classification task.

    Read more in the :ref:`User Guide <roc_metrics>`.

    Parameters
    ----------

    estimator : estimator instance
        Trained classifier.

    X : {array-like, sparse matrix}, shape (n_samples, n_features)
        Input values.

    y : array-like, shape (n_samples,)
        Target values.

    pos_label : int or str, default=None
        The label of the positive class.
        When ``pos_label=None``, if y_true is in {-1, 1} or {0, 1},
        ``pos_label`` is set to 1, otherwise an error will be raised.

    sample_weight : array-like of shape = [n_samples], optional (default=None)
        Sample weights.

    drop_intermediate : boolean, optional (default=True)
        Whether to drop some suboptimal thresholds which would not appear
        on a plotted ROC curve. This is useful in order to create lighter
        ROC curves.

    response_method : 'predict_proba' or 'decision_function' optional \
    (default='predict_proba')
        Method to call estimator to get target scores

    ax : matplotlib axes, optional (default=None)
        axes object to plot on

    Returns
    -------
    viz : :class:`sklearn.metrics.plot.RocCurveViz`
        object that stores computed values
    """
    viz = RocCurveViz(estimator,
                      X,
                      y,
                      sample_weight=sample_weight,
                      pos_label=pos_label,
                      drop_intermediate=drop_intermediate,
                      response_method=response_method)
    viz.plot(ax=ax)
    return viz
