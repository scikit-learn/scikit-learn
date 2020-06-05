from .. import brier_score_loss
from ...utils import check_matplotlib_support
from ...base import is_classifier
from ...calibration import calibration_curve


class CalibrationDisplay:
    """Calibration visualization.

    It is recommend to use :func:`~sklearn.metrics.plot_calibration_curve`
    to create a visualizer. All parameters are stored as attributes.

    Read more in the :ref:`User Guide <calibration>`.

    Parameters
    -----------
    prob_true : ndarray
        The proportion of samples whose class is the positive class (fraction
        of positives), in each bin.

    prob_pred : ndarray
        The mean predicted probability in each bin.

    y_prob : ndarray of shape (n_samples,)
        Probability estimates for the positive class.

    brier_value : int or None
        The Brier score value. If None, the Brier score is not shown.

    estimator_name : str, default=None
        Name of estimator. If None, then the estimator name is not shown.

    Attributes
    ----------
    line_ : matplotlib Artist
        Calibration curve.

    ax_ : matplotlib Axes
        Axes with calibration curve.

    figure_ : matplotlib Figure
        Figure containing the curve.
    """
    def __init__(self, prob_true, prob_pred, y_prob, *,
                 brier_value=None, estimator_name=None):
        self.prob_true = prob_true
        self.prob_pred = prob_pred
        self.y_prob = y_prob
        self.brier_value = brier_value
        self.estimator_name = estimator_name

    def plot(self, ax=None, *, name=None, ref_line=True, **kwargs):
        """Plot visualization.

        Extra keyword arguments will be passed to matplotlib's `plot`.

        Parameters
        ----------
        ax : Matplotlib Axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        name : str, default=None
            Name of calibration curve for labeling. If `None`, use the
            name of the estimator.

        ref_line : bool, default=True
            If `True`, plots a reference line representing a perfectly
            calibrated classifier.

        **kwargs : dict
            Keyword arguments to be passed to :func:`matplotlib.pyplot.plot`.

        Returns
        -------
        display : :class:`~sklearn.metrics.CalibrationDisplay`
            Object that stores computed values.
        """
        check_matplotlib_support("CalibrationDisplay.plot")
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots()

        name = self.estimator_name if name is None else name

        line_kwargs = {}
        if self.brier_value is not None and name is not None:
            line_kwargs["label"] = \
                f"{name} (Brier: {self.brier_value:0.3f})"
        elif self.brier_value is not None:
            line_kwargs["label"] = f"Brier: {self.brier_value:0.3f}"
        elif name is not None:
            line_kwargs["label"] = name
        line_kwargs.update(**kwargs)

        existing_ref_line = ('Perfectly calibrated' in
                             ax.get_legend_handles_labels()[1])
        if ref_line and not existing_ref_line:
            ax.plot([0, 1], [0, 1], "k:", label="Perfectly calibrated")
        self.line_ = ax.plot(self.prob_pred, self.prob_true, "s-",
                             **line_kwargs)[0]

        if "label" in line_kwargs:
            ax.legend(loc="lower right")

        ax.set(xlabel="Mean predicted probability",
               ylabel="Fraction of positives")

        self.ax_ = ax
        self.figure_ = ax.figure
        return self


def plot_calibration_curve(estimator, X, y, *,
                           n_bins=5, strategy='uniform',
                           name=None, ref_line=True, brier_score=True,
                           ax=None, **kwargs):
    """Plot calibration curve for binary classifiers.

    Extra keyword arguments will be passed to matplotlib's `plot`.

    Read more in the :ref:`User Guide <calibration>`.

    Parameters
    ----------
    estimator : estimator instance
        Fitted classifier or a fitted :class:`~sklearn.pipeline.Pipeline` in
        which the last estimator is a classifier.
        To calculate probability estimates :term:`predict_proba` will be used
        in priority. Otherwise, if no :term:`predict_proba` method exists,
        :term:`decision_function` will be used and the output confidence
        scores will be min-max scaled to the range [0,1].

    X : {array-like, sparse matrix} of shape (n_samples, n_features)
        Input values.

    y : array-like of shape (n_samples,)
        Binary target values.

    n_bins : int, default=5
        Number of bins to discretize the [0, 1] interval into when calculating
        the calibration curve.

    strategy : {'uniform', 'quantile'}, default='uniform'
        Strategy used to define the widths of the bins.

        `'uniform'`: The bins have identical widths.
        `'quantile'`: The bins have the same number of samples and depend on
        `estimator.predict_proba(X)`.

    name : str, default=None
        Name for labeling curve. If `None`, the name of the estimator is used.

    brier_score: bool, default=True
        If `True`, include Brier score in legend.

    ref_line : bool, default=True
        If `True`, plots a reference line representing a perfectly calibrated
        classifier.

    ax : matplotlib axes, default=None
        Axes object to plot on. If `None`, a new figure and axes is created.

    **kwargs : dict
        Keyword arguments to be passed to :func:`matplotlib.pyplot.plot`.

    Returns
    -------
    display : :class:`~sklearn.metrics.CalibrationDisplay`.
        Object that stores computed values.
    """
    check_matplotlib_support("plot_calibration_curve")
    binary_error = "Only binary classification is supported."

    if not is_classifier(estimator):
        raise ValueError("The estimator parameter should be a fitted binary "
                         "classifier")

    predict_proba = getattr(estimator, 'predict_proba', None)
    decision_function = getattr(estimator, 'decision_function', None)
    prediction_method = predict_proba or decision_function
    if prediction_method is None:
        msg = ("Neither response method 'predict_proba' nor "
               "'decision_function' are defined in "
               f"{estimator.__class__.__name__}")
        raise ValueError(msg)

    y_prob = prediction_method(X)

    if predict_proba is None:
        y_prob = (y_prob - y_prob.min()) / (y_prob.max() - y_prob.min())
    else:
        if not len(estimator.classes_) == 2:
            raise ValueError(binary_error)
        if y_prob.ndim != 1:
            if y_prob.shape[1] != 2:
                raise ValueError(binary_error)
            else:
                y_prob = y_prob[:, 1]

    prob_true, prob_pred = calibration_curve(
        y, y_prob, n_bins=n_bins, strategy=strategy
    )
    if brier_score:
        pos_label = estimator.classes_[1]
        brier_value = brier_score_loss(y, y_prob, pos_label=pos_label)
    else:
        brier_value = None
    name = name if name is not None else estimator.__class__.__name__
    viz = CalibrationDisplay(
        prob_true=prob_true, prob_pred=prob_pred, y_prob=y_prob,
        brier_value=brier_value, estimator_name=name
    )
    return viz.plot(ax=ax, name=name, ref_line=ref_line, **kwargs)
