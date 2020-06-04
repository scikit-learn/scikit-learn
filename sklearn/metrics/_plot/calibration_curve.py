from .base import _check_classifer_response_method

from ...utils import check_matplotlib_support
from ...utils.validation import _deprecate_positional_args
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
        The proportion of samples whose class is the positive class, in each
        bin (fraction of positives).

    prob_pred : ndarray
        The mean predicted probability in each bin.

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
    def __init__(self, prob_true, prob_pred, *,
                 estimator_name=None):
        self.prob_true = prob_true
        self.prob_pred = prob_pred
        self.estimator_name = estimator_name

    @_deprecate_positional_args
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
            Keyword arguments to be passed to matplotlib's `plot`.

        Returns
        -------
        display : :class:`~sklearn.metrics.CalibrationDisplay`
            Object that stores computed values.
        """
        check_matplotlib_support("CalifrationDisplay.plot")
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots()

        name = self.estimator_name if name is None else name

        line_kwargs = {}
        if name is not None:
            line_kwargs["label"] = name
        line_kwargs.update(**kwargs)

        existing_ref_line = ('Perfectly calibrated' in
                             ax.get_legend_handles_labels()[1])
        if ref_line and not existing_ref_line:
            ax.plot([0, 1], [0, 1], "k:", label="Perfectly calibrated")
        self.line_ = ax.plot(self.prob_pred, self.prob_true, "s-",
                             **line_kwargs)

        if "label" in line_kwargs:
            ax.legend(loc="lower right")

        ax.set(xlabel="Mean predicted probability",
               ylabel="Fraction of positives")

        self.axs_ = ax
        self.figure_ = ax.figure
        return self


@_deprecate_positional_args
def plot_calibration_curve(estimator, X, y, *,
                           normalize=False, n_bins=5, strategy='uniform',
                           name=None, ref_line=True, ax=None, **kwargs):
    """Plot calibration curve for binary classifiers.

    Extra keyword arguments will be passed to matplotlib's `plot`.

    Read more in the :ref:`User Guide <calibration>`.

    Parameters
    ----------
    estimator : estimator instance
        Fitted classifier or a fitted :class:`~sklearn.pipeline.Pipeline`
        in which the last estimator is a classifier. The classifier must
        have a `predict_proba` method; set `probability=True` for
        :class:`~sklearn.svm.SVC` and :class:`~sklearn.svm.NuSVC`
        (see: :ref:`User Guide <scores_probabilities>`) or use
        :class:`~sklearn.calibration.CalibratedClassifierCV`.

    X : {array-like, sparse matrix} of shape (n_samples, n_features)
        Input values.

    y : array-like of shape (n_samples,)
        Binary target values.

    normalize : bool, default=False
        Whether the output of `estimator.predict_proba(X)` needs to be
        normalized into the [0, 1] interval (i.e. is not a proper probability).
        If `True`, the smallest value is linearly mapped onto 0 and the largest
        one onto 1.

    n_bins : int, default=5
        Number of bins to discretize the [0, 1] interval into when calculating
        the calibration curve.

    strategy : {‘uniform’, ‘quantile’}, default=’uniform’
        Strategy used to define the widths of the bins.

        uniform
            The bins have identical widths.
        quantile
            The bins have the same number of samples and depend on `y_prob`.

    name : str, default=None
        Name for labeling curve. If `None`, the name of the estimator is used.

    ref_line : bool, default=True
        If `True`, plots a reference line representing a perfectly calibrated
        classifier.

    ax : matplotlib axes, default=None
        Axes object to plot on. If `None`, a new figure and axes is created.

    **kwargs : dict
        Keyword arguments to be passed to :func:`matplotlib.pyplot.plot`.

    Returns
    -------
    display : :class:`~sklearn.metrics.CalibrationDisplay`
        Object that stores computed values.
    """
    check_matplotlib_support("plot_calibration_curve")

    classification_error = ("{} should be a binary classifier".format(
        estimator.__class__.__name__))
    if not is_classifier(estimator):
        raise ValueError(classification_error)

    prediction_method = _check_classifer_response_method(
        estimator, response_method='predict_proba'
    )

    y_prob = prediction_method(X)

    if y_prob.ndim != 1:
        if y_prob.shape[1] != 2:
            raise ValueError(classification_error)
        else:
            y_prob = y_prob[:, 1]

    prob_true, prob_pred = calibration_curve(
        y, y_prob, normalize=normalize, n_bins=n_bins, strategy=strategy
    )
    name = name if name is not None else estimator.__class__.__name__
    viz = CalibrationDisplay(
        prob_true=prob_true, prob_pred=prob_pred, estimator_name=name
    )
    return viz.plot(ax=ax, name=name, ref_line=ref_line, **kwargs)
