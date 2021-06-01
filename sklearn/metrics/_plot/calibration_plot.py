from ...base import is_classifier
from ...utils import check_matplotlib_support
from .base import _check_classifier_response_method


class CalibrationDisplay:
    """Calibration curve (also known as reliability diagram) visualization.

    It is recommended to use
    :func:`~sklearn.calibration.CalibrationDisplay.from_estimator` or
    :func:`~sklearn.calibration.CalibrationDisplay.from_predictions`
    to create a `CalibrationDisplay`. All parameters are stored as attributes.

    Read more about calibration in the :ref:`User Guide <calibration>` and
    more about the scikit-learn visualization API in :ref:`visualizations`.

    .. versionadded:: 1.0

    Parameters
    -----------
    prob_true : ndarray
        The proportion of samples whose class is the positive class (fraction
        of positives), in each bin.

    prob_pred : ndarray
        The mean predicted probability in each bin.

    y_prob : ndarray of shape (n_samples,)
        Probability estimates for the positive class, for each sample.

    name : str, default=None
        Name for labeling curve.

    Attributes
    ----------
    line_ : matplotlib Artist
        Calibration curve.

    ax_ : matplotlib Axes
        Axes with calibration curve.

    figure_ : matplotlib Figure
        Figure containing the curve.

    Examples
    --------
    >>> from sklearn.datasets import make_classification
    >>> from sklearn.model_selection import train_test_split
    >>> from sklearn.linear_model import LogisticRegression
    >>> from sklearn.calibration import calibration_curve
    >>> from sklearn.metrics import CalibrationDisplay
    >>> X, y = make_classification(random_state=0)
    >>> X_train, X_test, y_train, y_test = train_test_split(
    ...     X, y, random_state=0)
    >>> clf = LogisticRegression(random_state=0)
    >>> clf.fit(X_train, y_train)
    LogisticRegression(random_state=0)
    >>> y_prob = clf.predict_proba(X_test)[:, 1]
    >>> prob_true, prob_pred = calibration_curve(y_test, y_prob, n_bins=10)
    >>> disp = CalibrationDisplay(prob_true, prob_pred, y_prob)
    >>> disp.plot() # doctest: +SKIP
    """
    def __init__(self, prob_true, prob_pred, y_prob, *, name=None):
        self.prob_true = prob_true
        self.prob_pred = prob_pred
        self.y_prob = y_prob
        self.name = name

    def plot(self, *, ax=None, name=None, ref_line=True, **kwargs):
        """Plot visualization.

        Extra keyword arguments will be passed to
        :func:`matplotlib.pyplot.plot`.

        Parameters
        ----------
        ax : Matplotlib Axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        name : str, default=None
            Name for labeling curve.

        ref_line : bool, default=True
            If `True`, plots a reference line representing a perfectly
            calibrated classifier.

        **kwargs : dict
            Keyword arguments to be passed to :func:`matplotlib.pyplot.plot`.

        Returns
        -------
        display : :class:`~sklearn.calibration.CalibrationDisplay`
            Object that stores computed values.
        """
        check_matplotlib_support("CalibrationDisplay.plot")
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots()

        name = self.name if name is None else name
        self.name = name

        line_kwargs = {}
        if name is not None:
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

    @classmethod
    def from_estimator(cls, estimator, X, y, *,
                       n_bins=5, strategy='uniform', name=None,
                       ref_line=True, ax=None, **kwargs):
        """Plot calibration curve, also known as reliability diagrams, for
        binary classifiers, using an estimator and data.

        The average predicted probability for each bin is plotted on the x-axis
        and the fraction of positive classes in each bin is plotted on the
        y-axis.

        Extra keyword arguments will be passed to
        :func:`matplotlib.pyplot.plot`.

        Read more about calibration in the :ref:`User Guide <calibration>` and
        more about the scikit-learn visualization API in :ref:`visualizations`.

        .. versionadded:: 1.0

        Parameters
        ----------
        estimator : estimator instance
            Fitted classifier or a fitted :class:`~sklearn.pipeline.Pipeline`
            in which the last estimator is a classifier. The classifier must
            have a :term:`predict_proba` method.

        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Input values.

        y : array-like of shape (n_samples,)
            Binary target values.

        n_bins : int, default=5
            Number of bins to discretize the [0, 1] interval into when
            calculating the calibration curve. A bigger number requires more
            data.

        strategy : {'uniform', 'quantile'}, default='uniform'
            Strategy used to define the widths of the bins.

            - `'uniform'`: The bins have identical widths.
            - `'quantile'`: The bins have the same number of samples and depend
              on predicted probabilities.

        name : str, default=None
            Name for labeling curve. If `None`, the name of the estimator is
            used.

        ref_line : bool, default=True
            If `True`, plots a reference line representing a perfectly
            calibrated classifier.

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        **kwargs : dict
            Keyword arguments to be passed to :func:`matplotlib.pyplot.plot`.

        Returns
        -------
        display : :class:`~sklearn.metrics.CalibrationDisplay`.
            Object that stores computed values.

        Examples
        --------
        >>> import matplotlib.pyplot as plt
        >>> from sklearn.datasets import make_classification
        >>> from sklearn.model_selection import train_test_split
        >>> from sklearn.linear_model import LogisticRegression
        >>> from sklearn.metrics import CalibrationDisplay
        >>> X, y = make_classification(random_state=0)
        >>> X_train, X_test, y_train, y_test = train_test_split(
        ...     X, y, random_state=0)
        >>> clf = LogisticRegression(random_state=0)
        >>> clf.fit(X_train, y_train)
        LogisticRegression(random_state=0)
        >>> disp = CalibrationDisplay.from_estimator(clf, X_test, y_test)
        >>> plt.show()
        """
        method_name = f"{cls.__name__}.from_estimator"
        check_matplotlib_support(method_name)

        if not is_classifier(estimator):
            raise ValueError("'estimator' should be a fitted classifier.")

        prediction_method = _check_classifier_response_method(
            estimator, response_method='predict_proba'
        )
        y_prob = prediction_method(X)

        binary_error = "Only binary classification is supported."
        if not len(estimator.classes_) == 2:
            raise ValueError(binary_error)
        if y_prob.ndim == 1:
            raise ValueError("'estimator.predict_proba' needs to return a 2d "
                             "array.")
        else:
            if y_prob.shape[1] != 2:
                raise ValueError(binary_error)
            else:
                y_prob = y_prob[:, 1]

        name = name if name is not None else estimator.__class__.__name__
        return cls.from_predictions(
            y, y_prob, n_bins=n_bins, strategy=strategy, name=name,
            ref_line=ref_line, ax=ax, **kwargs
        )

    @classmethod
    def from_predictions(cls, y_true, y_prob, *,
                         n_bins=5, strategy='uniform', name=None,
                         ref_line=True, ax=None, **kwargs):
        """Plot calibration curve, also known as reliability diagrams, for
        binary classifiers, using true and predicted labels.

        The average predicted probability for each bin is plotted on the x-axis
        and the fraction of positive classes in each bin is plotted on the
        y-axis.

        Extra keyword arguments will be passed to
        :func:`matplotlib.pyplot.plot`.

        Read more about calibration in the :ref:`User Guide <calibration>` and
        more about the scikit-learn visualization API in :ref:`visualizations`.

        .. versionadded:: 1.0

        Parameters
        ----------
        y_true : array-like of shape (n_samples,)
            True labels.

        y_prob : array-like of shape (n_samples,)
            The predicted probabilities of the positive class.

        n_bins : int, default=5
            Number of bins to discretize the [0, 1] interval into when
            calculating the calibration curve. A bigger number requires more
            data.

        strategy : {'uniform', 'quantile'}, default='uniform'
            Strategy used to define the widths of the bins.

            - `'uniform'`: The bins have identical widths.
            - `'quantile'`: The bins have the same number of samples and depend
              on predicted probabilities.

        name : str, default=None
            Name for labeling curve.

        ref_line : bool, default=True
            If `True`, plots a reference line representing a perfectly
            calibrated classifier.

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        **kwargs : dict
            Keyword arguments to be passed to :func:`matplotlib.pyplot.plot`.

        Returns
        -------
        display : :class:`~sklearn.metrics.CalibrationDisplay`.
            Object that stores computed values.

        Examples
        --------
        >>> import matplotlib.pyplot as plt
        >>> from sklearn.datasets import make_classification
        >>> from sklearn.model_selection import train_test_split
        >>> from sklearn.linear_model import LogisticRegression
        >>> from sklearn.metrics import CalibrationDisplay
        >>> X, y = make_classification(random_state=0)
        >>> X_train, X_test, y_train, y_test = train_test_split(
        ...     X, y, random_state=0)
        >>> clf = LogisticRegression(random_state=0)
        >>> clf.fit(X_train, y_train)
        LogisticRegression(random_state=0)
        >>> y_prob = clf.predict_proba(X_test)[:, 1]
        >>> disp = CalibrationDisplay.from_predictions(y_test, y_prob)
        >>> plt.show()
        """
        method_name = f"{cls.__name__}.from_estimator"
        check_matplotlib_support(method_name)
        from ...calibration import calibration_curve
        prob_true, prob_pred = calibration_curve(
            y_true, y_prob, n_bins=n_bins, strategy=strategy
        )

        disp = cls(prob_true=prob_true, prob_pred=prob_pred, y_prob=y_prob,
                   name=name)
        return disp.plot(ax=ax, ref_line=ref_line, **kwargs)
