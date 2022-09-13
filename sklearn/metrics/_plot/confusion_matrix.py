import numpy as np

from .. import confusion_matrix
from ...utils import check_matplotlib_support
from ...utils import deprecated
from ...utils.multiclass import unique_labels
from ...utils._plot import plot_heatmap
from ...base import is_classifier


class ConfusionMatrixDisplay:
    """Confusion Matrix visualization.

    It is recommend to use
    :func:`~sklearn.metrics.ConfusionMatrixDisplay.from_estimator` or
    :func:`~sklearn.metrics.ConfusionMatrixDisplay.from_predictions` to
    create a :class:`ConfusionMatrixDisplay`. All parameters are stored as
    attributes.

    Read more in the :ref:`User Guide <visualizations>`.

    Parameters
    ----------
    confusion_matrix : ndarray of shape (n_classes, n_classes)
        Confusion matrix.

    display_labels : ndarray of shape (n_classes,), default=None
        Display labels for plot. If None, display labels are set from 0 to
        `n_classes - 1`.

    Attributes
    ----------
    im_ : matplotlib AxesImage
        Image representing the confusion matrix.

    text_ : ndarray of shape (n_classes, n_classes), dtype=matplotlib Text, \
            or None
        Array of matplotlib axes. `None` if `include_values` is false.

    ax_ : matplotlib Axes
        Axes with confusion matrix.

    figure_ : matplotlib Figure
        Figure containing the confusion matrix.

    See Also
    --------
    confusion_matrix : Compute Confusion Matrix to evaluate the accuracy of a
        classification.
    ConfusionMatrixDisplay.from_estimator : Plot the confusion matrix
        given an estimator, the data, and the label.
    ConfusionMatrixDisplay.from_predictions : Plot the confusion matrix
        given the true and predicted labels.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from sklearn.datasets import make_classification
    >>> from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
    >>> from sklearn.model_selection import train_test_split
    >>> from sklearn.svm import SVC
    >>> X, y = make_classification(random_state=0)
    >>> X_train, X_test, y_train, y_test = train_test_split(X, y,
    ...                                                     random_state=0)
    >>> clf = SVC(random_state=0)
    >>> clf.fit(X_train, y_train)
    SVC(random_state=0)
    >>> predictions = clf.predict(X_test)
    >>> cm = confusion_matrix(y_test, predictions, labels=clf.classes_)
    >>> disp = ConfusionMatrixDisplay(confusion_matrix=cm,
    ...                               display_labels=clf.classes_)
    >>> disp.plot()
    <...>
    >>> plt.show()
    """

    def __init__(self, confusion_matrix, *, display_labels=None):
        self.confusion_matrix = confusion_matrix
        self.display_labels = display_labels

    def plot(
        self,
        *,
        include_values=True,
        cmap="viridis",
        xticks_rotation="horizontal",
        values_format=None,
        ax=None,
        colorbar=True,
        im_kw=None,
        text_kw=None,
    ):
        """Plot visualization.

        Parameters
        ----------
        include_values : bool, default=True
            Includes values in confusion matrix.

        cmap : str or matplotlib Colormap, default='viridis'
            Colormap recognized by matplotlib.

        xticks_rotation : {'vertical', 'horizontal'} or float, \
                         default='horizontal'
            Rotation of xtick labels.

        values_format : str, default=None
            Format specification for values in confusion matrix. If `None`,
            the format specification is 'd' or '.2g' whichever is shorter.

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        colorbar : bool, default=True
            Whether or not to add a colorbar to the plot.

        im_kw : dict, default=None
            Dict with keywords passed to `matplotlib.pyplot.imshow` call.

        text_kw : dict, default=None
            Dict with keywords passed to `matplotlib.pyplot.text` call.

            .. versionadded:: 1.2

        Returns
        -------
        display : :class:`~sklearn.metrics.ConfusionMatrixDisplay`
        """
        method_name = f"{self.__class__.__name__}.plot"
        check_matplotlib_support(method_name)
        cm = self.confusion_matrix
        n_classes = cm.shape[0]
        if self.display_labels is None:
            display_labels = np.arange(n_classes)
        else:
            display_labels = self.display_labels
        fig, ax, im, text = plot_heatmap(
            cm,
            ylabel="True label",
            xlabel="Predicted label",
            yticklabels=display_labels,
            xticklabels=display_labels,
            xticks_rotation=xticks_rotation,
            ax=ax,
            cmap=cmap,
            include_values=include_values,
            values_format=values_format,
            colorbar=colorbar,
            im_kw=im_kw,
            text_kw=text_kw,
        )

        ax.set_ylim((n_classes - 0.5, -0.5))

        self.figure_ = fig
        self.ax_ = ax
        self.im_ = im
        self.text_ = text
        return self

    @classmethod
    def from_estimator(
        cls,
        estimator,
        X,
        y,
        *,
        labels=None,
        sample_weight=None,
        normalize=None,
        display_labels=None,
        include_values=True,
        xticks_rotation="horizontal",
        values_format=None,
        cmap="viridis",
        ax=None,
        colorbar=True,
        im_kw=None,
        text_kw=None,
    ):
        """Plot Confusion Matrix given an estimator and some data.

        Read more in the :ref:`User Guide <confusion_matrix>`.

        .. versionadded:: 1.0

        Parameters
        ----------
        estimator : estimator instance
            Fitted classifier or a fitted :class:`~sklearn.pipeline.Pipeline`
            in which the last estimator is a classifier.

        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Input values.

        y : array-like of shape (n_samples,)
            Target values.

        labels : array-like of shape (n_classes,), default=None
            List of labels to index the confusion matrix. This may be used to
            reorder or select a subset of labels. If `None` is given, those
            that appear at least once in `y_true` or `y_pred` are used in
            sorted order.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights.

        normalize : {'true', 'pred', 'all'}, default=None
            Either to normalize the counts display in the matrix:

            - if `'true'`, the confusion matrix is normalized over the true
              conditions (e.g. rows);
            - if `'pred'`, the confusion matrix is normalized over the
              predicted conditions (e.g. columns);
            - if `'all'`, the confusion matrix is normalized by the total
              number of samples;
            - if `None` (default), the confusion matrix will not be normalized.

        display_labels : array-like of shape (n_classes,), default=None
            Target names used for plotting. By default, `labels` will be used
            if it is defined, otherwise the unique labels of `y_true` and
            `y_pred` will be used.

        include_values : bool, default=True
            Includes values in confusion matrix.

        xticks_rotation : {'vertical', 'horizontal'} or float, \
                default='horizontal'
            Rotation of xtick labels.

        values_format : str, default=None
            Format specification for values in confusion matrix. If `None`, the
            format specification is 'd' or '.2g' whichever is shorter.

        cmap : str or matplotlib Colormap, default='viridis'
            Colormap recognized by matplotlib.

        ax : matplotlib Axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        colorbar : bool, default=True
            Whether or not to add a colorbar to the plot.

        im_kw : dict, default=None
            Dict with keywords passed to `matplotlib.pyplot.imshow` call.

        text_kw : dict, default=None
            Dict with keywords passed to `matplotlib.pyplot.text` call.

            .. versionadded:: 1.2

        Returns
        -------
        display : :class:`~sklearn.metrics.ConfusionMatrixDisplay`

        See Also
        --------
        ConfusionMatrixDisplay.from_predictions : Plot the confusion matrix
            given the true and predicted labels.

        Examples
        --------
        >>> import matplotlib.pyplot as plt
        >>> from sklearn.datasets import make_classification
        >>> from sklearn.metrics import ConfusionMatrixDisplay
        >>> from sklearn.model_selection import train_test_split
        >>> from sklearn.svm import SVC
        >>> X, y = make_classification(random_state=0)
        >>> X_train, X_test, y_train, y_test = train_test_split(
        ...         X, y, random_state=0)
        >>> clf = SVC(random_state=0)
        >>> clf.fit(X_train, y_train)
        SVC(random_state=0)
        >>> ConfusionMatrixDisplay.from_estimator(
        ...     clf, X_test, y_test)
        <...>
        >>> plt.show()
        """
        method_name = f"{cls.__name__}.from_estimator"
        check_matplotlib_support(method_name)
        if not is_classifier(estimator):
            raise ValueError(f"{method_name} only supports classifiers")
        y_pred = estimator.predict(X)

        return cls.from_predictions(
            y,
            y_pred,
            sample_weight=sample_weight,
            labels=labels,
            normalize=normalize,
            display_labels=display_labels,
            include_values=include_values,
            cmap=cmap,
            ax=ax,
            xticks_rotation=xticks_rotation,
            values_format=values_format,
            colorbar=colorbar,
            im_kw=im_kw,
            text_kw=text_kw,
        )

    @classmethod
    def from_predictions(
        cls,
        y_true,
        y_pred,
        *,
        labels=None,
        sample_weight=None,
        normalize=None,
        display_labels=None,
        include_values=True,
        xticks_rotation="horizontal",
        values_format=None,
        cmap="viridis",
        ax=None,
        colorbar=True,
        im_kw=None,
        text_kw=None,
    ):
        """Plot Confusion Matrix given true and predicted labels.

        Read more in the :ref:`User Guide <confusion_matrix>`.

        .. versionadded:: 1.0

        Parameters
        ----------
        y_true : array-like of shape (n_samples,)
            True labels.

        y_pred : array-like of shape (n_samples,)
            The predicted labels given by the method `predict` of an
            classifier.

        labels : array-like of shape (n_classes,), default=None
            List of labels to index the confusion matrix. This may be used to
            reorder or select a subset of labels. If `None` is given, those
            that appear at least once in `y_true` or `y_pred` are used in
            sorted order.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights.

        normalize : {'true', 'pred', 'all'}, default=None
            Either to normalize the counts display in the matrix:

            - if `'true'`, the confusion matrix is normalized over the true
              conditions (e.g. rows);
            - if `'pred'`, the confusion matrix is normalized over the
              predicted conditions (e.g. columns);
            - if `'all'`, the confusion matrix is normalized by the total
              number of samples;
            - if `None` (default), the confusion matrix will not be normalized.

        display_labels : array-like of shape (n_classes,), default=None
            Target names used for plotting. By default, `labels` will be used
            if it is defined, otherwise the unique labels of `y_true` and
            `y_pred` will be used.

        include_values : bool, default=True
            Includes values in confusion matrix.

        xticks_rotation : {'vertical', 'horizontal'} or float, \
                default='horizontal'
            Rotation of xtick labels.

        values_format : str, default=None
            Format specification for values in confusion matrix. If `None`, the
            format specification is 'd' or '.2g' whichever is shorter.

        cmap : str or matplotlib Colormap, default='viridis'
            Colormap recognized by matplotlib.

        ax : matplotlib Axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        colorbar : bool, default=True
            Whether or not to add a colorbar to the plot.

        im_kw : dict, default=None
            Dict with keywords passed to `matplotlib.pyplot.imshow` call.

        text_kw : dict, default=None
            Dict with keywords passed to `matplotlib.pyplot.text` call.

            .. versionadded:: 1.2

        Returns
        -------
        display : :class:`~sklearn.metrics.ConfusionMatrixDisplay`

        See Also
        --------
        ConfusionMatrixDisplay.from_estimator : Plot the confusion matrix
            given an estimator, the data, and the label.

        Examples
        --------
        >>> import matplotlib.pyplot as plt
        >>> from sklearn.datasets import make_classification
        >>> from sklearn.metrics import ConfusionMatrixDisplay
        >>> from sklearn.model_selection import train_test_split
        >>> from sklearn.svm import SVC
        >>> X, y = make_classification(random_state=0)
        >>> X_train, X_test, y_train, y_test = train_test_split(
        ...         X, y, random_state=0)
        >>> clf = SVC(random_state=0)
        >>> clf.fit(X_train, y_train)
        SVC(random_state=0)
        >>> y_pred = clf.predict(X_test)
        >>> ConfusionMatrixDisplay.from_predictions(
        ...    y_test, y_pred)
        <...>
        >>> plt.show()
        """
        check_matplotlib_support(f"{cls.__name__}.from_predictions")

        if display_labels is None:
            if labels is None:
                display_labels = unique_labels(y_true, y_pred)
            else:
                display_labels = labels

        cm = confusion_matrix(
            y_true,
            y_pred,
            sample_weight=sample_weight,
            labels=labels,
            normalize=normalize,
        )

        disp = cls(confusion_matrix=cm, display_labels=display_labels)

        return disp.plot(
            include_values=include_values,
            cmap=cmap,
            ax=ax,
            xticks_rotation=xticks_rotation,
            values_format=values_format,
            colorbar=colorbar,
            im_kw=im_kw,
            text_kw=text_kw,
        )


@deprecated(
    "Function `plot_confusion_matrix` is deprecated in 1.0 and will be "
    "removed in 1.2. Use one of the class methods: "
    "ConfusionMatrixDisplay.from_predictions or "
    "ConfusionMatrixDisplay.from_estimator."
)
def plot_confusion_matrix(
    estimator,
    X,
    y_true,
    *,
    labels=None,
    sample_weight=None,
    normalize=None,
    display_labels=None,
    include_values=True,
    xticks_rotation="horizontal",
    values_format=None,
    cmap="viridis",
    ax=None,
    colorbar=True,
):
    """Plot Confusion Matrix.

    `plot_confusion_matrix` is deprecated in 1.0 and will be removed in
    1.2. Use one of the following class methods:
    :func:`~sklearn.metrics.ConfusionMatrixDisplay.from_predictions` or
    :func:`~sklearn.metrics.ConfusionMatrixDisplay.from_estimator`.

    Read more in the :ref:`User Guide <confusion_matrix>`.

    Parameters
    ----------
    estimator : estimator instance
        Fitted classifier or a fitted :class:`~sklearn.pipeline.Pipeline`
        in which the last estimator is a classifier.

    X : {array-like, sparse matrix} of shape (n_samples, n_features)
        Input values.

    y_true : array-like of shape (n_samples,)
        Target values.

    labels : array-like of shape (n_classes,), default=None
        List of labels to index the matrix. This may be used to reorder or
        select a subset of labels. If `None` is given, those that appear at
        least once in `y_true` or `y_pred` are used in sorted order.

    sample_weight : array-like of shape (n_samples,), default=None
        Sample weights.

    normalize : {'true', 'pred', 'all'}, default=None
        Either to normalize the counts display in the matrix:

            - if `'true'`, the confusion matrix is normalized over the true
              conditions (e.g. rows);
            - if `'pred'`, the confusion matrix is normalized over the
              predicted conditions (e.g. columns);
            - if `'all'`, the confusion matrix is normalized by the total
              number of samples;
            - if `None` (default), the confusion matrix will not be normalized.

    display_labels : array-like of shape (n_classes,), default=None
        Target names used for plotting. By default, `labels` will be used if
        it is defined, otherwise the unique labels of `y_true` and `y_pred`
        will be used.

    include_values : bool, default=True
        Includes values in confusion matrix.

    xticks_rotation : {'vertical', 'horizontal'} or float, \
                        default='horizontal'
        Rotation of xtick labels.

    values_format : str, default=None
        Format specification for values in confusion matrix. If `None`,
        the format specification is 'd' or '.2g' whichever is shorter.

    cmap : str or matplotlib Colormap, default='viridis'
        Colormap recognized by matplotlib.

    ax : matplotlib Axes, default=None
        Axes object to plot on. If `None`, a new figure and axes is
        created.

    colorbar : bool, default=True
        Whether or not to add a colorbar to the plot.

        .. versionadded:: 0.24

    Returns
    -------
    display : :class:`~sklearn.metrics.ConfusionMatrixDisplay`
        Object that stores computed values.

    See Also
    --------
    confusion_matrix : Compute Confusion Matrix to evaluate the accuracy of a
        classification.
    ConfusionMatrixDisplay : Confusion Matrix visualization.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from sklearn.datasets import make_classification
    >>> from sklearn.metrics import plot_confusion_matrix
    >>> from sklearn.model_selection import train_test_split
    >>> from sklearn.svm import SVC
    >>> X, y = make_classification(random_state=0)
    >>> X_train, X_test, y_train, y_test = train_test_split(
    ...         X, y, random_state=0)
    >>> clf = SVC(random_state=0)
    >>> clf.fit(X_train, y_train)
    SVC(random_state=0)
    >>> plot_confusion_matrix(clf, X_test, y_test)  # doctest: +SKIP
    >>> plt.show()
    """
    check_matplotlib_support("plot_confusion_matrix")

    if not is_classifier(estimator):
        raise ValueError("plot_confusion_matrix only supports classifiers")

    y_pred = estimator.predict(X)
    cm = confusion_matrix(
        y_true, y_pred, sample_weight=sample_weight, labels=labels, normalize=normalize
    )

    if display_labels is None:
        if labels is None:
            display_labels = unique_labels(y_true, y_pred)
        else:
            display_labels = labels

    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=display_labels)
    return disp.plot(
        include_values=include_values,
        cmap=cmap,
        ax=ax,
        xticks_rotation=xticks_rotation,
        values_format=values_format,
        colorbar=colorbar,
    )
