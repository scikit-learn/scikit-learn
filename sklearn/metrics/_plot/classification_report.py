# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np

from sklearn.base import is_classifier
from sklearn.metrics import classification_report
from sklearn.utils._optional_dependencies import check_matplotlib_support
from sklearn.utils._plotting import _validate_style_kwargs

# Score columns share the 0-1 colormap; support is a count, shown as text only.
_SCORE_COLUMNS = ("precision", "recall", "f1-score")


class ClassificationReportDisplay:
    """Classification report visualization.

    A heatmap of per-class ``precision``, ``recall`` and ``f1-score`` (all on a
    shared 0-1 color scale), with ``support`` shown as a text-only column. Rows
    are the individual classes followed by the ``macro avg`` and ``weighted avg``
    summary rows. The scalar ``accuracy`` entry of the report is not part of the
    grid.

    It is recommended to use
    :func:`~sklearn.metrics.ClassificationReportDisplay.from_estimator` or
    :func:`~sklearn.metrics.ClassificationReportDisplay.from_predictions` to
    create a :class:`ClassificationReportDisplay`. To visualize an existing
    report, pass the output of
    :func:`~sklearn.metrics.classification_report` (with ``output_dict=True``)
    to the constructor.

    For general information regarding `scikit-learn` visualization tools, see
    the :ref:`Visualization Guide <visualizations>`.

    Parameters
    ----------
    report : dict
        Classification report as returned by
        :func:`~sklearn.metrics.classification_report` with ``output_dict=True``.

    display_labels : array-like of shape (n_rows,), default=None
        Row labels for the plot. If `None`, the keys of `report` (excluding
        `accuracy`) are used.

    Attributes
    ----------
    im_ : matplotlib AxesImage
        Image representing the precision/recall/f1-score matrix.

    text_ : ndarray of matplotlib Text, or None
        Array of the annotation artists. `None` if `include_values` is False.

    display_labels : ndarray
        Resolved row labels.

    ax_ : matplotlib Axes
        Axes with the classification report.

    figure_ : matplotlib Figure
        Figure containing the classification report.

    See Also
    --------
    classification_report : Build a text/dict report of the main classification
        metrics.
    ClassificationReportDisplay.from_estimator : Plot the classification report
        given an estimator, the data, and the labels.
    ClassificationReportDisplay.from_predictions : Plot the classification report
        given the true and predicted labels.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from sklearn.datasets import make_classification
    >>> from sklearn.metrics import classification_report
    >>> from sklearn.metrics import ClassificationReportDisplay
    >>> from sklearn.model_selection import train_test_split
    >>> from sklearn.svm import SVC
    >>> X, y = make_classification(random_state=0)
    >>> X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    >>> clf = SVC(random_state=0).fit(X_train, y_train)
    >>> report = classification_report(y_test, clf.predict(X_test),
    ...                                output_dict=True)
    >>> disp = ClassificationReportDisplay(report)
    >>> disp.plot()
    <...>
    >>> plt.show()
    """

    def __init__(self, report, *, display_labels=None):
        self.report = report
        self.display_labels = display_labels

    def _rows(self):
        """Grid rows: every report entry except the scalar ``accuracy``."""
        return [key for key in self.report if key != "accuracy"]

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
            Whether to annotate the score cells and show the `support` column.

        cmap : str or matplotlib Colormap, default='viridis'
            Colormap recognized by matplotlib.

        xticks_rotation : {'vertical', 'horizontal'} or float, \
                default='horizontal'
            Rotation of the column (xtick) labels.

        values_format : str, default=None
            Format specification for the score values. If `None`, `'.2f'` is
            used.

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is created.

        colorbar : bool, default=True
            Whether or not to add a colorbar to the plot.

        im_kw : dict, default=None
            Dict with keywords passed to `matplotlib.pyplot.imshow`.

        text_kw : dict, default=None
            Dict with keywords passed to `matplotlib.pyplot.text`.

        Returns
        -------
        display : :class:`~sklearn.metrics.ClassificationReportDisplay`
            Object that stores the computed values.
        """
        check_matplotlib_support("ClassificationReportDisplay.plot")
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        rows = self._rows()
        if self.display_labels is None:
            self.display_labels = np.asarray(rows, dtype=object)
        n_rows = len(rows)
        n_scores = len(_SCORE_COLUMNS)

        scores = np.array(
            [[self.report[row][col] for col in _SCORE_COLUMNS] for row in rows],
            dtype=float,
        )

        default_im_kw = dict(interpolation="nearest", cmap=cmap, vmin=0.0, vmax=1.0)
        im_kw = _validate_style_kwargs(default_im_kw, im_kw or {})
        text_kw = text_kw or {}

        self.im_ = ax.imshow(scores, **im_kw)
        self.text_ = None
        cmap_min, cmap_max = self.im_.cmap(0.0), self.im_.cmap(1.0)

        if include_values:
            self.text_ = np.empty((n_rows, n_scores + 1), dtype=object)
            for i in range(n_rows):
                for j in range(n_scores):
                    color = cmap_max if scores[i, j] < 0.5 else cmap_min
                    text = format(scores[i, j], values_format or ".2f")
                    default_text_kw = dict(ha="center", va="center", color=color)
                    kw = _validate_style_kwargs(default_text_kw, text_kw)
                    self.text_[i, j] = ax.text(j, i, text, **kw)
                # support: a count, drawn on the (blank) axes background as text
                support = self.report[rows[i]]["support"]
                default_text_kw = dict(ha="center", va="center", color="black")
                kw = _validate_style_kwargs(default_text_kw, text_kw)
                self.text_[i, n_scores] = ax.text(
                    n_scores, i, format(int(support), "d"), **kw
                )
            col_labels = list(_SCORE_COLUMNS) + ["support"]
        else:
            col_labels = list(_SCORE_COLUMNS)

        n_cols = len(col_labels)
        if colorbar:
            fig.colorbar(self.im_, ax=ax)
        ax.set(
            xticks=np.arange(n_cols),
            yticks=np.arange(n_rows),
            xticklabels=col_labels,
            yticklabels=self.display_labels,
            ylabel="Class",
            xlabel="Metric",
        )
        ax.set_xlim(-0.5, n_cols - 0.5)
        ax.set_ylim(n_rows - 0.5, -0.5)
        plt.setp(ax.get_xticklabels(), rotation=xticks_rotation)

        self.figure_ = fig
        self.ax_ = ax
        return self

    @classmethod
    def from_estimator(
        cls,
        estimator,
        X,
        y,
        *,
        labels=None,
        target_names=None,
        sample_weight=None,
        zero_division="warn",
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
        """Plot the classification report given an estimator and some data.

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
            Optional list of label indices to include in the report.

        target_names : array-like of shape (n_classes,), default=None
            Optional display names matching the labels (same order).

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights.

        zero_division : {"warn", 0.0, 1.0, np.nan}, default="warn"
            Value to use when there is a zero division. Passed through to
            :func:`~sklearn.metrics.classification_report`.

        display_labels : array-like of shape (n_rows,), default=None
            Row labels for the plot. If `None`, the report keys are used.

        include_values : bool, default=True
            Whether to annotate the score cells and show the `support` column.

        xticks_rotation : {'vertical', 'horizontal'} or float, \
                default='horizontal'
            Rotation of the column (xtick) labels.

        values_format : str, default=None
            Format specification for the score values. If `None`, `'.2f'` is
            used.

        cmap : str or matplotlib Colormap, default='viridis'
            Colormap recognized by matplotlib.

        ax : matplotlib Axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is created.

        colorbar : bool, default=True
            Whether or not to add a colorbar to the plot.

        im_kw : dict, default=None
            Dict with keywords passed to `matplotlib.pyplot.imshow`.

        text_kw : dict, default=None
            Dict with keywords passed to `matplotlib.pyplot.text`.

        Returns
        -------
        display : :class:`~sklearn.metrics.ClassificationReportDisplay`

        See Also
        --------
        ClassificationReportDisplay.from_predictions : Plot the classification
            report given the true and predicted labels.

        Examples
        --------
        >>> import matplotlib.pyplot as plt
        >>> from sklearn.datasets import make_classification
        >>> from sklearn.metrics import ClassificationReportDisplay
        >>> from sklearn.model_selection import train_test_split
        >>> from sklearn.svm import SVC
        >>> X, y = make_classification(random_state=0)
        >>> X_train, X_test, y_train, y_test = train_test_split(
        ...     X, y, random_state=0)
        >>> clf = SVC(random_state=0).fit(X_train, y_train)
        >>> ClassificationReportDisplay.from_estimator(clf, X_test, y_test)
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
            labels=labels,
            target_names=target_names,
            sample_weight=sample_weight,
            zero_division=zero_division,
            display_labels=display_labels,
            include_values=include_values,
            xticks_rotation=xticks_rotation,
            values_format=values_format,
            cmap=cmap,
            ax=ax,
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
        target_names=None,
        sample_weight=None,
        zero_division="warn",
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
        """Plot the classification report given true and predicted labels.

        Parameters
        ----------
        y_true : array-like of shape (n_samples,)
            True labels.

        y_pred : array-like of shape (n_samples,)
            The predicted labels given by the method `predict` of a classifier.

        labels : array-like of shape (n_classes,), default=None
            Optional list of label indices to include in the report.

        target_names : array-like of shape (n_classes,), default=None
            Optional display names matching the labels (same order).

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights.

        zero_division : {"warn", 0.0, 1.0, np.nan}, default="warn"
            Value to use when there is a zero division. Passed through to
            :func:`~sklearn.metrics.classification_report`.

        display_labels : array-like of shape (n_rows,), default=None
            Row labels for the plot. If `None`, the report keys are used.

        include_values : bool, default=True
            Whether to annotate the score cells and show the `support` column.

        xticks_rotation : {'vertical', 'horizontal'} or float, \
                default='horizontal'
            Rotation of the column (xtick) labels.

        values_format : str, default=None
            Format specification for the score values. If `None`, `'.2f'` is
            used.

        cmap : str or matplotlib Colormap, default='viridis'
            Colormap recognized by matplotlib.

        ax : matplotlib Axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is created.

        colorbar : bool, default=True
            Whether or not to add a colorbar to the plot.

        im_kw : dict, default=None
            Dict with keywords passed to `matplotlib.pyplot.imshow`.

        text_kw : dict, default=None
            Dict with keywords passed to `matplotlib.pyplot.text`.

        Returns
        -------
        display : :class:`~sklearn.metrics.ClassificationReportDisplay`

        See Also
        --------
        ClassificationReportDisplay.from_estimator : Plot the classification
            report given an estimator, the data, and the labels.

        Examples
        --------
        >>> import matplotlib.pyplot as plt
        >>> from sklearn.datasets import make_classification
        >>> from sklearn.metrics import ClassificationReportDisplay
        >>> from sklearn.model_selection import train_test_split
        >>> from sklearn.svm import SVC
        >>> X, y = make_classification(random_state=0)
        >>> X_train, X_test, y_train, y_test = train_test_split(
        ...     X, y, random_state=0)
        >>> clf = SVC(random_state=0).fit(X_train, y_train)
        >>> y_pred = clf.predict(X_test)
        >>> ClassificationReportDisplay.from_predictions(y_test, y_pred)
        <...>
        >>> plt.show()
        """
        check_matplotlib_support(f"{cls.__name__}.from_predictions")

        report = classification_report(
            y_true,
            y_pred,
            labels=labels,
            target_names=target_names,
            sample_weight=sample_weight,
            zero_division=zero_division,
            output_dict=True,
        )

        disp = cls(report, display_labels=display_labels)
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
