# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np

from sklearn.base import is_classifier
from sklearn.utils.multiclass import unique_labels
from sklearn.metrics import confusion_matrix

import plotly.graph_objects as go
import plotly.express as px


class ConfusionMatrixSankeyDisplay:
    """Confusion Matrix Sankey visualization.

    It is recommended to use
    :func:`~sklearn.metrics.ConfusionMatrixSankeyDisplay.from_estimator` or
    :func:`~sklearn.metrics.ConfusionMatrixSankeyDisplay.from_predictions` to
    create a :class:`ConfusionMatrixSankeyDisplay`. All parameters are stored
    as attributes.

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
    figure_ : plotly.graph_objects.Figure
        Figure containing the Sankey diagram.

    See Also
    --------
    confusion_matrix : Compute Confusion Matrix to evaluate the accuracy of a
        classification.
    ConfusionMatrixDisplay: Traditional confusion matrix heatmap visualization.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from sklearn.datasets import make_classification
    >>> from sklearn.model_selection import train_test_split
    >>> from sklearn.ensemble import RandomForestClassifier
    >>> from sklearn.metrics import ConfusionMatrixSankeyDisplay
    >>> X, y = make_classification(random_state=0)
    >>> X_train, X_test, y_train, y_test
        = train_test_split(X, y, random_state=0)
    >>> clf = RandomForestClassifier(random_state=0)
    >>> clf.fit(X_train, y_train)
    >>> disp = ConfusionMatrixSankeyDisplay.from_estimator(clf, X_test, y_test)
    ...
    >>> disp.plot()
    <...>
    """

    def __init__(self, confusion_matrix, *, display_labels=None,
                 save_path="confusion_matrix_sankey.html"):
        self.cm = confusion_matrix
        self.display_labels = (display_labels if display_labels is not None
                               else np.arange(len(confusion_matrix)))
        self.save_path = save_path

    def plot(
            self,
            *,
            cmap=None,
            values_format='d',
            width=800,
            height=600,
            **kwargs
    ):
        """Plot visualization.

        Parameters
        ----------
        cmap : str, default=None
            Color scale name from plotly.express.colors.qualitative. If None,
            uses default Plotly color sequence.

        values_format : str, default='d'
            Format specification for values in Sankey flows. If format contains
            '%', percentage will be shown automatically.

        width : int, default=800
            Figure width in pixels.

        height : int, default=600
            Figure height in pixels.

        **kwargs : dict
            Additional keyword arguments passed to plotly.graph_objects.Sankey.

        Returns
        -------
        self : object
            Returns self.
        """

        n_classes = len(self.display_labels)
        labels = [str(label) for label in self.display_labels]

        # Create color sequence
        if cmap is None:
            colors = px.colors.qualitative.Plotly
        else:
            colors = getattr(px.colors.qualitative, cmap)

        # Create Sankey nodes (true on left, predicted on right)
        nodes = dict(
            label=labels + labels,
            x=[0] * n_classes + [1] * n_classes,  # Positions
            y=np.linspace(0, 1, n_classes).tolist() * 2,
            color=[colors[i % len(colors)] for i in range(n_classes)] * 2
        )

        # Create Sankey links
        sources, targets, values, link_colors = [], [], [], []
        for i in range(n_classes):
            for j in range(n_classes):
                if self.cm[i, j] > 0:
                    sources.append(i)
                    targets.append(n_classes + j)
                    values.append(self.cm[i, j])
                    link_colors.append(
                        f"rgba({self._hex_to_rgb(colors[i % len(colors)])}"
                        f",0.5)")

        fig = go.Figure(go.Sankey(
            node=nodes,
            link=dict(
                source=sources,
                target=targets,
                value=values,
                color=link_colors,
                hovertemplate='True: %{source.label}<br>Predicted:'
                              ' %{target.label}<br>Count: %{value}' +
                              ('<br>Percentage: %{value:.1f}%'
                               if values_format.endswith('%') else '')
            ),
            arrangement='snap'
        ))

        fig.write_html(self.save_path)
        print(f"figure saved to：{self.save_path}")
        return self

    def _hex_to_rgb(self, hex_color):
        """Convert hexadecimal color code to RGB string.

        Parameters
        ----------
        hex_color : str
            Hexadecimal color code (with or without leading #)

        Returns
        -------
        str
            Comma-separated RGB values (e.g., '255,0,0')
        """
        hex_color = hex_color.lstrip('#')
        return ','.join(str(int(hex_color[i:i + 2], 16)) for i in (0, 2, 4))

    @classmethod
    def from_estimator(
            cls,
            estimator,
            X,
            y,
            **kwargs
    ):
        """Plot Confusion Matrix Sankey given an estimator and some data.

        Parameters
        ----------
        estimator : estimator instance
            Fitted classifier or a fitted :class:`~sklearn.pipeline.Pipeline`
            in which the last estimator is a classifier.

        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Input values.

        y : array-like of shape (n_samples,)
            Target values.

        **kwargs : dict
            Additional keyword arguments passed to the `plot` method.

        Returns
        -------
        display : :class:`~sklearn.metrics.ConfusionMatrixSankeyDisplay`
        """
        if not is_classifier(estimator):
            raise ValueError("Estimator must be a classifier")
        return cls.from_predictions(y, estimator.predict(X), **kwargs)

    @classmethod
    def from_predictions(
            cls,
            y_true,
            y_pred,
            *,
            labels=None,
            normalize=None,
            display_labels=None,
            **kwargs
    ):
        """Plot Confusion Matrix Sankey given true and predicted labels.

        Parameters
        ----------
        y_true : array-like of shape (n_samples,)
            True labels.

        y_pred : array-like of shape (n_samples,)
            Predicted labels.

        labels : array-like of shape (n_classes,), default=None
            List of labels to index the confusion matrix.

        normalize : {'true', 'pred', 'all'}, default=None
            Normalization of confusion matrix counts:
            - 'true': normalizes over true labels
            - 'pred': normalizes over predicted labels
            - 'all': normalizes over total samples

        **kwargs : dict
            Additional keyword arguments passed to the `plot` method.

        Returns
        -------
        display : :class:`~sklearn.metrics.ConfusionMatrixSankeyDisplay`
        """
        cm = confusion_matrix(y_true, y_pred,
                              labels=labels, normalize=normalize)

        # Scale normalized values for percentage display
        if normalize:
            cm = (cm * 100).astype(int)
            kwargs['values_format'] = '.1f%'
        display_labels = (
            display_labels if display_labels is not None 
            else labels if labels is not None 
            else unique_labels(y_true, y_pred))
        return cls(cm, display_labels=display_labels)
