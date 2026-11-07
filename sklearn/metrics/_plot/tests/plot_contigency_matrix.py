from itertools import product

import numpy as np

from ...cluster._supervised import contingency_matrix
from ....utils import check_matplotlib_support
from ....utils.multiclass import unique_labels
from ....utils.validation import _deprecate_positional_args
from ....base import is_classifier


class contingencyMatrixDisplay:
    @_deprecate_positional_args
    def __init__(self, contingency_matrix, *, display_labels=None):
        self.contingency_matrix = contingency_matrix
        self.display_labels = display_labels

    @_deprecate_positional_args
    def plot(self, *, include_values=True, cmap='viridis',
             xticks_rotation='horizontal', values_format=None,
             ax=None, colorbar=True):
        check_matplotlib_support("contingencyMatrixDisplay.plot")
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        cm = self.contingency_matrix
        n_classes = cm.shape[0]
        self.im_ = ax.imshow(cm, interpolation='nearest', cmap=cmap)
        self.text_ = None
        cmap_min, cmap_max = self.im_.cmap(0), self.im_.cmap(256)

        if include_values:
            self.text_ = np.empty_like(cm, dtype=object)

            # print text with appropriate color depending on background
            thresh = (cm.max() + cm.min()) / 2.0

            for i, j in product(range(n_classes), range(n_classes)):
                color = cmap_max if cm[i, j] < thresh else cmap_min

                if values_format is None:
                    text_cm = format(cm[i, j], '.2g')
                    if cm.dtype.kind != 'f':
                        text_d = format(cm[i, j], 'd')
                        if len(text_d) < len(text_cm):
                            text_cm = text_d
                else:
                    text_cm = format(cm[i, j], values_format)

                self.text_[i, j] = ax.text(
                    j, i, text_cm,
                    ha="center", va="center",
                    color=color)

        if self.display_labels is None:
            display_labels = np.arange(n_classes)
        else:
            display_labels = self.display_labels
        if colorbar:
            fig.colorbar(self.im_, ax=ax)
        ax.set(xticks=np.arange(n_classes),
               yticks=np.arange(n_classes),
               xticklabels=display_labels,
               yticklabels=display_labels,
               ylabel="True label",
               xlabel="Predicted label")

        ax.set_ylim((n_classes - 0.5, -0.5))
        plt.setp(ax.get_xticklabels(), rotation=xticks_rotation)

        self.figure_ = fig
        self.ax_ = ax
        return self


@_deprecate_positional_args
def plot_contingency_matrix(estimator, X, y_true, *, labels=None,
                            sample_weight=None, normalize=None,
                            display_labels=None, include_values=True,
                            xticks_rotation='horizontal',
                            values_format=None,
                            cmap='inferno', ax=None, colorbar=True):
    check_matplotlib_support("plot_contingency_matrix")

    if not is_classifier(estimator):
        raise ValueError("plot_contingency_matrix only supports classifiers")

    y_pred = estimator.predict(X)
    cm = contingency_matrix(y_true,
                            y_pred,
                            eps=None, sparse=False)

    if display_labels is None:
        if labels is None:
            display_labels = unique_labels(y_true, y_pred)
        else:
            display_labels = labels

    disp = contingencyMatrixDisplay(contingency_matrix=cm,
                                    display_labels=display_labels)
    return disp.plot(include_values=include_values,
                     cmap=cmap, ax=ax, xticks_rotation=xticks_rotation,
                     values_format=values_format, colorbar=colorbar)
