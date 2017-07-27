import numpy as np
from sklearn.metrics import confusion_matrix
from sklearn.plot import plot_heatmap
from sklearn.utils.multiclass import unique_labels


def plot_confusion_matrix(y_true, y_pred, classes=None, sample_weight=None,
                          normalize=False,
                          xlabel="Predicted Label", ylabel="True Label",
                          title='Confusion matrix', cmap=None, vmin=None,
                          vmax=None, ax=None, fmt="{:.2f}",
                          xtickrotation=45, norm=None):
    """Plot the confusion matrix as a heatmap. A confusion matrix is computed
    using `y_true`, `y_pred` and `sample_weights` arguments. Normalization
    can be applied by setting `normalize=True`.

    Parameters
    ----------
    y_true : array, shape = [n_samples]
        Ground truth (correct) target values.

    y_pred : array, shape = [n_samples]
        Estimated targets as returned by a classifier.

    classes : list of strings, optional (default=None)
        The list of names of classes represented in the two-dimensional input
        array. If not passed in function call, the classes will be infered
        from y_true and y_pred

    sample_weight : array-like of shape = [n_samples], optional (default=None)
        Sample weights used to calculate the confusion matrix

    normalize : boolean, optional (default=False)
        If True, the confusion matrix will be normalized by row.

    xlabel : string, optional (default="Predicted Label")
        Label for the x-axis.

    ylabel : string, optional (default="True Label")
        Label for the y-axis.

    title : string, optional (default="Confusion matrix")
        Title for the heatmap.

    cmap : string or colormap, optional (default=None)
        Matpotlib colormap to use. If None, plt.cm.hot will be used.

    vmin : int, float or None, optional (default=None)
        Minimum clipping value. This argument will be passed on to the
        pcolormesh function from matplotlib used to generate the heatmap.

    vmax : int, float or None, optional (default=None)
        Maximum clipping value. This argument will be passed on to the
        pcolormesh function from matplotlib used to generate the heatmap.

    ax : axes object or None, optional (default=None)
        Matplotlib axes object to plot into. If None, the current axes are
        used.

    fmt : string, optional (default="{:.2f}")
        Format string to convert value to text. This will be ignored if
        normalize argument is False.

    xtickrotation : float, optional (default=45)
        Rotation of the xticklabels.

    norm : matplotlib normalizer, optional (default=None)
        Normalizer passed to pcolormesh function from matplotlib used to
        generate the heatmap.
    """

    import matplotlib.pyplot as plt

    unique_y = unique_labels(y_true, y_pred)

    if classes is None:
        classes = unique_y
    else:
        if len(classes) != len(unique_y):
            raise ValueError("y_true and y_pred contain %d unique classes,"
                             "which is not the same as %d"
                             "classes found in `classes=%s` paramter" %
                             (len(classes), len(unique_y), unique_y))

    values = confusion_matrix(y_true, y_pred, sample_weight=sample_weight)

    if normalize:
        values = values.astype('float') / values.sum(axis=1)[:, np.newaxis]

    fmt = fmt if normalize else '{:d}'

    img = plot_heatmap(values, xticklabels=classes, yticklabels=classes,
                       cmap=cmap, xlabel=xlabel, ylabel=ylabel, vmin=vmin,
                       vmax=vmax, ax=ax, fmt=fmt, xtickrotation=xtickrotation,
                       norm=norm)

    plt.title(title)

    return img
