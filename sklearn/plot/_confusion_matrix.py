import numpy as np
from sklearn.metrics import confusion_matrix
from sklearn.plot import plot_heatmap


def plot_confusion_matrix(y_true, y_pred, classes, sample_weight=None,
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

    classes : list of strings
        The list of classes represented in the two-dimensional input array.

    sample_weight : array-like of shape = [n_samples], optional
        Sample weights used to calculate the confusion matrix

    normalize : boolean, default=False
        If True, the confusion matrix will be normalized by row.

    xlabel : string, default="Predicted Label"
        Label for the x-axis.

    ylabel : string, default="True Label"
        Label for the y-axis.

    title : string, default="Confusion matrix"
        Title for the heatmap.

    cmap : string or colormap
        Matpotlib colormap to use.

    vmin : int, float or None
        Minimum clipping value. This argument will be passed on to the
        pcolormesh function from matplotlib used to generate the heatmap.

    vmax : int, float or None
        Maximum clipping value. This argument will be passed on to the
        pcolormesh function from matplotlib used to generate the heatmap.

    ax : axes object or None
        Matplotlib axes object to plot into. If None, the current axes are
        used.

    fmt : string, default="{:.2f}"
        Format string to convert value to text. This will be ignored if
        normalize argument is False.

    xtickrotation : float, default=45
        Rotation of the xticklabels.

    norm : matplotlib normalizer
        Normalizer passed to pcolormesh function from matplotlib used to
        generate the heatmap.
    """

    import matplotlib.pyplot as plt

    values = confusion_matrix(y_true, y_pred, sample_weight=sample_weight)

    if normalize:
        values = values.astype('float') / values.sum(axis=1)[:, np.newaxis]

    fmt = fmt if normalize else '{:d}'

    plot_heatmap(values, xticklabels=classes, yticklabels=classes, cmap=cmap,
                 xlabel=xlabel, ylabel=ylabel, vmin=vmin, vmax=vmax, ax=ax,
                 fmt=fmt, xtickrotation=xtickrotation, norm=norm)

    plt.title(title)
