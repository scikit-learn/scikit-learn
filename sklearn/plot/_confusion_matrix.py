import numpy as np
from sklearn.plot import plot_heatmap


def plot_confusion_matrix(values, classes, normalize=True,
                          xlabel="Predicted Label", ylabel="True Label",
                          title='Confusion matrix', cmap=None, vmin=None,
                          vmax=None, ax=None, fmt="{:.2f}",
                          xtickrotation=45, norm=None):
    """Print and plot the confusion matrix. Normalization can be applied by
    setting `normalize=True`.

    Parameters
    ----------
    values : ndarray
        Two-dimensional array to visualize.

    classes : list of strings
        The list of classes represented in the two-dimensional input array.

    normalize : boolean, default=True
        If True, the confusion matrix will be normalized by row.

    xlabel : string, default=""
        Label for the x-axis.

    ylabel : string, default=""
        Label for the y-axis.

    cmap : string or colormap
        Matpotlib colormap to use.

    vmin : int, float or None
        Minimum clipping value.

    vmax : int, float or None
        Maximum clipping value.

    ax : axes object or None
        Matplotlib axes object to plot into. If None, the current axes are
        used.

    fmt : string, default="{:.2f}"
        Format string to convert value to text. This will be ignored if
        normalize argument is False.

    xtickrotation : float, default=45
        Rotation of the xticklabels.

    norm : matplotlib normalizer
        Normalizer passed to pcolor
    """

    import matplotlib.pyplot as plt

    if normalize:
        values = values.astype('float') / values.sum(axis=1)[:, np.newaxis]

    print(title)
    print(values)

    fmt = fmt if normalize else '{:d}'

    plot_heatmap(values, xticklabels=classes, yticklabels=classes, cmap=cmap,
                 xlabel=xlabel, ylabel=ylabel, vmin=vmin, vmax=vmax, ax=ax,
                 fmt=fmt, xtickrotation=xtickrotation, norm=norm)

    plt.title(title)
