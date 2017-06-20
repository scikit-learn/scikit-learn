import numpy as np
from sklearn.plot import plot_heatmap


def plot_confusion_matrix(values, classes, normalize=False,
                          xlabel="Predicted Label", ylabel="True Label",
                          title='Confusion matrix', cmap=None, vmin=None,
                          vmax=None, ax=None, fmt="{:.2f}",
                          xtickrotation=45, norm=None):
    """Print and plot the confusion matrix as a heatmap. Normalization can be
    applied by setting `normalize=True`.

    Parameters
    ----------
    values : ndarray
        Two-dimensional array to visualize.

    classes : list of strings
        The list of classes represented in the two-dimensional input array.

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

    if normalize:
        values = values.astype('float') / values.sum(axis=1)[:, np.newaxis]

    print(title)
    print(values)

    fmt = fmt if normalize else '{:d}'

    plot_heatmap(values, xticklabels=classes, yticklabels=classes, cmap=cmap,
                 xlabel=xlabel, ylabel=ylabel, vmin=vmin, vmax=vmax, ax=ax,
                 fmt=fmt, xtickrotation=xtickrotation, norm=norm)

    plt.title(title)
