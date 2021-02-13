from itertools import product

import numpy as np

from ..utils import check_matplotlib_support
from ..utils.validation import _deprecate_positional_args


@_deprecate_positional_args
def plot_heatmap(
    data,
    ylabel,
    xlabel,
    ytick_labels,
    xtick_labels,
    *,
    xticks_rotation='horizontal',
    ax=None,
    cmap='viridis',
    include_values=True,
    values_format=None,
    colorbar=True,
):
    """Plot heatmap visualization.

    Parameters
    ----------
    data : array-like
        Data to plot the heatmap.

    ylabel : str
        Label for y-axis.

    xlabel : str
        Label for x-axis.

    ytick_labels : array-like of str
        Labels for y-axis ticks.

    xtick_labels : array-like of str
        Labels for x-axis ticks.

    xticks_rotation : {'vertical', 'horizontal'} or float, \
                     default='horizontal'
        Rotation of xtick labels.

    ax : matplotlib axes, default=None
        Axes object to plot on. If `None`, a new figure and axes is
        created.

    cmap : str or matplotlib Colormap, default='viridis'
        Colormap recognized by matplotlib.

    include_values : bool, default=True
        Includes values in confusion matrix.

    values_format : str, default=None
        Format specification for values in confusion matrix. If `None`,
        the format specification is 'd' or '.2g' whichever is shorter.

    colorbar : bool, default=True
        Whether or not to add a colorbar to the plot.

    Returns
    -------
    fig : matplotlib Figure
        Figure containing the heatmap.

    ax : matplotlib Axes
        Axes with the heatmap.

    im : matplotlib AxesImage
        Image of the heatmap.

    text : ndarray of shape (n_classes, n_classes), dtype=matplotlib Text, \
            or None
        Array of matplotlib axes. `None` if `include_values` is false.
    """
    check_matplotlib_support("utils.plot.plot_heatmap")
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    n_classes = data.shape[0]
    im = ax.imshow(data, interpolation='nearest', cmap=cmap)
    text = None
    cmap_min, cmap_max = im.cmap(0), im.cmap(256)

    if include_values:
        text = np.empty_like(data, dtype=object)

        # print text with appropriate color depending on background
        thresh = (data.max() + data.min()) / 2.0

        for i, j in product(range(n_classes), range(n_classes)):
            color = cmap_max if data[i, j] < thresh else cmap_min

            if values_format is None:
                text_data = format(data[i, j], '.2g')
                if data.dtype.kind != 'f':
                    text_d = format(data[i, j], 'd')
                    if len(text_d) < len(text_data):
                        text_data = text_d
            else:
                text_data = format(data[i, j], values_format)
            text[i, j] = ax.text(
                j, i, text_data,
                ha="center", va="center",
                color=color)

    if colorbar:
        fig.colorbar(im, ax=ax)
    ax.set(xticks=np.arange(n_classes),
           yticks=np.arange(n_classes),
           xticklabels=xtick_labels,
           yticklabels=ytick_labels,
           ylabel=ylabel,
           xlabel=xlabel)

    plt.setp(ax.get_xticklabels(), rotation=xticks_rotation)

    return fig, ax, im, text
