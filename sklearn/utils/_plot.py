import numpy as np

from ..utils import check_matplotlib_support


def plot_heatmap(
    data,
    *,
    ylabel,
    xlabel,
    yticklabels,
    xticklabels,
    xticks_rotation="horizontal",
    ax=None,
    cmap="viridis",
    include_values=True,
    values_format=None,
    colorbar=True,
    im_kwargs=None,
):
    """Plot a heatmap.

    .. versionadded:: 1.1

    Parameters
    ----------
    data : ndarray of shape (n_rows, n_columns)
        Data to plot the heatmap.

    ylabel : str
        Label for y-axis.

    xlabel : str
        Label for x-axis.

    yticklabels : array-like of shape (n_rows,), dtype=str
        Labels for y-axis ticks.

    xticklabels : array-like of shape (n_columns,), dtype=str
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

    im_kwargs : dict, default=None
        Dict with keywords passed to `matplotlib.pyplot.imshow` call.

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
        Array of matplotlib axes. `None` if `include_values` is `False`.
    """
    check_matplotlib_support("utils.plot.plot_heatmap")
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    default_im_kwargs = dict(interpolation="nearest", cmap=cmap)
    im_kwargs = im_kwargs or {}
    im_kwargs = {**default_im_kwargs, **im_kwargs}

    im = ax.imshow(data, **im_kwargs)
    text = None
    cmap_min, cmap_max = im.cmap(0), im.cmap(1.0)

    if include_values:
        text = np.empty_like(data, dtype=object)

        # print text with appropriate color depending on background
        thresh = (data.max() + data.min()) / 2.0

        for flat_index in range(data.size):
            row, col = np.unravel_index(flat_index, data.shape)
            color = cmap_max if data[row, col] < thresh else cmap_min

            if values_format is None:
                text_data = format(data[row, col], ".2g")
                if data.dtype.kind != "f":
                    text_d = format(data[row, col], "d")
                    if len(text_d) < len(text_data):
                        text_data = text_d
            else:
                text_data = format(data[row, col], values_format)
            text[row, col] = ax.text(
                col, row, text_data, ha="center", va="center", color=color
            )

    if colorbar:
        fig.colorbar(im, ax=ax)
    ax.set(
        xticks=np.arange(len(xticklabels)),
        yticks=np.arange(len(yticklabels)),
        xticklabels=xticklabels,
        yticklabels=yticklabels,
        ylabel=ylabel,
        xlabel=xlabel,
    )

    plt.setp(ax.get_xticklabels(), rotation=xticks_rotation)

    return fig, ax, im, text
