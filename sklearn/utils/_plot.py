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
    im_kw=None,
    text_kw=None,
):
    """Plot a heatmap.

    .. versionadded:: 1.2

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

    im_kw : dict, default=None
        Dict with keywords passed to `matplotlib.pyplot.imshow` call.

    text_kw : dict, default=None
            Dict with keywords passed to `matplotlib.pyplot.text` call.

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

    default_im_kw = dict(interpolation="nearest", cmap=cmap)
    im_kw = im_kw or {}
    im_kw = {**default_im_kw, **im_kw}
    text_kw = text_kw or {}

    im = ax.imshow(data, **im_kw)
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

            default_text_kwargs = dict(ha="center", va="center", color=color)
            text_kw = text_kw or {}
            text_kwargs = {**default_text_kwargs, **text_kw}
            text[row, col] = ax.text(col, row, text_data, **text_kwargs)

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
