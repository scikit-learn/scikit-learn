import numpy as np


def plot_heatmap(values, xlabel="", ylabel="", xticklabels=None,
                 yticklabels=None, cmap=None, vmin=None, vmax=None, ax=None,
                 fmt="{:.2f}"):
    """Plot a matrix as heatmap with explicit numbers.

    Parameters
    ----------
    values : ndarray
        Two-dimensional array to visualize.

    xlabel : string, default=""
        Label for the x-axis.

    ylabel : string, default=""
        Label for the y-axis.

    xticklabels : list of string or None, default=None
        Tick labels for the x-axis.

    yticklabels : list of string or None, default=None
        Tick labels for the y-axis

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
        Format string to convert value to text.
    """
    import matplotlib.pyplot as plt
    if ax is None:
        ax = plt.gca()
    img = ax.pcolor(values, cmap=cmap, vmin=None, vmax=None)
    # this will allow us to access the pixel values:
    img.update_scalarmappable()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if xticklabels is None:
        xticklabels = [""] * values.shape[1]
    if yticklabels is None:
        yticklabels = [""] * values.shape[0]

    ax.xaxis.set_ticks_position('bottom')

    # +.5 makes the ticks centered on the pixels
    ax.set_xticks(np.arange(values.shape[1]) + .5)
    ax.set_xticklabels(xticklabels, ha="center")
    ax.set_yticks(np.arange(values.shape[0]) + .5)
    ax.set_yticklabels(yticklabels, va="center")
    ax.set_aspect(1)

    for p, color, value in zip(img.get_paths(), img.get_facecolors(),
                               img.get_array()):
        x, y = p.vertices[:-2, :].mean(0)
        if np.mean(color[:3]) > 0.5:
            # pixel bright: use black for number
            c = 'k'
        else:
            c = 'w'
        ax.text(x, y, fmt.format(value), color=c, ha="center", va="center")
    return ax
