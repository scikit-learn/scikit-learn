import numpy as np
from sklearn.plot import plot_heatmap


def plot_gridsearch_results(cv_results, metric='mean_test_score',
                            xlabel=None, ylabel=None,
                            title='Grid Search Results', cmap=None,
                            vmin=None, vmax=None, ax=None, fmt="{:.2f}",
                            xtickrotation=45, norm=None):
    """Print and plot the confusion matrix as a heatmap. Normalization can be
    applied by setting `normalize=True`.

    Parameters
    ----------
    cv_results : dict of numpy (masked) ndarrays
        The cv_results_ attribute of the GridSearchCV object.

    xlabel : string, default=None
        Label for the x-axis. If None, the first key of the param_grid will
        be used as the xlabel.

    ylabel : string, default=None
        Label for the y-axis. If None, the second key of the param_grid will
        be used as the ylabel.

    metric : string, default="mean_test_score"
        The metric from the GridSearchCV results to display.

    title : string, default="Grid Search Results"
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
        Format string to convert value to text.

    xtickrotation : float, default=45
        Rotation of the xticklabels.

    norm : matplotlib normalizer
        Normalizer passed to pcolormesh function from matplotlib used to
        generate the heatmap.
    """

    import matplotlib.pyplot as plt

    params = sorted(cv_results['params'][0].keys())

    if len(params) == 1:
        # plot a line chart
        pass
    elif len(params) == 2:
        parameter1_values = np.unique(cv_results['param_%s' % params[0]])
        parameter2_values = np.unique(cv_results['param_%s' % params[1]])

        scores = cv_results[metric].reshape(len(parameter1_values),
                                            len(parameter2_values))

        xlabel = params[0] if xlabel is None else xlabel
        ylabel = params[1] if ylabel is None else ylabel

        plot_heatmap(scores, cmap=plt.cm.hot, xlabel=xlabel, ylabel=ylabel,
                     xticklabels=parameter1_values,
                     yticklabels=parameter2_values,
                     ax=ax, norm=norm)

        plt.title(title)
        plt.show()
    else:
        # print error statement
        pass
