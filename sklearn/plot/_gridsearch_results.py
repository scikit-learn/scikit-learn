import numpy as np
from sklearn.plot import plot_heatmap


def _plot_1D_results(cv_results, params, ax, xlabel, ylabel, title,
                     fmt, xtickrotation):
    import matplotlib.pyplot as plt

    param = params[0]
    param_range = sorted(cv_results['param_%s' % param])
    train_scores_mean = cv_results['mean_train_score']
    train_scores_std = cv_results['std_train_score']
    test_scores_mean = cv_results['mean_test_score']
    test_scores_std = cv_results['std_test_score']

    lw = 2
    x_vales = range(len(param_range))
    ax.plot(x_vales, train_scores_mean,
            label="Training score",
            color="darkorange", lw=lw)
    ax.fill_between(x_vales, train_scores_mean - train_scores_std,
                    train_scores_mean + train_scores_std, alpha=0.2,
                    color="darkorange", lw=lw)

    img = ax.plot(x_vales, test_scores_mean,
                  label="Cross-validation score",
                  color="navy", lw=lw)
    ax.fill_between(x_vales, test_scores_mean - test_scores_std,
                    test_scores_mean + test_scores_std, alpha=0.2,
                    color="navy", lw=lw)
    ax.set_xticks(x_vales)
    ax.set_xticklabels([fmt.format(x) for x in param_range],
                       rotation=xtickrotation)

    xlabel = params[0] if xlabel is None else xlabel
    ax.set_xlabel(param)
    ax.legend()
    ylabel = "Score" if ylabel is None else ylabel
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    plt.draw()
    return img


def _plot_2D_results(cv_results, params, metric, ax, xlabel,
                     ylabel, title, cmap, vmin, vmax, fmt,
                     xtickrotation, norm):
    import matplotlib.pyplot as plt

    parameter1_values = np.unique(cv_results['param_%s' % params[0]])
    parameter2_values = np.unique(cv_results['param_%s' % params[1]])

    scores = cv_results[metric].reshape(len(parameter1_values),
                                        len(parameter2_values))

    xlabel = params[0] if xlabel is None else xlabel
    ylabel = params[1] if ylabel is None else ylabel

    cmap = cmap if cmap is not None else plt.cm.hot

    img = plot_heatmap(scores, xlabel=xlabel, ylabel=ylabel,
                       xticklabels=parameter1_values,
                       yticklabels=parameter2_values,
                       title=title, cmap=cmap,
                       vmin=vmin, vmax=vmax, fmt=fmt, ax=ax,
                       xtickrotation=xtickrotation, norm=norm)
    return img


def plot_gridsearch_results(cv_results, metric='mean_test_score',
                            xlabel=None, ylabel=None,
                            title='Grid Search Results', cmap=None,
                            vmin=None, vmax=None, ax=None, fmt="{:.2f}",
                            xtickrotation=45, norm=None):
    """Plot the grid search results as a line chart for 1D search and heatmap
    for a 2D search. This function will not work if grid-search has more than
    2 parameters in the search space.

    Parameters
    ----------
    cv_results : dict of numpy (masked) ndarrays
        The cv_results_ attribute of the GridSearchCV object.

    xlabel : string, optional (default=None)
        Label for the x-axis. If None, the first key of the param_grid will
        be used as the xlabel.

    ylabel : string, optional (default=None)
        Label for the y-axis. If None, the second key of the param_grid will
        be used as the ylabel.

    metric : string, optional (default="mean_test_score")
        The metric from the GridSearchCV results to display. This is ignored
        if only 1 parameter is used in grid search.

    title : string, optional (default="Grid Search Results")
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
        Format string to convert value to text.

    xtickrotation : float, optional (default=45)
        Rotation of the xticklabels.

    norm : matplotlib normalizer, optional (default=None)
        Normalizer passed to pcolormesh function from matplotlib used to
        generate the heatmap. This is ignored if only 1 parameter is used in
        grid search.
    """
    import matplotlib.pyplot as plt

    params = sorted(cv_results['params'][0].keys())
    nparams = len(params)

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    if nparams == 1:
        img = _plot_1D_results(cv_results, params, ax, xlabel, ylabel, title,
                               fmt, xtickrotation)

    elif nparams == 2:
        img = _plot_2D_results(cv_results, params, metric, ax, xlabel,
                               ylabel, title, cmap, vmin, vmax, fmt,
                               xtickrotation, norm)

    else:
        raise ValueError('Plot function supports upto 2 parameters in grid'
                         'search, got %d.' % nparams)

    return img
