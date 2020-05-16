import numpy as np

from ...utils import check_matplotlib_support


def plot_gridsearch_results(cv_results, params=None, metric='mean_test_score',
                            xlabel=None, ylabel=None,
                            title='Grid Search Results', cmap=None,
                            vmin=None, vmax=None, ax=None, fmt="{:.2f}",
                            xtickrotation=45, norm=None):
    """Plot grid search results.

    The results are plotted as a line chart for 1D search and as a heatmap
    for a 2D search. This function will not work if grid-search has more than
    2 parameters in the search space.

    Parameters
    ----------
    cv_results : dict of numpy (masked) ndarrays
        The cv_results_ attribute of the GridSearchCV object.

    params : string or list of strings (default=None)
        The parameters from GridSearchCV results to display.
        If None, all parameters will be plot. Currently, there is support only
        for 1D and 2D parameters visualization.

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
    cv_params = list(cv_results['params'][0].keys())

    if params is None:
        if len(cv_params) > 2:
            raise ValueError("Plot function supports upto 2 parameters in grid"
                             "search, got %d." % len(cv_params))
        params = cv_params
    assert params is not None

    if len(params) == 0:
        raise ValueError("Expecting a non-null set of params. "
                         "Nothing to do with empty set of params")

    if not isinstance(params, (list, tuple)):
        assert isinstance(params, str)
        params = [params]

    if len(params) > 2:
        raise ValueError("Plot function supports upto 2 parameters in grid"
                         "search, got %d." % len(params))

    if not all(p in cv_params for p in params):
        msg = "Expected to have 'params' from {!s}, instead got {!s}"
        raise ValueError(msg.format(cv_params, params))

    gs = GridSearchDisplay(cv_results, params)
    return gs.plot(metric=metric, xlabel=xlabel, ylabel=ylabel,
                   title=title, cmap=cmap,
                   vmin=vmin, vmax=vmax, ax=ax, fmt=fmt,
                   xtickrotation=xtickrotation, norm=norm)


class GridSearchDisplay:
    """
    GridSearch display
    """
    def __init__(self, cv_results, params):
        self.cv_results = cv_results
        self.params = params

    def plot(self, metric='mean_test_score',
             xlabel=None, ylabel=None,
             title='Grid Search Results', cmap=None,
             vmin=None, vmax=None, ax=None, fmt="{:.2f}",
             xtickrotation=45, norm=None):

        check_matplotlib_support("GridSearchDisplay.plot")
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        img = _plot_gridsearch_results(
            self.cv_results, params=self.params,
            metric=metric, xlabel=xlabel, ylabel=ylabel,
            title=title, cmap=cmap,
            vmin=vmin, vmax=vmax, ax=ax, fmt=fmt,
            xtickrotation=xtickrotation, norm=norm)

        self.im_ = img

        # fig.colorbar(self.im_, ax=ax)
        # ax.set(xticks=np.arange(n_classes),
        #        yticks=np.arange(n_classes),
        #        xticklabels=self.display_labels,
        #        yticklabels=self.display_labels,
        #        ylabel="True label",
        #        xlabel="Predicted label")

        # ax.set_ylim((n_classes - 0.5, -0.5))
        plt.setp(ax.get_xticklabels(), rotation=xtickrotation)

        self.figure_ = fig
        self.ax_ = ax

        return self


def _plot_gridsearch_results(cv_results, params,
                             metric='mean_test_score',
                             xlabel=None, ylabel=None,
                             title='Grid Search Results', cmap=None,
                             vmin=None, vmax=None, ax=None, fmt="{:.2f}",
                             xtickrotation=45, norm=None):
    import matplotlib.pyplot as plt

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


def _plot_1D_results(cv_results, params, ax, xlabel, ylabel, title,
                     fmt, xtickrotation):
    param = params[0]
    param_range = sorted(cv_results['param_%s' % param])
    test_scores_mean = cv_results['mean_test_score']
    test_scores_std = cv_results['std_test_score']

    lw = 2
    x_vales = range(len(param_range))

    if 'mean_train_score' in cv_results and 'std_train_score' in cv_results:
        train_scores_mean = cv_results['mean_train_score']
        train_scores_std = cv_results['std_train_score']
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

    img = _plot_heatmap(scores, xlabel=xlabel, ylabel=ylabel,
                        xticklabels=parameter1_values,
                        yticklabels=parameter2_values,
                        title=title, cmap=cmap,
                        vmin=vmin, vmax=vmax, fmt=fmt, ax=ax,
                        xtickrotation=xtickrotation, norm=norm)
    return img


def _plot_heatmap(values, xlabel="", ylabel="", xticklabels=None,
                  yticklabels=None, title=None, cmap=None, vmin=None,
                  vmax=None, ax=None, fmt="{:.2f}", xtickrotation=45,
                  norm=None):
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

    title : string or None, default=None
        Title of the chart

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

    xtickrotation : float, default=45
        Rotation of the xticklabels.

    norm : matplotlib normalizer
        Normalizer passed to pcolor
    """

    import matplotlib.pyplot as plt
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    img = ax.pcolormesh(values, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)

    # this will allow us to access the pixel values:
    img.update_scalarmappable()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ax.set_xlim(0, values.shape[1])
    ax.set_ylim(0, values.shape[0])

    if xticklabels is None:
        xticklabels = [""] * values.shape[1]
    if yticklabels is None:
        yticklabels = [""] * values.shape[0]

    # +.5 makes the ticks centered on the pixels
    ax.set_xticks(np.arange(values.shape[1]) + .5)
    ax.set_xticklabels(xticklabels, ha="center", rotation=xtickrotation)
    ax.set_yticks(np.arange(values.shape[0]) + .5)
    ax.set_yticklabels(yticklabels, va="center")
    ax.set_aspect(1)

    for p, color, value in zip(img.get_paths(), img.get_facecolors(),
                               img.get_array()):
        x, y = p.vertices[:-2, :].mean(0)

        # adjusting x and y for alignment:
        x = x - 1./6
        y = y + 1./6

        if np.mean(color[:3]) > 0.5:
            # pixel bright: use black for number
            c = 'k'
        else:
            c = 'w'
        ax.text(x, y, fmt.format(value), color=c, ha="center", va="center")

    # Invert the y-axis so that the matrix looks like a diagonal matrix and
    # not anti-diagonal matrix
    ax.invert_yaxis()

    # set title if not none:
    if title is not None:
        ax.set_title(title)

    ax.figure.colorbar(img, ax=ax)

    return img
