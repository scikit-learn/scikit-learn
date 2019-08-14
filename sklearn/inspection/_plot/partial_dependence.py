from itertools import count
import warnings

import numpy as np

from ...utils import check_matplotlib_support


class PartialDependenceDisplay:
    """Partial Dependence visualization.

    If is recommended to use
    :func:`~sklearn.inspection.plot_partial_dependence` to create a
    :class:`~sklearn.inspection.PartialDependenceDisplay`. All parameters are
    stored as attributes.

    Read more in
    :ref:`sphx_glr_auto_examples_plot_roc_curve_visualization_api.py` and the
    :ref:`User Guide <visualizations>`.

        .. versionadded:: 0.22

    Parameters
    ----------
    pd_results : list of (ndarray, ndarray)
        Results of `sklearn.inspection.partial_dependence` for ``features``.

    features : list of {(int, ), (int, int)}
        Indicies of features for a given plot. A tuple of one int will plot
        a partial dependence curve of one feature. A tuple of two ints will
        plot a two-way partial dependence curve as a contour plot.

    feature_names : list of str
        Feature names corrsponding to the indicies in ``features``.

    target_idx : int
        - In a multiclass setting, specifies the class for which the PDPs
          should be computed. Note that for binary classification, the
          positive class (index 1) is always used.
        - In a multioutput setting, specifies the task for which the PDPs
          should be computed
        Ignored in binary classification or classical regression settings.

    pdp_lim : dict
        Global min and max average predictions. `pdp_lim[1]` is the global min
        and max for single partial dependence curves. `pdp_lim[2]` is the
        global min and max for two-way partial dependence curves.

    deciles : dict
        Deciles for feature indicies in ``features``.

    Attributes
    ----------
    bounding_ax_ : matplotlib Axes or None
        If `ax` is an axes or None, the `bounding_ax_` is the axes where the
        grid of partial dependence plots are drawn. If `ax` is  list of axes or
        a numpy array of axes, `bounding_ax_` is None.

    axes_ : ndarray of matplotlib Axes
        If `ax` is an axes or None, `axes_[i, j]` is the axes on the ith row
        and jth column. If `ax` is a list of axes, `axes_[i]` is the ith item
        in `ax`. Elements that are None corresponds to a nonexisting axes in
        that position.

    lines_ : ndarray of matplotlib Artists
        If `ax` is an axes or None, `line_[i, j]` is the partial dependence
        curve on the ith row and jth column. If `ax` is a list of axes,
        `lines_[i]` is the partial dependence curve corresponding to the ith
        item in `ax`. Elements that are None corresponds to a nonexisting axes
        or an axes that does not include a line plot.

    contours_ : ndarray of matplotlib Artists
        If `ax` is an axes or None, `contours_[i, j]` is the partial dependence
        plot on the ith row and jth column. If `ax` is a list of axes,
        `contours_[i]` is the partial dependence plot corresponding to the ith
        item in `ax`. Elements that are None corresponds to a nonexisting axes
        or an axes that does not include a contour plot.

    figure_ : matplotlib Figure
        Figure containing partial dependence plots.

    """
    def __init__(self, pd_results, features, feature_names, target_idx,
                 pdp_lim, deciles):
        self.pd_results = pd_results
        self.features = features
        self.feature_names = feature_names
        self.target_idx = target_idx
        self.pdp_lim = pdp_lim
        self.deciles = deciles

    def plot(self, ax=None, n_cols=3, line_kw=None, contour_kw=None, fig=None):
        """Plot partial dependence plots.

        Parameters
        ----------
        ax : Matplotlib axes, list of Matplotlib axes, ndarray of Matplotlib \
        Axes, default=None
            Axes to plot the partial dependence curves.
            - If a single axes is given, it is treated as a bounding axes and
              a grid of partial depdendence plots will be drawn on that top of
              it.
            - If a list of axes or a ndarray of axes are passed in, the partial
              dependence plots will be drawn on those axes.
            - By default, a single bounding axes is created and treated as the
              single axes case.

        n_cols : int, default=3
            The maximum number of columns in the grid plot.

        line_kw : dict or None, default=None
            Dict with keywords passed to the ``matplotlib.pyplot.plot`` call.
            For one-way partial dependence plots.

        contour_kw : dict or None, default=None
            Dict with keywords passed to the ``matplotlib.pyplot.contourf``
            call for two-way partial dependence plots.

        fig : Matplotlib figure object or None, default=None
            A figure object onto which the plots will be drawn, after the
            figure has been cleared. By default, a new one is created.

            .. deprecated:: 0.22

        Returns
        -------
        display: :class:`~sklearn.inspection.PartialDependenceDisplay`
        """

        check_matplotlib_support("plot_partial_dependence")
        import matplotlib.pyplot as plt  # noqa
        from matplotlib import transforms  # noqa
        from matplotlib.ticker import MaxNLocator  # noqa
        from matplotlib.ticker import ScalarFormatter  # noqa
        from matplotlib.gridspec import GridSpecFromSubplotSpec  # noqa

        if line_kw is None:
            line_kw = {}
        if contour_kw is None:
            contour_kw = {}

        contour_kw = {**{"alpha": 0.75}, **contour_kw}

        if fig is not None:
            warnings.warn("The fig parameter is deprecated in version "
                          "0.22 and will be removed in version 0.24",
                          DeprecationWarning)
            fig.clear()

        if ax is None:
            if fig is None:
                _, ax = plt.subplots()
            else:
                ax = fig.gca()

        n_features = len(self.features)

        if isinstance(ax, list):
            if len(ax) != n_features:
                raise ValueError("Expected len(ax) == len(features), "
                                 "got len(ax) = {}".format(len(ax)))
            self.bounding_ax_ = None
            self.figure_ = ax[0].figure
            self.axes_ = np.array(ax)
            self.lines_ = np.empty(n_features, dtype=np.object)
            self.contours_ = np.empty(n_features, dtype=np.object)

        elif isinstance(ax, np.ndarray):
            self.bounding_ax_ = None
            self.axes_ = ax
            self.figure_ = ax.ravel()[0].figure
            self.lines_ = np.empty_like(ax, dtype=np.object)
            self.contours_ = np.empty_like(ax, dtype=np.object)
        else:  # single axes
            self.bounding_ax_ = ax
            self.figure_ = ax.figure
            ax.set_axis_off()

            n_cols = min(n_cols, n_features)
            n_rows = int(np.ceil(n_features / float(n_cols)))

            self.axes_ = np.empty((n_rows, n_cols), dtype=np.object)
            self.lines_ = np.empty((n_rows, n_cols), dtype=np.object)
            self.contours_ = np.empty((n_rows, n_cols), dtype=np.object)

            axes_ravel = self.axes_.ravel()

            gs = GridSpecFromSubplotSpec(n_rows, n_cols,
                                         subplot_spec=ax.get_subplotspec())
            for i, spec in zip(range(n_features), gs):
                axes_ravel[i] = self.figure_.add_subplot(spec)

        # create contour levels for two-way plots
        if 2 in self.pdp_lim:
            Z_level = np.linspace(*self.pdp_lim[2], num=8)
        lines_ravel = self.lines_.ravel()
        contours_ravel = self.contours_.ravel()

        for i, axi, fx, (avg_preds, values) in zip(count(),
                                                   self.axes_.ravel(),
                                                   self.features,
                                                   self.pd_results):
            if len(values) == 1:
                lines_ravel[i] = axi.plot(values[0],
                                          avg_preds[self.target_idx].ravel(),
                                          **line_kw)[0]
            else:
                # contour plot
                XX, YY = np.meshgrid(values[0], values[1])
                Z = avg_preds[self.target_idx].T
                CS = axi.contour(XX, YY, Z, levels=Z_level, linewidths=0.5,
                                 colors='k')
                contours_ravel[i] = axi.contourf(XX, YY, Z, levels=Z_level,
                                                 vmax=Z_level[-1],
                                                 vmin=Z_level[0],
                                                 **contour_kw)
                axi.clabel(CS, fmt='%2.2f', colors='k', fontsize=10,
                           inline=True)

            trans = transforms.blended_transform_factory(axi.transData,
                                                         axi.transAxes)
            ylim = axi.get_ylim()
            axi.vlines(self.deciles[fx[0]], 0, 0.05, transform=trans,
                       color='k')
            axi.set_xlabel(self.feature_names[fx[0]])
            axi.set_ylim(ylim)

            if len(values) == 1:
                axi.set_ylabel('Partial dependence')
                axi.set_ylim(self.pdp_lim[1])
            else:
                # contour plot
                trans = transforms.blended_transform_factory(axi.transAxes,
                                                             axi.transData)
                xlim = axi.get_xlim()
                axi.hlines(self.deciles[fx[1]], 0, 0.05, transform=trans,
                           color='k')
                # hline erases xlim
                axi.set_ylabel(self.feature_names[fx[1]])
                axi.set_xlim(xlim)

        return self
