"""Partial dependence display"""
from itertools import count
from typing import Iterable

import numpy as np

from ..utils import check_matplotlib_support


class PartialDependenceDisplay:
    def __init__(self, pd_results, features, feature_names, target_idx,
                 pdp_lim, deciles):
        self.pd_results = pd_results
        self.features = features
        self.feature_names = feature_names
        self.target_idx = target_idx
        self.pdp_lim = pdp_lim
        self.deciles = deciles

    def plot(self, ax=None, n_cols=3, line_kw=None, contour_kw=None,
             ax_set_kw=None):
        check_matplotlib_support("plot_partial_dependence")
        import matplotlib.pyplot as plt  # noqa
        from matplotlib import transforms  # noqa
        from matplotlib.ticker import MaxNLocator  # noqa
        from matplotlib.ticker import ScalarFormatter  # noqa
        from matplotlib.gridspec import GridSpecFromSubplotSpec  # noqa

        # create contour levels for two-way plots
        if 2 in self.pdp_lim:
            Z_level = np.linspace(*self.pdp_lim[2], num=8)

        if ax is None:
            _, ax = plt.subplots()

        # check axes is the correct length

        if line_kw is None:
            line_kw = {}
        if contour_kw is None:
            contour_kw = {}
        if ax_set_kw is None:
            ax_set_kw = {}

        if isinstance(ax, Iterable):
            self.bounding_ax_ = None
            self.axes_ = dict(enumerate(ax))
            self.figure_ = ax[0].figure
            self.lines_ = {}
            self.contours_ = {}

            def _get_axes(i):
                return ax[i]

            def _store_line(line, i):
                self.lines_[i] = line

            def _store_contour(contour, i):
                self.contours_[i] = contour
        else:
            self.bounding_ax_ = ax
            self.axes_ = {}
            self.figure_ = ax.figure
            self.lines_ = {}
            self.contours_ = {}

            ax.set_axis_off()
            n_cols = min(n_cols, len(self.features))
            n_rows = int(np.ceil(len(self.features) / float(n_cols)))

            gs = GridSpecFromSubplotSpec(n_rows, n_cols,
                                         subplot_spec=ax.get_subplotspec())

            def _get_axes(i):
                row_idx, col_idx = i // n_cols, i % n_cols
                axi = self.figure_.add_subplot(gs[row_idx, col_idx])
                self.axes_[row_idx, col_idx] = axi
                return axi

            def _store_line(line, i):
                row_idx, col_idx = i // n_cols, i % n_cols
                self.lines_[row_idx, col_idx] = line

            def _store_contour(contour, i):
                row_idx, col_idx = i // n_cols, i % n_cols
                self.contours_[row_idx, col_idx] = contour

        for i, fx, name, (avg_preds, values) in zip(count(), self.features,
                                                    self.feature_names,
                                                    self.pd_results):
            axi = _get_axes(i)

            if len(values) == 1:
                line = axi.plot(values[0], avg_preds[self.target_idx].ravel(),
                                **line_kw)[0]
                _store_line(line, i)
            else:
                # contour plot
                assert len(values) == 2
                XX, YY = np.meshgrid(values[0], values[1])
                Z = avg_preds[self.target_idx].T
                CS = axi.contour(XX, YY, Z, levels=Z_level, linewidths=0.5,
                                 colors='k')
                contour = axi.contourf(XX, YY, Z, levels=Z_level,
                                       vmax=Z_level[-1], vmin=Z_level[0],
                                       alpha=0.75, **contour_kw)
                _store_contour(contour, i)
                axi.clabel(CS, fmt='%2.2f', colors='k', fontsize=10,
                           inline=True)

            trans = transforms.blended_transform_factory(axi.transData,
                                                         axi.transAxes)
            ylim = axi.get_ylim()
            axi.vlines(self.deciles[fx[0]], [0], 0.05, transform=trans,
                       color='k')
            axi.set_xlabel(name[0])
            axi.set_ylim(ylim)

            if len(values) == 1:
                axi.set_ylabel('Partial dependence')
                axi.set_ylim(self.pdp_lim[1])
            else:
                trans = transforms.blended_transform_factory(axi.transAxes,
                                                             axi.transData)
                xlim = axi.get_xlim()
                axi.hlines(self.deciles[fx[1]], [0], 0.05, transform=trans,
                           color='k')
                # hline erases xlim
                axi.set_ylabel(name[1])
                axi.set_xlim(xlim)

            axi.set(**ax_set_kw)

        return self
