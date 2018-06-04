"""
:mod:`~matplotlib.gridspec` is a module which specifies the location
of the subplot in the figure.

    `GridSpec`
        specifies the geometry of the grid that a subplot will be
        placed. The number of rows and number of columns of the grid
        need to be set. Optionally, the subplot layout parameters
        (e.g., left, right, etc.) can be tuned.

    `SubplotSpec`
        specifies the location of the subplot in the given `GridSpec`.

"""

from __future__ import absolute_import, division, print_function

import six

import copy
import logging
import warnings

import numpy as np

import matplotlib as mpl
from matplotlib import _pylab_helpers, tight_layout, rcParams
from matplotlib.transforms import Bbox
import matplotlib._layoutbox as layoutbox
from matplotlib.cbook import mplDeprecation

_log = logging.getLogger(__name__)


class GridSpecBase(object):
    """
    A base class of GridSpec that specifies the geometry of the grid
    that a subplot will be placed.
    """

    def __init__(self, nrows, ncols, height_ratios=None, width_ratios=None):
        """
        The number of rows and number of columns of the grid need to
        be set. Optionally, the ratio of heights and widths of rows and
        columns can be specified.
        """
        self._nrows, self._ncols = nrows, ncols
        self.set_height_ratios(height_ratios)
        self.set_width_ratios(width_ratios)

    def get_geometry(self):
        'get the geometry of the grid, e.g., 2,3'
        return self._nrows, self._ncols

    def get_subplot_params(self, figure=None, fig=None):
        pass

    def new_subplotspec(self, loc, rowspan=1, colspan=1):
        """
        create and return a SuplotSpec instance.
        """
        loc1, loc2 = loc
        subplotspec = self[loc1:loc1+rowspan, loc2:loc2+colspan]
        return subplotspec

    def set_width_ratios(self, width_ratios):
        if width_ratios is not None and len(width_ratios) != self._ncols:
            raise ValueError('Expected the given number of width ratios to '
                             'match the number of columns of the grid')
        self._col_width_ratios = width_ratios

    def get_width_ratios(self):
        return self._col_width_ratios

    def set_height_ratios(self, height_ratios):
        if height_ratios is not None and len(height_ratios) != self._nrows:
            raise ValueError('Expected the given number of height ratios to '
                             'match the number of rows of the grid')
        self._row_height_ratios = height_ratios

    def get_height_ratios(self):
        return self._row_height_ratios

    def get_grid_positions(self, fig, raw=False):
        """
        return lists of bottom and top position of rows, left and
        right positions of columns.

        If raw=True, then these are all in units relative to the container
        with no margins.  (used for constrained_layout).
        """
        nrows, ncols = self.get_geometry()

        if raw:
            left = 0.
            right = 1.
            bottom = 0.
            top = 1.
            wspace = 0.
            hspace = 0.
        else:
            subplot_params = self.get_subplot_params(fig)
            left = subplot_params.left
            right = subplot_params.right
            bottom = subplot_params.bottom
            top = subplot_params.top
            wspace = subplot_params.wspace
            hspace = subplot_params.hspace
        tot_width = right - left
        tot_height = top - bottom

        # calculate accumulated heights of columns
        cell_h = tot_height / (nrows + hspace*(nrows-1))
        sep_h = hspace * cell_h
        if self._row_height_ratios is not None:
            norm = cell_h * nrows / sum(self._row_height_ratios)
            cell_heights = [r * norm for r in self._row_height_ratios]
        else:
            cell_heights = [cell_h] * nrows
        sep_heights = [0] + ([sep_h] * (nrows-1))
        cell_hs = np.cumsum(np.column_stack([sep_heights, cell_heights]).flat)

        # calculate accumulated widths of rows
        cell_w = tot_width / (ncols + wspace*(ncols-1))
        sep_w = wspace * cell_w
        if self._col_width_ratios is not None:
            norm = cell_w * ncols / sum(self._col_width_ratios)
            cell_widths = [r * norm for r in self._col_width_ratios]
        else:
            cell_widths = [cell_w] * ncols
        sep_widths = [0] + ([sep_w] * (ncols-1))
        cell_ws = np.cumsum(np.column_stack([sep_widths, cell_widths]).flat)

        fig_tops, fig_bottoms = (top - cell_hs).reshape((-1, 2)).T
        fig_lefts, fig_rights = (left + cell_ws).reshape((-1, 2)).T
        return fig_bottoms, fig_tops, fig_lefts, fig_rights

    def __getitem__(self, key):
        """Create and return a SuplotSpec instance.
        """
        nrows, ncols = self.get_geometry()

        def _normalize(key, size):  # Includes last index.
            if isinstance(key, slice):
                start, stop, _ = key.indices(size)
                if stop > start:
                    return start, stop - 1
            else:
                if key < 0:
                    key += size
                if 0 <= key < size:
                    return key, key
            raise IndexError("invalid index")

        if isinstance(key, tuple):
            try:
                k1, k2 = key
            except ValueError:
                raise ValueError("unrecognized subplot spec")
            num1, num2 = np.ravel_multi_index(
                [_normalize(k1, nrows), _normalize(k2, ncols)], (nrows, ncols))
        else:  # Single key
            num1, num2 = _normalize(key, nrows * ncols)

        return SubplotSpec(self, num1, num2)


class GridSpec(GridSpecBase):
    """
    A class that specifies the geometry of the grid that a subplot
    will be placed. The location of grid is determined by similar way
    as the SubplotParams.
    """

    def __init__(self, nrows, ncols, figure=None,
                 left=None, bottom=None, right=None, top=None,
                 wspace=None, hspace=None,
                 width_ratios=None, height_ratios=None):
        """
        The number of rows and number of columns of the grid need to be set.
        Optionally, the subplot layout parameters (e.g., left, right, etc.)
        can be tuned.

        Parameters
        ----------
        nrows : int
            Number of rows in grid.

        ncols : int
            Number or columns in grid.

        Notes
        -----
        See `~.figure.SubplotParams` for descriptions of the layout parameters.
        """
        self.left = left
        self.bottom = bottom
        self.right = right
        self.top = top
        self.wspace = wspace
        self.hspace = hspace
        self.figure = figure

        GridSpecBase.__init__(self, nrows, ncols,
                              width_ratios=width_ratios,
                              height_ratios=height_ratios)

        if (self.figure is None) or not self.figure.get_constrained_layout():
            self._layoutbox = None
        else:
            self.figure.init_layoutbox()
            self._layoutbox = layoutbox.LayoutBox(
                parent=self.figure._layoutbox,
                name='gridspec' + layoutbox.seq_id(),
                artist=self)
        # by default the layoutbox for a gridsepc will fill a figure.
        # but this can change below if the gridspec is created from a
        # subplotspec. (GridSpecFromSubplotSpec)

    _AllowedKeys = ["left", "bottom", "right", "top", "wspace", "hspace"]

    def __getstate__(self):
        state = self.__dict__
        try:
            state.pop('_layoutbox')
        except KeyError:
            pass
        return state

    def __setstate__(self, state):
        self.__dict__ = state
        # layoutboxes don't survive pickling...
        self._layoutbox = None

    def update(self, **kwargs):
        """
        Update the current values.  If any kwarg is None, default to
        the current value, if set, otherwise to rc.
        """

        for k, v in six.iteritems(kwargs):
            if k in self._AllowedKeys:
                setattr(self, k, v)
            else:
                raise AttributeError("%s is unknown keyword" % (k,))

        for figmanager in six.itervalues(_pylab_helpers.Gcf.figs):
            for ax in figmanager.canvas.figure.axes:
                # copied from Figure.subplots_adjust
                if not isinstance(ax, mpl.axes.SubplotBase):
                    # Check if sharing a subplots axis
                    if isinstance(ax._sharex, mpl.axes.SubplotBase):
                        if ax._sharex.get_subplotspec().get_gridspec() == self:
                            ax._sharex.update_params()
                            ax._set_position(ax._sharex.figbox)
                    elif isinstance(ax._sharey, mpl.axes.SubplotBase):
                        if ax._sharey.get_subplotspec().get_gridspec() == self:
                            ax._sharey.update_params()
                            ax._set_position(ax._sharey.figbox)
                else:
                    ss = ax.get_subplotspec().get_topmost_subplotspec()
                    if ss.get_gridspec() == self:
                        ax.update_params()
                        ax._set_position(ax.figbox)

    def get_subplot_params(self, figure=None, fig=None):
        """
        Return a dictionary of subplot layout parameters. The default
        parameters are from rcParams unless a figure attribute is set.
        """
        if fig is not None:
            warnings.warn("the 'fig' kwarg is deprecated "
                          "use 'figure' instead", mplDeprecation)
        if figure is None:
            figure = fig

        if figure is None:
            kw = {k: rcParams["figure.subplot."+k] for k in self._AllowedKeys}
            subplotpars = mpl.figure.SubplotParams(**kw)
        else:
            subplotpars = copy.copy(figure.subplotpars)

        update_kw = {k: getattr(self, k) for k in self._AllowedKeys}
        subplotpars.update(**update_kw)

        return subplotpars

    def locally_modified_subplot_params(self):
        return [k for k in self._AllowedKeys if getattr(self, k)]

    def tight_layout(self, figure, renderer=None,
                     pad=1.08, h_pad=None, w_pad=None, rect=None):
        """
        Adjust subplot parameters to give specified padding.

        Parameters
        ----------

        pad : float
            Padding between the figure edge and the edges of subplots, as a
            fraction of the font-size.
        h_pad, w_pad : float, optional
            Padding (height/width) between edges of adjacent subplots.
            Defaults to ``pad_inches``.
        rect : tuple of 4 floats, optional
            (left, bottom, right, top) rectangle in normalized figure
            coordinates that the whole subplots area (including labels) will
            fit into.  Default is (0, 0, 1, 1).
        """

        subplotspec_list = tight_layout.get_subplotspec_list(
            figure.axes, grid_spec=self)
        if None in subplotspec_list:
            warnings.warn("This figure includes Axes that are not compatible "
                          "with tight_layout, so results might be incorrect.")

        if renderer is None:
            renderer = tight_layout.get_renderer(figure)

        kwargs = tight_layout.get_tight_layout_figure(
            figure, figure.axes, subplotspec_list, renderer,
            pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
        self.update(**kwargs)


class GridSpecFromSubplotSpec(GridSpecBase):
    """
    GridSpec whose subplot layout parameters are inherited from the
    location specified by a given SubplotSpec.
    """
    def __init__(self, nrows, ncols,
                 subplot_spec,
                 wspace=None, hspace=None,
                 height_ratios=None, width_ratios=None):
        """
        The number of rows and number of columns of the grid need to
        be set. An instance of SubplotSpec is also needed to be set
        from which the layout parameters will be inherited. The wspace
        and hspace of the layout can be optionally specified or the
        default values (from the figure or rcParams) will be used.
        """
        self._wspace = wspace
        self._hspace = hspace
        self._subplot_spec = subplot_spec
        GridSpecBase.__init__(self, nrows, ncols,
                              width_ratios=width_ratios,
                              height_ratios=height_ratios)
        # do the layoutboxes
        subspeclb = subplot_spec._layoutbox
        if subspeclb is None:
            self._layoutbox = None
        else:
            # OK, this is needed to divide the figure.
            self._layoutbox = subspeclb.layout_from_subplotspec(
                    subplot_spec,
                    name=subspeclb.name + '.gridspec' + layoutbox.seq_id(),
                    artist=self)

    def get_subplot_params(self, figure=None, fig=None):
        """Return a dictionary of subplot layout parameters.
        """
        if fig is not None:
            warnings.warn("the 'fig' kwarg is deprecated "
                          "use 'figure' instead", mplDeprecation)
        if figure is None:
            figure = fig

        hspace = (self._hspace if self._hspace is not None
                  else figure.subplotpars.hspace if figure is not None
                  else rcParams["figure.subplot.hspace"])
        wspace = (self._wspace if self._wspace is not None
                  else figure.subplotpars.wspace if figure is not None
                  else rcParams["figure.subplot.wspace"])

        figbox = self._subplot_spec.get_position(figure)
        left, bottom, right, top = figbox.extents

        return mpl.figure.SubplotParams(left=left, right=right,
                                        bottom=bottom, top=top,
                                        wspace=wspace, hspace=hspace)

    def get_topmost_subplotspec(self):
        """Get the topmost SubplotSpec instance associated with the subplot."""
        return self._subplot_spec.get_topmost_subplotspec()


class SubplotSpec(object):
    """Specifies the location of the subplot in the given `GridSpec`.
    """

    def __init__(self, gridspec, num1, num2=None):
        """
        The subplot will occupy the num1-th cell of the given
        gridspec.  If num2 is provided, the subplot will span between
        num1-th cell and num2-th cell.

        The index starts from 0.
        """
        self._gridspec = gridspec
        self.num1 = num1
        self.num2 = num2
        if gridspec._layoutbox is not None:
            glb = gridspec._layoutbox
            # So note that here we don't assign any layout yet,
            # just make the layoutbox that will conatin all items
            # associated w/ this axis.  This can include other axes like
            # a colorbar or a legend.
            self._layoutbox = layoutbox.LayoutBox(
                    parent=glb,
                    name=glb.name + '.ss' + layoutbox.seq_id(),
                    artist=self)
        else:
            self._layoutbox = None

    def __getstate__(self):
        state = self.__dict__
        try:
            state.pop('_layoutbox')
        except KeyError:
            pass
        return state

    def __setstate__(self, state):
        self.__dict__ = state
        # layoutboxes don't survive pickling...
        self._layoutbox = None

    def get_gridspec(self):
        return self._gridspec

    def get_geometry(self):
        """
        Get the subplot geometry (``n_rows, n_cols, start, stop``).

        start and stop are the index of the start and stop of the
        subplot.
        """
        rows, cols = self.get_gridspec().get_geometry()
        return rows, cols, self.num1, self.num2

    def get_rows_columns(self):
        """
        Get the subplot row and column numbers:
        (``n_rows, n_cols, row_start, row_stop, col_start, col_stop``)
        """
        gridspec = self.get_gridspec()
        nrows, ncols = gridspec.get_geometry()
        row_start, col_start = divmod(self.num1, ncols)
        if self.num2 is not None:
            row_stop, col_stop = divmod(self.num2, ncols)
        else:
            row_stop = row_start
            col_stop = col_start
        return nrows, ncols, row_start, row_stop, col_start, col_stop

    def get_position(self, figure, return_all=False):
        """Update the subplot position from ``figure.subplotpars``.
        """
        gridspec = self.get_gridspec()
        nrows, ncols = gridspec.get_geometry()
        rows, cols = np.unravel_index(
            [self.num1] if self.num2 is None else [self.num1, self.num2],
            (nrows, ncols))
        fig_bottoms, fig_tops, fig_lefts, fig_rights = \
            gridspec.get_grid_positions(figure)

        fig_bottom = fig_bottoms[rows].min()
        fig_top = fig_tops[rows].max()
        fig_left = fig_lefts[cols].min()
        fig_right = fig_rights[cols].max()
        figbox = Bbox.from_extents(fig_left, fig_bottom, fig_right, fig_top)

        if return_all:
            return figbox, rows[0], cols[0], nrows, ncols
        else:
            return figbox

    def get_topmost_subplotspec(self):
        'get the topmost SubplotSpec instance associated with the subplot'
        gridspec = self.get_gridspec()
        if hasattr(gridspec, "get_topmost_subplotspec"):
            return gridspec.get_topmost_subplotspec()
        else:
            return self

    def __eq__(self, other):
        # other may not even have the attributes we are checking.
        return ((self._gridspec, self.num1, self.num2)
                == (getattr(other, "_gridspec", object()),
                    getattr(other, "num1", object()),
                    getattr(other, "num2", object())))

    if six.PY2:
        def __ne__(self, other):
            return not self == other

    def __hash__(self):
        return hash((self._gridspec, self.num1, self.num2))
