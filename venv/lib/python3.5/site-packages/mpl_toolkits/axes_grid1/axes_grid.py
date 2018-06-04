from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import matplotlib.axes as maxes
import matplotlib.cbook as cbook
import matplotlib.ticker as ticker
from matplotlib.gridspec import SubplotSpec

from .axes_divider import Size, SubplotDivider, LocatableAxes, Divider
from .colorbar import Colorbar


def _extend_axes_pad(value):
    # Check whether a list/tuple/array or scalar has been passed
    ret = value
    if not hasattr(ret, "__getitem__"):
        ret = (value, value)
    return ret


def _tick_only(ax, bottom_on, left_on):
    bottom_off = not bottom_on
    left_off = not left_on
    # [l.set_visible(bottom_off) for l in ax.get_xticklabels()]
    # [l.set_visible(left_off) for l in ax.get_yticklabels()]
    # ax.xaxis.label.set_visible(bottom_off)
    # ax.yaxis.label.set_visible(left_off)
    ax.axis["bottom"].toggle(ticklabels=bottom_off, label=bottom_off)
    ax.axis["left"].toggle(ticklabels=left_off, label=left_off)


class CbarAxesBase(object):

    def colorbar(self, mappable, **kwargs):
        locator = kwargs.pop("locator", None)

        if locator is None:
            if "ticks" not in kwargs:
                kwargs["ticks"] = ticker.MaxNLocator(5)
        if locator is not None:
            if "ticks" in kwargs:
                raise ValueError("Either *locator* or *ticks* need" +
                                 " to be given, not both")
            else:
                kwargs["ticks"] = locator

        self._hold = True
        if self.orientation in ["top", "bottom"]:
            orientation = "horizontal"
        else:
            orientation = "vertical"

        cb = Colorbar(self, mappable, orientation=orientation, **kwargs)
        self._config_axes()

        def on_changed(m):
            cb.set_cmap(m.get_cmap())
            cb.set_clim(m.get_clim())
            cb.update_bruteforce(m)

        self.cbid = mappable.callbacksSM.connect('changed', on_changed)
        mappable.colorbar = cb

        self.locator = cb.cbar_axis.get_major_locator()

        return cb

    def _config_axes(self):
        '''
        Make an axes patch and outline.
        '''
        ax = self
        ax.set_navigate(False)

        ax.axis[:].toggle(all=False)
        b = self._default_label_on
        ax.axis[self.orientation].toggle(all=b)

        # for axis in ax.axis.values():
        #     axis.major_ticks.set_visible(False)
        #     axis.minor_ticks.set_visible(False)
        #     axis.major_ticklabels.set_visible(False)
        #     axis.minor_ticklabels.set_visible(False)
        #     axis.label.set_visible(False)

        # axis = ax.axis[self.orientation]
        # axis.major_ticks.set_visible(True)
        # axis.minor_ticks.set_visible(True)

        #axis.major_ticklabels.set_size(
        #    int(axis.major_ticklabels.get_size()*.9))
        #axis.major_tick_pad = 3

        # axis.major_ticklabels.set_visible(b)
        # axis.minor_ticklabels.set_visible(b)
        # axis.label.set_visible(b)

    def toggle_label(self, b):
        self._default_label_on = b
        axis = self.axis[self.orientation]
        axis.toggle(ticklabels=b, label=b)
        #axis.major_ticklabels.set_visible(b)
        #axis.minor_ticklabels.set_visible(b)
        #axis.label.set_visible(b)


class CbarAxes(CbarAxesBase, LocatableAxes):
    def __init__(self, *kl, **kwargs):
        orientation = kwargs.pop("orientation", None)
        if orientation is None:
            raise ValueError("orientation must be specified")
        self.orientation = orientation
        self._default_label_on = True
        self.locator = None

        super(LocatableAxes, self).__init__(*kl, **kwargs)

    def cla(self):
        super(LocatableAxes, self).cla()
        self._config_axes()


class Grid(object):
    """
    A class that creates a grid of Axes. In matplotlib, the axes
    location (and size) is specified in the normalized figure
    coordinates. This may not be ideal for images that needs to be
    displayed with a given aspect ratio.  For example, displaying
    images of a same size with some fixed padding between them cannot
    be easily done in matplotlib. AxesGrid is used in such case.
    """

    _defaultLocatableAxesClass = LocatableAxes

    def __init__(self, fig,
                 rect,
                 nrows_ncols,
                 ngrids=None,
                 direction="row",
                 axes_pad=0.02,
                 add_all=True,
                 share_all=False,
                 share_x=True,
                 share_y=True,
                 #aspect=True,
                 label_mode="L",
                 axes_class=None,
                 ):
        """
        Build an :class:`Grid` instance with a grid nrows*ncols
        :class:`~matplotlib.axes.Axes` in
        :class:`~matplotlib.figure.Figure` *fig* with
        *rect=[left, bottom, width, height]* (in
        :class:`~matplotlib.figure.Figure` coordinates) or
        the subplot position code (e.g., "121").

        Optional keyword arguments:

          ================  ========  =========================================
          Keyword           Default   Description
          ================  ========  =========================================
          direction         "row"     [ "row" | "column" ]
          axes_pad          0.02      float| pad between axes given in inches
                                      or tuple-like of floats,
                                      (horizontal padding, vertical padding)
          add_all           True      bool
          share_all         False     bool
          share_x           True      bool
          share_y           True      bool
          label_mode        "L"       [ "L" | "1" | "all" ]
          axes_class        None      a type object which must be a subclass
                                      of :class:`~matplotlib.axes.Axes`
          ================  ========  =========================================
        """
        self._nrows, self._ncols = nrows_ncols

        if ngrids is None:
            ngrids = self._nrows * self._ncols
        else:
            if (ngrids > self._nrows * self._ncols) or (ngrids <= 0):
                raise Exception("")

        self.ngrids = ngrids

        self._init_axes_pad(axes_pad)

        if direction not in ["column", "row"]:
            raise Exception("")

        self._direction = direction

        if axes_class is None:
            axes_class = self._defaultLocatableAxesClass
            axes_class_args = {}
        else:
            if (type(axes_class)) == type and \
                   issubclass(axes_class,
                              self._defaultLocatableAxesClass.Axes):
                axes_class_args = {}
            else:
                axes_class, axes_class_args = axes_class

        self.axes_all = []
        self.axes_column = [[] for _ in range(self._ncols)]
        self.axes_row = [[] for _ in range(self._nrows)]

        h = []
        v = []
        if isinstance(rect, six.string_types) or cbook.is_numlike(rect):
            self._divider = SubplotDivider(fig, rect, horizontal=h, vertical=v,
                                           aspect=False)
        elif isinstance(rect, SubplotSpec):
            self._divider = SubplotDivider(fig, rect, horizontal=h, vertical=v,
                                           aspect=False)
        elif len(rect) == 3:
            kw = dict(horizontal=h, vertical=v, aspect=False)
            self._divider = SubplotDivider(fig, *rect, **kw)
        elif len(rect) == 4:
            self._divider = Divider(fig, rect, horizontal=h, vertical=v,
                                    aspect=False)
        else:
            raise Exception("")

        rect = self._divider.get_position()

        # reference axes
        self._column_refax = [None for _ in range(self._ncols)]
        self._row_refax = [None for _ in range(self._nrows)]
        self._refax = None

        for i in range(self.ngrids):

            col, row = self._get_col_row(i)

            if share_all:
                sharex = self._refax
                sharey = self._refax
            else:
                if share_x:
                    sharex = self._column_refax[col]
                else:
                    sharex = None

                if share_y:
                    sharey = self._row_refax[row]
                else:
                    sharey = None

            ax = axes_class(fig, rect, sharex=sharex, sharey=sharey,
                            **axes_class_args)

            if share_all:
                if self._refax is None:
                    self._refax = ax
            else:
                if sharex is None:
                    self._column_refax[col] = ax
                if sharey is None:
                    self._row_refax[row] = ax

            self.axes_all.append(ax)
            self.axes_column[col].append(ax)
            self.axes_row[row].append(ax)

        self.axes_llc = self.axes_column[0][-1]

        self._update_locators()

        if add_all:
            for ax in self.axes_all:
                fig.add_axes(ax)

        self.set_label_mode(label_mode)

    def _init_axes_pad(self, axes_pad):
        axes_pad = _extend_axes_pad(axes_pad)
        self._axes_pad = axes_pad

        self._horiz_pad_size = Size.Fixed(axes_pad[0])
        self._vert_pad_size = Size.Fixed(axes_pad[1])

    def _update_locators(self):

        h = []

        h_ax_pos = []

        for _ in self._column_refax:
            #if h: h.append(Size.Fixed(self._axes_pad))
            if h:
                h.append(self._horiz_pad_size)

            h_ax_pos.append(len(h))

            sz = Size.Scaled(1)
            h.append(sz)

        v = []

        v_ax_pos = []
        for _ in self._row_refax[::-1]:
            #if v: v.append(Size.Fixed(self._axes_pad))
            if v:
                v.append(self._vert_pad_size)

            v_ax_pos.append(len(v))
            sz = Size.Scaled(1)
            v.append(sz)

        for i in range(self.ngrids):
            col, row = self._get_col_row(i)
            locator = self._divider.new_locator(nx=h_ax_pos[col],
                                ny=v_ax_pos[self._nrows - 1 - row])
            self.axes_all[i].set_axes_locator(locator)

        self._divider.set_horizontal(h)
        self._divider.set_vertical(v)

    def _get_col_row(self, n):
        if self._direction == "column":
            col, row = divmod(n, self._nrows)
        else:
            row, col = divmod(n, self._ncols)

        return col, row

    # Good to propagate __len__ if we have __getitem__
    def __len__(self):
        return len(self.axes_all)

    def __getitem__(self, i):
        return self.axes_all[i]

    def get_geometry(self):
        """
        get geometry of the grid. Returns a tuple of two integer,
        representing number of rows and number of columns.
        """
        return self._nrows, self._ncols

    def set_axes_pad(self, axes_pad):
        "set axes_pad"
        self._axes_pad = axes_pad

        # These two lines actually differ from ones in _init_axes_pad
        self._horiz_pad_size.fixed_size = axes_pad[0]
        self._vert_pad_size.fixed_size = axes_pad[1]

    def get_axes_pad(self):
        """
        get axes_pad

        Returns
        -------
        tuple
            Padding in inches, (horizontal pad, vertical pad)
        """
        return self._axes_pad

    def set_aspect(self, aspect):
        "set aspect"
        self._divider.set_aspect(aspect)

    def get_aspect(self):
        "get aspect"
        return self._divider.get_aspect()

    def set_label_mode(self, mode):
        "set label_mode"
        if mode == "all":
            for ax in self.axes_all:
                _tick_only(ax, False, False)
        elif mode == "L":
            # left-most axes
            for ax in self.axes_column[0][:-1]:
                _tick_only(ax, bottom_on=True, left_on=False)
            # lower-left axes
            ax = self.axes_column[0][-1]
            _tick_only(ax, bottom_on=False, left_on=False)

            for col in self.axes_column[1:]:
                # axes with no labels
                for ax in col[:-1]:
                    _tick_only(ax, bottom_on=True, left_on=True)

                # bottom
                ax = col[-1]
                _tick_only(ax, bottom_on=False, left_on=True)

        elif mode == "1":
            for ax in self.axes_all:
                _tick_only(ax, bottom_on=True, left_on=True)

            ax = self.axes_llc
            _tick_only(ax, bottom_on=False, left_on=False)

    def get_divider(self):
        return self._divider

    def set_axes_locator(self, locator):
        self._divider.set_locator(locator)

    def get_axes_locator(self):
        return self._divider.get_locator()

    def get_vsize_hsize(self):

        return self._divider.get_vsize_hsize()
#         from axes_size import AddList

#         vsize = AddList(self._divider.get_vertical())
#         hsize = AddList(self._divider.get_horizontal())

#         return vsize, hsize


class ImageGrid(Grid):
    """
    A class that creates a grid of Axes. In matplotlib, the axes
    location (and size) is specified in the normalized figure
    coordinates. This may not be ideal for images that needs to be
    displayed with a given aspect ratio.  For example, displaying
    images of a same size with some fixed padding between them cannot
    be easily done in matplotlib. ImageGrid is used in such case.
    """

    _defaultCbarAxesClass = CbarAxes

    def __init__(self, fig,
                 rect,
                 nrows_ncols,
                 ngrids=None,
                 direction="row",
                 axes_pad=0.02,
                 add_all=True,
                 share_all=False,
                 aspect=True,
                 label_mode="L",
                 cbar_mode=None,
                 cbar_location="right",
                 cbar_pad=None,
                 cbar_size="5%",
                 cbar_set_cax=True,
                 axes_class=None,
                 ):
        """
        Build an :class:`ImageGrid` instance with a grid nrows*ncols
        :class:`~matplotlib.axes.Axes` in
        :class:`~matplotlib.figure.Figure` *fig* with
        *rect=[left, bottom, width, height]* (in
        :class:`~matplotlib.figure.Figure` coordinates) or
        the subplot position code (e.g., "121").

        Optional keyword arguments:

          ================  ========  =========================================
          Keyword           Default   Description
          ================  ========  =========================================
          direction         "row"     [ "row" | "column" ]
          axes_pad          0.02      float| pad between axes given in inches
                                      or tuple-like of floats,
                                      (horizontal padding, vertical padding)
          add_all           True      bool
          share_all         False     bool
          aspect            True      bool
          label_mode        "L"       [ "L" | "1" | "all" ]
          cbar_mode         None      [ "each" | "single" | "edge" ]
          cbar_location     "right"   [ "left" | "right" | "bottom" | "top" ]
          cbar_pad          None
          cbar_size         "5%"
          cbar_set_cax      True      bool
          axes_class        None      a type object which must be a subclass
                                      of axes_grid's subclass of
                                      :class:`~matplotlib.axes.Axes`
          ================  ========  =========================================

        *cbar_set_cax* : if True, each axes in the grid has a cax
          attribute that is bind to associated cbar_axes.
        """
        self._nrows, self._ncols = nrows_ncols

        if ngrids is None:
            ngrids = self._nrows * self._ncols
        else:
            if not 0 <= ngrids < self._nrows * self._ncols:
                raise Exception

        self.ngrids = ngrids

        axes_pad = _extend_axes_pad(axes_pad)
        self._axes_pad = axes_pad

        self._colorbar_mode = cbar_mode
        self._colorbar_location = cbar_location
        if cbar_pad is None:
            # horizontal or vertical arrangement?
            if cbar_location in ("left", "right"):
                self._colorbar_pad = axes_pad[0]
            else:
                self._colorbar_pad = axes_pad[1]
        else:
            self._colorbar_pad = cbar_pad

        self._colorbar_size = cbar_size

        self._init_axes_pad(axes_pad)

        if direction not in ["column", "row"]:
            raise Exception("")

        self._direction = direction

        if axes_class is None:
            axes_class = self._defaultLocatableAxesClass
            axes_class_args = {}
        else:
            if isinstance(axes_class, maxes.Axes):
                axes_class_args = {}
            else:
                axes_class, axes_class_args = axes_class

        self.axes_all = []
        self.axes_column = [[] for _ in range(self._ncols)]
        self.axes_row = [[] for _ in range(self._nrows)]

        self.cbar_axes = []

        h = []
        v = []
        if isinstance(rect, six.string_types) or cbook.is_numlike(rect):
            self._divider = SubplotDivider(fig, rect, horizontal=h, vertical=v,
                                           aspect=aspect)
        elif isinstance(rect, SubplotSpec):
            self._divider = SubplotDivider(fig, rect, horizontal=h, vertical=v,
                                           aspect=aspect)
        elif len(rect) == 3:
            kw = dict(horizontal=h, vertical=v, aspect=aspect)
            self._divider = SubplotDivider(fig, *rect, **kw)
        elif len(rect) == 4:
            self._divider = Divider(fig, rect, horizontal=h, vertical=v,
                                    aspect=aspect)
        else:
            raise Exception("")

        rect = self._divider.get_position()

        # reference axes
        self._column_refax = [None for _ in range(self._ncols)]
        self._row_refax = [None for _ in range(self._nrows)]
        self._refax = None

        for i in range(self.ngrids):

            col, row = self._get_col_row(i)

            if share_all:
                if self.axes_all:
                    sharex = self.axes_all[0]
                    sharey = self.axes_all[0]
                else:
                    sharex = None
                    sharey = None
            else:
                sharex = self._column_refax[col]
                sharey = self._row_refax[row]

            ax = axes_class(fig, rect, sharex=sharex, sharey=sharey,
                            **axes_class_args)

            self.axes_all.append(ax)
            self.axes_column[col].append(ax)
            self.axes_row[row].append(ax)

            if share_all:
                if self._refax is None:
                    self._refax = ax
            if sharex is None:
                self._column_refax[col] = ax
            if sharey is None:
                self._row_refax[row] = ax

            cax = self._defaultCbarAxesClass(fig, rect,
                                        orientation=self._colorbar_location)
            self.cbar_axes.append(cax)

        self.axes_llc = self.axes_column[0][-1]

        self._update_locators()

        if add_all:
            for ax in self.axes_all+self.cbar_axes:
                fig.add_axes(ax)

        if cbar_set_cax:
            if self._colorbar_mode == "single":
                for ax in self.axes_all:
                    ax.cax = self.cbar_axes[0]
            elif self._colorbar_mode == "edge":
                for index, ax in enumerate(self.axes_all):
                    col, row = self._get_col_row(index)
                    if self._colorbar_location in ("left", "right"):
                        ax.cax = self.cbar_axes[row]
                    else:
                        ax.cax = self.cbar_axes[col]
            else:
                for ax, cax in zip(self.axes_all, self.cbar_axes):
                    ax.cax = cax

        self.set_label_mode(label_mode)

    def _update_locators(self):

        h = []
        v = []

        h_ax_pos = []
        h_cb_pos = []
        if (self._colorbar_mode == "single" and
             self._colorbar_location in ('left', 'bottom')):
            if self._colorbar_location == "left":
                #sz = Size.Fraction(Size.AxesX(self.axes_llc), self._nrows)
                sz = Size.Fraction(self._nrows, Size.AxesX(self.axes_llc))
                h.append(Size.from_any(self._colorbar_size, sz))
                h.append(Size.from_any(self._colorbar_pad, sz))
                locator = self._divider.new_locator(nx=0, ny=0, ny1=-1)
            elif self._colorbar_location == "bottom":
                #sz = Size.Fraction(Size.AxesY(self.axes_llc), self._ncols)
                sz = Size.Fraction(self._ncols, Size.AxesY(self.axes_llc))
                v.append(Size.from_any(self._colorbar_size, sz))
                v.append(Size.from_any(self._colorbar_pad, sz))
                locator = self._divider.new_locator(nx=0, nx1=-1, ny=0)
            for i in range(self.ngrids):
                self.cbar_axes[i].set_visible(False)
            self.cbar_axes[0].set_axes_locator(locator)
            self.cbar_axes[0].set_visible(True)

        for col, ax in enumerate(self.axes_row[0]):
            if h:
                h.append(self._horiz_pad_size)  # Size.Fixed(self._axes_pad))

            if ax:
                sz = Size.AxesX(ax, aspect="axes", ref_ax=self.axes_all[0])
            else:
                sz = Size.AxesX(self.axes_all[0],
                                aspect="axes", ref_ax=self.axes_all[0])

            if (self._colorbar_mode == "each" or
                    (self._colorbar_mode == 'edge' and
                        col == 0)) and self._colorbar_location == "left":
                h_cb_pos.append(len(h))
                h.append(Size.from_any(self._colorbar_size, sz))
                h.append(Size.from_any(self._colorbar_pad, sz))

            h_ax_pos.append(len(h))

            h.append(sz)

            if ((self._colorbar_mode == "each" or
                    (self._colorbar_mode == 'edge' and
                        col == self._ncols - 1)) and
                    self._colorbar_location == "right"):
                h.append(Size.from_any(self._colorbar_pad, sz))
                h_cb_pos.append(len(h))
                h.append(Size.from_any(self._colorbar_size, sz))

        v_ax_pos = []
        v_cb_pos = []
        for row, ax in enumerate(self.axes_column[0][::-1]):
            if v:
                v.append(self._vert_pad_size)  # Size.Fixed(self._axes_pad))

            if ax:
                sz = Size.AxesY(ax, aspect="axes", ref_ax=self.axes_all[0])
            else:
                sz = Size.AxesY(self.axes_all[0],
                                aspect="axes", ref_ax=self.axes_all[0])

            if (self._colorbar_mode == "each" or
                    (self._colorbar_mode == 'edge' and
                        row == 0)) and self._colorbar_location == "bottom":
                v_cb_pos.append(len(v))
                v.append(Size.from_any(self._colorbar_size, sz))
                v.append(Size.from_any(self._colorbar_pad, sz))

            v_ax_pos.append(len(v))
            v.append(sz)

            if ((self._colorbar_mode == "each" or
                    (self._colorbar_mode == 'edge' and
                        row == self._nrows - 1)) and
                        self._colorbar_location == "top"):
                v.append(Size.from_any(self._colorbar_pad, sz))
                v_cb_pos.append(len(v))
                v.append(Size.from_any(self._colorbar_size, sz))

        for i in range(self.ngrids):
            col, row = self._get_col_row(i)
            #locator = self._divider.new_locator(nx=4*col,
            #                                    ny=2*(self._nrows - row - 1))
            locator = self._divider.new_locator(nx=h_ax_pos[col],
                                                ny=v_ax_pos[self._nrows-1-row])
            self.axes_all[i].set_axes_locator(locator)

            if self._colorbar_mode == "each":
                if self._colorbar_location in ("right", "left"):
                    locator = self._divider.new_locator(
                        nx=h_cb_pos[col], ny=v_ax_pos[self._nrows - 1 - row])

                elif self._colorbar_location in ("top", "bottom"):
                    locator = self._divider.new_locator(
                        nx=h_ax_pos[col], ny=v_cb_pos[self._nrows - 1 - row])

                self.cbar_axes[i].set_axes_locator(locator)
            elif self._colorbar_mode == 'edge':
                if ((self._colorbar_location == 'left' and col == 0) or
                        (self._colorbar_location == 'right'
                         and col == self._ncols-1)):
                    locator = self._divider.new_locator(
                        nx=h_cb_pos[0], ny=v_ax_pos[self._nrows -1 - row])
                    self.cbar_axes[row].set_axes_locator(locator)
                elif ((self._colorbar_location == 'bottom' and
                       row == self._nrows - 1) or
                        (self._colorbar_location == 'top' and row == 0)):
                    locator = self._divider.new_locator(nx=h_ax_pos[col],
                                                        ny=v_cb_pos[0])
                    self.cbar_axes[col].set_axes_locator(locator)

        if self._colorbar_mode == "single":
            if self._colorbar_location == "right":
                #sz = Size.Fraction(Size.AxesX(self.axes_llc), self._nrows)
                sz = Size.Fraction(self._nrows, Size.AxesX(self.axes_llc))
                h.append(Size.from_any(self._colorbar_pad, sz))
                h.append(Size.from_any(self._colorbar_size, sz))
                locator = self._divider.new_locator(nx=-2, ny=0, ny1=-1)
            elif self._colorbar_location == "top":
                #sz = Size.Fraction(Size.AxesY(self.axes_llc), self._ncols)
                sz = Size.Fraction(self._ncols, Size.AxesY(self.axes_llc))
                v.append(Size.from_any(self._colorbar_pad, sz))
                v.append(Size.from_any(self._colorbar_size, sz))
                locator = self._divider.new_locator(nx=0, nx1=-1, ny=-2)
            if self._colorbar_location in ("right", "top"):
                for i in range(self.ngrids):
                    self.cbar_axes[i].set_visible(False)
                self.cbar_axes[0].set_axes_locator(locator)
                self.cbar_axes[0].set_visible(True)
        elif self._colorbar_mode == "each":
            for i in range(self.ngrids):
                self.cbar_axes[i].set_visible(True)
        elif self._colorbar_mode == "edge":
            if self._colorbar_location in ('right', 'left'):
                count = self._nrows
            else:
                count = self._ncols
            for i in range(count):
                self.cbar_axes[i].set_visible(True)
            for j in range(i + 1, self.ngrids):
                self.cbar_axes[j].set_visible(False)
        else:
            for i in range(self.ngrids):
                self.cbar_axes[i].set_visible(False)
                self.cbar_axes[i].set_position([1., 1., 0.001, 0.001],
                                               which="active")

        self._divider.set_horizontal(h)
        self._divider.set_vertical(v)


AxesGrid = ImageGrid

