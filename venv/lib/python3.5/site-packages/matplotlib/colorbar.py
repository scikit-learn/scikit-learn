'''
Colorbar toolkit with two classes and a function:

    :class:`ColorbarBase`
        the base class with full colorbar drawing functionality.
        It can be used as-is to make a colorbar for a given colormap;
        a mappable object (e.g., image) is not needed.

    :class:`Colorbar`
        the derived class for use with images or contour plots.

    :func:`make_axes`
        a function for resizing an axes and adding a second axes
        suitable for a colorbar

The :meth:`~matplotlib.figure.Figure.colorbar` method uses :func:`make_axes`
and :class:`Colorbar`; the :func:`~matplotlib.pyplot.colorbar` function
is a thin wrapper over :meth:`~matplotlib.figure.Figure.colorbar`.

'''
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves import xrange, zip

import warnings

import numpy as np

import matplotlib as mpl
import matplotlib.artist as martist
import matplotlib.cbook as cbook
import matplotlib.collections as collections
import matplotlib.colors as colors
import matplotlib.contour as contour
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.path as mpath
import matplotlib.ticker as ticker
import matplotlib.transforms as mtransforms
import matplotlib._layoutbox as layoutbox
import matplotlib._constrained_layout as constrained_layout
from matplotlib import docstring

make_axes_kw_doc = '''

    ============= ====================================================
    Property      Description
    ============= ====================================================
    *orientation* vertical or horizontal
    *fraction*    0.15; fraction of original axes to use for colorbar
    *pad*         0.05 if vertical, 0.15 if horizontal; fraction
                  of original axes between colorbar and new image axes
    *shrink*      1.0; fraction by which to multiply the size of the colorbar
    *aspect*      20; ratio of long to short dimensions
    *anchor*      (0.0, 0.5) if vertical; (0.5, 1.0) if horizontal;
                  the anchor point of the colorbar axes
    *panchor*     (1.0, 0.5) if vertical; (0.5, 0.0) if horizontal;
                  the anchor point of the colorbar parent axes. If
                  False, the parent axes' anchor will be unchanged
    ============= ====================================================

'''

colormap_kw_doc = '''

    ============  ====================================================
    Property      Description
    ============  ====================================================
    *extend*      [ 'neither' | 'both' | 'min' | 'max' ]
                  If not 'neither', make pointed end(s) for out-of-
                  range values.  These are set for a given colormap
                  using the colormap set_under and set_over methods.
    *extendfrac*  [ *None* | 'auto' | length | lengths ]
                  If set to *None*, both the minimum and maximum
                  triangular colorbar extensions with have a length of
                  5% of the interior colorbar length (this is the
                  default setting). If set to 'auto', makes the
                  triangular colorbar extensions the same lengths as
                  the interior boxes (when *spacing* is set to
                  'uniform') or the same lengths as the respective
                  adjacent interior boxes (when *spacing* is set to
                  'proportional'). If a scalar, indicates the length
                  of both the minimum and maximum triangular colorbar
                  extensions as a fraction of the interior colorbar
                  length. A two-element sequence of fractions may also
                  be given, indicating the lengths of the minimum and
                  maximum colorbar extensions respectively as a
                  fraction of the interior colorbar length.
    *extendrect*  bool
                  If *False* the minimum and maximum colorbar extensions
                  will be triangular (the default). If *True* the
                  extensions will be rectangular.
    *spacing*     [ 'uniform' | 'proportional' ]
                  Uniform spacing gives each discrete color the same
                  space; proportional makes the space proportional to
                  the data interval.
    *ticks*       [ None | list of ticks | Locator object ]
                  If None, ticks are determined automatically from the
                  input.
    *format*      [ None | format string | Formatter object ]
                  If None, the
                  :class:`~matplotlib.ticker.ScalarFormatter` is used.
                  If a format string is given, e.g., '%.3f', that is
                  used. An alternative
                  :class:`~matplotlib.ticker.Formatter` object may be
                  given instead.
    *drawedges*   bool
                  Whether to draw lines at color boundaries.
    ============  ====================================================

    The following will probably be useful only in the context of
    indexed colors (that is, when the mappable has norm=NoNorm()),
    or other unusual circumstances.

    ============   ===================================================
    Property       Description
    ============   ===================================================
    *boundaries*   None or a sequence
    *values*       None or a sequence which must be of length 1 less
                   than the sequence of *boundaries*. For each region
                   delimited by adjacent entries in *boundaries*, the
                   color mapped to the corresponding value in values
                   will be used.
    ============   ===================================================

'''

colorbar_doc = '''

Add a colorbar to a plot.

Function signatures for the :mod:`~matplotlib.pyplot` interface; all
but the first are also method signatures for the
:meth:`~matplotlib.figure.Figure.colorbar` method::

  colorbar(**kwargs)
  colorbar(mappable, **kwargs)
  colorbar(mappable, cax=cax, **kwargs)
  colorbar(mappable, ax=ax, **kwargs)

Parameters
----------
mappable :
    The :class:`~matplotlib.image.Image`,
    :class:`~matplotlib.contour.ContourSet`, etc. to
    which the colorbar applies; this argument is mandatory for the Figure
    :meth:`~matplotlib.figure.Figure.colorbar` method but optional for the
    pyplot :func:`~matplotlib.pyplot.colorbar` function, which sets the
    default to the current image.

cax : :class:`~matplotlib.axes.Axes` object, optional
    Axis into which the colorbar will be drawn

ax : :class:`~matplotlib.axes.Axes`, list of Axes, optional
    Parent axes from which space for a new colorbar axes will be stolen.
    If a list of axes is given they will all be resized to make room for the
    colorbar axes.

use_gridspec : bool, optional
    If *cax* is ``None``, a new *cax* is created as an instance of
    Axes. If *ax* is an instance of Subplot and *use_gridspec* is ``True``,
    *cax* is created as an instance of Subplot using the
    grid_spec module.


Returns
-------
:class:`~matplotlib.colorbar.Colorbar` instance
    See also its base class, :class:`~matplotlib.colorbar.ColorbarBase`.
    Call the :meth:`~matplotlib.colorbar.ColorbarBase.set_label` method
    to label the colorbar.

Notes
-----
Additional keyword arguments are of two kinds:

  axes properties:
%s
  colorbar properties:
%s

If *mappable* is a :class:`~matplotlib.contours.ContourSet`, its *extend*
kwarg is included automatically.

The *shrink* kwarg provides a simple way to scale the colorbar with respect
to the axes. Note that if *cax* is specified it determines the size of the
colorbar and *shrink* and *aspect* kwargs are ignored.

For more precise control, you can manually specify the positions of
the axes objects in which the mappable and the colorbar are drawn.  In
this case, do not use any of the axes properties kwargs.

It is known that some vector graphics viewer (svg and pdf) renders white gaps
between segments of the colorbar. This is due to bugs in the viewers not
matplotlib. As a workaround the colorbar can be rendered with overlapping
segments::

    cbar = colorbar()
    cbar.solids.set_edgecolor("face")
    draw()

However this has negative consequences in other circumstances. Particularly
with semi transparent images (alpha < 1) and colorbar extensions and is not
enabled by default see (issue #1188).

''' % (make_axes_kw_doc, colormap_kw_doc)

docstring.interpd.update(colorbar_doc=colorbar_doc)


def _set_ticks_on_axis_warn(*args, **kw):
    # a top level function which gets put in at the axes'
    # set_xticks set_yticks by _patch_ax
    warnings.warn("Use the colorbar set_ticks() method instead.")


class ColorbarBase(cm.ScalarMappable):
    '''
    Draw a colorbar in an existing axes.

    This is a base class for the :class:`Colorbar` class, which is the
    basis for the :func:`~matplotlib.pyplot.colorbar` function and the
    :meth:`~matplotlib.figure.Figure.colorbar` method, which are the
    usual ways of creating a colorbar.

    It is also useful by itself for showing a colormap.  If the *cmap*
    kwarg is given but *boundaries* and *values* are left as None,
    then the colormap will be displayed on a 0-1 scale. To show the
    under- and over-value colors, specify the *norm* as::

        colors.Normalize(clip=False)

    To show the colors versus index instead of on the 0-1 scale,
    use::

        norm=colors.NoNorm.

    Useful public methods are :meth:`set_label` and :meth:`add_lines`.

    Attributes
    ----------
    ax : Axes
        The `Axes` instance in which the colorbar is drawn.

    lines : list
        A list of `LineCollection` if lines were drawn, otherwise
        an empty list.

    dividers : LineCollection
        A LineCollection if *drawedges* is ``True``, otherwise ``None``.
    '''
    _slice_dict = {'neither': slice(0, None),
                   'both': slice(1, -1),
                   'min': slice(1, None),
                   'max': slice(0, -1)}

    n_rasterize = 50  # rasterize solids if number of colors >= n_rasterize

    def __init__(self, ax, cmap=None,
                 norm=None,
                 alpha=None,
                 values=None,
                 boundaries=None,
                 orientation='vertical',
                 ticklocation='auto',
                 extend='neither',
                 spacing='uniform',  # uniform or proportional
                 ticks=None,
                 format=None,
                 drawedges=False,
                 filled=True,
                 extendfrac=None,
                 extendrect=False,
                 label='',
                 ):
        #: The axes that this colorbar lives in.
        self.ax = ax
        self._patch_ax()
        if cmap is None:
            cmap = cm.get_cmap()
        if norm is None:
            norm = colors.Normalize()
        self.alpha = alpha
        cm.ScalarMappable.__init__(self, cmap=cmap, norm=norm)
        self.values = values
        self.boundaries = boundaries
        self.extend = extend
        self._inside = self._slice_dict[extend]
        self.spacing = spacing
        self.orientation = orientation
        self.drawedges = drawedges
        self.filled = filled
        self.extendfrac = extendfrac
        self.extendrect = extendrect
        self.solids = None
        self.lines = list()
        self.outline = None
        self.patch = None
        self.dividers = None

        if ticklocation == 'auto':
            ticklocation = 'bottom' if orientation == 'horizontal' else 'right'
        self.ticklocation = ticklocation

        self.set_label(label)
        if cbook.iterable(ticks):
            self.locator = ticker.FixedLocator(ticks, nbins=len(ticks))
        else:
            self.locator = ticks    # Handle default in _ticker()
        if format is None:
            if isinstance(self.norm, colors.LogNorm):
                self.formatter = ticker.LogFormatterSciNotation()
            elif isinstance(self.norm, colors.SymLogNorm):
                self.formatter = ticker.LogFormatterSciNotation(
                                        linthresh=self.norm.linthresh)
            else:
                self.formatter = ticker.ScalarFormatter()
        elif isinstance(format, six.string_types):
            self.formatter = ticker.FormatStrFormatter(format)
        else:
            self.formatter = format  # Assume it is a Formatter
        # The rest is in a method so we can recalculate when clim changes.
        self.config_axis()
        self.draw_all()

    def _extend_lower(self):
        """Returns whether the lower limit is open ended."""
        return self.extend in ('both', 'min')

    def _extend_upper(self):
        """Returns whether the uper limit is open ended."""
        return self.extend in ('both', 'max')

    def _patch_ax(self):
        # bind some methods to the axes to warn users
        # against using those methods.
        self.ax.set_xticks = _set_ticks_on_axis_warn
        self.ax.set_yticks = _set_ticks_on_axis_warn

    def draw_all(self):
        '''
        Calculate any free parameters based on the current cmap and norm,
        and do all the drawing.
        '''

        self._process_values()
        self._find_range()
        X, Y = self._mesh()
        C = self._values[:, np.newaxis]
        self._config_axes(X, Y)
        if self.filled:
            self._add_solids(X, Y, C)

    def config_axis(self):
        ax = self.ax
        if self.orientation == 'vertical':
            ax.xaxis.set_ticks([])
            # location is either one of 'bottom' or 'top'
            ax.yaxis.set_label_position(self.ticklocation)
            ax.yaxis.set_ticks_position(self.ticklocation)
        else:
            ax.yaxis.set_ticks([])
            # location is either one of 'left' or 'right'
            ax.xaxis.set_label_position(self.ticklocation)
            ax.xaxis.set_ticks_position(self.ticklocation)

        self._set_label()

    def update_ticks(self):
        """
        Force the update of the ticks and ticklabels. This must be
        called whenever the tick locator and/or tick formatter changes.
        """
        ax = self.ax
        ticks, ticklabels, offset_string = self._ticker()
        if self.orientation == 'vertical':
            ax.yaxis.set_ticks(ticks)
            ax.set_yticklabels(ticklabels)
            ax.yaxis.get_major_formatter().set_offset_string(offset_string)

        else:
            ax.xaxis.set_ticks(ticks)
            ax.set_xticklabels(ticklabels)
            ax.xaxis.get_major_formatter().set_offset_string(offset_string)

    def set_ticks(self, ticks, update_ticks=True):
        """
        Set tick locations.

        Parameters
        ----------
        ticks : {None, sequence, :class:`~matplotlib.ticker.Locator` instance}
            If None, a default Locator will be used.

        update_ticks : {True, False}, optional
            If True, tick locations are updated immediately.  If False,
            use :meth:`update_ticks` to manually update the ticks.

        """
        if cbook.iterable(ticks):
            self.locator = ticker.FixedLocator(ticks, nbins=len(ticks))
        else:
            self.locator = ticks

        if update_ticks:
            self.update_ticks()
        self.stale = True

    def get_ticks(self, minor=False):
        """Return the x ticks as a list of locations"""
        return self._tick_data_values

    def set_ticklabels(self, ticklabels, update_ticks=True):
        """
        set tick labels. Tick labels are updated immediately unless
        update_ticks is *False*. To manually update the ticks, call
        *update_ticks* method explicitly.
        """
        if isinstance(self.locator, ticker.FixedLocator):
            self.formatter = ticker.FixedFormatter(ticklabels)
            if update_ticks:
                self.update_ticks()
        else:
            warnings.warn("set_ticks() must have been called.")
        self.stale = True

    def _config_axes(self, X, Y):
        '''
        Make an axes patch and outline.
        '''
        ax = self.ax
        ax.set_frame_on(False)
        ax.set_navigate(False)
        xy = self._outline(X, Y)
        ax.update_datalim(xy)
        ax.set_xlim(*ax.dataLim.intervalx)
        ax.set_ylim(*ax.dataLim.intervaly)
        if self.outline is not None:
            self.outline.remove()
        self.outline = mpatches.Polygon(
            xy, edgecolor=mpl.rcParams['axes.edgecolor'],
            facecolor='none',
            linewidth=mpl.rcParams['axes.linewidth'],
            closed=True,
            zorder=2)
        ax.add_artist(self.outline)
        self.outline.set_clip_box(None)
        self.outline.set_clip_path(None)
        c = mpl.rcParams['axes.facecolor']
        if self.patch is not None:
            self.patch.remove()
        self.patch = mpatches.Polygon(xy, edgecolor=c,
                                      facecolor=c,
                                      linewidth=0.01,
                                      zorder=-1)
        ax.add_artist(self.patch)

        self.update_ticks()

    def _set_label(self):
        if self.orientation == 'vertical':
            self.ax.set_ylabel(self._label, **self._labelkw)
        else:
            self.ax.set_xlabel(self._label, **self._labelkw)
        self.stale = True

    def set_label(self, label, **kw):
        '''
        Label the long axis of the colorbar
        '''
        self._label = '%s' % (label, )
        self._labelkw = kw
        self._set_label()

    def _outline(self, X, Y):
        '''
        Return *x*, *y* arrays of colorbar bounding polygon,
        taking orientation into account.
        '''
        N = X.shape[0]
        ii = [0, 1, N - 2, N - 1, 2 * N - 1, 2 * N - 2, N + 1, N, 0]
        x = np.take(np.ravel(np.transpose(X)), ii)
        y = np.take(np.ravel(np.transpose(Y)), ii)
        x = x.reshape((len(x), 1))
        y = y.reshape((len(y), 1))
        if self.orientation == 'horizontal':
            return np.hstack((y, x))
        return np.hstack((x, y))

    def _edges(self, X, Y):
        '''
        Return the separator line segments; helper for _add_solids.
        '''
        N = X.shape[0]
        # Using the non-array form of these line segments is much
        # simpler than making them into arrays.
        if self.orientation == 'vertical':
            return [list(zip(X[i], Y[i])) for i in xrange(1, N - 1)]
        else:
            return [list(zip(Y[i], X[i])) for i in xrange(1, N - 1)]

    def _add_solids(self, X, Y, C):
        '''
        Draw the colors using :meth:`~matplotlib.axes.Axes.pcolormesh`;
        optionally add separators.
        '''
        if self.orientation == 'vertical':
            args = (X, Y, C)
        else:
            args = (np.transpose(Y), np.transpose(X), np.transpose(C))
        kw = dict(cmap=self.cmap,
                  norm=self.norm,
                  alpha=self.alpha,
                  edgecolors='None')
        # Save, set, and restore hold state to keep pcolor from
        # clearing the axes. Ordinarily this will not be needed,
        # since the axes object should already have hold set.
        _hold = self.ax._hold
        self.ax._hold = True
        col = self.ax.pcolormesh(*args, **kw)
        self.ax._hold = _hold
        #self.add_observer(col) # We should observe, not be observed...

        if self.solids is not None:
            self.solids.remove()
        self.solids = col
        if self.dividers is not None:
            self.dividers.remove()
            self.dividers = None
        if self.drawedges:
            linewidths = (0.5 * mpl.rcParams['axes.linewidth'],)
            self.dividers = collections.LineCollection(
                    self._edges(X, Y),
                    colors=(mpl.rcParams['axes.edgecolor'],),
                    linewidths=linewidths)
            self.ax.add_collection(self.dividers)
        elif len(self._y) >= self.n_rasterize:
            self.solids.set_rasterized(True)

    def add_lines(self, levels, colors, linewidths, erase=True):
        '''
        Draw lines on the colorbar.

        *colors* and *linewidths* must be scalars or
        sequences the same length as *levels*.

        Set *erase* to False to add lines without first
        removing any previously added lines.
        '''
        y = self._locate(levels)
        igood = (y < 1.001) & (y > -0.001)
        y = y[igood]
        if cbook.iterable(colors):
            colors = np.asarray(colors)[igood]
        if cbook.iterable(linewidths):
            linewidths = np.asarray(linewidths)[igood]
        N = len(y)
        x = np.array([0.0, 1.0])
        X, Y = np.meshgrid(x, y)
        if self.orientation == 'vertical':
            xy = [list(zip(X[i], Y[i])) for i in xrange(N)]
        else:
            xy = [list(zip(Y[i], X[i])) for i in xrange(N)]
        col = collections.LineCollection(xy, linewidths=linewidths)

        if erase and self.lines:
            for lc in self.lines:
                lc.remove()
            self.lines = []
        self.lines.append(col)
        col.set_color(colors)
        self.ax.add_collection(col)
        self.stale = True

    def _ticker(self):
        '''
        Return the sequence of ticks (colorbar data locations),
        ticklabels (strings), and the corresponding offset string.
        '''
        locator = self.locator
        formatter = self.formatter
        if locator is None:
            if self.boundaries is None:
                if isinstance(self.norm, colors.NoNorm):
                    nv = len(self._values)
                    base = 1 + int(nv / 10)
                    locator = ticker.IndexLocator(base=base, offset=0)
                elif isinstance(self.norm, colors.BoundaryNorm):
                    b = self.norm.boundaries
                    locator = ticker.FixedLocator(b, nbins=10)
                elif isinstance(self.norm, colors.LogNorm):
                    locator = ticker.LogLocator(subs='all')
                elif isinstance(self.norm, colors.SymLogNorm):
                    # The subs setting here should be replaced
                    # by logic in the locator.
                    locator = ticker.SymmetricalLogLocator(
                                      subs=np.arange(1, 10),
                                      linthresh=self.norm.linthresh,
                                      base=10)
                else:
                    if mpl.rcParams['_internal.classic_mode']:
                        locator = ticker.MaxNLocator()
                    else:
                        locator = ticker.AutoLocator()
            else:
                b = self._boundaries[self._inside]
                locator = ticker.FixedLocator(b, nbins=10)
        if isinstance(self.norm, colors.NoNorm) and self.boundaries is None:
            intv = self._values[0], self._values[-1]
        else:
            intv = self.vmin, self.vmax
        locator.create_dummy_axis(minpos=intv[0])
        formatter.create_dummy_axis(minpos=intv[0])
        locator.set_view_interval(*intv)
        locator.set_data_interval(*intv)
        formatter.set_view_interval(*intv)
        formatter.set_data_interval(*intv)

        b = np.array(locator())
        if isinstance(locator, ticker.LogLocator):
            eps = 1e-10
            b = b[(b <= intv[1] * (1 + eps)) & (b >= intv[0] * (1 - eps))]
        else:
            eps = (intv[1] - intv[0]) * 1e-10
            b = b[(b <= intv[1] + eps) & (b >= intv[0] - eps)]
        self._tick_data_values = b
        ticks = self._locate(b)
        formatter.set_locs(b)
        ticklabels = [formatter(t, i) for i, t in enumerate(b)]
        offset_string = formatter.get_offset()
        return ticks, ticklabels, offset_string

    def _process_values(self, b=None):
        '''
        Set the :attr:`_boundaries` and :attr:`_values` attributes
        based on the input boundaries and values.  Input boundaries
        can be *self.boundaries* or the argument *b*.
        '''
        if b is None:
            b = self.boundaries
        if b is not None:
            self._boundaries = np.asarray(b, dtype=float)
            if self.values is None:
                self._values = 0.5 * (self._boundaries[:-1]
                                      + self._boundaries[1:])
                if isinstance(self.norm, colors.NoNorm):
                    self._values = (self._values + 0.00001).astype(np.int16)
                return
            self._values = np.array(self.values)
            return
        if self.values is not None:
            self._values = np.array(self.values)
            if self.boundaries is None:
                b = np.zeros(len(self.values) + 1, 'd')
                b[1:-1] = 0.5 * (self._values[:-1] - self._values[1:])
                b[0] = 2.0 * b[1] - b[2]
                b[-1] = 2.0 * b[-2] - b[-3]
                self._boundaries = b
                return
            self._boundaries = np.array(self.boundaries)
            return
        # Neither boundaries nor values are specified;
        # make reasonable ones based on cmap and norm.
        if isinstance(self.norm, colors.NoNorm):
            b = self._uniform_y(self.cmap.N + 1) * self.cmap.N - 0.5
            v = np.zeros((len(b) - 1,), dtype=np.int16)
            v[self._inside] = np.arange(self.cmap.N, dtype=np.int16)
            if self._extend_lower():
                v[0] = -1
            if self._extend_upper():
                v[-1] = self.cmap.N
            self._boundaries = b
            self._values = v
            return
        elif isinstance(self.norm, colors.BoundaryNorm):
            b = list(self.norm.boundaries)
            if self._extend_lower():
                b = [b[0] - 1] + b
            if self._extend_upper():
                b = b + [b[-1] + 1]
            b = np.array(b)
            v = np.zeros((len(b) - 1,), dtype=float)
            bi = self.norm.boundaries
            v[self._inside] = 0.5 * (bi[:-1] + bi[1:])
            if self._extend_lower():
                v[0] = b[0] - 1
            if self._extend_upper():
                v[-1] = b[-1] + 1
            self._boundaries = b
            self._values = v
            return
        else:
            if not self.norm.scaled():
                self.norm.vmin = 0
                self.norm.vmax = 1

            self.norm.vmin, self.norm.vmax = mtransforms.nonsingular(
                self.norm.vmin,
                self.norm.vmax,
                expander=0.1)

            b = self.norm.inverse(self._uniform_y(self.cmap.N + 1))

            if isinstance(self.norm, colors.LogNorm):
                # If using a lognorm, ensure extensions don't go negative
                if self._extend_lower():
                    b[0] = 0.9 * b[0]
                if self._extend_upper():
                    b[-1] = 1.1 * b[-1]
            else:
                if self._extend_lower():
                    b[0] = b[0] - 1
                if self._extend_upper():
                    b[-1] = b[-1] + 1
        self._process_values(b)

    def _find_range(self):
        '''
        Set :attr:`vmin` and :attr:`vmax` attributes to the first and
        last boundary excluding extended end boundaries.
        '''
        b = self._boundaries[self._inside]
        self.vmin = b[0]
        self.vmax = b[-1]

    def _central_N(self):
        '''number of boundaries **before** extension of ends'''
        nb = len(self._boundaries)
        if self.extend == 'both':
            nb -= 2
        elif self.extend in ('min', 'max'):
            nb -= 1
        return nb

    def _extended_N(self):
        '''
        Based on the colormap and extend variable, return the
        number of boundaries.
        '''
        N = self.cmap.N + 1
        if self.extend == 'both':
            N += 2
        elif self.extend in ('min', 'max'):
            N += 1
        return N

    def _get_extension_lengths(self, frac, automin, automax, default=0.05):
        '''
        Get the lengths of colorbar extensions.

        A helper method for _uniform_y and _proportional_y.
        '''
        # Set the default value.
        extendlength = np.array([default, default])
        if isinstance(frac, six.string_types):
            if frac.lower() == 'auto':
                # Use the provided values when 'auto' is required.
                extendlength[0] = automin
                extendlength[1] = automax
            else:
                # Any other string is invalid.
                raise ValueError('invalid value for extendfrac')
        elif frac is not None:
            try:
                # Try to set min and max extension fractions directly.
                extendlength[:] = frac
                # If frac is a sequence containing None then NaN may
                # be encountered. This is an error.
                if np.isnan(extendlength).any():
                    raise ValueError()
            except (TypeError, ValueError):
                # Raise an error on encountering an invalid value for frac.
                raise ValueError('invalid value for extendfrac')
        return extendlength

    def _uniform_y(self, N):
        '''
        Return colorbar data coordinates for *N* uniformly
        spaced boundaries, plus ends if required.
        '''
        if self.extend == 'neither':
            y = np.linspace(0, 1, N)
        else:
            automin = automax = 1. / (N - 1.)
            extendlength = self._get_extension_lengths(self.extendfrac,
                                                       automin, automax,
                                                       default=0.05)
            if self.extend == 'both':
                y = np.zeros(N + 2, 'd')
                y[0] = 0. - extendlength[0]
                y[-1] = 1. + extendlength[1]
            elif self.extend == 'min':
                y = np.zeros(N + 1, 'd')
                y[0] = 0. - extendlength[0]
            else:
                y = np.zeros(N + 1, 'd')
                y[-1] = 1. + extendlength[1]
            y[self._inside] = np.linspace(0, 1, N)
        return y

    def _proportional_y(self):
        '''
        Return colorbar data coordinates for the boundaries of
        a proportional colorbar.
        '''
        if isinstance(self.norm, colors.BoundaryNorm):
            y = (self._boundaries - self._boundaries[0])
            y = y / (self._boundaries[-1] - self._boundaries[0])
        else:
            y = self.norm(self._boundaries.copy())
            y = np.ma.filled(y, np.nan)
        if self.extend == 'min':
            # Exclude leftmost interval of y.
            clen = y[-1] - y[1]
            automin = (y[2] - y[1]) / clen
            automax = (y[-1] - y[-2]) / clen
        elif self.extend == 'max':
            # Exclude rightmost interval in y.
            clen = y[-2] - y[0]
            automin = (y[1] - y[0]) / clen
            automax = (y[-2] - y[-3]) / clen
        elif self.extend == 'both':
            # Exclude leftmost and rightmost intervals in y.
            clen = y[-2] - y[1]
            automin = (y[2] - y[1]) / clen
            automax = (y[-2] - y[-3]) / clen
        if self.extend in ('both', 'min', 'max'):
            extendlength = self._get_extension_lengths(self.extendfrac,
                                                       automin, automax,
                                                       default=0.05)
        if self.extend in ('both', 'min'):
            y[0] = 0. - extendlength[0]
        if self.extend in ('both', 'max'):
            y[-1] = 1. + extendlength[1]
        yi = y[self._inside]
        norm = colors.Normalize(yi[0], yi[-1])
        y[self._inside] = np.ma.filled(norm(yi), np.nan)
        return y

    def _mesh(self):
        '''
        Return X,Y, the coordinate arrays for the colorbar pcolormesh.
        These are suitable for a vertical colorbar; swapping and
        transposition for a horizontal colorbar are done outside
        this function.
        '''
        x = np.array([0.0, 1.0])
        if self.spacing == 'uniform':
            y = self._uniform_y(self._central_N())
        else:
            y = self._proportional_y()
        self._y = y
        X, Y = np.meshgrid(x, y)
        if self._extend_lower() and not self.extendrect:
            X[0, :] = 0.5
        if self._extend_upper() and not self.extendrect:
            X[-1, :] = 0.5
        return X, Y

    def _locate(self, x):
        '''
        Given a set of color data values, return their
        corresponding colorbar data coordinates.
        '''
        if isinstance(self.norm, (colors.NoNorm, colors.BoundaryNorm)):
            b = self._boundaries
            xn = x
        else:
            # Do calculations using normalized coordinates so
            # as to make the interpolation more accurate.
            b = self.norm(self._boundaries, clip=False).filled()
            xn = self.norm(x, clip=False).filled()

        # The rest is linear interpolation with extrapolation at ends.
        ii = np.searchsorted(b, xn)
        i0 = ii - 1
        itop = (ii == len(b))
        ibot = (ii == 0)
        i0[itop] -= 1
        ii[itop] -= 1
        i0[ibot] += 1
        ii[ibot] += 1

        db = np.take(b, ii) - np.take(b, i0)
        y = self._y
        dy = np.take(y, ii) - np.take(y, i0)
        z = np.take(y, i0) + (xn - np.take(b, i0)) * dy / db
        return z

    def set_alpha(self, alpha):
        self.alpha = alpha

    def remove(self):
        """
        Remove this colorbar from the figure
        """

        fig = self.ax.figure
        fig.delaxes(self.ax)


class Colorbar(ColorbarBase):
    """
    This class connects a :class:`ColorbarBase` to a
    :class:`~matplotlib.cm.ScalarMappable` such as a
    :class:`~matplotlib.image.AxesImage` generated via
    :meth:`~matplotlib.axes.Axes.imshow`.

    It is not intended to be instantiated directly; instead,
    use :meth:`~matplotlib.figure.Figure.colorbar` or
    :func:`~matplotlib.pyplot.colorbar` to make your colorbar.

    """
    def __init__(self, ax, mappable, **kw):
        # Ensure the given mappable's norm has appropriate vmin and vmax set
        # even if mappable.draw has not yet been called.
        mappable.autoscale_None()

        self.mappable = mappable
        kw['cmap'] = cmap = mappable.cmap
        kw['norm'] = norm = mappable.norm

        if isinstance(mappable, contour.ContourSet):
            CS = mappable
            kw['alpha'] = mappable.get_alpha()
            kw['boundaries'] = CS._levels
            kw['values'] = CS.cvalues
            kw['extend'] = CS.extend
            #kw['ticks'] = CS._levels
            kw.setdefault('ticks', ticker.FixedLocator(CS.levels, nbins=10))
            kw['filled'] = CS.filled
            ColorbarBase.__init__(self, ax, **kw)
            if not CS.filled:
                self.add_lines(CS)
        else:
            if getattr(cmap, 'colorbar_extend', False) is not False:
                kw.setdefault('extend', cmap.colorbar_extend)

            if isinstance(mappable, martist.Artist):
                kw['alpha'] = mappable.get_alpha()

            ColorbarBase.__init__(self, ax, **kw)

    def on_mappable_changed(self, mappable):
        """
        Updates this colorbar to match the mappable's properties.

        Typically this is automatically registered as an event handler
        by :func:`colorbar_factory` and should not be called manually.

        """
        self.set_cmap(mappable.get_cmap())
        self.set_clim(mappable.get_clim())
        self.update_normal(mappable)

    def add_lines(self, CS, erase=True):
        '''
        Add the lines from a non-filled
        :class:`~matplotlib.contour.ContourSet` to the colorbar.

        Set *erase* to False if these lines should be added to
        any pre-existing lines.
        '''
        if not isinstance(CS, contour.ContourSet) or CS.filled:
            raise ValueError('add_lines is only for a ContourSet of lines')
        tcolors = [c[0] for c in CS.tcolors]
        tlinewidths = [t[0] for t in CS.tlinewidths]
        # The following was an attempt to get the colorbar lines
        # to follow subsequent changes in the contour lines,
        # but more work is needed: specifically, a careful
        # look at event sequences, and at how
        # to make one object track another automatically.
        #tcolors = [col.get_colors()[0] for col in CS.collections]
        #tlinewidths = [col.get_linewidth()[0] for lw in CS.collections]
        ColorbarBase.add_lines(self, CS.levels, tcolors, tlinewidths,
                               erase=erase)

    def update_normal(self, mappable):
        '''
        update solid, lines, etc. Unlike update_bruteforce, it does
        not clear the axes.  This is meant to be called when the image
        or contour plot to which this colorbar belongs is changed.
        '''
        self.draw_all()
        if isinstance(self.mappable, contour.ContourSet):
            CS = self.mappable
            if not CS.filled:
                self.add_lines(CS)
        self.stale = True

    def update_bruteforce(self, mappable):
        '''
        Destroy and rebuild the colorbar.  This is
        intended to become obsolete, and will probably be
        deprecated and then removed.  It is not called when
        the pyplot.colorbar function or the Figure.colorbar
        method are used to create the colorbar.

        '''
        # We are using an ugly brute-force method: clearing and
        # redrawing the whole thing.  The problem is that if any
        # properties have been changed by methods other than the
        # colorbar methods, those changes will be lost.
        self.ax.cla()
        # clearing the axes will delete outline, patch, solids, and lines:
        self.outline = None
        self.patch = None
        self.solids = None
        self.lines = list()
        self.dividers = None
        self.set_alpha(mappable.get_alpha())
        self.cmap = mappable.cmap
        self.norm = mappable.norm
        self.config_axis()
        self.draw_all()
        if isinstance(self.mappable, contour.ContourSet):
            CS = self.mappable
            if not CS.filled:
                self.add_lines(CS)
            #if self.lines is not None:
            #    tcolors = [c[0] for c in CS.tcolors]
            #    self.lines.set_color(tcolors)
        #Fixme? Recalculate boundaries, ticks if vmin, vmax have changed.
        #Fixme: Some refactoring may be needed; we should not
        # be recalculating everything if there was a simple alpha
        # change.

    def remove(self):
        """
        Remove this colorbar from the figure.  If the colorbar was created with
        ``use_gridspec=True`` then restore the gridspec to its previous value.
        """

        ColorbarBase.remove(self)
        self.mappable.callbacksSM.disconnect(self.mappable.colorbar_cid)
        self.mappable.colorbar = None
        self.mappable.colorbar_cid = None

        try:
            ax = self.mappable.axes
        except AttributeError:
            return

        try:
            gs = ax.get_subplotspec().get_gridspec()
            subplotspec = gs.get_topmost_subplotspec()
        except AttributeError:
            # use_gridspec was False
            pos = ax.get_position(original=True)
            ax._set_position(pos)
        else:
            # use_gridspec was True
            ax.set_subplotspec(subplotspec)


@docstring.Substitution(make_axes_kw_doc)
def make_axes(parents, location=None, orientation=None, fraction=0.15,
              shrink=1.0, aspect=20, **kw):
    '''
    Resize and reposition parent axes, and return a child
    axes suitable for a colorbar.

    Keyword arguments may include the following (with defaults):

        location : [None|'left'|'right'|'top'|'bottom']
            The position, relative to **parents**, where the colorbar axes
            should be created. If None, the value will either come from the
            given ``orientation``, else it will default to 'right'.

        orientation :  [None|'vertical'|'horizontal']
            The orientation of the colorbar. Typically, this keyword shouldn't
            be used, as it can be derived from the ``location`` keyword.

    %s

    Returns (cax, kw), the child axes and the reduced kw dictionary to be
    passed when creating the colorbar instance.
    '''
    locations = ["left", "right", "top", "bottom"]
    if orientation is not None and location is not None:
        raise TypeError('position and orientation are mutually exclusive. '
                        'Consider setting the position to any of {}'
                        .format(', '.join(locations)))

    # provide a default location
    if location is None and orientation is None:
        location = 'right'

    # allow the user to not specify the location by specifying the
    # orientation instead
    if location is None:
        location = 'right' if orientation == 'vertical' else 'bottom'

    if location not in locations:
        raise ValueError('Invalid colorbar location. Must be one '
                         'of %s' % ', '.join(locations))

    default_location_settings = {'left':   {'anchor': (1.0, 0.5),
                                            'panchor': (0.0, 0.5),
                                            'pad': 0.10,
                                            'orientation': 'vertical'},
                                 'right':  {'anchor': (0.0, 0.5),
                                            'panchor': (1.0, 0.5),
                                            'pad': 0.05,
                                            'orientation': 'vertical'},
                                 'top':    {'anchor': (0.5, 0.0),
                                            'panchor': (0.5, 1.0),
                                            'pad': 0.05,
                                            'orientation': 'horizontal'},
                                 'bottom': {'anchor': (0.5, 1.0),
                                            'panchor': (0.5, 0.0),
                                            'pad': 0.15,  # backwards compat
                                            'orientation': 'horizontal'},
                                 }

    loc_settings = default_location_settings[location]

    # put appropriate values into the kw dict for passing back to
    # the Colorbar class
    kw['orientation'] = loc_settings['orientation']
    kw['ticklocation'] = location

    anchor = kw.pop('anchor', loc_settings['anchor'])
    parent_anchor = kw.pop('panchor', loc_settings['panchor'])

    parents_iterable = cbook.iterable(parents)
    # turn parents into a list if it is not already. We do this w/ np
    # because `plt.subplots` can return an ndarray and is natural to
    # pass to `colorbar`.
    parents = np.atleast_1d(parents).ravel()

    # check if using constrained_layout:
    try:
        gs = parents[0].get_subplotspec().get_gridspec()
        using_constrained_layout = (gs._layoutbox is not None)
    except AttributeError:
        using_constrained_layout = False

    # defaults are not appropriate for constrained_layout:
    pad0 = loc_settings['pad']
    if using_constrained_layout:
        pad0 = 0.02
    pad = kw.pop('pad', pad0)

    fig = parents[0].get_figure()
    if not all(fig is ax.get_figure() for ax in parents):
        raise ValueError('Unable to create a colorbar axes as not all '
                         'parents share the same figure.')

    # take a bounding box around all of the given axes
    parents_bbox = mtransforms.Bbox.union(
        [ax.get_position(original=True).frozen() for ax in parents])

    pb = parents_bbox
    if location in ('left', 'right'):
        if location == 'left':
            pbcb, _, pb1 = pb.splitx(fraction, fraction + pad)
        else:
            pb1, _, pbcb = pb.splitx(1 - fraction - pad, 1 - fraction)
        pbcb = pbcb.shrunk(1.0, shrink).anchored(anchor, pbcb)
    else:
        if location == 'bottom':
            pbcb, _, pb1 = pb.splity(fraction, fraction + pad)
        else:
            pb1, _, pbcb = pb.splity(1 - fraction - pad, 1 - fraction)
        pbcb = pbcb.shrunk(shrink, 1.0).anchored(anchor, pbcb)

        # define the aspect ratio in terms of y's per x rather than x's per y
        aspect = 1.0 / aspect

    # define a transform which takes us from old axes coordinates to
    # new axes coordinates
    shrinking_trans = mtransforms.BboxTransform(parents_bbox, pb1)

    # transform each of the axes in parents using the new transform
    for ax in parents:
        new_posn = shrinking_trans.transform(ax.get_position())
        new_posn = mtransforms.Bbox(new_posn)
        ax._set_position(new_posn)
        if parent_anchor is not False:
            ax.set_anchor(parent_anchor)

    cax = fig.add_axes(pbcb)

    # OK, now make a layoutbox for the cb axis.  Later, we will use this
    # to make the colorbar fit nicely.
    if not using_constrained_layout:
        # no layout boxes:
        lb = None
        lbpos = None
        # and we need to set the aspect ratio by hand...
        cax.set_aspect(aspect, anchor=anchor, adjustable='box')
    else:
        if not parents_iterable:
            # this is a single axis...
            ax = parents[0]
            lb, lbpos = constrained_layout.layoutcolorbarsingle(
                    ax, cax, shrink, aspect, location, pad=pad)
        else:  # there is more than one parent, so lets use gridspec
            # the colorbar will be a sibling of this gridspec, so the
            # parent is the same parent as the gridspec.  Either the figure,
            # or a subplotspec.

            lb, lbpos = constrained_layout.layoutcolorbargridspec(
                    parents, cax, shrink, aspect, location, pad)

    cax._layoutbox = lb
    cax._poslayoutbox = lbpos

    return cax, kw


@docstring.Substitution(make_axes_kw_doc)
def make_axes_gridspec(parent, **kw):
    '''
    Resize and reposition a parent axes, and return a child axes
    suitable for a colorbar. This function is similar to
    make_axes. Prmary differences are

     * *make_axes_gridspec* only handles the *orientation* keyword
       and cannot handle the "location" keyword.

     * *make_axes_gridspec* should only be used with a subplot parent.

     * *make_axes* creates an instance of Axes. *make_axes_gridspec*
        creates an instance of Subplot.

     * *make_axes* updates the position of the
        parent. *make_axes_gridspec* replaces the grid_spec attribute
        of the parent with a new one.

    While this function is meant to be compatible with *make_axes*,
    there could be some minor differences.

    Keyword arguments may include the following (with defaults):

        *orientation*
            'vertical' or 'horizontal'

    %s

    All but the first of these are stripped from the input kw set.

    Returns (cax, kw), the child axes and the reduced kw dictionary to be
    passed when creating the colorbar instance.
    '''

    orientation = kw.setdefault('orientation', 'vertical')
    kw['ticklocation'] = 'auto'

    fraction = kw.pop('fraction', 0.15)
    shrink = kw.pop('shrink', 1.0)
    aspect = kw.pop('aspect', 20)

    x1 = 1 - fraction

    # for shrinking
    pad_s = (1 - shrink) * 0.5
    wh_ratios = [pad_s, shrink, pad_s]

    # we need to none the tree of layoutboxes because
    # constrained_layout can't remove and replace the tree
    # hierarchy w/o a seg fault.
    gs = parent.get_subplotspec().get_gridspec()
    layoutbox.nonetree(gs._layoutbox)
    gs_from_subplotspec = gridspec.GridSpecFromSubplotSpec
    if orientation == 'vertical':
        pad = kw.pop('pad', 0.05)
        wh_space = 2 * pad / (1 - pad)
        gs = gs_from_subplotspec(1, 2,
                                 subplot_spec=parent.get_subplotspec(),
                                 wspace=wh_space,
                                 width_ratios=[x1 - pad, fraction])
        gs2 = gs_from_subplotspec(3, 1,
                                  subplot_spec=gs[1],
                                  hspace=0.,
                                  height_ratios=wh_ratios)
        anchor = (0.0, 0.5)
        panchor = (1.0, 0.5)
    else:
        pad = kw.pop('pad', 0.15)
        wh_space = 2 * pad / (1 - pad)
        gs = gs_from_subplotspec(2, 1,
                                 subplot_spec=parent.get_subplotspec(),
                                 hspace=wh_space,
                                 height_ratios=[x1 - pad, fraction])
        gs2 = gs_from_subplotspec(1, 3,
                                  subplot_spec=gs[1],
                                  wspace=0.,
                                  width_ratios=wh_ratios)
        aspect = 1 / aspect
        anchor = (0.5, 1.0)
        panchor = (0.5, 0.0)

    parent.set_subplotspec(gs[0])
    parent.update_params()
    parent._set_position(parent.figbox)
    parent.set_anchor(panchor)

    fig = parent.get_figure()
    cax = fig.add_subplot(gs2[1])
    cax.set_aspect(aspect, anchor=anchor, adjustable='box')
    return cax, kw


class ColorbarPatch(Colorbar):
    """
    A Colorbar which is created using :class:`~matplotlib.patches.Patch`
    rather than the default :func:`~matplotlib.axes.pcolor`.

    It uses a list of Patch instances instead of a
    :class:`~matplotlib.collections.PatchCollection` because the
    latter does not allow the hatch pattern to vary among the
    members of the collection.
    """
    def __init__(self, ax, mappable, **kw):
        # we do not want to override the behaviour of solids
        # so add a new attribute which will be a list of the
        # colored patches in the colorbar
        self.solids_patches = []
        Colorbar.__init__(self, ax, mappable, **kw)

    def _add_solids(self, X, Y, C):
        """
        Draw the colors using :class:`~matplotlib.patches.Patch`;
        optionally add separators.
        """
        # Save, set, and restore hold state to keep pcolor from
        # clearing the axes. Ordinarily this will not be needed,
        # since the axes object should already have hold set.
        _hold = self.ax._hold
        self.ax._hold = True

        kw = {'alpha': self.alpha, }

        n_segments = len(C)

        # ensure there are sufficient hatches
        hatches = self.mappable.hatches * n_segments

        patches = []
        for i in xrange(len(X) - 1):
            val = C[i][0]
            hatch = hatches[i]

            xy = np.array([[X[i][0], Y[i][0]],
                           [X[i][1], Y[i][0]],
                           [X[i + 1][1], Y[i + 1][0]],
                           [X[i + 1][0], Y[i + 1][1]]])

            if self.orientation == 'horizontal':
                # if horizontal swap the xs and ys
                xy = xy[..., ::-1]

            patch = mpatches.PathPatch(mpath.Path(xy),
                                       facecolor=self.cmap(self.norm(val)),
                                       hatch=hatch, linewidth=0,
                                       antialiased=False, **kw)
            self.ax.add_patch(patch)
            patches.append(patch)

        if self.solids_patches:
            for solid in self.solids_patches:
                solid.remove()

        self.solids_patches = patches

        if self.dividers is not None:
            self.dividers.remove()
            self.dividers = None

        if self.drawedges:
            self.dividers = collections.LineCollection(
                    self._edges(X, Y),
                    colors=(mpl.rcParams['axes.edgecolor'],),
                    linewidths=(0.5 * mpl.rcParams['axes.linewidth'],))
            self.ax.add_collection(self.dividers)

        self.ax._hold = _hold


def colorbar_factory(cax, mappable, **kwargs):
    """
    Creates a colorbar on the given axes for the given mappable.

    Typically, for automatic colorbar placement given only a mappable use
    :meth:`~matplotlib.figure.Figure.colorbar`.

    """
    # if the given mappable is a contourset with any hatching, use
    # ColorbarPatch else use Colorbar
    if (isinstance(mappable, contour.ContourSet)
            and any([hatch is not None for hatch in mappable.hatches])):
        cb = ColorbarPatch(cax, mappable, **kwargs)
    else:
        cb = Colorbar(cax, mappable, **kwargs)

    cid = mappable.callbacksSM.connect('changed', cb.on_mappable_changed)
    mappable.colorbar = cb
    mappable.colorbar_cid = cid

    return cb
