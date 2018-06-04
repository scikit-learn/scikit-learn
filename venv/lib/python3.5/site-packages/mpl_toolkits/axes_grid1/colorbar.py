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

import numpy as np
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib import docstring
import matplotlib.ticker as ticker
import matplotlib.cbook as cbook
import matplotlib.collections as collections
import matplotlib.contour as contour
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.transforms import Bbox


make_axes_kw_doc = '''

    ============= ====================================================
    Property      Description
    ============= ====================================================
    *orientation* vertical or horizontal
    *fraction*    0.15; fraction of original axes to use for colorbar
    *pad*         0.05 if vertical, 0.15 if horizontal; fraction
                  of original axes between colorbar and new image axes
    *shrink*      1.0; fraction by which to shrink the colorbar
    *aspect*      20; ratio of long to short dimensions
    ============= ====================================================

'''

colormap_kw_doc = '''

    ===========   ====================================================
    Property      Description
    ===========   ====================================================
    *extend*      [ 'neither' | 'both' | 'min' | 'max' ]
                  If not 'neither', make pointed end(s) for out-of-
                  range values.  These are set for a given colormap
                  using the colormap set_under and set_over methods.
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
    ===========   ====================================================

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

arguments:

  *mappable*
    the :class:`~matplotlib.image.Image`,
    :class:`~matplotlib.contour.ContourSet`, etc. to
    which the colorbar applies; this argument is mandatory for the
    :meth:`~matplotlib.figure.Figure.colorbar` method but optional for the
    :func:`~matplotlib.pyplot.colorbar` function, which sets the
    default to the current image.

keyword arguments:

  *cax*
    None | axes object into which the colorbar will be drawn
  *ax*
    None | parent axes object from which space for a new
    colorbar axes will be stolen


Additional keyword arguments are of two kinds:

  axes properties:
%s
  colorbar properties:
%s

If *mappable* is a :class:`~matplotlib.contours.ContourSet`, its *extend*
kwarg is included automatically.

Note that the *shrink* kwarg provides a simple way to keep a vertical
colorbar, for example, from being taller than the axes of the mappable
to which the colorbar is attached; but it is a manual method requiring
some trial and error. If the colorbar is too tall (or a horizontal
colorbar is too wide) use a smaller value of *shrink*.

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

However this has negative consequences in other circumstances. Particularly with
semi transparent images (alpha < 1) and colorbar extensions and is not enabled
by default see (issue #1188).

returns:
    :class:`~matplotlib.colorbar.Colorbar` instance; see also its base class,
    :class:`~matplotlib.colorbar.ColorbarBase`.  Call the
    :meth:`~matplotlib.colorbar.ColorbarBase.set_label` method
    to label the colorbar.


The transData of the *cax* is adjusted so that the limits in the
longest axis actually corresponds to the limits in colorbar range. On
the other hand, the shortest axis has a data limits of [1,2], whose
unconventional value is to prevent underflow when log scale is used.
''' % (make_axes_kw_doc, colormap_kw_doc)

docstring.interpd.update(colorbar_doc=colorbar_doc)


class CbarAxesLocator(object):
    """
    CbarAxesLocator is a axes_locator for colorbar axes. It adjust the
    position of the axes to make a room for extended ends, i.e., the
    extended ends are located outside the axes area.
    """

    def __init__(self, locator=None, extend="neither", orientation="vertical"):
        """
        *locator* : the bbox returned from the locator is used as a
            initial axes location. If None, axes.bbox is used.

        *extend* : same as in ColorbarBase
        *orientation* : same as in ColorbarBase

        """
        self._locator = locator
        self.extesion_fraction = 0.05
        self.extend = extend
        self.orientation = orientation

    def get_original_position(self, axes, renderer):
        """
        get the original position of the axes.
        """
        if self._locator is None:
            bbox = axes.get_position(original=True)
        else:
            bbox = self._locator(axes, renderer)
        return bbox

    def get_end_vertices(self):
        """
        return a tuple of two vertices for the colorbar extended ends.
        The first vertices is for the minimum end, and the second is for
        the maximum end.
        """
        # Note that concatenating two vertices needs to make a
        # vertices for the frame.
        extesion_fraction = self.extesion_fraction

        corx = extesion_fraction*2.
        cory = 1./(1. - corx)
        x1, y1, w, h = 0, 0, 1, 1
        x2, y2 = x1 + w, y1 + h
        dw, dh = w*extesion_fraction, h*extesion_fraction*cory

        if self.extend in ["min", "both"]:
            bottom = [(x1, y1),
                      (x1+w/2., y1-dh),
                      (x2, y1)]
        else:
            bottom = [(x1, y1),
                      (x2, y1)]

        if self.extend in ["max", "both"]:
            top = [(x2, y2),
                   (x1+w/2., y2+dh),
                   (x1, y2)]
        else:
            top = [(x2, y2),
                   (x1, y2)]

        if self.orientation == "horizontal":
            bottom = [(y,x) for (x,y) in bottom]
            top = [(y,x) for (x,y) in top]

        return bottom, top


    def get_path_patch(self):
        """
        get the path for axes patch
        """
        end1, end2 = self.get_end_vertices()

        verts = [] + end1 + end2 + end1[:1]

        return Path(verts)


    def get_path_ends(self):
        """
        get the paths for extended ends
        """

        end1, end2 = self.get_end_vertices()

        return Path(end1), Path(end2)


    def __call__(self, axes, renderer):
        """
        Return the adjusted position of the axes
        """
        bbox0 = self.get_original_position(axes, renderer)
        bbox = bbox0

        x1, y1, w, h = bbox.bounds
        extesion_fraction = self.extesion_fraction
        dw, dh = w*extesion_fraction, h*extesion_fraction

        if self.extend in ["min", "both"]:
            if self.orientation == "horizontal":
                x1 = x1 + dw
            else:
                y1 = y1+dh

        if self.extend in ["max", "both"]:
            if self.orientation == "horizontal":
                w = w-2*dw
            else:
                h = h-2*dh

        return Bbox.from_bounds(x1, y1, w, h)



class ColorbarBase(cm.ScalarMappable):
    '''
    Draw a colorbar in an existing axes.

    This is a base class for the :class:`Colorbar` class, which is the
    basis for the :func:`~matplotlib.pyplot.colorbar` method and pylab
    function.

    It is also useful by itself for showing a colormap.  If the *cmap*
    kwarg is given but *boundaries* and *values* are left as None,
    then the colormap will be displayed on a 0-1 scale. To show the
    under- and over-value colors, specify the *norm* as::

        colors.Normalize(clip=False)

    To show the colors versus index instead of on the 0-1 scale,
    use::

        norm=colors.NoNorm.

    Useful attributes:

        :attr:`ax`
            the Axes instance in which the colorbar is drawn

        :attr:`lines`
            a LineCollection if lines were drawn, otherwise None

        :attr:`dividers`
            a LineCollection if *drawedges* is True, otherwise None

    Useful public methods are :meth:`set_label` and :meth:`add_lines`.

    '''

    def __init__(self, ax, cmap=None,
                           norm=None,
                           alpha=1.0,
                           values=None,
                           boundaries=None,
                           orientation='vertical',
                           extend='neither',
                           spacing='uniform',  # uniform or proportional
                           ticks=None,
                           format=None,
                           drawedges=False,
                           filled=True,
                           ):
        self.ax = ax

        if cmap is None: cmap = cm.get_cmap()
        if norm is None: norm = colors.Normalize()
        self.alpha = alpha
        cm.ScalarMappable.__init__(self, cmap=cmap, norm=norm)
        self.values = values
        self.boundaries = boundaries
        self.extend = extend
        self.spacing = spacing
        self.orientation = orientation
        self.drawedges = drawedges
        self.filled = filled

        # artists
        self.solids = None
        self.lines = None
        self.dividers = None
        self.extension_patch1 = None
        self.extension_patch2 = None

        if orientation == "vertical":
            self.cbar_axis = self.ax.yaxis
        else:
            self.cbar_axis = self.ax.xaxis


        if format is None:
            if isinstance(self.norm, colors.LogNorm):
                # change both axis for proper aspect
                self.ax.set_xscale("log")
                self.ax.set_yscale("log")
                self.cbar_axis.set_minor_locator(ticker.NullLocator())
                formatter = ticker.LogFormatter()
            else:
                formatter = None
        elif isinstance(format, six.string_types):
            formatter = ticker.FormatStrFormatter(format)
        else:
            formatter = format  # Assume it is a Formatter

        if formatter is None:
            formatter = self.cbar_axis.get_major_formatter()
        else:
            self.cbar_axis.set_major_formatter(formatter)

        if cbook.iterable(ticks):
            self.cbar_axis.set_ticks(ticks)
        elif ticks is not None:
            self.cbar_axis.set_major_locator(ticks)
        else:
            self._select_locator(formatter)


        self._config_axes()

        self.update_artists()

        self.set_label_text('')


    def _get_colorbar_limits(self):
        """
        initial limits for colorbar range. The returned min, max values
        will be used to create colorbar solid(?) and etc.
        """
        if self.boundaries is not None:
            C = self.boundaries
            if self.extend in ["min", "both"]:
                C = C[1:]

            if self.extend in ["max", "both"]:
                C = C[:-1]
            return min(C), max(C)
        else:
            return self.get_clim()


    def _config_axes(self):
        '''
        Adjust the properties of the axes to be adequate for colorbar display.
        '''
        ax = self.ax

        axes_locator = CbarAxesLocator(ax.get_axes_locator(),
                                       extend=self.extend,
                                       orientation=self.orientation)
        ax.set_axes_locator(axes_locator)

        # override the get_data_ratio for the aspect works.
        def _f():
            return 1.
        ax.get_data_ratio = _f
        ax.get_data_ratio_log = _f

        ax.set_frame_on(True)
        ax.set_navigate(False)

        self.ax.set_autoscalex_on(False)
        self.ax.set_autoscaley_on(False)

        if self.orientation == 'horizontal':
            ax.xaxis.set_label_position('bottom')
            ax.set_yticks([])
        else:
            ax.set_xticks([])
            ax.yaxis.set_label_position('right')
            ax.yaxis.set_ticks_position('right')



    def update_artists(self):
        """
        Update the colorbar associated artists, *filled* and
        *ends*. Note that *lines* are not updated.  This needs to be
        called whenever clim of associated image changes.
        """
        self._process_values()
        self._add_ends()

        X, Y = self._mesh()
        if self.filled:
            C = self._values[:,np.newaxis]
            self._add_solids(X, Y, C)

        ax = self.ax
        vmin, vmax = self._get_colorbar_limits()
        if self.orientation == 'horizontal':
            ax.set_ylim(1, 2)
            ax.set_xlim(vmin, vmax)
        else:
            ax.set_xlim(1, 2)
            ax.set_ylim(vmin, vmax)


    def _add_ends(self):
        """
        Create patches from extended ends and add them to the axes.
        """

        del self.extension_patch1
        del self.extension_patch2

        path1, path2 = self.ax.get_axes_locator().get_path_ends()
        fc=mpl.rcParams['axes.facecolor']
        ec=mpl.rcParams['axes.edgecolor']
        linewidths=0.5*mpl.rcParams['axes.linewidth']
        self.extension_patch1 = PathPatch(path1,
                                          fc=fc, ec=ec, lw=linewidths,
                                          zorder=2.,
                                          transform=self.ax.transAxes,
                                          clip_on=False)
        self.extension_patch2 = PathPatch(path2,
                                          fc=fc, ec=ec, lw=linewidths,
                                          zorder=2.,
                                          transform=self.ax.transAxes,
                                          clip_on=False)
        self.ax.add_artist(self.extension_patch1)
        self.ax.add_artist(self.extension_patch2)



    def _set_label_text(self):
        """
        set label.
        """
        self.cbar_axis.set_label_text(self._label, **self._labelkw)

    def set_label_text(self, label, **kw):
        '''
        Label the long axis of the colorbar
        '''
        self._label = label
        self._labelkw = kw
        self._set_label_text()


    def _edges(self, X, Y):
        '''
        Return the separator line segments; helper for _add_solids.
        '''
        N = X.shape[0]
        # Using the non-array form of these line segments is much
        # simpler than making them into arrays.
        if self.orientation == 'vertical':
            return [list(zip(X[i], Y[i])) for i in xrange(1, N-1)]
        else:
            return [list(zip(Y[i], X[i])) for i in xrange(1, N-1)]

    def _add_solids(self, X, Y, C):
        '''
        Draw the colors using :meth:`~matplotlib.axes.Axes.pcolormesh`;
        optionally add separators.
        '''
        ## Change to pcolorfast after fixing bugs in some backends...

        if self.extend in ["min", "both"]:
            cc = self.to_rgba([C[0][0]])
            self.extension_patch1.set_fc(cc[0])
            X, Y, C = X[1:], Y[1:], C[1:]

        if self.extend in ["max", "both"]:
            cc = self.to_rgba([C[-1][0]])
            self.extension_patch2.set_fc(cc[0])
            X, Y, C = X[:-1], Y[:-1], C[:-1]

        if self.orientation == 'vertical':
            args = (X, Y, C)
        else:
            args = (np.transpose(Y), np.transpose(X), np.transpose(C))
        kw = {'cmap':self.cmap, 'norm':self.norm,
              'shading':'flat', 'alpha':self.alpha,
              }

        del self.solids
        del self.dividers

        col = self.ax.pcolormesh(*args, **kw)

        self.solids = col
        if self.drawedges:
            self.dividers = collections.LineCollection(self._edges(X,Y),
                              colors=(mpl.rcParams['axes.edgecolor'],),
                              linewidths=(0.5*mpl.rcParams['axes.linewidth'],),
                              )
            self.ax.add_collection(self.dividers)
        else:
            self.dividers = None

    def add_lines(self, levels, colors, linewidths):
        '''
        Draw lines on the colorbar. It deletes preexisting lines.
        '''
        del self.lines

        N = len(levels)
        x = np.array([1.0, 2.0])
        X, Y = np.meshgrid(x,levels)
        if self.orientation == 'vertical':
            xy = [list(zip(X[i], Y[i])) for i in xrange(N)]
        else:
            xy = [list(zip(Y[i], X[i])) for i in xrange(N)]
        col = collections.LineCollection(xy, linewidths=linewidths,
                                         )
        self.lines = col
        col.set_color(colors)
        self.ax.add_collection(col)


    def _select_locator(self, formatter):
        '''
        select a suitable locator
        '''
        if self.boundaries is None:
            if isinstance(self.norm, colors.NoNorm):
                nv = len(self._values)
                base = 1 + int(nv/10)
                locator = ticker.IndexLocator(base=base, offset=0)
            elif isinstance(self.norm, colors.BoundaryNorm):
                b = self.norm.boundaries
                locator = ticker.FixedLocator(b, nbins=10)
            elif isinstance(self.norm, colors.LogNorm):
                locator = ticker.LogLocator()
            else:
                locator = ticker.MaxNLocator(nbins=5)
        else:
            b = self._boundaries[self._inside]
            locator = ticker.FixedLocator(b) #, nbins=10)

        self.cbar_axis.set_major_locator(locator)


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
                self._values = 0.5*(self._boundaries[:-1]
                                        + self._boundaries[1:])
                if isinstance(self.norm, colors.NoNorm):
                    self._values = (self._values + 0.00001).astype(np.int16)
                return
            self._values = np.array(self.values)
            return
        if self.values is not None:
            self._values = np.array(self.values)
            if self.boundaries is None:
                b = np.zeros(len(self.values)+1, 'd')
                b[1:-1] = 0.5*(self._values[:-1] - self._values[1:])
                b[0] = 2.0*b[1] - b[2]
                b[-1] = 2.0*b[-2] - b[-3]
                self._boundaries = b
                return
            self._boundaries = np.array(self.boundaries)
            return
        # Neither boundaries nor values are specified;
        # make reasonable ones based on cmap and norm.
        if isinstance(self.norm, colors.NoNorm):
            b = self._uniform_y(self.cmap.N+1) * self.cmap.N - 0.5
            v = np.zeros((len(b)-1,), dtype=np.int16)
            v = np.arange(self.cmap.N, dtype=np.int16)
            self._boundaries = b
            self._values = v
            return
        elif isinstance(self.norm, colors.BoundaryNorm):
            b = np.array(self.norm.boundaries)
            v = np.zeros((len(b)-1,), dtype=float)
            bi = self.norm.boundaries
            v = 0.5*(bi[:-1] + bi[1:])
            self._boundaries = b
            self._values = v
            return
        else:
            b = self._uniform_y(self.cmap.N+1)

        self._process_values(b)


    def _uniform_y(self, N):
        '''
        Return colorbar data coordinates for *N* uniformly
        spaced boundaries.
        '''
        vmin, vmax = self._get_colorbar_limits()
        if isinstance(self.norm, colors.LogNorm):
            y = np.logspace(np.log10(vmin), np.log10(vmax), N)
        else:
            y = np.linspace(vmin, vmax, N)
        return y

    def _mesh(self):
        '''
        Return X,Y, the coordinate arrays for the colorbar pcolormesh.
        These are suitable for a vertical colorbar; swapping and
        transposition for a horizontal colorbar are done outside
        this function.
        '''
        x = np.array([1.0, 2.0])
        if self.spacing == 'uniform':
            y = self._uniform_y(len(self._boundaries))
        else:
            y = self._boundaries
        self._y = y

        X, Y = np.meshgrid(x,y)
        return X, Y


    def set_alpha(self, alpha):
        """
        set alpha value.
        """
        self.alpha = alpha


class Colorbar(ColorbarBase):
    def __init__(self, ax, mappable, **kw):
        mappable.autoscale_None() # Ensure mappable.norm.vmin, vmax
                             # are set when colorbar is called,
                             # even if mappable.draw has not yet
                             # been called.  This will not change
                             # vmin, vmax if they are already set.
        self.mappable = mappable
        kw['cmap'] = mappable.cmap
        kw['norm'] = mappable.norm
        kw['alpha'] = mappable.get_alpha()
        if isinstance(mappable, contour.ContourSet):
            CS = mappable
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
            ColorbarBase.__init__(self, ax, **kw)


    def add_lines(self, CS):
        '''
        Add the lines from a non-filled
        :class:`~matplotlib.contour.ContourSet` to the colorbar.
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
        ColorbarBase.add_lines(self, CS.levels, tcolors, tlinewidths)

    def update_bruteforce(self, mappable):
        """
        Update the colorbar artists to reflect the change of the
        associated mappable.
        """
        self.update_artists()

        if isinstance(mappable, contour.ContourSet):
            if not mappable.filled:
                self.add_lines(mappable)

@docstring.Substitution(make_axes_kw_doc)
def make_axes(parent, **kw):
    '''
    Resize and reposition a parent axes, and return a child
    axes suitable for a colorbar::

        cax, kw = make_axes(parent, **kw)

    Keyword arguments may include the following (with defaults):

        *orientation*
            'vertical'  or 'horizontal'

    %s

    All but the first of these are stripped from the input kw set.

    Returns (cax, kw), the child axes and the reduced kw dictionary.
    '''
    orientation = kw.setdefault('orientation', 'vertical')
    fraction = kw.pop('fraction', 0.15)
    shrink = kw.pop('shrink', 1.0)
    aspect = kw.pop('aspect', 20)
    #pb = transforms.PBox(parent.get_position())
    pb = parent.get_position(original=True).frozen()
    if orientation == 'vertical':
        pad = kw.pop('pad', 0.05)
        x1 = 1.0-fraction
        pb1, pbx, pbcb = pb.splitx(x1-pad, x1)
        pbcb = pbcb.shrunk(1.0, shrink).anchored('C', pbcb)
        anchor = (0.0, 0.5)
        panchor = (1.0, 0.5)
    else:
        pad = kw.pop('pad', 0.15)
        pbcb, pbx, pb1 = pb.splity(fraction, fraction+pad)
        pbcb = pbcb.shrunk(shrink, 1.0).anchored('C', pbcb)
        aspect = 1.0/aspect
        anchor = (0.5, 1.0)
        panchor = (0.5, 0.0)
    parent.set_position(pb1)
    parent.set_anchor(panchor)
    fig = parent.get_figure()
    cax = fig.add_axes(pbcb)
    cax.set_aspect(aspect, anchor=anchor, adjustable='box')
    return cax, kw


def colorbar(mappable, cax=None, ax=None, **kw):
    """
    Create a colorbar for a ScalarMappable instance.

    Documentation for the pylab thin wrapper:
    %(colorbar_doc)s
    """
    import matplotlib.pyplot as plt
    if ax is None:
        ax = plt.gca()
    if cax is None:
        cax, kw = make_axes(ax, **kw)
    cax._hold = True
    cb = Colorbar(cax, mappable, **kw)

    def on_changed(m):
        cb.set_cmap(m.get_cmap())
        cb.set_clim(m.get_clim())
        cb.update_bruteforce(m)

    cbid = mappable.callbacksSM.connect('changed', on_changed)
    mappable.colorbar = cb
    ax.figure.sca(ax)
    return cb
