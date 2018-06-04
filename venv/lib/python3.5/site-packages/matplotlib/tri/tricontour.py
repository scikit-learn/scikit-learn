from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

from matplotlib.contour import ContourSet
from matplotlib.tri.triangulation import Triangulation
import matplotlib._tri as _tri
import numpy as np


class TriContourSet(ContourSet):
    """
    Create and store a set of contour lines or filled regions for
    a triangular grid.

    User-callable method: clabel

    Useful attributes:
      ax:
        the axes object in which the contours are drawn
      collections:
        a silent_list of LineCollections or PolyCollections
      levels:
        contour levels
      layers:
        same as levels for line contours; half-way between
        levels for filled contours.  See _process_colors method.
    """
    def __init__(self, ax, *args, **kwargs):
        """
        Draw triangular grid contour lines or filled regions,
        depending on whether keyword arg 'filled' is False
        (default) or True.

        The first argument of the initializer must be an axes
        object.  The remaining arguments and keyword arguments
        are described in the docstring of `tricontour`.
        """
        ContourSet.__init__(self, ax, *args, **kwargs)

    def _process_args(self, *args, **kwargs):
        """
        Process args and kwargs.
        """
        if isinstance(args[0], TriContourSet):
            C = args[0].cppContourGenerator
            if self.levels is None:
                self.levels = args[0].levels
        else:
            tri, z = self._contour_args(args, kwargs)
            C = _tri.TriContourGenerator(tri.get_cpp_triangulation(), z)
            self._mins = [tri.x.min(), tri.y.min()]
            self._maxs = [tri.x.max(), tri.y.max()]

        self.cppContourGenerator = C
        return kwargs

    def _get_allsegs_and_allkinds(self):
        """
        Create and return allsegs and allkinds by calling underlying C code.
        """
        allsegs = []
        if self.filled:
            lowers, uppers = self._get_lowers_and_uppers()
            allkinds = []
            for lower, upper in zip(lowers, uppers):
                segs, kinds = self.cppContourGenerator.create_filled_contour(
                    lower, upper)
                allsegs.append([segs])
                allkinds.append([kinds])
        else:
            allkinds = None
            for level in self.levels:
                segs = self.cppContourGenerator.create_contour(level)
                allsegs.append(segs)
        return allsegs, allkinds

    def _contour_args(self, args, kwargs):
        if self.filled:
            fn = 'contourf'
        else:
            fn = 'contour'
        tri, args, kwargs = Triangulation.get_from_args_and_kwargs(*args,
                                                                   **kwargs)
        z = np.asarray(args[0])
        if z.shape != tri.x.shape:
            raise ValueError('z array must have same length as triangulation x'
                             ' and y arrays')
        self.zmax = z.max()
        self.zmin = z.min()
        if self.logscale and self.zmin <= 0:
            raise ValueError('Cannot %s log of negative values.' % fn)
        self._contour_level_args(z, args[1:])
        return (tri, z)


def tricontour(ax, *args, **kwargs):
    """
    Draw contours on an unstructured triangular grid.
    :func:`~matplotlib.pyplot.tricontour` and
    :func:`~matplotlib.pyplot.tricontourf` draw contour lines and
    filled contours, respectively.  Except as noted, function
    signatures and return values are the same for both versions.

    The triangulation can be specified in one of two ways; either::

        tricontour(triangulation, ...)

    where triangulation is a :class:`matplotlib.tri.Triangulation`
    object, or

    ::

        tricontour(x, y, ...)
        tricontour(x, y, triangles, ...)
        tricontour(x, y, triangles=triangles, ...)
        tricontour(x, y, mask=mask, ...)
        tricontour(x, y, triangles, mask=mask, ...)

    in which case a Triangulation object will be created.  See
    :class:`~matplotlib.tri.Triangulation` for a explanation of
    these possibilities.

    The remaining arguments may be::

        tricontour(..., Z)

    where *Z* is the array of values to contour, one per point
    in the triangulation.  The level values are chosen
    automatically.

    ::

        tricontour(..., Z, N)

    contour up to *N+1* automatically chosen contour levels
    (*N* intervals).

    ::

        tricontour(..., Z, V)

    draw contour lines at the values specified in sequence *V*,
    which must be in increasing order.

    ::

        tricontourf(..., Z, V)

    fill the (len(*V*)-1) regions between the values in *V*,
    which must be in increasing order.

    ::

        tricontour(Z, **kwargs)

    Use keyword args to control colors, linewidth, origin, cmap ... see
    below for more details.

    ``C = tricontour(...)`` returns a
    :class:`~matplotlib.contour.TriContourSet` object.

    Optional keyword arguments:

        *colors*: [ *None* | string | (mpl_colors) ]
        If *None*, the colormap specified by cmap will be used.

        If a string, like 'r' or 'red', all levels will be plotted in this
        color.

        If a tuple of matplotlib color args (string, float, rgb, etc),
        different levels will be plotted in different colors in the order
        specified.

        *alpha*: float
        The alpha blending value

        *cmap*: [ *None* | Colormap ]
        A cm :class:`~matplotlib.colors.Colormap` instance or
        *None*. If *cmap* is *None* and *colors* is *None*, a
        default Colormap is used.

        *norm*: [ *None* | Normalize ]
        A :class:`matplotlib.colors.Normalize` instance for
        scaling data values to colors. If *norm* is *None* and
        *colors* is *None*, the default linear scaling is used.

        *levels* [level0, level1, ..., leveln]
        A list of floating point numbers indicating the level
        curves to draw, in increasing order; e.g., to draw just
        the zero contour pass ``levels=[0]``

        *origin*: [ *None* | 'upper' | 'lower' | 'image' ]
        If *None*, the first value of *Z* will correspond to the
        lower left corner, location (0,0). If 'image', the rc
        value for ``image.origin`` will be used.

        This keyword is not active if *X* and *Y* are specified in
        the call to contour.

        *extent*: [ *None* | (x0,x1,y0,y1) ]

        If *origin* is not *None*, then *extent* is interpreted as
        in :func:`matplotlib.pyplot.imshow`: it gives the outer
        pixel boundaries. In this case, the position of Z[0,0]
        is the center of the pixel, not a corner. If *origin* is
        *None*, then (*x0*, *y0*) is the position of Z[0,0], and
        (*x1*, *y1*) is the position of Z[-1,-1].

        This keyword is not active if *X* and *Y* are specified in
        the call to contour.

        *locator*: [ *None* | ticker.Locator subclass ]
        If *locator* is None, the default
        :class:`~matplotlib.ticker.MaxNLocator` is used. The
        locator is used to determine the contour levels if they
        are not given explicitly via the *V* argument.

        *extend*: [ 'neither' | 'both' | 'min' | 'max' ]
        Unless this is 'neither', contour levels are automatically
        added to one or both ends of the range so that all data
        are included. These added ranges are then mapped to the
        special colormap values which default to the ends of the
        colormap range, but can be set via
        :meth:`matplotlib.colors.Colormap.set_under` and
        :meth:`matplotlib.colors.Colormap.set_over` methods.

        *xunits*, *yunits*: [ *None* | registered units ]
        Override axis units by specifying an instance of a
        :class:`matplotlib.units.ConversionInterface`.


    tricontour-only keyword arguments:

        *linewidths*: [ *None* | number | tuple of numbers ]
        If *linewidths* is *None*, the default width in
        ``lines.linewidth`` in ``matplotlibrc`` is used.

        If a number, all levels will be plotted with this linewidth.

        If a tuple, different levels will be plotted with different
        linewidths in the order specified

        *linestyles*: [ *None* | 'solid' | 'dashed' | 'dashdot' | 'dotted' ]
        If *linestyles* is *None*, the 'solid' is used.

        *linestyles* can also be an iterable of the above strings
        specifying a set of linestyles to be used. If this
        iterable is shorter than the number of contour levels
        it will be repeated as necessary.

        If contour is using a monochrome colormap and the contour
        level is less than 0, then the linestyle specified
        in ``contour.negative_linestyle`` in ``matplotlibrc``
        will be used.

    tricontourf-only keyword arguments:

        *antialiased*: bool
        enable antialiasing

    Note: tricontourf fills intervals that are closed at the top; that
    is, for boundaries *z1* and *z2*, the filled region is::

        z1 < z <= z2

    There is one exception: if the lowest boundary coincides with
    the minimum value of the *z* array, then that minimum value
    will be included in the lowest interval.
    """
    if not ax._hold:
        ax.cla()
    kwargs['filled'] = False
    return TriContourSet(ax, *args, **kwargs)


def tricontourf(ax, *args, **kwargs):
    if not ax._hold:
        ax.cla()
    kwargs['filled'] = True
    return TriContourSet(ax, *args, **kwargs)
tricontourf.__doc__ = tricontour.__doc__
