"""
Support for plotting vector fields.

Presently this contains Quiver and Barb. Quiver plots an arrow in the
direction of the vector, with the size of the arrow related to the
magnitude of the vector.

Barbs are like quiver in that they point along a vector, but
the magnitude of the vector is given schematically by the presence of barbs
or flags on the barb.

This will also become a home for things such as standard
deviation ellipses, which can and will be derived very easily from
the Quiver code.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import math
import weakref

import numpy as np
from numpy import ma
import matplotlib.collections as mcollections
import matplotlib.transforms as transforms
import matplotlib.text as mtext
import matplotlib.artist as martist
from matplotlib.artist import allow_rasterization
from matplotlib import docstring
import matplotlib.font_manager as font_manager
from matplotlib.cbook import delete_masked_points
from matplotlib.patches import CirclePolygon


_quiver_doc = """
Plot a 2-D field of arrows.

Call signatures::

  quiver(U, V, **kw)
  quiver(U, V, C, **kw)
  quiver(X, Y, U, V, **kw)
  quiver(X, Y, U, V, C, **kw)

*U* and *V* are the arrow data, *X* and *Y* set the location of the
arrows, and *C* sets the color of the arrows. These arguments may be 1-D or
2-D arrays or sequences.

If *X* and *Y* are absent, they will be generated as a uniform grid.
If *U* and *V* are 2-D arrays and *X* and *Y* are 1-D, and if ``len(X)`` and
``len(Y)`` match the column and row dimensions of *U*, then *X* and *Y* will be
expanded with :func:`numpy.meshgrid`.

The default settings auto-scales the length of the arrows to a reasonable size.
To change this behavior see the *scale* and *scale_units* kwargs.

The defaults give a slightly swept-back arrow; to make the head a
triangle, make *headaxislength* the same as *headlength*. To make the
arrow more pointed, reduce *headwidth* or increase *headlength* and
*headaxislength*. To make the head smaller relative to the shaft,
scale down all the head parameters. You will probably do best to leave
minshaft alone.

*linewidths* and *edgecolors* can be used to customize the arrow
outlines.

Parameters
----------
X : 1D or 2D array, sequence, optional
    The x coordinates of the arrow locations
Y : 1D or 2D array, sequence, optional
    The y coordinates of the arrow locations
U : 1D or 2D array or masked array, sequence
    The x components of the arrow vectors
V : 1D or 2D array or masked array, sequence
    The y components of the arrow vectors
C : 1D or 2D array, sequence, optional
    The arrow colors
units : [ 'width' | 'height' | 'dots' | 'inches' | 'x' | 'y' | 'xy' ]
    The arrow dimensions (except for *length*) are measured in multiples of
    this unit.

    'width' or 'height': the width or height of the axis

    'dots' or 'inches': pixels or inches, based on the figure dpi

    'x', 'y', or 'xy': respectively *X*, *Y*, or :math:`\\sqrt{X^2 + Y^2}`
    in data units

    The arrows scale differently depending on the units.  For
    'x' or 'y', the arrows get larger as one zooms in; for other
    units, the arrow size is independent of the zoom state.  For
    'width or 'height', the arrow size increases with the width and
    height of the axes, respectively, when the window is resized;
    for 'dots' or 'inches', resizing does not change the arrows.
angles : [ 'uv' | 'xy' ], array, optional
    Method for determining the angle of the arrows. Default is 'uv'.

    'uv': the arrow axis aspect ratio is 1 so that
    if *U*==*V* the orientation of the arrow on the plot is 45 degrees
    counter-clockwise from the horizontal axis (positive to the right).

    'xy': arrows point from (x,y) to (x+u, y+v).
    Use this for plotting a gradient field, for example.

    Alternatively, arbitrary angles may be specified as an array
    of values in degrees, counter-clockwise from the horizontal axis.

    Note: inverting a data axis will correspondingly invert the
    arrows only with ``angles='xy'``.
scale : None, float, optional
    Number of data units per arrow length unit, e.g., m/s per plot width; a
    smaller scale parameter makes the arrow longer. Default is *None*.

    If *None*, a simple autoscaling algorithm is used, based on the average
    vector length and the number of vectors. The arrow length unit is given by
    the *scale_units* parameter
scale_units : [ 'width' | 'height' | 'dots' | 'inches' | 'x' | 'y' | 'xy' ], \
None, optional
    If the *scale* kwarg is *None*, the arrow length unit. Default is *None*.

    e.g. *scale_units* is 'inches', *scale* is 2.0, and
    ``(u,v) = (1,0)``, then the vector will be 0.5 inches long.

    If *scale_units* is 'width'/'height', then the vector will be half the
    width/height of the axes.

    If *scale_units* is 'x' then the vector will be 0.5 x-axis
    units. To plot vectors in the x-y plane, with u and v having
    the same units as x and y, use
    ``angles='xy', scale_units='xy', scale=1``.
width : scalar, optional
    Shaft width in arrow units; default depends on choice of units,
    above, and number of vectors; a typical starting value is about
    0.005 times the width of the plot.
headwidth : scalar, optional
    Head width as multiple of shaft width, default is 3
headlength : scalar, optional
    Head length as multiple of shaft width, default is 5
headaxislength : scalar, optional
    Head length at shaft intersection, default is 4.5
minshaft : scalar, optional
    Length below which arrow scales, in units of head length. Do not
    set this to less than 1, or small arrows will look terrible!
    Default is 1
minlength : scalar, optional
    Minimum length as a multiple of shaft width; if an arrow length
    is less than this, plot a dot (hexagon) of this diameter instead.
    Default is 1.
pivot : [ 'tail' | 'mid' | 'middle' | 'tip' ], optional
    The part of the arrow that is at the grid point; the arrow rotates
    about this point, hence the name *pivot*.
color : [ color | color sequence ], optional
    This is a synonym for the
    :class:`~matplotlib.collections.PolyCollection` facecolor kwarg.
    If *C* has been set, *color* has no effect.

Notes
-----
Additional :class:`~matplotlib.collections.PolyCollection`
keyword arguments:

%(PolyCollection)s

See Also
--------
quiverkey : Add a key to a quiver plot
""" % docstring.interpd.params

_quiverkey_doc = """
Add a key to a quiver plot.

Call signature::

  quiverkey(Q, X, Y, U, label, **kw)

Arguments:

  *Q*:
    The Quiver instance returned by a call to quiver.

  *X*, *Y*:
    The location of the key; additional explanation follows.

  *U*:
    The length of the key

  *label*:
    A string with the length and units of the key

Keyword arguments:

  *angle* = 0
    The angle of the key arrow. Measured in degrees anti-clockwise from the
    x-axis.

  *coordinates* = [ 'axes' | 'figure' | 'data' | 'inches' ]
    Coordinate system and units for *X*, *Y*: 'axes' and 'figure' are
    normalized coordinate systems with 0,0 in the lower left and 1,1
    in the upper right; 'data' are the axes data coordinates (used for
    the locations of the vectors in the quiver plot itself); 'inches'
    is position in the figure in inches, with 0,0 at the lower left
    corner.

  *color*:
    overrides face and edge colors from *Q*.

  *labelpos* = [ 'N' | 'S' | 'E' | 'W' ]
    Position the label above, below, to the right, to the left of the
    arrow, respectively.

  *labelsep*:
    Distance in inches between the arrow and the label.  Default is
    0.1

  *labelcolor*:
    defaults to default :class:`~matplotlib.text.Text` color.

  *fontproperties*:
    A dictionary with keyword arguments accepted by the
    :class:`~matplotlib.font_manager.FontProperties` initializer:
    *family*, *style*, *variant*, *size*, *weight*

Any additional keyword arguments are used to override vector
properties taken from *Q*.

The positioning of the key depends on *X*, *Y*, *coordinates*, and
*labelpos*.  If *labelpos* is 'N' or 'S', *X*, *Y* give the position
of the middle of the key arrow.  If *labelpos* is 'E', *X*, *Y*
positions the head, and if *labelpos* is 'W', *X*, *Y* positions the
tail; in either of these two cases, *X*, *Y* is somewhere in the
middle of the arrow+label key object.
"""


class QuiverKey(martist.Artist):
    """ Labelled arrow for use as a quiver plot scale key."""
    halign = {'N': 'center', 'S': 'center', 'E': 'left', 'W': 'right'}
    valign = {'N': 'bottom', 'S': 'top', 'E': 'center', 'W': 'center'}
    pivot = {'N': 'middle', 'S': 'middle', 'E': 'tip', 'W': 'tail'}

    def __init__(self, Q, X, Y, U, label, **kw):
        martist.Artist.__init__(self)
        self.Q = Q
        self.X = X
        self.Y = Y
        self.U = U
        self.angle = kw.pop('angle', 0)
        self.coord = kw.pop('coordinates', 'axes')
        self.color = kw.pop('color', None)
        self.label = label
        self._labelsep_inches = kw.pop('labelsep', 0.1)
        self.labelsep = (self._labelsep_inches * Q.ax.figure.dpi)

        # try to prevent closure over the real self
        weak_self = weakref.ref(self)

        def on_dpi_change(fig):
            self_weakref = weak_self()
            if self_weakref is not None:
                self_weakref.labelsep = (self_weakref._labelsep_inches*fig.dpi)
                self_weakref._initialized = False  # simple brute force update
                                                   # works because _init is
                                                   # called at the start of
                                                   # draw.

        self._cid = Q.ax.figure.callbacks.connect('dpi_changed',
                                                  on_dpi_change)

        self.labelpos = kw.pop('labelpos', 'N')
        self.labelcolor = kw.pop('labelcolor', None)
        self.fontproperties = kw.pop('fontproperties', dict())
        self.kw = kw
        _fp = self.fontproperties
        # boxprops = dict(facecolor='red')
        self.text = mtext.Text(
                        text=label,  # bbox=boxprops,
                        horizontalalignment=self.halign[self.labelpos],
                        verticalalignment=self.valign[self.labelpos],
                        fontproperties=font_manager.FontProperties(**_fp))

        if self.labelcolor is not None:
            self.text.set_color(self.labelcolor)
        self._initialized = False
        self.zorder = Q.zorder + 0.1

    def remove(self):
        """
        Overload the remove method
        """
        self.Q.ax.figure.callbacks.disconnect(self._cid)
        self._cid = None
        # pass the remove call up the stack
        martist.Artist.remove(self)

    __init__.__doc__ = _quiverkey_doc

    def _init(self):
        if True:  # not self._initialized:
            if not self.Q._initialized:
                self.Q._init()
            self._set_transform()
            _pivot = self.Q.pivot
            self.Q.pivot = self.pivot[self.labelpos]
            # Hack: save and restore the Umask
            _mask = self.Q.Umask
            self.Q.Umask = ma.nomask
            self.verts = self.Q._make_verts(np.array([self.U]),
                                            np.zeros((1,)),
                                            self.angle)
            self.Q.Umask = _mask
            self.Q.pivot = _pivot
            kw = self.Q.polykw
            kw.update(self.kw)
            self.vector = mcollections.PolyCollection(
                                        self.verts,
                                        offsets=[(self.X, self.Y)],
                                        transOffset=self.get_transform(),
                                        **kw)
            if self.color is not None:
                self.vector.set_color(self.color)
            self.vector.set_transform(self.Q.get_transform())
            self.vector.set_figure(self.get_figure())
            self._initialized = True

    def _text_x(self, x):
        if self.labelpos == 'E':
            return x + self.labelsep
        elif self.labelpos == 'W':
            return x - self.labelsep
        else:
            return x

    def _text_y(self, y):
        if self.labelpos == 'N':
            return y + self.labelsep
        elif self.labelpos == 'S':
            return y - self.labelsep
        else:
            return y

    @allow_rasterization
    def draw(self, renderer):
        self._init()
        self.vector.draw(renderer)
        x, y = self.get_transform().transform_point((self.X, self.Y))
        self.text.set_x(self._text_x(x))
        self.text.set_y(self._text_y(y))
        self.text.draw(renderer)
        self.stale = False

    def _set_transform(self):
        if self.coord == 'data':
            self.set_transform(self.Q.ax.transData)
        elif self.coord == 'axes':
            self.set_transform(self.Q.ax.transAxes)
        elif self.coord == 'figure':
            self.set_transform(self.Q.ax.figure.transFigure)
        elif self.coord == 'inches':
            self.set_transform(self.Q.ax.figure.dpi_scale_trans)
        else:
            raise ValueError('unrecognized coordinates')

    def set_figure(self, fig):
        martist.Artist.set_figure(self, fig)
        self.text.set_figure(fig)

    def contains(self, mouseevent):
        # Maybe the dictionary should allow one to
        # distinguish between a text hit and a vector hit.
        if (self.text.contains(mouseevent)[0] or
                self.vector.contains(mouseevent)[0]):
            return True, {}
        return False, {}

    quiverkey_doc = _quiverkey_doc


# This is a helper function that parses out the various combination of
# arguments for doing colored vector plots.  Pulling it out here
# allows both Quiver and Barbs to use it
def _parse_args(*args):
    X, Y, U, V, C = [None] * 5
    args = list(args)

    # The use of atleast_1d allows for handling scalar arguments while also
    # keeping masked arrays
    if len(args) == 3 or len(args) == 5:
        C = np.atleast_1d(args.pop(-1))
    V = np.atleast_1d(args.pop(-1))
    U = np.atleast_1d(args.pop(-1))
    if U.ndim == 1:
        nr, nc = 1, U.shape[0]
    else:
        nr, nc = U.shape
    if len(args) == 2:  # remaining after removing U,V,C
        X, Y = [np.array(a).ravel() for a in args]
        if len(X) == nc and len(Y) == nr:
            X, Y = [a.ravel() for a in np.meshgrid(X, Y)]
    else:
        indexgrid = np.meshgrid(np.arange(nc), np.arange(nr))
        X, Y = [np.ravel(a) for a in indexgrid]
    return X, Y, U, V, C


def _check_consistent_shapes(*arrays):
    all_shapes = set(a.shape for a in arrays)
    if len(all_shapes) != 1:
        raise ValueError('The shapes of the passed in arrays do not match.')


class Quiver(mcollections.PolyCollection):
    """
    Specialized PolyCollection for arrows.

    The only API method is set_UVC(), which can be used
    to change the size, orientation, and color of the
    arrows; their locations are fixed when the class is
    instantiated.  Possibly this method will be useful
    in animations.

    Much of the work in this class is done in the draw()
    method so that as much information as possible is available
    about the plot.  In subsequent draw() calls, recalculation
    is limited to things that might have changed, so there
    should be no performance penalty from putting the calculations
    in the draw() method.
    """

    _PIVOT_VALS = ('tail', 'mid', 'middle', 'tip')

    @docstring.Substitution(_quiver_doc)
    def __init__(self, ax, *args, **kw):
        """
        The constructor takes one required argument, an Axes
        instance, followed by the args and kwargs described
        by the following pylab interface documentation:
        %s
        """
        self.ax = ax
        X, Y, U, V, C = _parse_args(*args)
        self.X = X
        self.Y = Y
        self.XY = np.hstack((X[:, np.newaxis], Y[:, np.newaxis]))
        self.N = len(X)
        self.scale = kw.pop('scale', None)
        self.headwidth = kw.pop('headwidth', 3)
        self.headlength = float(kw.pop('headlength', 5))
        self.headaxislength = kw.pop('headaxislength', 4.5)
        self.minshaft = kw.pop('minshaft', 1)
        self.minlength = kw.pop('minlength', 1)
        self.units = kw.pop('units', 'width')
        self.scale_units = kw.pop('scale_units', None)
        self.angles = kw.pop('angles', 'uv')
        self.width = kw.pop('width', None)
        self.color = kw.pop('color', 'k')

        pivot = kw.pop('pivot', 'tail').lower()
        # validate pivot
        if pivot not in self._PIVOT_VALS:
            raise ValueError(
                'pivot must be one of {keys}, you passed {inp}'.format(
                      keys=self._PIVOT_VALS, inp=pivot))
        # normalize to 'middle'
        if pivot == 'mid':
            pivot = 'middle'
        self.pivot = pivot

        self.transform = kw.pop('transform', ax.transData)
        kw.setdefault('facecolors', self.color)
        kw.setdefault('linewidths', (0,))
        mcollections.PolyCollection.__init__(self, [], offsets=self.XY,
                                             transOffset=self.transform,
                                             closed=False,
                                             **kw)
        self.polykw = kw
        self.set_UVC(U, V, C)
        self._initialized = False

        self.keyvec = None
        self.keytext = None

        # try to prevent closure over the real self
        weak_self = weakref.ref(self)

        def on_dpi_change(fig):
            self_weakref = weak_self()
            if self_weakref is not None:
                self_weakref._new_UV = True  # vertices depend on width, span
                                             # which in turn depend on dpi
                self_weakref._initialized = False  # simple brute force update
                                                   # works because _init is
                                                   # called at the start of
                                                   # draw.

        self._cid = self.ax.figure.callbacks.connect('dpi_changed',
                                                     on_dpi_change)

    def remove(self):
        """
        Overload the remove method
        """
        # disconnect the call back
        self.ax.figure.callbacks.disconnect(self._cid)
        self._cid = None
        # pass the remove call up the stack
        mcollections.PolyCollection.remove(self)

    def _init(self):
        """
        Initialization delayed until first draw;
        allow time for axes setup.
        """
        # It seems that there are not enough event notifications
        # available to have this work on an as-needed basis at present.
        if True:  # not self._initialized:
            trans = self._set_transform()
            ax = self.ax
            sx, sy = trans.inverted().transform_point(
                                            (ax.bbox.width, ax.bbox.height))
            self.span = sx
            if self.width is None:
                sn = np.clip(math.sqrt(self.N), 8, 25)
                self.width = 0.06 * self.span / sn

            # _make_verts sets self.scale if not already specified
            if not self._initialized and self.scale is None:
                self._make_verts(self.U, self.V, self.angles)

            self._initialized = True

    def get_datalim(self, transData):
        trans = self.get_transform()
        transOffset = self.get_offset_transform()
        full_transform = (trans - transData) + (transOffset - transData)
        XY = full_transform.transform(self.XY)
        bbox = transforms.Bbox.null()
        bbox.update_from_data_xy(XY, ignore=True)
        return bbox

    @allow_rasterization
    def draw(self, renderer):
        self._init()
        verts = self._make_verts(self.U, self.V, self.angles)
        self.set_verts(verts, closed=False)
        self._new_UV = False
        mcollections.PolyCollection.draw(self, renderer)
        self.stale = False

    def set_UVC(self, U, V, C=None):
        # We need to ensure we have a copy, not a reference
        # to an array that might change before draw().
        U = ma.masked_invalid(U, copy=True).ravel()
        V = ma.masked_invalid(V, copy=True).ravel()
        mask = ma.mask_or(U.mask, V.mask, copy=False, shrink=True)
        if C is not None:
            C = ma.masked_invalid(C, copy=True).ravel()
            mask = ma.mask_or(mask, C.mask, copy=False, shrink=True)
            if mask is ma.nomask:
                C = C.filled()
            else:
                C = ma.array(C, mask=mask, copy=False)
        self.U = U.filled(1)
        self.V = V.filled(1)
        self.Umask = mask
        if C is not None:
            self.set_array(C)
        self._new_UV = True
        self.stale = True

    def _dots_per_unit(self, units):
        """
        Return a scale factor for converting from units to pixels
        """
        ax = self.ax
        if units in ('x', 'y', 'xy'):
            if units == 'x':
                dx0 = ax.viewLim.width
                dx1 = ax.bbox.width
            elif units == 'y':
                dx0 = ax.viewLim.height
                dx1 = ax.bbox.height
            else:  # 'xy' is assumed
                dxx0 = ax.viewLim.width
                dxx1 = ax.bbox.width
                dyy0 = ax.viewLim.height
                dyy1 = ax.bbox.height
                dx1 = np.hypot(dxx1, dyy1)
                dx0 = np.hypot(dxx0, dyy0)
            dx = dx1 / dx0
        else:
            if units == 'width':
                dx = ax.bbox.width
            elif units == 'height':
                dx = ax.bbox.height
            elif units == 'dots':
                dx = 1.0
            elif units == 'inches':
                dx = ax.figure.dpi
            else:
                raise ValueError('unrecognized units')
        return dx

    def _set_transform(self):
        """
        Sets the PolygonCollection transform to go
        from arrow width units to pixels.
        """
        dx = self._dots_per_unit(self.units)
        self._trans_scale = dx  # pixels per arrow width unit
        trans = transforms.Affine2D().scale(dx)
        self.set_transform(trans)
        return trans

    def _angles_lengths(self, U, V, eps=1):
        xy = self.ax.transData.transform(self.XY)
        uv = np.hstack((U[:, np.newaxis], V[:, np.newaxis]))
        xyp = self.ax.transData.transform(self.XY + eps * uv)
        dxy = xyp - xy
        angles = np.arctan2(dxy[:, 1], dxy[:, 0])
        lengths = np.hypot(*dxy.T) / eps
        return angles, lengths

    def _make_verts(self, U, V, angles):
        uv = (U + V * 1j)
        str_angles = angles if isinstance(angles, six.string_types) else ''
        if str_angles == 'xy' and self.scale_units == 'xy':
            # Here eps is 1 so that if we get U, V by diffing
            # the X, Y arrays, the vectors will connect the
            # points, regardless of the axis scaling (including log).
            angles, lengths = self._angles_lengths(U, V, eps=1)
        elif str_angles == 'xy' or self.scale_units == 'xy':
            # Calculate eps based on the extents of the plot
            # so that we don't end up with roundoff error from
            # adding a small number to a large.
            eps = np.abs(self.ax.dataLim.extents).max() * 0.001
            angles, lengths = self._angles_lengths(U, V, eps=eps)
        if str_angles and self.scale_units == 'xy':
            a = lengths
        else:
            a = np.abs(uv)
        if self.scale is None:
            sn = max(10, math.sqrt(self.N))
            if self.Umask is not ma.nomask:
                amean = a[~self.Umask].mean()
            else:
                amean = a.mean()
            # crude auto-scaling
            # scale is typical arrow length as a multiple of the arrow width
            scale = 1.8 * amean * sn / self.span
        if self.scale_units is None:
            if self.scale is None:
                self.scale = scale
            widthu_per_lenu = 1.0
        else:
            if self.scale_units == 'xy':
                dx = 1
            else:
                dx = self._dots_per_unit(self.scale_units)
            widthu_per_lenu = dx / self._trans_scale
            if self.scale is None:
                self.scale = scale * widthu_per_lenu
        length = a * (widthu_per_lenu / (self.scale * self.width))
        X, Y = self._h_arrows(length)
        if str_angles == 'xy':
            theta = angles
        elif str_angles == 'uv':
            theta = np.angle(uv)
        else:
            theta = ma.masked_invalid(np.deg2rad(angles)).filled(0)
        theta = theta.reshape((-1, 1))  # for broadcasting
        xy = (X + Y * 1j) * np.exp(1j * theta) * self.width
        xy = xy[:, :, np.newaxis]
        XY = np.concatenate((xy.real, xy.imag), axis=2)
        if self.Umask is not ma.nomask:
            XY = ma.array(XY)
            XY[self.Umask] = ma.masked
            # This might be handled more efficiently with nans, given
            # that nans will end up in the paths anyway.

        return XY

    def _h_arrows(self, length):
        """ length is in arrow width units """
        # It might be possible to streamline the code
        # and speed it up a bit by using complex (x,y)
        # instead of separate arrays; but any gain would be slight.
        minsh = self.minshaft * self.headlength
        N = len(length)
        length = length.reshape(N, 1)
        # This number is chosen based on when pixel values overflow in Agg
        # causing rendering errors
        # length = np.minimum(length, 2 ** 16)
        np.clip(length, 0, 2 ** 16, out=length)
        # x, y: normal horizontal arrow
        x = np.array([0, -self.headaxislength,
                      -self.headlength, 0],
                     np.float64)
        x = x + np.array([0, 1, 1, 1]) * length
        y = 0.5 * np.array([1, 1, self.headwidth, 0], np.float64)
        y = np.repeat(y[np.newaxis, :], N, axis=0)
        # x0, y0: arrow without shaft, for short vectors
        x0 = np.array([0, minsh - self.headaxislength,
                       minsh - self.headlength, minsh], np.float64)
        y0 = 0.5 * np.array([1, 1, self.headwidth, 0], np.float64)
        ii = [0, 1, 2, 3, 2, 1, 0, 0]
        X = x.take(ii, 1)
        Y = y.take(ii, 1)
        Y[:, 3:-1] *= -1
        X0 = x0.take(ii)
        Y0 = y0.take(ii)
        Y0[3:-1] *= -1
        shrink = length / minsh if minsh != 0. else 0.
        X0 = shrink * X0[np.newaxis, :]
        Y0 = shrink * Y0[np.newaxis, :]
        short = np.repeat(length < minsh, 8, axis=1)
        # Now select X0, Y0 if short, otherwise X, Y
        np.copyto(X, X0, where=short)
        np.copyto(Y, Y0, where=short)
        if self.pivot == 'middle':
            X -= 0.5 * X[:, 3, np.newaxis]
        elif self.pivot == 'tip':
            X = X - X[:, 3, np.newaxis]   # numpy bug? using -= does not
                                          # work here unless we multiply
                                          # by a float first, as with 'mid'.
        elif self.pivot != 'tail':
            raise ValueError(("Quiver.pivot must have value in {{'middle', "
                              "'tip', 'tail'}} not {0}").format(self.pivot))

        tooshort = length < self.minlength
        if tooshort.any():
            # Use a heptagonal dot:
            th = np.arange(0, 8, 1, np.float64) * (np.pi / 3.0)
            x1 = np.cos(th) * self.minlength * 0.5
            y1 = np.sin(th) * self.minlength * 0.5
            X1 = np.repeat(x1[np.newaxis, :], N, axis=0)
            Y1 = np.repeat(y1[np.newaxis, :], N, axis=0)
            tooshort = np.repeat(tooshort, 8, 1)
            np.copyto(X, X1, where=tooshort)
            np.copyto(Y, Y1, where=tooshort)
        # Mask handling is deferred to the caller, _make_verts.
        return X, Y

    quiver_doc = _quiver_doc


_barbs_doc = r"""
Plot a 2-D field of barbs.

Call signatures::

  barb(U, V, **kw)
  barb(U, V, C, **kw)
  barb(X, Y, U, V, **kw)
  barb(X, Y, U, V, C, **kw)

Arguments:

  *X*, *Y*:
    The x and y coordinates of the barb locations
    (default is head of barb; see *pivot* kwarg)

  *U*, *V*:
    Give the x and y components of the barb shaft

  *C*:
    An optional array used to map colors to the barbs

All arguments may be 1-D or 2-D arrays or sequences. If *X* and *Y*
are absent, they will be generated as a uniform grid.  If *U* and *V*
are 2-D arrays but *X* and *Y* are 1-D, and if ``len(X)`` and ``len(Y)``
match the column and row dimensions of *U*, then *X* and *Y* will be
expanded with :func:`numpy.meshgrid`.

*U*, *V*, *C* may be masked arrays, but masked *X*, *Y* are not
supported at present.

Keyword arguments:

  *length*:
    Length of the barb in points; the other parts of the barb
    are scaled against this.
    Default is 7.

  *pivot*: [ 'tip' | 'middle' | float ]
    The part of the arrow that is at the grid point; the arrow rotates
    about this point, hence the name *pivot*.  Default is 'tip'. Can
    also be a number, which shifts the start of the barb that many
    points from the origin.

  *barbcolor*: [ color | color sequence ]
    Specifies the color all parts of the barb except any flags.  This
    parameter is analogous to the *edgecolor* parameter for polygons,
    which can be used instead. However this parameter will override
    facecolor.

  *flagcolor*: [ color | color sequence ]
    Specifies the color of any flags on the barb.  This parameter is
    analogous to the *facecolor* parameter for polygons, which can be
    used instead. However this parameter will override facecolor.  If
    this is not set (and *C* has not either) then *flagcolor* will be
    set to match *barbcolor* so that the barb has a uniform color. If
    *C* has been set, *flagcolor* has no effect.

  *sizes*:
    A dictionary of coefficients specifying the ratio of a given
    feature to the length of the barb. Only those values one wishes to
    override need to be included.  These features include:

        - 'spacing' - space between features (flags, full/half barbs)

        - 'height' - height (distance from shaft to top) of a flag or
          full barb

        - 'width' - width of a flag, twice the width of a full barb

        - 'emptybarb' - radius of the circle used for low magnitudes

  *fill_empty*:
    A flag on whether the empty barbs (circles) that are drawn should
    be filled with the flag color.  If they are not filled, they will
    be drawn such that no color is applied to the center.  Default is
    False

  *rounding*:
    A flag to indicate whether the vector magnitude should be rounded
    when allocating barb components.  If True, the magnitude is
    rounded to the nearest multiple of the half-barb increment.  If
    False, the magnitude is simply truncated to the next lowest
    multiple.  Default is True

  *barb_increments*:
    A dictionary of increments specifying values to associate with
    different parts of the barb. Only those values one wishes to
    override need to be included.

        - 'half' - half barbs (Default is 5)

        - 'full' - full barbs (Default is 10)

        - 'flag' - flags (default is 50)

  *flip_barb*:
    Either a single boolean flag or an array of booleans.  Single
    boolean indicates whether the lines and flags should point
    opposite to normal for all barbs.  An array (which should be the
    same size as the other data arrays) indicates whether to flip for
    each individual barb.  Normal behavior is for the barbs and lines
    to point right (comes from wind barbs having these features point
    towards low pressure in the Northern Hemisphere.)  Default is
    False

Barbs are traditionally used in meteorology as a way to plot the speed
and direction of wind observations, but can technically be used to
plot any two dimensional vector quantity.  As opposed to arrows, which
give vector magnitude by the length of the arrow, the barbs give more
quantitative information about the vector magnitude by putting slanted
lines or a triangle for various increments in magnitude, as show
schematically below::

 :     /\    \\
 :    /  \    \\
 :   /    \    \    \\
 :  /      \    \    \\
 : ------------------------------

.. note the double \\ at the end of each line to make the figure
.. render correctly

The largest increment is given by a triangle (or "flag"). After those
come full lines (barbs). The smallest increment is a half line.  There
is only, of course, ever at most 1 half line.  If the magnitude is
small and only needs a single half-line and no full lines or
triangles, the half-line is offset from the end of the barb so that it
can be easily distinguished from barbs with a single full line.  The
magnitude for the barb shown above would nominally be 65, using the
standard increments of 50, 10, and 5.

linewidths and edgecolors can be used to customize the barb.
Additional :class:`~matplotlib.collections.PolyCollection` keyword
arguments:

%(PolyCollection)s
""" % docstring.interpd.params

docstring.interpd.update(barbs_doc=_barbs_doc)


class Barbs(mcollections.PolyCollection):
    '''
    Specialized PolyCollection for barbs.

    The only API method is :meth:`set_UVC`, which can be used to
    change the size, orientation, and color of the arrows.  Locations
    are changed using the :meth:`set_offsets` collection method.
    Possibly this method will be useful in animations.

    There is one internal function :meth:`_find_tails` which finds
    exactly what should be put on the barb given the vector magnitude.
    From there :meth:`_make_barbs` is used to find the vertices of the
    polygon to represent the barb based on this information.
    '''
    # This may be an abuse of polygons here to render what is essentially maybe
    # 1 triangle and a series of lines.  It works fine as far as I can tell
    # however.
    @docstring.interpd
    def __init__(self, ax, *args, **kw):
        """
        The constructor takes one required argument, an Axes
        instance, followed by the args and kwargs described
        by the following pylab interface documentation:
        %(barbs_doc)s
        """
        self._pivot = kw.pop('pivot', 'tip')
        self._length = kw.pop('length', 7)
        barbcolor = kw.pop('barbcolor', None)
        flagcolor = kw.pop('flagcolor', None)
        self.sizes = kw.pop('sizes', dict())
        self.fill_empty = kw.pop('fill_empty', False)
        self.barb_increments = kw.pop('barb_increments', dict())
        self.rounding = kw.pop('rounding', True)
        self.flip = kw.pop('flip_barb', False)
        transform = kw.pop('transform', ax.transData)

        # Flagcolor and barbcolor provide convenience parameters for
        # setting the facecolor and edgecolor, respectively, of the barb
        # polygon.  We also work here to make the flag the same color as the
        # rest of the barb by default

        if None in (barbcolor, flagcolor):
            kw['edgecolors'] = 'face'
            if flagcolor:
                kw['facecolors'] = flagcolor
            elif barbcolor:
                kw['facecolors'] = barbcolor
            else:
                # Set to facecolor passed in or default to black
                kw.setdefault('facecolors', 'k')
        else:
            kw['edgecolors'] = barbcolor
            kw['facecolors'] = flagcolor

        # Explicitly set a line width if we're not given one, otherwise
        # polygons are not outlined and we get no barbs
        if 'linewidth' not in kw and 'lw' not in kw:
            kw['linewidth'] = 1

        # Parse out the data arrays from the various configurations supported
        x, y, u, v, c = _parse_args(*args)
        self.x = x
        self.y = y
        xy = np.hstack((x[:, np.newaxis], y[:, np.newaxis]))

        # Make a collection
        barb_size = self._length ** 2 / 4  # Empirically determined
        mcollections.PolyCollection.__init__(self, [], (barb_size,),
                                             offsets=xy,
                                             transOffset=transform, **kw)
        self.set_transform(transforms.IdentityTransform())

        self.set_UVC(u, v, c)

    def _find_tails(self, mag, rounding=True, half=5, full=10, flag=50):
        '''
        Find how many of each of the tail pieces is necessary.  Flag
        specifies the increment for a flag, barb for a full barb, and half for
        half a barb. Mag should be the magnitude of a vector (i.e., >= 0).

        This returns a tuple of:

            (*number of flags*, *number of barbs*, *half_flag*, *empty_flag*)

        *half_flag* is a boolean whether half of a barb is needed,
        since there should only ever be one half on a given
        barb. *empty_flag* flag is an array of flags to easily tell if
        a barb is empty (too low to plot any barbs/flags.
        '''

        # If rounding, round to the nearest multiple of half, the smallest
        # increment
        if rounding:
            mag = half * (mag / half + 0.5).astype(int)

        num_flags = np.floor(mag / flag).astype(int)
        mag = np.mod(mag, flag)

        num_barb = np.floor(mag / full).astype(int)
        mag = np.mod(mag, full)

        half_flag = mag >= half
        empty_flag = ~(half_flag | (num_flags > 0) | (num_barb > 0))

        return num_flags, num_barb, half_flag, empty_flag

    def _make_barbs(self, u, v, nflags, nbarbs, half_barb, empty_flag, length,
                    pivot, sizes, fill_empty, flip):
        '''
        This function actually creates the wind barbs.  *u* and *v*
        are components of the vector in the *x* and *y* directions,
        respectively.

        *nflags*, *nbarbs*, and *half_barb*, empty_flag* are,
        *respectively, the number of flags, number of barbs, flag for
        *half a barb, and flag for empty barb, ostensibly obtained
        *from :meth:`_find_tails`.

        *length* is the length of the barb staff in points.

        *pivot* specifies the point on the barb around which the
        entire barb should be rotated.  Right now, valid options are
        'tip' and 'middle'. Can also be a number, which shifts the start
        of the barb that many points from the origin.

        *sizes* is a dictionary of coefficients specifying the ratio
        of a given feature to the length of the barb. These features
        include:

            - *spacing*: space between features (flags, full/half
               barbs)

            - *height*: distance from shaft of top of a flag or full
               barb

            - *width* - width of a flag, twice the width of a full barb

            - *emptybarb* - radius of the circle used for low
               magnitudes

        *fill_empty* specifies whether the circle representing an
        empty barb should be filled or not (this changes the drawing
        of the polygon).

        *flip* is a flag indicating whether the features should be flipped to
        the other side of the barb (useful for winds in the southern
        hemisphere).

        This function returns list of arrays of vertices, defining a polygon
        for each of the wind barbs.  These polygons have been rotated to
        properly align with the vector direction.
        '''

        # These control the spacing and size of barb elements relative to the
        # length of the shaft
        spacing = length * sizes.get('spacing', 0.125)
        full_height = length * sizes.get('height', 0.4)
        full_width = length * sizes.get('width', 0.25)
        empty_rad = length * sizes.get('emptybarb', 0.15)

        # Controls y point where to pivot the barb.
        pivot_points = dict(tip=0.0, middle=-length / 2.)

        # Check for flip
        if flip:
            full_height = -full_height

        endx = 0.0
        try:
            endy = float(pivot)
        except ValueError:
            endy = pivot_points[pivot.lower()]

        # Get the appropriate angle for the vector components.  The offset is
        # due to the way the barb is initially drawn, going down the y-axis.
        # This makes sense in a meteorological mode of thinking since there 0
        # degrees corresponds to north (the y-axis traditionally)
        angles = -(ma.arctan2(v, u) + np.pi / 2)

        # Used for low magnitude.  We just get the vertices, so if we make it
        # out here, it can be reused.  The center set here should put the
        # center of the circle at the location(offset), rather than at the
        # same point as the barb pivot; this seems more sensible.
        circ = CirclePolygon((0, 0), radius=empty_rad).get_verts()
        if fill_empty:
            empty_barb = circ
        else:
            # If we don't want the empty one filled, we make a degenerate
            # polygon that wraps back over itself
            empty_barb = np.concatenate((circ, circ[::-1]))

        barb_list = []
        for index, angle in np.ndenumerate(angles):
            # If the vector magnitude is too weak to draw anything, plot an
            # empty circle instead
            if empty_flag[index]:
                # We can skip the transform since the circle has no preferred
                # orientation
                barb_list.append(empty_barb)
                continue

            poly_verts = [(endx, endy)]
            offset = length

            # Add vertices for each flag
            for i in range(nflags[index]):
                # The spacing that works for the barbs is a little to much for
                # the flags, but this only occurs when we have more than 1
                # flag.
                if offset != length:
                    offset += spacing / 2.
                poly_verts.extend(
                    [[endx, endy + offset],
                     [endx + full_height, endy - full_width / 2 + offset],
                     [endx, endy - full_width + offset]])

                offset -= full_width + spacing

            # Add vertices for each barb.  These really are lines, but works
            # great adding 3 vertices that basically pull the polygon out and
            # back down the line
            for i in range(nbarbs[index]):
                poly_verts.extend(
                    [(endx, endy + offset),
                     (endx + full_height, endy + offset + full_width / 2),
                     (endx, endy + offset)])

                offset -= spacing

            # Add the vertices for half a barb, if needed
            if half_barb[index]:
                # If the half barb is the first on the staff, traditionally it
                # is offset from the end to make it easy to distinguish from a
                # barb with a full one
                if offset == length:
                    poly_verts.append((endx, endy + offset))
                    offset -= 1.5 * spacing
                poly_verts.extend(
                    [(endx, endy + offset),
                     (endx + full_height / 2, endy + offset + full_width / 4),
                     (endx, endy + offset)])

            # Rotate the barb according the angle. Making the barb first and
            # then rotating it made the math for drawing the barb really easy.
            # Also, the transform framework makes doing the rotation simple.
            poly_verts = transforms.Affine2D().rotate(-angle).transform(
                poly_verts)
            barb_list.append(poly_verts)

        return barb_list

    def set_UVC(self, U, V, C=None):
        self.u = ma.masked_invalid(U, copy=False).ravel()
        self.v = ma.masked_invalid(V, copy=False).ravel()
        if C is not None:
            c = ma.masked_invalid(C, copy=False).ravel()
            x, y, u, v, c = delete_masked_points(self.x.ravel(),
                                                 self.y.ravel(),
                                                 self.u, self.v, c)
            _check_consistent_shapes(x, y, u, v, c)
        else:
            x, y, u, v = delete_masked_points(self.x.ravel(), self.y.ravel(),
                                              self.u, self.v)
            _check_consistent_shapes(x, y, u, v)

        magnitude = np.hypot(u, v)
        flags, barbs, halves, empty = self._find_tails(magnitude,
                                                       self.rounding,
                                                       **self.barb_increments)

        # Get the vertices for each of the barbs

        plot_barbs = self._make_barbs(u, v, flags, barbs, halves, empty,
                                      self._length, self._pivot, self.sizes,
                                      self.fill_empty, self.flip)
        self.set_verts(plot_barbs)

        # Set the color array
        if C is not None:
            self.set_array(c)

        # Update the offsets in case the masked data changed
        xy = np.hstack((x[:, np.newaxis], y[:, np.newaxis]))
        self._offsets = xy
        self.stale = True

    def set_offsets(self, xy):
        """
        Set the offsets for the barb polygons.  This saves the offsets passed
        in and actually sets version masked as appropriate for the existing
        U/V data. *offsets* should be a sequence.

        ACCEPTS: sequence of pairs of floats
        """
        self.x = xy[:, 0]
        self.y = xy[:, 1]
        x, y, u, v = delete_masked_points(self.x.ravel(), self.y.ravel(),
                                          self.u, self.v)
        _check_consistent_shapes(x, y, u, v)
        xy = np.hstack((x[:, np.newaxis], y[:, np.newaxis]))
        mcollections.PolyCollection.set_offsets(self, xy)
        self.stale = True

    set_offsets.__doc__ = mcollections.PolyCollection.set_offsets.__doc__

    barbs_doc = _barbs_doc
