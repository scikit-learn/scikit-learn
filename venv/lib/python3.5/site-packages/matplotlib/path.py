r"""
A module for dealing with the polylines used throughout Matplotlib.

The primary class for polyline handling in Matplotlib is `Path`.  Almost all
vector drawing makes use of `Path`\s somewhere in the drawing pipeline.

Whilst a `Path` instance itself cannot be drawn, some `.Artist` subclasses,
such as `.PathPatch` and `.PathCollection`, can be used for convenient `Path`
visualisation.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

from weakref import WeakValueDictionary

import numpy as np

from . import _path, rcParams
from .cbook import (_to_unmasked_float_array, simple_linear_interpolation,
                    maxdict)


class Path(object):
    """
    :class:`Path` represents a series of possibly disconnected,
    possibly closed, line and curve segments.

    The underlying storage is made up of two parallel numpy arrays:
      - *vertices*: an Nx2 float array of vertices
      - *codes*: an N-length uint8 array of vertex types

    These two arrays always have the same length in the first
    dimension.  For example, to represent a cubic curve, you must
    provide three vertices as well as three codes ``CURVE3``.

    The code types are:

       - ``STOP``   :  1 vertex (ignored)
           A marker for the end of the entire path (currently not
           required and ignored)

       - ``MOVETO`` :  1 vertex
            Pick up the pen and move to the given vertex.

       - ``LINETO`` :  1 vertex
            Draw a line from the current position to the given vertex.

       - ``CURVE3`` :  1 control point, 1 endpoint
          Draw a quadratic Bezier curve from the current position,
          with the given control point, to the given end point.

       - ``CURVE4`` :  2 control points, 1 endpoint
          Draw a cubic Bezier curve from the current position, with
          the given control points, to the given end point.

       - ``CLOSEPOLY`` : 1 vertex (ignored)
          Draw a line segment to the start point of the current
          polyline.

    Users of Path objects should not access the vertices and codes
    arrays directly.  Instead, they should use :meth:`iter_segments`
    or :meth:`cleaned` to get the vertex/code pairs.  This is important,
    since many :class:`Path` objects, as an optimization, do not store a
    *codes* at all, but have a default one provided for them by
    :meth:`iter_segments`.

    Some behavior of Path objects can be controlled by rcParams. See
    the rcParams whose keys contain 'path.'.

    .. note::

        The vertices and codes arrays should be treated as
        immutable -- there are a number of optimizations and assumptions
        made up front in the constructor that will not change when the
        data changes.

    """

    # Path codes
    STOP = 0         # 1 vertex
    MOVETO = 1       # 1 vertex
    LINETO = 2       # 1 vertex
    CURVE3 = 3       # 2 vertices
    CURVE4 = 4       # 3 vertices
    CLOSEPOLY = 79   # 1 vertex

    #: A dictionary mapping Path codes to the number of vertices that the
    #: code expects.
    NUM_VERTICES_FOR_CODE = {STOP: 1,
                             MOVETO: 1,
                             LINETO: 1,
                             CURVE3: 2,
                             CURVE4: 3,
                             CLOSEPOLY: 1}

    code_type = np.uint8

    def __init__(self, vertices, codes=None, _interpolation_steps=1,
                 closed=False, readonly=False):
        """
        Create a new path with the given vertices and codes.

        Parameters
        ----------
        vertices : array_like
            The ``(n, 2)`` float array, masked array or sequence of pairs
            representing the vertices of the path.

            If *vertices* contains masked values, they will be converted
            to NaNs which are then handled correctly by the Agg
            PathIterator and other consumers of path data, such as
            :meth:`iter_segments`.
        codes : {None, array_like}, optional
            n-length array integers representing the codes of the path.
            If not None, codes must be the same length as vertices.
            If None, *vertices* will be treated as a series of line segments.
        _interpolation_steps : int, optional
            Used as a hint to certain projections, such as Polar, that this
            path should be linearly interpolated immediately before drawing.
            This attribute is primarily an implementation detail and is not
            intended for public use.
        closed : bool, optional
            If *codes* is None and closed is True, vertices will be treated as
            line segments of a closed polygon.
        readonly : bool, optional
            Makes the path behave in an immutable way and sets the vertices
            and codes as read-only arrays.
        """
        vertices = _to_unmasked_float_array(vertices)
        if (vertices.ndim != 2) or (vertices.shape[1] != 2):
            raise ValueError(
                "'vertices' must be a 2D list or array with shape Nx2")

        if codes is not None:
            codes = np.asarray(codes, self.code_type)
            if (codes.ndim != 1) or len(codes) != len(vertices):
                raise ValueError("'codes' must be a 1D list or array with the "
                                 "same length of 'vertices'")
            if len(codes) and codes[0] != self.MOVETO:
                raise ValueError("The first element of 'code' must be equal "
                                 "to 'MOVETO' ({})".format(self.MOVETO))
        elif closed:
            codes = np.empty(len(vertices), dtype=self.code_type)
            codes[0] = self.MOVETO
            codes[1:-1] = self.LINETO
            codes[-1] = self.CLOSEPOLY

        self._vertices = vertices
        self._codes = codes
        self._interpolation_steps = _interpolation_steps
        self._update_values()

        if readonly:
            self._vertices.flags.writeable = False
            if self._codes is not None:
                self._codes.flags.writeable = False
            self._readonly = True
        else:
            self._readonly = False

    @classmethod
    def _fast_from_codes_and_verts(cls, verts, codes, internals=None):
        """
        Creates a Path instance without the expense of calling the constructor

        Parameters
        ----------
        verts : numpy array
        codes : numpy array
        internals : dict or None
            The attributes that the resulting path should have.
            Allowed keys are ``readonly``, ``should_simplify``,
            ``simplify_threshold``, ``has_nonfinite`` and
            ``interpolation_steps``.

        """
        internals = internals or {}
        pth = cls.__new__(cls)
        pth._vertices = _to_unmasked_float_array(verts)
        pth._codes = codes
        pth._readonly = internals.pop('readonly', False)
        pth.should_simplify = internals.pop('should_simplify', True)
        pth.simplify_threshold = (
            internals.pop('simplify_threshold',
                          rcParams['path.simplify_threshold'])
        )
        pth._has_nonfinite = internals.pop('has_nonfinite', False)
        pth._interpolation_steps = internals.pop('interpolation_steps', 1)
        if internals:
            raise ValueError('Unexpected internals provided to '
                             '_fast_from_codes_and_verts: '
                             '{0}'.format('\n *'.join(internals)))
        return pth

    def _update_values(self):
        self._simplify_threshold = rcParams['path.simplify_threshold']
        self._should_simplify = (
            self._simplify_threshold > 0 and
            rcParams['path.simplify'] and
            len(self._vertices) >= 128 and
            (self._codes is None or np.all(self._codes <= Path.LINETO))
        )
        self._has_nonfinite = not np.isfinite(self._vertices).all()

    @property
    def vertices(self):
        """
        The list of vertices in the `Path` as an Nx2 numpy array.
        """
        return self._vertices

    @vertices.setter
    def vertices(self, vertices):
        if self._readonly:
            raise AttributeError("Can't set vertices on a readonly Path")
        self._vertices = vertices
        self._update_values()

    @property
    def codes(self):
        """
        The list of codes in the `Path` as a 1-D numpy array.  Each
        code is one of `STOP`, `MOVETO`, `LINETO`, `CURVE3`, `CURVE4`
        or `CLOSEPOLY`.  For codes that correspond to more than one
        vertex (`CURVE3` and `CURVE4`), that code will be repeated so
        that the length of `self.vertices` and `self.codes` is always
        the same.
        """
        return self._codes

    @codes.setter
    def codes(self, codes):
        if self._readonly:
            raise AttributeError("Can't set codes on a readonly Path")
        self._codes = codes
        self._update_values()

    @property
    def simplify_threshold(self):
        """
        The fraction of a pixel difference below which vertices will
        be simplified out.
        """
        return self._simplify_threshold

    @simplify_threshold.setter
    def simplify_threshold(self, threshold):
        self._simplify_threshold = threshold

    @property
    def has_nonfinite(self):
        """
        `True` if the vertices array has nonfinite values.
        """
        return self._has_nonfinite

    @property
    def should_simplify(self):
        """
        `True` if the vertices array should be simplified.
        """
        return self._should_simplify

    @should_simplify.setter
    def should_simplify(self, should_simplify):
        self._should_simplify = should_simplify

    @property
    def readonly(self):
        """
        `True` if the `Path` is read-only.
        """
        return self._readonly

    def __copy__(self):
        """
        Returns a shallow copy of the `Path`, which will share the
        vertices and codes with the source `Path`.
        """
        import copy
        return copy.copy(self)

    copy = __copy__

    def __deepcopy__(self, memo=None):
        """
        Returns a deepcopy of the `Path`.  The `Path` will not be
        readonly, even if the source `Path` is.
        """
        try:
            codes = self.codes.copy()
        except AttributeError:
            codes = None
        return self.__class__(
            self.vertices.copy(), codes,
            _interpolation_steps=self._interpolation_steps)

    deepcopy = __deepcopy__

    @classmethod
    def make_compound_path_from_polys(cls, XY):
        """
        Make a compound path object to draw a number
        of polygons with equal numbers of sides XY is a (numpolys x
        numsides x 2) numpy array of vertices.  Return object is a
        :class:`Path`

        .. plot:: gallery/api/histogram_path.py

        """

        # for each poly: 1 for the MOVETO, (numsides-1) for the LINETO, 1 for
        # the CLOSEPOLY; the vert for the closepoly is ignored but we still
        # need it to keep the codes aligned with the vertices
        numpolys, numsides, two = XY.shape
        if two != 2:
            raise ValueError("The third dimension of 'XY' must be 2")
        stride = numsides + 1
        nverts = numpolys * stride
        verts = np.zeros((nverts, 2))
        codes = np.ones(nverts, int) * cls.LINETO
        codes[0::stride] = cls.MOVETO
        codes[numsides::stride] = cls.CLOSEPOLY
        for i in range(numsides):
            verts[i::stride] = XY[:, i]

        return cls(verts, codes)

    @classmethod
    def make_compound_path(cls, *args):
        """Make a compound path from a list of Path objects."""
        # Handle an empty list in args (i.e. no args).
        if not args:
            return Path(np.empty([0, 2], dtype=np.float32))

        lengths = [len(x) for x in args]
        total_length = sum(lengths)

        vertices = np.vstack([x.vertices for x in args])
        vertices.reshape((total_length, 2))

        codes = np.empty(total_length, dtype=cls.code_type)
        i = 0
        for path in args:
            if path.codes is None:
                codes[i] = cls.MOVETO
                codes[i + 1:i + len(path.vertices)] = cls.LINETO
            else:
                codes[i:i + len(path.codes)] = path.codes
            i += len(path.vertices)

        return cls(vertices, codes)

    def __repr__(self):
        return "Path(%r, %r)" % (self.vertices, self.codes)

    def __len__(self):
        return len(self.vertices)

    def iter_segments(self, transform=None, remove_nans=True, clip=None,
                      snap=False, stroke_width=1.0, simplify=None,
                      curves=True, sketch=None):
        """
        Iterates over all of the curve segments in the path.  Each
        iteration returns a 2-tuple (*vertices*, *code*), where
        *vertices* is a sequence of 1 - 3 coordinate pairs, and *code* is
        one of the :class:`Path` codes.

        Additionally, this method can provide a number of standard
        cleanups and conversions to the path.

        Parameters
        ----------
        transform : None or :class:`~matplotlib.transforms.Transform` instance
            If not None, the given affine transformation will
            be applied to the path.
        remove_nans : {False, True}, optional
            If True, will remove all NaNs from the path and
            insert MOVETO commands to skip over them.
        clip : None or sequence, optional
            If not None, must be a four-tuple (x1, y1, x2, y2)
            defining a rectangle in which to clip the path.
        snap : None or bool, optional
            If None, auto-snap to pixels, to reduce
            fuzziness of rectilinear lines.  If True, force snapping, and
            if False, don't snap.
        stroke_width : float, optional
            The width of the stroke being drawn.  Needed
             as a hint for the snapping algorithm.
        simplify : None or bool, optional
            If True, perform simplification, to remove
             vertices that do not affect the appearance of the path.  If
             False, perform no simplification.  If None, use the
             should_simplify member variable.  See also the rcParams
             path.simplify and path.simplify_threshold.
        curves : {True, False}, optional
            If True, curve segments will be returned as curve
            segments.  If False, all curves will be converted to line
            segments.
        sketch : None or sequence, optional
            If not None, must be a 3-tuple of the form
            (scale, length, randomness), representing the sketch
            parameters.
        """
        if not len(self):
            return

        cleaned = self.cleaned(transform=transform,
                               remove_nans=remove_nans, clip=clip,
                               snap=snap, stroke_width=stroke_width,
                               simplify=simplify, curves=curves,
                               sketch=sketch)
        vertices = cleaned.vertices
        codes = cleaned.codes
        len_vertices = vertices.shape[0]

        # Cache these object lookups for performance in the loop.
        NUM_VERTICES_FOR_CODE = self.NUM_VERTICES_FOR_CODE
        STOP = self.STOP

        i = 0
        while i < len_vertices:
            code = codes[i]
            if code == STOP:
                return
            else:
                num_vertices = NUM_VERTICES_FOR_CODE[code]
                curr_vertices = vertices[i:i+num_vertices].flatten()
                yield curr_vertices, code
                i += num_vertices

    def cleaned(self, transform=None, remove_nans=False, clip=None,
                quantize=False, simplify=False, curves=False,
                stroke_width=1.0, snap=False, sketch=None):
        """
        Cleans up the path according to the parameters returning a new
        Path instance.

        .. seealso::

            See :meth:`iter_segments` for details of the keyword arguments.

        Returns
        -------
        Path instance with cleaned up vertices and codes.

        """
        vertices, codes = _path.cleanup_path(self, transform,
                                             remove_nans, clip,
                                             snap, stroke_width,
                                             simplify, curves, sketch)
        internals = {'should_simplify': self.should_simplify and not simplify,
                     'has_nonfinite': self.has_nonfinite and not remove_nans,
                     'simplify_threshold': self.simplify_threshold,
                     'interpolation_steps': self._interpolation_steps}
        return Path._fast_from_codes_and_verts(vertices, codes, internals)

    def transformed(self, transform):
        """
        Return a transformed copy of the path.

        .. seealso::

            :class:`matplotlib.transforms.TransformedPath`
                A specialized path class that will cache the
                transformed result and automatically update when the
                transform changes.
        """
        return Path(transform.transform(self.vertices), self.codes,
                    self._interpolation_steps)

    def contains_point(self, point, transform=None, radius=0.0):
        """
        Returns whether the (closed) path contains the given point.

        If *transform* is not ``None``, the path will be transformed before
        performing the test.

        *radius* allows the path to be made slightly larger or smaller.
        """
        if transform is not None:
            transform = transform.frozen()
        # `point_in_path` does not handle nonlinear transforms, so we
        # transform the path ourselves.  If `transform` is affine, letting
        # `point_in_path` handle the transform avoids allocating an extra
        # buffer.
        if transform and not transform.is_affine:
            self = transform.transform_path(self)
            transform = None
        return _path.point_in_path(point[0], point[1], radius, self, transform)

    def contains_points(self, points, transform=None, radius=0.0):
        """
        Returns a bool array which is ``True`` if the (closed) path contains
        the corresponding point.

        If *transform* is not ``None``, the path will be transformed before
        performing the test.

        *radius* allows the path to be made slightly larger or smaller.
        """
        if transform is not None:
            transform = transform.frozen()
        result = _path.points_in_path(points, radius, self, transform)
        return result.astype('bool')

    def contains_path(self, path, transform=None):
        """
        Returns whether this (closed) path completely contains the given path.

        If *transform* is not ``None``, the path will be transformed before
        performing the test.
        """
        if transform is not None:
            transform = transform.frozen()
        return _path.path_in_path(self, None, path, transform)

    def get_extents(self, transform=None):
        """
        Returns the extents (*xmin*, *ymin*, *xmax*, *ymax*) of the
        path.

        Unlike computing the extents on the *vertices* alone, this
        algorithm will take into account the curves and deal with
        control points appropriately.
        """
        from .transforms import Bbox
        path = self
        if transform is not None:
            transform = transform.frozen()
            if not transform.is_affine:
                path = self.transformed(transform)
                transform = None
        return Bbox(_path.get_path_extents(path, transform))

    def intersects_path(self, other, filled=True):
        """
        Returns *True* if this path intersects another given path.

        *filled*, when True, treats the paths as if they were filled.
        That is, if one path completely encloses the other,
        :meth:`intersects_path` will return True.
        """
        return _path.path_intersects_path(self, other, filled)

    def intersects_bbox(self, bbox, filled=True):
        """
        Returns *True* if this path intersects a given
        :class:`~matplotlib.transforms.Bbox`.

        *filled*, when True, treats the path as if it was filled.
        That is, if the path completely encloses the bounding box,
        :meth:`intersects_bbox` will return True.

        The bounding box is always considered filled.
        """
        return _path.path_intersects_rectangle(self,
            bbox.x0, bbox.y0, bbox.x1, bbox.y1, filled)

    def interpolated(self, steps):
        """
        Returns a new path resampled to length N x steps.  Does not
        currently handle interpolating curves.
        """
        if steps == 1:
            return self

        vertices = simple_linear_interpolation(self.vertices, steps)
        codes = self.codes
        if codes is not None:
            new_codes = Path.LINETO * np.ones(((len(codes) - 1) * steps + 1, ))
            new_codes[0::steps] = codes
        else:
            new_codes = None
        return Path(vertices, new_codes)

    def to_polygons(self, transform=None, width=0, height=0, closed_only=True):
        """
        Convert this path to a list of polygons or polylines.  Each
        polygon/polyline is an Nx2 array of vertices.  In other words,
        each polygon has no ``MOVETO`` instructions or curves.  This
        is useful for displaying in backends that do not support
        compound paths or Bezier curves, such as GDK.

        If *width* and *height* are both non-zero then the lines will
        be simplified so that vertices outside of (0, 0), (width,
        height) will be clipped.

        If *closed_only* is `True` (default), only closed polygons,
        with the last point being the same as the first point, will be
        returned.  Any unclosed polylines in the path will be
        explicitly closed.  If *closed_only* is `False`, any unclosed
        polygons in the path will be returned as unclosed polygons,
        and the closed polygons will be returned explicitly closed by
        setting the last point to the same as the first point.
        """
        if len(self.vertices) == 0:
            return []

        if transform is not None:
            transform = transform.frozen()

        if self.codes is None and (width == 0 or height == 0):
            vertices = self.vertices
            if closed_only:
                if len(vertices) < 3:
                    return []
                elif np.any(vertices[0] != vertices[-1]):
                    vertices = list(vertices) + [vertices[0]]

            if transform is None:
                return [vertices]
            else:
                return [transform.transform(vertices)]

        # Deal with the case where there are curves and/or multiple
        # subpaths (using extension code)
        return _path.convert_path_to_polygons(
            self, transform, width, height, closed_only)

    _unit_rectangle = None

    @classmethod
    def unit_rectangle(cls):
        """
        Return a :class:`Path` instance of the unit rectangle
        from (0, 0) to (1, 1).
        """
        if cls._unit_rectangle is None:
            cls._unit_rectangle = \
                cls([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0],
                     [0.0, 0.0]],
                    [cls.MOVETO, cls.LINETO, cls.LINETO, cls.LINETO,
                     cls.CLOSEPOLY],
                    readonly=True)
        return cls._unit_rectangle

    _unit_regular_polygons = WeakValueDictionary()

    @classmethod
    def unit_regular_polygon(cls, numVertices):
        """
        Return a :class:`Path` instance for a unit regular
        polygon with the given *numVertices* and radius of 1.0,
        centered at (0, 0).
        """
        if numVertices <= 16:
            path = cls._unit_regular_polygons.get(numVertices)
        else:
            path = None
        if path is None:
            theta = (2*np.pi/numVertices *
                     np.arange(numVertices + 1).reshape((numVertices + 1, 1)))
            # This initial rotation is to make sure the polygon always
            # "points-up"
            theta += np.pi / 2.0
            verts = np.concatenate((np.cos(theta), np.sin(theta)), 1)
            codes = np.empty((numVertices + 1,))
            codes[0] = cls.MOVETO
            codes[1:-1] = cls.LINETO
            codes[-1] = cls.CLOSEPOLY
            path = cls(verts, codes, readonly=True)
            if numVertices <= 16:
                cls._unit_regular_polygons[numVertices] = path
        return path

    _unit_regular_stars = WeakValueDictionary()

    @classmethod
    def unit_regular_star(cls, numVertices, innerCircle=0.5):
        """
        Return a :class:`Path` for a unit regular star
        with the given numVertices and radius of 1.0, centered at (0,
        0).
        """
        if numVertices <= 16:
            path = cls._unit_regular_stars.get((numVertices, innerCircle))
        else:
            path = None
        if path is None:
            ns2 = numVertices * 2
            theta = (2*np.pi/ns2 * np.arange(ns2 + 1))
            # This initial rotation is to make sure the polygon always
            # "points-up"
            theta += np.pi / 2.0
            r = np.ones(ns2 + 1)
            r[1::2] = innerCircle
            verts = np.vstack((r*np.cos(theta), r*np.sin(theta))).transpose()
            codes = np.empty((ns2 + 1,))
            codes[0] = cls.MOVETO
            codes[1:-1] = cls.LINETO
            codes[-1] = cls.CLOSEPOLY
            path = cls(verts, codes, readonly=True)
            if numVertices <= 16:
                cls._unit_regular_stars[(numVertices, innerCircle)] = path
        return path

    @classmethod
    def unit_regular_asterisk(cls, numVertices):
        """
        Return a :class:`Path` for a unit regular
        asterisk with the given numVertices and radius of 1.0,
        centered at (0, 0).
        """
        return cls.unit_regular_star(numVertices, 0.0)

    _unit_circle = None

    @classmethod
    def unit_circle(cls):
        """
        Return the readonly :class:`Path` of the unit circle.

        For most cases, :func:`Path.circle` will be what you want.

        """
        if cls._unit_circle is None:
            cls._unit_circle = cls.circle(center=(0, 0), radius=1,
                                          readonly=True)
        return cls._unit_circle

    @classmethod
    def circle(cls, center=(0., 0.), radius=1., readonly=False):
        """
        Return a Path representing a circle of a given radius and center.

        Parameters
        ----------
        center : pair of floats
            The center of the circle. Default ``(0, 0)``.
        radius : float
            The radius of the circle. Default is 1.
        readonly : bool
            Whether the created path should have the "readonly" argument
            set when creating the Path instance.

        Notes
        -----
        The circle is approximated using cubic Bezier curves.  This
        uses 8 splines around the circle using the approach presented
        here:

          Lancaster, Don.  `Approximating a Circle or an Ellipse Using Four
          Bezier Cubic Splines <http://www.tinaja.com/glib/ellipse4.pdf>`_.

        """
        MAGIC = 0.2652031
        SQRTHALF = np.sqrt(0.5)
        MAGIC45 = SQRTHALF * MAGIC

        vertices = np.array([[0.0, -1.0],

                             [MAGIC, -1.0],
                             [SQRTHALF-MAGIC45, -SQRTHALF-MAGIC45],
                             [SQRTHALF, -SQRTHALF],

                             [SQRTHALF+MAGIC45, -SQRTHALF+MAGIC45],
                             [1.0, -MAGIC],
                             [1.0, 0.0],

                             [1.0, MAGIC],
                             [SQRTHALF+MAGIC45, SQRTHALF-MAGIC45],
                             [SQRTHALF, SQRTHALF],

                             [SQRTHALF-MAGIC45, SQRTHALF+MAGIC45],
                             [MAGIC, 1.0],
                             [0.0, 1.0],

                             [-MAGIC, 1.0],
                             [-SQRTHALF+MAGIC45, SQRTHALF+MAGIC45],
                             [-SQRTHALF, SQRTHALF],

                             [-SQRTHALF-MAGIC45, SQRTHALF-MAGIC45],
                             [-1.0, MAGIC],
                             [-1.0, 0.0],

                             [-1.0, -MAGIC],
                             [-SQRTHALF-MAGIC45, -SQRTHALF+MAGIC45],
                             [-SQRTHALF, -SQRTHALF],

                             [-SQRTHALF+MAGIC45, -SQRTHALF-MAGIC45],
                             [-MAGIC, -1.0],
                             [0.0, -1.0],

                             [0.0, -1.0]],
                            dtype=float)

        codes = [cls.CURVE4] * 26
        codes[0] = cls.MOVETO
        codes[-1] = cls.CLOSEPOLY
        return Path(vertices * radius + center, codes, readonly=readonly)

    _unit_circle_righthalf = None

    @classmethod
    def unit_circle_righthalf(cls):
        """
        Return a :class:`Path` of the right half
        of a unit circle. The circle is approximated using cubic Bezier
        curves.  This uses 4 splines around the circle using the approach
        presented here:

          Lancaster, Don.  `Approximating a Circle or an Ellipse Using Four
          Bezier Cubic Splines <http://www.tinaja.com/glib/ellipse4.pdf>`_.
        """
        if cls._unit_circle_righthalf is None:
            MAGIC = 0.2652031
            SQRTHALF = np.sqrt(0.5)
            MAGIC45 = SQRTHALF * MAGIC

            vertices = np.array(
                [[0.0, -1.0],

                 [MAGIC, -1.0],
                 [SQRTHALF-MAGIC45, -SQRTHALF-MAGIC45],
                 [SQRTHALF, -SQRTHALF],

                 [SQRTHALF+MAGIC45, -SQRTHALF+MAGIC45],
                 [1.0, -MAGIC],
                 [1.0, 0.0],

                 [1.0, MAGIC],
                 [SQRTHALF+MAGIC45, SQRTHALF-MAGIC45],
                 [SQRTHALF, SQRTHALF],

                 [SQRTHALF-MAGIC45, SQRTHALF+MAGIC45],
                 [MAGIC, 1.0],
                 [0.0, 1.0],

                 [0.0, -1.0]],

                float)

            codes = cls.CURVE4 * np.ones(14)
            codes[0] = cls.MOVETO
            codes[-1] = cls.CLOSEPOLY

            cls._unit_circle_righthalf = cls(vertices, codes, readonly=True)
        return cls._unit_circle_righthalf

    @classmethod
    def arc(cls, theta1, theta2, n=None, is_wedge=False):
        """
        Return an arc on the unit circle from angle
        *theta1* to angle *theta2* (in degrees).

        *theta2* is unwrapped to produce the shortest arc within 360 degrees.
        That is, if *theta2* > *theta1* + 360, the arc will be from *theta1* to
        *theta2* - 360 and not a full circle plus some extra overlap.

        If *n* is provided, it is the number of spline segments to make.
        If *n* is not provided, the number of spline segments is
        determined based on the delta between *theta1* and *theta2*.

           Masionobe, L.  2003.  `Drawing an elliptical arc using
           polylines, quadratic or cubic Bezier curves
           <http://www.spaceroots.org/documents/ellipse/index.html>`_.
        """
        halfpi = np.pi * 0.5

        eta1 = theta1
        eta2 = theta2 - 360 * np.floor((theta2 - theta1) / 360)
        # Ensure 2pi range is not flattened to 0 due to floating-point errors,
        # but don't try to expand existing 0 range.
        if theta2 != theta1 and eta2 <= eta1:
            eta2 += 360
        eta1, eta2 = np.deg2rad([eta1, eta2])

        # number of curve segments to make
        if n is None:
            n = int(2 ** np.ceil((eta2 - eta1) / halfpi))
        if n < 1:
            raise ValueError("n must be >= 1 or None")

        deta = (eta2 - eta1) / n
        t = np.tan(0.5 * deta)
        alpha = np.sin(deta) * (np.sqrt(4.0 + 3.0 * t * t) - 1) / 3.0

        steps = np.linspace(eta1, eta2, n + 1, True)
        cos_eta = np.cos(steps)
        sin_eta = np.sin(steps)

        xA = cos_eta[:-1]
        yA = sin_eta[:-1]
        xA_dot = -yA
        yA_dot = xA

        xB = cos_eta[1:]
        yB = sin_eta[1:]
        xB_dot = -yB
        yB_dot = xB

        if is_wedge:
            length = n * 3 + 4
            vertices = np.zeros((length, 2), float)
            codes = cls.CURVE4 * np.ones((length, ), cls.code_type)
            vertices[1] = [xA[0], yA[0]]
            codes[0:2] = [cls.MOVETO, cls.LINETO]
            codes[-2:] = [cls.LINETO, cls.CLOSEPOLY]
            vertex_offset = 2
            end = length - 2
        else:
            length = n * 3 + 1
            vertices = np.empty((length, 2), float)
            codes = cls.CURVE4 * np.ones((length, ), cls.code_type)
            vertices[0] = [xA[0], yA[0]]
            codes[0] = cls.MOVETO
            vertex_offset = 1
            end = length

        vertices[vertex_offset:end:3, 0] = xA + alpha * xA_dot
        vertices[vertex_offset:end:3, 1] = yA + alpha * yA_dot
        vertices[vertex_offset+1:end:3, 0] = xB - alpha * xB_dot
        vertices[vertex_offset+1:end:3, 1] = yB - alpha * yB_dot
        vertices[vertex_offset+2:end:3, 0] = xB
        vertices[vertex_offset+2:end:3, 1] = yB

        return cls(vertices, codes, readonly=True)

    @classmethod
    def wedge(cls, theta1, theta2, n=None):
        """
        Return a wedge of the unit circle from angle
        *theta1* to angle *theta2* (in degrees).

        *theta2* is unwrapped to produce the shortest wedge within 360 degrees.
        That is, if *theta2* > *theta1* + 360, the wedge will be from *theta1*
        to *theta2* - 360 and not a full circle plus some extra overlap.

        If *n* is provided, it is the number of spline segments to make.
        If *n* is not provided, the number of spline segments is
        determined based on the delta between *theta1* and *theta2*.
        """
        return cls.arc(theta1, theta2, n, True)

    _hatch_dict = maxdict(8)

    @classmethod
    def hatch(cls, hatchpattern, density=6):
        """
        Given a hatch specifier, *hatchpattern*, generates a Path that
        can be used in a repeated hatching pattern.  *density* is the
        number of lines per unit square.
        """
        from matplotlib.hatch import get_path

        if hatchpattern is None:
            return None

        hatch_path = cls._hatch_dict.get((hatchpattern, density))
        if hatch_path is not None:
            return hatch_path

        hatch_path = get_path(hatchpattern, density)
        cls._hatch_dict[(hatchpattern, density)] = hatch_path
        return hatch_path

    def clip_to_bbox(self, bbox, inside=True):
        """
        Clip the path to the given bounding box.

        The path must be made up of one or more closed polygons.  This
        algorithm will not behave correctly for unclosed paths.

        If *inside* is `True`, clip to the inside of the box, otherwise
        to the outside of the box.
        """
        # Use make_compound_path_from_polys
        verts = _path.clip_path_to_rect(self, bbox, inside)
        paths = [Path(poly) for poly in verts]
        return self.make_compound_path(*paths)


def get_path_collection_extents(
        master_transform, paths, transforms, offsets, offset_transform):
    """
    Given a sequence of :class:`Path` objects,
    :class:`~matplotlib.transforms.Transform` objects and offsets, as
    found in a :class:`~matplotlib.collections.PathCollection`,
    returns the bounding box that encapsulates all of them.

    *master_transform* is a global transformation to apply to all paths

    *paths* is a sequence of :class:`Path` instances.

    *transforms* is a sequence of
    :class:`~matplotlib.transforms.Affine2D` instances.

    *offsets* is a sequence of (x, y) offsets (or an Nx2 array)

    *offset_transform* is a :class:`~matplotlib.transforms.Affine2D`
    to apply to the offsets before applying the offset to the path.

    The way that *paths*, *transforms* and *offsets* are combined
    follows the same method as for collections.  Each is iterated over
    independently, so if you have 3 paths, 2 transforms and 1 offset,
    their combinations are as follows:

        (A, A, A), (B, B, A), (C, A, A)
    """
    from .transforms import Bbox
    if len(paths) == 0:
        raise ValueError("No paths provided")
    return Bbox.from_extents(*_path.get_path_collection_extents(
        master_transform, paths, np.atleast_3d(transforms),
        offsets, offset_transform))


def get_paths_extents(paths, transforms=[]):
    """
    Given a sequence of :class:`Path` objects and optional
    :class:`~matplotlib.transforms.Transform` objects, returns the
    bounding box that encapsulates all of them.

    *paths* is a sequence of :class:`Path` instances.

    *transforms* is an optional sequence of
    :class:`~matplotlib.transforms.Affine2D` instances to apply to
    each path.
    """
    from .transforms import Bbox, Affine2D
    if len(paths) == 0:
        raise ValueError("No paths provided")
    return Bbox.from_extents(*_path.get_path_collection_extents(
        Affine2D(), paths, transforms, [], Affine2D()))
