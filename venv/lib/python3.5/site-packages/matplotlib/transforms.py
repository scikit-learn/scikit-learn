"""
matplotlib includes a framework for arbitrary geometric
transformations that is used determine the final position of all
elements drawn on the canvas.

Transforms are composed into trees of :class:`TransformNode` objects
whose actual value depends on their children.  When the contents of
children change, their parents are automatically invalidated.  The
next time an invalidated transform is accessed, it is recomputed to
reflect those changes.  This invalidation/caching approach prevents
unnecessary recomputations of transforms, and contributes to better
interactive performance.

For example, here is a graph of the transform tree used to plot data
to the graph:

.. image:: ../_static/transforms.png

The framework can be used for both affine and non-affine
transformations.  However, for speed, we want use the backend
renderers to perform affine transformations whenever possible.
Therefore, it is possible to perform just the affine or non-affine
part of a transformation on a set of data.  The affine is always
assumed to occur after the non-affine.  For any transform::

  full transform == non-affine part + affine part

The backends are not expected to handle non-affine transformations
themselves.
"""

# Note: There are a number of places in the code where we use `np.min` or
# `np.minimum` instead of the builtin `min`, and likewise for `max`.  This is
# done so that `nan`s are propagated, instead of being silently dropped.

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import numpy as np
from matplotlib._path import (affine_transform, count_bboxes_overlapping_bbox,
    update_path_extents)
from numpy.linalg import inv

import re
import weakref
import warnings

from . import cbook
from .path import Path

DEBUG = False


def _indent_str(obj):  # textwrap.indent(str(obj), 4) on Py3.
    return re.sub("(^|\n)", r"\1    ", str(obj))


class TransformNode(object):
    """
    :class:`TransformNode` is the base class for anything that
    participates in the transform tree and needs to invalidate its
    parents or be invalidated.  This includes classes that are not
    really transforms, such as bounding boxes, since some transforms
    depend on bounding boxes to compute their values.
    """
    _gid = 0

    # Invalidation may affect only the affine part.  If the
    # invalidation was "affine-only", the _invalid member is set to
    # INVALID_AFFINE_ONLY
    INVALID_NON_AFFINE = 1
    INVALID_AFFINE = 2
    INVALID = INVALID_NON_AFFINE | INVALID_AFFINE

    # Some metadata about the transform, used to determine whether an
    # invalidation is affine-only
    is_affine = False
    is_bbox = False

    pass_through = False
    """
    If pass_through is True, all ancestors will always be
    invalidated, even if 'self' is already invalid.
    """

    def __init__(self, shorthand_name=None):
        """
        Creates a new :class:`TransformNode`.

        Parameters
        ----------
        shorthand_name : str
            A string representing the "name" of the transform. The name carries
            no significance other than to improve the readability of
            ``str(transform)`` when DEBUG=True.
        """
        self._parents = {}

        # TransformNodes start out as invalid until their values are
        # computed for the first time.
        self._invalid = 1
        self._shorthand_name = shorthand_name or ''

    if DEBUG:
        def __str__(self):
            # either just return the name of this TransformNode, or it's repr
            return self._shorthand_name or repr(self)

    def __getstate__(self):
        d = self.__dict__.copy()
        # turn the dictionary with weak values into a normal dictionary
        d['_parents'] = dict((k, v()) for (k, v) in
                             six.iteritems(self._parents))
        return d

    def __setstate__(self, data_dict):
        self.__dict__ = data_dict
        # turn the normal dictionary back into a dictionary with weak
        # values
        self._parents = dict((k, weakref.ref(v)) for (k, v) in
                             six.iteritems(self._parents) if v is not None)

    def __copy__(self, *args):
        raise NotImplementedError(
            "TransformNode instances can not be copied. "
            "Consider using frozen() instead.")
    __deepcopy__ = __copy__

    def invalidate(self):
        """
        Invalidate this :class:`TransformNode` and triggers an
        invalidation of its ancestors.  Should be called any
        time the transform changes.
        """
        value = self.INVALID
        if self.is_affine:
            value = self.INVALID_AFFINE
        return self._invalidate_internal(value, invalidating_node=self)

    def _invalidate_internal(self, value, invalidating_node):
        """
        Called by :meth:`invalidate` and subsequently ascends the transform
        stack calling each TransformNode's _invalidate_internal method.
        """
        # determine if this call will be an extension to the invalidation
        # status. If not, then a shortcut means that we needn't invoke an
        # invalidation up the transform stack as it will already have been
        # invalidated.

        # N.B This makes the invalidation sticky, once a transform has been
        # invalidated as NON_AFFINE, then it will always be invalidated as
        # NON_AFFINE even when triggered with a AFFINE_ONLY invalidation.
        # In most cases this is not a problem (i.e. for interactive panning and
        # zooming) and the only side effect will be on performance.
        status_changed = self._invalid < value

        if self.pass_through or status_changed:
            self._invalid = value

            for parent in list(six.itervalues(self._parents)):
                # Dereference the weak reference
                parent = parent()
                if parent is not None:
                    parent._invalidate_internal(
                        value=value, invalidating_node=self)

    def set_children(self, *children):
        """
        Set the children of the transform, to let the invalidation
        system know which transforms can invalidate this transform.
        Should be called from the constructor of any transforms that
        depend on other transforms.
        """
        # Parents are stored as weak references, so that if the
        # parents are destroyed, references from the children won't
        # keep them alive.
        for child in children:
            child._parents[id(self)] = weakref.ref(self)

    if DEBUG:
        _set_children = set_children

        def set_children(self, *children):
            self._set_children(*children)
            self._children = children
        set_children.__doc__ = _set_children.__doc__

    def frozen(self):
        """
        Returns a frozen copy of this transform node.  The frozen copy
        will not update when its children change.  Useful for storing
        a previously known state of a transform where
        ``copy.deepcopy()`` might normally be used.
        """
        return self

    if DEBUG:
        def write_graphviz(self, fobj, highlight=[]):
            """
            For debugging purposes.

            Writes the transform tree rooted at 'self' to a graphviz "dot"
            format file.  This file can be run through the "dot" utility
            to produce a graph of the transform tree.

            Affine transforms are marked in blue.  Bounding boxes are
            marked in yellow.

            *fobj*: A Python file-like object

            Once the "dot" file has been created, it can be turned into a
            png easily with::

                $> dot -Tpng -o $OUTPUT_FILE $DOT_FILE

            """
            seen = set()

            def recurse(root):
                if root in seen:
                    return
                seen.add(root)
                props = {}
                label = root.__class__.__name__
                if root._invalid:
                    label = '[%s]' % label
                if root in highlight:
                    props['style'] = 'bold'
                props['shape'] = 'box'
                props['label'] = '"%s"' % label
                props = ' '.join(['%s=%s' % (key, val)
                                  for key, val
                                  in six.iteritems(props)])

                fobj.write('%s [%s];\n' %
                           (hash(root), props))

                if hasattr(root, '_children'):
                    for child in root._children:
                        name = '?'
                        for key, val in six.iteritems(root.__dict__):
                            if val is child:
                                name = key
                                break
                        fobj.write('"%s" -> "%s" [label="%s", fontsize=10];\n'
                                    % (hash(root),
                                    hash(child),
                                    name))
                        recurse(child)

            fobj.write("digraph G {\n")
            recurse(self)
            fobj.write("}\n")


class BboxBase(TransformNode):
    """
    This is the base class of all bounding boxes, and provides
    read-only access to its data.  A mutable bounding box is provided
    by the :class:`Bbox` class.

    The canonical representation is as two points, with no
    restrictions on their ordering.  Convenience properties are
    provided to get the left, bottom, right and top edges and width
    and height, but these are not stored explicitly.
    """
    is_bbox = True
    is_affine = True

    if DEBUG:
        def _check(points):
            if isinstance(points, np.ma.MaskedArray):
                warnings.warn("Bbox bounds are a masked array.")
            points = np.asarray(points)
            if (points[1, 0] - points[0, 0] == 0 or
                points[1, 1] - points[0, 1] == 0):
                warnings.warn("Singular Bbox.")
        _check = staticmethod(_check)

    def frozen(self):
        return Bbox(self.get_points().copy())
    frozen.__doc__ = TransformNode.__doc__

    def __array__(self, *args, **kwargs):
        return self.get_points()

    def is_unit(self):
        """
        Returns True if the :class:`Bbox` is the unit bounding box
        from (0, 0) to (1, 1).
        """
        return list(self.get_points().flatten()) == [0., 0., 1., 1.]

    @property
    def x0(self):
        """
        :attr:`x0` is the first of the pair of *x* coordinates that
        define the bounding box. :attr:`x0` is not guaranteed to be less than
        :attr:`x1`.  If you require that, use :attr:`xmin`.
        """
        return self.get_points()[0, 0]

    @property
    def y0(self):
        """
        :attr:`y0` is the first of the pair of *y* coordinates that
        define the bounding box. :attr:`y0` is not guaranteed to be less than
        :attr:`y1`.  If you require that, use :attr:`ymin`.
        """
        return self.get_points()[0, 1]

    @property
    def x1(self):
        """
        :attr:`x1` is the second of the pair of *x* coordinates that
        define the bounding box. :attr:`x1` is not guaranteed to be greater
        than :attr:`x0`.  If you require that, use :attr:`xmax`.
        """
        return self.get_points()[1, 0]

    @property
    def y1(self):
        """
        :attr:`y1` is the second of the pair of *y* coordinates that
        define the bounding box. :attr:`y1` is not guaranteed to be greater
        than :attr:`y0`.  If you require that, use :attr:`ymax`.
        """
        return self.get_points()[1, 1]

    @property
    def p0(self):
        """
        :attr:`p0` is the first pair of (*x*, *y*) coordinates that
        define the bounding box.  It is not guaranteed to be the bottom-left
        corner.  For that, use :attr:`min`.
        """
        return self.get_points()[0]

    @property
    def p1(self):
        """
        :attr:`p1` is the second pair of (*x*, *y*) coordinates that
        define the bounding box.  It is not guaranteed to be the top-right
        corner.  For that, use :attr:`max`.
        """
        return self.get_points()[1]

    @property
    def xmin(self):
        """
        :attr:`xmin` is the left edge of the bounding box.
        """
        return np.min(self.get_points()[:, 0])

    @property
    def ymin(self):
        """
        :attr:`ymin` is the bottom edge of the bounding box.
        """
        return np.min(self.get_points()[:, 1])

    @property
    def xmax(self):
        """
        :attr:`xmax` is the right edge of the bounding box.
        """
        return np.max(self.get_points()[:, 0])

    @property
    def ymax(self):
        """
        :attr:`ymax` is the top edge of the bounding box.
        """
        return np.max(self.get_points()[:, 1])

    @property
    def min(self):
        """
        :attr:`min` is the bottom-left corner of the bounding box.
        """
        return np.min(self.get_points(), axis=0)

    @property
    def max(self):
        """
        :attr:`max` is the top-right corner of the bounding box.
        """
        return np.max(self.get_points(), axis=0)

    @property
    def intervalx(self):
        """
        :attr:`intervalx` is the pair of *x* coordinates that define
        the bounding box. It is not guaranteed to be sorted from left to right.
        """
        return self.get_points()[:, 0]

    @property
    def intervaly(self):
        """
        :attr:`intervaly` is the pair of *y* coordinates that define
        the bounding box.  It is not guaranteed to be sorted from bottom to
        top.
        """
        return self.get_points()[:, 1]

    @property
    def width(self):
        """
        The width of the bounding box.  It may be negative if
        :attr:`x1` < :attr:`x0`.
        """
        points = self.get_points()
        return points[1, 0] - points[0, 0]

    @property
    def height(self):
        """
        The height of the bounding box.  It may be negative if
        :attr:`y1` < :attr:`y0`.
        """
        points = self.get_points()
        return points[1, 1] - points[0, 1]

    @property
    def size(self):
        """
        The width and height of the bounding box.  May be negative,
        in the same way as :attr:`width` and :attr:`height`.
        """
        points = self.get_points()
        return points[1] - points[0]

    @property
    def bounds(self):
        """
        Returns (:attr:`x0`, :attr:`y0`, :attr:`width`,
        :attr:`height`).
        """
        x0, y0, x1, y1 = self.get_points().flatten()
        return (x0, y0, x1 - x0, y1 - y0)

    @property
    def extents(self):
        """
        Returns (:attr:`x0`, :attr:`y0`, :attr:`x1`,
        :attr:`y1`).
        """
        return self.get_points().flatten().copy()

    def get_points(self):
        raise NotImplementedError

    def containsx(self, x):
        """
        Returns whether *x* is in the closed (:attr:`x0`, :attr:`x1`) interval.
        """
        x0, x1 = self.intervalx
        return x0 <= x <= x1 or x0 >= x >= x1

    def containsy(self, y):
        """
        Returns whether *y* is in the closed (:attr:`y0`, :attr:`y1`) interval.
        """
        y0, y1 = self.intervaly
        return y0 <= y <= y1 or y0 >= y >= y1

    def contains(self, x, y):
        """
        Returns whether ``(x, y)`` is in the bounding box or on its edge.
        """
        return self.containsx(x) and self.containsy(y)

    def overlaps(self, other):
        """
        Returns whether this bounding box overlaps with the other bounding box.

        Parameters
        ----------
        other : BboxBase
        """
        ax1, ay1, ax2, ay2 = self.extents
        bx1, by1, bx2, by2 = other.extents
        if ax2 < ax1:
            ax2, ax1 = ax1, ax2
        if ay2 < ay1:
            ay2, ay1 = ay1, ay2
        if bx2 < bx1:
            bx2, bx1 = bx1, bx2
        if by2 < by1:
            by2, by1 = by1, by2
        return ax1 <= bx2 and bx1 <= ax2 and ay1 <= by2 and by1 <= ay2

    def fully_containsx(self, x):
        """
        Returns whether *x* is in the open (:attr:`x0`, :attr:`x1`) interval.
        """
        x0, x1 = self.intervalx
        return x0 < x < x1 or x0 > x > x1

    def fully_containsy(self, y):
        """
        Returns whether *y* is in the open (:attr:`y0`, :attr:`y1`) interval.
        """
        y0, y1 = self.intervaly
        return y0 < y < y1 or y0 > y > y1

    def fully_contains(self, x, y):
        """
        Returns whether ``x, y`` is in the bounding box, but not on its edge.
        """
        return self.fully_containsx(x) and self.fully_containsy(y)

    def fully_overlaps(self, other):
        """
        Returns whether this bounding box overlaps with the other bounding box,
        not including the edges.

        Parameters
        ----------
        other : BboxBase
        """
        ax1, ay1, ax2, ay2 = self.extents
        bx1, by1, bx2, by2 = other.extents
        if ax2 < ax1:
            ax2, ax1 = ax1, ax2
        if ay2 < ay1:
            ay2, ay1 = ay1, ay2
        if bx2 < bx1:
            bx2, bx1 = bx1, bx2
        if by2 < by1:
            by2, by1 = by1, by2
        return ax1 < bx2 and bx1 < ax2 and ay1 < by2 and by1 < ay2

    def transformed(self, transform):
        """
        Return a new :class:`Bbox` object, statically transformed by
        the given transform.
        """
        pts = self.get_points()
        ll, ul, lr = transform.transform(np.array([pts[0],
            [pts[0, 0], pts[1, 1]], [pts[1, 0], pts[0, 1]]]))
        return Bbox([ll, [lr[0], ul[1]]])

    def inverse_transformed(self, transform):
        """
        Return a new :class:`Bbox` object, statically transformed by
        the inverse of the given transform.
        """
        return self.transformed(transform.inverted())

    coefs = {'C':  (0.5, 0.5),
             'SW': (0, 0),
             'S':  (0.5, 0),
             'SE': (1.0, 0),
             'E':  (1.0, 0.5),
             'NE': (1.0, 1.0),
             'N':  (0.5, 1.0),
             'NW': (0, 1.0),
             'W':  (0, 0.5)}

    def anchored(self, c, container=None):
        """
        Return a copy of the :class:`Bbox`, shifted to position *c*
        within a container.

        Parameters
        ----------
        c :
            May be either:

            * A sequence (*cx*, *cy*) where *cx* and *cy* range from 0
              to 1, where 0 is left or bottom and 1 is right or top

            * a string:
              - 'C' for centered
              - 'S' for bottom-center
              - 'SE' for bottom-left
              - 'E' for left
              - etc.

        container : Bbox, optional
            The box within which the :class:`Bbox` is positioned; it defaults
            to the initial :class:`Bbox`.
        """
        if container is None:
            container = self
        l, b, w, h = container.bounds
        if isinstance(c, six.string_types):
            cx, cy = self.coefs[c]
        else:
            cx, cy = c
        L, B, W, H = self.bounds
        return Bbox(self._points +
                    [(l + cx * (w - W)) - L,
                     (b + cy * (h - H)) - B])

    def shrunk(self, mx, my):
        """
        Return a copy of the :class:`Bbox`, shrunk by the factor *mx*
        in the *x* direction and the factor *my* in the *y* direction.
        The lower left corner of the box remains unchanged.  Normally
        *mx* and *my* will be less than 1, but this is not enforced.
        """
        w, h = self.size
        return Bbox([self._points[0],
                    self._points[0] + [mx * w, my * h]])

    def shrunk_to_aspect(self, box_aspect, container=None, fig_aspect=1.0):
        """
        Return a copy of the :class:`Bbox`, shrunk so that it is as
        large as it can be while having the desired aspect ratio,
        *box_aspect*.  If the box coordinates are relative---that
        is, fractions of a larger box such as a figure---then the
        physical aspect ratio of that figure is specified with
        *fig_aspect*, so that *box_aspect* can also be given as a
        ratio of the absolute dimensions, not the relative dimensions.
        """
        if box_aspect <= 0 or fig_aspect <= 0:
            raise ValueError("'box_aspect' and 'fig_aspect' must be positive")
        if container is None:
            container = self
        w, h = container.size
        H = w * box_aspect / fig_aspect
        if H <= h:
            W = w
        else:
            W = h * fig_aspect / box_aspect
            H = h
        return Bbox([self._points[0],
                     self._points[0] + (W, H)])

    def splitx(self, *args):
        """
        e.g., ``bbox.splitx(f1, f2, ...)``

        Returns a list of new :class:`Bbox` objects formed by
        splitting the original one with vertical lines at fractional
        positions *f1*, *f2*, ...
        """
        xf = [0] + list(args) + [1]
        x0, y0, x1, y1 = self.extents
        w = x1 - x0
        return [Bbox([[x0 + xf0 * w, y0], [x0 + xf1 * w, y1]])
                for xf0, xf1 in zip(xf[:-1], xf[1:])]

    def splity(self, *args):
        """
        e.g., ``bbox.splitx(f1, f2, ...)``

        Returns a list of new :class:`Bbox` objects formed by
        splitting the original one with horizontal lines at fractional
        positions *f1*, *f2*, ...
        """
        yf = [0] + list(args) + [1]
        x0, y0, x1, y1 = self.extents
        h = y1 - y0
        return [Bbox([[x0, y0 + yf0 * h], [x1, y0 + yf1 * h]])
                for yf0, yf1 in zip(yf[:-1], yf[1:])]

    def count_contains(self, vertices):
        """
        Count the number of vertices contained in the :class:`Bbox`.
        Any vertices with a non-finite x or y value are ignored.

        Parameters
        ----------
        vertices : Nx2 Numpy array.
        """
        if len(vertices) == 0:
            return 0
        vertices = np.asarray(vertices)
        with np.errstate(invalid='ignore'):
            return (((self.min < vertices) &
                     (vertices < self.max)).all(axis=1).sum())

    def count_overlaps(self, bboxes):
        """
        Count the number of bounding boxes that overlap this one.

        Parameters
        ----------
        bboxes : sequence of :class:`BboxBase` objects
        """
        return count_bboxes_overlapping_bbox(
            self, np.atleast_3d([np.array(x) for x in bboxes]))

    def expanded(self, sw, sh):
        """
        Return a new :class:`Bbox` which is this :class:`Bbox`
        expanded around its center by the given factors *sw* and
        *sh*.
        """
        width = self.width
        height = self.height
        deltaw = (sw * width - width) / 2.0
        deltah = (sh * height - height) / 2.0
        a = np.array([[-deltaw, -deltah], [deltaw, deltah]])
        return Bbox(self._points + a)

    def padded(self, p):
        """
        Return a new :class:`Bbox` that is padded on all four sides by
        the given value.
        """
        points = self.get_points()
        return Bbox(points + [[-p, -p], [p, p]])

    def translated(self, tx, ty):
        """
        Return a copy of the :class:`Bbox`, statically translated by
        *tx* and *ty*.
        """
        return Bbox(self._points + (tx, ty))

    def corners(self):
        """
        Return an array of points which are the four corners of this
        rectangle.  For example, if this :class:`Bbox` is defined by
        the points (*a*, *b*) and (*c*, *d*), :meth:`corners` returns
        (*a*, *b*), (*a*, *d*), (*c*, *b*) and (*c*, *d*).
        """
        l, b, r, t = self.get_points().flatten()
        return np.array([[l, b], [l, t], [r, b], [r, t]])

    def rotated(self, radians):
        """
        Return a new bounding box that bounds a rotated version of
        this bounding box by the given radians.  The new bounding box
        is still aligned with the axes, of course.
        """
        corners = self.corners()
        corners_rotated = Affine2D().rotate(radians).transform(corners)
        bbox = Bbox.unit()
        bbox.update_from_data_xy(corners_rotated, ignore=True)
        return bbox

    @staticmethod
    def union(bboxes):
        """
        Return a :class:`Bbox` that contains all of the given bboxes.
        """
        if not len(bboxes):
            raise ValueError("'bboxes' cannot be empty")
        x0 = np.min([bbox.xmin for bbox in bboxes])
        x1 = np.max([bbox.xmax for bbox in bboxes])
        y0 = np.min([bbox.ymin for bbox in bboxes])
        y1 = np.max([bbox.ymax for bbox in bboxes])
        return Bbox([[x0, y0], [x1, y1]])

    @staticmethod
    def intersection(bbox1, bbox2):
        """
        Return the intersection of the two bboxes or None
        if they do not intersect.
        """
        x0 = np.maximum(bbox1.xmin, bbox2.xmin)
        x1 = np.minimum(bbox1.xmax, bbox2.xmax)
        y0 = np.maximum(bbox1.ymin, bbox2.ymin)
        y1 = np.minimum(bbox1.ymax, bbox2.ymax)
        return Bbox([[x0, y0], [x1, y1]]) if x0 <= x1 and y0 <= y1 else None


class Bbox(BboxBase):
    """
    A mutable bounding box.
    """

    def __init__(self, points, **kwargs):
        """
        Parameters
        ----------
        points : ndarray
            A 2x2 numpy array of the form ``[[x0, y0], [x1, y1]]``.

        Notes
        -----
        If you need to create a :class:`Bbox` object from another form
        of data, consider the static methods :meth:`unit`,
        :meth:`from_bounds` and :meth:`from_extents`.
        """
        BboxBase.__init__(self, **kwargs)
        points = np.asarray(points, float)
        if points.shape != (2, 2):
            raise ValueError('Bbox points must be of the form '
                             '"[[x0, y0], [x1, y1]]".')
        self._points = points
        self._minpos = np.array([np.inf, np.inf])
        self._ignore = True
        # it is helpful in some contexts to know if the bbox is a
        # default or has been mutated; we store the orig points to
        # support the mutated methods
        self._points_orig = self._points.copy()
    if DEBUG:
        ___init__ = __init__

        def __init__(self, points, **kwargs):
            self._check(points)
            self.___init__(points, **kwargs)

        def invalidate(self):
            self._check(self._points)
            TransformNode.invalidate(self)

    @staticmethod
    def unit():
        """
        (staticmethod) Create a new unit :class:`Bbox` from (0, 0) to
        (1, 1).
        """
        return Bbox(np.array([[0.0, 0.0], [1.0, 1.0]], float))

    @staticmethod
    def null():
        """
        (staticmethod) Create a new null :class:`Bbox` from (inf, inf) to
        (-inf, -inf).
        """
        return Bbox(np.array([[np.inf, np.inf], [-np.inf, -np.inf]], float))

    @staticmethod
    def from_bounds(x0, y0, width, height):
        """
        (staticmethod) Create a new :class:`Bbox` from *x0*, *y0*,
        *width* and *height*.

        *width* and *height* may be negative.
        """
        return Bbox.from_extents(x0, y0, x0 + width, y0 + height)

    @staticmethod
    def from_extents(*args):
        """
        (staticmethod) Create a new Bbox from *left*, *bottom*,
        *right* and *top*.

        The *y*-axis increases upwards.
        """
        points = np.array(args, dtype=float).reshape(2, 2)
        return Bbox(points)

    def __format__(self, fmt):
        return (
            'Bbox(x0={0.x0:{1}}, y0={0.y0:{1}}, x1={0.x1:{1}}, y1={0.y1:{1}})'.
            format(self, fmt))

    def __str__(self):
        return format(self, '')

    def __repr__(self):
        return 'Bbox([[{0.x0}, {0.y0}], [{0.x1}, {0.y1}]])'.format(self)

    def ignore(self, value):
        """
        Set whether the existing bounds of the box should be ignored
        by subsequent calls to :meth:`update_from_data_xy`.

        value : bool
           - When ``True``, subsequent calls to :meth:`update_from_data_xy`
             will ignore the existing bounds of the :class:`Bbox`.

           - When ``False``, subsequent calls to :meth:`update_from_data_xy`
             will include the existing bounds of the :class:`Bbox`.
        """
        self._ignore = value

    def update_from_path(self, path, ignore=None, updatex=True, updatey=True):
        """
        Update the bounds of the :class:`Bbox` based on the passed in
        data.  After updating, the bounds will have positive *width*
        and *height*; *x0* and *y0* will be the minimal values.

        Parameters
        ----------
        path : :class:`~matplotlib.path.Path`

        ignore : bool, optional
           - when ``True``, ignore the existing bounds of the :class:`Bbox`.
           - when ``False``, include the existing bounds of the :class:`Bbox`.
           - when ``None``, use the last value passed to :meth:`ignore`.

        updatex, updatey : bool, optional
            When ``True``, update the x/y values.
        """
        if ignore is None:
            ignore = self._ignore

        if path.vertices.size == 0:
            return

        points, minpos, changed = update_path_extents(
            path, None, self._points, self._minpos, ignore)

        if changed:
            self.invalidate()
            if updatex:
                self._points[:, 0] = points[:, 0]
                self._minpos[0] = minpos[0]
            if updatey:
                self._points[:, 1] = points[:, 1]
                self._minpos[1] = minpos[1]

    def update_from_data_xy(self, xy, ignore=None, updatex=True, updatey=True):
        """
        Update the bounds of the :class:`Bbox` based on the passed in
        data.  After updating, the bounds will have positive *width*
        and *height*; *x0* and *y0* will be the minimal values.

        Parameters
        ----------
        xy : ndarray
            A numpy array of 2D points.

        ignore : bool, optional
           - When ``True``, ignore the existing bounds of the :class:`Bbox`.
           - When ``False``, include the existing bounds of the :class:`Bbox`.
           - When ``None``, use the last value passed to :meth:`ignore`.

        updatex, updatey : bool, optional
            When ``True``, update the x/y values.
        """
        if len(xy) == 0:
            return

        path = Path(xy)
        self.update_from_path(path, ignore=ignore,
                                    updatex=updatex, updatey=updatey)

    @BboxBase.x0.setter
    def x0(self, val):
        self._points[0, 0] = val
        self.invalidate()

    @BboxBase.y0.setter
    def y0(self, val):
        self._points[0, 1] = val
        self.invalidate()

    @BboxBase.x1.setter
    def x1(self, val):
        self._points[1, 0] = val
        self.invalidate()

    @BboxBase.y1.setter
    def y1(self, val):
        self._points[1, 1] = val
        self.invalidate()

    @BboxBase.p0.setter
    def p0(self, val):
        self._points[0] = val
        self.invalidate()

    @BboxBase.p1.setter
    def p1(self, val):
        self._points[1] = val
        self.invalidate()

    @BboxBase.intervalx.setter
    def intervalx(self, interval):
        self._points[:, 0] = interval
        self.invalidate()

    @BboxBase.intervaly.setter
    def intervaly(self, interval):
        self._points[:, 1] = interval
        self.invalidate()

    @BboxBase.bounds.setter
    def bounds(self, bounds):
        l, b, w, h = bounds
        points = np.array([[l, b], [l + w, b + h]], float)
        if np.any(self._points != points):
            self._points = points
            self.invalidate()

    @property
    def minpos(self):
        return self._minpos

    @property
    def minposx(self):
        return self._minpos[0]

    @property
    def minposy(self):
        return self._minpos[1]

    def get_points(self):
        """
        Get the points of the bounding box directly as a numpy array
        of the form: ``[[x0, y0], [x1, y1]]``.
        """
        self._invalid = 0
        return self._points

    def set_points(self, points):
        """
        Set the points of the bounding box directly from a numpy array
        of the form: ``[[x0, y0], [x1, y1]]``.  No error checking is
        performed, as this method is mainly for internal use.
        """
        if np.any(self._points != points):
            self._points = points
            self.invalidate()

    def set(self, other):
        """
        Set this bounding box from the "frozen" bounds of another
        :class:`Bbox`.
        """
        if np.any(self._points != other.get_points()):
            self._points = other.get_points()
            self.invalidate()

    def mutated(self):
        'Return whether the bbox has changed since init.'
        return self.mutatedx() or self.mutatedy()

    def mutatedx(self):
        'Return whether the x-limits have changed since init.'
        return (self._points[0, 0] != self._points_orig[0, 0] or
                self._points[1, 0] != self._points_orig[1, 0])

    def mutatedy(self):
        'Return whether the y-limits have changed since init.'
        return (self._points[0, 1] != self._points_orig[0, 1] or
                self._points[1, 1] != self._points_orig[1, 1])


class TransformedBbox(BboxBase):
    """
    A :class:`Bbox` that is automatically transformed by a given
    transform.  When either the child bounding box or transform
    changes, the bounds of this bbox will update accordingly.
    """
    def __init__(self, bbox, transform, **kwargs):
        """
        Parameters
        ----------
        bbox : :class:`Bbox`

        transform : :class:`Transform`
        """
        if not bbox.is_bbox:
            raise ValueError("'bbox' is not a bbox")
        if not isinstance(transform, Transform):
            raise ValueError("'transform' must be an instance of "
                             "'matplotlib.transform.Transform'")
        if transform.input_dims != 2 or transform.output_dims != 2:
            raise ValueError(
                "The input and output dimensions of 'transform' must be 2")

        BboxBase.__init__(self, **kwargs)
        self._bbox = bbox
        self._transform = transform
        self.set_children(bbox, transform)
        self._points = None

    def __str__(self):
        return ("{}(\n"
                    "{},\n"
                    "{})"
                .format(type(self).__name__,
                        _indent_str(self._bbox),
                        _indent_str(self._transform)))

    __repr__ = __str__

    def get_points(self):
        if self._invalid:
            p = self._bbox.get_points()
            # Transform all four points, then make a new bounding box
            # from the result, taking care to make the orientation the
            # same.
            points = self._transform.transform(
                [[p[0, 0], p[0, 1]],
                 [p[1, 0], p[0, 1]],
                 [p[0, 0], p[1, 1]],
                 [p[1, 0], p[1, 1]]])
            points = np.ma.filled(points, 0.0)

            xs = min(points[:, 0]), max(points[:, 0])
            if p[0, 0] > p[1, 0]:
                xs = xs[::-1]

            ys = min(points[:, 1]), max(points[:, 1])
            if p[0, 1] > p[1, 1]:
                ys = ys[::-1]

            self._points = np.array([
                [xs[0], ys[0]],
                [xs[1], ys[1]]
            ])

            self._invalid = 0
        return self._points
    get_points.__doc__ = Bbox.get_points.__doc__

    if DEBUG:
        _get_points = get_points

        def get_points(self):
            points = self._get_points()
            self._check(points)
            return points


class LockableBbox(BboxBase):
    """
    A :class:`Bbox` where some elements may be locked at certain values.

    When the child bounding box changes, the bounds of this bbox will update
    accordingly with the exception of the locked elements.
    """
    def __init__(self, bbox, x0=None, y0=None, x1=None, y1=None, **kwargs):
        """
        Parameters
        ----------
        bbox : Bbox
            The child bounding box to wrap.

        x0 : float or None
            The locked value for x0, or None to leave unlocked.

        y0 : float or None
            The locked value for y0, or None to leave unlocked.

        x1 : float or None
            The locked value for x1, or None to leave unlocked.

        y1 : float or None
            The locked value for y1, or None to leave unlocked.

        """
        if not bbox.is_bbox:
            raise ValueError("'bbox' is not a bbox")

        BboxBase.__init__(self, **kwargs)
        self._bbox = bbox
        self.set_children(bbox)
        self._points = None
        fp = [x0, y0, x1, y1]
        mask = [val is None for val in fp]
        self._locked_points = np.ma.array(fp, float, mask=mask).reshape((2, 2))

    def __str__(self):
        return ("{}(\n"
                    "{},\n"
                    "{})"
                .format(type(self).__name__,
                        _indent_str(self._bbox),
                        _indent_str(self._locked_points)))

    __repr__ = __str__

    def get_points(self):
        if self._invalid:
            points = self._bbox.get_points()
            self._points = np.where(self._locked_points.mask,
                                    points,
                                    self._locked_points)
            self._invalid = 0
        return self._points
    get_points.__doc__ = Bbox.get_points.__doc__

    if DEBUG:
        _get_points = get_points

        def get_points(self):
            points = self._get_points()
            self._check(points)
            return points

    @property
    def locked_x0(self):
        """
        float or None: The value used for the locked x0.
        """
        if self._locked_points.mask[0, 0]:
            return None
        else:
            return self._locked_points[0, 0]

    @locked_x0.setter
    def locked_x0(self, x0):
        self._locked_points.mask[0, 0] = x0 is None
        self._locked_points.data[0, 0] = x0
        self.invalidate()

    @property
    def locked_y0(self):
        """
        float or None: The value used for the locked y0.
        """
        if self._locked_points.mask[0, 1]:
            return None
        else:
            return self._locked_points[0, 1]

    @locked_y0.setter
    def locked_y0(self, y0):
        self._locked_points.mask[0, 1] = y0 is None
        self._locked_points.data[0, 1] = y0
        self.invalidate()

    @property
    def locked_x1(self):
        """
        float or None: The value used for the locked x1.
        """
        if self._locked_points.mask[1, 0]:
            return None
        else:
            return self._locked_points[1, 0]

    @locked_x1.setter
    def locked_x1(self, x1):
        self._locked_points.mask[1, 0] = x1 is None
        self._locked_points.data[1, 0] = x1
        self.invalidate()

    @property
    def locked_y1(self):
        """
        float or None: The value used for the locked y1.
        """
        if self._locked_points.mask[1, 1]:
            return None
        else:
            return self._locked_points[1, 1]

    @locked_y1.setter
    def locked_y1(self, y1):
        self._locked_points.mask[1, 1] = y1 is None
        self._locked_points.data[1, 1] = y1
        self.invalidate()


class Transform(TransformNode):
    """
    The base class of all :class:`TransformNode` instances that
    actually perform a transformation.

    All non-affine transformations should be subclasses of this class.
    New affine transformations should be subclasses of
    :class:`Affine2D`.

    Subclasses of this class should override the following members (at
    minimum):

      - :attr:`input_dims`
      - :attr:`output_dims`
      - :meth:`transform`
      - :attr:`is_separable`
      - :attr:`has_inverse`
      - :meth:`inverted` (if :attr:`has_inverse` is True)

    If the transform needs to do something non-standard with
    :class:`matplotlib.path.Path` objects, such as adding curves
    where there were once line segments, it should override:

      - :meth:`transform_path`
    """
    input_dims = None
    """
    The number of input dimensions of this transform.
    Must be overridden (with integers) in the subclass.
    """

    output_dims = None
    """
    The number of output dimensions of this transform.
    Must be overridden (with integers) in the subclass.
    """

    has_inverse = False
    """True if this transform has a corresponding inverse transform."""

    is_separable = False
    """True if this transform is separable in the x- and y- dimensions."""

    def __add__(self, other):
        """
        Composes two transforms together such that *self* is followed
        by *other*.
        """
        if isinstance(other, Transform):
            return composite_transform_factory(self, other)
        raise TypeError(
            "Can not add Transform to object of type '%s'" % type(other))

    def __radd__(self, other):
        """
        Composes two transforms together such that *self* is followed
        by *other*.
        """
        if isinstance(other, Transform):
            return composite_transform_factory(other, self)
        raise TypeError(
            "Can not add Transform to object of type '%s'" % type(other))

    # Equality is based on object identity for `Transform`s (so we don't
    # override `__eq__`), but some subclasses, such as TransformWrapper &
    # AffineBase, override this behavior.

    if six.PY2:
        def __ne__(self, other):
            return not (self == other)

    def _iter_break_from_left_to_right(self):
        """
        Returns an iterator breaking down this transform stack from left to
        right recursively. If self == ((A, N), A) then the result will be an
        iterator which yields I : ((A, N), A), followed by A : (N, A),
        followed by (A, N) : (A), but not ((A, N), A) : I.

        This is equivalent to flattening the stack then yielding
        ``flat_stack[:i], flat_stack[i:]`` where i=0..(n-1).

        """
        yield IdentityTransform(), self

    @property
    def depth(self):
        """
        Returns the number of transforms which have been chained
        together to form this Transform instance.

        .. note::

            For the special case of a Composite transform, the maximum depth
            of the two is returned.

        """
        return 1

    def contains_branch(self, other):
        """
        Return whether the given transform is a sub-tree of this transform.

        This routine uses transform equality to identify sub-trees, therefore
        in many situations it is object id which will be used.

        For the case where the given transform represents the whole
        of this transform, returns True.

        """
        if self.depth < other.depth:
            return False

        # check that a subtree is equal to other (starting from self)
        for _, sub_tree in self._iter_break_from_left_to_right():
            if sub_tree == other:
                return True
        return False

    def contains_branch_seperately(self, other_transform):
        """
        Returns whether the given branch is a sub-tree of this transform on
        each separate dimension.

        A common use for this method is to identify if a transform is a blended
        transform containing an axes' data transform. e.g.::

            x_isdata, y_isdata = trans.contains_branch_seperately(ax.transData)

        """
        if self.output_dims != 2:
            raise ValueError('contains_branch_seperately only supports '
                             'transforms with 2 output dimensions')
        # for a non-blended transform each separate dimension is the same, so
        # just return the appropriate shape.
        return [self.contains_branch(other_transform)] * 2

    def __sub__(self, other):
        """
        Returns a transform stack which goes all the way down self's transform
        stack, and then ascends back up other's stack. If it can, this is
        optimised::

            # normally
            A - B == a + b.inverted()

            # sometimes, when A contains the tree B there is no need to
            # descend all the way down to the base of A (via B), instead we
            # can just stop at B.

            (A + B) - (B)^-1 == A

            # similarly, when B contains tree A, we can avoid decending A at
            # all, basically:
            A - (A + B) == ((B + A) - A).inverted() or B^-1

        For clarity, the result of ``(A + B) - B + B == (A + B)``.

        """
        # we only know how to do this operation if other is a Transform.
        if not isinstance(other, Transform):
            return NotImplemented

        for remainder, sub_tree in self._iter_break_from_left_to_right():
            if sub_tree == other:
                return remainder

        for remainder, sub_tree in other._iter_break_from_left_to_right():
            if sub_tree == self:
                if not remainder.has_inverse:
                    raise ValueError("The shortcut cannot be computed since "
                     "other's transform includes a non-invertable component.")
                return remainder.inverted()

        # if we have got this far, then there was no shortcut possible
        if other.has_inverse:
            return self + other.inverted()
        else:
            raise ValueError('It is not possible to compute transA - transB '
                             'since transB cannot be inverted and there is no '
                             'shortcut possible.')

    def __array__(self, *args, **kwargs):
        """
        Array interface to get at this Transform's affine matrix.
        """
        return self.get_affine().get_matrix()

    def transform(self, values):
        """
        Performs the transformation on the given array of values.

        Accepts a numpy array of shape (N x :attr:`input_dims`) and
        returns a numpy array of shape (N x :attr:`output_dims`).

        Alternatively, accepts a numpy array of length :attr:`input_dims`
        and returns a numpy array of length :attr:`output_dims`.
        """
        # Ensure that values is a 2d array (but remember whether
        # we started with a 1d or 2d array).
        values = np.asanyarray(values)
        ndim = values.ndim
        values = values.reshape((-1, self.input_dims))

        # Transform the values
        res = self.transform_affine(self.transform_non_affine(values))

        # Convert the result back to the shape of the input values.
        if ndim == 0:
            assert not np.ma.is_masked(res)  # just to be on the safe side
            return res[0, 0]
        if ndim == 1:
            return res.reshape(-1)
        elif ndim == 2:
            return res
        raise ValueError(
            "Input values must have shape (N x {dims}) "
            "or ({dims}).".format(dims=self.input_dims))

    def transform_affine(self, values):
        """
        Performs only the affine part of this transformation on the
        given array of values.

        ``transform(values)`` is always equivalent to
        ``transform_affine(transform_non_affine(values))``.

        In non-affine transformations, this is generally a no-op.  In
        affine transformations, this is equivalent to
        ``transform(values)``.

        Accepts a numpy array of shape (N x :attr:`input_dims`) and
        returns a numpy array of shape (N x :attr:`output_dims`).

        Alternatively, accepts a numpy array of length :attr:`input_dims`
        and returns a numpy array of length :attr:`output_dims`.
        """
        return self.get_affine().transform(values)

    def transform_non_affine(self, values):
        """
        Performs only the non-affine part of the transformation.

        ``transform(values)`` is always equivalent to
        ``transform_affine(transform_non_affine(values))``.

        In non-affine transformations, this is generally equivalent to
        ``transform(values)``.  In affine transformations, this is
        always a no-op.

        Accepts a numpy array of shape (N x :attr:`input_dims`) and
        returns a numpy array of shape (N x :attr:`output_dims`).

        Alternatively, accepts a numpy array of length :attr:`input_dims`
        and returns a numpy array of length :attr:`output_dims`.
        """
        return values

    def transform_bbox(self, bbox):
        """
        Transform the given bounding box.

        Note, for smarter transforms including caching (a common
        requirement for matplotlib figures), see :class:`TransformedBbox`.
        """
        return Bbox(self.transform(bbox.get_points()))

    def get_affine(self):
        """
        Get the affine part of this transform.
        """
        return IdentityTransform()

    def get_matrix(self):
        """
        Get the Affine transformation array for the affine part
        of this transform.

        """
        return self.get_affine().get_matrix()

    def transform_point(self, point):
        """
        A convenience function that returns the transformed copy of a
        single point.

        The point is given as a sequence of length :attr:`input_dims`.
        The transformed point is returned as a sequence of length
        :attr:`output_dims`.
        """
        if len(point) != self.input_dims:
            raise ValueError("The length of 'point' must be 'self.input_dims'")
        return self.transform(np.asarray([point]))[0]

    def transform_path(self, path):
        """
        Returns a transformed path.

        *path*: a :class:`~matplotlib.path.Path` instance.

        In some cases, this transform may insert curves into the path
        that began as line segments.
        """
        return self.transform_path_affine(self.transform_path_non_affine(path))

    def transform_path_affine(self, path):
        """
        Returns a path, transformed only by the affine part of
        this transform.

        *path*: a :class:`~matplotlib.path.Path` instance.

        ``transform_path(path)`` is equivalent to
        ``transform_path_affine(transform_path_non_affine(values))``.
        """
        return self.get_affine().transform_path_affine(path)

    def transform_path_non_affine(self, path):
        """
        Returns a path, transformed only by the non-affine
        part of this transform.

        *path*: a :class:`~matplotlib.path.Path` instance.

        ``transform_path(path)`` is equivalent to
        ``transform_path_affine(transform_path_non_affine(values))``.
        """
        x = self.transform_non_affine(path.vertices)
        return Path._fast_from_codes_and_verts(x, path.codes,
                {'interpolation_steps': path._interpolation_steps,
                 'should_simplify': path.should_simplify})

    def transform_angles(self, angles, pts, radians=False, pushoff=1e-5):
        """
        Performs transformation on a set of angles anchored at
        specific locations.

        The *angles* must be a column vector (i.e., numpy array).

        The *pts* must be a two-column numpy array of x,y positions
        (angle transforms currently only work in 2D).  This array must
        have the same number of rows as *angles*.

        *radians* indicates whether or not input angles are given in
         radians (True) or degrees (False; the default).

        *pushoff* is the distance to move away from *pts* for
         determining transformed angles (see discussion of method
         below).

        The transformed angles are returned in an array with the same
        size as *angles*.

        The generic version of this method uses a very generic
        algorithm that transforms *pts*, as well as locations very
        close to *pts*, to find the angle in the transformed system.
        """
        # Must be 2D
        if self.input_dims != 2 or self.output_dims != 2:
            raise NotImplementedError('Only defined in 2D')

        if pts.shape[1] != 2:
            raise ValueError("'pts' must be array with 2 columns for x,y")

        if angles.ndim != 1 or angles.shape[0] != pts.shape[0]:
            raise ValueError("'angles' must be a column vector and have same "
                             "number of rows as 'pts'")

        # Convert to radians if desired
        if not radians:
            angles = angles / 180.0 * np.pi

        # Move a short distance away
        pts2 = pts + pushoff * np.c_[np.cos(angles), np.sin(angles)]

        # Transform both sets of points
        tpts = self.transform(pts)
        tpts2 = self.transform(pts2)

        # Calculate transformed angles
        d = tpts2 - tpts
        a = np.arctan2(d[:, 1], d[:, 0])

        # Convert back to degrees if desired
        if not radians:
            a = np.rad2deg(a)

        return a

    def inverted(self):
        """
        Return the corresponding inverse transformation.

        The return value of this method should be treated as
        temporary.  An update to *self* does not cause a corresponding
        update to its inverted copy.

        ``x === self.inverted().transform(self.transform(x))``
        """
        raise NotImplementedError()

    def __repr__(self):
        return str(self)


class TransformWrapper(Transform):
    """
    A helper class that holds a single child transform and acts
    equivalently to it.

    This is useful if a node of the transform tree must be replaced at
    run time with a transform of a different type.  This class allows
    that replacement to correctly trigger invalidation.

    Note that :class:`TransformWrapper` instances must have the same
    input and output dimensions during their entire lifetime, so the
    child transform may only be replaced with another child transform
    of the same dimensions.
    """
    pass_through = True

    def __init__(self, child):
        """
        *child*: A class:`Transform` instance.  This child may later
        be replaced with :meth:`set`.
        """
        if not isinstance(child, Transform):
            raise ValueError("'child' must be an instance of "
                             "'matplotlib.transform.Transform'")
        self._init(child)
        self.set_children(child)

    def _init(self, child):
        Transform.__init__(self)
        self.input_dims = child.input_dims
        self.output_dims = child.output_dims
        self._set(child)
        self._invalid = 0

    def __eq__(self, other):
        return self._child.__eq__(other)

    # NOTE: Transform.__[gs]etstate__ should be sufficient when using only
    # Python 3.4+.
    def __getstate__(self):
        # only store the child information and parents
        return {
            'child': self._child,
            'input_dims': self.input_dims,
            'output_dims': self.output_dims,
            # turn the weak-values dictionary into a normal dictionary
            'parents': dict((k, v()) for (k, v) in
                            six.iteritems(self._parents))
        }

    def __setstate__(self, state):
        # re-initialise the TransformWrapper with the state's child
        self._init(state['child'])
        # The child may not be unpickled yet, so restore its information.
        self.input_dims = state['input_dims']
        self.output_dims = state['output_dims']
        # turn the normal dictionary back into a dictionary with weak
        # values
        self._parents = dict((k, weakref.ref(v)) for (k, v) in
                             six.iteritems(state['parents']) if v is not None)

    def __str__(self):
        return ("{}(\n"
                    "{})"
                .format(type(self).__name__,
                        _indent_str(self._child)))

    def frozen(self):
        return self._child.frozen()
    frozen.__doc__ = Transform.frozen.__doc__

    def _set(self, child):
        self._child = child

        self.transform = child.transform
        self.transform_affine = child.transform_affine
        self.transform_non_affine = child.transform_non_affine
        self.transform_path = child.transform_path
        self.transform_path_affine = child.transform_path_affine
        self.transform_path_non_affine = child.transform_path_non_affine
        self.get_affine = child.get_affine
        self.inverted = child.inverted
        self.get_matrix = child.get_matrix

        # note we do not wrap other properties here since the transform's
        # child can be changed with WrappedTransform.set and so checking
        # is_affine and other such properties may be dangerous.

    def set(self, child):
        """
        Replace the current child of this transform with another one.

        The new child must have the same number of input and output
        dimensions as the current child.
        """
        if (child.input_dims != self.input_dims or
                child.output_dims != self.output_dims):
            raise ValueError(
                "The new child must have the same number of input and output "
                "dimensions as the current child")

        self.set_children(child)
        self._set(child)

        self._invalid = 0
        self.invalidate()
        self._invalid = 0

    def _get_is_affine(self):
        return self._child.is_affine
    is_affine = property(_get_is_affine)

    def _get_is_separable(self):
        return self._child.is_separable
    is_separable = property(_get_is_separable)

    def _get_has_inverse(self):
        return self._child.has_inverse
    has_inverse = property(_get_has_inverse)


class AffineBase(Transform):
    """
    The base class of all affine transformations of any number of
    dimensions.
    """
    is_affine = True

    def __init__(self, *args, **kwargs):
        Transform.__init__(self, *args, **kwargs)
        self._inverted = None

    def __array__(self, *args, **kwargs):
        # optimises the access of the transform matrix vs the superclass
        return self.get_matrix()

    @staticmethod
    def _concat(a, b):
        """
        Concatenates two transformation matrices (represented as numpy
        arrays) together.
        """
        return np.dot(b, a)

    def __eq__(self, other):
        if getattr(other, "is_affine", False):
            return np.all(self.get_matrix() == other.get_matrix())
        return NotImplemented

    def transform(self, values):
        return self.transform_affine(values)
    transform.__doc__ = Transform.transform.__doc__

    def transform_affine(self, values):
        raise NotImplementedError('Affine subclasses should override this '
                                  'method.')
    transform_affine.__doc__ = Transform.transform_affine.__doc__

    def transform_non_affine(self, points):
        return points
    transform_non_affine.__doc__ = Transform.transform_non_affine.__doc__

    def transform_path(self, path):
        return self.transform_path_affine(path)
    transform_path.__doc__ = Transform.transform_path.__doc__

    def transform_path_affine(self, path):
        return Path(self.transform_affine(path.vertices),
                    path.codes, path._interpolation_steps)
    transform_path_affine.__doc__ = Transform.transform_path_affine.__doc__

    def transform_path_non_affine(self, path):
        return path
    transform_path_non_affine.__doc__ = Transform.transform_path_non_affine.__doc__

    def get_affine(self):
        return self
    get_affine.__doc__ = Transform.get_affine.__doc__


class Affine2DBase(AffineBase):
    """
    The base class of all 2D affine transformations.

    2D affine transformations are performed using a 3x3 numpy array::

        a c e
        b d f
        0 0 1

    This class provides the read-only interface.  For a mutable 2D
    affine transformation, use :class:`Affine2D`.

    Subclasses of this class will generally only need to override a
    constructor and :meth:`get_matrix` that generates a custom 3x3 matrix.
    """
    has_inverse = True

    input_dims = 2
    output_dims = 2

    def frozen(self):
        return Affine2D(self.get_matrix().copy())
    frozen.__doc__ = AffineBase.frozen.__doc__

    def _get_is_separable(self):
        mtx = self.get_matrix()
        return mtx[0, 1] == 0.0 and mtx[1, 0] == 0.0
    is_separable = property(_get_is_separable)

    def to_values(self):
        """
        Return the values of the matrix as a sequence (a,b,c,d,e,f)
        """
        mtx = self.get_matrix()
        return tuple(mtx[:2].swapaxes(0, 1).flatten())

    @staticmethod
    def matrix_from_values(a, b, c, d, e, f):
        """
        (staticmethod) Create a new transformation matrix as a 3x3
        numpy array of the form::

          a c e
          b d f
          0 0 1
        """
        return np.array([[a, c, e], [b, d, f], [0.0, 0.0, 1.0]], float)

    def transform_affine(self, points):
        mtx = self.get_matrix()
        if isinstance(points, np.ma.MaskedArray):
            tpoints = affine_transform(points.data, mtx)
            return np.ma.MaskedArray(tpoints, mask=np.ma.getmask(points))
        return affine_transform(points, mtx)

    def transform_point(self, point):
        mtx = self.get_matrix()
        return affine_transform([point], mtx)[0]
    transform_point.__doc__ = AffineBase.transform_point.__doc__

    if DEBUG:
        _transform_affine = transform_affine

        def transform_affine(self, points):
            # The major speed trap here is just converting to the
            # points to an array in the first place.  If we can use
            # more arrays upstream, that should help here.
            if not isinstance(points, (np.ma.MaskedArray, np.ndarray)):
                warnings.warn(
                    ('A non-numpy array of type %s was passed in for ' +
                     'transformation.  Please correct this.')
                    % type(points))
            return self._transform_affine(points)
    transform_affine.__doc__ = AffineBase.transform_affine.__doc__

    def inverted(self):
        if self._inverted is None or self._invalid:
            mtx = self.get_matrix()
            shorthand_name = None
            if self._shorthand_name:
                shorthand_name = '(%s)-1' % self._shorthand_name
            self._inverted = Affine2D(inv(mtx), shorthand_name=shorthand_name)
            self._invalid = 0
        return self._inverted
    inverted.__doc__ = AffineBase.inverted.__doc__


class Affine2D(Affine2DBase):
    """
    A mutable 2D affine transformation.
    """

    def __init__(self, matrix=None, **kwargs):
        """
        Initialize an Affine transform from a 3x3 numpy float array::

          a c e
          b d f
          0 0 1

        If *matrix* is None, initialize with the identity transform.
        """
        Affine2DBase.__init__(self, **kwargs)
        if matrix is None:
            # A bit faster than np.identity(3).
            matrix = IdentityTransform._mtx.copy()
        self._mtx = matrix
        self._invalid = 0

    def __str__(self):
        return ("{}(\n"
                    "{})"
                .format(type(self).__name__,
                        _indent_str(self._mtx)))

    @staticmethod
    def from_values(a, b, c, d, e, f):
        """
        (staticmethod) Create a new Affine2D instance from the given
        values::

          a c e
          b d f
          0 0 1

        .
        """
        return Affine2D(
            np.array([a, c, e, b, d, f, 0.0, 0.0, 1.0], float).reshape((3, 3)))

    def get_matrix(self):
        """
        Get the underlying transformation matrix as a 3x3 numpy array::

          a c e
          b d f
          0 0 1

        .
        """
        self._invalid = 0
        return self._mtx

    def set_matrix(self, mtx):
        """
        Set the underlying transformation matrix from a 3x3 numpy array::

          a c e
          b d f
          0 0 1

        .
        """
        self._mtx = mtx
        self.invalidate()

    def set(self, other):
        """
        Set this transformation from the frozen copy of another
        :class:`Affine2DBase` object.
        """
        if not isinstance(other, Affine2DBase):
            raise ValueError("'other' must be an instance of "
                             "'matplotlib.transform.Affine2DBase'")
        self._mtx = other.get_matrix()
        self.invalidate()

    @staticmethod
    def identity():
        """
        (staticmethod) Return a new :class:`Affine2D` object that is
        the identity transform.

        Unless this transform will be mutated later on, consider using
        the faster :class:`IdentityTransform` class instead.
        """
        return Affine2D()

    def clear(self):
        """
        Reset the underlying matrix to the identity transform.
        """
        # A bit faster than np.identity(3).
        self._mtx = IdentityTransform._mtx.copy()
        self.invalidate()
        return self

    def rotate(self, theta):
        """
        Add a rotation (in radians) to this transform in place.

        Returns *self*, so this method can easily be chained with more
        calls to :meth:`rotate`, :meth:`rotate_deg`, :meth:`translate`
        and :meth:`scale`.
        """
        a = np.cos(theta)
        b = np.sin(theta)
        rotate_mtx = np.array([[a, -b, 0.0], [b, a, 0.0], [0.0, 0.0, 1.0]],
                              float)
        self._mtx = np.dot(rotate_mtx, self._mtx)
        self.invalidate()
        return self

    def rotate_deg(self, degrees):
        """
        Add a rotation (in degrees) to this transform in place.

        Returns *self*, so this method can easily be chained with more
        calls to :meth:`rotate`, :meth:`rotate_deg`, :meth:`translate`
        and :meth:`scale`.
        """
        return self.rotate(np.deg2rad(degrees))

    def rotate_around(self, x, y, theta):
        """
        Add a rotation (in radians) around the point (x, y) in place.

        Returns *self*, so this method can easily be chained with more
        calls to :meth:`rotate`, :meth:`rotate_deg`, :meth:`translate`
        and :meth:`scale`.
        """
        return self.translate(-x, -y).rotate(theta).translate(x, y)

    def rotate_deg_around(self, x, y, degrees):
        """
        Add a rotation (in degrees) around the point (x, y) in place.

        Returns *self*, so this method can easily be chained with more
        calls to :meth:`rotate`, :meth:`rotate_deg`, :meth:`translate`
        and :meth:`scale`.
        """
        # Cast to float to avoid wraparound issues with uint8's
        x, y = float(x), float(y)
        return self.translate(-x, -y).rotate_deg(degrees).translate(x, y)

    def translate(self, tx, ty):
        """
        Adds a translation in place.

        Returns *self*, so this method can easily be chained with more
        calls to :meth:`rotate`, :meth:`rotate_deg`, :meth:`translate`
        and :meth:`scale`.
        """
        translate_mtx = np.array(
            [[1.0, 0.0, tx], [0.0, 1.0, ty], [0.0, 0.0, 1.0]], float)
        self._mtx = np.dot(translate_mtx, self._mtx)
        self.invalidate()
        return self

    def scale(self, sx, sy=None):
        """
        Adds a scale in place.

        If *sy* is None, the same scale is applied in both the *x*- and
        *y*-directions.

        Returns *self*, so this method can easily be chained with more
        calls to :meth:`rotate`, :meth:`rotate_deg`, :meth:`translate`
        and :meth:`scale`.
        """
        if sy is None:
            sy = sx
        scale_mtx = np.array(
            [[sx, 0.0, 0.0], [0.0, sy, 0.0], [0.0, 0.0, 1.0]], float)
        self._mtx = np.dot(scale_mtx, self._mtx)
        self.invalidate()
        return self

    def skew(self, xShear, yShear):
        """
        Adds a skew in place.

        *xShear* and *yShear* are the shear angles along the *x*- and
        *y*-axes, respectively, in radians.

        Returns *self*, so this method can easily be chained with more
        calls to :meth:`rotate`, :meth:`rotate_deg`, :meth:`translate`
        and :meth:`scale`.
        """
        rotX = np.tan(xShear)
        rotY = np.tan(yShear)
        skew_mtx = np.array(
            [[1.0, rotX, 0.0], [rotY, 1.0, 0.0], [0.0, 0.0, 1.0]], float)
        self._mtx = np.dot(skew_mtx, self._mtx)
        self.invalidate()
        return self

    def skew_deg(self, xShear, yShear):
        """
        Adds a skew in place.

        *xShear* and *yShear* are the shear angles along the *x*- and
        *y*-axes, respectively, in degrees.

        Returns *self*, so this method can easily be chained with more
        calls to :meth:`rotate`, :meth:`rotate_deg`, :meth:`translate`
        and :meth:`scale`.
        """
        return self.skew(np.deg2rad(xShear), np.deg2rad(yShear))

    def _get_is_separable(self):
        mtx = self.get_matrix()
        return mtx[0, 1] == 0.0 and mtx[1, 0] == 0.0
    is_separable = property(_get_is_separable)


class IdentityTransform(Affine2DBase):
    """
    A special class that does one thing, the identity transform, in a
    fast way.
    """
    _mtx = np.identity(3)

    def frozen(self):
        return self
    frozen.__doc__ = Affine2DBase.frozen.__doc__

    def __str__(self):
        return ("{}()"
                .format(type(self).__name__))

    def get_matrix(self):
        return self._mtx
    get_matrix.__doc__ = Affine2DBase.get_matrix.__doc__

    def transform(self, points):
        return np.asanyarray(points)
    transform.__doc__ = Affine2DBase.transform.__doc__

    transform_affine = transform
    transform_affine.__doc__ = Affine2DBase.transform_affine.__doc__

    transform_non_affine = transform
    transform_non_affine.__doc__ = Affine2DBase.transform_non_affine.__doc__

    def transform_path(self, path):
        return path
    transform_path.__doc__ = Affine2DBase.transform_path.__doc__

    transform_path_affine = transform_path
    transform_path_affine.__doc__ = Affine2DBase.transform_path_affine.__doc__

    transform_path_non_affine = transform_path
    transform_path_non_affine.__doc__ = Affine2DBase.transform_path_non_affine.__doc__

    def get_affine(self):
        return self
    get_affine.__doc__ = Affine2DBase.get_affine.__doc__

    inverted = get_affine
    inverted.__doc__ = Affine2DBase.inverted.__doc__


class BlendedGenericTransform(Transform):
    """
    A "blended" transform uses one transform for the *x*-direction, and
    another transform for the *y*-direction.

    This "generic" version can handle any given child transform in the
    *x*- and *y*-directions.
    """
    input_dims = 2
    output_dims = 2
    is_separable = True
    pass_through = True

    def __init__(self, x_transform, y_transform, **kwargs):
        """
        Create a new "blended" transform using *x_transform* to
        transform the *x*-axis and *y_transform* to transform the
        *y*-axis.

        You will generally not call this constructor directly but use
        the :func:`blended_transform_factory` function instead, which
        can determine automatically which kind of blended transform to
        create.
        """
        # Here we ask: "Does it blend?"

        Transform.__init__(self, **kwargs)
        self._x = x_transform
        self._y = y_transform
        self.set_children(x_transform, y_transform)
        self._affine = None

    def __eq__(self, other):
        # Note, this is an exact copy of BlendedAffine2D.__eq__
        if isinstance(other, (BlendedAffine2D, BlendedGenericTransform)):
            return (self._x == other._x) and (self._y == other._y)
        elif self._x == self._y:
            return self._x == other
        else:
            return NotImplemented

    def contains_branch_seperately(self, transform):
        # Note, this is an exact copy of BlendedAffine2D.contains_branch_seperately
        return self._x.contains_branch(transform), self._y.contains_branch(transform)

    @property
    def depth(self):
        return max(self._x.depth, self._y.depth)

    def contains_branch(self, other):
        # a blended transform cannot possibly contain a branch from two different transforms.
        return False

    def _get_is_affine(self):
        return self._x.is_affine and self._y.is_affine
    is_affine = property(_get_is_affine)

    def _get_has_inverse(self):
        return self._x.has_inverse and self._y.has_inverse
    has_inverse = property(_get_has_inverse)

    def frozen(self):
        return blended_transform_factory(self._x.frozen(), self._y.frozen())
    frozen.__doc__ = Transform.frozen.__doc__

    def __str__(self):
        return ("{}(\n"
                    "{},\n"
                    "{})"
                .format(type(self).__name__,
                        _indent_str(self._x),
                        _indent_str(self._y)))

    def transform_non_affine(self, points):
        if self._x.is_affine and self._y.is_affine:
            return points
        x = self._x
        y = self._y

        if x == y and x.input_dims == 2:
            return x.transform_non_affine(points)

        if x.input_dims == 2:
            x_points = x.transform_non_affine(points)[:, 0:1]
        else:
            x_points = x.transform_non_affine(points[:, 0])
            x_points = x_points.reshape((len(x_points), 1))

        if y.input_dims == 2:
            y_points = y.transform_non_affine(points)[:, 1:]
        else:
            y_points = y.transform_non_affine(points[:, 1])
            y_points = y_points.reshape((len(y_points), 1))

        if (isinstance(x_points, np.ma.MaskedArray) or
                isinstance(y_points, np.ma.MaskedArray)):
            return np.ma.concatenate((x_points, y_points), 1)
        else:
            return np.concatenate((x_points, y_points), 1)
    transform_non_affine.__doc__ = Transform.transform_non_affine.__doc__

    def inverted(self):
        return BlendedGenericTransform(self._x.inverted(), self._y.inverted())
    inverted.__doc__ = Transform.inverted.__doc__

    def get_affine(self):
        if self._invalid or self._affine is None:
            if self._x == self._y:
                self._affine = self._x.get_affine()
            else:
                x_mtx = self._x.get_affine().get_matrix()
                y_mtx = self._y.get_affine().get_matrix()
                # This works because we already know the transforms are
                # separable, though normally one would want to set b and
                # c to zero.
                mtx = np.vstack((x_mtx[0], y_mtx[1], [0.0, 0.0, 1.0]))
                self._affine = Affine2D(mtx)
            self._invalid = 0
        return self._affine
    get_affine.__doc__ = Transform.get_affine.__doc__


class BlendedAffine2D(Affine2DBase):
    """
    A "blended" transform uses one transform for the *x*-direction, and
    another transform for the *y*-direction.

    This version is an optimization for the case where both child
    transforms are of type :class:`Affine2DBase`.
    """
    is_separable = True

    def __init__(self, x_transform, y_transform, **kwargs):
        """
        Create a new "blended" transform using *x_transform* to
        transform the *x*-axis and *y_transform* to transform the
        *y*-axis.

        Both *x_transform* and *y_transform* must be 2D affine
        transforms.

        You will generally not call this constructor directly but use
        the :func:`blended_transform_factory` function instead, which
        can determine automatically which kind of blended transform to
        create.
        """
        is_affine = x_transform.is_affine and y_transform.is_affine
        is_separable = x_transform.is_separable and y_transform.is_separable
        is_correct = is_affine and is_separable
        if not is_correct:
            raise ValueError("Both *x_transform* and *y_transform* must be 2D "
                             "affine transforms")

        Transform.__init__(self, **kwargs)
        self._x = x_transform
        self._y = y_transform
        self.set_children(x_transform, y_transform)

        Affine2DBase.__init__(self)
        self._mtx = None

    def __eq__(self, other):
        # Note, this is an exact copy of BlendedGenericTransform.__eq__
        if isinstance(other, (BlendedAffine2D, BlendedGenericTransform)):
            return (self._x == other._x) and (self._y == other._y)
        elif self._x == self._y:
            return self._x == other
        else:
            return NotImplemented

    def contains_branch_seperately(self, transform):
        # Note, this is an exact copy of BlendedTransform.contains_branch_seperately
        return self._x.contains_branch(transform), self._y.contains_branch(transform)

    def __str__(self):
        return ("{}(\n"
                    "{},\n"
                    "{})"
                .format(type(self).__name__,
                        _indent_str(self._x),
                        _indent_str(self._y)))

    def get_matrix(self):
        if self._invalid:
            if self._x == self._y:
                self._mtx = self._x.get_matrix()
            else:
                x_mtx = self._x.get_matrix()
                y_mtx = self._y.get_matrix()
                # This works because we already know the transforms are
                # separable, though normally one would want to set b and
                # c to zero.
                self._mtx = np.vstack((x_mtx[0], y_mtx[1], [0.0, 0.0, 1.0]))
            self._inverted = None
            self._invalid = 0
        return self._mtx
    get_matrix.__doc__ = Affine2DBase.get_matrix.__doc__


def blended_transform_factory(x_transform, y_transform):
    """
    Create a new "blended" transform using *x_transform* to transform
    the *x*-axis and *y_transform* to transform the *y*-axis.

    A faster version of the blended transform is returned for the case
    where both child transforms are affine.
    """
    if (isinstance(x_transform, Affine2DBase)
        and isinstance(y_transform, Affine2DBase)):
        return BlendedAffine2D(x_transform, y_transform)
    return BlendedGenericTransform(x_transform, y_transform)


class CompositeGenericTransform(Transform):
    """
    A composite transform formed by applying transform *a* then
    transform *b*.

    This "generic" version can handle any two arbitrary
    transformations.
    """
    pass_through = True

    def __init__(self, a, b, **kwargs):
        """
        Create a new composite transform that is the result of
        applying transform *a* then transform *b*.

        You will generally not call this constructor directly but use
        the :func:`composite_transform_factory` function instead,
        which can automatically choose the best kind of composite
        transform instance to create.
        """
        if a.output_dims != b.input_dims:
            raise ValueError("The output dimension of 'a' must be equal to "
                             "the input dimensions of 'b'")
        self.input_dims = a.input_dims
        self.output_dims = b.output_dims

        Transform.__init__(self, **kwargs)
        self._a = a
        self._b = b
        self.set_children(a, b)

    is_affine = property(lambda self: self._a.is_affine and self._b.is_affine)

    def frozen(self):
        self._invalid = 0
        frozen = composite_transform_factory(self._a.frozen(), self._b.frozen())
        if not isinstance(frozen, CompositeGenericTransform):
            return frozen.frozen()
        return frozen
    frozen.__doc__ = Transform.frozen.__doc__

    def _invalidate_internal(self, value, invalidating_node):
        # In some cases for a composite transform, an invalidating call to AFFINE_ONLY needs
        # to be extended to invalidate the NON_AFFINE part too. These cases are when the right
        # hand transform is non-affine and either:
        # (a) the left hand transform is non affine
        # (b) it is the left hand node which has triggered the invalidation
        if value == Transform.INVALID_AFFINE \
            and not self._b.is_affine \
            and (not self._a.is_affine or invalidating_node is self._a):

            value = Transform.INVALID

        Transform._invalidate_internal(self, value=value,
                                       invalidating_node=invalidating_node)

    def __eq__(self, other):
        if isinstance(other, (CompositeGenericTransform, CompositeAffine2D)):
            return self is other or (self._a == other._a and self._b == other._b)
        else:
            return False

    def _iter_break_from_left_to_right(self):
        for lh_compliment, rh_compliment in self._a._iter_break_from_left_to_right():
            yield lh_compliment, rh_compliment + self._b
        for lh_compliment, rh_compliment in self._b._iter_break_from_left_to_right():
            yield self._a + lh_compliment, rh_compliment

    @property
    def depth(self):
        return self._a.depth + self._b.depth

    def _get_is_affine(self):
        return self._a.is_affine and self._b.is_affine
    is_affine = property(_get_is_affine)

    def _get_is_separable(self):
        return self._a.is_separable and self._b.is_separable
    is_separable = property(_get_is_separable)

    def __str__(self):
        return ("{}(\n"
                    "{},\n"
                    "{})"
                .format(type(self).__name__,
                        _indent_str(self._a),
                        _indent_str(self._b)))

    def transform_affine(self, points):
        return self.get_affine().transform(points)
    transform_affine.__doc__ = Transform.transform_affine.__doc__

    def transform_non_affine(self, points):
        if self._a.is_affine and self._b.is_affine:
            return points
        elif not self._a.is_affine and self._b.is_affine:
            return self._a.transform_non_affine(points)
        else:
            return self._b.transform_non_affine(
                                self._a.transform(points))
    transform_non_affine.__doc__ = Transform.transform_non_affine.__doc__

    def transform_path_non_affine(self, path):
        if self._a.is_affine and self._b.is_affine:
            return path
        elif not self._a.is_affine and self._b.is_affine:
            return self._a.transform_path_non_affine(path)
        else:
            return self._b.transform_path_non_affine(
                                    self._a.transform_path(path))
    transform_path_non_affine.__doc__ = Transform.transform_path_non_affine.__doc__

    def get_affine(self):
        if not self._b.is_affine:
            return self._b.get_affine()
        else:
            return Affine2D(np.dot(self._b.get_affine().get_matrix(),
                                self._a.get_affine().get_matrix()))
    get_affine.__doc__ = Transform.get_affine.__doc__

    def inverted(self):
        return CompositeGenericTransform(self._b.inverted(), self._a.inverted())
    inverted.__doc__ = Transform.inverted.__doc__

    def _get_has_inverse(self):
        return self._a.has_inverse and self._b.has_inverse
    has_inverse = property(_get_has_inverse)


class CompositeAffine2D(Affine2DBase):
    """
    A composite transform formed by applying transform *a* then transform *b*.

    This version is an optimization that handles the case where both *a*
    and *b* are 2D affines.
    """
    def __init__(self, a, b, **kwargs):
        """
        Create a new composite transform that is the result of
        applying transform *a* then transform *b*.

        Both *a* and *b* must be instances of :class:`Affine2DBase`.

        You will generally not call this constructor directly but use
        the :func:`composite_transform_factory` function instead,
        which can automatically choose the best kind of composite
        transform instance to create.
        """
        if not a.is_affine or not b.is_affine:
            raise ValueError("'a' and 'b' must be affine transforms")
        if a.output_dims != b.input_dims:
            raise ValueError("The output dimension of 'a' must be equal to "
                             "the input dimensions of 'b'")
        self.input_dims = a.input_dims
        self.output_dims = b.output_dims

        Affine2DBase.__init__(self, **kwargs)
        self._a = a
        self._b = b
        self.set_children(a, b)
        self._mtx = None

    @property
    def depth(self):
        return self._a.depth + self._b.depth

    def _iter_break_from_left_to_right(self):
        for lh_compliment, rh_compliment in self._a._iter_break_from_left_to_right():
            yield lh_compliment, rh_compliment + self._b
        for lh_compliment, rh_compliment in self._b._iter_break_from_left_to_right():
            yield self._a + lh_compliment, rh_compliment

    def __str__(self):
        return ("{}(\n"
                    "{},\n"
                    "{})"
                .format(type(self).__name__,
                        _indent_str(self._a),
                        _indent_str(self._b)))

    def get_matrix(self):
        if self._invalid:
            self._mtx = np.dot(
                self._b.get_matrix(),
                self._a.get_matrix())
            self._inverted = None
            self._invalid = 0
        return self._mtx
    get_matrix.__doc__ = Affine2DBase.get_matrix.__doc__


def composite_transform_factory(a, b):
    """
    Create a new composite transform that is the result of applying
    transform a then transform b.

    Shortcut versions of the blended transform are provided for the
    case where both child transforms are affine, or one or the other
    is the identity transform.

    Composite transforms may also be created using the '+' operator,
    e.g.::

      c = a + b
    """
    # check to see if any of a or b are IdentityTransforms. We use
    # isinstance here to guarantee that the transforms will *always*
    # be IdentityTransforms. Since TransformWrappers are mutable,
    # use of equality here would be wrong.
    if isinstance(a, IdentityTransform):
        return b
    elif isinstance(b, IdentityTransform):
        return a
    elif isinstance(a, Affine2D) and isinstance(b, Affine2D):
        return CompositeAffine2D(a, b)
    return CompositeGenericTransform(a, b)


class BboxTransform(Affine2DBase):
    """
    :class:`BboxTransform` linearly transforms points from one
    :class:`Bbox` to another :class:`Bbox`.
    """
    is_separable = True

    def __init__(self, boxin, boxout, **kwargs):
        """
        Create a new :class:`BboxTransform` that linearly transforms
        points from *boxin* to *boxout*.
        """
        if not boxin.is_bbox or not boxout.is_bbox:
            raise ValueError("'boxin' and 'boxout' must be bbox")

        Affine2DBase.__init__(self, **kwargs)
        self._boxin = boxin
        self._boxout = boxout
        self.set_children(boxin, boxout)
        self._mtx = None
        self._inverted = None

    def __str__(self):
        return ("{}(\n"
                    "{},\n"
                    "{})"
                .format(type(self).__name__,
                        _indent_str(self._boxin),
                        _indent_str(self._boxout)))

    def get_matrix(self):
        if self._invalid:
            inl, inb, inw, inh = self._boxin.bounds
            outl, outb, outw, outh = self._boxout.bounds
            x_scale = outw / inw
            y_scale = outh / inh
            if DEBUG and (x_scale == 0 or y_scale == 0):
                raise ValueError("Transforming from or to a singular bounding box.")
            self._mtx = np.array([[x_scale, 0.0    , (-inl*x_scale+outl)],
                                  [0.0    , y_scale, (-inb*y_scale+outb)],
                                  [0.0    , 0.0    , 1.0        ]],
                                 float)
            self._inverted = None
            self._invalid = 0
        return self._mtx
    get_matrix.__doc__ = Affine2DBase.get_matrix.__doc__


class BboxTransformTo(Affine2DBase):
    """
    :class:`BboxTransformTo` is a transformation that linearly
    transforms points from the unit bounding box to a given
    :class:`Bbox`.
    """
    is_separable = True

    def __init__(self, boxout, **kwargs):
        """
        Create a new :class:`BboxTransformTo` that linearly transforms
        points from the unit bounding box to *boxout*.
        """
        if not boxout.is_bbox:
            raise ValueError("'boxout' must be bbox")

        Affine2DBase.__init__(self, **kwargs)
        self._boxout = boxout
        self.set_children(boxout)
        self._mtx = None
        self._inverted = None

    def __str__(self):
        return ("{}(\n"
                    "{})"
                .format(type(self).__name__,
                        _indent_str(self._boxout)))

    def get_matrix(self):
        if self._invalid:
            outl, outb, outw, outh = self._boxout.bounds
            if DEBUG and (outw == 0 or outh == 0):
                raise ValueError("Transforming to a singular bounding box.")
            self._mtx = np.array([[outw,  0.0, outl],
                                  [ 0.0, outh, outb],
                                  [ 0.0,  0.0,  1.0]],
                                  float)
            self._inverted = None
            self._invalid = 0
        return self._mtx
    get_matrix.__doc__ = Affine2DBase.get_matrix.__doc__


class BboxTransformToMaxOnly(BboxTransformTo):
    """
    :class:`BboxTransformTo` is a transformation that linearly
    transforms points from the unit bounding box to a given
    :class:`Bbox` with a fixed upper left of (0, 0).
    """
    def get_matrix(self):
        if self._invalid:
            xmax, ymax = self._boxout.max
            if DEBUG and (xmax == 0 or ymax == 0):
                raise ValueError("Transforming to a singular bounding box.")
            self._mtx = np.array([[xmax,  0.0, 0.0],
                                  [ 0.0, ymax, 0.0],
                                  [ 0.0,  0.0, 1.0]],
                                 float)
            self._inverted = None
            self._invalid = 0
        return self._mtx
    get_matrix.__doc__ = Affine2DBase.get_matrix.__doc__


class BboxTransformFrom(Affine2DBase):
    """
    :class:`BboxTransformFrom` linearly transforms points from a given
    :class:`Bbox` to the unit bounding box.
    """
    is_separable = True

    def __init__(self, boxin, **kwargs):
        if not boxin.is_bbox:
            raise ValueError("'boxin' must be bbox")

        Affine2DBase.__init__(self, **kwargs)
        self._boxin = boxin
        self.set_children(boxin)
        self._mtx = None
        self._inverted = None

    def __str__(self):
        return ("{}(\n"
                    "{})"
                .format(type(self).__name__,
                        _indent_str(self._boxin)))

    def get_matrix(self):
        if self._invalid:
            inl, inb, inw, inh = self._boxin.bounds
            if DEBUG and (inw == 0 or inh == 0):
                raise ValueError("Transforming from a singular bounding box.")
            x_scale = 1.0 / inw
            y_scale = 1.0 / inh
            self._mtx = np.array([[x_scale, 0.0    , (-inl*x_scale)],
                                  [0.0    , y_scale, (-inb*y_scale)],
                                  [0.0    , 0.0    , 1.0        ]],
                                 float)
            self._inverted = None
            self._invalid = 0
        return self._mtx
    get_matrix.__doc__ = Affine2DBase.get_matrix.__doc__


class ScaledTranslation(Affine2DBase):
    """
    A transformation that translates by *xt* and *yt*, after *xt* and *yt*
    have been transformad by the given transform *scale_trans*.
    """
    def __init__(self, xt, yt, scale_trans, **kwargs):
        Affine2DBase.__init__(self, **kwargs)
        self._t = (xt, yt)
        self._scale_trans = scale_trans
        self.set_children(scale_trans)
        self._mtx = None
        self._inverted = None

    def __str__(self):
        return ("{}(\n"
                    "{})"
                .format(type(self).__name__,
                        _indent_str(self._t)))

    def get_matrix(self):
        if self._invalid:
            xt, yt = self._scale_trans.transform_point(self._t)
            self._mtx = np.array([[1.0, 0.0, xt],
                                  [0.0, 1.0, yt],
                                  [0.0, 0.0, 1.0]],
                                 float)
            self._invalid = 0
            self._inverted = None
        return self._mtx
    get_matrix.__doc__ = Affine2DBase.get_matrix.__doc__


class TransformedPath(TransformNode):
    """
    A :class:`TransformedPath` caches a non-affine transformed copy of
    the :class:`~matplotlib.path.Path`.  This cached copy is
    automatically updated when the non-affine part of the transform
    changes.

    .. note::

        Paths are considered immutable by this class. Any update to the
        path's vertices/codes will not trigger a transform recomputation.

    """
    def __init__(self, path, transform):
        """
        Create a new :class:`TransformedPath` from the given
        :class:`~matplotlib.path.Path` and :class:`Transform`.
        """
        if not isinstance(transform, Transform):
            raise ValueError("'transform' must be an instance of "
                             "'matplotlib.transform.Transform'")
        TransformNode.__init__(self)

        self._path = path
        self._transform = transform
        self.set_children(transform)
        self._transformed_path = None
        self._transformed_points = None

    def _revalidate(self):
        # only recompute if the invalidation includes the non_affine part of the transform
        if ((self._invalid & self.INVALID_NON_AFFINE == self.INVALID_NON_AFFINE)
            or self._transformed_path is None):
            self._transformed_path = \
                self._transform.transform_path_non_affine(self._path)
            self._transformed_points = \
                Path._fast_from_codes_and_verts(
                    self._transform.transform_non_affine(self._path.vertices),
                    None,
                    {'interpolation_steps': self._path._interpolation_steps,
                     'should_simplify': self._path.should_simplify})
        self._invalid = 0

    def get_transformed_points_and_affine(self):
        """
        Return a copy of the child path, with the non-affine part of
        the transform already applied, along with the affine part of
        the path necessary to complete the transformation.  Unlike
        :meth:`get_transformed_path_and_affine`, no interpolation will
        be performed.
        """
        self._revalidate()
        return self._transformed_points, self.get_affine()

    def get_transformed_path_and_affine(self):
        """
        Return a copy of the child path, with the non-affine part of
        the transform already applied, along with the affine part of
        the path necessary to complete the transformation.
        """
        self._revalidate()
        return self._transformed_path, self.get_affine()

    def get_fully_transformed_path(self):
        """
        Return a fully-transformed copy of the child path.
        """
        self._revalidate()
        return self._transform.transform_path_affine(self._transformed_path)

    def get_affine(self):
        return self._transform.get_affine()


class TransformedPatchPath(TransformedPath):
    """
    A :class:`TransformedPatchPath` caches a non-affine transformed copy of
    the :class:`~matplotlib.path.Patch`. This cached copy is automatically
    updated when the non-affine part of the transform or the patch changes.
    """
    def __init__(self, patch):
        """
        Create a new :class:`TransformedPatchPath` from the given
        :class:`~matplotlib.path.Patch`.
        """
        TransformNode.__init__(self)

        transform = patch.get_transform()
        self._patch = patch
        self._transform = transform
        self.set_children(transform)
        self._path = patch.get_path()
        self._transformed_path = None
        self._transformed_points = None

    def _revalidate(self):
        patch_path = self._patch.get_path()
        # Only recompute if the invalidation includes the non_affine part of
        # the transform, or the Patch's Path has changed.
        if (self._transformed_path is None or self._path != patch_path or
                (self._invalid & self.INVALID_NON_AFFINE ==
                    self.INVALID_NON_AFFINE)):
            self._path = patch_path
            self._transformed_path = \
                self._transform.transform_path_non_affine(patch_path)
            self._transformed_points = \
                Path._fast_from_codes_and_verts(
                    self._transform.transform_non_affine(patch_path.vertices),
                    None,
                    {'interpolation_steps': patch_path._interpolation_steps,
                     'should_simplify': patch_path.should_simplify})
        self._invalid = 0


def nonsingular(vmin, vmax, expander=0.001, tiny=1e-15, increasing=True):
    """
    Modify the endpoints of a range as needed to avoid singularities.

    Parameters
    ----------
    vmin, vmax : float
        The initial endpoints.
    expander : float, optional, default: 0.001
        Fractional amount by which *vmin* and *vmax* are expanded if
        the original interval is too small, based on *tiny*.
    tiny : float, optional, default: 1e-15
        Threshold for the ratio of the interval to the maximum absolute
        value of its endpoints.  If the interval is smaller than
        this, it will be expanded.  This value should be around
        1e-15 or larger; otherwise the interval will be approaching
        the double precision resolution limit.
    increasing : bool, optional, default: True
        If True, swap *vmin*, *vmax* if *vmin* > *vmax*.

    Returns
    -------
    vmin, vmax : float
        Endpoints, expanded and/or swapped if necessary.
        If either input is inf or NaN, or if both inputs are 0 or very
        close to zero, it returns -*expander*, *expander*.
    """

    if (not np.isfinite(vmin)) or (not np.isfinite(vmax)):
        return -expander, expander

    swapped = False
    if vmax < vmin:
        vmin, vmax = vmax, vmin
        swapped = True

    maxabsvalue = max(abs(vmin), abs(vmax))
    if maxabsvalue < (1e6 / tiny) * np.finfo(float).tiny:
        vmin = -expander
        vmax = expander

    elif vmax - vmin <= maxabsvalue * tiny:
        if vmax == 0 and vmin == 0:
            vmin = -expander
            vmax = expander
        else:
            vmin -= expander*abs(vmin)
            vmax += expander*abs(vmax)

    if swapped and not increasing:
        vmin, vmax = vmax, vmin
    return vmin, vmax


def interval_contains(interval, val):
    """
    Check, inclusively, whether an interval includes a given value.

    Parameters
    ----------
    interval : sequence of scalar
        A 2-length sequence, endpoints that define the interval.
    val : scalar
        Value to check is within interval.

    Returns
    -------
    bool
        Returns true if given val is within the interval.
    """
    a, b = interval
    return a <= val <= b or a >= val >= b


def interval_contains_open(interval, val):
    """
    Check, excluding endpoints, whether an interval includes a given value.

    Parameters
    ----------
    interval : sequence of scalar
        A 2-length sequence, endpoints that define the interval.
    val : scalar
        Value to check is within interval.

    Returns
    -------
    bool
        Returns true if given val is within the interval.
    """
    a, b = interval
    return a < val < b or a > val > b


def offset_copy(trans, fig=None, x=0.0, y=0.0, units='inches'):
    """
    Return a new transform with an added offset.

    Parameters
    ----------
    trans : :class:`Transform` instance
        Any transform, to which offset will be applied.
    fig : :class:`~matplotlib.figure.Figure`, optional, default: None
        Current figure. It can be None if *units* are 'dots'.
    x, y : float, optional, default: 0.0
        Specifies the offset to apply.
    units : {'inches', 'points', 'dots'}, optional
        Units of the offset.

    Returns
    -------
    trans : :class:`Transform` instance
        Transform with applied offset.
    """
    if units == 'dots':
        return trans + Affine2D().translate(x, y)
    if fig is None:
        raise ValueError('For units of inches or points a fig kwarg is needed')
    if units == 'points':
        x /= 72.0
        y /= 72.0
    elif not units == 'inches':
        raise ValueError('units must be dots, points, or inches')
    return trans + ScaledTranslation(x, y, fig.dpi_scale_trans)
