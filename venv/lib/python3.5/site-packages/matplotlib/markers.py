"""
This module contains functions to handle markers.  Used by both the
marker functionality of `~matplotlib.axes.Axes.plot` and
`~matplotlib.axes.Axes.scatter`.

All possible markers are defined here:

============================== ===============================================
marker                         description
============================== ===============================================
`"."`                          point
`","`                          pixel
`"o"`                          circle
`"v"`                          triangle_down
`"^"`                          triangle_up
`"<"`                          triangle_left
`">"`                          triangle_right
`"1"`                          tri_down
`"2"`                          tri_up
`"3"`                          tri_left
`"4"`                          tri_right
`"8"`                          octagon
`"s"`                          square
`"p"`                          pentagon
`"P"`                          plus (filled)
`"*"`                          star
`"h"`                          hexagon1
`"H"`                          hexagon2
`"+"`                          plus
`"x"`                          x
`"X"`                          x (filled)
`"D"`                          diamond
`"d"`                          thin_diamond
`"|"`                          vline
`"_"`                          hline
TICKLEFT                       tickleft
TICKRIGHT                      tickright
TICKUP                         tickup
TICKDOWN                       tickdown
CARETLEFT                      caretleft (centered at tip)
CARETRIGHT                     caretright (centered at tip)
CARETUP                        caretup (centered at tip)
CARETDOWN                      caretdown (centered at tip)
CARETLEFTBASE                  caretleft (centered at base)
CARETRIGHTBASE                 caretright (centered at base)
CARETUPBASE                    caretup (centered at base)
`"None"`, `" "` or `""`        nothing
``'$...$'``                    render the string using mathtext.
`verts`                        a list of (x, y) pairs used for Path vertices.
                               The center of the marker is located at (0,0) and
                               the size is normalized.
path                           a `~matplotlib.path.Path` instance.
(`numsides`, `style`, `angle`) The marker can also be a tuple (`numsides`,
                               `style`, `angle`), which will create a custom,
                               regular symbol.

                               `numsides`:
                                   the number of sides

                               `style`:
                                   the style of the regular symbol:

                                   0
                                     a regular polygon
                                   1
                                     a star-like symbol
                                   2
                                     an asterisk
                                   3
                                     a circle (`numsides` and `angle` is
                                     ignored)

                               `angle`:
                                   the angle of rotation of the symbol
============================== ===============================================

For backward compatibility, the form (`verts`, 0) is also accepted,
but it is equivalent to just `verts` for giving a raw set of vertices
that define the shape.

`None` is the default which means 'nothing', however this table is
referred to from other docs for the valid inputs from marker inputs and in
those cases `None` still means 'default'.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves import xrange

from collections import Sized

import numpy as np

from . import rcParams
from .cbook import is_math_text, is_numlike
from .path import Path
from .transforms import IdentityTransform, Affine2D

# special-purpose marker identifiers:
(TICKLEFT, TICKRIGHT, TICKUP, TICKDOWN,
 CARETLEFT, CARETRIGHT, CARETUP, CARETDOWN,
 CARETLEFTBASE, CARETRIGHTBASE, CARETUPBASE, CARETDOWNBASE) = xrange(12)

_empty_path = Path(np.empty((0, 2)))


class MarkerStyle(object):

    markers = {
        '.': 'point',
        ',': 'pixel',
        'o': 'circle',
        'v': 'triangle_down',
        '^': 'triangle_up',
        '<': 'triangle_left',
        '>': 'triangle_right',
        '1': 'tri_down',
        '2': 'tri_up',
        '3': 'tri_left',
        '4': 'tri_right',
        '8': 'octagon',
        's': 'square',
        'p': 'pentagon',
        '*': 'star',
        'h': 'hexagon1',
        'H': 'hexagon2',
        '+': 'plus',
        'x': 'x',
        'D': 'diamond',
        'd': 'thin_diamond',
        '|': 'vline',
        '_': 'hline',
        'P': 'plus_filled',
        'X': 'x_filled',
        TICKLEFT: 'tickleft',
        TICKRIGHT: 'tickright',
        TICKUP: 'tickup',
        TICKDOWN: 'tickdown',
        CARETLEFT: 'caretleft',
        CARETRIGHT: 'caretright',
        CARETUP: 'caretup',
        CARETDOWN: 'caretdown',
        CARETLEFTBASE: 'caretleftbase',
        CARETRIGHTBASE: 'caretrightbase',
        CARETUPBASE: 'caretupbase',
        CARETDOWNBASE: 'caretdownbase',
        "None": 'nothing',
        None: 'nothing',
        ' ': 'nothing',
        '': 'nothing'
    }

    # Just used for informational purposes.  is_filled()
    # is calculated in the _set_* functions.
    filled_markers = (
        'o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd',
        'P', 'X')

    fillstyles = ('full', 'left', 'right', 'bottom', 'top', 'none')
    _half_fillstyles = ('left', 'right', 'bottom', 'top')

    # TODO: Is this ever used as a non-constant?
    _point_size_reduction = 0.5

    def __init__(self, marker=None, fillstyle=None):
        """
        MarkerStyle

        Attributes
        ----------
        markers : list of known marks

        fillstyles : list of known fillstyles

        filled_markers : list of known filled markers.

        Parameters
        ----------
        marker : string or array_like, optional, default: None
            See the descriptions of possible markers in the module docstring.

        fillstyle : string, optional, default: 'full'
            'full', 'left", 'right', 'bottom', 'top', 'none'
        """
        self._marker_function = None
        self.set_fillstyle(fillstyle)
        self.set_marker(marker)

    def __getstate__(self):
        d = self.__dict__.copy()
        d.pop('_marker_function')
        return d

    def __setstate__(self, statedict):
        self.__dict__ = statedict
        self.set_marker(self._marker)

    def _recache(self):
        if self._marker_function is None:
            return
        self._path = _empty_path
        self._transform = IdentityTransform()
        self._alt_path = None
        self._alt_transform = None
        self._snap_threshold = None
        self._joinstyle = 'round'
        self._capstyle = 'butt'
        self._filled = True
        self._marker_function()

    if six.PY3:
        def __bool__(self):
            return bool(len(self._path.vertices))
    else:
        def __nonzero__(self):
            return bool(len(self._path.vertices))

    def is_filled(self):
        return self._filled

    def get_fillstyle(self):
        return self._fillstyle

    def set_fillstyle(self, fillstyle):
        """
        Sets fillstyle

        Parameters
        ----------
        fillstyle : string amongst known fillstyles
        """
        if fillstyle is None:
            fillstyle = rcParams['markers.fillstyle']
        if fillstyle not in self.fillstyles:
            raise ValueError("Unrecognized fillstyle %s"
                             % ' '.join(self.fillstyles))
        self._fillstyle = fillstyle
        self._recache()

    def get_joinstyle(self):
        return self._joinstyle

    def get_capstyle(self):
        return self._capstyle

    def get_marker(self):
        return self._marker

    def set_marker(self, marker):
        if (isinstance(marker, np.ndarray) and marker.ndim == 2 and
                marker.shape[1] == 2):
            self._marker_function = self._set_vertices
        elif (isinstance(marker, Sized) and len(marker) in (2, 3) and
                marker[1] in (0, 1, 2, 3)):
            self._marker_function = self._set_tuple_marker
        elif (not isinstance(marker, (np.ndarray, list)) and
              marker in self.markers):
            self._marker_function = getattr(
                self, '_set_' + self.markers[marker])
        elif isinstance(marker, six.string_types) and is_math_text(marker):
            self._marker_function = self._set_mathtext_path
        elif isinstance(marker, Path):
            self._marker_function = self._set_path_marker
        else:
            try:
                Path(marker)
                self._marker_function = self._set_vertices
            except ValueError:
                raise ValueError('Unrecognized marker style'
                                 ' {0}'.format(marker))

        self._marker = marker
        self._recache()

    def get_path(self):
        return self._path

    def get_transform(self):
        return self._transform.frozen()

    def get_alt_path(self):
        return self._alt_path

    def get_alt_transform(self):
        return self._alt_transform.frozen()

    def get_snap_threshold(self):
        return self._snap_threshold

    def _set_nothing(self):
        self._filled = False

    def _set_custom_marker(self, path):
        verts = path.vertices
        rescale = max(np.max(np.abs(verts[:, 0])),
                      np.max(np.abs(verts[:, 1])))
        self._transform = Affine2D().scale(0.5 / rescale)
        self._path = path

    def _set_path_marker(self):
        self._set_custom_marker(self._marker)

    def _set_vertices(self):
        verts = self._marker
        marker = Path(verts)
        self._set_custom_marker(marker)

    def _set_tuple_marker(self):
        marker = self._marker
        if is_numlike(marker[0]):
            if len(marker) == 2:
                numsides, rotation = marker[0], 0.0
            elif len(marker) == 3:
                numsides, rotation = marker[0], marker[2]
            symstyle = marker[1]
            if symstyle == 0:
                self._path = Path.unit_regular_polygon(numsides)
                self._joinstyle = 'miter'
            elif symstyle == 1:
                self._path = Path.unit_regular_star(numsides)
                self._joinstyle = 'bevel'
            elif symstyle == 2:
                self._path = Path.unit_regular_asterisk(numsides)
                self._filled = False
                self._joinstyle = 'bevel'
            elif symstyle == 3:
                self._path = Path.unit_circle()
            self._transform = Affine2D().scale(0.5).rotate_deg(rotation)
        else:
            verts = np.asarray(marker[0])
            path = Path(verts)
            self._set_custom_marker(path)

    def _set_mathtext_path(self):
        """
        Draws mathtext markers '$...$' using TextPath object.

        Submitted by tcb
        """
        from matplotlib.text import TextPath
        from matplotlib.font_manager import FontProperties

        # again, the properties could be initialised just once outside
        # this function
        # Font size is irrelevant here, it will be rescaled based on
        # the drawn size later
        props = FontProperties(size=1.0)
        text = TextPath(xy=(0, 0), s=self.get_marker(), fontproperties=props,
                        usetex=rcParams['text.usetex'])
        if len(text.vertices) == 0:
            return

        xmin, ymin = text.vertices.min(axis=0)
        xmax, ymax = text.vertices.max(axis=0)
        width = xmax - xmin
        height = ymax - ymin
        max_dim = max(width, height)
        self._transform = Affine2D() \
            .translate(-xmin + 0.5 * -width, -ymin + 0.5 * -height) \
            .scale(1.0 / max_dim)
        self._path = text
        self._snap = False

    def _half_fill(self):
        fs = self.get_fillstyle()
        result = fs in self._half_fillstyles
        return result

    def _set_circle(self, reduction=1.0):
        self._transform = Affine2D().scale(0.5 * reduction)
        self._snap_threshold = np.inf
        fs = self.get_fillstyle()
        if not self._half_fill():
            self._path = Path.unit_circle()
        else:
            # build a right-half circle
            if fs == 'bottom':
                rotate = 270.
            elif fs == 'top':
                rotate = 90.
            elif fs == 'left':
                rotate = 180.
            else:
                rotate = 0.

            self._path = self._alt_path = Path.unit_circle_righthalf()
            self._transform.rotate_deg(rotate)
            self._alt_transform = self._transform.frozen().rotate_deg(180.)

    def _set_pixel(self):
        self._path = Path.unit_rectangle()
        # Ideally, you'd want -0.5, -0.5 here, but then the snapping
        # algorithm in the Agg backend will round this to a 2x2
        # rectangle from (-1, -1) to (1, 1).  By offsetting it
        # slightly, we can force it to be (0, 0) to (1, 1), which both
        # makes it only be a single pixel and places it correctly
        # aligned to 1-width stroking (i.e. the ticks).  This hack is
        # the best of a number of bad alternatives, mainly because the
        # backends are not aware of what marker is actually being used
        # beyond just its path data.
        self._transform = Affine2D().translate(-0.49999, -0.49999)
        self._snap_threshold = None

    def _set_point(self):
        self._set_circle(reduction=self._point_size_reduction)

    _triangle_path = Path(
        [[0.0, 1.0], [-1.0, -1.0], [1.0, -1.0], [0.0, 1.0]],
        [Path.MOVETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY])
    # Going down halfway looks to small.  Golden ratio is too far.
    _triangle_path_u = Path(
        [[0.0, 1.0], [-3 / 5., -1 / 5.], [3 / 5., -1 / 5.], [0.0, 1.0]],
        [Path.MOVETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY])
    _triangle_path_d = Path(
        [[-3 / 5., -1 / 5.], [3 / 5., -1 / 5.], [1.0, -1.0], [-1.0, -1.0],
         [-3 / 5., -1 / 5.]],
        [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY])
    _triangle_path_l = Path(
        [[0.0, 1.0], [0.0, -1.0], [-1.0, -1.0], [0.0, 1.0]],
        [Path.MOVETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY])
    _triangle_path_r = Path(
        [[0.0, 1.0], [0.0, -1.0], [1.0, -1.0], [0.0, 1.0]],
        [Path.MOVETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY])

    def _set_triangle(self, rot, skip):
        self._transform = Affine2D().scale(0.5, 0.5).rotate_deg(rot)
        self._snap_threshold = 5.0
        fs = self.get_fillstyle()

        if not self._half_fill():
            self._path = self._triangle_path
        else:
            mpaths = [self._triangle_path_u,
                      self._triangle_path_l,
                      self._triangle_path_d,
                      self._triangle_path_r]

            if fs == 'top':
                self._path = mpaths[(0 + skip) % 4]
                self._alt_path = mpaths[(2 + skip) % 4]
            elif fs == 'bottom':
                self._path = mpaths[(2 + skip) % 4]
                self._alt_path = mpaths[(0 + skip) % 4]
            elif fs == 'left':
                self._path = mpaths[(1 + skip) % 4]
                self._alt_path = mpaths[(3 + skip) % 4]
            else:
                self._path = mpaths[(3 + skip) % 4]
                self._alt_path = mpaths[(1 + skip) % 4]

            self._alt_transform = self._transform

        self._joinstyle = 'miter'

    def _set_triangle_up(self):
        return self._set_triangle(0.0, 0)

    def _set_triangle_down(self):
        return self._set_triangle(180.0, 2)

    def _set_triangle_left(self):
        return self._set_triangle(90.0, 3)

    def _set_triangle_right(self):
        return self._set_triangle(270.0, 1)

    def _set_square(self):
        self._transform = Affine2D().translate(-0.5, -0.5)
        self._snap_threshold = 2.0
        fs = self.get_fillstyle()
        if not self._half_fill():
            self._path = Path.unit_rectangle()
        else:
            # build a bottom filled square out of two rectangles, one
            # filled.  Use the rotation to support left, right, bottom
            # or top
            if fs == 'bottom':
                rotate = 0.
            elif fs == 'top':
                rotate = 180.
            elif fs == 'left':
                rotate = 270.
            else:
                rotate = 90.

            self._path = Path([[0.0, 0.0], [1.0, 0.0], [1.0, 0.5],
                               [0.0, 0.5], [0.0, 0.0]])
            self._alt_path = Path([[0.0, 0.5], [1.0, 0.5], [1.0, 1.0],
                                   [0.0, 1.0], [0.0, 0.5]])
            self._transform.rotate_deg(rotate)
            self._alt_transform = self._transform

        self._joinstyle = 'miter'

    def _set_diamond(self):
        self._transform = Affine2D().translate(-0.5, -0.5).rotate_deg(45)
        self._snap_threshold = 5.0
        fs = self.get_fillstyle()
        if not self._half_fill():
            self._path = Path.unit_rectangle()
        else:
            self._path = Path([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 0.0]])
            self._alt_path = Path([[0.0, 0.0], [0.0, 1.0],
                                   [1.0, 1.0], [0.0, 0.0]])

            if fs == 'bottom':
                rotate = 270.
            elif fs == 'top':
                rotate = 90.
            elif fs == 'left':
                rotate = 180.
            else:
                rotate = 0.

            self._transform.rotate_deg(rotate)
            self._alt_transform = self._transform

        self._joinstyle = 'miter'

    def _set_thin_diamond(self):
        self._set_diamond()
        self._transform.scale(0.6, 1.0)

    def _set_pentagon(self):
        self._transform = Affine2D().scale(0.5)
        self._snap_threshold = 5.0

        polypath = Path.unit_regular_polygon(5)
        fs = self.get_fillstyle()

        if not self._half_fill():
            self._path = polypath
        else:
            verts = polypath.vertices

            y = (1 + np.sqrt(5)) / 4.
            top = Path([verts[0], verts[1], verts[4], verts[0]])
            bottom = Path([verts[1], verts[2], verts[3], verts[4], verts[1]])
            left = Path([verts[0], verts[1], verts[2], [0, -y], verts[0]])
            right = Path([verts[0], verts[4], verts[3], [0, -y], verts[0]])

            if fs == 'top':
                mpath, mpath_alt = top, bottom
            elif fs == 'bottom':
                mpath, mpath_alt = bottom, top
            elif fs == 'left':
                mpath, mpath_alt = left, right
            else:
                mpath, mpath_alt = right, left
            self._path = mpath
            self._alt_path = mpath_alt
            self._alt_transform = self._transform

        self._joinstyle = 'miter'

    def _set_star(self):
        self._transform = Affine2D().scale(0.5)
        self._snap_threshold = 5.0

        fs = self.get_fillstyle()
        polypath = Path.unit_regular_star(5, innerCircle=0.381966)

        if not self._half_fill():
            self._path = polypath
        else:
            verts = polypath.vertices

            top = Path(np.vstack((verts[0:4, :], verts[7:10, :], verts[0])))
            bottom = Path(np.vstack((verts[3:8, :], verts[3])))
            left = Path(np.vstack((verts[0:6, :], verts[0])))
            right = Path(np.vstack((verts[0], verts[5:10, :], verts[0])))

            if fs == 'top':
                mpath, mpath_alt = top, bottom
            elif fs == 'bottom':
                mpath, mpath_alt = bottom, top
            elif fs == 'left':
                mpath, mpath_alt = left, right
            else:
                mpath, mpath_alt = right, left
            self._path = mpath
            self._alt_path = mpath_alt
            self._alt_transform = self._transform

        self._joinstyle = 'bevel'

    def _set_hexagon1(self):
        self._transform = Affine2D().scale(0.5)
        self._snap_threshold = None

        fs = self.get_fillstyle()
        polypath = Path.unit_regular_polygon(6)

        if not self._half_fill():
            self._path = polypath
        else:
            verts = polypath.vertices

            # not drawing inside lines
            x = np.abs(np.cos(5 * np.pi / 6.))
            top = Path(np.vstack(([-x, 0], verts[(1, 0, 5), :], [x, 0])))
            bottom = Path(np.vstack(([-x, 0], verts[2:5, :], [x, 0])))
            left = Path(verts[(0, 1, 2, 3), :])
            right = Path(verts[(0, 5, 4, 3), :])

            if fs == 'top':
                mpath, mpath_alt = top, bottom
            elif fs == 'bottom':
                mpath, mpath_alt = bottom, top
            elif fs == 'left':
                mpath, mpath_alt = left, right
            else:
                mpath, mpath_alt = right, left

            self._path = mpath
            self._alt_path = mpath_alt
            self._alt_transform = self._transform

        self._joinstyle = 'miter'

    def _set_hexagon2(self):
        self._transform = Affine2D().scale(0.5).rotate_deg(30)
        self._snap_threshold = None

        fs = self.get_fillstyle()
        polypath = Path.unit_regular_polygon(6)

        if not self._half_fill():
            self._path = polypath
        else:
            verts = polypath.vertices

            # not drawing inside lines
            x, y = np.sqrt(3) / 4, 3 / 4.
            top = Path(verts[(1, 0, 5, 4, 1), :])
            bottom = Path(verts[(1, 2, 3, 4), :])
            left = Path(np.vstack(([x, y], verts[(0, 1, 2), :],
                                   [-x, -y], [x, y])))
            right = Path(np.vstack(([x, y], verts[(5, 4, 3), :], [-x, -y])))

            if fs == 'top':
                mpath, mpath_alt = top, bottom
            elif fs == 'bottom':
                mpath, mpath_alt = bottom, top
            elif fs == 'left':
                mpath, mpath_alt = left, right
            else:
                mpath, mpath_alt = right, left

            self._path = mpath
            self._alt_path = mpath_alt
            self._alt_transform = self._transform

        self._joinstyle = 'miter'

    def _set_octagon(self):
        self._transform = Affine2D().scale(0.5)
        self._snap_threshold = 5.0

        fs = self.get_fillstyle()
        polypath = Path.unit_regular_polygon(8)

        if not self._half_fill():
            self._transform.rotate_deg(22.5)
            self._path = polypath
        else:
            x = np.sqrt(2.) / 4.
            half = Path([[0, -1], [0, 1], [-x, 1], [-1, x],
                         [-1, -x], [-x, -1], [0, -1]])

            if fs == 'bottom':
                rotate = 90.
            elif fs == 'top':
                rotate = 270.
            elif fs == 'right':
                rotate = 180.
            else:
                rotate = 0.

            self._transform.rotate_deg(rotate)
            self._path = self._alt_path = half
            self._alt_transform = self._transform.frozen().rotate_deg(180.0)

        self._joinstyle = 'miter'

    _line_marker_path = Path([[0.0, -1.0], [0.0, 1.0]])

    def _set_vline(self):
        self._transform = Affine2D().scale(0.5)
        self._snap_threshold = 1.0
        self._filled = False
        self._path = self._line_marker_path

    def _set_hline(self):
        self._set_vline()
        self._transform = self._transform.rotate_deg(90)

    _tickhoriz_path = Path([[0.0, 0.0], [1.0, 0.0]])

    def _set_tickleft(self):
        self._transform = Affine2D().scale(-1.0, 1.0)
        self._snap_threshold = 1.0
        self._filled = False
        self._path = self._tickhoriz_path

    def _set_tickright(self):
        self._transform = Affine2D().scale(1.0, 1.0)
        self._snap_threshold = 1.0
        self._filled = False
        self._path = self._tickhoriz_path

    _tickvert_path = Path([[-0.0, 0.0], [-0.0, 1.0]])

    def _set_tickup(self):
        self._transform = Affine2D().scale(1.0, 1.0)
        self._snap_threshold = 1.0
        self._filled = False
        self._path = self._tickvert_path

    def _set_tickdown(self):
        self._transform = Affine2D().scale(1.0, -1.0)
        self._snap_threshold = 1.0
        self._filled = False
        self._path = self._tickvert_path

    _tri_path = Path([[0.0, 0.0], [0.0, -1.0],
                      [0.0, 0.0], [0.8, 0.5],
                      [0.0, 0.0], [-0.8, 0.5]],
                     [Path.MOVETO, Path.LINETO,
                      Path.MOVETO, Path.LINETO,
                      Path.MOVETO, Path.LINETO])

    def _set_tri_down(self):
        self._transform = Affine2D().scale(0.5)
        self._snap_threshold = 5.0
        self._filled = False
        self._path = self._tri_path

    def _set_tri_up(self):
        self._set_tri_down()
        self._transform = self._transform.rotate_deg(180)

    def _set_tri_left(self):
        self._set_tri_down()
        self._transform = self._transform.rotate_deg(270)

    def _set_tri_right(self):
        self._set_tri_down()
        self._transform = self._transform.rotate_deg(90)

    _caret_path = Path([[-1.0, 1.5], [0.0, 0.0], [1.0, 1.5]])

    def _set_caretdown(self):
        self._transform = Affine2D().scale(0.5)
        self._snap_threshold = 3.0
        self._filled = False
        self._path = self._caret_path
        self._joinstyle = 'miter'

    def _set_caretup(self):
        self._set_caretdown()
        self._transform = self._transform.rotate_deg(180)

    def _set_caretleft(self):
        self._set_caretdown()
        self._transform = self._transform.rotate_deg(270)

    def _set_caretright(self):
        self._set_caretdown()
        self._transform = self._transform.rotate_deg(90)

    _caret_path_base = Path([[-1.0, 0.0], [0.0, -1.5], [1.0, 0]])

    def _set_caretdownbase(self):
        self._set_caretdown()
        self._path = self._caret_path_base

    def _set_caretupbase(self):
        self._set_caretdownbase()
        self._transform = self._transform.rotate_deg(180)

    def _set_caretleftbase(self):
        self._set_caretdownbase()
        self._transform = self._transform.rotate_deg(270)

    def _set_caretrightbase(self):
        self._set_caretdownbase()
        self._transform = self._transform.rotate_deg(90)

    _plus_path = Path([[-1.0, 0.0], [1.0, 0.0],
                       [0.0, -1.0], [0.0, 1.0]],
                      [Path.MOVETO, Path.LINETO,
                       Path.MOVETO, Path.LINETO])

    def _set_plus(self):
        self._transform = Affine2D().scale(0.5)
        self._snap_threshold = 1.0
        self._filled = False
        self._path = self._plus_path

    _x_path = Path([[-1.0, -1.0], [1.0, 1.0],
                    [-1.0, 1.0], [1.0, -1.0]],
                   [Path.MOVETO, Path.LINETO,
                    Path.MOVETO, Path.LINETO])

    def _set_x(self):
        self._transform = Affine2D().scale(0.5)
        self._snap_threshold = 3.0
        self._filled = False
        self._path = self._x_path

    _plus_filled_path = Path([(1/3, 0), (2/3, 0), (2/3, 1/3),
                              (1, 1/3), (1, 2/3), (2/3, 2/3),
                              (2/3, 1), (1/3, 1), (1/3, 2/3),
                              (0, 2/3), (0, 1/3), (1/3, 1/3),
                              (1/3, 0)],
                             [Path.MOVETO, Path.LINETO, Path.LINETO,
                              Path.LINETO, Path.LINETO, Path.LINETO,
                              Path.LINETO, Path.LINETO, Path.LINETO,
                              Path.LINETO, Path.LINETO, Path.LINETO,
                              Path.CLOSEPOLY])

    _plus_filled_path_t = Path([(1, 1/2), (1, 2/3), (2/3, 2/3),
                                (2/3, 1), (1/3, 1), (1/3, 2/3),
                                (0, 2/3), (0, 1/2), (1, 1/2)],
                               [Path.MOVETO, Path.LINETO, Path.LINETO,
                                Path.LINETO, Path.LINETO, Path.LINETO,
                                Path.LINETO, Path.LINETO,
                                Path.CLOSEPOLY])

    def _set_plus_filled(self):
        self._transform = Affine2D().translate(-0.5, -0.5)
        self._snap_threshold = 5.0
        self._joinstyle = 'miter'
        fs = self.get_fillstyle()
        if not self._half_fill():
            self._path = self._plus_filled_path
        else:
            # Rotate top half path to support all partitions
            if fs == 'top':
                rotate, rotate_alt = 0, 180
            elif fs == 'bottom':
                rotate, rotate_alt = 180, 0
            elif fs == 'left':
                rotate, rotate_alt = 90, 270
            else:
                rotate, rotate_alt = 270, 90

            self._path = self._plus_filled_path_t
            self._alt_path = self._plus_filled_path_t
            self._alt_transform = Affine2D().translate(-0.5, -0.5)
            self._transform.rotate_deg(rotate)
            self._alt_transform.rotate_deg(rotate_alt)

    _x_filled_path = Path([(0.25, 0), (0.5, 0.25), (0.75, 0), (1, 0.25),
                           (0.75, 0.5), (1, 0.75), (0.75, 1), (0.5, 0.75),
                           (0.25, 1), (0, 0.75), (0.25, 0.5), (0, 0.25),
                           (0.25, 0)],
                          [Path.MOVETO, Path.LINETO, Path.LINETO,
                           Path.LINETO, Path.LINETO, Path.LINETO,
                           Path.LINETO, Path.LINETO, Path.LINETO,
                           Path.LINETO, Path.LINETO, Path.LINETO,
                           Path.CLOSEPOLY])

    _x_filled_path_t = Path([(0.75, 0.5), (1, 0.75), (0.75, 1),
                             (0.5, 0.75), (0.25, 1), (0, 0.75),
                             (0.25, 0.5), (0.75, 0.5)],
                            [Path.MOVETO, Path.LINETO, Path.LINETO,
                             Path.LINETO, Path.LINETO, Path.LINETO,
                             Path.LINETO, Path.CLOSEPOLY])

    def _set_x_filled(self):
        self._transform = Affine2D().translate(-0.5, -0.5)
        self._snap_threshold = 5.0
        self._joinstyle = 'miter'
        fs = self.get_fillstyle()
        if not self._half_fill():
            self._path = self._x_filled_path
        else:
            # Rotate top half path to support all partitions
            if fs == 'top':
                rotate, rotate_alt = 0, 180
            elif fs == 'bottom':
                rotate, rotate_alt = 180, 0
            elif fs == 'left':
                rotate, rotate_alt = 90, 270
            else:
                rotate, rotate_alt = 270, 90

            self._path = self._x_filled_path_t
            self._alt_path = self._x_filled_path_t
            self._alt_transform = Affine2D().translate(-0.5, -0.5)
            self._transform.rotate_deg(rotate)
            self._alt_transform.rotate_deg(rotate_alt)
