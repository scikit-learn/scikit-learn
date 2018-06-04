from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

from matplotlib.patches import _Style, FancyArrowPatch
from matplotlib.transforms import IdentityTransform
from matplotlib.path import Path
import numpy as np

class _FancyAxislineStyle(object):
    class SimpleArrow(FancyArrowPatch):
        """
        The artist class that will be returned for SimpleArrow style.
        """
        _ARROW_STYLE = "->"

        def __init__(self, axis_artist, line_path, transform,
                     line_mutation_scale):
            self._axis_artist = axis_artist
            self._line_transform = transform
            self._line_path = line_path
            self._line_mutation_scale = line_mutation_scale

            FancyArrowPatch.__init__(self,
                                     path=self._line_path,
                                     arrowstyle=self._ARROW_STYLE,
                                     arrow_transmuter=None,
                                     patchA=None,
                                     patchB=None,
                                     shrinkA=0.,
                                     shrinkB=0.,
                                     mutation_scale=line_mutation_scale,
                                     mutation_aspect=None,
                                     transform=IdentityTransform(),
                                     )

        def set_line_mutation_scale(self, scale):
            self.set_mutation_scale(scale*self._line_mutation_scale)

        def _extend_path(self, path, mutation_size=10):
            """
            Extend the path to make a room for drawing arrow.
            """
            from matplotlib.bezier import get_cos_sin

            x0, y0 = path.vertices[-2]
            x1, y1 = path.vertices[-1]
            cost, sint = get_cos_sin(x0, y0, x1, y1)

            d = mutation_size * 1.
            x2, y2 = x1 + cost*d, y1+sint*d

            if path.codes is None:
                _path = Path(np.concatenate([path.vertices, [[x2, y2]]]))
            else:
                _path = Path(np.concatenate([path.vertices, [[x2, y2]]]),
                             np.concatenate([path.codes, [Path.LINETO]]))

            return _path

        def set_path(self, path):
            self._line_path = path

        def draw(self, renderer):
            """
            Draw the axis line.
             1) transform the path to the display coordinate.
             2) extend the path to make a room for arrow
             3) update the path of the FancyArrowPatch.
             4) draw
            """
            path_in_disp = self._line_transform.transform_path(self._line_path)
            mutation_size = self.get_mutation_scale() #line_mutation_scale()
            extented_path = self._extend_path(path_in_disp,
                                              mutation_size=mutation_size)

            self._path_original = extented_path
            FancyArrowPatch.draw(self, renderer)

    class FilledArrow(SimpleArrow):
        """
        The artist class that will be returned for SimpleArrow style.
        """
        _ARROW_STYLE = "-|>"


class AxislineStyle(_Style):
    """
    :class:`AxislineStyle` is a container class which defines style classes
    for AxisArtists.

    An instance of any axisline style class is an callable object,
    whose call signature is ::

       __call__(self, axis_artist, path, transform)

    When called, this should return a mpl artist with following
    methods implemented. ::

      def set_path(self, path):
          # set the path for axisline.

      def set_line_mutation_scale(self, scale):
          # set the scale

      def draw(self, renderer):
          # draw


    """

    _style_list = {}


    class _Base(object):
        # The derived classes are required to be able to be initialized
        # w/o arguments, i.e., all its argument (except self) must have
        # the default values.

        def __init__(self):
            """
            initialization.
            """
            super(AxislineStyle._Base, self).__init__()




        def __call__(self, axis_artist, transform):
            """
            Given the AxisArtist instance, and transform for the path
            (set_path method), return the mpl artist for drawing the axis line.
            """

            return self.new_line(axis_artist, transform)


    class SimpleArrow(_Base):
        """
        A simple arrow.
        """

        ArrowAxisClass = _FancyAxislineStyle.SimpleArrow

        def __init__(self, size=1):
            """
             *size*
                size of the arrow as a fraction of the ticklabel size.
            """

            self.size = size
            super(AxislineStyle.SimpleArrow, self).__init__()

        def new_line(self, axis_artist, transform):

            linepath = Path([(0,0), (0, 1)])
            axisline = self.ArrowAxisClass(axis_artist, linepath, transform,
                                           line_mutation_scale=self.size)
            return axisline


    _style_list["->"] = SimpleArrow

    class FilledArrow(SimpleArrow):
        ArrowAxisClass = _FancyAxislineStyle.FilledArrow

    _style_list["-|>"] = FilledArrow
