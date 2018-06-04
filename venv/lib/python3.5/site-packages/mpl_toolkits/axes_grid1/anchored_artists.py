from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import six

from matplotlib import docstring
from matplotlib.offsetbox import (AnchoredOffsetbox, AuxTransformBox,
                                  DrawingArea, TextArea, VPacker)
from matplotlib.patches import Rectangle, Ellipse


__all__ = ['AnchoredDrawingArea', 'AnchoredAuxTransformBox',
           'AnchoredEllipse', 'AnchoredSizeBar']


class AnchoredDrawingArea(AnchoredOffsetbox):
    @docstring.dedent
    def __init__(self, width, height, xdescent, ydescent,
                 loc, pad=0.4, borderpad=0.5, prop=None, frameon=True,
                 **kwargs):
        """
        An anchored container with a fixed size and fillable DrawingArea.

        Artists added to the *drawing_area* will have their coordinates
        interpreted as pixels. Any transformations set on the artists will be
        overridden.

        Parameters
        ----------
        width, height : int or float
            width and height of the container, in pixels.

        xdescent, ydescent : int or float
            descent of the container in the x- and y- direction, in pixels.

        loc : int
            Location of this artist. Valid location codes are::

                'upper right'  : 1,
                'upper left'   : 2,
                'lower left'   : 3,
                'lower right'  : 4,
                'right'        : 5,
                'center left'  : 6,
                'center right' : 7,
                'lower center' : 8,
                'upper center' : 9,
                'center'       : 10

        pad : int or float, optional
            Padding around the child objects, in fraction of the font
            size. Defaults to 0.4.

        borderpad : int or float, optional
            Border padding, in fraction of the font size.
            Defaults to 0.5.

        prop : `matplotlib.font_manager.FontProperties`, optional
            Font property used as a reference for paddings.

        frameon : bool, optional
            If True, draw a box around this artists. Defaults to True.

        **kwargs :
            Keyworded arguments to pass to
            :class:`matplotlib.offsetbox.AnchoredOffsetbox`.

        Attributes
        ----------
        drawing_area : `matplotlib.offsetbox.DrawingArea`
            A container for artists to display.

        Examples
        --------
        To display blue and red circles of different sizes in the upper right
        of an axes *ax*:

        >>> ada = AnchoredDrawingArea(20, 20, 0, 0, loc=1, frameon=False)
        >>> ada.drawing_area.add_artist(Circle((10, 10), 10, fc="b"))
        >>> ada.drawing_area.add_artist(Circle((30, 10), 5, fc="r"))
        >>> ax.add_artist(ada)
        """
        self.da = DrawingArea(width, height, xdescent, ydescent)
        self.drawing_area = self.da

        super(AnchoredDrawingArea, self).__init__(
            loc, pad=pad, borderpad=borderpad, child=self.da, prop=None,
            frameon=frameon, **kwargs
        )


class AnchoredAuxTransformBox(AnchoredOffsetbox):
    @docstring.dedent
    def __init__(self, transform, loc,
                 pad=0.4, borderpad=0.5, prop=None, frameon=True, **kwargs):
        """
        An anchored container with transformed coordinates.

        Artists added to the *drawing_area* are scaled according to the
        coordinates of the transformation used. The dimensions of this artist
        will scale to contain the artists added.

        Parameters
        ----------
        transform : `matplotlib.transforms.Transform`
            The transformation object for the coordinate system in use, i.e.,
            :attr:`matplotlib.axes.Axes.transData`.

        loc : int
            Location of this artist. Valid location codes are::

                'upper right'  : 1,
                'upper left'   : 2,
                'lower left'   : 3,
                'lower right'  : 4,
                'right'        : 5,
                'center left'  : 6,
                'center right' : 7,
                'lower center' : 8,
                'upper center' : 9,
                'center'       : 10

        pad : int or float, optional
            Padding around the child objects, in fraction of the font
            size. Defaults to 0.4.

        borderpad : int or float, optional
            Border padding, in fraction of the font size.
            Defaults to 0.5.

        prop : `matplotlib.font_manager.FontProperties`, optional
            Font property used as a reference for paddings.

        frameon : bool, optional
            If True, draw a box around this artists. Defaults to True.

        **kwargs :
            Keyworded arguments to pass to
            :class:`matplotlib.offsetbox.AnchoredOffsetbox`.

        Attributes
        ----------
        drawing_area : `matplotlib.offsetbox.AuxTransformBox`
            A container for artists to display.

        Examples
        --------
        To display an ellipse in the upper left, with a width of 0.1 and
        height of 0.4 in data coordinates:

        >>> box = AnchoredAuxTransformBox(ax.transData, loc=2)
        >>> el = Ellipse((0,0), width=0.1, height=0.4, angle=30)
        >>> box.drawing_area.add_artist(el)
        >>> ax.add_artist(box)
        """
        self.drawing_area = AuxTransformBox(transform)

        AnchoredOffsetbox.__init__(self, loc, pad=pad, borderpad=borderpad,
                                   child=self.drawing_area,
                                   prop=prop,
                                   frameon=frameon,
                                   **kwargs)


class AnchoredEllipse(AnchoredOffsetbox):
    @docstring.dedent
    def __init__(self, transform, width, height, angle, loc,
                 pad=0.1, borderpad=0.1, prop=None, frameon=True, **kwargs):
        """
        Draw an anchored ellipse of a given size.

        Parameters
        ----------
        transform : `matplotlib.transforms.Transform`
            The transformation object for the coordinate system in use, i.e.,
            :attr:`matplotlib.axes.Axes.transData`.

        width, height : int or float
            Width and height of the ellipse, given in coordinates of
            *transform*.

        angle : int or float
            Rotation of the ellipse, in degrees, anti-clockwise.

        loc : int
            Location of this size bar. Valid location codes are::

                'upper right'  : 1,
                'upper left'   : 2,
                'lower left'   : 3,
                'lower right'  : 4,
                'right'        : 5,
                'center left'  : 6,
                'center right' : 7,
                'lower center' : 8,
                'upper center' : 9,
                'center'       : 10

        pad : int or float, optional
            Padding around the ellipse, in fraction of the font size. Defaults
            to 0.1.

        borderpad : int or float, optional
            Border padding, in fraction of the font size. Defaults to 0.1.

        frameon : bool, optional
            If True, draw a box around the ellipse. Defaults to True.

        prop : `matplotlib.font_manager.FontProperties`, optional
            Font property used as a reference for paddings.

        **kwargs :
            Keyworded arguments to pass to
            :class:`matplotlib.offsetbox.AnchoredOffsetbox`.

        Attributes
        ----------
        ellipse : `matplotlib.patches.Ellipse`
            Ellipse patch drawn.
        """
        self._box = AuxTransformBox(transform)
        self.ellipse = Ellipse((0, 0), width, height, angle)
        self._box.add_artist(self.ellipse)

        AnchoredOffsetbox.__init__(self, loc, pad=pad, borderpad=borderpad,
                                   child=self._box,
                                   prop=prop,
                                   frameon=frameon, **kwargs)


class AnchoredSizeBar(AnchoredOffsetbox):
    @docstring.dedent
    def __init__(self, transform, size, label, loc,
                 pad=0.1, borderpad=0.1, sep=2,
                 frameon=True, size_vertical=0, color='black',
                 label_top=False, fontproperties=None, fill_bar=None,
                 **kwargs):
        """
        Draw a horizontal scale bar with a center-aligned label underneath.

        Parameters
        ----------
        transform : `matplotlib.transforms.Transform`
            The transformation object for the coordinate system in use, i.e.,
            :attr:`matplotlib.axes.Axes.transData`.

        size : int or float
            Horizontal length of the size bar, given in coordinates of
            *transform*.

        label : str
            Label to display.

        loc : int
            Location of this size bar. Valid location codes are::

                'upper right'  : 1,
                'upper left'   : 2,
                'lower left'   : 3,
                'lower right'  : 4,
                'right'        : 5,
                'center left'  : 6,
                'center right' : 7,
                'lower center' : 8,
                'upper center' : 9,
                'center'       : 10

        pad : int or float, optional
            Padding around the label and size bar, in fraction of the font
            size. Defaults to 0.1.

        borderpad : int or float, optional
            Border padding, in fraction of the font size.
            Defaults to 0.1.

        sep : int or float, optional
            Separation between the label and the size bar, in points.
            Defaults to 2.

        frameon : bool, optional
            If True, draw a box around the horizontal bar and label.
            Defaults to True.

        size_vertical : int or float, optional
            Vertical length of the size bar, given in coordinates of
            *transform*. Defaults to 0.

        color : str, optional
            Color for the size bar and label.
            Defaults to black.

        label_top : bool, optional
            If True, the label will be over the size bar.
            Defaults to False.

        fontproperties : `matplotlib.font_manager.FontProperties`, optional
            Font properties for the label text.

        fill_bar : bool, optional
            If True and if size_vertical is nonzero, the size bar will
            be filled in with the color specified by the size bar.
            Defaults to True if `size_vertical` is greater than
            zero and False otherwise.

        **kwargs :
            Keyworded arguments to pass to
            :class:`matplotlib.offsetbox.AnchoredOffsetbox`.

        Attributes
        ----------
        size_bar : `matplotlib.offsetbox.AuxTransformBox`
            Container for the size bar.

        txt_label : `matplotlib.offsetbox.TextArea`
            Container for the label of the size bar.

        Notes
        -----
        If *prop* is passed as a keyworded argument, but *fontproperties* is
        not, then *prop* is be assumed to be the intended *fontproperties*.
        Using both *prop* and *fontproperties* is not supported.

        Examples
        --------
        >>> import matplotlib.pyplot as plt
        >>> import numpy as np
        >>> from mpl_toolkits.axes_grid1.anchored_artists import \
AnchoredSizeBar
        >>> fig, ax = plt.subplots()
        >>> ax.imshow(np.random.random((10,10)))
        >>> bar = AnchoredSizeBar(ax.transData, 3, '3 data units', 4)
        >>> ax.add_artist(bar)
        >>> fig.show()

        Using all the optional parameters

        >>> import matplotlib.font_manager as fm
        >>> fontprops = fm.FontProperties(size=14, family='monospace')
        >>> bar = AnchoredSizeBar(ax.transData, 3, '3 units', 4, pad=0.5, \
sep=5, borderpad=0.5, frameon=False, \
size_vertical=0.5, color='white', \
fontproperties=fontprops)
        """
        if fill_bar is None:
            fill_bar = size_vertical > 0

        self.size_bar = AuxTransformBox(transform)
        self.size_bar.add_artist(Rectangle((0, 0), size, size_vertical,
                                           fill=fill_bar, facecolor=color,
                                           edgecolor=color))

        if fontproperties is None and 'prop' in kwargs:
            fontproperties = kwargs.pop('prop')

        if fontproperties is None:
            textprops = {'color': color}
        else:
            textprops = {'color': color, 'fontproperties': fontproperties}

        self.txt_label = TextArea(
            label,
            minimumdescent=False,
            textprops=textprops)

        if label_top:
            _box_children = [self.txt_label, self.size_bar]
        else:
            _box_children = [self.size_bar, self.txt_label]

        self._box = VPacker(children=_box_children,
                            align="center",
                            pad=0, sep=sep)

        AnchoredOffsetbox.__init__(self, loc, pad=pad, borderpad=borderpad,
                                   child=self._box,
                                   prop=fontproperties,
                                   frameon=frameon, **kwargs)
