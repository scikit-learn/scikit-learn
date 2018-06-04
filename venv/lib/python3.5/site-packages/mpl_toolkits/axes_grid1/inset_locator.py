"""
A collection of functions and objects for creating or placing inset axes.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from matplotlib import docstring
import six
from matplotlib.offsetbox import AnchoredOffsetbox
from matplotlib.patches import Patch, Rectangle
from matplotlib.path import Path
from matplotlib.transforms import Bbox, BboxTransformTo
from matplotlib.transforms import IdentityTransform, TransformedBbox

from . import axes_size as Size
from .parasite_axes import HostAxes


class InsetPosition(object):
    @docstring.dedent_interpd
    def __init__(self, parent, lbwh):
        """
        An object for positioning an inset axes.

        This is created by specifying the normalized coordinates in the axes,
        instead of the figure.

        Parameters
        ----------
        parent : `matplotlib.axes.Axes`
            Axes to use for normalizing coordinates.

        lbwh : iterable of four floats
            The left edge, bottom edge, width, and height of the inset axes, in
            units of the normalized coordinate of the *parent* axes.

        See Also
        --------
        :meth:`matplotlib.axes.Axes.set_axes_locator`

        Examples
        --------
        The following bounds the inset axes to a box with 20%% of the parent
        axes's height and 40%% of the width. The size of the axes specified
        ([0, 0, 1, 1]) ensures that the axes completely fills the bounding box:

        >>> parent_axes = plt.gca()
        >>> ax_ins = plt.axes([0, 0, 1, 1])
        >>> ip = InsetPosition(ax, [0.5, 0.1, 0.4, 0.2])
        >>> ax_ins.set_axes_locator(ip)
        """
        self.parent = parent
        self.lbwh = lbwh

    def __call__(self, ax, renderer):
        bbox_parent = self.parent.get_position(original=False)
        trans = BboxTransformTo(bbox_parent)
        bbox_inset = Bbox.from_bounds(*self.lbwh)
        bb = TransformedBbox(bbox_inset, trans)
        return bb


class AnchoredLocatorBase(AnchoredOffsetbox):
    def __init__(self, bbox_to_anchor, offsetbox, loc,
                 borderpad=0.5, bbox_transform=None):
        super(AnchoredLocatorBase, self).__init__(
            loc, pad=0., child=None, borderpad=borderpad,
            bbox_to_anchor=bbox_to_anchor, bbox_transform=bbox_transform
        )

    def draw(self, renderer):
        raise RuntimeError("No draw method should be called")

    def __call__(self, ax, renderer):
        self.axes = ax

        fontsize = renderer.points_to_pixels(self.prop.get_size_in_points())
        self._update_offset_func(renderer, fontsize)

        width, height, xdescent, ydescent = self.get_extent(renderer)

        px, py = self.get_offset(width, height, 0, 0, renderer)
        bbox_canvas = Bbox.from_bounds(px, py, width, height)
        tr = ax.figure.transFigure.inverted()
        bb = TransformedBbox(bbox_canvas, tr)

        return bb


class AnchoredSizeLocator(AnchoredLocatorBase):
    def __init__(self, bbox_to_anchor, x_size, y_size, loc,
                 borderpad=0.5, bbox_transform=None):

        super(AnchoredSizeLocator, self).__init__(
            bbox_to_anchor, None, loc,
            borderpad=borderpad, bbox_transform=bbox_transform
        )

        self.x_size = Size.from_any(x_size)
        self.y_size = Size.from_any(y_size)

    def get_extent(self, renderer):
        x, y, w, h = self.get_bbox_to_anchor().bounds

        dpi = renderer.points_to_pixels(72.)

        r, a = self.x_size.get_size(renderer)
        width = w*r + a*dpi

        r, a = self.y_size.get_size(renderer)
        height = h*r + a*dpi
        xd, yd = 0, 0

        fontsize = renderer.points_to_pixels(self.prop.get_size_in_points())
        pad = self.pad * fontsize

        return width+2*pad, height+2*pad, xd+pad, yd+pad


class AnchoredZoomLocator(AnchoredLocatorBase):
    def __init__(self, parent_axes, zoom, loc,
                 borderpad=0.5,
                 bbox_to_anchor=None,
                 bbox_transform=None):

        self.parent_axes = parent_axes
        self.zoom = zoom

        if bbox_to_anchor is None:
            bbox_to_anchor = parent_axes.bbox

        super(AnchoredZoomLocator, self).__init__(
            bbox_to_anchor, None, loc, borderpad=borderpad,
            bbox_transform=bbox_transform)

    def get_extent(self, renderer):
        bb = TransformedBbox(self.axes.viewLim,
                             self.parent_axes.transData)

        x, y, w, h = bb.bounds
        fontsize = renderer.points_to_pixels(self.prop.get_size_in_points())
        pad = self.pad * fontsize

        return abs(w*self.zoom)+2*pad, abs(h*self.zoom)+2*pad, pad, pad


class BboxPatch(Patch):
    @docstring.dedent_interpd
    def __init__(self, bbox, **kwargs):
        """
        Patch showing the shape bounded by a Bbox.

        Parameters
        ----------
        bbox : `matplotlib.transforms.Bbox`
            Bbox to use for the extents of this patch.

        **kwargs
            Patch properties. Valid arguments include:
            %(Patch)s
        """
        if "transform" in kwargs:
            raise ValueError("transform should not be set")

        kwargs["transform"] = IdentityTransform()
        Patch.__init__(self, **kwargs)
        self.bbox = bbox

    def get_path(self):
        x0, y0, x1, y1 = self.bbox.extents

        verts = [(x0, y0),
                 (x1, y0),
                 (x1, y1),
                 (x0, y1),
                 (x0, y0),
                 (0, 0)]

        codes = [Path.MOVETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.CLOSEPOLY]

        return Path(verts, codes)
    get_path.__doc__ = Patch.get_path.__doc__


class BboxConnector(Patch):
    @staticmethod
    def get_bbox_edge_pos(bbox, loc):
        """
        Helper function to obtain the location of a corner of a bbox

        Parameters
        ----------
        bbox : `matplotlib.transforms.Bbox`

        loc : {1, 2, 3, 4}
            Corner of *bbox*. Valid values are::

                'upper right'  : 1,
                'upper left'   : 2,
                'lower left'   : 3,
                'lower right'  : 4

        Returns
        -------
        x, y : float
            Coordinates of the corner specified by *loc*.
        """
        x0, y0, x1, y1 = bbox.extents
        if loc == 1:
            return x1, y1
        elif loc == 2:
            return x0, y1
        elif loc == 3:
            return x0, y0
        elif loc == 4:
            return x1, y0

    @staticmethod
    def connect_bbox(bbox1, bbox2, loc1, loc2=None):
        """
        Helper function to obtain a Path from one bbox to another.

        Parameters
        ----------
        bbox1, bbox2 : `matplotlib.transforms.Bbox`
            Bounding boxes to connect.

        loc1 : {1, 2, 3, 4}
            Corner of *bbox1* to use. Valid values are::

                'upper right'  : 1,
                'upper left'   : 2,
                'lower left'   : 3,
                'lower right'  : 4

        loc2 : {1, 2, 3, 4}, optional
            Corner of *bbox2* to use. If None, defaults to *loc1*.
            Valid values are::

                'upper right'  : 1,
                'upper left'   : 2,
                'lower left'   : 3,
                'lower right'  : 4

        Returns
        -------
        path : `matplotlib.path.Path`
            A line segment from the *loc1* corner of *bbox1* to the *loc2*
            corner of *bbox2*.
        """
        if isinstance(bbox1, Rectangle):
            transform = bbox1.get_transfrom()
            bbox1 = Bbox.from_bounds(0, 0, 1, 1)
            bbox1 = TransformedBbox(bbox1, transform)

        if isinstance(bbox2, Rectangle):
            transform = bbox2.get_transform()
            bbox2 = Bbox.from_bounds(0, 0, 1, 1)
            bbox2 = TransformedBbox(bbox2, transform)

        if loc2 is None:
            loc2 = loc1

        x1, y1 = BboxConnector.get_bbox_edge_pos(bbox1, loc1)
        x2, y2 = BboxConnector.get_bbox_edge_pos(bbox2, loc2)

        verts = [[x1, y1], [x2, y2]]
        codes = [Path.MOVETO, Path.LINETO]

        return Path(verts, codes)

    @docstring.dedent_interpd
    def __init__(self, bbox1, bbox2, loc1, loc2=None, **kwargs):
        """
        Connect two bboxes with a straight line.

        Parameters
        ----------
        bbox1, bbox2 : `matplotlib.transforms.Bbox`
            Bounding boxes to connect.

        loc1 : {1, 2, 3, 4}
            Corner of *bbox1* to draw the line. Valid values are::

                'upper right'  : 1,
                'upper left'   : 2,
                'lower left'   : 3,
                'lower right'  : 4

        loc2 : {1, 2, 3, 4}, optional
            Corner of *bbox2* to draw the line. If None, defaults to *loc1*.
            Valid values are::

                'upper right'  : 1,
                'upper left'   : 2,
                'lower left'   : 3,
                'lower right'  : 4

        **kwargs
            Patch properties for the line drawn. Valid arguments include:
            %(Patch)s
        """
        if "transform" in kwargs:
            raise ValueError("transform should not be set")

        kwargs["transform"] = IdentityTransform()
        Patch.__init__(self, fill=False, **kwargs)
        self.bbox1 = bbox1
        self.bbox2 = bbox2
        self.loc1 = loc1
        self.loc2 = loc2

    def get_path(self):
        return self.connect_bbox(self.bbox1, self.bbox2,
                                 self.loc1, self.loc2)
    get_path.__doc__ = Patch.get_path.__doc__


class BboxConnectorPatch(BboxConnector):
    @docstring.dedent_interpd
    def __init__(self, bbox1, bbox2, loc1a, loc2a, loc1b, loc2b, **kwargs):
        """
        Connect two bboxes with a quadrilateral.

        The quadrilateral is specified by two lines that start and end at corners
        of the bboxes. The four sides of the quadrilateral are defined by the two
        lines given, the line between the two corners specified in *bbox1* and the
        line between the two corners specified in *bbox2*.

        Parameters
        ----------
        bbox1, bbox2 : `matplotlib.transforms.Bbox`
            Bounding boxes to connect.

        loc1a, loc2a : {1, 2, 3, 4}
            Corners of *bbox1* and *bbox2* to draw the first line.
            Valid values are::

                'upper right'  : 1,
                'upper left'   : 2,
                'lower left'   : 3,
                'lower right'  : 4

        loc1b, loc2b : {1, 2, 3, 4}
            Corners of *bbox1* and *bbox2* to draw the second line.
            Valid values are::

                'upper right'  : 1,
                'upper left'   : 2,
                'lower left'   : 3,
                'lower right'  : 4

        **kwargs
            Patch properties for the line drawn:
            %(Patch)s
        """
        if "transform" in kwargs:
            raise ValueError("transform should not be set")
        BboxConnector.__init__(self, bbox1, bbox2, loc1a, loc2a, **kwargs)
        self.loc1b = loc1b
        self.loc2b = loc2b

    def get_path(self):
        path1 = self.connect_bbox(self.bbox1, self.bbox2, self.loc1, self.loc2)
        path2 = self.connect_bbox(self.bbox2, self.bbox1,
                                  self.loc2b, self.loc1b)
        path_merged = (list(path1.vertices) +
                       list(path2.vertices) +
                       [path1.vertices[0]])
        return Path(path_merged)
    get_path.__doc__ = BboxConnector.get_path.__doc__


def _add_inset_axes(parent_axes, inset_axes):
    """Helper function to add an inset axes and disable navigation in it"""
    parent_axes.figure.add_axes(inset_axes)
    inset_axes.set_navigate(False)


@docstring.dedent_interpd
def inset_axes(parent_axes, width, height, loc=1,
               bbox_to_anchor=None, bbox_transform=None,
               axes_class=None,
               axes_kwargs=None,
               borderpad=0.5):
    """
    Create an inset axes with a given width and height.

    Both sizes used can be specified either in inches or percentage of the
    parent axes.

    Parameters
    ----------
    parent_axes : `matplotlib.axes.Axes`
        Axes to place the inset axes.

    width, height : float or str
        Size of the inset axes to create.

    loc : int or string, optional, default to 1
        Location to place the inset axes. The valid locations are::

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

    bbox_to_anchor : tuple or `matplotlib.transforms.BboxBase`, optional
        Bbox that the inset axes will be anchored. Can be a tuple of
        [left, bottom, width, height], or a tuple of [left, bottom].

    bbox_transform : `matplotlib.transforms.Transform`, optional
        Transformation for the bbox. if None, `parent_axes.transAxes` is used.

    axes_class : `matplotlib.axes.Axes` type, optional
        If specified, the inset axes created with be created with this class's
        constructor.

    axes_kwargs : dict, optional
        Keyworded arguments to pass to the constructor of the inset axes.
        Valid arguments include:
        %(Axes)s

    borderpad : float, optional
        Padding between inset axes and the bbox_to_anchor. Defaults to 0.5.

    Returns
    -------
    inset_axes : `axes_class`
        Inset axes object created.
    """

    if axes_class is None:
        axes_class = HostAxes

    if axes_kwargs is None:
        inset_axes = axes_class(parent_axes.figure, parent_axes.get_position())
    else:
        inset_axes = axes_class(parent_axes.figure, parent_axes.get_position(),
                                **axes_kwargs)

    if bbox_to_anchor is None:
        bbox_to_anchor = parent_axes.bbox

    axes_locator = AnchoredSizeLocator(bbox_to_anchor,
                                       width, height,
                                       loc=loc,
                                       bbox_transform=bbox_transform,
                                       borderpad=borderpad)

    inset_axes.set_axes_locator(axes_locator)

    _add_inset_axes(parent_axes, inset_axes)

    return inset_axes


@docstring.dedent_interpd
def zoomed_inset_axes(parent_axes, zoom, loc=1,
                      bbox_to_anchor=None, bbox_transform=None,
                      axes_class=None,
                      axes_kwargs=None,
                      borderpad=0.5):
    """
    Create an anchored inset axes by scaling a parent axes.

    Parameters
    ----------
    parent_axes : `matplotlib.axes.Axes`
        Axes to place the inset axes.

    zoom : float
        Scaling factor of the data axes. *zoom* > 1 will enlargen the
        coordinates (i.e., "zoomed in"), while *zoom* < 1 will shrink the
        coordinates (i.e., "zoomed out").

    loc : int or string, optional, default to 1
        Location to place the inset axes. The valid locations are::

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

    bbox_to_anchor : tuple or `matplotlib.transforms.BboxBase`, optional
        Bbox that the inset axes will be anchored. Can be a tuple of
        [left, bottom, width, height], or a tuple of [left, bottom].

    bbox_transform : `matplotlib.transforms.Transform`, optional
        Transformation for the bbox. if None, `parent_axes.transAxes` is used.

    axes_class : `matplotlib.axes.Axes` type, optional
        If specified, the inset axes created with be created with this class's
        constructor.

    axes_kwargs : dict, optional
        Keyworded arguments to pass to the constructor of the inset axes.
        Valid arguments include:
        %(Axes)s

    borderpad : float, optional
        Padding between inset axes and the bbox_to_anchor. Defaults to 0.5.

    Returns
    -------
    inset_axes : `axes_class`
        Inset axes object created.
    """

    if axes_class is None:
        axes_class = HostAxes

    if axes_kwargs is None:
        inset_axes = axes_class(parent_axes.figure, parent_axes.get_position())
    else:
        inset_axes = axes_class(parent_axes.figure, parent_axes.get_position(),
                                **axes_kwargs)

    axes_locator = AnchoredZoomLocator(parent_axes, zoom=zoom, loc=loc,
                                       bbox_to_anchor=bbox_to_anchor,
                                       bbox_transform=bbox_transform,
                                       borderpad=borderpad)
    inset_axes.set_axes_locator(axes_locator)

    _add_inset_axes(parent_axes, inset_axes)

    return inset_axes


@docstring.dedent_interpd
def mark_inset(parent_axes, inset_axes, loc1, loc2, **kwargs):
    """
    Draw a box to mark the location of an area represented by an inset axes.

    This function draws a box in *parent_axes* at the bounding box of
    *inset_axes*, and shows a connection with the inset axes by drawing lines
    at the corners, giving a "zoomed in" effect.

    Parameters
    ----------
    parent_axes : `matplotlib.axes.Axes`
        Axes which contains the area of the inset axes.

    inset_axes : `matplotlib.axes.Axes`
        The inset axes.

    loc1, loc2 : {1, 2, 3, 4}
        Corners to use for connecting the inset axes and the area in the
        parent axes.

    **kwargs
        Patch properties for the lines and box drawn:
        %(Patch)s

    Returns
    -------
    pp : `matplotlib.patches.Patch`
        The patch drawn to represent the area of the inset axes.

    p1, p2 : `matplotlib.patches.Patch`
        The patches connecting two corners of the inset axes and its area.
    """
    rect = TransformedBbox(inset_axes.viewLim, parent_axes.transData)

    fill = kwargs.pop("fill", False)
    pp = BboxPatch(rect, fill=fill, **kwargs)
    parent_axes.add_patch(pp)

    p1 = BboxConnector(inset_axes.bbox, rect, loc1=loc1, **kwargs)
    inset_axes.add_patch(p1)
    p1.set_clip_on(False)
    p2 = BboxConnector(inset_axes.bbox, rect, loc1=loc2, **kwargs)
    inset_axes.add_patch(p2)
    p2.set_clip_on(False)

    return pp, p1, p2
