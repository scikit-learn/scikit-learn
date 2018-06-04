"""
The OffsetBox is a simple container artist. The child artist are meant
to be drawn at a relative position to its parent.  The [VH]Packer,
DrawingArea and TextArea are derived from the OffsetBox.

The [VH]Packer automatically adjust the relative postisions of their
children, which should be instances of the OffsetBox. This is used to
align similar artists together, e.g., in legend.

The DrawingArea can contain any Artist as a child. The
DrawingArea has a fixed width and height. The position of children
relative to the parent is fixed.  The TextArea is contains a single
Text instance. The width and height of the TextArea instance is the
width and height of the its child text.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves import xrange, zip

import warnings
import matplotlib.transforms as mtransforms
import matplotlib.artist as martist
import matplotlib.text as mtext
import matplotlib.path as mpath
import numpy as np
from matplotlib.transforms import Bbox, BboxBase, TransformedBbox

from matplotlib.font_manager import FontProperties
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib import rcParams

from matplotlib import docstring

from matplotlib.image import BboxImage

from matplotlib.patches import bbox_artist as mbbox_artist
from matplotlib.text import _AnnotationBase


DEBUG = False


# for debugging use
def bbox_artist(*args, **kwargs):
    if DEBUG:
        mbbox_artist(*args, **kwargs)

# _get_packed_offsets() and _get_aligned_offsets() are coded assuming
# that we are packing boxes horizontally. But same function will be
# used with vertical packing.


def _get_packed_offsets(wd_list, total, sep, mode="fixed"):
    """
    Geiven a list of (width, xdescent) of each boxes, calculate the
    total width and the x-offset positions of each items according to
    *mode*. xdescent is analogous to the usual descent, but along the
    x-direction. xdescent values are currently ignored.

    *wd_list* : list of (width, xdescent) of boxes to be packed.
    *sep* : spacing between boxes
    *total* : Intended total length. None if not used.
    *mode* : packing mode. 'fixed', 'expand', or 'equal'.
    """

    w_list, d_list = zip(*wd_list)
    # d_list is currently not used.

    if mode == "fixed":
        offsets_ = np.cumsum([0] + [w + sep for w in w_list])
        offsets = offsets_[:-1]
        if total is None:
            total = offsets_[-1] - sep
        return total, offsets

    elif mode == "expand":
        # This is a bit of a hack to avoid a TypeError when *total*
        # is None and used in conjugation with tight layout.
        if total is None:
            total = 1
        if len(w_list) > 1:
            sep = (total - sum(w_list)) / (len(w_list) - 1)
        else:
            sep = 0
        offsets_ = np.cumsum([0] + [w + sep for w in w_list])
        offsets = offsets_[:-1]
        return total, offsets

    elif mode == "equal":
        maxh = max(w_list)
        if total is None:
            total = (maxh + sep) * len(w_list)
        else:
            sep = total / len(w_list) - maxh
        offsets = (maxh + sep) * np.arange(len(w_list))
        return total, offsets

    else:
        raise ValueError("Unknown mode : %s" % (mode,))


def _get_aligned_offsets(hd_list, height, align="baseline"):
    """
    Given a list of (height, descent) of each boxes, align the boxes
    with *align* and calculate the y-offsets of each boxes.
    total width and the offset positions of each items according to
    *mode*. xdescent is analogous to the usual descent, but along the
    x-direction. xdescent values are currently ignored.

    *hd_list* : list of (width, xdescent) of boxes to be aligned.
    *sep* : spacing between boxes
    *height* : Intended total length. None if not used.
    *align* : align mode. 'baseline', 'top', 'bottom', or 'center'.
    """

    if height is None:
        height = max(h for h, d in hd_list)

    if align == "baseline":
        height_descent = max(h - d for h, d in hd_list)
        descent = max(d for h, d in hd_list)
        height = height_descent + descent
        offsets = [0. for h, d in hd_list]
    elif align in ["left", "top"]:
        descent = 0.
        offsets = [d for h, d in hd_list]
    elif align in ["right", "bottom"]:
        descent = 0.
        offsets = [height - h + d for h, d in hd_list]
    elif align == "center":
        descent = 0.
        offsets = [(height - h) * .5 + d for h, d in hd_list]
    else:
        raise ValueError("Unknown Align mode : %s" % (align,))

    return height, descent, offsets


class OffsetBox(martist.Artist):
    """
    The OffsetBox is a simple container artist. The child artist are meant
    to be drawn at a relative position to its parent.
    """
    def __init__(self, *args, **kwargs):

        super(OffsetBox, self).__init__(*args, **kwargs)

        # Clipping has not been implemented in the OffesetBox family, so
        # disable the clip flag for consistency. It can always be turned back
        # on to zero effect.
        self.set_clip_on(False)

        self._children = []
        self._offset = (0, 0)

    def __getstate__(self):
        state = martist.Artist.__getstate__(self)

        # pickle cannot save instancemethods, so handle them here
        from .cbook import _InstanceMethodPickler
        import inspect

        offset = state['_offset']
        if inspect.ismethod(offset):
            state['_offset'] = _InstanceMethodPickler(offset)
        return state

    def __setstate__(self, state):
        self.__dict__ = state
        from .cbook import _InstanceMethodPickler
        if isinstance(self._offset, _InstanceMethodPickler):
            self._offset = self._offset.get_instancemethod()
        self.stale = True

    def set_figure(self, fig):
        """
        Set the figure

        accepts a class:`~matplotlib.figure.Figure` instance
        """
        martist.Artist.set_figure(self, fig)
        for c in self.get_children():
            c.set_figure(fig)

    @martist.Artist.axes.setter
    def axes(self, ax):
        # TODO deal with this better
        martist.Artist.axes.fset(self, ax)
        for c in self.get_children():
            if c is not None:
                c.axes = ax

    def contains(self, mouseevent):
        for c in self.get_children():
            a, b = c.contains(mouseevent)
            if a:
                return a, b
        return False, {}

    def set_offset(self, xy):
        """
        Set the offset

        accepts x, y, tuple, or a callable object.
        """
        self._offset = xy
        self.stale = True

    def get_offset(self, width, height, xdescent, ydescent, renderer):
        """
        Get the offset

        accepts extent of the box
        """
        return (self._offset(width, height, xdescent, ydescent, renderer)
                if callable(self._offset)
                else self._offset)

    def set_width(self, width):
        """
        Set the width

        accepts float
        """
        self.width = width
        self.stale = True

    def set_height(self, height):
        """
        Set the height

        accepts float
        """
        self.height = height
        self.stale = True

    def get_visible_children(self):
        """
        Return a list of visible artists it contains.
        """
        return [c for c in self._children if c.get_visible()]

    def get_children(self):
        """
        Return a list of artists it contains.
        """
        return self._children

    def get_extent_offsets(self, renderer):
        raise Exception("")

    def get_extent(self, renderer):
        """
        Return with, height, xdescent, ydescent of box
        """
        w, h, xd, yd, offsets = self.get_extent_offsets(renderer)
        return w, h, xd, yd

    def get_window_extent(self, renderer):
        '''
        get the bounding box in display space.
        '''
        w, h, xd, yd, offsets = self.get_extent_offsets(renderer)
        px, py = self.get_offset(w, h, xd, yd, renderer)
        return mtransforms.Bbox.from_bounds(px - xd, py - yd, w, h)

    def draw(self, renderer):
        """
        Update the location of children if necessary and draw them
        to the given *renderer*.
        """

        width, height, xdescent, ydescent, offsets = self.get_extent_offsets(
                                                        renderer)

        px, py = self.get_offset(width, height, xdescent, ydescent, renderer)

        for c, (ox, oy) in zip(self.get_visible_children(), offsets):
            c.set_offset((px + ox, py + oy))
            c.draw(renderer)

        bbox_artist(self, renderer, fill=False, props=dict(pad=0.))
        self.stale = False


class PackerBase(OffsetBox):
    def __init__(self, pad=None, sep=None, width=None, height=None,
                 align=None, mode=None,
                 children=None):
        """
        Parameters
        ----------
        pad : float, optional
            Boundary pad.

        sep : float, optional
            Spacing between items.

        width : float, optional

        height : float, optional
           Width and height of the container box, calculated if
           `None`.

        align : str, optional
            Alignment of boxes. Can be one of ``top``, ``bottom``,
            ``left``, ``right``, ``center`` and ``baseline``

        mode : str, optional
            Packing mode.

        Notes
        -----
        *pad* and *sep* need to given in points and will be scale with
        the renderer dpi, while *width* and *height* need to be in
        pixels.
        """
        super(PackerBase, self).__init__()

        self.height = height
        self.width = width
        self.sep = sep
        self.pad = pad
        self.mode = mode
        self.align = align

        self._children = children


class VPacker(PackerBase):
    """
    The VPacker has its children packed vertically. It automatically
    adjust the relative positions of children in the drawing time.
    """
    def __init__(self, pad=None, sep=None, width=None, height=None,
                 align="baseline", mode="fixed",
                 children=None):
        """
        Parameters
        ----------
        pad : float, optional
            Boundary pad.

        sep : float, optional
            Spacing between items.

        width : float, optional

        height : float, optional

            width and height of the container box, calculated if
            `None`.

        align : str, optional
            Alignment of boxes.

        mode : str, optional
            Packing mode.

        Notes
        -----
        *pad* and *sep* need to given in points and will be scale with
        the renderer dpi, while *width* and *height* need to be in
        pixels.
        """
        super(VPacker, self).__init__(pad, sep, width, height,
                                      align, mode,
                                      children)

    def get_extent_offsets(self, renderer):
        """
        update offset of childrens and return the extents of the box
        """

        dpicor = renderer.points_to_pixels(1.)
        pad = self.pad * dpicor
        sep = self.sep * dpicor

        if self.width is not None:
            for c in self.get_visible_children():
                if isinstance(c, PackerBase) and c.mode == "expand":
                    c.set_width(self.width)

        whd_list = [c.get_extent(renderer)
                    for c in self.get_visible_children()]
        whd_list = [(w, h, xd, (h - yd)) for w, h, xd, yd in whd_list]

        wd_list = [(w, xd) for w, h, xd, yd in whd_list]
        width, xdescent, xoffsets = _get_aligned_offsets(wd_list,
                                                         self.width,
                                                         self.align)

        pack_list = [(h, yd) for w, h, xd, yd in whd_list]
        height, yoffsets_ = _get_packed_offsets(pack_list, self.height,
                                                sep, self.mode)

        yoffsets = yoffsets_ + [yd for w, h, xd, yd in whd_list]
        ydescent = height - yoffsets[0]
        yoffsets = height - yoffsets

        yoffsets = yoffsets - ydescent

        return width + 2 * pad, height + 2 * pad, \
               xdescent + pad, ydescent + pad, \
               list(zip(xoffsets, yoffsets))


class HPacker(PackerBase):
    """
    The HPacker has its children packed horizontally. It automatically
    adjusts the relative positions of children at draw time.
    """
    def __init__(self, pad=None, sep=None, width=None, height=None,
                 align="baseline", mode="fixed",
                 children=None):
        """
        Parameters
        ----------
        pad : float, optional
            Boundary pad.

        sep : float, optional
            Spacing between items.

        width : float, optional

        height : float, optional
           Width and height of the container box, calculated if
           `None`.

        align : str
           Alignment of boxes.

        mode : str
           Packing mode.

        Notes
        -----
        *pad* and *sep* need to given in points and will be scale with
        the renderer dpi, while *width* and *height* need to be in
        pixels.
        """
        super(HPacker, self).__init__(pad, sep, width, height,
                                      align, mode, children)

    def get_extent_offsets(self, renderer):
        """
        update offset of children and return the extents of the box
        """
        dpicor = renderer.points_to_pixels(1.)
        pad = self.pad * dpicor
        sep = self.sep * dpicor

        whd_list = [c.get_extent(renderer)
                    for c in self.get_visible_children()]

        if not whd_list:
            return 2 * pad, 2 * pad, pad, pad, []

        if self.height is None:
            height_descent = max(h - yd for w, h, xd, yd in whd_list)
            ydescent = max(yd for w, h, xd, yd in whd_list)
            height = height_descent + ydescent
        else:
            height = self.height - 2 * pad  # width w/o pad

        hd_list = [(h, yd) for w, h, xd, yd in whd_list]
        height, ydescent, yoffsets = _get_aligned_offsets(hd_list,
                                                          self.height,
                                                          self.align)

        pack_list = [(w, xd) for w, h, xd, yd in whd_list]

        width, xoffsets_ = _get_packed_offsets(pack_list, self.width,
                                               sep, self.mode)

        xoffsets = xoffsets_ + [xd for w, h, xd, yd in whd_list]

        xdescent = whd_list[0][2]
        xoffsets = xoffsets - xdescent

        return width + 2 * pad, height + 2 * pad, \
               xdescent + pad, ydescent + pad, \
               list(zip(xoffsets, yoffsets))


class PaddedBox(OffsetBox):
    def __init__(self, child, pad=None, draw_frame=False, patch_attrs=None):
        """
        *pad* : boundary pad

        .. note::
          *pad* need to given in points and will be
          scale with the renderer dpi, while *width* and *height*
          need to be in pixels.
        """

        super(PaddedBox, self).__init__()

        self.pad = pad
        self._children = [child]

        self.patch = FancyBboxPatch(
            xy=(0.0, 0.0), width=1., height=1.,
            facecolor='w', edgecolor='k',
            mutation_scale=1,  # self.prop.get_size_in_points(),
            snap=True
            )

        self.patch.set_boxstyle("square", pad=0)

        if patch_attrs is not None:
            self.patch.update(patch_attrs)

        self._drawFrame = draw_frame

    def get_extent_offsets(self, renderer):
        """
        update offset of childrens and return the extents of the box
        """

        dpicor = renderer.points_to_pixels(1.)
        pad = self.pad * dpicor

        w, h, xd, yd = self._children[0].get_extent(renderer)

        return w + 2 * pad, h + 2 * pad, \
               xd + pad, yd + pad, \
               [(0, 0)]

    def draw(self, renderer):
        """
        Update the location of children if necessary and draw them
        to the given *renderer*.
        """

        width, height, xdescent, ydescent, offsets = self.get_extent_offsets(
                                                        renderer)

        px, py = self.get_offset(width, height, xdescent, ydescent, renderer)

        for c, (ox, oy) in zip(self.get_visible_children(), offsets):
            c.set_offset((px + ox, py + oy))

        self.draw_frame(renderer)

        for c in self.get_visible_children():
            c.draw(renderer)

        #bbox_artist(self, renderer, fill=False, props=dict(pad=0.))
        self.stale = False

    def update_frame(self, bbox, fontsize=None):
        self.patch.set_bounds(bbox.x0, bbox.y0,
                              bbox.width, bbox.height)

        if fontsize:
            self.patch.set_mutation_scale(fontsize)
        self.stale = True

    def draw_frame(self, renderer):
        # update the location and size of the legend
        bbox = self.get_window_extent(renderer)
        self.update_frame(bbox)

        if self._drawFrame:
            self.patch.draw(renderer)


class DrawingArea(OffsetBox):
    """
    The DrawingArea can contain any Artist as a child. The DrawingArea
    has a fixed width and height. The position of children relative to
    the parent is fixed. The children can be clipped at the
    boundaries of the parent.
    """

    def __init__(self, width, height, xdescent=0.,
                 ydescent=0., clip=False):
        """
        *width*, *height* : width and height of the container box.
        *xdescent*, *ydescent* : descent of the box in x- and y-direction.
        *clip* : Whether to clip the children
        """

        super(DrawingArea, self).__init__()

        self.width = width
        self.height = height
        self.xdescent = xdescent
        self.ydescent = ydescent
        self._clip_children = clip

        self.offset_transform = mtransforms.Affine2D()
        self.offset_transform.clear()
        self.offset_transform.translate(0, 0)

        self.dpi_transform = mtransforms.Affine2D()

    @property
    def clip_children(self):
        """
        If the children of this DrawingArea should be clipped
        by DrawingArea bounding box.
        """
        return self._clip_children

    @clip_children.setter
    def clip_children(self, val):
        self._clip_children = bool(val)
        self.stale = True

    def get_transform(self):
        """
        Return the :class:`~matplotlib.transforms.Transform` applied
        to the children
        """
        return self.dpi_transform + self.offset_transform

    def set_transform(self, t):
        """
        set_transform is ignored.
        """
        pass

    def set_offset(self, xy):
        """
        set offset of the container.

        Accept : tuple of x,y coordinate in display units.
        """
        self._offset = xy

        self.offset_transform.clear()
        self.offset_transform.translate(xy[0], xy[1])
        self.stale = True

    def get_offset(self):
        """
        return offset of the container.
        """
        return self._offset

    def get_window_extent(self, renderer):
        '''
        get the bounding box in display space.
        '''
        w, h, xd, yd = self.get_extent(renderer)
        ox, oy = self.get_offset()  # w, h, xd, yd)

        return mtransforms.Bbox.from_bounds(ox - xd, oy - yd, w, h)

    def get_extent(self, renderer):
        """
        Return with, height, xdescent, ydescent of box
        """

        dpi_cor = renderer.points_to_pixels(1.)
        return self.width * dpi_cor, self.height * dpi_cor, \
               self.xdescent * dpi_cor, self.ydescent * dpi_cor

    def add_artist(self, a):
        'Add any :class:`~matplotlib.artist.Artist` to the container box'
        self._children.append(a)
        if not a.is_transform_set():
            a.set_transform(self.get_transform())
        if self.axes is not None:
            a.axes = self.axes
        fig = self.figure
        if fig is not None:
            a.set_figure(fig)

    def draw(self, renderer):
        """
        Draw the children
        """

        dpi_cor = renderer.points_to_pixels(1.)
        self.dpi_transform.clear()
        self.dpi_transform.scale(dpi_cor, dpi_cor)

        # At this point the DrawingArea has a transform
        # to the display space so the path created is
        # good for clipping children
        tpath = mtransforms.TransformedPath(
            mpath.Path([[0, 0], [0, self.height],
                        [self.width, self.height],
                        [self.width, 0]]),
            self.get_transform())
        for c in self._children:
            if self._clip_children and not (c.clipbox or c._clippath):
                c.set_clip_path(tpath)
            c.draw(renderer)

        bbox_artist(self, renderer, fill=False, props=dict(pad=0.))
        self.stale = False


class TextArea(OffsetBox):
    """
    The TextArea is contains a single Text instance. The text is
    placed at (0,0) with baseline+left alignment. The width and height
    of the TextArea instance is the width and height of the its child
    text.
    """
    def __init__(self, s,
                 textprops=None,
                 multilinebaseline=None,
                 minimumdescent=True,
                 ):
        """
        Parameters
        ----------
        s : str
            a string to be displayed.

        textprops : `~matplotlib.font_manager.FontProperties`, optional

        multilinebaseline : bool, optional
            If `True`, baseline for multiline text is adjusted so that
            it is (approximatedly) center-aligned with singleline
            text.

        minimumdescent : bool, optional
            If `True`, the box has a minimum descent of "p".
        """
        if textprops is None:
            textprops = {}

        if "va" not in textprops:
            textprops["va"] = "baseline"

        self._text = mtext.Text(0, 0, s, **textprops)

        OffsetBox.__init__(self)

        self._children = [self._text]

        self.offset_transform = mtransforms.Affine2D()
        self.offset_transform.clear()
        self.offset_transform.translate(0, 0)
        self._baseline_transform = mtransforms.Affine2D()
        self._text.set_transform(self.offset_transform +
                                 self._baseline_transform)

        self._multilinebaseline = multilinebaseline
        self._minimumdescent = minimumdescent

    def set_text(self, s):
        "Set the text of this area as a string."
        self._text.set_text(s)
        self.stale = True

    def get_text(self):
        "Returns the string representation of this area's text"
        return self._text.get_text()

    def set_multilinebaseline(self, t):
        """
        Set multilinebaseline .

        If True, baseline for multiline text is
        adjusted so that it is (approximatedly) center-aligned with
        singleline text.
        """
        self._multilinebaseline = t
        self.stale = True

    def get_multilinebaseline(self):
        """
        get multilinebaseline .
        """
        return self._multilinebaseline

    def set_minimumdescent(self, t):
        """
        Set minimumdescent .

        If True, extent of the single line text is adjusted so that
        it has minimum descent of "p"
        """
        self._minimumdescent = t
        self.stale = True

    def get_minimumdescent(self):
        """
        get minimumdescent.
        """
        return self._minimumdescent

    def set_transform(self, t):
        """
        set_transform is ignored.
        """
        pass

    def set_offset(self, xy):
        """
        set offset of the container.

        Accept : tuple of x,y coordinates in display units.
        """
        self._offset = xy

        self.offset_transform.clear()
        self.offset_transform.translate(xy[0], xy[1])
        self.stale = True

    def get_offset(self):
        """
        return offset of the container.
        """
        return self._offset

    def get_window_extent(self, renderer):
        '''
        get the bounding box in display space.
        '''
        w, h, xd, yd = self.get_extent(renderer)
        ox, oy = self.get_offset()  # w, h, xd, yd)
        return mtransforms.Bbox.from_bounds(ox - xd, oy - yd, w, h)

    def get_extent(self, renderer):
        clean_line, ismath = self._text.is_math_text(self._text._text)
        _, h_, d_ = renderer.get_text_width_height_descent(
            "lp", self._text._fontproperties, ismath=False)

        bbox, info, d = self._text._get_layout(renderer)
        w, h = bbox.width, bbox.height

        line = info[-1][0]  # last line

        self._baseline_transform.clear()

        if len(info) > 1 and self._multilinebaseline:
            d_new = 0.5 * h - 0.5 * (h_ - d_)
            self._baseline_transform.translate(0, d - d_new)
            d = d_new

        else:  # single line

            h_d = max(h_ - d_, h - d)

            if self.get_minimumdescent():
                ## to have a minimum descent, #i.e., "l" and "p" have same
                ## descents.
                d = max(d, d_)
            #else:
            #    d = d

            h = h_d + d

        return w, h, 0., d

    def draw(self, renderer):
        """
        Draw the children
        """

        self._text.draw(renderer)

        bbox_artist(self, renderer, fill=False, props=dict(pad=0.))
        self.stale = False


class AuxTransformBox(OffsetBox):
    """
    Offset Box with the aux_transform . Its children will be
    transformed with the aux_transform first then will be
    offseted. The absolute coordinate of the aux_transform is meaning
    as it will be automatically adjust so that the left-lower corner
    of the bounding box of children will be set to (0,0) before the
    offset transform.

    It is similar to drawing area, except that the extent of the box
    is not predetermined but calculated from the window extent of its
    children. Furthermore, the extent of the children will be
    calculated in the transformed coordinate.
    """
    def __init__(self, aux_transform):
        self.aux_transform = aux_transform
        OffsetBox.__init__(self)

        self.offset_transform = mtransforms.Affine2D()
        self.offset_transform.clear()
        self.offset_transform.translate(0, 0)

        # ref_offset_transform is used to make the offset_transform is
        # always reference to the lower-left corner of the bbox of its
        # children.
        self.ref_offset_transform = mtransforms.Affine2D()
        self.ref_offset_transform.clear()

    def add_artist(self, a):
        'Add any :class:`~matplotlib.artist.Artist` to the container box'
        self._children.append(a)
        a.set_transform(self.get_transform())
        self.stale = True

    def get_transform(self):
        """
        Return the :class:`~matplotlib.transforms.Transform` applied
        to the children
        """
        return self.aux_transform + \
               self.ref_offset_transform + \
               self.offset_transform

    def set_transform(self, t):
        """
        set_transform is ignored.
        """
        pass

    def set_offset(self, xy):
        """
        set offset of the container.

        Accept : tuple of x,y coordinate in display units.
        """
        self._offset = xy

        self.offset_transform.clear()
        self.offset_transform.translate(xy[0], xy[1])
        self.stale = True

    def get_offset(self):
        """
        return offset of the container.
        """
        return self._offset

    def get_window_extent(self, renderer):
        '''
        get the bounding box in display space.
        '''
        w, h, xd, yd = self.get_extent(renderer)
        ox, oy = self.get_offset()  # w, h, xd, yd)
        return mtransforms.Bbox.from_bounds(ox - xd, oy - yd, w, h)

    def get_extent(self, renderer):

        # clear the offset transforms
        _off = self.offset_transform.to_values()  # to be restored later
        self.ref_offset_transform.clear()
        self.offset_transform.clear()

        # calculate the extent
        bboxes = [c.get_window_extent(renderer) for c in self._children]
        ub = mtransforms.Bbox.union(bboxes)

        # adjust ref_offset_tansform
        self.ref_offset_transform.translate(-ub.x0, -ub.y0)

        # restor offset transform
        mtx = self.offset_transform.matrix_from_values(*_off)
        self.offset_transform.set_matrix(mtx)

        return ub.width, ub.height, 0., 0.

    def draw(self, renderer):
        """
        Draw the children
        """

        for c in self._children:
            c.draw(renderer)

        bbox_artist(self, renderer, fill=False, props=dict(pad=0.))
        self.stale = False


class AnchoredOffsetbox(OffsetBox):
    """
    An offset box placed according to the legend location
    loc. AnchoredOffsetbox has a single child. When multiple children
    is needed, use other OffsetBox class to enclose them.  By default,
    the offset box is anchored against its parent axes. You may
    explicitly specify the bbox_to_anchor.
    """
    zorder = 5  # zorder of the legend

    # Location codes
    codes = {'upper right': 1,
             'upper left': 2,
             'lower left': 3,
             'lower right': 4,
             'right': 5,
             'center left': 6,
             'center right': 7,
             'lower center': 8,
             'upper center': 9,
             'center': 10,
             }

    def __init__(self, loc,
                 pad=0.4, borderpad=0.5,
                 child=None, prop=None, frameon=True,
                 bbox_to_anchor=None,
                 bbox_transform=None,
                 **kwargs):
        """
        loc is a string or an integer specifying the legend location.
        The valid location codes are::

        'upper right'  : 1,
        'upper left'   : 2,
        'lower left'   : 3,
        'lower right'  : 4,
        'right'        : 5, (same as 'center right', for back-compatibility)
        'center left'  : 6,
        'center right' : 7,
        'lower center' : 8,
        'upper center' : 9,
        'center'       : 10,

        pad : pad around the child for drawing a frame. given in
          fraction of fontsize.

        borderpad : pad between offsetbox frame and the bbox_to_anchor,

        child : OffsetBox instance that will be anchored.

        prop : font property. This is only used as a reference for paddings.

        frameon : draw a frame box if True.

        bbox_to_anchor : bbox to anchor. Use self.axes.bbox if None.

        bbox_transform : with which the bbox_to_anchor will be transformed.

        """
        super(AnchoredOffsetbox, self).__init__(**kwargs)

        self.set_bbox_to_anchor(bbox_to_anchor, bbox_transform)
        self.set_child(child)

        if isinstance(loc, six.string_types):
            try:
                loc = self.codes[loc]
            except KeyError:
                raise ValueError('Unrecognized location "%s". Valid '
                                 'locations are\n\t%s\n'
                                 % (loc, '\n\t'.join(self.codes)))

        self.loc = loc
        self.borderpad = borderpad
        self.pad = pad

        if prop is None:
            self.prop = FontProperties(size=rcParams["legend.fontsize"])
        elif isinstance(prop, dict):
            self.prop = FontProperties(**prop)
            if "size" not in prop:
                self.prop.set_size(rcParams["legend.fontsize"])
        else:
            self.prop = prop

        self.patch = FancyBboxPatch(
            xy=(0.0, 0.0), width=1., height=1.,
            facecolor='w', edgecolor='k',
            mutation_scale=self.prop.get_size_in_points(),
            snap=True
            )
        self.patch.set_boxstyle("square", pad=0)
        self._drawFrame = frameon

    def set_child(self, child):
        "set the child to be anchored"
        self._child = child
        if child is not None:
            child.axes = self.axes
        self.stale = True

    def get_child(self):
        "return the child"
        return self._child

    def get_children(self):
        "return the list of children"
        return [self._child]

    def get_extent(self, renderer):
        """
        return the extent of the artist. The extent of the child
        added with the pad is returned
        """
        w, h, xd, yd = self.get_child().get_extent(renderer)
        fontsize = renderer.points_to_pixels(self.prop.get_size_in_points())
        pad = self.pad * fontsize

        return w + 2 * pad, h + 2 * pad, xd + pad, yd + pad

    def get_bbox_to_anchor(self):
        """
        return the bbox that the legend will be anchored
        """
        if self._bbox_to_anchor is None:
            return self.axes.bbox
        else:
            transform = self._bbox_to_anchor_transform
            if transform is None:
                return self._bbox_to_anchor
            else:
                return TransformedBbox(self._bbox_to_anchor,
                                       transform)

    def set_bbox_to_anchor(self, bbox, transform=None):
        """
        set the bbox that the child will be anchored.

        *bbox* can be a Bbox instance, a list of [left, bottom, width,
        height], or a list of [left, bottom] where the width and
        height will be assumed to be zero. The bbox will be
        transformed to display coordinate by the given transform.
        """
        if bbox is None or isinstance(bbox, BboxBase):
            self._bbox_to_anchor = bbox
        else:
            try:
                l = len(bbox)
            except TypeError:
                raise ValueError("Invalid argument for bbox : %s" % str(bbox))

            if l == 2:
                bbox = [bbox[0], bbox[1], 0, 0]

            self._bbox_to_anchor = Bbox.from_bounds(*bbox)

        self._bbox_to_anchor_transform = transform
        self.stale = True

    def get_window_extent(self, renderer):
        '''
        get the bounding box in display space.
        '''
        self._update_offset_func(renderer)
        w, h, xd, yd = self.get_extent(renderer)
        ox, oy = self.get_offset(w, h, xd, yd, renderer)
        return Bbox.from_bounds(ox - xd, oy - yd, w, h)

    def _update_offset_func(self, renderer, fontsize=None):
        """
        Update the offset func which depends on the dpi of the
        renderer (because of the padding).
        """
        if fontsize is None:
            fontsize = renderer.points_to_pixels(
                            self.prop.get_size_in_points())

        def _offset(w, h, xd, yd, renderer, fontsize=fontsize, self=self):
            bbox = Bbox.from_bounds(0, 0, w, h)
            borderpad = self.borderpad * fontsize
            bbox_to_anchor = self.get_bbox_to_anchor()

            x0, y0 = self._get_anchored_bbox(self.loc,
                                             bbox,
                                             bbox_to_anchor,
                                             borderpad)
            return x0 + xd, y0 + yd

        self.set_offset(_offset)

    def update_frame(self, bbox, fontsize=None):
        self.patch.set_bounds(bbox.x0, bbox.y0,
                              bbox.width, bbox.height)

        if fontsize:
            self.patch.set_mutation_scale(fontsize)

    def draw(self, renderer):
        "draw the artist"

        if not self.get_visible():
            return

        fontsize = renderer.points_to_pixels(self.prop.get_size_in_points())
        self._update_offset_func(renderer, fontsize)

        if self._drawFrame:
            # update the location and size of the legend
            bbox = self.get_window_extent(renderer)
            self.update_frame(bbox, fontsize)
            self.patch.draw(renderer)

        width, height, xdescent, ydescent = self.get_extent(renderer)

        px, py = self.get_offset(width, height, xdescent, ydescent, renderer)

        self.get_child().set_offset((px, py))
        self.get_child().draw(renderer)
        self.stale = False

    def _get_anchored_bbox(self, loc, bbox, parentbbox, borderpad):
        """
        return the position of the bbox anchored at the parentbbox
        with the loc code, with the borderpad.
        """
        assert loc in range(1, 11)  # called only internally

        BEST, UR, UL, LL, LR, R, CL, CR, LC, UC, C = xrange(11)

        anchor_coefs = {UR: "NE",
                        UL: "NW",
                        LL: "SW",
                        LR: "SE",
                        R: "E",
                        CL: "W",
                        CR: "E",
                        LC: "S",
                        UC: "N",
                        C: "C"}

        c = anchor_coefs[loc]

        container = parentbbox.padded(-borderpad)
        anchored_box = bbox.anchored(c, container=container)
        return anchored_box.x0, anchored_box.y0


class AnchoredText(AnchoredOffsetbox):
    """
    AnchoredOffsetbox with Text.
    """

    def __init__(self, s, loc, pad=0.4, borderpad=0.5, prop=None, **kwargs):
        """
        Parameters
        ----------
        s : string
            Text.

        loc : str
            Location code.

        pad : float, optional
            Pad between the text and the frame as fraction of the font
            size.

        borderpad : float, optional
            Pad between the frame and the axes (or *bbox_to_anchor*).

        prop : `matplotlib.font_manager.FontProperties`
            Font properties.

        Notes
        -----
        Other keyword parameters of `AnchoredOffsetbox` are also
        allowed.
        """

        if prop is None:
            prop = {}
        badkwargs = {'ha', 'horizontalalignment', 'va', 'verticalalignment'}
        if badkwargs & set(prop):
            warnings.warn("Mixing horizontalalignment or verticalalignment "
                          "with AnchoredText is not supported.")

        self.txt = TextArea(s, textprops=prop, minimumdescent=False)
        fp = self.txt._text.get_fontproperties()
        super(AnchoredText, self).__init__(
            loc, pad=pad, borderpad=borderpad, child=self.txt, prop=fp,
            **kwargs)


class OffsetImage(OffsetBox):
    def __init__(self, arr,
                 zoom=1,
                 cmap=None,
                 norm=None,
                 interpolation=None,
                 origin=None,
                 filternorm=1,
                 filterrad=4.0,
                 resample=False,
                 dpi_cor=True,
                 **kwargs
                 ):

        OffsetBox.__init__(self)
        self._dpi_cor = dpi_cor

        self.image = BboxImage(bbox=self.get_window_extent,
                               cmap=cmap,
                               norm=norm,
                               interpolation=interpolation,
                               origin=origin,
                               filternorm=filternorm,
                               filterrad=filterrad,
                               resample=resample,
                               **kwargs
                               )

        self._children = [self.image]

        self.set_zoom(zoom)
        self.set_data(arr)

    def set_data(self, arr):
        self._data = np.asarray(arr)
        self.image.set_data(self._data)
        self.stale = True

    def get_data(self):
        return self._data

    def set_zoom(self, zoom):
        self._zoom = zoom
        self.stale = True

    def get_zoom(self):
        return self._zoom

#     def set_axes(self, axes):
#         self.image.set_axes(axes)
#         martist.Artist.set_axes(self, axes)

#     def set_offset(self, xy):
#         """
#         set offset of the container.

#         Accept : tuple of x,y coordinate in display units.
#         """
#         self._offset = xy

#         self.offset_transform.clear()
#         self.offset_transform.translate(xy[0], xy[1])

    def get_offset(self):
        """
        return offset of the container.
        """
        return self._offset

    def get_children(self):
        return [self.image]

    def get_window_extent(self, renderer):
        '''
        get the bounding box in display space.
        '''
        w, h, xd, yd = self.get_extent(renderer)
        ox, oy = self.get_offset()
        return mtransforms.Bbox.from_bounds(ox - xd, oy - yd, w, h)

    def get_extent(self, renderer):
        if self._dpi_cor:  # True, do correction
            dpi_cor = renderer.points_to_pixels(1.)
        else:
            dpi_cor = 1.

        zoom = self.get_zoom()
        data = self.get_data()
        ny, nx = data.shape[:2]
        w, h = dpi_cor * nx * zoom, dpi_cor * ny * zoom

        return w, h, 0, 0

    def draw(self, renderer):
        """
        Draw the children
        """
        self.image.draw(renderer)
        # bbox_artist(self, renderer, fill=False, props=dict(pad=0.))
        self.stale = False


class AnnotationBbox(martist.Artist, _AnnotationBase):
    """
    Annotation-like class, but with offsetbox instead of Text.
    """
    zorder = 3

    def __str__(self):
        return "AnnotationBbox(%g,%g)" % (self.xy[0], self.xy[1])

    @docstring.dedent_interpd
    def __init__(self, offsetbox, xy,
                 xybox=None,
                 xycoords='data',
                 boxcoords=None,
                 frameon=True, pad=0.4,  # BboxPatch
                 annotation_clip=None,
                 box_alignment=(0.5, 0.5),
                 bboxprops=None,
                 arrowprops=None,
                 fontsize=None,
                 **kwargs):
        """
        *offsetbox* : OffsetBox instance

        *xycoords* : same as Annotation but can be a tuple of two
           strings which are interpreted as x and y coordinates.

        *boxcoords* : similar to textcoords as Annotation but can be a
           tuple of two strings which are interpreted as x and y
           coordinates.

        *box_alignment* : a tuple of two floats for a vertical and
           horizontal alignment of the offset box w.r.t. the *boxcoords*.
           The lower-left corner is (0.0) and upper-right corner is (1.1).

        other parameters are identical to that of Annotation.
        """

        martist.Artist.__init__(self, **kwargs)
        _AnnotationBase.__init__(self,
                                 xy,
                                 xycoords=xycoords,
                                 annotation_clip=annotation_clip)

        self.offsetbox = offsetbox

        self.arrowprops = arrowprops

        self.set_fontsize(fontsize)

        if xybox is None:
            self.xybox = xy
        else:
            self.xybox = xybox

        if boxcoords is None:
            self.boxcoords = xycoords
        else:
            self.boxcoords = boxcoords

        if arrowprops is not None:
            self._arrow_relpos = self.arrowprops.pop("relpos", (0.5, 0.5))
            self.arrow_patch = FancyArrowPatch((0, 0), (1, 1),
                                               **self.arrowprops)
        else:
            self._arrow_relpos = None
            self.arrow_patch = None

        #self._fw, self._fh = 0., 0. # for alignment
        self._box_alignment = box_alignment

        # frame
        self.patch = FancyBboxPatch(
            xy=(0.0, 0.0), width=1., height=1.,
            facecolor='w', edgecolor='k',
            mutation_scale=self.prop.get_size_in_points(),
            snap=True
            )
        self.patch.set_boxstyle("square", pad=pad)
        if bboxprops:
            self.patch.set(**bboxprops)
        self._drawFrame = frameon

    @property
    def xyann(self):
        return self.xybox

    @xyann.setter
    def xyann(self, xyann):
        self.xybox = xyann
        self.stale = True

    @property
    def anncoords(self):
        return self.boxcoords

    @anncoords.setter
    def anncoords(self, coords):
        self.boxcoords = coords
        self.stale = True

    def contains(self, event):
        t, tinfo = self.offsetbox.contains(event)
        #if self.arrow_patch is not None:
        #    a,ainfo=self.arrow_patch.contains(event)
        #    t = t or a

        # self.arrow_patch is currently not checked as this can be a line - JJ

        return t, tinfo

    def get_children(self):
        children = [self.offsetbox, self.patch]
        if self.arrow_patch:
            children.append(self.arrow_patch)
        return children

    def set_figure(self, fig):

        if self.arrow_patch is not None:
            self.arrow_patch.set_figure(fig)
        self.offsetbox.set_figure(fig)
        martist.Artist.set_figure(self, fig)

    def set_fontsize(self, s=None):
        """
        set fontsize in points
        """
        if s is None:
            s = rcParams["legend.fontsize"]

        self.prop = FontProperties(size=s)
        self.stale = True

    def get_fontsize(self, s=None):
        """
        return fontsize in points
        """
        return self.prop.get_size_in_points()

    def update_positions(self, renderer):
        """
        Update the pixel positions of the annotated point and the text.
        """
        xy_pixel = self._get_position_xy(renderer)
        self._update_position_xybox(renderer, xy_pixel)

        mutation_scale = renderer.points_to_pixels(self.get_fontsize())
        self.patch.set_mutation_scale(mutation_scale)

        if self.arrow_patch:
            self.arrow_patch.set_mutation_scale(mutation_scale)

    def _update_position_xybox(self, renderer, xy_pixel):
        """
        Update the pixel positions of the annotation text and the arrow
        patch.
        """

        x, y = self.xybox
        if isinstance(self.boxcoords, tuple):
            xcoord, ycoord = self.boxcoords
            x1, y1 = self._get_xy(renderer, x, y, xcoord)
            x2, y2 = self._get_xy(renderer, x, y, ycoord)
            ox0, oy0 = x1, y2
        else:
            ox0, oy0 = self._get_xy(renderer, x, y, self.boxcoords)

        w, h, xd, yd = self.offsetbox.get_extent(renderer)

        _fw, _fh = self._box_alignment
        self.offsetbox.set_offset((ox0 - _fw * w + xd, oy0 - _fh * h + yd))

        # update patch position
        bbox = self.offsetbox.get_window_extent(renderer)
        #self.offsetbox.set_offset((ox0-_fw*w, oy0-_fh*h))
        self.patch.set_bounds(bbox.x0, bbox.y0,
                              bbox.width, bbox.height)

        x, y = xy_pixel

        ox1, oy1 = x, y

        if self.arrowprops:
            x0, y0 = x, y

            d = self.arrowprops.copy()

            # Use FancyArrowPatch if self.arrowprops has "arrowstyle" key.

            # adjust the starting point of the arrow relative to
            # the textbox.
            # TODO : Rotation needs to be accounted.
            relpos = self._arrow_relpos

            ox0 = bbox.x0 + bbox.width * relpos[0]
            oy0 = bbox.y0 + bbox.height * relpos[1]

            # The arrow will be drawn from (ox0, oy0) to (ox1,
            # oy1). It will be first clipped by patchA and patchB.
            # Then it will be shrunk by shrinkA and shrinkB
            # (in points). If patch A is not set, self.bbox_patch
            # is used.

            self.arrow_patch.set_positions((ox0, oy0), (ox1, oy1))
            fs = self.prop.get_size_in_points()
            mutation_scale = d.pop("mutation_scale", fs)
            mutation_scale = renderer.points_to_pixels(mutation_scale)
            self.arrow_patch.set_mutation_scale(mutation_scale)

            patchA = d.pop("patchA", self.patch)
            self.arrow_patch.set_patchA(patchA)

    def draw(self, renderer):
        """
        Draw the :class:`Annotation` object to the given *renderer*.
        """

        if renderer is not None:
            self._renderer = renderer
        if not self.get_visible():
            return

        xy_pixel = self._get_position_xy(renderer)

        if not self._check_xy(renderer, xy_pixel):
            return

        self.update_positions(renderer)

        if self.arrow_patch is not None:
            if self.arrow_patch.figure is None and self.figure is not None:
                self.arrow_patch.figure = self.figure
            self.arrow_patch.draw(renderer)

        if self._drawFrame:
            self.patch.draw(renderer)

        self.offsetbox.draw(renderer)
        self.stale = False


class DraggableBase(object):
    """
    helper code for a draggable artist (legend, offsetbox)
    The derived class must override following two method.

      def save_offset(self):
          pass

      def update_offset(self, dx, dy):
          pass

    *save_offset* is called when the object is picked for dragging and it
    is meant to save reference position of the artist.

    *update_offset* is called during the dragging. dx and dy is the pixel
     offset from the point where the mouse drag started.

    Optionally you may override following two methods.

      def artist_picker(self, artist, evt):
          return self.ref_artist.contains(evt)

      def finalize_offset(self):
          pass

    *artist_picker* is a picker method that will be
     used. *finalize_offset* is called when the mouse is released. In
     current implementation of DraggableLegend and DraggableAnnotation,
     *update_offset* places the artists simply in display
     coordinates. And *finalize_offset* recalculate their position in
     the normalized axes coordinate and set a relavant attribute.

    """
    def __init__(self, ref_artist, use_blit=False):
        self.ref_artist = ref_artist
        self.got_artist = False

        self.canvas = self.ref_artist.figure.canvas
        self._use_blit = use_blit and self.canvas.supports_blit

        c2 = self.canvas.mpl_connect('pick_event', self.on_pick)
        c3 = self.canvas.mpl_connect('button_release_event', self.on_release)

        ref_artist.set_picker(self.artist_picker)
        self.cids = [c2, c3]

    def on_motion(self, evt):
        if self.got_artist:
            dx = evt.x - self.mouse_x
            dy = evt.y - self.mouse_y
            self.update_offset(dx, dy)
            self.canvas.draw()

    def on_motion_blit(self, evt):
        if self.got_artist:
            dx = evt.x - self.mouse_x
            dy = evt.y - self.mouse_y
            self.update_offset(dx, dy)
            self.canvas.restore_region(self.background)
            self.ref_artist.draw(self.ref_artist.figure._cachedRenderer)
            self.canvas.blit(self.ref_artist.figure.bbox)

    def on_pick(self, evt):
        if evt.artist == self.ref_artist:

            self.mouse_x = evt.mouseevent.x
            self.mouse_y = evt.mouseevent.y
            self.got_artist = True

            if self._use_blit:
                self.ref_artist.set_animated(True)
                self.canvas.draw()
                self.background = self.canvas.copy_from_bbox(
                                    self.ref_artist.figure.bbox)
                self.ref_artist.draw(self.ref_artist.figure._cachedRenderer)
                self.canvas.blit(self.ref_artist.figure.bbox)
                self._c1 = self.canvas.mpl_connect('motion_notify_event',
                                                   self.on_motion_blit)
            else:
                self._c1 = self.canvas.mpl_connect('motion_notify_event',
                                                   self.on_motion)
            self.save_offset()

    def on_release(self, event):
        if self.got_artist:
            self.finalize_offset()
            self.got_artist = False
            self.canvas.mpl_disconnect(self._c1)

            if self._use_blit:
                self.ref_artist.set_animated(False)

    def disconnect(self):
        """disconnect the callbacks"""
        for cid in self.cids:
            self.canvas.mpl_disconnect(cid)
        try:
            c1 = self._c1
        except AttributeError:
            pass
        else:
            self.canvas.mpl_disconnect(c1)

    def artist_picker(self, artist, evt):
        return self.ref_artist.contains(evt)

    def save_offset(self):
        pass

    def update_offset(self, dx, dy):
        pass

    def finalize_offset(self):
        pass


class DraggableOffsetBox(DraggableBase):
    def __init__(self, ref_artist, offsetbox, use_blit=False):
        DraggableBase.__init__(self, ref_artist, use_blit=use_blit)
        self.offsetbox = offsetbox

    def save_offset(self):
        offsetbox = self.offsetbox
        renderer = offsetbox.figure._cachedRenderer
        w, h, xd, yd = offsetbox.get_extent(renderer)
        offset = offsetbox.get_offset(w, h, xd, yd, renderer)
        self.offsetbox_x, self.offsetbox_y = offset
        self.offsetbox.set_offset(offset)

    def update_offset(self, dx, dy):
        loc_in_canvas = self.offsetbox_x + dx, self.offsetbox_y + dy
        self.offsetbox.set_offset(loc_in_canvas)

    def get_loc_in_canvas(self):

        offsetbox = self.offsetbox
        renderer = offsetbox.figure._cachedRenderer
        w, h, xd, yd = offsetbox.get_extent(renderer)
        ox, oy = offsetbox._offset
        loc_in_canvas = (ox - xd, oy - yd)

        return loc_in_canvas


class DraggableAnnotation(DraggableBase):
    def __init__(self, annotation, use_blit=False):
        DraggableBase.__init__(self, annotation, use_blit=use_blit)
        self.annotation = annotation

    def save_offset(self):
        ann = self.annotation
        self.ox, self.oy = ann.get_transform().transform(ann.xyann)

    def update_offset(self, dx, dy):
        ann = self.annotation
        ann.xyann = ann.get_transform().inverted().transform(
            (self.ox + dx, self.oy + dy))


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    fig = plt.figure(1)
    fig.clf()
    ax = plt.subplot(121)

    #txt = ax.text(0.5, 0.5, "Test", size=30, ha="center", color="w")
    kwargs = dict()

    a = np.arange(256).reshape(16, 16) / 256.
    myimage = OffsetImage(a,
                          zoom=2,
                          norm=None,
                          origin=None,
                          **kwargs
                          )
    ax.add_artist(myimage)

    myimage.set_offset((100, 100))

    myimage2 = OffsetImage(a,
                           zoom=2,
                           norm=None,
                           origin=None,
                           **kwargs
                           )
    ann = AnnotationBbox(myimage2, (0.5, 0.5),
                         xybox=(30, 30),
                         xycoords='data',
                         boxcoords="offset points",
                         frameon=True, pad=0.4,  # BboxPatch
                         bboxprops=dict(boxstyle="round", fc="y"),
                         fontsize=None,
                         arrowprops=dict(arrowstyle="->"),
                         )

    ax.add_artist(ann)

    plt.draw()
    plt.show()
