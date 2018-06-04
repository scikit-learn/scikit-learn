"""
Classes for including text in a figure.
"""
from __future__ import absolute_import, division, print_function

import six
from six.moves import zip

import contextlib
import logging
import math
import warnings
import weakref

import numpy as np

from . import artist, cbook, docstring, rcParams
from .artist import Artist
from .font_manager import FontProperties
from .lines import Line2D
from .patches import FancyArrowPatch, FancyBboxPatch, Rectangle
from .textpath import TextPath  # Unused, but imported by others.
from .transforms import (
    Affine2D, Bbox, BboxBase, BboxTransformTo, IdentityTransform, Transform)


_log = logging.getLogger(__name__)


def _process_text_args(override, fontdict=None, **kwargs):
    "Return an override dict.  See :func:`~pyplot.text' docstring for info"

    if fontdict is not None:
        override.update(fontdict)

    override.update(kwargs)
    return override


@contextlib.contextmanager
def _wrap_text(textobj):
    """Temporarily inserts newlines to the text if the wrap option is enabled.
    """
    if textobj.get_wrap():
        old_text = textobj.get_text()
        try:
            textobj.set_text(textobj._get_wrapped_text())
            yield textobj
        finally:
            textobj.set_text(old_text)
    else:
        yield textobj


# Extracted from Text's method to serve as a function
def get_rotation(rotation):
    """
    Return the text angle as float. The returned
    angle is between 0 and 360 deg.

    *rotation* may be 'horizontal', 'vertical', or a numeric value in degrees.
    """
    try:
        angle = float(rotation)
    except (ValueError, TypeError):
        isString = isinstance(rotation, six.string_types)
        if ((isString and rotation == 'horizontal') or rotation is None):
            angle = 0.
        elif (isString and rotation == 'vertical'):
            angle = 90.
        else:
            raise ValueError("rotation is {0} expected either 'horizontal'"
                             " 'vertical', numeric value or"
                             "None".format(rotation))

    return angle % 360


def _get_textbox(text, renderer):
    """
    Calculate the bounding box of the text. Unlike
    :meth:`matplotlib.text.Text.get_extents` method, The bbox size of
    the text before the rotation is calculated.
    """
    # TODO : This function may move into the Text class as a method. As a
    # matter of fact, The information from the _get_textbox function
    # should be available during the Text._get_layout() call, which is
    # called within the _get_textbox. So, it would better to move this
    # function as a method with some refactoring of _get_layout method.

    projected_xs = []
    projected_ys = []

    theta = np.deg2rad(text.get_rotation())
    tr = Affine2D().rotate(-theta)

    _, parts, d = text._get_layout(renderer)

    for t, wh, x, y in parts:
        w, h = wh

        xt1, yt1 = tr.transform_point((x, y))
        yt1 -= d
        xt2, yt2 = xt1 + w, yt1 + h

        projected_xs.extend([xt1, xt2])
        projected_ys.extend([yt1, yt2])

    xt_box, yt_box = min(projected_xs), min(projected_ys)
    w_box, h_box = max(projected_xs) - xt_box, max(projected_ys) - yt_box

    tr = Affine2D().rotate(theta)

    x_box, y_box = tr.transform_point((xt_box, yt_box))

    return x_box, y_box, w_box, h_box


class Text(Artist):
    """
    Handle storing and drawing of text in window or data coordinates.
    """
    zorder = 3
    _cached = cbook.maxdict(50)

    def __repr__(self):
        return "Text(%g,%g,%s)" % (self._x, self._y, repr(self._text))

    def __init__(self,
                 x=0, y=0, text='',
                 color=None,           # defaults to rc params
                 verticalalignment='baseline',
                 horizontalalignment='left',
                 multialignment=None,
                 fontproperties=None,  # defaults to FontProperties()
                 rotation=None,
                 linespacing=None,
                 rotation_mode=None,
                 usetex=None,          # defaults to rcParams['text.usetex']
                 wrap=False,
                 **kwargs
                 ):
        """
        Create a :class:`~matplotlib.text.Text` instance at *x*, *y*
        with string *text*.

        Valid kwargs are
        %(Text)s
        """

        Artist.__init__(self)
        self._x, self._y = x, y

        if color is None:
            color = rcParams['text.color']
        if fontproperties is None:
            fontproperties = FontProperties()
        elif isinstance(fontproperties, six.string_types):
            fontproperties = FontProperties(fontproperties)

        self.set_text(text)
        self.set_color(color)
        self.set_usetex(usetex)
        self.set_wrap(wrap)
        self._verticalalignment = verticalalignment
        self._horizontalalignment = horizontalalignment
        self._multialignment = multialignment
        self._rotation = rotation
        self._fontproperties = fontproperties
        self._bbox_patch = None  # a FancyBboxPatch instance
        self._renderer = None
        if linespacing is None:
            linespacing = 1.2   # Maybe use rcParam later.
        self._linespacing = linespacing
        self.set_rotation_mode(rotation_mode)
        self.update(kwargs)

    def update(self, kwargs):
        """
        Update properties from a dictionary.
        """
        # Update bbox last, as it depends on font properties.
        sentinel = object()  # bbox can be None, so use another sentinel.
        bbox = kwargs.pop("bbox", sentinel)
        super(Text, self).update(kwargs)
        if bbox is not sentinel:
            self.set_bbox(bbox)

    def __getstate__(self):
        d = super(Text, self).__getstate__()
        # remove the cached _renderer (if it exists)
        d['_renderer'] = None
        return d

    def contains(self, mouseevent):
        """Test whether the mouse event occurred in the patch.

        In the case of text, a hit is true anywhere in the
        axis-aligned bounding-box containing the text.

        Returns True or False.
        """
        if callable(self._contains):
            return self._contains(self, mouseevent)

        if not self.get_visible() or self._renderer is None:
            return False, {}

        l, b, w, h = self.get_window_extent().bounds
        r, t = l + w, b + h

        x, y = mouseevent.x, mouseevent.y
        inside = (l <= x <= r and b <= y <= t)
        cattr = {}

        # if the text has a surrounding patch, also check containment for it,
        # and merge the results with the results for the text.
        if self._bbox_patch:
            patch_inside, patch_cattr = self._bbox_patch.contains(mouseevent)
            inside = inside or patch_inside
            cattr["bbox_patch"] = patch_cattr

        return inside, cattr

    def _get_xy_display(self):
        'get the (possibly unit converted) transformed x, y in display coords'
        x, y = self.get_unitless_position()
        return self.get_transform().transform_point((x, y))

    def _get_multialignment(self):
        if self._multialignment is not None:
            return self._multialignment
        else:
            return self._horizontalalignment

    def get_rotation(self):
        'return the text angle as float in degrees'
        return get_rotation(self._rotation)  # string_or_number -> number

    def set_rotation_mode(self, m):
        """
        Set text rotation mode.

        .. ACCEPTS: [ None | "default" | "anchor" ]

        Parameters
        ----------
        m : ``None`` or ``"default"`` or ``"anchor"``
            If ``None`` or ``"default"``, the text will be first rotated, then
            aligned according to their horizontal and vertical alignments.  If
            ``"anchor"``, then alignment occurs before rotation.
        """
        if m is None or m in ["anchor", "default"]:
            self._rotation_mode = m
        else:
            raise ValueError("Unknown rotation_mode : %s" % repr(m))
        self.stale = True

    def get_rotation_mode(self):
        "get text rotation mode"
        return self._rotation_mode

    def update_from(self, other):
        'Copy properties from other to self'
        Artist.update_from(self, other)
        self._color = other._color
        self._multialignment = other._multialignment
        self._verticalalignment = other._verticalalignment
        self._horizontalalignment = other._horizontalalignment
        self._fontproperties = other._fontproperties.copy()
        self._rotation = other._rotation
        self._picker = other._picker
        self._linespacing = other._linespacing
        self.stale = True

    def _get_layout(self, renderer):
        """
        return the extent (bbox) of the text together with
        multiple-alignment information. Note that it returns an extent
        of a rotated text when necessary.
        """
        key = self.get_prop_tup(renderer=renderer)
        if key in self._cached:
            return self._cached[key]

        horizLayout = []

        thisx, thisy = 0.0, 0.0
        xmin, ymin = 0.0, 0.0
        width, height = 0.0, 0.0
        lines = self.get_text().split('\n')

        whs = np.zeros((len(lines), 2))
        horizLayout = np.zeros((len(lines), 4))

        # Find full vertical extent of font,
        # including ascenders and descenders:
        tmp, lp_h, lp_bl = renderer.get_text_width_height_descent('lp',
                                                         self._fontproperties,
                                                         ismath=False)
        offsety = (lp_h - lp_bl) * self._linespacing

        baseline = 0
        for i, line in enumerate(lines):
            clean_line, ismath = self.is_math_text(line, self.get_usetex())
            if clean_line:
                w, h, d = renderer.get_text_width_height_descent(clean_line,
                                                        self._fontproperties,
                                                        ismath=ismath)
            else:
                w, h, d = 0, 0, 0

            # For multiline text, increase the line spacing when the
            # text net-height(excluding baseline) is larger than that
            # of a "l" (e.g., use of superscripts), which seems
            # what TeX does.
            h = max(h, lp_h)
            d = max(d, lp_bl)

            whs[i] = w, h

            baseline = (h - d) - thisy
            thisy -= max(offsety, (h - d) * self._linespacing)
            horizLayout[i] = thisx, thisy, w, h
            thisy -= d
            width = max(width, w)
            descent = d

        ymin = horizLayout[-1][1]
        ymax = horizLayout[0][1] + horizLayout[0][3]
        height = ymax - ymin
        xmax = xmin + width

        # get the rotation matrix
        M = Affine2D().rotate_deg(self.get_rotation())

        offsetLayout = np.zeros((len(lines), 2))
        offsetLayout[:] = horizLayout[:, 0:2]
        # now offset the individual text lines within the box
        if len(lines) > 1:  # do the multiline aligment
            malign = self._get_multialignment()
            if malign == 'center':
                offsetLayout[:, 0] += width / 2.0 - horizLayout[:, 2] / 2.0
            elif malign == 'right':
                offsetLayout[:, 0] += width - horizLayout[:, 2]

        # the corners of the unrotated bounding box
        cornersHoriz = np.array(
            [(xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin)], float)
        cornersHoriz[:, 1] -= descent

        # now rotate the bbox
        cornersRotated = M.transform(cornersHoriz)

        txs = cornersRotated[:, 0]
        tys = cornersRotated[:, 1]

        # compute the bounds of the rotated box
        xmin, xmax = txs.min(), txs.max()
        ymin, ymax = tys.min(), tys.max()
        width = xmax - xmin
        height = ymax - ymin

        # Now move the box to the target position offset the display
        # bbox by alignment
        halign = self._horizontalalignment
        valign = self._verticalalignment

        rotation_mode = self.get_rotation_mode()
        if rotation_mode != "anchor":
            # compute the text location in display coords and the offsets
            # necessary to align the bbox with that location
            if halign == 'center':
                offsetx = (xmin + width / 2.0)
            elif halign == 'right':
                offsetx = (xmin + width)
            else:
                offsetx = xmin

            if valign == 'center':
                offsety = (ymin + height / 2.0)
            elif valign == 'top':
                offsety = (ymin + height)
            elif valign == 'baseline':
                offsety = (ymin + height) - baseline
            elif valign == 'center_baseline':
                offsety = ymin + height - baseline / 2.0
            else:
                offsety = ymin
        else:
            xmin1, ymin1 = cornersHoriz[0]
            xmax1, ymax1 = cornersHoriz[2]

            if halign == 'center':
                offsetx = (xmin1 + xmax1) / 2.0
            elif halign == 'right':
                offsetx = xmax1
            else:
                offsetx = xmin1

            if valign == 'center':
                offsety = (ymin1 + ymax1) / 2.0
            elif valign == 'top':
                offsety = ymax1
            elif valign == 'baseline':
                offsety = ymax1 - baseline
            elif valign == 'center_baseline':
                offsety = (ymin1 + ymax1 - baseline) / 2.0
            else:
                offsety = ymin1

            offsetx, offsety = M.transform_point((offsetx, offsety))

        xmin -= offsetx
        ymin -= offsety

        bbox = Bbox.from_bounds(xmin, ymin, width, height)

        # now rotate the positions around the first x,y position
        xys = M.transform(offsetLayout)
        xys -= (offsetx, offsety)

        xs, ys = xys[:, 0], xys[:, 1]

        ret = bbox, list(zip(lines, whs, xs, ys)), descent
        self._cached[key] = ret
        return ret

    def set_bbox(self, rectprops):
        """
        Draw a bounding box around self.  rectprops are any settable
        properties for a FancyBboxPatch, e.g., facecolor='red', alpha=0.5.

          t.set_bbox(dict(facecolor='red', alpha=0.5))

        The default boxstyle is 'square'. The mutation
        scale of the FancyBboxPatch is set to the fontsize.

        ACCEPTS: FancyBboxPatch prop dict
        """

        if rectprops is not None:
            props = rectprops.copy()
            boxstyle = props.pop("boxstyle", None)
            pad = props.pop("pad", None)
            if boxstyle is None:
                boxstyle = "square"
                if pad is None:
                    pad = 4  # points
                pad /= self.get_size()  # to fraction of font size
            else:
                if pad is None:
                    pad = 0.3

            # boxstyle could be a callable or a string
            if (isinstance(boxstyle, six.string_types)
                    and "pad" not in boxstyle):
                boxstyle += ",pad=%0.2f" % pad

            bbox_transmuter = props.pop("bbox_transmuter", None)

            self._bbox_patch = FancyBboxPatch(
                                    (0., 0.),
                                    1., 1.,
                                    boxstyle=boxstyle,
                                    bbox_transmuter=bbox_transmuter,
                                    transform=IdentityTransform(),
                                    **props)
        else:
            self._bbox_patch = None

        self._update_clip_properties()

    def get_bbox_patch(self):
        """
        Return the bbox Patch object. Returns None if the
        FancyBboxPatch is not made.
        """
        return self._bbox_patch

    def update_bbox_position_size(self, renderer):
        """
        Update the location and the size of the bbox. This method
        should be used when the position and size of the bbox needs to
        be updated before actually drawing the bbox.
        """

        if self._bbox_patch:

            trans = self.get_transform()

            # don't use self.get_unitless_position here, which refers to text
            # position in Text, and dash position in TextWithDash:
            posx = float(self.convert_xunits(self._x))
            posy = float(self.convert_yunits(self._y))

            posx, posy = trans.transform_point((posx, posy))

            x_box, y_box, w_box, h_box = _get_textbox(self, renderer)
            self._bbox_patch.set_bounds(0., 0., w_box, h_box)
            theta = np.deg2rad(self.get_rotation())
            tr = Affine2D().rotate(theta)
            tr = tr.translate(posx + x_box, posy + y_box)
            self._bbox_patch.set_transform(tr)
            fontsize_in_pixel = renderer.points_to_pixels(self.get_size())
            self._bbox_patch.set_mutation_scale(fontsize_in_pixel)

    def _draw_bbox(self, renderer, posx, posy):

        """ Update the location and the size of the bbox
        (FancyBboxPatch), and draw
        """

        x_box, y_box, w_box, h_box = _get_textbox(self, renderer)
        self._bbox_patch.set_bounds(0., 0., w_box, h_box)
        theta = np.deg2rad(self.get_rotation())
        tr = Affine2D().rotate(theta)
        tr = tr.translate(posx + x_box, posy + y_box)
        self._bbox_patch.set_transform(tr)
        fontsize_in_pixel = renderer.points_to_pixels(self.get_size())
        self._bbox_patch.set_mutation_scale(fontsize_in_pixel)
        self._bbox_patch.draw(renderer)

    def _update_clip_properties(self):
        clipprops = dict(clip_box=self.clipbox,
                         clip_path=self._clippath,
                         clip_on=self._clipon)

        if self._bbox_patch:
            bbox = self._bbox_patch.update(clipprops)

    def set_clip_box(self, clipbox):
        """
        Set the artist's clip :class:`~matplotlib.transforms.Bbox`.

        ACCEPTS: a :class:`matplotlib.transforms.Bbox` instance
        """
        super(Text, self).set_clip_box(clipbox)
        self._update_clip_properties()

    def set_clip_path(self, path, transform=None):
        """
        Set the artist's clip path, which may be:

          * a :class:`~matplotlib.patches.Patch` (or subclass) instance

          * a :class:`~matplotlib.path.Path` instance, in which case
             an optional :class:`~matplotlib.transforms.Transform`
             instance may be provided, which will be applied to the
             path before using it for clipping.

          * *None*, to remove the clipping path

        For efficiency, if the path happens to be an axis-aligned
        rectangle, this method will set the clipping box to the
        corresponding rectangle and set the clipping path to *None*.

        ACCEPTS: [ (:class:`~matplotlib.path.Path`,
        :class:`~matplotlib.transforms.Transform`) |
        :class:`~matplotlib.patches.Patch` | None ]
        """
        super(Text, self).set_clip_path(path, transform)
        self._update_clip_properties()

    def set_clip_on(self, b):
        """
        Set whether artist uses clipping.

        When False, artists will be visible outside of the axes, which can lead
        to unexpected results.

        Parameters
        ----------
        b : bool
            .. ACCEPTS: bool
        """
        super(Text, self).set_clip_on(b)
        self._update_clip_properties()

    def get_wrap(self):
        """Returns the wrapping state for the text."""
        return self._wrap

    def set_wrap(self, wrap):
        """Sets the wrapping state for the text.

        Parameters
        ----------
        wrap : bool
            .. ACCEPTS: bool
        """
        self._wrap = wrap

    def _get_wrap_line_width(self):
        """
        Returns the maximum line width for wrapping text based on the
        current orientation.
        """
        x0, y0 = self.get_transform().transform(self.get_position())
        figure_box = self.get_figure().get_window_extent()

        # Calculate available width based on text alignment
        alignment = self.get_horizontalalignment()
        self.set_rotation_mode('anchor')
        rotation = self.get_rotation()

        left = self._get_dist_to_box(rotation, x0, y0, figure_box)
        right = self._get_dist_to_box(
            (180 + rotation) % 360,
            x0,
            y0,
            figure_box)

        if alignment == 'left':
            line_width = left
        elif alignment == 'right':
            line_width = right
        else:
            line_width = 2 * min(left, right)

        return line_width

    def _get_dist_to_box(self, rotation, x0, y0, figure_box):
        """
        Returns the distance from the given points, to the boundaries
        of a rotated box in pixels.
        """
        if rotation > 270:
            quad = rotation - 270
            h1 = y0 / math.cos(math.radians(quad))
            h2 = (figure_box.x1 - x0) / math.cos(math.radians(90 - quad))
        elif rotation > 180:
            quad = rotation - 180
            h1 = x0 / math.cos(math.radians(quad))
            h2 = y0 / math.cos(math.radians(90 - quad))
        elif rotation > 90:
            quad = rotation - 90
            h1 = (figure_box.y1 - y0) / math.cos(math.radians(quad))
            h2 = x0 / math.cos(math.radians(90 - quad))
        else:
            h1 = (figure_box.x1 - x0) / math.cos(math.radians(rotation))
            h2 = (figure_box.y1 - y0) / math.cos(math.radians(90 - rotation))

        return min(h1, h2)

    def _get_rendered_text_width(self, text):
        """
        Returns the width of a given text string, in pixels.
        """
        w, h, d = self._renderer.get_text_width_height_descent(
            text,
            self.get_fontproperties(),
            False)
        return math.ceil(w)

    def _get_wrapped_text(self):
        """
        Return a copy of the text with new lines added, so that
        the text is wrapped relative to the parent figure.
        """
        # Not fit to handle breaking up latex syntax correctly, so
        # ignore latex for now.
        if self.get_usetex():
            return self.get_text()

        # Build the line incrementally, for a more accurate measure of length
        line_width = self._get_wrap_line_width()
        wrapped_str = ""
        line = ""

        for word in self.get_text().split(' '):
            # New lines in the user's test need to force a split, so that it's
            # not using the longest current line width in the line being built
            sub_words = word.split('\n')
            for i in range(len(sub_words)):
                current_width = self._get_rendered_text_width(
                    line + ' ' + sub_words[i])

                # Split long lines, and each newline found in the current word
                if current_width > line_width or i > 0:
                    wrapped_str += line + '\n'
                    line = ""

                if line == "":
                    line = sub_words[i]
                else:
                    line += ' ' + sub_words[i]

        return wrapped_str + line

    @artist.allow_rasterization
    def draw(self, renderer):
        """
        Draws the :class:`Text` object to the given *renderer*.
        """
        if renderer is not None:
            self._renderer = renderer
        if not self.get_visible():
            return
        if self.get_text() == '':
            return

        renderer.open_group('text', self.get_gid())

        with _wrap_text(self) as textobj:
            bbox, info, descent = textobj._get_layout(renderer)
            trans = textobj.get_transform()

            # don't use textobj.get_position here, which refers to text
            # position in Text, and dash position in TextWithDash:
            posx = float(textobj.convert_xunits(textobj._x))
            posy = float(textobj.convert_yunits(textobj._y))
            posx, posy = trans.transform_point((posx, posy))
            if not np.isfinite(posx) or not np.isfinite(posy):
                _log.warning("posx and posy should be finite values")
                return
            canvasw, canvash = renderer.get_canvas_width_height()

            # draw the FancyBboxPatch
            if textobj._bbox_patch:
                textobj._draw_bbox(renderer, posx, posy)

            gc = renderer.new_gc()
            gc.set_foreground(textobj.get_color())
            gc.set_alpha(textobj.get_alpha())
            gc.set_url(textobj._url)
            textobj._set_gc_clip(gc)

            angle = textobj.get_rotation()

            for line, wh, x, y in info:

                mtext = textobj if len(info) == 1 else None
                x = x + posx
                y = y + posy
                if renderer.flipy():
                    y = canvash - y
                clean_line, ismath = textobj.is_math_text(line,
                                                          self.get_usetex())

                if textobj.get_path_effects():
                    from matplotlib.patheffects import PathEffectRenderer
                    textrenderer = PathEffectRenderer(
                                        textobj.get_path_effects(), renderer)
                else:
                    textrenderer = renderer

                if textobj.get_usetex():
                    textrenderer.draw_tex(gc, x, y, clean_line,
                                          textobj._fontproperties, angle,
                                          mtext=mtext)
                else:
                    textrenderer.draw_text(gc, x, y, clean_line,
                                           textobj._fontproperties, angle,
                                           ismath=ismath, mtext=mtext)

        gc.restore()
        renderer.close_group('text')
        self.stale = False

    def get_color(self):
        "Return the color of the text"
        return self._color

    def get_fontproperties(self):
        "Return the :class:`~font_manager.FontProperties` object"
        return self._fontproperties

    def get_font_properties(self):
        'alias for get_fontproperties'
        return self.get_fontproperties()

    def get_family(self):
        "Return the list of font families used for font lookup"
        return self._fontproperties.get_family()

    def get_fontfamily(self):
        'alias for get_family'
        return self.get_family()

    def get_name(self):
        "Return the font name as string"
        return self._fontproperties.get_name()

    def get_style(self):
        "Return the font style as string"
        return self._fontproperties.get_style()

    def get_size(self):
        "Return the font size as integer"
        return self._fontproperties.get_size_in_points()

    def get_variant(self):
        "Return the font variant as a string"
        return self._fontproperties.get_variant()

    def get_fontvariant(self):
        'alias for get_variant'
        return self.get_variant()

    def get_weight(self):
        "Get the font weight as string or number"
        return self._fontproperties.get_weight()

    def get_fontname(self):
        'alias for get_name'
        return self.get_name()

    def get_fontstyle(self):
        'alias for get_style'
        return self.get_style()

    def get_fontsize(self):
        'alias for get_size'
        return self.get_size()

    def get_fontweight(self):
        'alias for get_weight'
        return self.get_weight()

    def get_stretch(self):
        'Get the font stretch as a string or number'
        return self._fontproperties.get_stretch()

    def get_fontstretch(self):
        'alias for get_stretch'
        return self.get_stretch()

    def get_ha(self):
        'alias for get_horizontalalignment'
        return self.get_horizontalalignment()

    def get_horizontalalignment(self):
        """
        Return the horizontal alignment as string.  Will be one of
        'left', 'center' or 'right'.
        """
        return self._horizontalalignment

    def get_unitless_position(self):
        "Return the unitless position of the text as a tuple (*x*, *y*)"
        # This will get the position with all unit information stripped away.
        # This is here for convenience since it is done in several locations.
        x = float(self.convert_xunits(self._x))
        y = float(self.convert_yunits(self._y))
        return x, y

    def get_position(self):
        "Return the position of the text as a tuple (*x*, *y*)"
        # This should return the same data (possible unitized) as was
        # specified with 'set_x' and 'set_y'.
        return self._x, self._y

    def get_prop_tup(self, renderer=None):
        """
        Return a hashable tuple of properties.

        Not intended to be human readable, but useful for backends who
        want to cache derived information about text (e.g., layouts) and
        need to know if the text has changed.
        """
        x, y = self.get_unitless_position()
        renderer = renderer or self._renderer
        return (x, y, self.get_text(), self._color,
                self._verticalalignment, self._horizontalalignment,
                hash(self._fontproperties),
                self._rotation, self._rotation_mode,
                self.figure.dpi, weakref.ref(renderer),
                self._linespacing
                )

    def get_text(self):
        "Get the text as string"
        return self._text

    def get_va(self):
        'alias for :meth:`getverticalalignment`'
        return self.get_verticalalignment()

    def get_verticalalignment(self):
        """
        Return the vertical alignment as string.  Will be one of
        'top', 'center', 'bottom' or 'baseline'.
        """
        return self._verticalalignment

    def get_window_extent(self, renderer=None, dpi=None):
        '''
        Return a :class:`~matplotlib.transforms.Bbox` object bounding
        the text, in display units.

        In addition to being used internally, this is useful for
        specifying clickable regions in a png file on a web page.

        *renderer* defaults to the _renderer attribute of the text
        object.  This is not assigned until the first execution of
        :meth:`draw`, so you must use this kwarg if you want
        to call :meth:`get_window_extent` prior to the first
        :meth:`draw`.  For getting web page regions, it is
        simpler to call the method after saving the figure.

        *dpi* defaults to self.figure.dpi; the renderer dpi is
        irrelevant.  For the web application, if figure.dpi is not
        the value used when saving the figure, then the value that
        was used must be specified as the *dpi* argument.
        '''
        #return _unit_box
        if not self.get_visible():
            return Bbox.unit()
        if dpi is not None:
            dpi_orig = self.figure.dpi
            self.figure.dpi = dpi
        if self.get_text() == '':
            tx, ty = self._get_xy_display()
            return Bbox.from_bounds(tx, ty, 0, 0)

        if renderer is not None:
            self._renderer = renderer
        if self._renderer is None:
            raise RuntimeError('Cannot get window extent w/o renderer')

        bbox, info, descent = self._get_layout(self._renderer)
        x, y = self.get_unitless_position()
        x, y = self.get_transform().transform_point((x, y))
        bbox = bbox.translated(x, y)
        if dpi is not None:
            self.figure.dpi = dpi_orig
        return bbox

    def set_backgroundcolor(self, color):
        """
        Set the background color of the text by updating the bbox.

        .. seealso::

            :meth:`set_bbox`
               To change the position of the bounding box.

        ACCEPTS: any matplotlib color
        """
        if self._bbox_patch is None:
            self.set_bbox(dict(facecolor=color, edgecolor=color))
        else:
            self._bbox_patch.update(dict(facecolor=color))

        self._update_clip_properties()
        self.stale = True

    def set_color(self, color):
        """
        Set the foreground color of the text

        ACCEPTS: any matplotlib color
        """
        # Make sure it is hashable, or get_prop_tup will fail.
        try:
            hash(color)
        except TypeError:
            color = tuple(color)
        self._color = color
        self.stale = True

    def set_ha(self, align):
        'alias for set_horizontalalignment'
        self.set_horizontalalignment(align)

    def set_horizontalalignment(self, align):
        """
        Set the horizontal alignment to one of

        ACCEPTS: [ 'center' | 'right' | 'left' ]
        """
        legal = ('center', 'right', 'left')
        if align not in legal:
            raise ValueError('Horizontal alignment must be one of %s' %
                             str(legal))
        self._horizontalalignment = align
        self.stale = True

    def set_ma(self, align):
        'alias for set_multialignment'
        self.set_multialignment(align)

    def set_multialignment(self, align):
        """
        Set the alignment for multiple lines layout.  The layout of the
        bounding box of all the lines is determined bu the horizontalalignment
        and verticalalignment properties, but the multiline text within that
        box can be

        ACCEPTS: ['left' | 'right' | 'center' ]
        """
        legal = ('center', 'right', 'left')
        if align not in legal:
            raise ValueError('Horizontal alignment must be one of %s' %
                             str(legal))
        self._multialignment = align
        self.stale = True

    def set_linespacing(self, spacing):
        """
        Set the line spacing as a multiple of the font size.
        Default is 1.2.

        ACCEPTS: float (multiple of font size)
        """
        self._linespacing = spacing
        self.stale = True

    def set_family(self, fontname):
        """
        Set the font family.  May be either a single string, or a list
        of strings in decreasing priority.  Each string may be either
        a real font name or a generic font class name.  If the latter,
        the specific font names will be looked up in the
        :file:`matplotlibrc` file.

        ACCEPTS: [FONTNAME | 'serif' | 'sans-serif' | 'cursive' | 'fantasy' |
                  'monospace' ]
        """
        self._fontproperties.set_family(fontname)
        self.stale = True

    def set_variant(self, variant):
        """
        Set the font variant, either 'normal' or 'small-caps'.

        ACCEPTS: [ 'normal' | 'small-caps' ]
        """
        self._fontproperties.set_variant(variant)
        self.stale = True

    def set_fontvariant(self, variant):
        'alias for set_variant'
        return self.set_variant(variant)

    def set_name(self, fontname):
        """alias for set_family"""
        return self.set_family(fontname)

    def set_fontname(self, fontname):
        """alias for set_family"""
        self.set_family(fontname)

    def set_style(self, fontstyle):
        """
        Set the font style.

        ACCEPTS: [ 'normal' | 'italic' | 'oblique']
        """
        self._fontproperties.set_style(fontstyle)
        self.stale = True

    def set_fontstyle(self, fontstyle):
        'alias for set_style'
        return self.set_style(fontstyle)

    def set_size(self, fontsize):
        """
        Set the font size.  May be either a size string, relative to
        the default font size, or an absolute font size in points.

        ACCEPTS: [size in points | 'xx-small' | 'x-small' | 'small' |
                  'medium' | 'large' | 'x-large' | 'xx-large' ]
        """
        self._fontproperties.set_size(fontsize)
        self.stale = True

    def set_fontsize(self, fontsize):
        'alias for set_size'
        return self.set_size(fontsize)

    def set_weight(self, weight):
        """
        Set the font weight.

        ACCEPTS: [a numeric value in range 0-1000 | 'ultralight' | 'light' |
                  'normal' | 'regular' | 'book' | 'medium' | 'roman' |
                  'semibold' | 'demibold' | 'demi' | 'bold' | 'heavy' |
                  'extra bold' | 'black' ]
        """
        self._fontproperties.set_weight(weight)
        self.stale = True

    def set_fontweight(self, weight):
        'alias for set_weight'
        return self.set_weight(weight)

    def set_stretch(self, stretch):
        """
        Set the font stretch (horizontal condensation or expansion).

        ACCEPTS: [a numeric value in range 0-1000 | 'ultra-condensed' |
                  'extra-condensed' | 'condensed' | 'semi-condensed' |
                  'normal' | 'semi-expanded' | 'expanded' | 'extra-expanded' |
                  'ultra-expanded' ]
        """
        self._fontproperties.set_stretch(stretch)
        self.stale = True

    def set_fontstretch(self, stretch):
        'alias for set_stretch'
        return self.set_stretch(stretch)

    def set_position(self, xy):
        """
        Set the (*x*, *y*) position of the text

        ACCEPTS: (x,y)
        """
        self.set_x(xy[0])
        self.set_y(xy[1])

    def set_x(self, x):
        """
        Set the *x* position of the text

        ACCEPTS: float
        """
        self._x = x
        self.stale = True

    def set_y(self, y):
        """
        Set the *y* position of the text

        ACCEPTS: float
        """
        self._y = y
        self.stale = True

    def set_rotation(self, s):
        """
        Set the rotation of the text

        ACCEPTS: [ angle in degrees | 'vertical' | 'horizontal' ]
        """
        self._rotation = s
        self.stale = True

    def set_va(self, align):
        'alias for set_verticalalignment'
        self.set_verticalalignment(align)

    def set_verticalalignment(self, align):
        """
        Set the vertical alignment

        ACCEPTS: [ 'center' | 'top' | 'bottom' | 'baseline' ]
        """
        legal = ('top', 'bottom', 'center', 'baseline')
        if align not in legal:
            raise ValueError('Vertical alignment must be one of %s' %
                             str(legal))

        self._verticalalignment = align
        self.stale = True

    def set_text(self, s):
        """
        Set the text string *s*

        It may contain newlines (``\\n``) or math in LaTeX syntax.

        ACCEPTS: string or anything printable with '%s' conversion.
        """
        self._text = '%s' % (s,)
        self.stale = True

    @staticmethod
    def is_math_text(s, usetex=None):
        """
        Returns a cleaned string and a boolean flag.
        The flag indicates if the given string *s* contains any mathtext,
        determined by counting unescaped dollar signs. If no mathtext
        is present, the cleaned string has its dollar signs unescaped.
        If usetex is on, the flag always has the value "TeX".
        """
        # Did we find an even number of non-escaped dollar signs?
        # If so, treat is as math text.
        if usetex is None:
            usetex = rcParams['text.usetex']
        if usetex:
            if s == ' ':
                s = r'\ '
            return s, 'TeX'

        if cbook.is_math_text(s):
            return s, True
        else:
            return s.replace(r'\$', '$'), False

    def set_fontproperties(self, fp):
        """
        Set the font properties that control the text.  *fp* must be a
        :class:`matplotlib.font_manager.FontProperties` object.

        ACCEPTS: a :class:`matplotlib.font_manager.FontProperties` instance
        """
        if isinstance(fp, six.string_types):
            fp = FontProperties(fp)
        self._fontproperties = fp.copy()
        self.stale = True

    def set_font_properties(self, fp):
        'alias for set_fontproperties'
        self.set_fontproperties(fp)

    def set_usetex(self, usetex):
        """
        Parameters
        ----------
        usetex : bool or None
            Whether to render using TeX, ``None`` means to use
            :rc:`text.usetex`.

            .. ACCEPTS: bool or None
        """
        if usetex is None:
            self._usetex = rcParams['text.usetex']
        else:
            self._usetex = bool(usetex)
        self.stale = True

    def get_usetex(self):
        """
        Return whether this `Text` object uses TeX for rendering.

        If the user has not manually set this value, it defaults to
        :rc:`text.usetex`.
        """
        if self._usetex is None:
            return rcParams['text.usetex']
        else:
            return self._usetex

docstring.interpd.update(Text=artist.kwdoc(Text))
docstring.dedent_interpd(Text.__init__)


class TextWithDash(Text):
    """
    This is basically a :class:`~matplotlib.text.Text` with a dash
    (drawn with a :class:`~matplotlib.lines.Line2D`) before/after
    it. It is intended to be a drop-in replacement for
    :class:`~matplotlib.text.Text`, and should behave identically to
    it when *dashlength* = 0.0.

    The dash always comes between the point specified by
    :meth:`~matplotlib.text.Text.set_position` and the text. When a
    dash exists, the text alignment arguments (*horizontalalignment*,
    *verticalalignment*) are ignored.

    *dashlength* is the length of the dash in canvas units.
    (default = 0.0).

    *dashdirection* is one of 0 or 1, where 0 draws the dash after the
    text and 1 before.  (default = 0).

    *dashrotation* specifies the rotation of the dash, and should
    generally stay *None*. In this case
    :meth:`~matplotlib.text.TextWithDash.get_dashrotation` returns
    :meth:`~matplotlib.text.Text.get_rotation`.  (i.e., the dash takes
    its rotation from the text's rotation). Because the text center is
    projected onto the dash, major deviations in the rotation cause
    what may be considered visually unappealing results.
    (default = *None*)

    *dashpad* is a padding length to add (or subtract) space
    between the text and the dash, in canvas units.
    (default = 3)

    *dashpush* "pushes" the dash and text away from the point
    specified by :meth:`~matplotlib.text.Text.set_position` by the
    amount in canvas units.  (default = 0)

    .. note::

        The alignment of the two objects is based on the bounding box
        of the :class:`~matplotlib.text.Text`, as obtained by
        :meth:`~matplotlib.artist.Artist.get_window_extent`.  This, in
        turn, appears to depend on the font metrics as given by the
        rendering backend. Hence the quality of the "centering" of the
        label text with respect to the dash varies depending on the
        backend used.

    .. note::

        I'm not sure that I got the
        :meth:`~matplotlib.text.TextWithDash.get_window_extent` right,
        or whether that's sufficient for providing the object bounding
        box.

    """
    __name__ = 'textwithdash'

    def __str__(self):
        return "TextWithDash(%g,%g,%s)" % (self._x, self._y, repr(self._text))

    def __init__(self,
                 x=0, y=0, text='',
                 color=None,          # defaults to rc params
                 verticalalignment='center',
                 horizontalalignment='center',
                 multialignment=None,
                 fontproperties=None,  # defaults to FontProperties()
                 rotation=None,
                 linespacing=None,
                 dashlength=0.0,
                 dashdirection=0,
                 dashrotation=None,
                 dashpad=3,
                 dashpush=0,
                 ):

        Text.__init__(self, x=x, y=y, text=text, color=color,
                      verticalalignment=verticalalignment,
                      horizontalalignment=horizontalalignment,
                      multialignment=multialignment,
                      fontproperties=fontproperties,
                      rotation=rotation,
                      linespacing=linespacing)

        # The position (x,y) values for text and dashline
        # are bogus as given in the instantiation; they will
        # be set correctly by update_coords() in draw()

        self.dashline = Line2D(xdata=(x, x),
                               ydata=(y, y),
                               color='k',
                               linestyle='-')

        self._dashx = float(x)
        self._dashy = float(y)
        self._dashlength = dashlength
        self._dashdirection = dashdirection
        self._dashrotation = dashrotation
        self._dashpad = dashpad
        self._dashpush = dashpush

        #self.set_bbox(dict(pad=0))

    def get_unitless_position(self):
        "Return the unitless position of the text as a tuple (*x*, *y*)"
        # This will get the position with all unit information stripped away.
        # This is here for convenience since it is done in several locations.
        x = float(self.convert_xunits(self._dashx))
        y = float(self.convert_yunits(self._dashy))
        return x, y

    def get_position(self):
        "Return the position of the text as a tuple (*x*, *y*)"
        # This should return the same data (possibly unitized) as was
        # specified with set_x and set_y
        return self._dashx, self._dashy

    def get_prop_tup(self, renderer=None):
        """
        Return a hashable tuple of properties.

        Not intended to be human readable, but useful for backends who
        want to cache derived information about text (e.g., layouts) and
        need to know if the text has changed.
        """
        props = [p for p in Text.get_prop_tup(self, renderer=renderer)]
        props.extend([self._x, self._y, self._dashlength,
                      self._dashdirection, self._dashrotation, self._dashpad,
                      self._dashpush])
        return tuple(props)

    def draw(self, renderer):
        """
        Draw the :class:`TextWithDash` object to the given *renderer*.
        """
        self.update_coords(renderer)
        Text.draw(self, renderer)
        if self.get_dashlength() > 0.0:
            self.dashline.draw(renderer)
        self.stale = False

    def update_coords(self, renderer):
        """
        Computes the actual *x*, *y* coordinates for text based on the
        input *x*, *y* and the *dashlength*. Since the rotation is
        with respect to the actual canvas's coordinates we need to map
        back and forth.
        """
        dashx, dashy = self.get_unitless_position()
        dashlength = self.get_dashlength()
        # Shortcircuit this process if we don't have a dash
        if dashlength == 0.0:
            self._x, self._y = dashx, dashy
            return

        dashrotation = self.get_dashrotation()
        dashdirection = self.get_dashdirection()
        dashpad = self.get_dashpad()
        dashpush = self.get_dashpush()

        angle = get_rotation(dashrotation)
        theta = np.pi * (angle / 180.0 + dashdirection - 1)
        cos_theta, sin_theta = np.cos(theta), np.sin(theta)

        transform = self.get_transform()

        # Compute the dash end points
        # The 'c' prefix is for canvas coordinates
        cxy = transform.transform_point((dashx, dashy))
        cd = np.array([cos_theta, sin_theta])
        c1 = cxy + dashpush * cd
        c2 = cxy + (dashpush + dashlength) * cd

        inverse = transform.inverted()
        (x1, y1) = inverse.transform_point(tuple(c1))
        (x2, y2) = inverse.transform_point(tuple(c2))
        self.dashline.set_data((x1, x2), (y1, y2))

        # We now need to extend this vector out to
        # the center of the text area.
        # The basic problem here is that we're "rotating"
        # two separate objects but want it to appear as
        # if they're rotated together.
        # This is made non-trivial because of the
        # interaction between text rotation and alignment -
        # text alignment is based on the bbox after rotation.
        # We reset/force both alignments to 'center'
        # so we can do something relatively reasonable.
        # There's probably a better way to do this by
        # embedding all this in the object's transformations,
        # but I don't grok the transformation stuff
        # well enough yet.
        we = Text.get_window_extent(self, renderer=renderer)
        w, h = we.width, we.height
        # Watch for zeros
        if sin_theta == 0.0:
            dx = w
            dy = 0.0
        elif cos_theta == 0.0:
            dx = 0.0
            dy = h
        else:
            tan_theta = sin_theta / cos_theta
            dx = w
            dy = w * tan_theta
            if dy > h or dy < -h:
                dy = h
                dx = h / tan_theta
        cwd = np.array([dx, dy]) / 2
        cwd *= 1 + dashpad / np.sqrt(np.dot(cwd, cwd))
        cw = c2 + (dashdirection * 2 - 1) * cwd

        newx, newy = inverse.transform_point(tuple(cw))
        self._x, self._y = newx, newy

        # Now set the window extent
        # I'm not at all sure this is the right way to do this.
        we = Text.get_window_extent(self, renderer=renderer)
        self._twd_window_extent = we.frozen()
        self._twd_window_extent.update_from_data_xy(np.array([c1]), False)

        # Finally, make text align center
        Text.set_horizontalalignment(self, 'center')
        Text.set_verticalalignment(self, 'center')

    def get_window_extent(self, renderer=None):
        '''
        Return a :class:`~matplotlib.transforms.Bbox` object bounding
        the text, in display units.

        In addition to being used internally, this is useful for
        specifying clickable regions in a png file on a web page.

        *renderer* defaults to the _renderer attribute of the text
        object.  This is not assigned until the first execution of
        :meth:`draw`, so you must use this kwarg if you want
        to call :meth:`get_window_extent` prior to the first
        :meth:`draw`.  For getting web page regions, it is
        simpler to call the method after saving the figure.
        '''
        self.update_coords(renderer)
        if self.get_dashlength() == 0.0:
            return Text.get_window_extent(self, renderer=renderer)
        else:
            return self._twd_window_extent

    def get_dashlength(self):
        """
        Get the length of the dash.
        """
        return self._dashlength

    def set_dashlength(self, dl):
        """
        Set the length of the dash.

        ACCEPTS: float (canvas units)
        """
        self._dashlength = dl
        self.stale = True

    def get_dashdirection(self):
        """
        Get the direction dash.  1 is before the text and 0 is after.
        """
        return self._dashdirection

    def set_dashdirection(self, dd):
        """
        Set the direction of the dash following the text.
        1 is before the text and 0 is after. The default
        is 0, which is what you'd want for the typical
        case of ticks below and on the left of the figure.

        ACCEPTS: int (1 is before, 0 is after)
        """
        self._dashdirection = dd
        self.stale = True

    def get_dashrotation(self):
        """
        Get the rotation of the dash in degrees.
        """
        if self._dashrotation is None:
            return self.get_rotation()
        else:
            return self._dashrotation

    def set_dashrotation(self, dr):
        """
        Set the rotation of the dash, in degrees

        ACCEPTS: float (degrees)
        """
        self._dashrotation = dr
        self.stale = True

    def get_dashpad(self):
        """
        Get the extra spacing between the dash and the text, in canvas units.
        """
        return self._dashpad

    def set_dashpad(self, dp):
        """
        Set the "pad" of the TextWithDash, which is the extra spacing
        between the dash and the text, in canvas units.

        ACCEPTS: float (canvas units)
        """
        self._dashpad = dp
        self.stale = True

    def get_dashpush(self):
        """
        Get the extra spacing between the dash and the specified text
        position, in canvas units.
        """
        return self._dashpush

    def set_dashpush(self, dp):
        """
        Set the "push" of the TextWithDash, which
        is the extra spacing between the beginning
        of the dash and the specified position.

        ACCEPTS: float (canvas units)
        """
        self._dashpush = dp
        self.stale = True

    def set_position(self, xy):
        """
        Set the (*x*, *y*) position of the :class:`TextWithDash`.

        ACCEPTS: (x, y)
        """
        self.set_x(xy[0])
        self.set_y(xy[1])

    def set_x(self, x):
        """
        Set the *x* position of the :class:`TextWithDash`.

        ACCEPTS: float
        """
        self._dashx = float(x)
        self.stale = True

    def set_y(self, y):
        """
        Set the *y* position of the :class:`TextWithDash`.

        ACCEPTS: float
        """
        self._dashy = float(y)
        self.stale = True

    def set_transform(self, t):
        """
        Set the :class:`matplotlib.transforms.Transform` instance used
        by this artist.

        ACCEPTS: a :class:`matplotlib.transforms.Transform` instance
        """
        Text.set_transform(self, t)
        self.dashline.set_transform(t)
        self.stale = True

    def get_figure(self):
        'return the figure instance the artist belongs to'
        return self.figure

    def set_figure(self, fig):
        """
        Set the figure instance the artist belong to.

        ACCEPTS: a :class:`matplotlib.figure.Figure` instance
        """
        Text.set_figure(self, fig)
        self.dashline.set_figure(fig)

docstring.interpd.update(TextWithDash=artist.kwdoc(TextWithDash))


class OffsetFrom(object):
    'Callable helper class for working with `Annotation`'
    def __init__(self, artist, ref_coord, unit="points"):
        '''
        Parameters
        ----------
        artist : `Artist`, `BboxBase`, or `Transform`
            The object to compute the offset from.

        ref_coord : length 2 sequence
            If `artist` is an `Artist` or `BboxBase`, this values is
            the location to of the offset origin in fractions of the
            `artist` bounding box.

            If `artist` is a transform, the offset origin is the
            transform applied to this value.

        unit : {'points, 'pixels'}
            The screen units to use (pixels or points) for the offset
            input.

        '''
        self._artist = artist
        self._ref_coord = ref_coord
        self.set_unit(unit)

    def set_unit(self, unit):
        '''
        The unit for input to the transform used by ``__call__``

        Parameters
        ----------
        unit : {'points', 'pixels'}
        '''
        if unit not in ["points", "pixels"]:
            raise ValueError("'unit' must be one of [ 'points' | 'pixels' ]")
        self._unit = unit

    def get_unit(self):
        'The unit for input to the transform used by ``__call__``'
        return self._unit

    def _get_scale(self, renderer):
        unit = self.get_unit()
        if unit == "pixels":
            return 1.
        else:
            return renderer.points_to_pixels(1.)

    def __call__(self, renderer):
        '''
        Return the offset transform.

        Parameters
        ----------
        renderer : `RendererBase`
            The renderer to use to compute the offset

        Returns
        -------
        transform : `Transform`
            Maps (x, y) in pixel or point units to screen units
            relative to the given artist.
        '''
        if isinstance(self._artist, Artist):
            bbox = self._artist.get_window_extent(renderer)
            l, b, w, h = bbox.bounds
            xf, yf = self._ref_coord
            x, y = l + w * xf, b + h * yf
        elif isinstance(self._artist, BboxBase):
            l, b, w, h = self._artist.bounds
            xf, yf = self._ref_coord
            x, y = l + w * xf, b + h * yf
        elif isinstance(self._artist, Transform):
            x, y = self._artist.transform_point(self._ref_coord)
        else:
            raise RuntimeError("unknown type")

        sc = self._get_scale(renderer)
        tr = Affine2D().scale(sc, sc).translate(x, y)

        return tr


class _AnnotationBase(object):
    def __init__(self,
                 xy,
                 xycoords='data',
                 annotation_clip=None):

        self.xy = xy
        self.xycoords = xycoords
        self.set_annotation_clip(annotation_clip)

        self._draggable = None

    def _get_xy(self, renderer, x, y, s):
        if isinstance(s, tuple):
            s1, s2 = s
        else:
            s1, s2 = s, s

        if s1 == 'data':
            x = float(self.convert_xunits(x))
        if s2 == 'data':
            y = float(self.convert_yunits(y))

        tr = self._get_xy_transform(renderer, s)
        x1, y1 = tr.transform_point((x, y))
        return x1, y1

    def _get_xy_transform(self, renderer, s):

        if isinstance(s, tuple):
            s1, s2 = s
            from matplotlib.transforms import blended_transform_factory
            tr1 = self._get_xy_transform(renderer, s1)
            tr2 = self._get_xy_transform(renderer, s2)
            tr = blended_transform_factory(tr1, tr2)
            return tr
        elif callable(s):
            tr = s(renderer)
            if isinstance(tr, BboxBase):
                return BboxTransformTo(tr)
            elif isinstance(tr, Transform):
                return tr
            else:
                raise RuntimeError("unknown return type ...")
        elif isinstance(s, Artist):
            bbox = s.get_window_extent(renderer)
            return BboxTransformTo(bbox)
        elif isinstance(s, BboxBase):
            return BboxTransformTo(s)
        elif isinstance(s, Transform):
            return s
        elif not isinstance(s, six.string_types):
            raise RuntimeError("unknown coordinate type : %s" % (s,))

        if s == 'data':
            return self.axes.transData
        elif s == 'polar':
            from matplotlib.projections import PolarAxes
            tr = PolarAxes.PolarTransform()
            trans = tr + self.axes.transData
            return trans

        s_ = s.split()
        if len(s_) != 2:
            raise ValueError("%s is not a recognized coordinate" % s)

        bbox0, xy0 = None, None

        bbox_name, unit = s_
        # if unit is offset-like
        if bbox_name == "figure":
            bbox0 = self.figure.bbox
        elif bbox_name == "axes":
            bbox0 = self.axes.bbox
        # elif bbox_name == "bbox":
        #     if bbox is None:
        #         raise RuntimeError("bbox is specified as a coordinate but "
        #                            "never set")
        #     bbox0 = self._get_bbox(renderer, bbox)

        if bbox0 is not None:
            xy0 = bbox0.bounds[:2]
        elif bbox_name == "offset":
            xy0 = self._get_ref_xy(renderer)

        if xy0 is not None:
            # reference x, y in display coordinate
            ref_x, ref_y = xy0
            from matplotlib.transforms import Affine2D
            if unit == "points":
                # dots per points
                dpp = self.figure.get_dpi() / 72.
                tr = Affine2D().scale(dpp, dpp)
            elif unit == "pixels":
                tr = Affine2D()
            elif unit == "fontsize":
                fontsize = self.get_size()
                dpp = fontsize * self.figure.get_dpi() / 72.
                tr = Affine2D().scale(dpp, dpp)
            elif unit == "fraction":
                w, h = bbox0.bounds[2:]
                tr = Affine2D().scale(w, h)
            else:
                raise ValueError("%s is not a recognized coordinate" % s)

            return tr.translate(ref_x, ref_y)

        else:
            raise ValueError("%s is not a recognized coordinate" % s)

    def _get_ref_xy(self, renderer):
        """
        return x, y (in display coordinate) that is to be used for a reference
        of any offset coordinate
        """

        if isinstance(self.xycoords, tuple):
            s1, s2 = self.xycoords
            if ((isinstance(s1, six.string_types)
                 and s1.split()[0] == "offset")
                    or (isinstance(s2, six.string_types)
                        and s2.split()[0] == "offset")):
                raise ValueError("xycoords should not be an offset coordinate")
            x, y = self.xy
            x1, y1 = self._get_xy(renderer, x, y, s1)
            x2, y2 = self._get_xy(renderer, x, y, s2)
            return x1, y2
        elif (isinstance(self.xycoords, six.string_types) and
              self.xycoords.split()[0] == "offset"):
            raise ValueError("xycoords should not be an offset coordinate")
        else:
            x, y = self.xy
            return self._get_xy(renderer, x, y, self.xycoords)
        #raise RuntimeError("must be defined by the derived class")

    # def _get_bbox(self, renderer):
    #     if hasattr(bbox, "bounds"):
    #         return bbox
    #     elif hasattr(bbox, "get_window_extent"):
    #         bbox = bbox.get_window_extent()
    #         return bbox
    #     else:
    #         raise ValueError("A bbox instance is expected but got %s" %
    #                          str(bbox))

    def set_annotation_clip(self, b):
        """
        set *annotation_clip* attribute.

          * True: the annotation will only be drawn when self.xy is inside
                  the axes.
          * False: the annotation will always be drawn regardless of its
                   position.
          * None: the self.xy will be checked only if *xycoords* is "data"
        """
        self._annotation_clip = b

    def get_annotation_clip(self):
        """
        Return *annotation_clip* attribute.
        See :meth:`set_annotation_clip` for the meaning of return values.
        """
        return self._annotation_clip

    def _get_position_xy(self, renderer):
        "Return the pixel position of the annotated point."
        x, y = self.xy
        return self._get_xy(renderer, x, y, self.xycoords)

    def _check_xy(self, renderer, xy_pixel):
        """
        given the xy pixel coordinate, check if the annotation need to
        be drawn.
        """

        b = self.get_annotation_clip()

        if b or (b is None and self.xycoords == "data"):
            # check if self.xy is inside the axes.
            if not self.axes.contains_point(xy_pixel):
                return False

        return True

    def draggable(self, state=None, use_blit=False):
        """
        Set the draggable state -- if state is

          * None : toggle the current state

          * True : turn draggable on

          * False : turn draggable off

        If draggable is on, you can drag the annotation on the canvas with
        the mouse.  The DraggableAnnotation helper instance is returned if
        draggable is on.
        """
        from matplotlib.offsetbox import DraggableAnnotation
        is_draggable = self._draggable is not None

        # if state is None we'll toggle
        if state is None:
            state = not is_draggable

        if state:
            if self._draggable is None:
                self._draggable = DraggableAnnotation(self, use_blit)
        else:
            if self._draggable is not None:
                self._draggable.disconnect()
            self._draggable = None

        return self._draggable


class Annotation(Text, _AnnotationBase):
    def __str__(self):
        return "Annotation(%g,%g,%s)" % (self.xy[0],
                                         self.xy[1],
                                         repr(self._text))

    @docstring.dedent_interpd
    def __init__(self, s, xy,
                 xytext=None,
                 xycoords='data',
                 textcoords=None,
                 arrowprops=None,
                 annotation_clip=None,
                 **kwargs):
        '''
        Annotate the point ``xy`` with text ``s``.

        Additional kwargs are passed to `~matplotlib.text.Text`.

        Parameters
        ----------

        s : str
            The text of the annotation

        xy : iterable
            Length 2 sequence specifying the *(x,y)* point to annotate

        xytext : iterable, optional
            Length 2 sequence specifying the *(x,y)* to place the text
            at.  If None, defaults to ``xy``.

        xycoords : str, Artist, Transform, callable or tuple, optional

            The coordinate system that ``xy`` is given in.

            For a `str` the allowed values are:

            =================   ===============================================
            Property            Description
            =================   ===============================================
            'figure points'     points from the lower left of the figure
            'figure pixels'     pixels from the lower left of the figure
            'figure fraction'   fraction of figure from lower left
            'axes points'       points from lower left corner of axes
            'axes pixels'       pixels from lower left corner of axes
            'axes fraction'     fraction of axes from lower left
            'data'              use the coordinate system of the object being
                                annotated (default)
            'polar'             *(theta,r)* if not native 'data' coordinates
            =================   ===============================================

            If a `~matplotlib.artist.Artist` object is passed in the units are
            fraction if it's bounding box.

            If a `~matplotlib.transforms.Transform` object is passed
            in use that to transform ``xy`` to screen coordinates

            If a callable it must take a
            `~matplotlib.backend_bases.RendererBase` object as input
            and return a `~matplotlib.transforms.Transform` or
            `~matplotlib.transforms.Bbox` object

            If a `tuple` must be length 2 tuple of str, `Artist`,
            `Transform` or callable objects.  The first transform is
            used for the *x* coordinate and the second for *y*.

            See :ref:`plotting-guide-annotation` for more details.

            Defaults to ``'data'``

        textcoords : str, `Artist`, `Transform`, callable or tuple, optional
            The coordinate system that ``xytext`` is given, which
            may be different than the coordinate system used for
            ``xy``.

            All ``xycoords`` values are valid as well as the following
            strings:

            =================   =========================================
            Property            Description
            =================   =========================================
            'offset points'     offset (in points) from the *xy* value
            'offset pixels'     offset (in pixels) from the *xy* value
            =================   =========================================

            defaults to the input of ``xycoords``

        arrowprops : dict, optional
            If not None, properties used to draw a
            `~matplotlib.patches.FancyArrowPatch` arrow between ``xy`` and
            ``xytext``.

            If `arrowprops` does not contain the key ``'arrowstyle'`` the
            allowed keys are:

            ==========   ======================================================
            Key          Description
            ==========   ======================================================
            width        the width of the arrow in points
            headwidth    the width of the base of the arrow head in points
            headlength   the length of the arrow head in points
            shrink       fraction of total length to 'shrink' from both ends
            ?            any key to :class:`matplotlib.patches.FancyArrowPatch`
            ==========   ======================================================

            If the `arrowprops` contains the key ``'arrowstyle'`` the
            above keys are forbidden.  The allowed values of
            ``'arrowstyle'`` are:

            ============   =============================================
            Name           Attrs
            ============   =============================================
            ``'-'``        None
            ``'->'``       head_length=0.4,head_width=0.2
            ``'-['``       widthB=1.0,lengthB=0.2,angleB=None
            ``'|-|'``      widthA=1.0,widthB=1.0
            ``'-|>'``      head_length=0.4,head_width=0.2
            ``'<-'``       head_length=0.4,head_width=0.2
            ``'<->'``      head_length=0.4,head_width=0.2
            ``'<|-'``      head_length=0.4,head_width=0.2
            ``'<|-|>'``    head_length=0.4,head_width=0.2
            ``'fancy'``    head_length=0.4,head_width=0.4,tail_width=0.4
            ``'simple'``   head_length=0.5,head_width=0.5,tail_width=0.2
            ``'wedge'``    tail_width=0.3,shrink_factor=0.5
            ============   =============================================

            Valid keys for `~matplotlib.patches.FancyArrowPatch` are:

            ===============  ==================================================
            Key              Description
            ===============  ==================================================
            arrowstyle       the arrow style
            connectionstyle  the connection style
            relpos           default is (0.5, 0.5)
            patchA           default is bounding box of the text
            patchB           default is None
            shrinkA          default is 2 points
            shrinkB          default is 2 points
            mutation_scale   default is text size (in points)
            mutation_aspect  default is 1.
            ?                any key for :class:`matplotlib.patches.PathPatch`
            ===============  ==================================================

            Defaults to None

        annotation_clip : bool, optional
            Controls the visibility of the annotation when it goes
            outside the axes area.

            If `True`, the annotation will only be drawn when the
            ``xy`` is inside the axes. If `False`, the annotation will
            always be drawn regardless of its position.

            The default is `None`, which behave as `True` only if
            *xycoords* is "data".

        Returns
        -------
        Annotation

        '''

        _AnnotationBase.__init__(self,
                                 xy,
                                 xycoords=xycoords,
                                 annotation_clip=annotation_clip)
        # warn about wonky input data
        if (xytext is None and
                textcoords is not None and
                textcoords != xycoords):
            warnings.warn("You have used the `textcoords` kwarg, but not "
                          "the `xytext` kwarg.  This can lead to surprising "
                          "results.")

        # clean up textcoords and assign default
        if textcoords is None:
            textcoords = self.xycoords
        self._textcoords = textcoords

        # cleanup xytext defaults
        if xytext is None:
            xytext = self.xy
        x, y = xytext

        Text.__init__(self, x, y, s, **kwargs)

        self.arrowprops = arrowprops

        self.arrow = None

        if arrowprops is not None:
            if "arrowstyle" in arrowprops:
                arrowprops = self.arrowprops.copy()
                self._arrow_relpos = arrowprops.pop("relpos", (0.5, 0.5))
            else:
                # modified YAArrow API to be used with FancyArrowPatch
                shapekeys = ('width', 'headwidth', 'headlength',
                             'shrink', 'frac')
                arrowprops = dict()
                for key, val in self.arrowprops.items():
                    if key not in shapekeys:
                        arrowprops[key] = val  # basic Patch properties
            self.arrow_patch = FancyArrowPatch((0, 0), (1, 1),
                                               **arrowprops)
        else:
            self.arrow_patch = None

    def contains(self, event):
        contains, tinfo = Text.contains(self, event)
        if self.arrow is not None:
            in_arrow, _ = self.arrow.contains(event)
            contains = contains or in_arrow
        if self.arrow_patch is not None:
            in_patch, _ = self.arrow_patch.contains(event)
            contains = contains or in_patch

        return contains, tinfo

    @property
    def xyann(self):
        return self.get_position()

    @xyann.setter
    def xyann(self, xytext):
        self.set_position(xytext)

    @property
    def anncoords(self):
        return self._textcoords

    @anncoords.setter
    def anncoords(self, coords):
        self._textcoords = coords

    def set_figure(self, fig):

        if self.arrow is not None:
            self.arrow.set_figure(fig)
        if self.arrow_patch is not None:
            self.arrow_patch.set_figure(fig)
        Artist.set_figure(self, fig)

    def update_positions(self, renderer):
        """"Update the pixel positions of the annotated point and the
        text.
        """
        xy_pixel = self._get_position_xy(renderer)
        self._update_position_xytext(renderer, xy_pixel)

    def _update_position_xytext(self, renderer, xy_pixel):
        """Update the pixel positions of the annotation text and the arrow
        patch.
        """
        # generate transformation,
        self.set_transform(self._get_xy_transform(renderer, self.anncoords))

        ox0, oy0 = self._get_xy_display()
        ox1, oy1 = xy_pixel

        if self.arrowprops is not None:
            x0, y0 = xy_pixel
            l, b, w, h = Text.get_window_extent(self, renderer).bounds
            r = l + w
            t = b + h
            xc = 0.5 * (l + r)
            yc = 0.5 * (b + t)

            d = self.arrowprops.copy()
            ms = d.pop("mutation_scale", self.get_size())
            self.arrow_patch.set_mutation_scale(ms)

            if "arrowstyle" not in d:
                # Approximately simulate the YAArrow.
                # Pop its kwargs:
                shrink = d.pop('shrink', 0.0)
                width = d.pop('width', 4)
                headwidth = d.pop('headwidth', 12)
                # Ignore frac--it is useless.
                frac = d.pop('frac', None)
                if frac is not None:
                    warnings.warn(
                        "'frac' option in 'arrowprops' is no longer supported;"
                        " use 'headlength' to set the head length in points.")
                headlength = d.pop('headlength', 12)

                # NB: ms is in pts
                stylekw = dict(head_length=headlength / ms,
                               head_width=headwidth / ms,
                               tail_width=width / ms)

                self.arrow_patch.set_arrowstyle('simple', **stylekw)

                # using YAArrow style:
                # pick the x,y corner of the text bbox closest to point
                # annotated
                xpos = ((l, 0), (xc, 0.5), (r, 1))
                ypos = ((b, 0), (yc, 0.5), (t, 1))

                _, (x, relposx) = min((abs(val[0] - x0), val) for val in xpos)
                _, (y, relposy) = min((abs(val[0] - y0), val) for val in ypos)

                self._arrow_relpos = (relposx, relposy)

                r = np.hypot((y - y0), (x - x0))
                shrink_pts = shrink * r / renderer.points_to_pixels(1)
                self.arrow_patch.shrinkA = shrink_pts
                self.arrow_patch.shrinkB = shrink_pts

            # adjust the starting point of the arrow relative to
            # the textbox.
            # TODO : Rotation needs to be accounted.
            relpos = self._arrow_relpos
            bbox = Text.get_window_extent(self, renderer)
            ox0 = bbox.x0 + bbox.width * relpos[0]
            oy0 = bbox.y0 + bbox.height * relpos[1]

            # The arrow will be drawn from (ox0, oy0) to (ox1,
            # oy1). It will be first clipped by patchA and patchB.
            # Then it will be shrunk by shrinkA and shrinkB
            # (in points). If patch A is not set, self.bbox_patch
            # is used.

            self.arrow_patch.set_positions((ox0, oy0), (ox1, oy1))

            if "patchA" in d:
                self.arrow_patch.set_patchA(d.pop("patchA"))
            else:
                if self._bbox_patch:
                    self.arrow_patch.set_patchA(self._bbox_patch)
                else:
                    pad = renderer.points_to_pixels(4)
                    if self.get_text() == "":
                        self.arrow_patch.set_patchA(None)
                        return

                    bbox = Text.get_window_extent(self, renderer)
                    l, b, w, h = bbox.bounds
                    l -= pad / 2.
                    b -= pad / 2.
                    w += pad
                    h += pad
                    r = Rectangle(xy=(l, b),
                                  width=w,
                                  height=h,
                                  )
                    r.set_transform(IdentityTransform())
                    r.set_clip_on(False)

                    self.arrow_patch.set_patchA(r)

    @artist.allow_rasterization
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

        self._update_position_xytext(renderer, xy_pixel)
        self.update_bbox_position_size(renderer)

        if self.arrow_patch is not None:   # FancyArrowPatch
            if self.arrow_patch.figure is None and self.figure is not None:
                self.arrow_patch.figure = self.figure
            self.arrow_patch.draw(renderer)

        # Draw text, including FancyBboxPatch, after FancyArrowPatch.
        # Otherwise, a wedge arrowstyle can land partly on top of the Bbox.
        Text.draw(self, renderer)

    def get_window_extent(self, renderer=None):
        '''
        Return a :class:`~matplotlib.transforms.Bbox` object bounding
        the text and arrow annotation, in display units.

        *renderer* defaults to the _renderer attribute of the text
        object.  This is not assigned until the first execution of
        :meth:`draw`, so you must use this kwarg if you want
        to call :meth:`get_window_extent` prior to the first
        :meth:`draw`.  For getting web page regions, it is
        simpler to call the method after saving the figure. The
        *dpi* used defaults to self.figure.dpi; the renderer dpi is
        irrelevant.

        '''
        if not self.get_visible():
            return Bbox.unit()
        arrow = self.arrow
        arrow_patch = self.arrow_patch

        text_bbox = Text.get_window_extent(self, renderer=renderer)
        bboxes = [text_bbox]

        if self.arrow is not None:
            bboxes.append(arrow.get_window_extent(renderer=renderer))
        elif self.arrow_patch is not None:
            bboxes.append(arrow_patch.get_window_extent(renderer=renderer))

        return Bbox.union(bboxes)


docstring.interpd.update(Annotation=Annotation.__init__.__doc__)
