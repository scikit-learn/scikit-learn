"""
axis_artist.py module provides axis-related artists. They are

 * axis line
 * tick lines
 * tick labels
 * axis label
 * grid lines

The main artist class is a AxisArtist and a GridlinesCollection. The
GridlinesCollection is responsible for drawing grid lines and the
AxisArtist is responsible for all other artists. The AxisArtist class
has attributes that are associated with each type of artists.

 * line : axis line
 * major_ticks : major tick lines
 * major_ticklabels : major tick labels
 * minor_ticks : minor tick lines
 * minor_ticklabels : minor tick labels
 * label : axis label

Typically, the AxisArtist associated with a axes will be accessed with
the *axis* dictionary of the axes, i.e., the AxisArtist for the bottom
axis is

  ax.axis["bottom"]

where *ax* is an instance of axes (mpl_toolkits.axislines.Axes).  Thus,
ax.axis["bottom"].line is an artist associated with the axis line, and
ax.axis["bottom"].major_ticks is an artist associated with the major tick
lines.

You can change the colors, fonts, line widths, etc. of these artists
by calling suitable set method. For example, to change the color of the major
ticks of the bottom axis to red,

  ax.axis["bottom"].major_ticks.set_color("r")

However, things like the locations of ticks, and their ticklabels need
to be changed from the side of the grid_helper.

axis_direction
--------------

AxisArtist, AxisLabel, TickLabels have *axis_direction* attribute,
which adjusts the location, angle, etc.,. The *axis_direction* must be
one of [left, right, bottom, top] and they follow the matplotlib
convention for the rectangle axis.

For example, for the *bottom* axis (the left and right is relative to
the direction of the increasing coordinate),

 * ticklabels and axislabel are on the right
 * ticklabels and axislabel have text angle of 0
 * ticklabels are baseline, center-aligned
 * axislabel is top, center-aligned


The text angles are actually relative to (90 + angle of the direction
to the ticklabel), which gives 0 for bottom axis.

                        left bottom right top
 ticklabels location    left right  right left
 axislabel location     left right  right left
 ticklabels angle       90    0      -90  180
 axislabel angle        180   0     0     180
 ticklabel va           center baseline center baseline
 axislabel va           center top      center bottom
 ticklabel ha           right  center   right  center
 axislabel ha           right  center   right  center


Ticks are by default direct opposite side of the ticklabels. To make
ticks to the same side of the ticklabels,

  ax.axis["bottom"].major_ticks.set_ticks_out(True)


Following attributes can be customized (use set_xxx method)

 * Ticks : ticksize, tick_out
 * TickLabels : pad
 * AxisLabel : pad

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

# FIXME :
# angles are given in data coordinate - need to convert it to canvas coordinate


import matplotlib.artist as martist
import matplotlib.text as mtext
import matplotlib.font_manager as font_manager

from matplotlib.path import Path
from matplotlib.transforms import (
    Affine2D, Bbox, IdentityTransform, ScaledTranslation, TransformedPath)
from matplotlib.collections import LineCollection

from matplotlib import rcParams

from matplotlib.artist import allow_rasterization

import warnings

import numpy as np


import matplotlib.lines as mlines
from .axisline_style import AxislineStyle


class BezierPath(mlines.Line2D):

    def __init__(self, path, *kl, **kw):
        mlines.Line2D.__init__(self, [], [], *kl, **kw)
        self._path = path
        self._invalid = False

    def recache(self):

        self._transformed_path = TransformedPath(self._path, self.get_transform())

        self._invalid = False

    def set_path(self, path):
        self._path = path
        self._invalid = True


    def draw(self, renderer):
        if self._invalid:
            self.recache()

        if not self._visible: return
        renderer.open_group('line2d')

        gc = renderer.new_gc()
        self._set_gc_clip(gc)

        gc.set_foreground(self._color)
        gc.set_antialiased(self._antialiased)
        gc.set_linewidth(self._linewidth)
        gc.set_alpha(self._alpha)
        if self.is_dashed():
            cap = self._dashcapstyle
            join = self._dashjoinstyle
        else:
            cap = self._solidcapstyle
            join = self._solidjoinstyle
        gc.set_joinstyle(join)
        gc.set_capstyle(cap)
        gc.set_dashes(self._dashOffset, self._dashSeq)

        if self._lineStyles[self._linestyle] != '_draw_nothing':
            tpath, affine = (
                self._transformed_path.get_transformed_path_and_affine())
            renderer.draw_path(gc, tpath, affine.frozen())

        gc.restore()
        renderer.close_group('line2d')



class UnimplementedException(Exception):
    pass

from matplotlib.artist import Artist

class AttributeCopier(object):
    def __init__(self, ref_artist, klass=Artist):
        self._klass = klass
        self._ref_artist = ref_artist
        super(AttributeCopier, self).__init__()

    def set_ref_artist(self, artist):
        self._ref_artist = artist

    def get_ref_artist(self):
        raise RuntimeError("get_ref_artist must overridden")
    #return self._ref_artist

    def get_attribute_from_ref_artist(self, attr_name, default_value):
        get_attr_method_name = "get_"+attr_name
        c = getattr(self._klass, get_attr_method_name)(self)
        if c == 'auto':
            ref_artist = self.get_ref_artist()
            if ref_artist:
                attr = getattr(ref_artist,
                               get_attr_method_name)()
                return attr
            else:
                return default_value

        return c


from matplotlib.lines import Line2D

class Ticks(Line2D, AttributeCopier):
    """
    Ticks are derived from Line2D, and note that ticks themselves
    are markers. Thus, you should use set_mec, set_mew, etc.

    To change the tick size (length), you need to use
    set_ticksize. To change the direction of the ticks (ticks are
    in opposite direction of ticklabels by default), use
    set_tick_out(False).
    """

    def __init__(self, ticksize, tick_out=False, **kwargs):
        self._ticksize = ticksize
        self.locs_angles_labels = []

        self.set_tick_out(tick_out)

        self._axis = kwargs.pop("axis", None)
        if self._axis is not None:
            if "color" not in kwargs:
                kwargs["color"] = "auto"
            if ("mew" not in kwargs) and ("markeredgewidth" not in kwargs):
                kwargs["markeredgewidth"] = "auto"

        Line2D.__init__(self, [0.], [0.], **kwargs)
        AttributeCopier.__init__(self, self._axis, klass=Line2D)
        self.set_snap(True)

    def get_ref_artist(self):
        #return self._ref_artist.get_ticklines()[0]
        return self._ref_artist.majorTicks[0].tick1line

    def get_color(self):
        return self.get_attribute_from_ref_artist("color", "k")

    def get_markeredgecolor(self):
        if self._markeredgecolor == 'auto':
            return self.get_color()
        else:
            return self._markeredgecolor

    def get_markeredgewidth(self):
        return self.get_attribute_from_ref_artist("markeredgewidth", .5)


    def set_tick_out(self, b):
        """
        set True if tick need to be rotated by 180 degree.
        """
        self._tick_out = b

    def get_tick_out(self):
        """
        Return True if the tick will be rotated by 180 degree.
        """
        return self._tick_out


    def set_ticksize(self, ticksize):
        """
        set length of the ticks in points.
        """
        self._ticksize = ticksize


    def get_ticksize(self):
        """
        Return length of the ticks in points.
        """
        return self._ticksize

    def set_locs_angles(self, locs_angles):
        self.locs_angles = locs_angles


    def _update(self, renderer):
        pass

    _tickvert_path = Path([[0., 0.], [1., 0.]])

    def draw(self, renderer):
        if not self.get_visible():
            return

        self._update(renderer) # update the tick

        size = self._ticksize
        path_trans = self.get_transform()

        # set gc : copied from lines.py
#         gc = renderer.new_gc()
#         self._set_gc_clip(gc)

#         gc.set_foreground(self.get_color())
#         gc.set_antialiased(self._antialiased)
#         gc.set_linewidth(self._linewidth)
#         gc.set_alpha(self._alpha)
#         if self.is_dashed():
#             cap = self._dashcapstyle
#             join = self._dashjoinstyle
#         else:
#             cap = self._solidcapstyle
#             join = self._solidjoinstyle
#         gc.set_joinstyle(join)
#         gc.set_capstyle(cap)
#         gc.set_snap(self.get_snap())


        gc = renderer.new_gc()
        gc.set_foreground(self.get_markeredgecolor())
        gc.set_linewidth(self.get_markeredgewidth())
        gc.set_alpha(self._alpha)

        offset = renderer.points_to_pixels(size)
        marker_scale = Affine2D().scale(offset, offset)

        if self.get_tick_out():
            add_angle = 180
        else:
            add_angle = 0

        marker_rotation = Affine2D()
        marker_transform = marker_scale + marker_rotation

        for loc, angle in self.locs_angles:
            marker_rotation.clear().rotate_deg(angle+add_angle)
            locs = path_trans.transform_non_affine([loc])
            if (self.axes and
                    not self.axes.viewLim.contains(locs[0][0], locs[0][1])):
                continue
            renderer.draw_markers(gc, self._tickvert_path, marker_transform,
                                  Path(locs), path_trans.get_affine())

        gc.restore()


class LabelBase(mtext.Text):
    """
    A base class for AxisLabel and TickLabels. The position and angle
    of the text are calculated by to offset_ref_angle,
    text_ref_angle, and offset_radius attributes.
    """

    def __init__(self, *kl, **kwargs):
        self.locs_angles_labels = []
        self._ref_angle = 0
        self._offset_radius = 0.

        super(LabelBase, self).__init__(*kl,
                                        **kwargs)

        self.set_rotation_mode("anchor")
        self._text_follow_ref_angle = True
        #self._offset_ref_angle = 0

    def _set_ref_angle(self, a):
        self._ref_angle = a

    def _get_ref_angle(self):
        return self._ref_angle

    def _get_text_ref_angle(self):
        if self._text_follow_ref_angle:
            return self._get_ref_angle()+90
        else:
            return 0 #self.get_ref_angle()

    def _get_offset_ref_angle(self):
        return self._get_ref_angle()

    def _set_offset_radius(self, offset_radius):
        self._offset_radius = offset_radius

    def _get_offset_radius(self):
        return self._offset_radius


    _get_opposite_direction = {"left":"right",
                               "right":"left",
                               "top":"bottom",
                               "bottom":"top"}.__getitem__


    def _update(self, renderer):
        pass

    def draw(self, renderer):
        if not self.get_visible(): return

        self._update(renderer)

        # save original and adjust some properties
        tr = self.get_transform()
        angle_orig = self.get_rotation()

        offset_tr = Affine2D()
        self.set_transform(tr+offset_tr)

        text_ref_angle = self._get_text_ref_angle()
        offset_ref_angle = self._get_offset_ref_angle()

        theta = (offset_ref_angle)/180.*np.pi
        dd = self._get_offset_radius()
        dx, dy = dd * np.cos(theta), dd * np.sin(theta)
        offset_tr.translate(dx, dy)
        self.set_rotation(text_ref_angle+angle_orig)
        super(LabelBase, self).draw(renderer)
        offset_tr.clear()


        # restore original properties
        self.set_transform(tr)
        self.set_rotation(angle_orig)


    def get_window_extent(self, renderer):

        self._update(renderer)

        # save original and adjust some properties
        tr = self.get_transform()
        angle_orig = self.get_rotation()

        offset_tr = Affine2D()
        self.set_transform(tr+offset_tr)

        text_ref_angle = self._get_text_ref_angle()
        offset_ref_angle = self._get_offset_ref_angle()

        theta = (offset_ref_angle)/180.*np.pi
        dd = self._get_offset_radius()
        dx, dy = dd * np.cos(theta), dd * np.sin(theta)
        offset_tr.translate(dx, dy)
        self.set_rotation(text_ref_angle+angle_orig)

        bbox = super(LabelBase, self).get_window_extent(renderer).frozen()

        offset_tr.clear()


        # restore original properties
        self.set_transform(tr)
        self.set_rotation(angle_orig)

        return bbox


class AxisLabel(LabelBase, AttributeCopier):
    """
    Axis Label. Derived from Text. The position of the text is updated
    in the fly, so changing text position has no effect. Otherwise, the
    properties can be changed as a normal Text.

    To change the pad between ticklabels and axis label, use set_pad.
    """

    def __init__(self, *kl, **kwargs):

        axis_direction = kwargs.pop("axis_direction", "bottom")
        self._axis = kwargs.pop("axis", None)
        #super(AxisLabel, self).__init__(*kl, **kwargs)
        LabelBase.__init__(self, *kl, **kwargs)
        AttributeCopier.__init__(self, self._axis, klass=LabelBase)

        self.set_axis_direction(axis_direction)
        self._pad = 5
        self._extra_pad = 0

    def set_pad(self, pad):
        """
        Set the pad in points. Note that the actual pad will be the
        sum of the internal pad and the external pad (that are set
        automatically by the AxisArtist), and it only set the internal
        pad
        """
        self._pad = pad

    def get_pad(self):
        """
        return pad in points. See set_pad for more details.
        """
        return self._pad


    def _set_external_pad(self, p):
        """
        Set external pad IN PIXELS. This is intended to be set by the
        AxisArtist, bot by user..
        """
        self._extra_pad = p

    def _get_external_pad(self):
        """
        Get external pad.
        """
        return self._extra_pad


    def get_ref_artist(self):
        return self._axis.get_label()


    def get_text(self):
        t = super(AxisLabel, self).get_text()
        if t == "__from_axes__":
            return self._axis.get_label().get_text()
        return self._text

    _default_alignments = dict(left=("bottom", "center"),
                               right=("top", "center"),
                               bottom=("top", "center"),
                               top=("bottom", "center"))



    def set_default_alignment(self, d):
        if d not in ["left", "right", "top", "bottom"]:
            raise ValueError('direction must be on of "left", "right", "top", "bottom"')

        va, ha = self._default_alignments[d]
        self.set_va(va)
        self.set_ha(ha)


    _default_angles = dict(left=180,
                           right=0,
                           bottom=0,
                           top=180)


    def set_default_angle(self, d):
        if d not in ["left", "right", "top", "bottom"]:
            raise ValueError('direction must be on of "left", "right", "top", "bottom"')

        self.set_rotation(self._default_angles[d])


    def set_axis_direction(self, d):
        """
        Adjust the text angle and text alignment of axis label
        according to the matplotlib convention.


        =====================    ========== ========= ========== ==========
        property                 left       bottom    right      top
        =====================    ========== ========= ========== ==========
        axislabel angle          180        0         0          180
        axislabel va             center     top       center     bottom
        axislabel ha             right      center    right      center
        =====================    ========== ========= ========== ==========

        Note that the text angles are actually relative to (90 + angle
        of the direction to the ticklabel), which gives 0 for bottom
        axis.

        """
        if d not in ["left", "right", "top", "bottom"]:
            raise ValueError('direction must be on of "left", "right", "top", "bottom"')

        self.set_default_alignment(d)
        self.set_default_angle(d)

    def get_color(self):
        return self.get_attribute_from_ref_artist("color", "k")

    def draw(self, renderer):
        if not self.get_visible():
            return

        pad = renderer.points_to_pixels(self.get_pad())
        r = self._get_external_pad() + pad
        self._set_offset_radius(r)

        super(AxisLabel, self).draw(renderer)


    def get_window_extent(self, renderer):

        if not self.get_visible():
            return

        pad = renderer.points_to_pixels(self.get_pad())
        r = self._get_external_pad() + pad
        self._set_offset_radius(r)

        bb = super(AxisLabel, self).get_window_extent(renderer)

        return bb


class TickLabels(AxisLabel, AttributeCopier): # mtext.Text
    """
    Tick Labels. While derived from Text, this single artist draws all
    ticklabels. As in AxisLabel, the position of the text is updated
    in the fly, so changing text position has no effect. Otherwise,
    the properties can be changed as a normal Text. Unlike the
    ticklabels of the mainline matplotlib, properties of single
    ticklabel alone cannot modified.

    To change the pad between ticks and ticklabels, use set_pad.
    """

    def __init__(self, **kwargs):

        axis_direction = kwargs.pop("axis_direction", "bottom")
        AxisLabel.__init__(self, **kwargs)
        self.set_axis_direction(axis_direction)
        #self._axis_direction = axis_direction
        self._axislabel_pad = 0
        #self._extra_pad = 0


    # attribute copier
    def get_ref_artist(self):
        return self._axis.get_ticklabels()[0]

    def set_axis_direction(self, label_direction):
        """
        Adjust the text angle and text alignment of ticklabels
        according to the matplotlib convention.

        The *label_direction* must be one of [left, right, bottom,
        top].

        =====================    ========== ========= ========== ==========
        property                 left       bottom    right      top
        =====================    ========== ========= ========== ==========
        ticklabels angle         90         0         -90        180
        ticklabel va             center     baseline  center     baseline
        ticklabel ha             right      center    right      center
        =====================    ========== ========= ========== ==========


        Note that the text angles are actually relative to (90 + angle
        of the direction to the ticklabel), which gives 0 for bottom
        axis.

        """

        if label_direction not in ["left", "right", "top", "bottom"]:
            raise ValueError('direction must be one of "left", "right", "top", "bottom"')

        self._axis_direction = label_direction
        self.set_default_alignment(label_direction)
        self.set_default_angle(label_direction)


    def invert_axis_direction(self):
        label_direction = self._get_opposite_direction(self._axis_direction)
        self.set_axis_direction(label_direction)

    def _get_ticklabels_offsets(self, renderer, label_direction):
        """
        Calculates the offsets of the ticklabels from the tick and
        their total heights. The offset only takes account the offset
        due to the vertical alignment of the ticklabels, i.e.,if axis
        direction is bottom and va is ;top', it will return 0. if va
        is 'baseline', it will return (height-descent).
        """
        whd_list = self.get_texts_widths_heights_descents(renderer)

        if not whd_list:
            return 0, 0

        r = 0
        va, ha = self.get_va(), self.get_ha()

        if label_direction == "left":
            pad = max(w for w, h, d in whd_list)
            if ha == "left":
                r = pad
            elif ha == "center":
                r = .5 * pad
        elif label_direction == "right":
            pad = max(w for w, h, d in whd_list)
            if ha == "right":
                r = pad
            elif ha == "center":
                r = .5 * pad
        elif label_direction == "bottom":
            pad = max(h for w, h, d in whd_list)
            if va == "bottom":
                r = pad
            elif va == "center":
                r =.5 * pad
            elif va == "baseline":
                max_ascent = max(h - d for w, h, d in whd_list)
                max_descent = max(d for w, h, d in whd_list)
                r  = max_ascent
                pad = max_ascent + max_descent
        elif label_direction == "top":
            pad = max(h for w, h, d in whd_list)
            if va == "top":
                r = pad
            elif va == "center":
                r =.5 * pad
            elif va == "baseline":
                max_ascent = max(h - d for w, h, d in whd_list)
                max_descent = max(d for w, h, d in whd_list)
                r  = max_descent
                pad = max_ascent + max_descent

        #tick_pad = renderer.points_to_pixels(self.get_pad())

        # r : offset

        # pad : total height of the ticklabels. This will be used to
        # calculate the pad for the axislabel.
        return r, pad



    _default_alignments = dict(left=("center", "right"),
                               right=("center", "left"),
                               bottom=("baseline", "center"),
                               top=("baseline", "center"))



    # set_default_alignments(self, d)

    _default_angles = dict(left=90,
                           right=-90,
                           bottom=0,
                           top=180)


    def draw(self, renderer):
        if not self.get_visible():
            self._axislabel_pad = self._get_external_pad()
            return

        r, total_width = self._get_ticklabels_offsets(renderer,
                                                      self._axis_direction)

        #self._set_external_pad(r+self._get_external_pad())
        pad = self._get_external_pad() + \
              renderer.points_to_pixels(self.get_pad())
        self._set_offset_radius(r+pad)

        #self._set_offset_radius(r)

        for (x, y), a, l in self._locs_angles_labels:
            if not l.strip(): continue
            self._set_ref_angle(a) #+ add_angle
            self.set_x(x)
            self.set_y(y)
            self.set_text(l)
            LabelBase.draw(self, renderer)

        self._axislabel_pad = total_width \
                              + pad # the value saved will be used to draw axislabel.


    def set_locs_angles_labels(self, locs_angles_labels):
        self._locs_angles_labels = locs_angles_labels

    def get_window_extents(self, renderer):

        if not self.get_visible():
            self._axislabel_pad = self._get_external_pad()
            return []

        bboxes = []

        r, total_width = self._get_ticklabels_offsets(renderer,
                                                     self._axis_direction)

        pad = self._get_external_pad() + \
              renderer.points_to_pixels(self.get_pad())
        self._set_offset_radius(r+pad)


        for (x, y), a, l in self._locs_angles_labels:
            self._set_ref_angle(a) #+ add_angle
            self.set_x(x)
            self.set_y(y)
            self.set_text(l)
            bb = LabelBase.get_window_extent(self, renderer)
            bboxes.append(bb)

        self._axislabel_pad = total_width \
                              + pad # the value saved will be used to draw axislabel.

        return bboxes


    def get_texts_widths_heights_descents(self, renderer):
        """
        return a list of width, height, descent for ticklabels.
        """
        whd_list = []
        for (x, y), a, l in self._locs_angles_labels:
            if not l.strip(): continue
            clean_line, ismath = self.is_math_text(l)
            whd = renderer.get_text_width_height_descent(
                clean_line, self._fontproperties, ismath=ismath)
            whd_list.append(whd)

        return whd_list


class GridlinesCollection(LineCollection):
    def __init__(self, *kl, **kwargs):
        """
        *which* : "major" or "minor"
        *axis* : "both", "x" or "y"
        """
        self._which = kwargs.pop("which", "major")
        self._axis = kwargs.pop("axis", "both")
        super(GridlinesCollection, self).__init__(*kl, **kwargs)
        self.set_grid_helper(None)

    def set_which(self, which):
        self._which = which

    def set_axis(self, axis):
        self._axis = axis

    def set_grid_helper(self, grid_helper):
        self._grid_helper = grid_helper

    def draw(self, renderer):
        if self._grid_helper is not None:
            self._grid_helper.update_lim(self.axes)
            gl = self._grid_helper.get_gridlines(self._which, self._axis)
            if gl:
                self.set_segments([np.transpose(l) for l in gl])
            else:
                self.set_segments([])
        super(GridlinesCollection, self).draw(renderer)




class AxisArtist(martist.Artist):
    """
    An artist which draws axis (a line along which the n-th axes coord
    is constant) line, ticks, ticklabels, and axis label.
    """

    ZORDER=2.5

    @property
    def LABELPAD(self):
        return self.label.get_pad()

    @LABELPAD.setter
    def LABELPAD(self, v):
        return self.label.set_pad(v)

    def __init__(self, axes,
                 helper,
                 offset=None,
                 axis_direction="bottom",
                 **kw):
        """
        *axes* : axes
        *helper* : an AxisArtistHelper instance.
        """
        #axes is also used to follow the axis attribute (tick color, etc).

        super(AxisArtist, self).__init__(**kw)

        self.axes = axes

        self._axis_artist_helper = helper

        if offset is None:
            offset = (0, 0)
        self.dpi_transform = Affine2D()
        self.offset_transform = ScaledTranslation(offset[0], offset[1],
                                                  self.dpi_transform)

        self._label_visible = True
        self._majortick_visible = True
        self._majorticklabel_visible = True
        self._minortick_visible = True
        self._minorticklabel_visible = True


        #if self._axis_artist_helper._loc in ["left", "right"]:
        if axis_direction in ["left", "right"]:
            axis_name = "ytick"
            self.axis = axes.yaxis
        else:
            axis_name = "xtick"
            self.axis = axes.xaxis


        self._axisline_style = None


        self._axis_direction = axis_direction


        self._init_line()
        self._init_ticks(axis_name, **kw)
        self._init_offsetText(axis_direction)
        self._init_label()

        self.set_zorder(self.ZORDER)

        self._rotate_label_along_line = False

        # axis direction
        self._tick_add_angle = 180.
        self._ticklabel_add_angle = 0.
        self._axislabel_add_angle = 0.
        self.set_axis_direction(axis_direction)


    # axis direction

    def set_axis_direction(self, axis_direction):
        """
        Adjust the direction, text angle, text alignment of
        ticklabels, labels following the matplotlib convention for
        the rectangle axes.

        The *axis_direction* must be one of [left, right, bottom,
        top].

        =====================    ========== ========= ========== ==========
        property                 left       bottom    right      top
        =====================    ========== ========= ========== ==========
        ticklabels location      "-"        "+"       "+"        "-"
        axislabel location       "-"        "+"       "+"        "-"
        ticklabels angle         90         0         -90        180
        ticklabel va             center     baseline  center     baseline
        ticklabel ha             right      center    right      center
        axislabel angle          180        0         0          180
        axislabel va             center     top       center     bottom
        axislabel ha             right      center    right      center
        =====================    ========== ========= ========== ==========


        Note that the direction "+" and "-" are relative to the direction of
        the increasing coordinate. Also, the text angles are actually
        relative to (90 + angle of the direction to the ticklabel),
        which gives 0 for bottom axis.

        """

        if axis_direction not in ["left", "right", "top", "bottom"]:
            raise ValueError('direction must be on of "left", "right", "top", "bottom"')
        self._axis_direction = axis_direction
        if axis_direction in ["left", "top"]:
            #self._set_tick_direction("+")
            self.set_ticklabel_direction("-")
            self.set_axislabel_direction("-")
        else:
            #self._set_tick_direction("-")
            self.set_ticklabel_direction("+")
            self.set_axislabel_direction("+")

        self.major_ticklabels.set_axis_direction(axis_direction)
        self.label.set_axis_direction(axis_direction)

    # def _set_tick_direction(self, d):
    #     if d not in ["+", "-"]:
    #         raise ValueError('direction must be on of "in", "out"')

    #     if d == "+":
    #         self._tick_add_angle = 0 #get_helper()._extremes=0, 10
    #     else:
    #         self._tick_add_angle = 180 #get_helper()._extremes=0, 10

    def set_ticklabel_direction(self, tick_direction):
        """
        Adjust the direction of the ticklabel.

         ACCEPTS: [ "+" | "-" ]

        Note that the label_direction '+' and '-' are relative to the
        direction of the increasing coordinate.
        """

        if tick_direction not in ["+", "-"]:
            raise ValueError('direction must be one of "+", "-"')

        if tick_direction == "-":
            self._ticklabel_add_angle = 180
        else:
            self._ticklabel_add_angle = 0

    def invert_ticklabel_direction(self):
        self._ticklabel_add_angle = (self._ticklabel_add_angle + 180) % 360
        self.major_ticklabels.invert_axis_direction()
        self.minor_ticklabels.invert_axis_direction()

    # def invert_ticks_direction(self):
    #     self.major_ticks.set_tick_out(not self.major_ticks.get_tick_out())
    #     self.minor_ticks.set_tick_out(not self.minor_ticks.get_tick_out())

    def set_axislabel_direction(self, label_direction):
        """
        Adjust the direction of the axislabel.

         ACCEPTS: [ "+" | "-" ]

        Note that the label_direction '+' and '-' are relative to the
        direction of the increasing coordinate.
        """
        if label_direction not in ["+", "-"]:
            raise ValueError('direction must be one of "+", "-"')

        if label_direction == "-":
            self._axislabel_add_angle = 180
        else:
            self._axislabel_add_angle = 0



    def get_transform(self):
        return self.axes.transAxes + self.offset_transform

    def get_helper(self):
        """
        Return axis artist helper instance.
        """
        return self._axis_artist_helper


    def set_axisline_style(self, axisline_style=None, **kw):
        """
        Set the axisline style.

        *axisline_style* can be a string with axisline style name with optional
         comma-separated attributes. Alternatively, the attrs can
         be provided as keywords.

         set_arrowstyle("->,size=1.5")
         set_arrowstyle("->", size=1.5)

        Old attrs simply are forgotten.

        Without argument (or with arrowstyle=None), return
        available styles as a list of strings.
        """

        if axisline_style==None:
            return AxislineStyle.pprint_styles()

        if isinstance(axisline_style, AxislineStyle._Base):
            self._axisline_style = axisline_style
        else:
            self._axisline_style = AxislineStyle(axisline_style, **kw)


        self._init_line()


    def get_axisline_style(self):
        """
        return the current axisline style.
        """
        return self._axisline_style

    def _init_line(self):
        """
        Initialize the *line* artist that is responsible to draw the axis line.
        """
        tran = self._axis_artist_helper.get_line_transform(self.axes) \
               + self.offset_transform

        axisline_style = self.get_axisline_style()
        if axisline_style is None:
            self.line = BezierPath(self._axis_artist_helper.get_line(self.axes),
                                   color=rcParams['axes.edgecolor'],
                                   linewidth=rcParams['axes.linewidth'],
                                   transform=tran)
        else:
            self.line = axisline_style(self, transform=tran)

    def _draw_line(self, renderer):
        self.line.set_path(self._axis_artist_helper.get_line(self.axes))
        if self.get_axisline_style() is not None:
            self.line.set_line_mutation_scale(self.major_ticklabels.get_size())
        self.line.draw(renderer)


    def _init_ticks(self, axis_name, **kw):

        trans=self._axis_artist_helper.get_tick_transform(self.axes) \
               + self.offset_transform


        major_tick_size = kw.get("major_tick_size",
                                 rcParams['%s.major.size'%axis_name])
        major_tick_pad = kw.get("major_tick_pad",
                                rcParams['%s.major.pad'%axis_name])
        minor_tick_size = kw.get("minor_tick_size",
                                 rcParams['%s.minor.size'%axis_name])
        minor_tick_pad = kw.get("minor_tick_pad",
                                rcParams['%s.minor.pad'%axis_name])

        self.major_ticks = Ticks(major_tick_size,
                                 axis=self.axis,
                                 transform=trans)
        self.minor_ticks = Ticks(minor_tick_size,
                                 axis=self.axis,
                                 transform=trans)

        if axis_name == "xaxis":
            size = rcParams['xtick.labelsize']
        else:
            size = rcParams['ytick.labelsize']


        fontprops = font_manager.FontProperties(size=size)

        self.major_ticklabels = TickLabels(size=size, axis=self.axis,
                                           axis_direction=self._axis_direction)
        self.minor_ticklabels = TickLabels(size=size, axis=self.axis,
                                           axis_direction=self._axis_direction)


        self.major_ticklabels.set(figure = self.axes.figure,
                                  transform=trans,
                                  fontproperties=fontprops)
        self.major_ticklabels.set_pad(major_tick_pad)

        self.minor_ticklabels.set(figure = self.axes.figure,
                                  transform=trans,
                                  fontproperties=fontprops)
        self.minor_ticklabels.set_pad(minor_tick_pad)



    def _get_tick_info(self, tick_iter):
        """
        return ticks_loc_angle, ticklabels_loc_angle_label

        ticks_loc_angle : list of locs and angles for ticks
        ticklabels_loc_angle_label : list of locs, angles and labels for tickslabels
        """
        ticks_loc_angle = []
        ticklabels_loc_angle_label = []

        tick_add_angle = self._tick_add_angle
        ticklabel_add_angle = self._ticklabel_add_angle

        for loc, angle_normal, angle_tangent, label in tick_iter:
            angle_label = angle_tangent  - 90
            angle_label += ticklabel_add_angle

            if np.cos((angle_label - angle_normal)/180.*np.pi) < 0.:
                angle_tick = angle_normal
            else:
                angle_tick = angle_normal + 180

            ticks_loc_angle.append([loc, angle_tick])
            ticklabels_loc_angle_label.append([loc, angle_label, label])

        return ticks_loc_angle, ticklabels_loc_angle_label


    def _update_ticks(self, renderer):


        # set extra pad for major and minor ticklabels:
        # use ticksize of majorticks even for minor ticks. not clear what is best.

        dpi_cor = renderer.points_to_pixels(1.)
        if self.major_ticks.get_visible() and self.major_ticks.get_tick_out():
            self.major_ticklabels._set_external_pad(self.major_ticks._ticksize*dpi_cor)
            self.minor_ticklabels._set_external_pad(self.major_ticks._ticksize*dpi_cor)
        else:
            self.major_ticklabels._set_external_pad(0)
            self.minor_ticklabels._set_external_pad(0)


        majortick_iter,  minortick_iter = \
                self._axis_artist_helper.get_tick_iterators(self.axes)

        tick_loc_angle, ticklabel_loc_angle_label \
                              = self._get_tick_info(majortick_iter)

        self.major_ticks.set_locs_angles(tick_loc_angle)
        self.major_ticklabels.set_locs_angles_labels(ticklabel_loc_angle_label)

        #self.major_ticks.draw(renderer)
        #self.major_ticklabels.draw(renderer)


        # minor ticks
        tick_loc_angle, ticklabel_loc_angle_label \
                              = self._get_tick_info(minortick_iter)

        self.minor_ticks.set_locs_angles(tick_loc_angle)
        self.minor_ticklabels.set_locs_angles_labels(ticklabel_loc_angle_label)

        #self.minor_ticks.draw(renderer)
        #self.minor_ticklabels.draw(renderer)


        #if (self.major_ticklabels.get_visible() or self.minor_ticklabels.get_visible()):
        #    self._draw_offsetText(renderer)

        return self.major_ticklabels.get_window_extents(renderer)


    def _draw_ticks(self, renderer):

        extents = self._update_ticks(renderer)

        self.major_ticks.draw(renderer)
        self.major_ticklabels.draw(renderer)

        self.minor_ticks.draw(renderer)
        self.minor_ticklabels.draw(renderer)


        if (self.major_ticklabels.get_visible() or self.minor_ticklabels.get_visible()):
            self._draw_offsetText(renderer)

        return extents

    def _draw_ticks2(self, renderer):


        # set extra pad for major and minor ticklabels:
        # use ticksize of majorticks even for minor ticks. not clear what is best.

        dpi_cor = renderer.points_to_pixels(1.)
        if self.major_ticks.get_visible() and self.major_ticks.get_tick_out():
            self.major_ticklabels._set_external_pad(self.major_ticks._ticksize*dpi_cor)
            self.minor_ticklabels._set_external_pad(self.major_ticks._ticksize*dpi_cor)
        else:
            self.major_ticklabels._set_external_pad(0)
            self.minor_ticklabels._set_external_pad(0)


        majortick_iter,  minortick_iter = \
                self._axis_artist_helper.get_tick_iterators(self.axes)

        tick_loc_angle, ticklabel_loc_angle_label \
                              = self._get_tick_info(majortick_iter)

        self.major_ticks.set_locs_angles(tick_loc_angle)
        self.major_ticklabels.set_locs_angles_labels(ticklabel_loc_angle_label)

        self.major_ticks.draw(renderer)
        self.major_ticklabels.draw(renderer)


        # minor ticks
        tick_loc_angle, ticklabel_loc_angle_label \
                              = self._get_tick_info(minortick_iter)

        self.minor_ticks.set_locs_angles(tick_loc_angle)
        self.minor_ticklabels.set_locs_angles_labels(ticklabel_loc_angle_label)

        self.minor_ticks.draw(renderer)
        self.minor_ticklabels.draw(renderer)


        if (self.major_ticklabels.get_visible() or self.minor_ticklabels.get_visible()):
            self._draw_offsetText(renderer)

        return self.major_ticklabels.get_window_extents(renderer)




    _offsetText_pos = dict(left=(0, 1, "bottom", "right"),
                           right=(1, 1, "bottom", "left"),
                           bottom=(1, 0, "top", "right"),
                           top=(1, 1, "bottom", "right"))

    def _init_offsetText(self, direction):

        x,y,va,ha = self._offsetText_pos[direction]

        self.offsetText = mtext.Annotation("",
                                           xy=(x,y), xycoords="axes fraction",
                                           xytext=(0,0), textcoords="offset points",
                                           #fontproperties = fp,
                                           color = rcParams['xtick.color'],
                                           verticalalignment=va,
                                           horizontalalignment=ha,
                                           )
        self.offsetText.set_transform(IdentityTransform())
        self.axes._set_artist_props(self.offsetText)


    def _update_offsetText(self):
        self.offsetText.set_text( self.axis.major.formatter.get_offset() )
        self.offsetText.set_size(self.major_ticklabels.get_size())
        offset = self.major_ticklabels.get_pad() + self.major_ticklabels.get_size() + 2.
        self.offsetText.xyann= (0, offset)


    def _draw_offsetText(self, renderer):
        self._update_offsetText()
        self.offsetText.draw(renderer)



    def _init_label(self, **kw):
        # x in axes coords, y in display coords (to be updated at draw
        # time by _update_label_positions)

        labelsize = kw.get("labelsize",
                           rcParams['axes.labelsize'])
        #labelcolor = kw.get("labelcolor",
        #                    rcParams['axes.labelcolor'])
        fontprops = font_manager.FontProperties(
            size=labelsize,
            weight=rcParams['axes.labelweight'])
        textprops = dict(fontproperties = fontprops)
                         #color = labelcolor)

        tr = self._axis_artist_helper.get_axislabel_transform(self.axes) \
             + self.offset_transform

        self.label = AxisLabel(0, 0, "__from_axes__",
                               color = "auto", #rcParams['axes.labelcolor'],
                               fontproperties=fontprops,
                               axis=self.axis,
                               transform=tr,
                               axis_direction=self._axis_direction,
                               )

        self.label.set_figure(self.axes.figure)

        labelpad = kw.get("labelpad", 5)
        self.label.set_pad(labelpad)


    def _update_label(self, renderer):

        if not self.label.get_visible():
            return

        fontprops = font_manager.FontProperties(
            size=rcParams['axes.labelsize'],
            weight=rcParams['axes.labelweight'])

        #pad_points = self.major_tick_pad

        #if abs(self._ticklabel_add_angle - self._axislabel_add_angle)%360 > 90:
        if self._ticklabel_add_angle != self._axislabel_add_angle:
            if (self.major_ticks.get_visible() and not self.major_ticks.get_tick_out()) \
               or \
               (self.minor_ticks.get_visible() and not self.major_ticks.get_tick_out()):
                axislabel_pad = self.major_ticks._ticksize
            else:
                axislabel_pad = 0
        else:
            axislabel_pad = max(self.major_ticklabels._axislabel_pad,
                                self.minor_ticklabels._axislabel_pad)


        #label_offset =  axislabel_pad + self.LABELPAD

        #self.label._set_offset_radius(label_offset)
        self.label._set_external_pad(axislabel_pad)

        xy, angle_tangent = self._axis_artist_helper.get_axislabel_pos_angle(self.axes)
        if xy is None: return

        angle_label = angle_tangent  - 90


        x, y = xy
        self.label._set_ref_angle(angle_label+self._axislabel_add_angle)
        self.label.set(x=x, y=y)


    def _draw_label(self, renderer):
        self._update_label(renderer)
        self.label.draw(renderer)

    def _draw_label2(self, renderer):

        if not self.label.get_visible():
            return

        fontprops = font_manager.FontProperties(
            size=rcParams['axes.labelsize'],
            weight=rcParams['axes.labelweight'])

        #pad_points = self.major_tick_pad

        #if abs(self._ticklabel_add_angle - self._axislabel_add_angle)%360 > 90:
        if self._ticklabel_add_angle != self._axislabel_add_angle:
            if (self.major_ticks.get_visible() and not self.major_ticks.get_tick_out()) \
               or \
               (self.minor_ticks.get_visible() and not self.major_ticks.get_tick_out()):
                axislabel_pad = self.major_ticks._ticksize
            else:
                axislabel_pad = 0
        else:
            axislabel_pad = max(self.major_ticklabels._axislabel_pad,
                                self.minor_ticklabels._axislabel_pad)

        #label_offset =  axislabel_pad + self.LABELPAD

        #self.label._set_offset_radius(label_offset)
        self.label._set_external_pad(axislabel_pad)

        xy, angle_tangent = self._axis_artist_helper.get_axislabel_pos_angle(self.axes)
        if xy is None: return

        angle_label = angle_tangent - 90

        x, y = xy
        self.label._set_ref_angle(angle_label+self._axislabel_add_angle)
        self.label.set(x=x, y=y)
        self.label.draw(renderer)



    def set_label(self, s):
        self.label.set_text(s)



    def get_tightbbox(self, renderer):
        if not self.get_visible(): return

        self._axis_artist_helper.update_lim(self.axes)

        dpi_cor = renderer.points_to_pixels(1.)
        self.dpi_transform.clear().scale(dpi_cor, dpi_cor)


        bb = []

        self._update_ticks(renderer)

        #if self.major_ticklabels.get_visible():
        bb.extend(self.major_ticklabels.get_window_extents(renderer))
        #if self.minor_ticklabels.get_visible():
        bb.extend(self.minor_ticklabels.get_window_extents(renderer))


        self._update_label(renderer)

        #if self.label.get_visible():
        bb.append(self.label.get_window_extent(renderer))
        bb.append(self.offsetText.get_window_extent(renderer))

        bb = [b for b in bb if b and (b.width!=0 or b.height!=0)]
        if bb:
            _bbox = Bbox.union(bb)
            return _bbox
        else:
            return None

        #self._draw_line(renderer)

        #self._draw_ticks(renderer)

        #self._draw_offsetText(renderer)
        #self._draw_label(renderer)



    @allow_rasterization
    def draw(self, renderer):
        'Draw the axis lines, tick lines and labels'

        if not self.get_visible(): return

        renderer.open_group(__name__)

        self._axis_artist_helper.update_lim(self.axes)

        dpi_cor = renderer.points_to_pixels(1.)
        self.dpi_transform.clear().scale(dpi_cor, dpi_cor)


        self._draw_ticks(renderer)

        self._draw_line(renderer)

        #self._draw_offsetText(renderer)
        self._draw_label(renderer)

        renderer.close_group(__name__)

    #def get_ticklabel_extents(self, renderer):
    #    pass

    def toggle(self, all=None, ticks=None, ticklabels=None, label=None):
        """
        Toggle visibility of ticks, ticklabels, and (axis) label.
        To turn all off, ::

          axis.toggle(all=False)

        To turn all off but ticks on ::

          axis.toggle(all=False, ticks=True)

        To turn all on but (axis) label off ::

          axis.toggle(all=True, label=False))

        """
        if all:
            _ticks, _ticklabels, _label = True, True, True
        elif all is not None:
            _ticks, _ticklabels, _label = False, False, False
        else:
            _ticks, _ticklabels, _label = None, None, None

        if ticks is not None:
            _ticks = ticks
        if ticklabels is not None:
            _ticklabels = ticklabels
        if label is not None:
            _label = label

        if _ticks is not None:
            self.major_ticks.set_visible(_ticks)
            self.minor_ticks.set_visible(_ticks)
        if _ticklabels is not None:
            self.major_ticklabels.set_visible(_ticklabels)
            self.minor_ticklabels.set_visible(_ticklabels)
        if _label is not None:
            self.label.set_visible(_label)
