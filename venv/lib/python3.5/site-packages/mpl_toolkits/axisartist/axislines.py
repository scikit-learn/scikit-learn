"""
Axislines includes modified implementation of the Axes class. The
biggest difference is that the artists responsible for drawing the axis spine,
ticks, ticklabels and axis labels are separated out from mpl's Axis
class. Originally, this change was motivated to support curvilinear
grid. Here are a few reasons that I came up with a new axes class:


 * "top" and "bottom" x-axis (or "left" and "right" y-axis) can have
   different ticks (tick locations and labels). This is not possible
   with the current mpl, although some twin axes trick can help.

 * Curvilinear grid.

 * angled ticks.

In the new axes class, xaxis and yaxis is set to not visible by
default, and new set of artist (AxisArtist) are defined to draw axis
line, ticks, ticklabels and axis label. Axes.axis attribute serves as
a dictionary of these artists, i.e., ax.axis["left"] is a AxisArtist
instance responsible to draw left y-axis. The default Axes.axis contains
"bottom", "left", "top" and "right".

AxisArtist can be considered as a container artist and
has following children artists which will draw ticks, labels, etc.

 * line
 * major_ticks, major_ticklabels
 * minor_ticks, minor_ticklabels
 * offsetText
 * label

Note that these are separate artists from Axis class of the
original mpl, thus most of tick-related command in the original mpl
won't work, although some effort has made to work with. For example,
color and markerwidth of the ax.axis["bottom"].major_ticks will follow
those of Axes.xaxis unless explicitly specified.

In addition to AxisArtist, the Axes will have *gridlines* attribute,
which obviously draws grid lines. The gridlines needs to be separated
from the axis as some gridlines can never pass any axis.

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import warnings

import numpy as np

from matplotlib import rcParams
import matplotlib.artist as martist
import matplotlib.axes as maxes
from matplotlib.path import Path
from matplotlib.transforms import Bbox
from .axisline_style import AxislineStyle
from .axis_artist import AxisArtist, GridlinesCollection


class AxisArtistHelper(object):
    """
    AxisArtistHelper should define
    following method with given APIs. Note that the first axes argument
    will be axes attribute of the caller artist.::


        # LINE (spinal line?)

        def get_line(self, axes):
            # path : Path
            return path

        def get_line_transform(self, axes):
            # ...
            # trans : transform
            return trans

        # LABEL

        def get_label_pos(self, axes):
            # x, y : position
            return (x, y), trans


        def get_label_offset_transform(self, \
                axes,
                pad_points, fontprops, renderer,
                bboxes,
                ):
            # va : vertical alignment
            # ha : horizontal alignment
            # a : angle
            return trans, va, ha, a

        # TICK

        def get_tick_transform(self, axes):
            return trans

        def get_tick_iterators(self, axes):
            # iter : iterable object that yields (c, angle, l) where
            # c, angle, l is position, tick angle, and label

            return iter_major, iter_minor


    """

    class _Base(object):
        """
        Base class for axis helper.
        """
        def __init__(self):
            """
            """
            self.delta1, self.delta2 = 0.00001, 0.00001

        def update_lim(self, axes):
            pass


    class Fixed(_Base):
        """
        Helper class for a fixed (in the axes coordinate) axis.
        """

        _default_passthru_pt = dict(left=(0, 0),
                                    right=(1, 0),
                                    bottom=(0, 0),
                                    top=(0, 1))

        def __init__(self,
                     loc, nth_coord=None,
                     ):
            """
            nth_coord = along which coordinate value varies
            in 2d, nth_coord = 0 ->  x axis, nth_coord = 1 -> y axis
            """

            self._loc = loc

            if loc not in ["left", "right", "bottom", "top"]:
                raise ValueError("%s" % loc)

            if nth_coord is None:
                if loc in ["left", "right"]:
                    nth_coord = 1
                elif loc in ["bottom", "top"]:
                    nth_coord = 0

            self.nth_coord = nth_coord

            super(AxisArtistHelper.Fixed, self).__init__()

            self.passthru_pt = self._default_passthru_pt[loc]



            _verts = np.array([[0., 0.],
                               [1., 1.]])
            fixed_coord = 1-nth_coord
            _verts[:,fixed_coord] = self.passthru_pt[fixed_coord]

            # axis line in transAxes
            self._path = Path(_verts)


        def get_nth_coord(self):
            return self.nth_coord

        # LINE

        def get_line(self, axes):
            return self._path

        def get_line_transform(self, axes):
            return axes.transAxes

        # LABEL

        def get_axislabel_transform(self, axes):
            return axes.transAxes

        def get_axislabel_pos_angle(self, axes):
            """
            label reference position in transAxes.

            get_label_transform() returns a transform of (transAxes+offset)
            """
            loc = self._loc
            pos, angle_tangent = dict(left=((0., 0.5), 90),
                                      right=((1., 0.5), 90),
                                      bottom=((0.5, 0.), 0),
                                      top=((0.5, 1.), 0))[loc]

            return pos, angle_tangent



        # TICK

        def get_tick_transform(self, axes):
            trans_tick = [axes.get_xaxis_transform(),
                          axes.get_yaxis_transform()][self.nth_coord]

            return trans_tick


    class Floating(_Base):
        def __init__(self, nth_coord,
                     value):

            self.nth_coord = nth_coord

            self._value = value

            super(AxisArtistHelper.Floating,
                  self).__init__()


        def get_nth_coord(self):
            return self.nth_coord

        def get_line(self, axes):
            raise RuntimeError("get_line method should be defined by the derived class")




class AxisArtistHelperRectlinear(object):

    class Fixed(AxisArtistHelper.Fixed):

        def __init__(self, axes, loc, nth_coord=None):
            """
            nth_coord = along which coordinate value varies
            in 2d, nth_coord = 0 ->  x axis, nth_coord = 1 -> y axis
            """
            super(AxisArtistHelperRectlinear.Fixed, self).__init__(
                loc, nth_coord)
            self.axis = [axes.xaxis, axes.yaxis][self.nth_coord]

        # TICK

        def get_tick_iterators(self, axes):
            """tick_loc, tick_angle, tick_label"""

            loc = self._loc

            if loc in ["bottom", "top"]:
                angle_normal, angle_tangent = 90, 0
            else:
                angle_normal, angle_tangent = 0, 90

            major = self.axis.major
            majorLocs = major.locator()
            major.formatter.set_locs(majorLocs)
            majorLabels = [major.formatter(val, i) for i, val in enumerate(majorLocs)]

            minor = self.axis.minor
            minorLocs = minor.locator()
            minor.formatter.set_locs(minorLocs)
            minorLabels = [minor.formatter(val, i) for i, val in enumerate(minorLocs)]

            trans_tick = self.get_tick_transform(axes)

            tr2ax = trans_tick + axes.transAxes.inverted()

            def _f(locs, labels):
                for x, l in zip(locs, labels):

                    c = list(self.passthru_pt) # copy
                    c[self.nth_coord] = x

                    # check if the tick point is inside axes
                    c2 = tr2ax.transform_point(c)
                    #delta=0.00001
                    if 0. -self.delta1<= c2[self.nth_coord] <= 1.+self.delta2:
                        yield c, angle_normal, angle_tangent, l

            return _f(majorLocs, majorLabels), _f(minorLocs, minorLabels)



    class Floating(AxisArtistHelper.Floating):
        def __init__(self, axes, nth_coord,
                     passingthrough_point, axis_direction="bottom"):
            super(AxisArtistHelperRectlinear.Floating, self).__init__(
                nth_coord, passingthrough_point)
            self._axis_direction = axis_direction
            self.axis = [axes.xaxis, axes.yaxis][self.nth_coord]

        def get_line(self, axes):
            _verts = np.array([[0., 0.],
                               [1., 1.]])

            fixed_coord = 1-self.nth_coord
            trans_passingthrough_point = axes.transData + axes.transAxes.inverted()
            p = trans_passingthrough_point.transform_point([self._value,
                                                            self._value])
            _verts[:,fixed_coord] = p[fixed_coord]

            return Path(_verts)

        def get_line_transform(self, axes):
            return axes.transAxes

        def get_axislabel_transform(self, axes):
            return axes.transAxes

        def get_axislabel_pos_angle(self, axes):
            """
            label reference position in transAxes.

            get_label_transform() returns a transform of (transAxes+offset)
            """
            loc = self._axis_direction
            #angle = dict(left=0,
            #             right=0,
            #             bottom=.5*np.pi,
            #             top=.5*np.pi)[loc]

            if self.nth_coord == 0:
                angle = 0
            else:
                angle = 90

            _verts = [0.5, 0.5]

            fixed_coord = 1-self.nth_coord
            trans_passingthrough_point = axes.transData + axes.transAxes.inverted()
            p = trans_passingthrough_point.transform_point([self._value,
                                                            self._value])
            _verts[fixed_coord] = p[fixed_coord]
            if not (0. <= _verts[fixed_coord] <= 1.):
                return None, None
            else:
                return _verts, angle



        def get_tick_transform(self, axes):
            return axes.transData


        def get_tick_iterators(self, axes):
            """tick_loc, tick_angle, tick_label"""

            loc = self._axis_direction

            if loc in ["bottom", "top"]:
                angle_normal, angle_tangent = 90, 0
            else:
                angle_normal, angle_tangent = 0, 90

            if self.nth_coord == 0:
                angle_normal, angle_tangent = 90, 0
            else:
                angle_normal, angle_tangent = 0, 90

            #angle = 90 - 90 * self.nth_coord

            major = self.axis.major
            majorLocs = major.locator()
            major.formatter.set_locs(majorLocs)
            majorLabels = [major.formatter(val, i) for i, val in enumerate(majorLocs)]

            minor = self.axis.minor
            minorLocs = minor.locator()
            minor.formatter.set_locs(minorLocs)
            minorLabels = [minor.formatter(val, i) for i, val in enumerate(minorLocs)]

            tr2ax = axes.transData + axes.transAxes.inverted()

            def _f(locs, labels):
                for x, l in zip(locs, labels):

                    c = [self._value, self._value]
                    c[self.nth_coord] = x
                    c1, c2 = tr2ax.transform_point(c)
                    if 0. <= c1 <= 1. and 0. <= c2 <= 1.:
                        if 0. - self.delta1 <= [c1, c2][self.nth_coord] <= 1. + self.delta2:
                            yield c, angle_normal, angle_tangent, l

            return _f(majorLocs, majorLabels), _f(minorLocs, minorLabels)





class GridHelperBase(object):

    def __init__(self):
        self._force_update = True
        self._old_limits = None
        super(GridHelperBase, self).__init__()


    def update_lim(self, axes):
        x1, x2 = axes.get_xlim()
        y1, y2 = axes.get_ylim()

        if self._force_update or self._old_limits != (x1, x2, y1, y2):
            self._update(x1, x2, y1, y2)
            self._force_update = False
            self._old_limits = (x1, x2, y1, y2)


    def _update(self, x1, x2, y1, y2):
        pass


    def invalidate(self):
        self._force_update = True

    def valid(self):
        return not self._force_update


    def get_gridlines(self, which, axis):
        """
        Return list of grid lines as a list of paths (list of points).

        *which* : "major" or "minor"
        *axis* : "both", "x" or "y"
        """
        return []

    def new_gridlines(self, ax):
        """
        Create and return a new GridlineCollection instance.

        *which* : "major" or "minor"
        *axis* : "both", "x" or "y"

        """
        gridlines = GridlinesCollection(None, transform=ax.transData,
                                        colors=rcParams['grid.color'],
                                        linestyles=rcParams['grid.linestyle'],
                                        linewidths=rcParams['grid.linewidth'])
        ax._set_artist_props(gridlines)
        gridlines.set_grid_helper(self)

        ax.axes._set_artist_props(gridlines)
        # gridlines.set_clip_path(self.axes.patch)
        # set_clip_path need to be deferred after Axes.cla is completed.
        # It is done inside the cla.

        return gridlines


class GridHelperRectlinear(GridHelperBase):


    def __init__(self, axes):

        super(GridHelperRectlinear, self).__init__()
        self.axes = axes



    def new_fixed_axis(self, loc,
                       nth_coord=None,
                       axis_direction=None,
                       offset=None,
                       axes=None,
                       ):

        if axes is None:
            warnings.warn("'new_fixed_axis' explicitly requires the axes keyword.")
            axes = self.axes

        _helper = AxisArtistHelperRectlinear.Fixed(axes, loc, nth_coord)

        if axis_direction is None:
            axis_direction = loc
        axisline = AxisArtist(axes, _helper, offset=offset,
                              axis_direction=axis_direction,
                              )

        return axisline


    def new_floating_axis(self, nth_coord, value,
                          axis_direction="bottom",
                          axes=None,
                          ):

        if axes is None:
            warnings.warn(
                "'new_floating_axis' explicitly requires the axes keyword.")
            axes = self.axes

        passthrough_point = (value, value)
        transform = axes.transData

        _helper = AxisArtistHelperRectlinear.Floating(
            axes, nth_coord, value, axis_direction)

        axisline = AxisArtist(axes, _helper)

        axisline.line.set_clip_on(True)
        axisline.line.set_clip_box(axisline.axes.bbox)
        return axisline


    def get_gridlines(self, which="major", axis="both"):
        """
        return list of gridline coordinates in data coordinates.

        *which* : "major" or "minor"
        *axis* : "both", "x" or "y"
        """

        gridlines = []


        if axis in ["both", "x"]:
            locs = []
            y1, y2 = self.axes.get_ylim()
            #if self.axes.xaxis._gridOnMajor:
            if which in ["both", "major"]:
                locs.extend(self.axes.xaxis.major.locator())
            #if self.axes.xaxis._gridOnMinor:
            if which in ["both", "minor"]:
                locs.extend(self.axes.xaxis.minor.locator())

            for x in locs:
                gridlines.append([[x, x], [y1, y2]])


        if axis in ["both", "y"]:
            x1, x2 = self.axes.get_xlim()
            locs = []
            if self.axes.yaxis._gridOnMajor:
            #if which in ["both", "major"]:
                locs.extend(self.axes.yaxis.major.locator())
            if self.axes.yaxis._gridOnMinor:
            #if which in ["both", "minor"]:
                locs.extend(self.axes.yaxis.minor.locator())

            for y in locs:
                gridlines.append([[x1, x2], [y, y]])

        return gridlines






class SimpleChainedObjects(object):
    def __init__(self, objects):
        self._objects = objects

    def __getattr__(self, k):
        _a = SimpleChainedObjects([getattr(a, k) for a in self._objects])
        return _a

    def __call__(self, *kl, **kwargs):
        for m in self._objects:
            m(*kl, **kwargs)


class Axes(maxes.Axes):

    class AxisDict(dict):
        def __init__(self, axes):
            self.axes = axes
            super(Axes.AxisDict, self).__init__()

        def __getitem__(self, k):
            if isinstance(k, tuple):
                r = SimpleChainedObjects([dict.__getitem__(self, k1) for k1 in k])
                return r
            elif isinstance(k, slice):
                if k.start == None and k.stop == None and k.step == None:
                    r = SimpleChainedObjects(list(six.itervalues(self)))
                    return r
                else:
                    raise ValueError("Unsupported slice")
            else:
                return dict.__getitem__(self, k)

        def __call__(self, *v, **kwargs):
            return maxes.Axes.axis(self.axes, *v, **kwargs)


    def __init__(self, *kl, **kw):


        helper = kw.pop("grid_helper", None)

        self._axisline_on = True

        if helper:
            self._grid_helper = helper
        else:
            self._grid_helper = GridHelperRectlinear(self)

        super(Axes, self).__init__(*kl, **kw)

        self.toggle_axisline(True)


    def toggle_axisline(self, b=None):
        if b is None:
            b = not self._axisline_on
        if b:
            self._axisline_on = True
            for s in self.spines.values():
                s.set_visible(False)
            self.xaxis.set_visible(False)
            self.yaxis.set_visible(False)
        else:
            self._axisline_on = False
            for s in self.spines.values():
                s.set_visible(True)
            self.xaxis.set_visible(True)
            self.yaxis.set_visible(True)


    def _init_axis(self):
        super(Axes, self)._init_axis()


    def _init_axis_artists(self, axes=None):
        if axes is None:
            axes = self

        self._axislines = self.AxisDict(self)
        new_fixed_axis = self.get_grid_helper().new_fixed_axis
        for loc in ["bottom", "top", "left", "right"]:
            self._axislines[loc] = new_fixed_axis(loc=loc, axes=axes,
                                                  axis_direction=loc)

        for axisline in [self._axislines["top"], self._axislines["right"]]:
            axisline.label.set_visible(False)
            axisline.major_ticklabels.set_visible(False)
            axisline.minor_ticklabels.set_visible(False)

    @property
    def axis(self):
        return self._axislines

    def new_gridlines(self, grid_helper=None):
        """
        Create and return a new GridlineCollection instance.

        *which* : "major" or "minor"
        *axis* : "both", "x" or "y"

        """
        if grid_helper is None:
            grid_helper = self.get_grid_helper()

        gridlines = grid_helper.new_gridlines(self)

        return gridlines


    def _init_gridlines(self, grid_helper=None):
        # It is done inside the cla.
        gridlines = self.new_gridlines(grid_helper)

        self.gridlines = gridlines

    def cla(self):
        # gridlines need to b created before cla() since cla calls grid()

        self._init_gridlines()
        super(Axes, self).cla()

        # the clip_path should be set after Axes.cla() since that's
        # when a patch is created.
        self.gridlines.set_clip_path(self.axes.patch)

        self._init_axis_artists()

    def get_grid_helper(self):
        return self._grid_helper


    def grid(self, b=None, which='major', axis="both", **kwargs):
        """
        Toggle the gridlines, and optionally set the properties of the lines.
        """
        # their are some discrepancy between the behavior of grid in
        # axes_grid and the original mpl's grid, because axes_grid
        # explicitly set the visibility of the gridlines.

        super(Axes, self).grid(b, which=which, axis=axis, **kwargs)
        if not self._axisline_on:
            return

        if b is None:

            if self.axes.xaxis._gridOnMinor or self.axes.xaxis._gridOnMajor or \
                   self.axes.yaxis._gridOnMinor or self.axes.yaxis._gridOnMajor:
                b=True
            else:
                b=False

        self.gridlines.set_which(which)
        self.gridlines.set_axis(axis)
        self.gridlines.set_visible(b)

        if len(kwargs):
            martist.setp(self.gridlines, **kwargs)

    def get_children(self):
        if self._axisline_on:
            children = list(six.itervalues(self._axislines)) + [self.gridlines]
        else:
            children = []
        children.extend(super(Axes, self).get_children())
        return children

    def invalidate_grid_helper(self):
        self._grid_helper.invalidate()


    def new_fixed_axis(self, loc, offset=None):
        gh = self.get_grid_helper()
        axis = gh.new_fixed_axis(loc,
                                 nth_coord=None,
                                 axis_direction=None,
                                 offset=offset,
                                 axes=self,
                                 )
        return axis


    def new_floating_axis(self, nth_coord, value,
                          axis_direction="bottom",
                          ):
        gh = self.get_grid_helper()
        axis = gh.new_floating_axis(nth_coord, value,
                                    axis_direction=axis_direction,
                                    axes=self)
        return axis



    def draw(self, renderer, inframe=False):

        if not self._axisline_on:
            super(Axes, self).draw(renderer, inframe)
            return

        orig_artists = self.artists
        self.artists = self.artists + list(self._axislines.values()) + [self.gridlines]

        super(Axes, self).draw(renderer, inframe)

        self.artists = orig_artists


    def get_tightbbox(self, renderer, call_axes_locator=True):

        bb0 = super(Axes, self).get_tightbbox(renderer, call_axes_locator)

        if not self._axisline_on:
            return bb0

        bb = [bb0]

        for axisline in list(six.itervalues(self._axislines)):
            if not axisline.get_visible():
                continue

            bb.append(axisline.get_tightbbox(renderer))
            # if axisline.label.get_visible():
            #     bb.append(axisline.label.get_window_extent(renderer))


            # if axisline.major_ticklabels.get_visible():
            #     bb.extend(axisline.major_ticklabels.get_window_extents(renderer))
            # if axisline.minor_ticklabels.get_visible():
            #     bb.extend(axisline.minor_ticklabels.get_window_extents(renderer))
            # if axisline.major_ticklabels.get_visible() or \
            #    axisline.minor_ticklabels.get_visible():
            #     bb.append(axisline.offsetText.get_window_extent(renderer))

        #bb.extend([c.get_window_extent(renderer) for c in artists \
        #           if c.get_visible()])

        _bbox = Bbox.union([b for b in bb if b and (b.width!=0 or b.height!=0)])

        return _bbox




Subplot = maxes.subplot_class_factory(Axes)

class AxesZero(Axes):
    def __init__(self, *kl, **kw):

        super(AxesZero, self).__init__(*kl, **kw)


    def _init_axis_artists(self):
        super(AxesZero, self)._init_axis_artists()

        new_floating_axis = self._grid_helper.new_floating_axis
        xaxis_zero = new_floating_axis(nth_coord=0,
                                       value=0.,
                                       axis_direction="bottom",
                                       axes=self)

        xaxis_zero.line.set_clip_path(self.patch)
        xaxis_zero.set_visible(False)
        self._axislines["xzero"] = xaxis_zero

        yaxis_zero = new_floating_axis(nth_coord=1,
                                       value=0.,
                                       axis_direction="left",
                                       axes=self)


        yaxis_zero.line.set_clip_path(self.patch)
        yaxis_zero.set_visible(False)
        self._axislines["yzero"] = yaxis_zero

SubplotZero = maxes.subplot_class_factory(AxesZero)
