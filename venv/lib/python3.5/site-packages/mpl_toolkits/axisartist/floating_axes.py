"""
An experimental support for curvilinear grid.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves import zip

# TODO :
# see if tick_iterator method can be simplified by reusing the parent method.

import numpy as np

from matplotlib.transforms import Affine2D, IdentityTransform
from . import grid_helper_curvelinear
from .axislines import AxisArtistHelper, GridHelperBase
from .axis_artist import AxisArtist
from .grid_finder import GridFinder


class FloatingAxisArtistHelper(grid_helper_curvelinear.FloatingAxisArtistHelper):
    pass


class FixedAxisArtistHelper(grid_helper_curvelinear.FloatingAxisArtistHelper):

    def __init__(self, grid_helper, side, nth_coord_ticks=None):
        """
        nth_coord = along which coordinate value varies.
         nth_coord = 0 ->  x axis, nth_coord = 1 -> y axis
        """

        value, nth_coord = grid_helper.get_data_boundary(side) # return v= 0 , nth=1, extremes of the other coordinate.
        super(FixedAxisArtistHelper, self).__init__(grid_helper,
                                                    nth_coord,
                                                    value,
                                                    axis_direction=side,
                                                    )
        #self.grid_helper = grid_helper
        if nth_coord_ticks is None:
            nth_coord_ticks = nth_coord
        self.nth_coord_ticks = nth_coord_ticks

        self.value = value
        self.grid_helper = grid_helper
        self._side = side


    def update_lim(self, axes):
        self.grid_helper.update_lim(axes)

        self.grid_info = self.grid_helper.grid_info



    def get_axislabel_pos_angle(self, axes):

        extremes = self.grid_info["extremes"]

        if self.nth_coord == 0:
            xx0 = self.value
            yy0 = (extremes[2]+extremes[3])/2.
            dxx, dyy = 0., abs(extremes[2]-extremes[3])/1000.
        elif self.nth_coord == 1:
            xx0 = (extremes[0]+extremes[1])/2.
            yy0 = self.value
            dxx, dyy = abs(extremes[0]-extremes[1])/1000., 0.

        grid_finder = self.grid_helper.grid_finder
        xx1, yy1 = grid_finder.transform_xy([xx0], [yy0])

        trans_passingthrough_point = axes.transData + axes.transAxes.inverted()
        p = trans_passingthrough_point.transform_point([xx1[0], yy1[0]])


        if (0. <= p[0] <= 1.) and (0. <= p[1] <= 1.):
            xx1c, yy1c = axes.transData.transform_point([xx1[0], yy1[0]])
            xx2, yy2 = grid_finder.transform_xy([xx0+dxx], [yy0+dyy])
            xx2c, yy2c = axes.transData.transform_point([xx2[0], yy2[0]])

            return (xx1c, yy1c), np.arctan2(yy2c-yy1c, xx2c-xx1c)/np.pi*180.
        else:
            return None, None



    def get_tick_transform(self, axes):
        return IdentityTransform() #axes.transData

    def get_tick_iterators(self, axes):
        """tick_loc, tick_angle, tick_label, (optionally) tick_label"""


        grid_finder = self.grid_helper.grid_finder

        lat_levs, lat_n, lat_factor = self.grid_info["lat_info"]
        lon_levs, lon_n, lon_factor = self.grid_info["lon_info"]

        lon_levs, lat_levs = np.asarray(lon_levs), np.asarray(lat_levs)
        if lat_factor is not None:
            yy0 = lat_levs / lat_factor
            dy = 0.001 / lat_factor
        else:
            yy0 = lat_levs
            dy = 0.001

        if lon_factor is not None:
            xx0 = lon_levs / lon_factor
            dx = 0.001 / lon_factor
        else:
            xx0 = lon_levs
            dx = 0.001

        _extremes = self.grid_helper._extremes
        xmin, xmax = sorted(_extremes[:2])
        ymin, ymax = sorted(_extremes[2:])
        if self.nth_coord == 0:
            mask = (ymin <= yy0) & (yy0 <= ymax)
            yy0 = yy0[mask]
        elif self.nth_coord == 1:
            mask = (xmin <= xx0) & (xx0 <= xmax)
            xx0 = xx0[mask]

        def transform_xy(x, y):
            x1, y1 = grid_finder.transform_xy(x, y)
            x2y2 = axes.transData.transform(np.array([x1, y1]).transpose())
            x2, y2 = x2y2.transpose()
            return x2, y2

        # find angles
        if self.nth_coord == 0:
            xx0 = np.empty_like(yy0)
            xx0.fill(self.value)

            #yy0_ = yy0.copy()

            xx1, yy1 = transform_xy(xx0, yy0)

            xx00 = xx0.astype(float, copy=True)
            xx00[xx0+dx>xmax] -= dx
            xx1a, yy1a = transform_xy(xx00, yy0)
            xx1b, yy1b = transform_xy(xx00+dx, yy0)

            yy00 = yy0.astype(float, copy=True)
            yy00[yy0+dy>ymax] -= dy
            xx2a, yy2a = transform_xy(xx0, yy00)
            xx2b, yy2b = transform_xy(xx0, yy00+dy)

            labels = self.grid_info["lat_labels"]
            labels = [l for l, m in zip(labels, mask) if m]

        elif self.nth_coord == 1:
            yy0 = np.empty_like(xx0)
            yy0.fill(self.value)

            #xx0_ = xx0.copy()
            xx1, yy1 = transform_xy(xx0, yy0)


            yy00 = yy0.astype(float, copy=True)
            yy00[yy0+dy>ymax] -= dy
            xx1a, yy1a = transform_xy(xx0, yy00)
            xx1b, yy1b = transform_xy(xx0, yy00+dy)

            xx00 = xx0.astype(float, copy=True)
            xx00[xx0+dx>xmax] -= dx
            xx2a, yy2a = transform_xy(xx00, yy0)
            xx2b, yy2b = transform_xy(xx00+dx, yy0)

            labels = self.grid_info["lon_labels"]
            labels = [l for l, m in zip(labels, mask) if m]


        def f1():
            dd = np.arctan2(yy1b-yy1a, xx1b-xx1a) # angle normal
            dd2 = np.arctan2(yy2b-yy2a, xx2b-xx2a) # angle tangent
            mm = ((yy1b-yy1a)==0.) & ((xx1b-xx1a)==0.) # mask where dd1 is not defined
            dd[mm] = dd2[mm] + np.pi / 2

            #dd += np.pi
            #dd = np.arctan2(xx2-xx1, angle_tangent-yy1)
            trans_tick = self.get_tick_transform(axes)
            tr2ax = trans_tick + axes.transAxes.inverted()
            for x, y, d, d2, lab in zip(xx1, yy1, dd, dd2, labels):
                c2 = tr2ax.transform_point((x, y))
                delta=0.00001
                if (0. -delta<= c2[0] <= 1.+delta) and \
                       (0. -delta<= c2[1] <= 1.+delta):
                    d1 = d/3.14159*180.
                    d2 = d2/3.14159*180.
                    #_mod = (d2-d1+180)%360
                    #if _mod < 180:
                    #    d1 += 180
                    ##_div, _mod = divmod(d2-d1, 360)
                    yield [x, y], d1, d2, lab
                    #, d2/3.14159*180.+da)

        return f1(), iter([])

    def get_line_transform(self, axes):
        return axes.transData

    def get_line(self, axes):

        self.update_lim(axes)
        from matplotlib.path import Path
        k, v = dict(left=("lon_lines0", 0),
                    right=("lon_lines0", 1),
                    bottom=("lat_lines0", 0),
                    top=("lat_lines0", 1))[self._side]

        xx, yy = self.grid_info[k][v]
        return Path(np.column_stack([xx, yy]))



from .grid_finder import ExtremeFinderSimple

class ExtremeFinderFixed(ExtremeFinderSimple):
    def __init__(self, extremes):
        self._extremes = extremes

    def __call__(self, transform_xy, x1, y1, x2, y2):
        """
        get extreme values.

        x1, y1, x2, y2 in image coordinates (0-based)
        nx, ny : number of division in each axis
        """
        #lon_min, lon_max, lat_min, lat_max = self._extremes
        return self._extremes



class GridHelperCurveLinear(grid_helper_curvelinear.GridHelperCurveLinear):

    def __init__(self, aux_trans, extremes,
                 grid_locator1=None,
                 grid_locator2=None,
                 tick_formatter1=None,
                 tick_formatter2=None):
        """
        aux_trans : a transform from the source (curved) coordinate to
        target (rectilinear) coordinate. An instance of MPL's Transform
        (inverse transform should be defined) or a tuple of two callable
        objects which defines the transform and its inverse. The callables
        need take two arguments of array of source coordinates and
        should return two target coordinates:
          e.g., x2, y2 = trans(x1, y1)
        """

        self._old_values = None

        self._extremes = extremes
        extreme_finder = ExtremeFinderFixed(extremes)

        super(GridHelperCurveLinear, self).__init__(aux_trans,
                                                    extreme_finder,
                                                    grid_locator1=grid_locator1,
                                                    grid_locator2=grid_locator2,
                                                    tick_formatter1=tick_formatter1,
                                                    tick_formatter2=tick_formatter2)


    # def update_grid_finder(self, aux_trans=None, **kw):

    #     if aux_trans is not None:
    #         self.grid_finder.update_transform(aux_trans)

    #     self.grid_finder.update(**kw)
    #     self.invalidate()


    # def _update(self, x1, x2, y1, y2):
    #     "bbox in 0-based image coordinates"
    #     # update wcsgrid

    #     if self.valid() and self._old_values == (x1, x2, y1, y2):
    #         return

    #     self._update_grid(x1, y1, x2, y2)

    #     self._old_values = (x1, x2, y1, y2)

    #     self._force_update = False


    def get_data_boundary(self, side):
        """
        return v= 0 , nth=1
        """
        lon1, lon2, lat1, lat2 = self._extremes
        return dict(left=(lon1, 0),
                    right=(lon2, 0),
                    bottom=(lat1, 1),
                    top=(lat2, 1))[side]


    def new_fixed_axis(self, loc,
                       nth_coord=None,
                       axis_direction=None,
                       offset=None,
                       axes=None):

        if axes is None:
            axes = self.axes

        if axis_direction is None:
            axis_direction = loc

        _helper = FixedAxisArtistHelper(self, loc,
                                        nth_coord_ticks=nth_coord)


        axisline = AxisArtist(axes, _helper, axis_direction=axis_direction)
        axisline.line.set_clip_on(True)
        axisline.line.set_clip_box(axisline.axes.bbox)


        return axisline


    # new_floating_axis will inherit the grid_helper's extremes.

    # def new_floating_axis(self, nth_coord,
    #                       value,
    #                       axes=None,
    #                       axis_direction="bottom"
    #                       ):

    #     axis = super(GridHelperCurveLinear,
    #                  self).new_floating_axis(nth_coord,
    #                                          value, axes=axes,
    #                                          axis_direction=axis_direction)

    #     # set extreme values of the axis helper
    #     if nth_coord == 1:
    #         axis.get_helper().set_extremes(*self._extremes[:2])
    #     elif nth_coord == 0:
    #         axis.get_helper().set_extremes(*self._extremes[2:])

    #     return axis


    def _update_grid(self, x1, y1, x2, y2):

        #self.grid_info = self.grid_finder.get_grid_info(x1, y1, x2, y2)

        if self.grid_info is None:
            self.grid_info = dict()

        grid_info = self.grid_info

        grid_finder = self.grid_finder
        extremes = grid_finder.extreme_finder(grid_finder.inv_transform_xy,
                                              x1, y1, x2, y2)

        lon_min, lon_max = sorted(extremes[:2])
        lat_min, lat_max = sorted(extremes[2:])
        lon_levs, lon_n, lon_factor = \
                  grid_finder.grid_locator1(lon_min, lon_max)
        lat_levs, lat_n, lat_factor = \
                  grid_finder.grid_locator2(lat_min, lat_max)
        grid_info["extremes"] = lon_min, lon_max, lat_min, lat_max #extremes

        grid_info["lon_info"] = lon_levs, lon_n, lon_factor
        grid_info["lat_info"] = lat_levs, lat_n, lat_factor

        grid_info["lon_labels"] = grid_finder.tick_formatter1("bottom",
                                                              lon_factor,
                                                              lon_levs)

        grid_info["lat_labels"] = grid_finder.tick_formatter2("bottom",
                                                              lat_factor,
                                                              lat_levs)

        if lon_factor is None:
            lon_values = np.asarray(lon_levs[:lon_n])
        else:
            lon_values = np.asarray(lon_levs[:lon_n]/lon_factor)
        if lat_factor is None:
            lat_values = np.asarray(lat_levs[:lat_n])
        else:
            lat_values = np.asarray(lat_levs[:lat_n]/lat_factor)

        lon_values0 = lon_values[(lon_min<lon_values) & (lon_values<lon_max)]
        lat_values0 = lat_values[(lat_min<lat_values) & (lat_values<lat_max)]
        lon_lines, lat_lines = grid_finder._get_raw_grid_lines(lon_values0,
                                                               lat_values0,
                                                               lon_min, lon_max,
                                                               lat_min, lat_max)


        grid_info["lon_lines"] = lon_lines
        grid_info["lat_lines"] = lat_lines


        lon_lines, lat_lines = grid_finder._get_raw_grid_lines(extremes[:2],
                                                               extremes[2:],
                                                               *extremes)
        #lon_min, lon_max,
        #                                                       lat_min, lat_max)


        grid_info["lon_lines0"] = lon_lines
        grid_info["lat_lines0"] = lat_lines



    def get_gridlines(self, which="major", axis="both"):
        grid_lines = []
        if axis in ["both", "x"]:
            for gl in self.grid_info["lon_lines"]:
                grid_lines.extend([gl])
        if axis in ["both", "y"]:
            for gl in self.grid_info["lat_lines"]:
                grid_lines.extend([gl])

        return grid_lines


    def get_boundary(self):
        """
        return Nx2 array of x,y coordinate of the boundary
        """
        x0, x1, y0, y1 = self._extremes
        tr = self._aux_trans
        xx = np.linspace(x0, x1, 100)
        yy0, yy1 = np.empty_like(xx), np.empty_like(xx)
        yy0.fill(y0)
        yy1.fill(y1)

        yy = np.linspace(y0, y1, 100)
        xx0, xx1 = np.empty_like(yy), np.empty_like(yy)
        xx0.fill(x0)
        xx1.fill(x1)

        xxx = np.concatenate([xx[:-1], xx1[:-1], xx[-1:0:-1], xx0])
        yyy = np.concatenate([yy0[:-1], yy[:-1], yy1[:-1], yy[::-1]])
        t = tr.transform(np.array([xxx, yyy]).transpose())

        return t












class FloatingAxesBase(object):


    def __init__(self, *kl, **kwargs):
        grid_helper = kwargs.get("grid_helper", None)
        if grid_helper is None:
            raise ValueError("FloatingAxes requires grid_helper argument")
        if not hasattr(grid_helper, "get_boundary"):
            raise ValueError("grid_helper must implement get_boundary method")

        self._axes_class_floating.__init__(self, *kl, **kwargs)

        self.set_aspect(1.)
        self.adjust_axes_lim()


    def _gen_axes_patch(self):
        """
        Returns the patch used to draw the background of the axes.  It
        is also used as the clipping path for any data elements on the
        axes.

        In the standard axes, this is a rectangle, but in other
        projections it may not be.

        .. note::
            Intended to be overridden by new projection types.
        """
        import matplotlib.patches as mpatches
        grid_helper = self.get_grid_helper()
        t = grid_helper.get_boundary()
        return mpatches.Polygon(t)

    def cla(self):
        self._axes_class_floating.cla(self)
        #HostAxes.cla(self)
        self.patch.set_transform(self.transData)


        patch = self._axes_class_floating._gen_axes_patch(self)
        patch.set_figure(self.figure)
        patch.set_visible(False)
        patch.set_transform(self.transAxes)

        self.patch.set_clip_path(patch)
        self.gridlines.set_clip_path(patch)

        self._original_patch = patch


    def adjust_axes_lim(self):

        #t = self.get_boundary()
        grid_helper = self.get_grid_helper()
        t = grid_helper.get_boundary()
        x, y = t[:,0], t[:,1]

        xmin, xmax = min(x), max(x)
        ymin, ymax = min(y), max(y)

        dx = (xmax-xmin)/100.
        dy = (ymax-ymin)/100.

        self.set_xlim(xmin-dx, xmax+dx)
        self.set_ylim(ymin-dy, ymax+dy)



_floatingaxes_classes = {}

def floatingaxes_class_factory(axes_class):

    new_class = _floatingaxes_classes.get(axes_class)
    if new_class is None:
        new_class = type(str("Floating %s" % (axes_class.__name__)),
                         (FloatingAxesBase, axes_class),
                         {'_axes_class_floating': axes_class})
        _floatingaxes_classes[axes_class] = new_class

    return new_class

from .axislines import Axes
from mpl_toolkits.axes_grid1.parasite_axes import host_axes_class_factory

FloatingAxes = floatingaxes_class_factory(host_axes_class_factory(Axes))


import matplotlib.axes as maxes
FloatingSubplot = maxes.subplot_class_factory(FloatingAxes)
