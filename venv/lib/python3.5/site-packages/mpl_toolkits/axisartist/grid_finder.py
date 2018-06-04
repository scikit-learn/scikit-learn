from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import numpy as np
from matplotlib.transforms import Bbox
from . import clip_path
clip_line_to_rect = clip_path.clip_line_to_rect

import matplotlib.ticker as mticker
from matplotlib.transforms import Transform

# extremes finder

class ExtremeFinderSimple(object):
    def __init__(self, nx, ny):
        self.nx, self.ny = nx, ny

    def __call__(self, transform_xy, x1, y1, x2, y2):
        """
        get extreme values.

        x1, y1, x2, y2 in image coordinates (0-based)
        nx, ny : number of division in each axis
        """
        x_, y_ = np.linspace(x1, x2, self.nx), np.linspace(y1, y2, self.ny)
        x, y = np.meshgrid(x_, y_)
        lon, lat = transform_xy(np.ravel(x), np.ravel(y))

        lon_min, lon_max = lon.min(), lon.max()
        lat_min, lat_max = lat.min(), lat.max()

        return self._add_pad(lon_min, lon_max, lat_min, lat_max)

    def _add_pad(self, lon_min, lon_max, lat_min, lat_max):
        """ a small amount of padding is added because the current
        clipping algorithms seems to fail when the gridline ends at
        the bbox boundary.
        """
        dlon = (lon_max - lon_min) / self.nx
        dlat = (lat_max - lat_min) / self.ny

        lon_min, lon_max = lon_min - dlon, lon_max + dlon
        lat_min, lat_max = lat_min - dlat, lat_max + dlat

        return lon_min, lon_max, lat_min, lat_max



class GridFinderBase(object):
    def __init__(self,
                 extreme_finder,
                 grid_locator1,
                 grid_locator2,
                 tick_formatter1=None,
                 tick_formatter2=None):
        """
        the transData of the axes to the world coordinate.
        locator1, locator2 : grid locator for 1st and 2nd axis.

        Derived must define "transform_xy, inv_transform_xy"
        (may use update_transform)
        """
        super(GridFinderBase, self).__init__()

        self.extreme_finder = extreme_finder
        self.grid_locator1 = grid_locator1
        self.grid_locator2 = grid_locator2
        self.tick_formatter1 = tick_formatter1
        self.tick_formatter2 = tick_formatter2

    def get_grid_info(self,
                      x1, y1, x2, y2):
        """
        lon_values, lat_values : list of grid values. if integer is given,
                           rough number of grids in each direction.
        """

        extremes = self.extreme_finder(self.inv_transform_xy, x1, y1, x2, y2)

        # min & max rage of lat (or lon) for each grid line will be drawn.
        # i.e., gridline of lon=0 will be drawn from lat_min to lat_max.

        lon_min, lon_max, lat_min, lat_max = extremes
        lon_levs, lon_n, lon_factor = \
                  self.grid_locator1(lon_min, lon_max)
        lat_levs, lat_n, lat_factor = \
                  self.grid_locator2(lat_min, lat_max)

        if lon_factor is None:
            lon_values = np.asarray(lon_levs[:lon_n])
        else:
            lon_values = np.asarray(lon_levs[:lon_n]/lon_factor)
        if lat_factor is None:
            lat_values = np.asarray(lat_levs[:lat_n])
        else:
            lat_values = np.asarray(lat_levs[:lat_n]/lat_factor)


        lon_lines, lat_lines = self._get_raw_grid_lines(lon_values,
                                                        lat_values,
                                                        lon_min, lon_max,
                                                        lat_min, lat_max)

        ddx = (x2-x1)*1.e-10
        ddy = (y2-y1)*1.e-10
        bb = Bbox.from_extents(x1-ddx, y1-ddy, x2+ddx, y2+ddy)

        grid_info = {}
        grid_info["extremes"] = extremes
        grid_info["lon_lines"] = lon_lines
        grid_info["lat_lines"] = lat_lines

        grid_info["lon"] = self._clip_grid_lines_and_find_ticks(lon_lines,
                                                                lon_values,
                                                                lon_levs,
                                                                bb)

        grid_info["lat"] = self._clip_grid_lines_and_find_ticks(lat_lines,
                                                                lat_values,
                                                                lat_levs,
                                                                bb)

        tck_labels = grid_info["lon"]["tick_labels"] = dict()
        for direction in ["left", "bottom", "right", "top"]:
            levs = grid_info["lon"]["tick_levels"][direction]
            tck_labels[direction] = self.tick_formatter1(direction,
                                                         lon_factor, levs)

        tck_labels = grid_info["lat"]["tick_labels"] = dict()
        for direction in ["left", "bottom", "right", "top"]:
            levs = grid_info["lat"]["tick_levels"][direction]
            tck_labels[direction] = self.tick_formatter2(direction,
                                                         lat_factor, levs)

        return grid_info


    def _get_raw_grid_lines(self,
                            lon_values, lat_values,
                            lon_min, lon_max, lat_min, lat_max):

        lons_i = np.linspace(lon_min, lon_max, 100) # for interpolation
        lats_i = np.linspace(lat_min, lat_max, 100)

        lon_lines = [self.transform_xy(np.zeros_like(lats_i) + lon, lats_i)
                     for lon in lon_values]
        lat_lines = [self.transform_xy(lons_i, np.zeros_like(lons_i) + lat)
                     for lat in lat_values]

        return lon_lines, lat_lines


    def _clip_grid_lines_and_find_ticks(self, lines, values, levs, bb):
        gi = dict()
        gi["values"] = []
        gi["levels"] = []
        gi["tick_levels"] = dict(left=[], bottom=[], right=[], top=[])
        gi["tick_locs"] = dict(left=[], bottom=[], right=[], top=[])
        gi["lines"] = []

        tck_levels = gi["tick_levels"]
        tck_locs = gi["tick_locs"]
        for (lx, ly), v, lev in zip(lines, values, levs):
            xy, tcks = clip_line_to_rect(lx, ly, bb)
            if not xy:
                continue
            gi["levels"].append(v)
            gi["lines"].append(xy)

            for tck, direction in zip(tcks,
                                      ["left", "bottom", "right", "top"]):
                for t in tck:
                    tck_levels[direction].append(lev)
                    tck_locs[direction].append(t)

        return gi


    def update_transform(self, aux_trans):
        if isinstance(aux_trans, Transform):
            def transform_xy(x, y):
                x, y = np.asarray(x), np.asarray(y)
                ll1 = np.concatenate((x[:,np.newaxis], y[:,np.newaxis]), 1)
                ll2 = aux_trans.transform(ll1)
                lon, lat = ll2[:,0], ll2[:,1]
                return lon, lat

            def inv_transform_xy(x, y):
                x, y = np.asarray(x), np.asarray(y)
                ll1 = np.concatenate((x[:,np.newaxis], y[:,np.newaxis]), 1)
                ll2 = aux_trans.inverted().transform(ll1)
                lon, lat = ll2[:,0], ll2[:,1]
                return lon, lat

        else:
            transform_xy, inv_transform_xy = aux_trans

        self.transform_xy = transform_xy
        self.inv_transform_xy = inv_transform_xy


    def update(self, **kw):
        for k in kw:
            if k in ["extreme_finder",
                     "grid_locator1",
                     "grid_locator2",
                     "tick_formatter1",
                     "tick_formatter2"]:
                setattr(self, k, kw[k])
            else:
                raise ValueError("unknown update property '%s'" % k)


class GridFinder(GridFinderBase):

    def __init__(self,
                 transform,
                 extreme_finder=None,
                 grid_locator1=None,
                 grid_locator2=None,
                 tick_formatter1=None,
                 tick_formatter2=None):
        """
        transform : transform from the image coordinate (which will be
        the transData of the axes to the world coordinate.

        or transform = (transform_xy, inv_transform_xy)

        locator1, locator2 : grid locator for 1st and 2nd axis.
        """
        if extreme_finder is None:
            extreme_finder = ExtremeFinderSimple(20, 20)
        if grid_locator1 is None:
            grid_locator1 = MaxNLocator()
        if grid_locator2 is None:
            grid_locator2 = MaxNLocator()
        if tick_formatter1 is None:
            tick_formatter1 = FormatterPrettyPrint()
        if tick_formatter2 is None:
            tick_formatter2 = FormatterPrettyPrint()
        super(GridFinder, self).__init__(
            extreme_finder,
            grid_locator1,
            grid_locator2,
            tick_formatter1,
            tick_formatter2)
        self.update_transform(transform)


class MaxNLocator(mticker.MaxNLocator):
    def __init__(self, nbins=10, steps=None,
                 trim=True,
                 integer=False,
                 symmetric=False,
                 prune=None):
        # trim argument has no effect. It has been left for API compatibility
        mticker.MaxNLocator.__init__(self, nbins, steps=steps,
                                     integer=integer,
                                     symmetric=symmetric, prune=prune)
        self.create_dummy_axis()
        self._factor = None

    def __call__(self, v1, v2):
        if self._factor is not None:
            self.set_bounds(v1*self._factor, v2*self._factor)
            locs = mticker.MaxNLocator.__call__(self)
            return np.array(locs), len(locs), self._factor
        else:
            self.set_bounds(v1, v2)
            locs = mticker.MaxNLocator.__call__(self)
            return np.array(locs), len(locs), None

    def set_factor(self, f):
        self._factor = f


class FixedLocator(object):
    def __init__(self, locs):
        self._locs = locs
        self._factor = None


    def __call__(self, v1, v2):
        if self._factor is None:
            v1, v2 = sorted([v1, v2])
        else:
            v1, v2 = sorted([v1*self._factor, v2*self._factor])
        locs = np.array([l for l in self._locs if ((v1 <= l) and (l <= v2))])
        return locs, len(locs), self._factor

    def set_factor(self, f):
        self._factor = f



# Tick Formatter

class FormatterPrettyPrint(object):
    def __init__(self, useMathText=True):
        self._fmt = mticker.ScalarFormatter(
            useMathText=useMathText, useOffset=False)
        self._fmt.create_dummy_axis()
        self._ignore_factor = True

    def __call__(self, direction, factor, values):
        if not self._ignore_factor:
            if factor is None:
                factor = 1.
            values = [v/factor for v in values]
        #values = [v for v in values]
        self._fmt.set_locs(values)
        return [self._fmt(v) for v in values]


class DictFormatter(object):
    def __init__(self, format_dict, formatter=None):
        """
        format_dict : dictionary for format strings to be used.
        formatter : fall-back formatter
        """
        super(DictFormatter, self).__init__()
        self._format_dict = format_dict
        self._fallback_formatter = formatter

    def __call__(self, direction, factor, values):
        """
        factor is ignored if value is found in the dictionary
        """

        if self._fallback_formatter:
            fallback_strings = self._fallback_formatter(
                direction, factor, values)
        else:
            fallback_strings = [""]*len(values)

        r = [self._format_dict.get(k, v) for k, v in zip(values,
                                                         fallback_strings)]
        return r
