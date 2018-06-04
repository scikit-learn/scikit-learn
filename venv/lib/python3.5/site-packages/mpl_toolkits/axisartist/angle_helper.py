from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import numpy as np
import math

from mpl_toolkits.axisartist.grid_finder import ExtremeFinderSimple

def select_step_degree(dv):

    degree_limits_ = [1.5, 3, 7, 13, 20, 40, 70, 120, 270, 520]
    degree_steps_  = [  1, 2, 5, 10, 15, 30, 45,  90, 180, 360]
    degree_factors = [1.] * len(degree_steps_)

    minsec_limits_ = [1.5, 2.5, 3.5, 8, 11, 18, 25, 45]
    minsec_steps_  = [1,   2,   3,   5, 10, 15, 20, 30]

    minute_limits_ = np.array(minsec_limits_) / 60
    minute_factors = [60.] * len(minute_limits_)

    second_limits_ = np.array(minsec_limits_) / 3600
    second_factors = [3600.] * len(second_limits_)

    degree_limits = np.concatenate([second_limits_,
                                    minute_limits_,
                                    degree_limits_])

    degree_steps = np.concatenate([minsec_steps_,
                                   minsec_steps_,
                                   degree_steps_])

    degree_factors = np.concatenate([second_factors,
                                     minute_factors,
                                     degree_factors])

    n = degree_limits.searchsorted(dv)
    step = degree_steps[n]
    factor = degree_factors[n]

    return step, factor



def select_step_hour(dv):

    hour_limits_ = [1.5, 2.5, 3.5, 5, 7, 10, 15, 21, 36]
    hour_steps_  = [1,   2  , 3,   4, 6,  8, 12, 18, 24]
    hour_factors = [1.] * len(hour_steps_)

    minsec_limits_ = [1.5, 2.5, 3.5, 4.5, 5.5, 8, 11, 14, 18, 25, 45]
    minsec_steps_  = [1,   2,   3,   4,   5,   6, 10, 12, 15, 20, 30]

    minute_limits_ = np.array(minsec_limits_) / 60
    minute_factors = [60.] * len(minute_limits_)

    second_limits_ = np.array(minsec_limits_) / 3600
    second_factors = [3600.] * len(second_limits_)

    hour_limits = np.concatenate([second_limits_,
                                  minute_limits_,
                                  hour_limits_])

    hour_steps = np.concatenate([minsec_steps_,
                                 minsec_steps_,
                                 hour_steps_])

    hour_factors = np.concatenate([second_factors,
                                   minute_factors,
                                   hour_factors])

    n = hour_limits.searchsorted(dv)
    step = hour_steps[n]
    factor = hour_factors[n]

    return step, factor


def select_step_sub(dv):

    # subarcsec or degree
    tmp = 10.**(int(math.log10(dv))-1.)

    factor = 1./tmp

    if 1.5*tmp >= dv:
        step = 1
    elif 3.*tmp >= dv:
        step = 2
    elif 7.*tmp >= dv:
        step = 5
    else:
        step = 1
        factor = 0.1*factor

    return step, factor


def select_step(v1, v2, nv, hour=False, include_last=True,
                threshold_factor=3600.):

    if v1 > v2:
        v1, v2 = v2, v1

    dv = (v2 - v1) / nv

    if hour:
        _select_step = select_step_hour
        cycle = 24.
    else:
        _select_step = select_step_degree
        cycle = 360.

    # for degree
    if dv > 1./threshold_factor:
        step, factor = _select_step(dv)
    else:
        step, factor = select_step_sub(dv*threshold_factor)

        factor = factor * threshold_factor


    f1, f2, fstep = v1*factor, v2*factor, step/factor
    levs = np.arange(np.floor(f1/step), np.ceil(f2/step)+0.5, dtype=int) * step

    # n : number of valid levels. If there is a cycle, e.g., [0, 90, 180,
    # 270, 360], the grid line needs to be extended from 0 to 360, so
    # we need to return the whole array. However, the last level (360)
    # needs to be ignored often. In this case, so we return n=4.

    n = len(levs)


    # we need to check the range of values
    # for example, -90 to 90, 0 to 360,

    if factor == 1. and (levs[-1] >= levs[0]+cycle): # check for cycle
        nv = int(cycle / step)
        if include_last:
            levs = levs[0] + np.arange(0, nv+1, 1) * step
        else:
            levs = levs[0] + np.arange(0, nv, 1) * step

        n = len(levs)

    return np.array(levs), n, factor


def select_step24(v1, v2, nv, include_last=True, threshold_factor=3600):
    v1, v2 = v1/15., v2/15.
    levs, n, factor =  select_step(v1, v2, nv, hour=True,
                                   include_last=include_last,
                                   threshold_factor=threshold_factor)
    return levs*15., n, factor

def select_step360(v1, v2, nv, include_last=True, threshold_factor=3600):
    return select_step(v1, v2, nv, hour=False,
                       include_last=include_last,
                       threshold_factor=threshold_factor)


class LocatorBase(object):
    def __init__(self, den, include_last=True):
        self.den = den
        self._include_last = include_last

    @property
    def nbins(self):
        return self.den

    @nbins.setter
    def nbins(self, v):
        self.den = v

    def set_params(self, nbins=None):
        if nbins is not None:
            self.den = int(nbins)


class LocatorHMS(LocatorBase):
    def __call__(self, v1, v2):
        return select_step24(v1, v2, self.den, self._include_last)

class LocatorHM(LocatorBase):
    def __call__(self, v1, v2):
        return select_step24(v1, v2, self.den, self._include_last,
                             threshold_factor=60)

class LocatorH(LocatorBase):
    def __call__(self, v1, v2):
        return select_step24(v1, v2, self.den, self._include_last,
                             threshold_factor=1)


class LocatorDMS(LocatorBase):
    def __call__(self, v1, v2):
        return select_step360(v1, v2, self.den, self._include_last)

class LocatorDM(LocatorBase):
    def __call__(self, v1, v2):
        return select_step360(v1, v2, self.den, self._include_last,
                              threshold_factor=60)

class LocatorD(LocatorBase):
    def __call__(self, v1, v2):
        return select_step360(v1, v2, self.den, self._include_last,
                              threshold_factor=1)


class FormatterDMS(object):
    deg_mark = "^{\circ}"
    min_mark = "^{\prime}"
    sec_mark = "^{\prime\prime}"

    fmt_d = "$%d" + deg_mark + "$"
    fmt_ds = r"$%d.%s" + deg_mark + "$"

    # %s for sign
    fmt_d_m = r"$%s%d" + deg_mark + "\,%02d" + min_mark + "$"
    fmt_d_ms = r"$%s%d" + deg_mark + "\,%02d.%s" + min_mark + "$"

    fmt_d_m_partial = "$%s%d" + deg_mark + "\,%02d" + min_mark + "\,"
    fmt_s_partial = "%02d" + sec_mark + "$"
    fmt_ss_partial = "%02d.%s" + sec_mark + "$"

    def _get_number_fraction(self, factor):
        ## check for fractional numbers
        number_fraction = None
        # check for 60

        for threshold in [1, 60, 3600]:
            if factor <= threshold:
                break

            d = factor // threshold
            int_log_d = int(np.floor(np.log10(d)))
            if 10**int_log_d == d and d != 1:
                number_fraction = int_log_d
                factor = factor // 10**int_log_d
                return factor, number_fraction

        return factor, number_fraction


    def __call__(self, direction, factor, values):
        if len(values) == 0:
            return []
        #ss = [[-1, 1][v>0] for v in values] #not py24 compliant
        values = np.asarray(values)
        ss = np.where(values>0, 1, -1)

        sign_map = {(-1, True):"-"}
        signs = [sign_map.get((s, v!=0), "") for s, v in zip(ss, values)]

        factor, number_fraction = self._get_number_fraction(factor)

        values = np.abs(values)

        if number_fraction is not None:
            values, frac_part = divmod(values, 10**number_fraction)
            frac_fmt = "%%0%dd" % (number_fraction,)
            frac_str = [frac_fmt % (f1,) for f1 in frac_part]

        if factor == 1:
            if number_fraction is None:
                return [self.fmt_d % (s*int(v),) for (s, v) in zip(ss, values)]
            else:
                return [self.fmt_ds % (s*int(v), f1)
                        for (s, v, f1) in zip(ss, values, frac_str)]
        elif factor == 60:
            deg_part, min_part = divmod(values, 60)
            if number_fraction is None:
                return [self.fmt_d_m % (s1, d1, m1)
                        for s1, d1, m1 in zip(signs, deg_part, min_part)]
            else:
                return [self.fmt_d_ms % (s, d1, m1, f1)
                        for s, d1, m1, f1 in zip(signs, deg_part, min_part, frac_str)]

        elif factor == 3600:
            if ss[-1] == -1:
                inverse_order = True
                values = values[::-1]
                signs = signs[::-1]
            else:
                inverse_order = False

            l_hm_old = ""
            r = []

            deg_part, min_part_ = divmod(values, 3600)
            min_part, sec_part = divmod(min_part_, 60)

            if number_fraction is None:
                sec_str = [self.fmt_s_partial % (s1,) for s1 in sec_part]
            else:
                sec_str = [self.fmt_ss_partial % (s1, f1) for s1, f1 in zip(sec_part, frac_str)]

            for s, d1, m1, s1 in zip(signs, deg_part, min_part, sec_str):
                l_hm = self.fmt_d_m_partial % (s, d1, m1)
                if l_hm != l_hm_old:
                    l_hm_old = l_hm
                    l = l_hm + s1 #l_s
                else:
                    l = "$" + s + s1
                r.append(l)

            if inverse_order:
                return r[::-1]
            else:
                return r

        else: # factor > 3600.
            return [r"$%s^{\circ}$" % (str(v),) for v in ss*values]


class FormatterHMS(FormatterDMS):
    deg_mark = "^\mathrm{h}"
    min_mark = "^\mathrm{m}"
    sec_mark = "^\mathrm{s}"

    fmt_d = "$%d" + deg_mark + "$"
    fmt_ds = r"$%d.%s" + deg_mark + "$"

    # %s for sign
    fmt_d_m = r"$%s%d" + deg_mark + "\,%02d" + min_mark+"$"
    fmt_d_ms = r"$%s%d" + deg_mark + "\,%02d.%s" + min_mark+"$"

    fmt_d_m_partial = "$%s%d" + deg_mark + "\,%02d" + min_mark + "\,"
    fmt_s_partial = "%02d" + sec_mark + "$"
    fmt_ss_partial = "%02d.%s" + sec_mark + "$"

    def __call__(self, direction, factor, values): # hour
        return FormatterDMS.__call__(self, direction, factor, np.asarray(values)/15.)





class ExtremeFinderCycle(ExtremeFinderSimple):
    """
    When there is a cycle, e.g., longitude goes from 0-360.
    """
    def __init__(self,
                 nx, ny,
                 lon_cycle = 360.,
                 lat_cycle = None,
                 lon_minmax = None,
                 lat_minmax = (-90, 90)
                 ):
        #self.transfrom_xy = transform_xy
        #self.inv_transfrom_xy = inv_transform_xy
        self.nx, self.ny = nx, ny
        self.lon_cycle, self.lat_cycle = lon_cycle, lat_cycle
        self.lon_minmax = lon_minmax
        self.lat_minmax = lat_minmax


    def __call__(self, transform_xy, x1, y1, x2, y2):
        """
        get extreme values.

        x1, y1, x2, y2 in image coordinates (0-based)
        nx, ny : number of divisions in each axis
        """
        x_, y_ = np.linspace(x1, x2, self.nx), np.linspace(y1, y2, self.ny)
        x, y = np.meshgrid(x_, y_)
        lon, lat = transform_xy(np.ravel(x), np.ravel(y))

        # iron out jumps, but algorithm should be improved.
        # This is just naive way of doing and my fail for some cases.
        # Consider replacing this with numpy.unwrap
        # We are ignoring invalid warnings. They are triggered when
        # comparing arrays with NaNs using > We are already handling
        # that correctly using np.nanmin and np.nanmax
        with np.errstate(invalid='ignore'):
            if self.lon_cycle is not None:
                lon0 = np.nanmin(lon)
                lon -= 360. * ((lon - lon0) > 180.)
            if self.lat_cycle is not None:
                lat0 = np.nanmin(lat)
                lat -= 360. * ((lat - lat0) > 180.)

        lon_min, lon_max = np.nanmin(lon), np.nanmax(lon)
        lat_min, lat_max = np.nanmin(lat), np.nanmax(lat)

        lon_min, lon_max, lat_min, lat_max = \
                 self._adjust_extremes(lon_min, lon_max, lat_min, lat_max)

        return lon_min, lon_max, lat_min, lat_max


    def _adjust_extremes(self, lon_min, lon_max, lat_min, lat_max):

        lon_min, lon_max, lat_min, lat_max = \
                 self._add_pad(lon_min, lon_max, lat_min, lat_max)

        # check cycle
        if self.lon_cycle:
            lon_max = min(lon_max, lon_min + self.lon_cycle)
        if self.lat_cycle:
            lat_max = min(lat_max, lat_min + self.lat_cycle)

        if self.lon_minmax is not None:
            min0 = self.lon_minmax[0]
            lon_min = max(min0, lon_min)
            max0 = self.lon_minmax[1]
            lon_max = min(max0, lon_max)

        if self.lat_minmax is not None:
            min0 = self.lat_minmax[0]
            lat_min = max(min0, lat_min)
            max0 = self.lat_minmax[1]
            lat_max = min(max0, lat_max)

        return lon_min, lon_max, lat_min, lat_max
