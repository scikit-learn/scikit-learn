from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

#from mpl_toolkits.axes_grid1.axes_rgb import *
from mpl_toolkits.axes_grid1.axes_rgb import make_rgb_axes, imshow_rgb, RGBAxesBase

#import mpl_toolkits.axes_grid1.axes_rgb as axes_rgb_orig
from .axislines import Axes

class RGBAxes(RGBAxesBase):
    _defaultAxesClass = Axes
