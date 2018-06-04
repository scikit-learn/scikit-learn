from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

from . import axes_size as Size
from .axes_divider import Divider, SubplotDivider, LocatableAxes, \
     make_axes_locatable
from .axes_grid import Grid, ImageGrid, AxesGrid
#from axes_divider import make_axes_locatable

from .parasite_axes import host_subplot, host_axes
