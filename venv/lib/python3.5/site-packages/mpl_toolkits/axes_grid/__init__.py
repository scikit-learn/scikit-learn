from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from . import axes_size as Size
from .axes_divider import Divider, SubplotDivider, LocatableAxes, \
     make_axes_locatable
from .axes_grid import Grid, ImageGrid, AxesGrid
#from axes_divider import make_axes_locatable
from matplotlib.cbook import warn_deprecated
warn_deprecated(since='2.1',
                name='mpl_toolkits.axes_grid',
                alternative='mpl_toolkits.axes_grid1 and'
                            ' mpl_toolkits.axisartist provies the same'
                            ' functionality',
                obj_type='module')
