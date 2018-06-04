from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import mpl_toolkits.axes_grid1.axes_grid as axes_grid_orig
from .axes_divider import LocatableAxes

class CbarAxes(axes_grid_orig.CbarAxesBase, LocatableAxes):
    def __init__(self, *kl, **kwargs):
        orientation=kwargs.pop("orientation", None)
        if orientation is None:
            raise ValueError("orientation must be specified")
        self.orientation = orientation
        self._default_label_on = False
        self.locator = None

        super(LocatableAxes, self).__init__(*kl, **kwargs)

    def cla(self):
        super(LocatableAxes, self).cla()
        self._config_axes()


class Grid(axes_grid_orig.Grid):
    _defaultLocatableAxesClass = LocatableAxes

class ImageGrid(axes_grid_orig.ImageGrid):
    _defaultLocatableAxesClass = LocatableAxes
    _defaultCbarAxesClass = CbarAxes

AxesGrid = ImageGrid
