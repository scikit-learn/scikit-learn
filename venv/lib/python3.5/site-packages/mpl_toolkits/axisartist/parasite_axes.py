from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from mpl_toolkits.axes_grid1.parasite_axes import (
    host_axes_class_factory, parasite_axes_class_factory,
    parasite_axes_auxtrans_class_factory, subplot_class_factory)

from .axislines import Axes


ParasiteAxes = parasite_axes_class_factory(Axes)

ParasiteAxesAuxTrans = \
    parasite_axes_auxtrans_class_factory(axes_class=ParasiteAxes)

HostAxes = host_axes_class_factory(axes_class=Axes)

SubplotHost = subplot_class_factory(HostAxes)
