from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves import map

from matplotlib.gridspec import GridSpec, SubplotSpec
from matplotlib import docstring
import matplotlib.artist as martist
from matplotlib.axes._axes import Axes

import matplotlib._layoutbox as layoutbox

import warnings
from matplotlib.cbook import mplDeprecation


class SubplotBase(object):
    """
    Base class for subplots, which are :class:`Axes` instances with
    additional methods to facilitate generating and manipulating a set
    of :class:`Axes` within a figure.
    """

    def __init__(self, fig, *args, **kwargs):
        """
        *fig* is a :class:`matplotlib.figure.Figure` instance.

        *args* is the tuple (*numRows*, *numCols*, *plotNum*), where
        the array of subplots in the figure has dimensions *numRows*,
        *numCols*, and where *plotNum* is the number of the subplot
        being created.  *plotNum* starts at 1 in the upper left
        corner and increases to the right.

        If *numRows* <= *numCols* <= *plotNum* < 10, *args* can be the
        decimal integer *numRows* * 100 + *numCols* * 10 + *plotNum*.
        """

        self.figure = fig

        if len(args) == 1:
            if isinstance(args[0], SubplotSpec):
                self._subplotspec = args[0]
            else:
                try:
                    s = str(int(args[0]))
                    rows, cols, num = map(int, s)
                except ValueError:
                    raise ValueError('Single argument to subplot must be '
                        'a 3-digit integer')
                self._subplotspec = GridSpec(rows, cols,
                                             figure=self.figure)[num - 1]
                # num - 1 for converting from MATLAB to python indexing
        elif len(args) == 3:
            rows, cols, num = args
            rows = int(rows)
            cols = int(cols)
            if isinstance(num, tuple) and len(num) == 2:
                num = [int(n) for n in num]
                self._subplotspec = GridSpec(
                        rows, cols,
                        figure=self.figure)[(num[0] - 1):num[1]]
            else:
                if num < 1 or num > rows*cols:
                    raise ValueError(
                        ("num must be 1 <= num <= {maxn}, not {num}"
                        ).format(maxn=rows*cols, num=num))
                self._subplotspec = GridSpec(
                        rows, cols, figure=self.figure)[int(num) - 1]
                # num - 1 for converting from MATLAB to python indexing
        else:
            raise ValueError('Illegal argument(s) to subplot: %s' % (args,))

        self.update_params()

        # _axes_class is set in the subplot_class_factory
        self._axes_class.__init__(self, fig, self.figbox, **kwargs)
        # add a layout box to this, for both the full axis, and the poss
        # of the axis.  We need both because the axes may become smaller
        # due to parasitic axes and hence no longer fill the subplotspec.
        if self._subplotspec._layoutbox is None:
            self._layoutbox = None
            self._poslayoutbox = None
        else:
            name = self._subplotspec._layoutbox.name + '.ax'
            name = name + layoutbox.seq_id()
            self._layoutbox = layoutbox.LayoutBox(
                    parent=self._subplotspec._layoutbox,
                    name=name,
                    artist=self)
            self._poslayoutbox = layoutbox.LayoutBox(
                    parent=self._layoutbox,
                    name=self._layoutbox.name+'.pos',
                    pos=True, subplot=True, artist=self)

    def __reduce__(self):
        # get the first axes class which does not
        # inherit from a subplotbase

        def not_subplotbase(c):
            return issubclass(c, Axes) and not issubclass(c, SubplotBase)

        axes_class = [c for c in self.__class__.mro()
                      if not_subplotbase(c)][0]
        r = [_PicklableSubplotClassConstructor(),
             (axes_class,),
             self.__getstate__()]
        return tuple(r)

    def get_geometry(self):
        """get the subplot geometry, e.g., 2,2,3"""
        rows, cols, num1, num2 = self.get_subplotspec().get_geometry()
        return rows, cols, num1 + 1  # for compatibility

    # COVERAGE NOTE: Never used internally or from examples
    def change_geometry(self, numrows, numcols, num):
        """change subplot geometry, e.g., from 1,1,1 to 2,2,3"""
        self._subplotspec = GridSpec(numrows, numcols,
                                     figure=self.figure)[num - 1]
        self.update_params()
        self.set_position(self.figbox)

    def get_subplotspec(self):
        """get the SubplotSpec instance associated with the subplot"""
        return self._subplotspec

    def set_subplotspec(self, subplotspec):
        """set the SubplotSpec instance associated with the subplot"""
        self._subplotspec = subplotspec

    def update_params(self):
        """update the subplot position from fig.subplotpars"""

        self.figbox, self.rowNum, self.colNum, self.numRows, self.numCols = \
            self.get_subplotspec().get_position(self.figure,
                                                return_all=True)

    def is_first_col(self):
        return self.colNum == 0

    def is_first_row(self):
        return self.rowNum == 0

    def is_last_row(self):
        return self.rowNum == self.numRows - 1

    def is_last_col(self):
        return self.colNum == self.numCols - 1

    # COVERAGE NOTE: Never used internally.
    def label_outer(self):
        """Only show "outer" labels and tick labels.

        x-labels are only kept for subplots on the last row; y-labels only for
        subplots on the first column.
        """
        lastrow = self.is_last_row()
        firstcol = self.is_first_col()
        if not lastrow:
            for label in self.get_xticklabels(which="both"):
                label.set_visible(False)
            self.get_xaxis().get_offset_text().set_visible(False)
            self.set_xlabel("")
        if not firstcol:
            for label in self.get_yticklabels(which="both"):
                label.set_visible(False)
            self.get_yaxis().get_offset_text().set_visible(False)
            self.set_ylabel("")

    def _make_twin_axes(self, *kl, **kwargs):
        """
        Make a twinx axes of self. This is used for twinx and twiny.
        """
        from matplotlib.projections import process_projection_requirements
        if 'sharex' in kwargs and 'sharey' in kwargs:
            # The following line is added in v2.2 to avoid breaking Seaborn,
            # which currently uses this internal API.
            if kwargs["sharex"] is not self and kwargs["sharey"] is not self:
                raise ValueError("Twinned Axes may share only one axis.")
        kl = (self.get_subplotspec(),) + kl
        projection_class, kwargs, key = process_projection_requirements(
            self.figure, *kl, **kwargs)

        ax2 = subplot_class_factory(projection_class)(self.figure,
                                                      *kl, **kwargs)
        self.figure.add_subplot(ax2)
        self.set_adjustable('datalim')
        ax2.set_adjustable('datalim')

        if self._layoutbox is not None and ax2._layoutbox is not None:
            # make the layout boxes be explicitly the same
            ax2._layoutbox.constrain_same(self._layoutbox)
            ax2._poslayoutbox.constrain_same(self._poslayoutbox)

        self._twinned_axes.join(self, ax2)
        return ax2

_subplot_classes = {}


def subplot_class_factory(axes_class=None):
    # This makes a new class that inherits from SubplotBase and the
    # given axes_class (which is assumed to be a subclass of Axes).
    # This is perhaps a little bit roundabout to make a new class on
    # the fly like this, but it means that a new Subplot class does
    # not have to be created for every type of Axes.
    if axes_class is None:
        axes_class = Axes

    new_class = _subplot_classes.get(axes_class)
    if new_class is None:
        new_class = type(str("%sSubplot") % (axes_class.__name__),
                         (SubplotBase, axes_class),
                         {'_axes_class': axes_class})
        _subplot_classes[axes_class] = new_class

    return new_class

# This is provided for backward compatibility
Subplot = subplot_class_factory()


class _PicklableSubplotClassConstructor(object):
    """
    This stub class exists to return the appropriate subplot
    class when __call__-ed with an axes class. This is purely to
    allow Pickling of Axes and Subplots.
    """
    def __call__(self, axes_class):
        # create a dummy object instance
        subplot_instance = _PicklableSubplotClassConstructor()
        subplot_class = subplot_class_factory(axes_class)
        # update the class to the desired subplot class
        subplot_instance.__class__ = subplot_class
        return subplot_instance


docstring.interpd.update(Axes=martist.kwdoc(Axes))
docstring.interpd.update(Subplot=martist.kwdoc(Axes))

"""
# this is some discarded code I was using to find the minimum positive
# data point for some log scaling fixes.  I realized there was a
# cleaner way to do it, but am keeping this around as an example for
# how to get the data out of the axes.  Might want to make something
# like this a method one day, or better yet make get_verts an Artist
# method

            minx, maxx = self.get_xlim()
            if minx<=0 or maxx<=0:
                # find the min pos value in the data
                xs = []
                for line in self.lines:
                    xs.extend(line.get_xdata(orig=False))
                for patch in self.patches:
                    xs.extend([x for x,y in patch.get_verts()])
                for collection in self.collections:
                    xs.extend([x for x,y in collection.get_verts()])
                posx = [x for x in xs if x>0]
                if len(posx):

                    minx = min(posx)
                    maxx = max(posx)
                    # warning, probably breaks inverted axis
                    self.set_xlim((0.1*minx, maxx))

"""
