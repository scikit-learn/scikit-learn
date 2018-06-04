from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import matplotlib.axes as maxes
from matplotlib.artist import Artist
from matplotlib.axis import XAxis, YAxis

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
                r = SimpleChainedObjects(
                    [super(Axes.AxisDict, self).__getitem__(k1) for k1 in k])
                return r
            elif isinstance(k, slice):
                if k.start is None and k.stop is None and k.step is None:
                    r = SimpleChainedObjects(list(six.itervalues(self)))
                    return r
                else:
                    raise ValueError("Unsupported slice")
            else:
                return dict.__getitem__(self, k)

        def __call__(self, *v, **kwargs):
            return maxes.Axes.axis(self.axes, *v, **kwargs)

    def __init__(self, *kl, **kw):
        super(Axes, self).__init__(*kl, **kw)

    def _init_axis_artists(self, axes=None):
        if axes is None:
            axes = self

        self._axislines = self.AxisDict(self)

        self._axislines["bottom"] = SimpleAxisArtist(self.xaxis, 1, self.spines["bottom"])
        self._axislines["top"] = SimpleAxisArtist(self.xaxis, 2, self.spines["top"])
        self._axislines["left"] = SimpleAxisArtist(self.yaxis, 1, self.spines["left"])
        self._axislines["right"] = SimpleAxisArtist(self.yaxis, 2, self.spines["right"])


    def _get_axislines(self):
        return self._axislines

    axis = property(_get_axislines)

    def cla(self):

        super(Axes, self).cla()
        self._init_axis_artists()


class SimpleAxisArtist(Artist):
    def __init__(self, axis, axisnum, spine):
        self._axis = axis
        self._axisnum = axisnum
        self.line = spine

        if isinstance(axis, XAxis):
            self._axis_direction = ["bottom", "top"][axisnum-1]
        elif isinstance(axis, YAxis):
            self._axis_direction = ["left", "right"][axisnum-1]
        else:
            raise ValueError("axis must be instance of XAxis or YAxis : %s is provided" % (axis,))
        Artist.__init__(self)


    def _get_major_ticks(self):
        tickline = "tick%dline" % self._axisnum
        return SimpleChainedObjects([getattr(tick, tickline)
                                     for tick in self._axis.get_major_ticks()])

    def _get_major_ticklabels(self):
        label = "label%d" % self._axisnum
        return SimpleChainedObjects([getattr(tick, label)
                                     for tick in self._axis.get_major_ticks()])

    def _get_label(self):
        return self._axis.label

    major_ticks = property(_get_major_ticks)
    major_ticklabels = property(_get_major_ticklabels)
    label = property(_get_label)

    def set_visible(self, b):
        self.toggle(all=b)
        self.line.set_visible(b)
        self._axis.set_visible(True)
        Artist.set_visible(self, b)

    def set_label(self, txt):
        self._axis.set_label_text(txt)

    def toggle(self, all=None, ticks=None, ticklabels=None, label=None):

        if all:
            _ticks, _ticklabels, _label = True, True, True
        elif all is not None:
            _ticks, _ticklabels, _label = False, False, False
        else:
            _ticks, _ticklabels, _label = None, None, None

        if ticks is not None:
            _ticks = ticks
        if ticklabels is not None:
            _ticklabels = ticklabels
        if label is not None:
            _label = label

        tickOn = "tick%dOn" % self._axisnum
        labelOn = "label%dOn" % self._axisnum

        if _ticks is not None:
            tickparam = {tickOn: _ticks}
            self._axis.set_tick_params(**tickparam)
        if _ticklabels is not None:
            tickparam = {labelOn: _ticklabels}
            self._axis.set_tick_params(**tickparam)

        if _label is not None:
            pos = self._axis.get_label_position()
            if (pos == self._axis_direction) and not _label:
                self._axis.label.set_visible(False)
            elif _label:
                self._axis.label.set_visible(True)
                self._axis.set_label_position(self._axis_direction)


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = Axes(fig, [0.1, 0.1, 0.8, 0.8])
    fig.add_axes(ax)
    ax.cla()
