"""
Stacked area plot for 1D arrays inspired by Douglas Y'barbo's stackoverflow
answer:
http://stackoverflow.com/questions/2225995/how-can-i-create-stacked-line-graph-with-matplotlib

(http://stackoverflow.com/users/66549/doug)

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves import xrange

from cycler import cycler
import numpy as np

__all__ = ['stackplot']


def stackplot(axes, x, *args, **kwargs):
    """
    Draws a stacked area plot.

    Parameters
    ----------
    x : 1d array of dimension N

    y : 2d array (dimension MxN), or sequence of 1d arrays (each dimension 1xN)

        The data is assumed to be unstacked. Each of the following
        calls is legal::

            stackplot(x, y)               # where y is MxN
            stackplot(x, y1, y2, y3, y4)  # where y1, y2, y3, y4, are all 1xNm

    baseline : ['zero' | 'sym' | 'wiggle' | 'weighted_wiggle']
        Method used to calculate the baseline:

        - ``'zero'``: Constant zero baseline, i.e. a simple stacked plot.
        - ``'sym'``:  Symmetric around zero and is sometimes called
          'ThemeRiver'.
        - ``'wiggle'``: Minimizes the sum of the squared slopes.
        - ``'weighted_wiggle'``: Does the same but weights to account for
          size of each layer. It is also called 'Streamgraph'-layout. More
          details can be found at http://leebyron.com/streamgraph/.

    labels : Length N sequence of strings
        Labels to assign to each data series.

    colors : Length N sequence of colors
        A list or tuple of colors. These will be cycled through and used to
        colour the stacked areas.

    **kwargs :
        All other keyword arguments are passed to `Axes.fill_between()`.


    Returns
    -------
    list of `.PolyCollection`
        A list of `.PolyCollection` instances, one for each element in the
        stacked area plot.
    """

    y = np.row_stack(args)

    labels = iter(kwargs.pop('labels', []))

    colors = kwargs.pop('colors', None)
    if colors is not None:
        axes.set_prop_cycle(cycler('color', colors))

    baseline = kwargs.pop('baseline', 'zero')
    # Assume data passed has not been 'stacked', so stack it here.
    # We'll need a float buffer for the upcoming calculations.
    stack = np.cumsum(y, axis=0, dtype=np.promote_types(y.dtype, np.float32))

    if baseline == 'zero':
        first_line = 0.

    elif baseline == 'sym':
        first_line = -np.sum(y, 0) * 0.5
        stack += first_line[None, :]

    elif baseline == 'wiggle':
        m = y.shape[0]
        first_line = (y * (m - 0.5 - np.arange(m)[:, None])).sum(0)
        first_line /= -m
        stack += first_line

    elif baseline == 'weighted_wiggle':
        m, n = y.shape
        total = np.sum(y, 0)
        # multiply by 1/total (or zero) to avoid infinities in the division:
        inv_total = np.zeros_like(total)
        mask = total > 0
        inv_total[mask] = 1.0 / total[mask]
        increase = np.hstack((y[:, 0:1], np.diff(y)))
        below_size = total - stack
        below_size += 0.5 * y
        move_up = below_size * inv_total
        move_up[:, 0] = 0.5
        center = (move_up - 0.5) * increase
        center = np.cumsum(center.sum(0))
        first_line = center - 0.5 * total
        stack += first_line

    else:
        errstr = "Baseline method %s not recognised. " % baseline
        errstr += "Expected 'zero', 'sym', 'wiggle' or 'weighted_wiggle'"
        raise ValueError(errstr)

    # Color between x = 0 and the first array.
    color = axes._get_lines.get_next_color()
    coll = axes.fill_between(x, first_line, stack[0, :],
                             facecolor=color, label=next(labels, None),
                             **kwargs)
    coll.sticky_edges.y[:] = [0]
    r = [coll]

    # Color between array i-1 and array i
    for i in xrange(len(y) - 1):
        color = axes._get_lines.get_next_color()
        r.append(axes.fill_between(x, stack[i, :], stack[i + 1, :],
                                   facecolor=color, label=next(labels, None),
                                   **kwargs))
    return r
