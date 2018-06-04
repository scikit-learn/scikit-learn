from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import numpy as np
from matplotlib.tri.triangulation import Triangulation


def triplot(ax, *args, **kwargs):
    """
    Draw a unstructured triangular grid as lines and/or markers.

    The triangulation to plot can be specified in one of two ways;
    either::

      triplot(triangulation, ...)

    where triangulation is a :class:`matplotlib.tri.Triangulation`
    object, or

    ::

      triplot(x, y, ...)
      triplot(x, y, triangles, ...)
      triplot(x, y, triangles=triangles, ...)
      triplot(x, y, mask=mask, ...)
      triplot(x, y, triangles, mask=mask, ...)

    in which case a Triangulation object will be created.  See
    :class:`~matplotlib.tri.Triangulation` for a explanation of these
    possibilities.

    The remaining args and kwargs are the same as for
    :meth:`~matplotlib.axes.Axes.plot`.

    Return a list of 2 :class:`~matplotlib.lines.Line2D` containing
    respectively:

        - the lines plotted for triangles edges
        - the markers plotted for triangles nodes
    """
    import matplotlib.axes

    tri, args, kwargs = Triangulation.get_from_args_and_kwargs(*args, **kwargs)
    x, y, edges = (tri.x, tri.y, tri.edges)

    # Decode plot format string, e.g., 'ro-'
    fmt = ""
    if len(args) > 0:
        fmt = args[0]
    linestyle, marker, color = matplotlib.axes._base._process_plot_format(fmt)

    # Insert plot format string into a copy of kwargs (kwargs values prevail).
    kw = kwargs.copy()
    for key, val in zip(('linestyle', 'marker', 'color'),
                        (linestyle, marker, color)):
        if val is not None:
            kw[key] = kwargs.get(key, val)

    # Draw lines without markers.
    # Note 1: If we drew markers here, most markers would be drawn more than
    #         once as they belong to several edges.
    # Note 2: We insert nan values in the flattened edges arrays rather than
    #         plotting directly (triang.x[edges].T, triang.y[edges].T)
    #         as it considerably speeds-up code execution.
    linestyle = kw['linestyle']
    kw_lines = kw.copy()
    kw_lines['marker'] = 'None'  # No marker to draw.
    kw_lines['zorder'] = kw.get('zorder', 1)  # Path default zorder is used.
    if (linestyle is not None) and (linestyle not in ['None', '', ' ']):
        tri_lines_x = np.insert(x[edges], 2, np.nan, axis=1)
        tri_lines_y = np.insert(y[edges], 2, np.nan, axis=1)
        tri_lines = ax.plot(tri_lines_x.ravel(), tri_lines_y.ravel(),
                            **kw_lines)
    else:
        tri_lines = ax.plot([], [], **kw_lines)

    # Draw markers separately.
    marker = kw['marker']
    kw_markers = kw.copy()
    kw_markers['linestyle'] = 'None'  # No line to draw.
    if (marker is not None) and (marker not in ['None', '', ' ']):
        tri_markers = ax.plot(x, y, **kw_markers)
    else:
        tri_markers = ax.plot([], [], **kw_markers)

    return tri_lines + tri_markers
