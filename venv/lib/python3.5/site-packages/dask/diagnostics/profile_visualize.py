from __future__ import absolute_import, division, print_function

import random
from bisect import bisect_left
from distutils.version import LooseVersion
from itertools import cycle
from operator import itemgetter, add

from ..utils import funcname, import_required
from ..core import istask
from ..compatibility import apply


_BOKEH_MISSING_MSG = "Diagnostics plots require `bokeh` to be installed"
_TOOLZ_MISSING_MSG = "Diagnostics plots require `toolz` to be installed"


def unquote(expr):
    if istask(expr):
        if expr[0] in (tuple, list, set):
            return expr[0](map(unquote, expr[1]))
        elif (expr[0] == dict and
              isinstance(expr[1], list) and
              isinstance(expr[1][0], list)):
            return dict(map(unquote, expr[1]))
    return expr


def pprint_task(task, keys, label_size=60):
    """Return a nicely formatted string for a task.

    Parameters
    ----------
    task:
        Value within dask graph to render as text
    keys: iterable
        List of keys within dask graph
    label_size: int (optional)
        Maximum size of output label, defaults to 60

    Examples
    --------
    >>> from operator import add, mul
    >>> dsk = {'a': 1,
    ...        'b': 2,
    ...        'c': (add, 'a', 'b'),
    ...        'd': (add, (mul, 'a', 'b'), 'c'),
    ...        'e': (sum, ['a', 'b', 5]),
    ...        'f': (add,),
    ...        'g': []}

    >>> pprint_task(dsk['c'], dsk)
    'add(_, _)'
    >>> pprint_task(dsk['d'], dsk)
    'add(mul(_, _), _)'
    >>> pprint_task(dsk['e'], dsk)
    'sum([_, _, *])'
    >>> pprint_task(dsk['f'], dsk)
    'add()'
    >>> pprint_task(dsk['g'], dsk)
    '[]'
    """
    if istask(task):
        func = task[0]
        if func is apply:
            head = funcname(task[1])
            tail = ')'
            args = unquote(task[2]) if len(task) > 2 else ()
            kwargs = unquote(task[3]) if len(task) > 3 else {}
        else:
            if hasattr(func, 'funcs'):
                head = '('.join(funcname(f) for f in func.funcs)
                tail = ')' * len(func.funcs)
            else:
                head = funcname(task[0])
                tail = ')'
            args = task[1:]
            kwargs = {}
        if args or kwargs:
            label_size2 = int((label_size - len(head) - len(tail)) //
                              (len(args) + len(kwargs)))
            pprint = lambda t: pprint_task(t, keys, label_size2)
        if args:
            if label_size2 > 5:
                args = ', '.join(pprint(t) for t in args)
            else:
                args = '...'
        else:
            args = ''
        if kwargs:
            if label_size2 > 5:
                kwargs = ', ' + ', '.join('{0}={1}'.format(k, pprint(v))
                                          for k, v in sorted(kwargs.items()))
            else:
                kwargs = ', ...'
        else:
            kwargs = ''
        return '{0}({1}{2}{3}'.format(head, args, kwargs, tail)
    elif isinstance(task, list):
        if not task:
            return '[]'
        elif len(task) > 3:
            result = pprint_task(task[:3], keys, label_size)
            return result[:-1] + ', ...]'
        else:
            label_size2 = int((label_size - 2 - 2 * len(task)) // len(task))
            args = ', '.join(pprint_task(t, keys, label_size2) for t in task)
            return '[{0}]'.format(args)
    else:
        try:
            if task in keys:
                return '_'
            else:
                return '*'
        except TypeError:
            return '*'


def get_colors(palette, funcs):
    """Get a dict mapping funcs to colors from palette.

    Parameters
    ----------
    palette : string
        Name of the bokeh palette to use, must be a member of
        bokeh.palettes.all_palettes.
    funcs : iterable
        Iterable of function names
    """
    palettes = import_required('bokeh.palettes', _BOKEH_MISSING_MSG)
    tz = import_required('toolz', _TOOLZ_MISSING_MSG)

    unique_funcs = list(sorted(tz.unique(funcs)))
    n_funcs = len(unique_funcs)
    palette_lookup = palettes.all_palettes[palette]
    keys = list(sorted(palette_lookup.keys()))
    index = keys[min(bisect_left(keys, n_funcs), len(keys) - 1)]
    palette = palette_lookup[index]
    # Some bokeh palettes repeat colors, we want just the unique set
    palette = list(tz.unique(palette))
    if len(palette) > n_funcs:
        # Consistently shuffle palette - prevents just using low-range
        random.Random(42).shuffle(palette)
    color_lookup = dict(zip(unique_funcs, cycle(palette)))
    return [color_lookup[n] for n in funcs]


def visualize(profilers, file_path=None, show=True, save=True, **kwargs):
    """Visualize the results of profiling in a bokeh plot.

    If multiple profilers are passed in, the plots are stacked vertically.

    Parameters
    ----------
    profilers : profiler or list
        Profiler or list of profilers.
    file_path : string, optional
        Name of the plot output file.
    show : boolean, optional
        If True (default), the plot is opened in a browser.
    save : boolean, optional
        If True (default), the plot is saved to disk.
    **kwargs
        Other keyword arguments, passed to bokeh.figure. These will override
        all defaults set by visualize.

    Returns
    -------
    The completed bokeh plot object.
    """
    bp = import_required('bokeh.plotting', _BOKEH_MISSING_MSG)
    import bokeh

    if LooseVersion(bokeh.__version__) >= "0.12.10":
        from bokeh.io import state
        in_notebook = state.curstate().notebook
    else:
        from bokeh.io import _state
        in_notebook = _state._notebook

    if not in_notebook:
        file_path = file_path or "profile.html"
        bp.output_file(file_path)

    if not isinstance(profilers, list):
        profilers = [profilers]
    figs = [prof._plot(**kwargs) for prof in profilers]
    # Stack the plots
    if len(figs) == 1:
        p = figs[0]
    else:
        top = figs[0]
        for f in figs[1:]:
            f.x_range = top.x_range
            f.title = None
            f.min_border_top = 20
            f.plot_height -= 30
        for f in figs[:-1]:
            f.xaxis.axis_label = None
            f.min_border_bottom = 20
            f.plot_height -= 30
        for f in figs:
            f.min_border_left = 75
            f.min_border_right = 75
        p = bp.gridplot([[f] for f in figs])
    if show:
        bp.show(p)
    if file_path and save:
        bp.save(p)
    return p


def _get_figure_keywords():
    bp = import_required('bokeh.plotting', _BOKEH_MISSING_MSG)
    o = bp.Figure.properties()
    o.add('tools')
    return o


def plot_tasks(results, dsk, palette='Viridis', label_size=60, **kwargs):
    """Visualize the results of profiling in a bokeh plot.

    Parameters
    ----------
    results : sequence
        Output of Profiler.results
    dsk : dict
        The dask graph being profiled.
    palette : string, optional
        Name of the bokeh palette to use, must be a member of
        bokeh.palettes.all_palettes.
    label_size: int (optional)
        Maximum size of output labels in plot, defaults to 60
    **kwargs
        Other keyword arguments, passed to bokeh.figure. These will override
        all defaults set by visualize.

    Returns
    -------
    The completed bokeh plot object.
    """
    bp = import_required('bokeh.plotting', _BOKEH_MISSING_MSG)
    from bokeh.models import HoverTool
    tz = import_required('toolz', _TOOLZ_MISSING_MSG)

    defaults = dict(title="Profile Results",
                    tools="hover,save,reset,xwheel_zoom,xpan",
                    toolbar_location='above',
                    plot_width=800, plot_height=300)
    defaults.update((k, v) for (k, v) in kwargs.items() if k in
                    _get_figure_keywords())

    if results:
        keys, tasks, starts, ends, ids = zip(*results)

        id_group = tz.groupby(itemgetter(4), results)
        timings = dict((k, [i.end_time - i.start_time for i in v]) for (k, v) in
                       id_group.items())
        id_lk = dict((t[0], n) for (n, t) in enumerate(sorted(timings.items(),
                     key=itemgetter(1), reverse=True)))

        left = min(starts)
        right = max(ends)

        p = bp.figure(y_range=[str(i) for i in range(len(id_lk))],
                      x_range=[0, right - left], **defaults)

        data = {}
        data['width'] = width = [e - s for (s, e) in zip(starts, ends)]
        data['x'] = [w / 2 + s - left for (w, s) in zip(width, starts)]
        data['y'] = [id_lk[i] + 1 for i in ids]
        data['function'] = funcs = [pprint_task(i, dsk, label_size) for i in tasks]
        data['color'] = get_colors(palette, funcs)
        data['key'] = [str(i) for i in keys]

        source = bp.ColumnDataSource(data=data)

        p.rect(source=source, x='x', y='y', height=1, width='width',
               color='color', line_color='gray')
    else:
        p = bp.figure(y_range=[str(i) for i in range(8)], x_range=[0, 10],
                      **defaults)
    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.yaxis.axis_label = "Worker ID"
    p.xaxis.axis_label = "Time (s)"

    hover = p.select(HoverTool)
    hover.tooltips = """
    <div>
        <span style="font-size: 14px; font-weight: bold;">Key:</span>&nbsp;
        <span style="font-size: 10px; font-family: Monaco, monospace;">@key</span>
    </div>
    <div>
        <span style="font-size: 14px; font-weight: bold;">Task:</span>&nbsp;
        <span style="font-size: 10px; font-family: Monaco, monospace;">@function</span>
    </div>
    """
    hover.point_policy = 'follow_mouse'

    return p


def plot_resources(results, palette='Viridis', **kwargs):
    """Plot resource usage in a bokeh plot.

    Parameters
    ----------
    results : sequence
        Output of ResourceProfiler.results
    palette : string, optional
        Name of the bokeh palette to use, must be a member of
        bokeh.palettes.all_palettes.
    **kwargs
        Other keyword arguments, passed to bokeh.figure. These will override
        all defaults set by plot_resources.

    Returns
    -------
    The completed bokeh plot object.
    """
    bp = import_required('bokeh.plotting', _BOKEH_MISSING_MSG)
    from bokeh import palettes
    from bokeh.models import LinearAxis, Range1d

    defaults = dict(title="Profile Results",
                    tools="save,reset,xwheel_zoom,xpan",
                    toolbar_location='above',
                    plot_width=800, plot_height=300)
    defaults.update((k, v) for (k, v) in kwargs.items() if k in
                    _get_figure_keywords())
    if results:
        t, mem, cpu = zip(*results)
        left, right = min(t), max(t)
        t = [i - left for i in t]
        p = bp.figure(y_range=fix_bounds(0, max(cpu), 100),
                      x_range=fix_bounds(0, right - left, 1),
                      **defaults)
    else:
        t = mem = cpu = []
        p = bp.figure(y_range=(0, 100), x_range=(0, 1), **defaults)
    colors = palettes.all_palettes[palette][6]
    p.line(t, cpu, color=colors[0], line_width=4, legend='% CPU')
    p.yaxis.axis_label = "% CPU"
    p.extra_y_ranges = {'memory': Range1d(*fix_bounds(min(mem) if mem else 0,
                                                      max(mem) if mem else 100,
                                                      100))}
    p.line(t, mem, color=colors[2], y_range_name='memory', line_width=4,
           legend='Memory')
    p.add_layout(LinearAxis(y_range_name='memory', axis_label='Memory (MB)'),
                 'right')
    p.xaxis.axis_label = "Time (s)"
    return p


def fix_bounds(start, end, min_span):
    """Adjust end point to ensure span of at least `min_span`"""
    return start, max(end, start + min_span)


def plot_cache(results, dsk, start_time, metric_name, palette='Viridis',
               label_size=60, **kwargs):
    """Visualize the results of profiling in a bokeh plot.

    Parameters
    ----------
    results : sequence
        Output of CacheProfiler.results
    dsk : dict
        The dask graph being profiled.
    start_time : float
        Start time of the profile.
    metric_name : string
        Metric used to measure cache size
    palette : string, optional
        Name of the bokeh palette to use, must be a member of
        bokeh.palettes.all_palettes.
    label_size: int (optional)
        Maximum size of output labels in plot, defaults to 60
    **kwargs
        Other keyword arguments, passed to bokeh.figure. These will override
        all defaults set by visualize.

    Returns
    -------
    The completed bokeh plot object.
    """
    bp = import_required('bokeh.plotting', _BOKEH_MISSING_MSG)
    from bokeh.models import HoverTool
    tz = import_required('toolz', _TOOLZ_MISSING_MSG)

    defaults = dict(title="Profile Results",
                    tools="hover,save,reset,wheel_zoom,xpan",
                    toolbar_location='above',
                    plot_width=800, plot_height=300)
    defaults.update((k, v) for (k, v) in kwargs.items() if k in
                    _get_figure_keywords())

    if results:
        starts, ends = list(zip(*results))[3:]
        tics = list(sorted(tz.unique(starts + ends)))
        groups = tz.groupby(lambda d: pprint_task(d[1], dsk, label_size), results)
        data = {}
        for k, vals in groups.items():
            cnts = dict.fromkeys(tics, 0)
            for v in vals:
                cnts[v.cache_time] += v.metric
                cnts[v.free_time] -= v.metric
            data[k] = [0] + list(tz.accumulate(add, tz.pluck(1, sorted(cnts.items()))))

        tics = [0] + [i - start_time for i in tics]
        p = bp.figure(x_range=[0, max(tics)], **defaults)

        for (key, val), color in zip(data.items(), get_colors(palette, data.keys())):
            p.line('x', 'y', line_color=color, line_width=3,
                   source=bp.ColumnDataSource({'x': tics, 'y': val,
                                               'label': [key for i in val]}))

    else:
        p = bp.figure(y_range=[0, 10], x_range=[0, 10], **defaults)
    p.yaxis.axis_label = "Cache Size ({0})".format(metric_name)
    p.xaxis.axis_label = "Time (s)"

    hover = p.select(HoverTool)
    hover.tooltips = """
    <div>
        <span style="font-size: 14px; font-weight: bold;">Task:</span>&nbsp;
        <span style="font-size: 10px; font-family: Monaco, monospace;">@label</span>
    </div>
    """
    return p
