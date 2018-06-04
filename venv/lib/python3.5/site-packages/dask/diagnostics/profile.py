from __future__ import absolute_import, division, print_function

from collections import namedtuple
from itertools import starmap
from timeit import default_timer
from time import sleep
from multiprocessing import Process, Pipe, current_process

from ..callbacks import Callback
from ..utils import import_required


# Stores execution data for each task
TaskData = namedtuple('TaskData', ('key', 'task', 'start_time',
                                   'end_time', 'worker_id'))


class Profiler(Callback):
    """A profiler for dask execution at the task level.

    Records the following information for each task:
        1. Key
        2. Task
        3. Start time in seconds since the epoch
        4. Finish time in seconds since the epoch
        5. Worker id

    Examples
    --------

    >>> from operator import add, mul
    >>> from dask.threaded import get
    >>> dsk = {'x': 1, 'y': (add, 'x', 10), 'z': (mul, 'y', 2)}
    >>> with Profiler() as prof:
    ...     get(dsk, 'z')
    22

    >>> prof.results  # doctest: +SKIP
    [('y', (add, 'x', 10), 1435352238.48039, 1435352238.480655, 140285575100160),
     ('z', (mul, 'y', 2), 1435352238.480657, 1435352238.480803, 140285566707456)]

    These results can be visualized in a bokeh plot using the ``visualize``
    method. Note that this requires bokeh to be installed.

    >>> prof.visualize() # doctest: +SKIP

    You can activate the profiler globally

    >>> prof.register()  # doctest: +SKIP

    If you use the profiler globally you will need to clear out old results
    manually.

    >>> prof.clear()

    """
    def __init__(self):
        self._results = {}
        self.results = []
        self._dsk = {}

    def __enter__(self):
        self.clear()
        return super(Profiler, self).__enter__()

    def _start(self, dsk):
        self._dsk.update(dsk)

    def _pretask(self, key, dsk, state):
        start = default_timer()
        self._results[key] = (key, dsk[key], start)

    def _posttask(self, key, value, dsk, state, id):
        end = default_timer()
        self._results[key] += (end, id)

    def _finish(self, dsk, state, failed):
        results = dict((k, v) for k, v in self._results.items() if len(v) == 5)
        self.results += list(starmap(TaskData, results.values()))
        self._results.clear()

    def _plot(self, **kwargs):
        from .profile_visualize import plot_tasks
        return plot_tasks(self.results, self._dsk, **kwargs)

    def visualize(self, **kwargs):
        """Visualize the profiling run in a bokeh plot.

        See also
        --------
        dask.diagnostics.profile_visualize.visualize
        """
        from .profile_visualize import visualize
        return visualize(self, **kwargs)

    def clear(self):
        """Clear out old results from profiler"""
        self._results.clear()
        del self.results[:]
        self._dsk = {}


ResourceData = namedtuple('ResourceData', ('time', 'mem', 'cpu'))


class ResourceProfiler(Callback):
    """A profiler for resource use.

    Records the following each timestep
        1. Time in seconds since the epoch
        2. Memory usage in MB
        3. % CPU usage

    Examples
    --------

    >>> from operator import add, mul
    >>> from dask.threaded import get
    >>> dsk = {'x': 1, 'y': (add, 'x', 10), 'z': (mul, 'y', 2)}
    >>> with ResourceProfiler() as prof:  # doctest: +SKIP
    ...     get(dsk, 'z')
    22

    These results can be visualized in a bokeh plot using the ``visualize``
    method. Note that this requires bokeh to be installed.

    >>> prof.visualize() # doctest: +SKIP

    You can activate the profiler globally

    >>> prof.register()  # doctest: +SKIP

    If you use the profiler globally you will need to clear out old results
    manually.

    >>> prof.clear()  # doctest: +SKIP

    Note that when used as a context manager data will be collected throughout
    the duration of the enclosed block. In contrast, when registered globally
    data will only be collected while a dask scheduler is active.
    """
    def __init__(self, dt=1):
        self._dt = dt
        self._entered = False
        self._tracker = None
        self.results = []

    def _is_running(self):
        return self._tracker is not None and self._tracker.is_alive()

    def _start_collect(self):
        if not self._is_running():
            self._tracker = _Tracker(self._dt)
            self._tracker.start()
        self._tracker.parent_conn.send('collect')

    def _stop_collect(self):
        if self._is_running():
            self._tracker.parent_conn.send('send_data')
            self.results.extend(starmap(ResourceData, self._tracker.parent_conn.recv()))

    def __enter__(self):
        self._entered = True
        self.clear()
        self._start_collect()
        return super(ResourceProfiler, self).__enter__()

    def __exit__(self, *args):
        self._entered = False
        self._stop_collect()
        self.close()
        super(ResourceProfiler, self).__exit__(*args)

    def _start(self, dsk):
        self._start_collect()

    def _finish(self, dsk, state, failed):
        if not self._entered:
            self._stop_collect()

    def close(self):
        """Shutdown the resource tracker process"""
        if self._is_running():
            self._tracker.shutdown()
            self._tracker = None

    __del__ = close

    def clear(self):
        self.results = []

    def _plot(self, **kwargs):
        from .profile_visualize import plot_resources
        return plot_resources(self.results, **kwargs)

    def visualize(self, **kwargs):
        """Visualize the profiling run in a bokeh plot.

        See also
        --------
        dask.diagnostics.profile_visualize.visualize
        """
        from .profile_visualize import visualize
        return visualize(self, **kwargs)


class _Tracker(Process):
    """Background process for tracking resource usage"""
    def __init__(self, dt=1):
        psutil = import_required("psutil", "Tracking resource usage requires "
                                           "`psutil` to be installed")
        Process.__init__(self)
        self.daemon = True
        self.dt = dt
        self.parent = psutil.Process(current_process().pid)
        self.parent_conn, self.child_conn = Pipe()

    def shutdown(self):
        if not self.parent_conn.closed:
            self.parent_conn.send('shutdown')
            self.parent_conn.close()
        self.join()

    def _update_pids(self, pid):
        return [self.parent] + [p for p in self.parent.children()
                                if p.pid != pid and p.status() != 'zombie']

    def run(self):
        pid = current_process()
        data = []
        while True:
            try:
                msg = self.child_conn.recv()
            except KeyboardInterrupt:
                continue
            if msg == 'shutdown':
                break
            elif msg == 'collect':
                ps = self._update_pids(pid)
                while not data or not self.child_conn.poll():
                    tic = default_timer()
                    mem = cpu = 0
                    for p in ps:
                        try:
                            mem2 = p.memory_info().rss
                            cpu2 = p.cpu_percent()
                        except Exception: # could be a few different exceptions
                            pass
                        else:
                            # Only increment if both were successful
                            mem += mem2
                            cpu += cpu2
                    data.append((tic, mem / 1e6, cpu))
                    sleep(self.dt)
            elif msg == 'send_data':
                self.child_conn.send(data)
                data = []
        self.child_conn.close()


CacheData = namedtuple('CacheData', ('key', 'task', 'metric', 'cache_time',
                                     'free_time'))


class CacheProfiler(Callback):
    """A profiler for dask execution at the scheduler cache level.

    Records the following information for each task:
        1. Key
        2. Task
        3. Size metric
        4. Cache entry time in seconds since the epoch
        5. Cache exit time in seconds since the epoch

    Examples
    --------

    >>> from operator import add, mul
    >>> from dask.threaded import get
    >>> dsk = {'x': 1, 'y': (add, 'x', 10), 'z': (mul, 'y', 2)}
    >>> with CacheProfiler() as prof:
    ...     get(dsk, 'z')
    22

    >>> prof.results    # doctest: +SKIP
    [CacheData('y', (add, 'x', 10), 1, 1435352238.48039, 1435352238.480655),
     CacheData('z', (mul, 'y', 2), 1, 1435352238.480657, 1435352238.480803)]

    The default is to count each task (``metric`` is 1 for all tasks). Other
    functions may used as a metric instead through the ``metric`` keyword. For
    example, the ``nbytes`` function found in ``cachey`` can be used to measure
    the number of bytes in the cache.

    >>> from cachey import nbytes    # doctest: +SKIP
    >>> with CacheProfiler(metric=nbytes) as prof:  # doctest: +SKIP
    ...     get(dsk, 'z')

    The profiling results can be visualized in a bokeh plot using the
    ``visualize`` method. Note that this requires bokeh to be installed.

    >>> prof.visualize() # doctest: +SKIP

    You can activate the profiler globally

    >>> prof.register()  # doctest: +SKIP

    If you use the profiler globally you will need to clear out old results
    manually.

    >>> prof.clear()

    """

    def __init__(self, metric=None, metric_name=None):
        self.clear()
        self._metric = metric if metric else lambda value: 1
        if metric_name:
            self._metric_name = metric_name
        elif metric:
            self._metric_name = metric.__name__
        else:
            self._metric_name = 'count'

    def __enter__(self):
        self.clear()
        return super(CacheProfiler, self).__enter__()

    def _start(self, dsk):
        self._dsk.update(dsk)
        if not self._start_time:
            self._start_time = default_timer()

    def _posttask(self, key, value, dsk, state, id):
        t = default_timer()
        self._cache[key] = (self._metric(value), t)
        for k in state['released'].intersection(self._cache):
            metric, start = self._cache.pop(k)
            self.results.append(CacheData(k, dsk[k], metric, start, t))

    def _finish(self, dsk, state, failed):
        t = default_timer()
        for k, (metric, start) in self._cache.items():
            self.results.append(CacheData(k, dsk[k], metric, start, t))
        self._cache.clear()

    def _plot(self, **kwargs):
        from .profile_visualize import plot_cache
        return plot_cache(self.results, self._dsk, self._start_time,
                          self._metric_name, **kwargs)

    def visualize(self, **kwargs):
        """Visualize the profiling run in a bokeh plot.

        See also
        --------
        dask.diagnostics.profile_visualize.visualize
        """
        from .profile_visualize import visualize
        return visualize(self, **kwargs)

    def clear(self):
        """Clear out old results from profiler"""
        self.results = []
        self._cache = {}
        self._dsk = {}
        self._start_time = None
