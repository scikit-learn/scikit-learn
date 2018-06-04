from __future__ import absolute_import, division, print_function

import sys
import threading
import time
from timeit import default_timer

from ..callbacks import Callback
from ..utils import ignoring


def format_time(t):
    """Format seconds into a human readable form.

    >>> format_time(10.4)
    '10.4s'
    >>> format_time(1000.4)
    '16min 40.4s'
    """
    m, s = divmod(t, 60)
    h, m = divmod(m, 60)
    if h:
        return '{0:2.0f}hr {1:2.0f}min {2:4.1f}s'.format(h, m, s)
    elif m:
        return '{0:2.0f}min {1:4.1f}s'.format(m, s)
    else:
        return '{0:4.1f}s'.format(s)


class ProgressBar(Callback):
    """A progress bar for dask.

    Parameters
    ----------
    minimum : int, optional
        Minimum time threshold in seconds before displaying a progress bar.
        Default is 0 (always display)
    width : int, optional
        Width of the bar
    dt : float, optional
        Update resolution in seconds, default is 0.1 seconds

    Examples
    --------

    Below we create a progress bar with a minimum threshold of 1 second before
    displaying. For cheap computations nothing is shown:

    >>> with ProgressBar(minimum=1.0):      # doctest: +SKIP
    ...     out = some_fast_computation.compute()

    But for expensive computations a full progress bar is displayed:

    >>> with ProgressBar(minimum=1.0):      # doctest: +SKIP
    ...     out = some_slow_computation.compute()
    [########################################] | 100% Completed | 10.4 s

    The duration of the last computation is available as an attribute

    >>> pbar = ProgressBar()
    >>> with pbar:                          # doctest: +SKIP
    ...     out = some_computation.compute()
    [########################################] | 100% Completed | 10.4 s
    >>> pbar.last_duration                  # doctest: +SKIP
    10.4

    You can also register a progress bar so that it displays for all
    computations:

    >>> pbar = ProgressBar()                # doctest: +SKIP
    >>> pbar.register()                     # doctest: +SKIP
    >>> some_slow_computation.compute()     # doctest: +SKIP
    [########################################] | 100% Completed | 10.4 s
    """

    def __init__(self, minimum=0, width=40, dt=0.1, out=None):
        if out is None:
            out = sys.stdout
        self._minimum = minimum
        self._width = width
        self._dt = dt
        self._file = out
        self.last_duration = 0

    def _start(self, dsk):
        self._state = None
        self._start_time = default_timer()
        # Start background thread
        self._running = True
        self._timer = threading.Thread(target=self._timer_func)
        self._timer.daemon = True
        self._timer.start()

    def _pretask(self, key, dsk, state):
        self._state = state
        self._file.flush()

    def _finish(self, dsk, state, errored):
        self._running = False
        self._timer.join()
        elapsed = default_timer() - self._start_time
        self.last_duration = elapsed
        if elapsed < self._minimum:
            return
        if not errored:
            self._draw_bar(1, elapsed)
        else:
            self._update_bar(elapsed)
        self._file.write('\n')
        self._file.flush()

    def _timer_func(self):
        """Background thread for updating the progress bar"""
        while self._running:
            elapsed = default_timer() - self._start_time
            if elapsed > self._minimum:
                self._update_bar(elapsed)
            time.sleep(self._dt)

    def _update_bar(self, elapsed):
        s = self._state
        if not s:
            self._draw_bar(0, elapsed)
            return
        ndone = len(s['finished'])
        ntasks = sum(len(s[k]) for k in ['ready', 'waiting', 'running']) + ndone
        self._draw_bar(ndone / ntasks if ntasks else 0, elapsed)

    def _draw_bar(self, frac, elapsed):
        bar = '#' * int(self._width * frac)
        percent = int(100 * frac)
        elapsed = format_time(elapsed)
        msg = '\r[{0:<{1}}] | {2}% Completed | {3}'.format(bar, self._width,
                                                           percent, elapsed)
        with ignoring(ValueError):
            self._file.write(msg)
            self._file.flush()
