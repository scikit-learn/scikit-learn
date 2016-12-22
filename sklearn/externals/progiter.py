# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import sys
import time
import datetime
import collections
from math import log10, floor
from sklearn.externals import six

# VT100 ANSI definitions
# https://en.wikipedia.org/wiki/ANSI_escape_code#CSI_codes
CLEARLINE_EL0 = '\33[0K'  # clear line to right
CLEARLINE_EL1 = '\33[1K'  # clear line to left
CLEARLINE_EL2 = '\33[2K'  # clear line
DECTCEM_HIDE = '\033[?25l'  # hide cursor
DECTCEM_SHOW = '\033[?25h'  # show cursor

WIN32 = sys.platform.startswith('win32')
WITH_ANSI = not WIN32

if WIN32:
    # Use time.clock in win32
    default_timer = time.clock
else:
    default_timer = time.time

if WITH_ANSI:
    CLEAR_BEFORE = '\r'
    AT_END = '\n'
    CLEAR_AFTER = ''
else:
    CLEAR_BEFORE = '\r' + CLEARLINE_EL2 + DECTCEM_HIDE
    CLEAR_AFTER = CLEARLINE_EL0
    AT_END = DECTCEM_SHOW + '\n'


def _infer_length(iterable):
    # use PEP 424 length hint if available
    # adapted from click implementation
    try:
        return len(iterable)
    except (AttributeError, TypeError):
        try:
            get_hint = type(iterable).__length_hint__
        except AttributeError:
            return None
        try:
            hint = get_hint(iterable)
        except TypeError:
            return None
        if (hint is NotImplemented or
             not isinstance(hint, int) or
             hint < 0):
            return None
        return hint


class ProgIter(object):
    """
    Attributes
    ----------

    iterable : sequence
        A python iterable

    label : int
        Maximum length of the process
            (estimated from iterable if not specified)

    label : str
        Message to print

    freq : int
        How many iterations to wait between messages.

    adjust : bool
        if True freq is adjusted based on time_thresh

    eta_window : int
        number of previous measurements to use in eta calculation

    clearline : bool
        if true messages are printed on the same line

    adjust : bool
        if True `freq` is adjusted based on time_thresh

    time_thresh : float
        desired amount of time to wait between messages if adjust is True
        otherwise does nothing

    stream : file
        defaults to sys.stdout

    enabled : bool
         if False nothing happens.

    verbose : int
        verbosity mode
        0 - no verbosity,
        1 - verbosity with clearline=True and adjust=True
        2 - verbosity without clearline=False and adjust=True
        3 - verbosity without clearline=False and adjust=False

    Examples
    ----------
    >>> from sklearn.externals.progiter import ProgIter
    >>> def is_prime(n):
    ...     return n >= 2 and not any(n % i == 0 for i in range(2, n))
    >>> for n in ProgIter(range(10000), verbose=2):
    >>>     # do some work
    >>>     is_prime(n)
    10000/10000... rate=13294.94 Hz, eta=0:00:00, total=0:00:00, wall=13:34 EST

    Notes
    ----------
    Either use ProgIter in a with statement or call prog.end() at the end of
    the computation if there is a possibility that the entire iterable may not
    be exhausted.
    """
    def __init__(self, iterable=None, label=None, length=None, freq=1,
                 eta_window=64, clearline=True, adjust=True, time_thresh=2.0,
                 enabled=True, verbose=None, stream=None):
        if label is None:
            label = ''
        if verbose is not None:
            if verbose <= 0:
                enabled = False
            elif verbose == 1:
                enabled, clearline, adjust = 1, 1, 1
            elif verbose == 2:
                enabled, clearline, adjust = 1, 0, 1
            elif verbose >= 3:
                enabled, clearline, adjust = 1, 0, 0
        if stream is None:
            stream = sys.stdout

        self.stream = stream
        self.iterable = iterable
        self.label = label
        self.length = length
        self.freq = freq
        self.enabled = enabled
        self.adjust = adjust
        self.eta_window = eta_window
        self.time_thresh = 1.0
        self.clearline = clearline
        self.extra = ''
        self.started = False
        self.finished = False

    def __call__(self, iterable):
        self.iterable = iterable
        return iter(self)

    def __enter__(self):
        return self

    def __exit__(self, type_, value, trace):
        if trace is not None:
            return False
        else:
            self.end()

    def __iter__(self):
        if not self.enabled:
            return iter(self.iterable)
        else:
            return self.iter_rate()

    def set_extra(self, extra):
        """
        specify a custom info appended to the end of the next message
        TODO: come up with a better name and rename
        """
        self.extra = extra

    def iter_rate(self):
        self.begin()
        # Wrap input iterable in a generator
        for self._iter_idx, item in enumerate(self.iterable, start=1):
            yield item
            if (self._iter_idx) % self.freq == 0:
                # update progress information every so often
                self.update_measurements()
                self.update_estimates()
                self.display_message()
        self.end()

    def mark(self):
        self.update_measurements()
        self.update_estimates()
        self.display_message()

    def reset_internals(self):
        # Prepare for iteration
        if self.length is None:
            self.length = _infer_length(self.iterable)
        self._est_seconds_left = None
        self._total_seconds = 0
        self._between_time = 0
        self._iter_idx = 0
        self._last_idx = -1
        # now time is actually not right now
        # now refers the the most recent measurment
        # last refers to the measurement before that
        self._now_idx = 0
        self._now_time = 0
        self._between_count = -1
        self._max_between_time = -1.0
        self._max_between_count = -1.0
        self._iters_per_second = 0.0

    def begin(self):
        """
        Initializes information used to measure progress
        """
        # Prepare for iteration
        if self.length is None:
            self.length = _infer_length(self.iterable)

        self.reset_internals()
        self._msg_fmtstr = self.build_message_template()

        self.tryflush()
        self.display_message()

        # Time progress was initialized
        self._start_time = default_timer()
        # Last time measures were udpated
        self._last_time  = self._start_time
        self._now_idx = self._iter_idx
        self._now_time = self._start_time

        # use last few times to compute a more stable average rate
        if self.eta_window is not None:
            self._measured_times = collections.deque(
                [], maxlen=self.eta_window)
            self._measured_times.append((self._iter_idx, self._start_time))

        self._cursor_at_newline = True
        self.started = True
        self.finished = False

    def end(self):
        if not self.enabled or self.finished:
            return
        # Write the final progress line if it was not written in the loop
        if self._iter_idx != self._now_idx:
            self.update_measurements()
            self.update_estimates()
            self._est_seconds_left = 0
            self.display_message()
        self.ensure_newline()
        self._cursor_at_newline = True
        self.finished = True

    def adjust_frequency(self):
        # Adjust frequency so the next print will not happen until
        # approximatly `time_thresh` seconds have passed as estimated by
        # iter_idx.
        eps = 1E-9
        self._max_between_time = max(self._max_between_time,
                                     self._between_time)
        self._max_between_time = max(self._max_between_time, eps)
        self._max_between_count = max(self._max_between_count,
                                      self._between_count)

        # If progress was uniform and all time estimates were
        # perfect this would be the new freq to achieve self.time_thresh
        new_freq = int(self.time_thresh * self._max_between_count /
                       self._max_between_time)
        new_freq = max(new_freq, 1)
        # But things are not perfect. So, don't make drastic changes
        factor = 1.5
        max_freq_change_up = max(256, int(self.freq * factor))
        max_freq_change_down = int(self.freq // factor)
        if (new_freq - self.freq) > max_freq_change_up:
            self.freq += max_freq_change_up
        elif (self.freq - new_freq) > max_freq_change_down:
            self.freq -= max_freq_change_down
        else:
            self.freq = new_freq

    def update_measurements(self):
        """
        update current measurements and estimated of time and progress
        """
        self._last_idx = self._now_idx
        self._last_time  = self._now_time

        self._now_idx = self._iter_idx
        self._now_time = default_timer()

        self._between_time = self._now_time - self._last_time
        self._between_count = self._now_idx - self._last_idx
        self._total_seconds = self._now_time - self._start_time

        # Record that measures were updated

    def update_estimates(self):
        # Estimate rate of progress
        if self.eta_window is None:
            self._iters_per_second = self._now_idx / self._total_seconds
        else:
            # Smooth out rate with a window
            self._measured_times.append((self._now_idx, self._now_time))
            prev_idx, prev_time = self._measured_times[0]
            self._iters_per_second =  ((self._now_idx - prev_idx) /
                                       (self._now_time - prev_time))

        if self.length is not None:
            # Estimate time remaining if length is given
            iters_left = self.length - self._now_idx
            est_eta = iters_left / self._iters_per_second
            self._est_seconds_left  = est_eta

        # Adjust frequency if printing too quickly
        # so progress doesnt slow down actual function
        if self.adjust and (self._between_time < self.time_thresh or
                            self._between_time > self.time_thresh * 2.0):
            self.adjust_frequency()

    def build_message_template(self):
        """ Defines the template for the progress line """
        tzname = time.tzname[0]
        if self.length <= 0:
            n_chrs = 4
        else:
            n_chrs = int(floor(log10(float(self.length))) + 1)
        msg_body = [
            (self.label),
            (' {iter_idx:' + str(n_chrs) + 'd}/'),
            ('?' if self.length <= 0 else six.text_type(self.length)),
            ('... '),
            ('rate={rate:4.2f} Hz,'),
            ('' if self.length == 0 else ' eta={eta},'),
            (' total={total},'),
            (' wall={wall} ' + tzname),
            (' {extra}'),
        ]
        if self.clearline:
            msg_body = [CLEAR_BEFORE] + msg_body + [CLEAR_AFTER]
        else:
            msg_body = msg_body + [AT_END]
        msg_fmtstr_time = ''.join(msg_body)
        return msg_fmtstr_time

    def format_message(self):
        """ formats the progress template with current values """
        if self._est_seconds_left is None:
            eta = '?'
        else:
            eta = six.text_type(datetime.timedelta(
                seconds=int(self._est_seconds_left)))
        total = six.text_type(datetime.timedelta(
            seconds=int(self._total_seconds)))
        msg = self._msg_fmtstr.format(
            iter_idx=self._now_idx,
            rate=self._iters_per_second,
            eta=eta, total=total,
            wall=time.strftime('%H:%M'),
            extra=self.extra,
        )
        return msg

    def ensure_newline(self):
        """
        use before any custom printing when using the progress iter to ensure
        your print statement starts on a new line instead of at the end of a
        progress line
        """
        if not self._cursor_at_newline:
            self.write(AT_END)
            self._cursor_at_newline = True

    def display_message(self):
        """ Writes current progress to the output stream """
        msg = self.format_message()
        self.write(msg)
        self.tryflush()
        self._cursor_at_newline = not self.clearline

    def tryflush(self):
        try:
            # flush sometimes causes issues in IPython notebooks
            self.stream.flush()
        except IOError:
            pass

    def write(self, msg):
        self.stream.write(msg)


class Timer(object):
    """
    Timer with-statment context object.
    """
    def __init__(self, msg='', verbose=True, newline=True):
        self.msg = msg
        self.verbose = verbose
        self.newline = newline
        self.tstart = -1
        self.ellapsed = -1

    def tic(self):
        if self.verbose:
            sys.stdout.flush()
            sys.stdout.write('\ntic(%r)' % self.msg)
            if self.newline:
                sys.stdout.write('\n')
            sys.stdout.flush()
        self.tstart = default_timer()

    def toc(self):
        ellapsed = (default_timer() - self.tstart)
        if self.verbose:
            sys.stdout.write('...toc(%r)=%.4fs\n' % (self.msg, ellapsed))
            sys.stdout.flush()
        return ellapsed

    start = tic
    stop = toc

    def __enter__(self):
        self.tic()
        return self

    def __exit__(self, type_, value, trace):
        self.ellapsed = self.toc()
        if trace is not None:
            return False


def test_progiter():
    from six.moves import cStringIO
    from sklearn.externals.progiter import ProgIter
    # Define a function that takes some time
    def is_prime(n):
        return n >= 2 and not any(n % i == 0 for i in range(2, n))
    N = 5000

    if False:
        stream = cStringIO()
        prog = ProgIter(range(N), clearline=False, stream=stream, freq=500,
                        adjust=False)
        stream.seek(0)
        print(stream.read())

        prog = ProgIter(range(N), clearline=False)
        for n in prog:
            was_prime = is_prime(n)
            prog.set_extra('n=%r, was_prime=%r' % (n, was_prime,))
            if (n + 1) % 128 == 0 and was_prime:
                prog.set_extra('n=%r, was_prime=%r EXTRA' % (n, was_prime,))
        stream.seek(0)
        print(stream.read())

    length = 1000
    N = 50000
    N0 = N - length
    print('N = %r' % (N,))
    print('N0 = %r' % (N0,))

    print('\n-----')
    print('Demo #0: progress can be disabled and incur essentially 0 overhead')
    print('However, the overhead of enabled progress is minimal and typically '
          'insignificant')
    print('this is verbosity mode verbose=0')
    iterable = (is_prime(n) for n in range(N0, N))
    with Timer('demo0'):
        piterable = ProgIter(iterable, length=length, label='demo0',
                             enabled=False)
        list(piterable)

    print('\n-----')
    print('Demo #1: progress is shown by default in the same line')
    print('this is verbosity mode verbose=1')
    iterable = (is_prime(n) for n in range(N0, N))
    with Timer('demo1'):
        piterable = ProgIter(iterable, length=length, label='demo1')
        list(piterable)

    # Default behavior adjusts frequency of progress reporting so
    # the performance of the loop is minimally impacted
    print('\n-----')
    print('Demo #2: clearline=False prints multiple lines.')
    print('Progress is only printed as needed')
    print('Notice the adjustment behavior of the print frequency')
    print('this is verbosity mode verbose=2')
    with Timer('demo2'):
        iterable = (is_prime(n) for n in range(N0, N))
        piterable = ProgIter(iterable, length=length, clearline=False,
                             label='demo2')
        list(piterable)
        # import utool as ut
        # print(ut.repr4(piterable.__dict__))

    print('\n-----')
    print('Demo #3: Adjustments can be turned off to give constant feedback')
    print('this is verbosity mode verbose=3')
    iterable = (is_prime(n) for n in range(N0, N))
    with Timer('demo3'):
        piterable = ProgIter(iterable, length=length, adjust=False,
                             clearline=False, freq=100, label='demo3')
        list(piterable)


def time_progiter_overhead():
    # Time the overhead of this function
    import timeit
    import textwrap
    setup = textwrap.dedent(
        '''
        from sklearn.externals.progiter import ProgIter
        import numpy as np
        import time
        from six.moves import cStringIO, range
        import utool as ut
        N = 500
        stream = cStringIO()
        rng = np.random.RandomState(42)
        ndims = 2
        vec1 = rng.rand(113, ndims)
        vec2 = rng.rand(71, ndims)

        def minimal_wraper1(iterable):
            for item in iterable:
                yield item

        def minimal_wraper2(iterable):
            for count, item in enumerate(iterable, start=1):
                yield item

        def minimal_wraper3(iterable):
            count = 0
            for item in iterable:
                yield item
                count += 1

        def minwrap4(iterable):
            for count, item in enumerate(iterable, start=1):
                yield item
                if count % 100:
                    pass

        def minwrap5(iterable):
            for count, item in enumerate(iterable, start=1):
                yield item
                if time.time() < 100:
                    pass
        '''
    )
    statements = {
        'baseline'       : '[{work} for n in range(N)]',
        'creation'       : 'ProgIter(range(N))',
        'minwrap1'       : '[{work} for n in minimal_wraper1(range(N))]',
        'minwrap2'       : '[{work} for n in minimal_wraper2(range(N))]',
        'minwrap3'       : '[{work} for n in minimal_wraper3(range(N))]',
        'minwrap4'       : '[{work} for n in minwrap4(range(N))]',
        'minwrap5'       : '[{work} for n in minwrap5(range(N))]',
        '(sk-disabled)'  : '[{work} for n in ProgIter(range(N), enabled=False, stream=stream)]',  # NOQA
        '(sk-plain)'     : '[{work} for n in ProgIter(range(N), stream=stream)]',  # NOQA
        '(sk-freq)'      : '[{work} for n in ProgIter(range(N), stream=stream, freq=100)]',  # NOQA
        '(sk-no-adjust)' : '[{work} for n in ProgIter(range(N), stream=stream, adjust=False, freq=200)]',  # NOQA
        '(sk-high-freq)' : '[{work} for n in ProgIter(range(N), stream=stream, adjust=False, freq=200)]',  # NOQA

        # '(ut-disabled)'  : '[{work} for n in ut.ProgIter(range(N), enabled=False, stream=stream)]',    # NOQA
        # '(ut-plain)'     : '[{work} for n in ut.ProgIter(range(N), stream=stream)]',  # NOQA
        # '(ut-freq)'      : '[{work} for n in ut.ProgIter(range(N), freq=100, stream=stream)]',  # NOQA
        # '(ut-no-adjust)' : '[{work} for n in ut.ProgIter(range(N), freq=200, adjust=False, stream=stream)]',  # NOQA
        # '(ut-high-freq)' : '[{work} for n in ut.ProgIter(range(N), stream=stream, adjust=False, freq=200)]',  # NOQA
    }
    # statements = {
    #     'calc_baseline': '[vec1.dot(vec2.T) for n in range(M)]',  # NOQA
    #     'calc_plain': '[vec1.dot(vec2.T) for n in ProgIter(range(M), stream=stream)]',  # NOQA
    #     'calc_plain_ut': '[vec1.dot(vec2.T) for n in ut.ProgIter(range(M), stream=stream)]',  # NOQA
    # }
    timeings = {}

    work_strs = [
        'None',
        'vec1.dot(vec2.T)',
        'n % 10 == 0',
    ]
    work = work_strs[0]
    # work = work_strs[1]

    number = 10000
    prog = ProgIter(label='timing', adjust=True)
    for key, stmt in prog(statements.items()):
        prog.set_extra(key)
        secs = timeit.timeit(stmt.format(work=work), setup, number=number)
        timeings[key] = secs / number

    import utool as ut
    print(ut.align(ut.repr4(timeings, precision=8), ':'))


if __name__ == '__main__':
    r"""
    CommandLine:
        python -m sklearn.externals.progiter
        python -m sklearn.externals.progiter --allexamples
    """
    test_progiter()
    # time_progiter_overhead()
