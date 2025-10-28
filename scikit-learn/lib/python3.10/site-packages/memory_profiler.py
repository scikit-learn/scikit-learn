"""Profile the memory usage of a Python program"""

# .. we'll use this to pass it to the child script ..
_CLEAN_GLOBALS = globals().copy()

__version__ = '0.61.0'

_CMD_USAGE = "python -m memory_profiler script_file.py"

from asyncio import iscoroutinefunction
from contextlib import contextmanager
from functools import partial, wraps
from types import coroutine
import builtins
import inspect
import linecache
import logging
import os
import io
import pdb
import subprocess
import sys
import time
import traceback
import warnings

if sys.platform == "win32":
    # any value except signal.CTRL_C_EVENT and signal.CTRL_BREAK_EVENT
    # can be used to kill a process unconditionally in Windows
    SIGKILL = -1
else:
    from signal import SIGKILL
import psutil


# TODO: provide alternative when multiprocessing is not available
try:
    from multiprocessing import Process, Pipe
except ImportError:
    from multiprocessing.dummy import Process, Pipe

try:
    from IPython.core.magic import Magics, line_cell_magic, magics_class
except ImportError:
    # ipython_version < '0.13'
    Magics = object
    line_cell_magic = lambda func: func
    magics_class = lambda cls: cls

_TWO_20 = float(2 ** 20)


# .. get available packages ..
try:
    import tracemalloc

    has_tracemalloc = True
except ImportError:
    has_tracemalloc = False


class MemitResult(object):
    """memit magic run details.

    Object based on IPython's TimeitResult
    """

    def __init__(self, mem_usage, baseline, repeat, timeout, interval,
                 include_children):
        self.mem_usage = mem_usage
        self.baseline = baseline
        self.repeat = repeat
        self.timeout = timeout
        self.interval = interval
        self.include_children = include_children

    def __str__(self):
        max_mem = max(self.mem_usage)
        inc = max_mem - self.baseline
        return 'peak memory: %.02f MiB, increment: %.02f MiB' % (max_mem, inc)

    def _repr_pretty_(self, p, cycle):
        msg = str(self)
        p.text(u'<MemitResult : ' + msg + u'>')


def _get_child_memory(process, meminfo_attr=None, memory_metric=0):
    """
    Returns a generator that yields memory for all child processes.
    """
    # Convert a pid to a process
    if isinstance(process, int):
        if process == -1: process = os.getpid()
        process = psutil.Process(process)

    if not meminfo_attr:
        # Use the psutil 2.0 attr if the older version isn't passed in.
        meminfo_attr = 'memory_info' if hasattr(process, 'memory_info') else 'get_memory_info'

    # Select the psutil function get the children similar to how we selected
    # the memory_info attr (a change from excepting the AttributeError).
    children_attr = 'children' if hasattr(process, 'children') else 'get_children'

    # Loop over the child processes and yield their memory
    try:
        for child in getattr(process, children_attr)(recursive=True):
            if isinstance(memory_metric, str):
                meminfo = getattr(child, meminfo_attr)()
                yield child.pid, getattr(meminfo, memory_metric) / _TWO_20
            else:
                yield child.pid, getattr(child, meminfo_attr)()[memory_metric] / _TWO_20
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        # https://github.com/fabianp/memory_profiler/issues/71
        yield (0, 0.0)


def _get_memory(pid, backend, timestamps=False, include_children=False, filename=None):
    # .. low function to get memory consumption ..
    if pid == -1:
        pid = os.getpid()

    def tracemalloc_tool():
        # .. cross-platform but but requires Python 3.4 or higher ..
        stat = next(filter(lambda item: str(item).startswith(filename),
                           tracemalloc.take_snapshot().statistics('filename')))
        mem = stat.size / _TWO_20
        if timestamps:
            return mem, time.time()
        else:
            return mem

    def ps_util_tool():
        # .. cross-platform but but requires psutil ..
        process = psutil.Process(pid)
        try:
            # avoid using get_memory_info since it does not exists
            # in psutil > 2.0 and accessing it will cause exception.
            meminfo_attr = 'memory_info' if hasattr(process, 'memory_info') \
                else 'get_memory_info'
            mem = getattr(process, meminfo_attr)()[0] / _TWO_20
            if include_children:
                mem +=  sum([mem for (pid, mem) in _get_child_memory(process, meminfo_attr)])
            if timestamps:
                return mem, time.time()
            else:
                return mem
        except psutil.AccessDenied:
            pass
            # continue and try to get this from ps

    def _ps_util_full_tool(memory_metric):

        # .. cross-platform but requires psutil > 4.0.0 ..
        process = psutil.Process(pid)
        try:
            if not hasattr(process, 'memory_full_info'):
                raise NotImplementedError("Backend `{}` requires psutil > 4.0.0".format(memory_metric))

            meminfo_attr = 'memory_full_info'
            meminfo = getattr(process, meminfo_attr)()

            if not hasattr(meminfo, memory_metric):
                raise NotImplementedError(
                    "Metric `{}` not available. For details, see:".format(memory_metric) +
                    "https://psutil.readthedocs.io/en/latest/index.html?highlight=memory_info#psutil.Process.memory_full_info")
            mem = getattr(meminfo, memory_metric) / _TWO_20

            if include_children:
                mem +=  sum([mem for (pid, mem) in _get_child_memory(process, meminfo_attr, memory_metric)])

            if timestamps:
                return mem, time.time()
            else:
                return mem
        
        except psutil.AccessDenied:
            pass
            # continue and try to get this from ps

    def posix_tool():
        # .. scary stuff ..
        if include_children:
            raise NotImplementedError((
                "The psutil module is required to monitor the "
                "memory usage of child processes."
            ))

        warnings.warn("psutil module not found. memory_profiler will be slow")
        # ..
        # .. memory usage in MiB ..
        # .. this should work on both Mac and Linux ..
        # .. subprocess.check_output appeared in 2.7, using Popen ..
        # .. for backwards compatibility ..
        out = subprocess.Popen(['ps', 'v', '-p', str(pid)],
                               stdout=subprocess.PIPE
                               ).communicate()[0].split(b'\n')
        try:
            vsz_index = out[0].split().index(b'RSS')
            mem = float(out[1].split()[vsz_index]) / 1024
            if timestamps:
                return mem, time.time()
            else:
                return mem
        except:
            if timestamps:
                return -1, time.time()
            else:
                return -1

    if backend == 'tracemalloc' and \
            (filename is None or filename == '<unknown>'):
        raise RuntimeError(
            'There is no access to source file of the profiled function'
        )

    tools = {'tracemalloc': tracemalloc_tool,
             'psutil': ps_util_tool,
             'psutil_pss': lambda: _ps_util_full_tool(memory_metric="pss"),
             'psutil_uss': lambda: _ps_util_full_tool(memory_metric="uss"),
             'posix': posix_tool}
    return tools[backend]()


class MemTimer(Process):
    """
    Fetch memory consumption from over a time interval
    """

    def __init__(self, monitor_pid, interval, pipe, backend, max_usage=False,
                 *args, **kw):
        self.monitor_pid = monitor_pid
        self.interval = interval
        self.pipe = pipe
        self.cont = True
        self.backend = backend
        self.max_usage = max_usage
        self.n_measurements = 1

        self.timestamps = kw.pop("timestamps", False)
        self.include_children = kw.pop("include_children", False)

        # get baseline memory usage
        self.mem_usage = [
            _get_memory(self.monitor_pid, self.backend, timestamps=self.timestamps,
                        include_children=self.include_children)]
        super(MemTimer, self).__init__(*args, **kw)

    def run(self):
        self.pipe.send(0)  # we're ready
        stop = False
        while True:
            cur_mem = _get_memory(
                self.monitor_pid, self.backend, timestamps=self.timestamps,
                include_children=self.include_children,)
            if not self.max_usage:
                self.mem_usage.append(cur_mem)
            else:
                self.mem_usage[0] = max(cur_mem, self.mem_usage[0])
            self.n_measurements += 1
            if stop:
                break
            stop = self.pipe.poll(self.interval)
            # do one more iteration

        self.pipe.send(self.mem_usage)
        self.pipe.send(self.n_measurements)


def memory_usage(proc=-1, interval=.1, timeout=None, timestamps=False,
                 include_children=False, multiprocess=False, max_usage=False,
                 retval=False, stream=None, backend=None, max_iterations=None):
    """
    Return the memory usage of a process or piece of code

    Parameters
    ----------
    proc : {int, string, tuple, subprocess.Popen}, optional
        The process to monitor. Can be given by an integer/string
        representing a PID, by a Popen object or by a tuple
        representing a Python function. The tuple contains three
        values (f, args, kw) and specifies to run the function
        f(*args, **kw).
        Set to -1 (default) for current process.

    interval : float, optional
        Interval at which measurements are collected.

    timeout : float, optional
        Maximum amount of time (in seconds) to wait before returning.

    max_usage : bool, optional
        Only return the maximum memory usage (default False)

    retval : bool, optional
        For profiling python functions. Save the return value of the profiled
        function. Return value of memory_usage becomes a tuple:
        (mem_usage, retval)

    timestamps : bool, optional
        if True, timestamps of memory usage measurement are collected as well.

    include_children : bool, optional
        if True, sum the memory of all forked processes as well

    multiprocess : bool, optional
        if True, track the memory usage of all forked processes.

    stream : File
        if stream is a File opened with write access, then results are written
        to this file instead of stored in memory and returned at the end of
        the subprocess. Useful for long-running processes.
        Implies timestamps=True.

    backend : str, optional
        Current supported backends: 'psutil', 'psutil_pss', 'psutil_uss', 'posix', 'tracemalloc'
        If `backend=None` the default is "psutil" which measures RSS aka "Resident Set Size". 
        For more information on "psutil_pss" (measuring PSS) and "psutil_uss" please refer to:
        https://psutil.readthedocs.io/en/latest/index.html?highlight=memory_info#psutil.Process.memory_full_info 

    max_iterations : int
        Limits the number of iterations (calls to the process being monitored). Relevant
        when the process is a python function.

    Returns
    -------
    mem_usage : list of floating-point values
        memory usage, in MiB. It's length is always < timeout / interval
        if max_usage is given, returns the two elements maximum memory and
        number of measurements effectuated
    ret : return value of the profiled function
        Only returned if retval is set to True
    """
    backend = choose_backend(backend)
    if stream is not None:
        timestamps = True

    if not max_usage:
        ret = []
    else:
        ret = -1

    if timeout is not None:
        max_iter = int(round(timeout / interval))
    elif isinstance(proc, int):
        # external process and no timeout
        max_iter = 1
    else:
        # for a Python function wait until it finishes
        max_iter = float('inf')
        if max_iterations is not None:
            max_iter = max_iterations

    if callable(proc):
        proc = (proc, (), {})
    if isinstance(proc, (list, tuple)):
        if len(proc) == 1:
            f, args, kw = (proc[0], (), {})
        elif len(proc) == 2:
            f, args, kw = (proc[0], proc[1], {})
        elif len(proc) == 3:
            f, args, kw = (proc[0], proc[1], proc[2])
        else:
            raise ValueError

        current_iter = 0
        while True:
            current_iter += 1
            child_conn, parent_conn = Pipe()  # this will store MemTimer's results
            p = MemTimer(os.getpid(), interval, child_conn, backend,
                         timestamps=timestamps,
                         max_usage=max_usage,
                         include_children=include_children)
            p.start()
            parent_conn.recv()  # wait until we start getting memory

            # When there is an exception in the "proc" - the (spawned) monitoring processes don't get killed.
            # Therefore, the whole process hangs indefinitely. Here, we are ensuring that the process gets killed!
            try:
                returned = f(*args, **kw)
                parent_conn.send(0)  # finish timing
                ret = parent_conn.recv()
                n_measurements = parent_conn.recv()
                if max_usage:
                    # Convert the one element list produced by MemTimer to a singular value
                    ret = ret[0]
                if retval:
                    ret = ret, returned
            except Exception:
                parent = psutil.Process(os.getpid())
                for child in parent.children(recursive=True):
                    os.kill(child.pid, SIGKILL)
                p.join(0)
                raise

            p.join(5 * interval)

            if (n_measurements > 4) or (current_iter == max_iter) or (interval < 1e-6):
                break
            interval /= 10.
    elif isinstance(proc, subprocess.Popen):
        # external process, launched from Python
        line_count = 0
        while True:
            if not max_usage:
                mem_usage = _get_memory(
                    proc.pid, backend, timestamps=timestamps,
                    include_children=include_children)

                if mem_usage and stream is not None:
                    stream.write("MEM {0:.6f} {1:.4f}\n".format(*mem_usage))

                    # Write children to the stream file
                    if multiprocess:
                        for idx, chldmem in _get_child_memory(proc.pid):
                            stream.write("CHLD {0} {1:.6f} {2:.4f}\n".format(idx, chldmem, time.time()))
                else:
                    # Create a nested list with the child memory
                    if multiprocess:
                        mem_usage = [mem_usage]
                        for _, chldmem in _get_child_memory(proc.pid):
                            mem_usage.append(chldmem)

                    # Append the memory usage to the return value
                    ret.append(mem_usage)
            else:
                ret = max(ret,
                          _get_memory(
                              proc.pid, backend, include_children=include_children))
            time.sleep(interval)
            line_count += 1
            # flush every 50 lines. Make 'tail -f' usable on profile file
            if line_count > 50:
                line_count = 0
                if stream is not None:
                    stream.flush()
            if timeout is not None:
                max_iter -= 1
                if max_iter == 0:
                    break
            if proc.poll() is not None:
                break
    else:
        # external process
        if max_iter == -1:
            max_iter = 1
        counter = 0
        while counter < max_iter:
            counter += 1
            if not max_usage:
                mem_usage = _get_memory(
                    proc, backend, timestamps=timestamps,
                    include_children=include_children)
                if stream is not None:
                    stream.write("MEM {0:.6f} {1:.4f}\n".format(*mem_usage))

                    # Write children to the stream file
                    if multiprocess:
                        for idx, chldmem in _get_child_memory(proc):
                            stream.write("CHLD {0} {1:.6f} {2:.4f}\n".format(idx, chldmem, time.time()))
                else:
                    # Create a nested list with the child memory
                    if multiprocess:
                        mem_usage = [mem_usage]
                        for _, chldmem in _get_child_memory(proc):
                            mem_usage.append(chldmem)

                    # Append the memory usage to the return value
                    ret.append(mem_usage)
            else:
                ret = max([ret,
                           _get_memory(proc, backend, include_children=include_children)
                           ])

            time.sleep(interval)
            # Flush every 50 lines.
            if counter % 50 == 0 and stream is not None:
                stream.flush()
    if stream:
        return None
    return ret


# ..
# .. utility functions for line-by-line ..


def _find_script(script_name):
    """ Find the script.

    If the input is not a file, then $PATH will be searched.
    """
    if os.path.isfile(script_name):
        return script_name
    path = os.getenv('PATH', os.defpath).split(os.pathsep)
    for folder in path:
        if not folder:
            continue
        fn = os.path.join(folder, script_name)
        if os.path.isfile(fn):
            return fn

    sys.stderr.write('Could not find script {0}\n'.format(script_name))
    raise SystemExit(1)


class _TimeStamperCM(object):
    """Time-stamping context manager."""

    def __init__(self, timestamps, filename, backend, timestamper=None, func=None,
                 include_children=False):
        self.timestamps = timestamps
        self.filename = filename
        self.backend = backend
        self.ts = timestamper
        self.func = func
        self.include_children = include_children

    def __enter__(self):
        if self.ts is not None:
            self.ts.current_stack_level += 1
            self.ts.stack[self.func].append(self.ts.current_stack_level)

        self.timestamps.append(
            _get_memory(os.getpid(), self.backend, timestamps=True,
                        include_children=self.include_children, filename=self.filename))

    def __exit__(self, *args):
        if self.ts is not None:
            self.ts.current_stack_level -= 1

        self.timestamps.append(
            _get_memory(os.getpid(), self.backend, timestamps=True,
                        include_children=self.include_children, filename=self.filename))


class TimeStamper:
    """ A profiler that just records start and end execution times for
    any decorated function.
    """

    def __init__(self, backend, include_children=False):
        self.functions = {}
        self.backend = backend
        self.include_children = include_children
        self.current_stack_level = -1
        self.stack = {}

    def __call__(self, func=None, precision=None):
        if func is not None:
            if not callable(func):
                raise ValueError("Value must be callable")

            self.add_function(func)
            f = self.wrap_function(func)
            f.__module__ = func.__module__
            f.__name__ = func.__name__
            f.__doc__ = func.__doc__
            f.__dict__.update(getattr(func, '__dict__', {}))
            return f
        else:
            def inner_partial(f):
                return self.__call__(f, precision=precision)

            return inner_partial

    def timestamp(self, name="<block>"):
        """Returns a context manager for timestamping a block of code."""
        # Make a fake function
        func = lambda x: x
        func.__module__ = ""
        func.__name__ = name
        self.add_function(func)
        timestamps = []
        self.functions[func].append(timestamps)
        # A new object is required each time, since there can be several
        # nested context managers.
        try:
            filename = inspect.getsourcefile(func)
        except TypeError:
            filename = '<unknown>'
        return _TimeStamperCM(
            timestamps,
            filename,
            self.backend,
            timestamper=self,
            func=func
        )

    def add_function(self, func):
        if func not in self.functions:
            self.functions[func] = []
            self.stack[func] = []

    def wrap_function(self, func):
        """ Wrap a function to timestamp it.
        """

        def f(*args, **kwds):
            # Start time
            try:
                filename = inspect.getsourcefile(func)
            except TypeError:
                filename = '<unknown>'
            timestamps = [
                _get_memory(os.getpid(), self.backend, timestamps=True,
                            include_children=self.include_children, filename=filename)]
            self.functions[func].append(timestamps)
            try:
                with self.call_on_stack(func, *args, **kwds) as result:
                    return result
            finally:
                # end time
                timestamps.append(_get_memory(os.getpid(), self.backend, timestamps=True,
                                              include_children=self.include_children,
                                              filename=filename))

        return f

    @contextmanager
    def call_on_stack(self, func, *args, **kwds):
        self.current_stack_level += 1
        self.stack[func].append(self.current_stack_level)

        yield func(*args, **kwds)

        self.current_stack_level -= 1

    def show_results(self, stream=None):
        if stream is None:
            stream = sys.stdout

        for func, timestamps in self.functions.items():
            function_name = "%s.%s" % (func.__module__, func.__name__)
            for ts, level in zip(timestamps, self.stack[func]):
                stream.write("FUNC %s %.4f %.4f %.4f %.4f %d\n" % (
                    (function_name,) + ts[0] + ts[1] + (level,)))


class CodeMap(dict):
    def __init__(self, include_children, backend):
        self.include_children = include_children
        self._toplevel = []
        self.backend = backend

    def add(self, code, toplevel_code=None):
        if code in self:
            return

        if toplevel_code is None:
            filename = code.co_filename
            if filename.endswith((".pyc", ".pyo")):
                filename = filename[:-1]
            if not os.path.exists(filename):
                print('ERROR: Could not find file ' + filename)
                if filename.startswith(("ipython-input", "<ipython-input")):
                    print(
                        "NOTE: %mprun can only be used on functions defined in"
                        " physical files, and not in the IPython environment.")
                return

            toplevel_code = code
            (sub_lines, start_line) = inspect.getsourcelines(code)
            linenos = range(start_line,
                            start_line + len(sub_lines))
            self._toplevel.append((filename, code, linenos))
            self[code] = {}
        else:
            self[code] = self[toplevel_code]

        for subcode in filter(inspect.iscode, code.co_consts):
            self.add(subcode, toplevel_code=toplevel_code)

    def trace(self, code, lineno, prev_lineno):
        memory = _get_memory(-1, self.backend, include_children=self.include_children,
                             filename=code.co_filename)
        prev_value = self[code].get(lineno, None)
        previous_memory = prev_value[1] if prev_value else 0
        previous_inc = prev_value[0] if prev_value else 0

        prev_line_value = self[code].get(prev_lineno, None) if prev_lineno else None
        prev_line_memory = prev_line_value[1] if prev_line_value else 0
        occ_count = self[code][lineno][2] + 1 if lineno in self[code] else 1
        self[code][lineno] = (
            previous_inc + (memory - prev_line_memory),
            max(memory, previous_memory),
            occ_count,
        )

    def items(self):
        """Iterate on the toplevel code blocks."""
        for (filename, code, linenos) in self._toplevel:
            measures = self[code]
            if not measures:
                continue  # skip if no measurement
            line_iterator = ((line, measures.get(line)) for line in linenos)
            yield (filename, line_iterator)


class LineProfiler(object):
    """ A profiler that records the amount of memory for each line """

    def __init__(self, **kw):
        include_children = kw.get('include_children', False)
        backend = kw.get('backend', 'psutil')
        self.code_map = CodeMap(
            include_children=include_children, backend=backend)
        self.enable_count = 0
        self.max_mem = kw.get('max_mem', None)
        self.prevlines = []
        self.backend = choose_backend(kw.get('backend', None))
        self.prev_lineno = None

    def __call__(self, func=None, precision=1):
        if func is not None:
            self.add_function(func)
            f = self.wrap_function(func)
            f.__module__ = func.__module__
            f.__name__ = func.__name__
            f.__doc__ = func.__doc__
            f.__dict__.update(getattr(func, '__dict__', {}))
            return f
        else:
            def inner_partial(f):
                return self.__call__(f, precision=precision)

            return inner_partial

    def add_function(self, func):
        """ Record line profiling information for the given Python function.
        """
        try:
            # func_code does not exist in Python3
            code = func.__code__
        except AttributeError:
            warnings.warn("Could not extract a code object for the object %r"
                          % func)
        else:
            self.code_map.add(code)

    @contextmanager
    def _count_ctxmgr(self):
        self.enable_by_count()
        try:
            yield
        finally:
            self.disable_by_count()

    def wrap_function(self, func):
        """ Wrap a function to profile it.
        """

        if iscoroutinefunction(func):
            @coroutine
            def f(*args, **kwargs):
                with self._count_ctxmgr():
                    res = yield from func(*args, **kwargs)
                    return res
        else:
            def f(*args, **kwds):
                with self._count_ctxmgr():
                    return func(*args, **kwds)

        return f

    def runctx(self, cmd, globals, locals):
        """ Profile a single executable statement in the given namespaces.
        """
        self.enable_by_count()
        try:
            exec(cmd, globals, locals)
        finally:
            self.disable_by_count()
        return self

    def enable_by_count(self):
        """ Enable the profiler if it hasn't been enabled before.
        """
        if self.enable_count == 0:
            self.enable()
        self.enable_count += 1

    def disable_by_count(self):
        """ Disable the profiler if the number of disable requests matches the
        number of enable requests.
        """
        if self.enable_count > 0:
            self.enable_count -= 1
            if self.enable_count == 0:
                self.disable()

    def trace_memory_usage(self, frame, event, arg):
        """Callback for sys.settrace"""
        if frame.f_code in self.code_map:
            if event == 'call':
                # "call" event just saves the lineno but not the memory
                self.prevlines.append(frame.f_lineno)
            elif event == 'line':
                # trace needs current line and previous line
                self.code_map.trace(frame.f_code, self.prevlines[-1], self.prev_lineno)
                # saving previous line
                self.prev_lineno = self.prevlines[-1]
                self.prevlines[-1] = frame.f_lineno
            elif event == 'return':
                lineno = self.prevlines.pop()
                self.code_map.trace(frame.f_code, lineno, self.prev_lineno)
                self.prev_lineno = lineno

        if self._original_trace_function is not None:
            self._original_trace_function(frame, event, arg)

        return self.trace_memory_usage

    def trace_max_mem(self, frame, event, arg):
        # run into PDB as soon as memory is higher than MAX_MEM
        if event in ('line', 'return') and frame.f_code in self.code_map:
            c = _get_memory(-1, self.backend, filename=frame.f_code.co_filename)
            if c >= self.max_mem:
                t = ('Current memory {0:.2f} MiB exceeded the '
                     'maximum of {1:.2f} MiB\n'.format(c, self.max_mem))
                sys.stdout.write(t)
                sys.stdout.write('Stepping into the debugger \n')
                frame.f_lineno -= 2
                p = pdb.Pdb()
                p.quitting = False
                p.stopframe = frame
                p.returnframe = None
                p.stoplineno = frame.f_lineno - 3
                p.botframe = None
                return p.trace_dispatch

        if self._original_trace_function is not None:
            (self._original_trace_function)(frame, event, arg)

        return self.trace_max_mem

    def __enter__(self):
        self.enable_by_count()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.disable_by_count()

    def enable(self):
        self._original_trace_function = sys.gettrace()
        if self.max_mem is not None:
            sys.settrace(self.trace_max_mem)
        else:
            sys.settrace(self.trace_memory_usage)

    def disable(self):
        sys.settrace(self._original_trace_function)


def show_results(prof, stream=None, precision=1):
    if stream is None:
        stream = sys.stdout
    template = '{0:>6} {1:>12} {2:>12}  {3:>10}   {4:<}'

    for (filename, lines) in prof.code_map.items():
        header = template.format('Line #', 'Mem usage', 'Increment', 'Occurrences',
                                 'Line Contents')

        stream.write(u'Filename: ' + filename + '\n\n')
        stream.write(header + u'\n')
        stream.write(u'=' * len(header) + '\n')

        all_lines = linecache.getlines(filename)

        float_format = u'{0}.{1}f'.format(precision + 4, precision)
        template_mem = u'{0:' + float_format + '} MiB'
        for (lineno, mem) in lines:
            if mem:
                inc = mem[0]
                total_mem = mem[1]
                total_mem = template_mem.format(total_mem)
                occurrences = mem[2]
                inc = template_mem.format(inc)
            else:
                total_mem = u''
                inc = u''
                occurrences = u''
            tmp = template.format(lineno, total_mem, inc, occurrences, all_lines[lineno - 1])
            stream.write(tmp)
        stream.write(u'\n\n')


def _func_exec(stmt, ns):
    # helper for magic_memit, just a function proxy for the exec
    # statement
    exec(stmt, ns)


@magics_class
class MemoryProfilerMagics(Magics):
    # A lprun-style %mprun magic for IPython.
    @line_cell_magic
    def mprun(self, parameter_s='', cell=None):
        """ Execute a statement under the line-by-line memory profiler from the
        memory_profiler module.

        Usage, in line mode:
          %mprun -f func1 -f func2 <statement>

        Usage, in cell mode:
          %%mprun -f func1 -f func2 [statement]
          code...
          code...

        In cell mode, the additional code lines are appended to the (possibly
        empty) statement in the first line. Cell mode allows you to easily
        profile multiline blocks without having to put them in a separate
        function.

        The given statement (which doesn't require quote marks) is run via the
        LineProfiler. Profiling is enabled for the functions specified by the -f
        options. The statistics will be shown side-by-side with the code through
        the pager once the statement has completed.

        Options:

        -f <function>: LineProfiler only profiles functions and methods it is told
        to profile.  This option tells the profiler about these functions. Multiple
        -f options may be used. The argument may be any expression that gives
        a Python function or method object. However, one must be careful to avoid
        spaces that may confuse the option parser. Additionally, functions defined
        in the interpreter at the In[] prompt or via %run currently cannot be
        displayed.  Write these functions out to a separate file and import them.

        One or more -f options are required to get any useful results.

        -T <filename>: dump the text-formatted statistics with the code
        side-by-side out to a text file.

        -r: return the LineProfiler object after it has completed profiling.

        -c: If present, add the memory usage of any children process to the report.
        """
        from io import StringIO
        from memory_profiler import show_results, LineProfiler

        # Local imports to avoid hard dependency.
        from distutils.version import LooseVersion
        import IPython
        ipython_version = LooseVersion(IPython.__version__)
        if ipython_version < '0.11':
            from IPython.genutils import page
            from IPython.ipstruct import Struct
            from IPython.ipapi import UsageError
        else:
            from IPython.core.page import page
            from IPython.utils.ipstruct import Struct
            from IPython.core.error import UsageError

        # Escape quote markers.
        opts_def = Struct(T=[''], f=[])
        parameter_s = parameter_s.replace('"', r'\"').replace("'", r"\'")
        opts, arg_str = self.parse_options(parameter_s, 'rf:T:c',
                                           list_all=True)
        opts.merge(opts_def)
        global_ns = self.shell.user_global_ns
        local_ns = self.shell.user_ns

        if cell is not None:
            arg_str += '\n' + cell

        # Get the requested functions.
        funcs = []
        for name in opts.f:
            try:
                funcs.append(eval(name, global_ns, local_ns))
            except Exception as e:
                raise UsageError('Could not find function %r.\n%s: %s' % (name,
                                                                          e.__class__.__name__,
                                                                          e))

        include_children = 'c' in opts
        profile = LineProfiler(include_children=include_children)
        for func in funcs:
            profile(func)

        # Add the profiler to the builtins for @profile.
        if 'profile' in builtins.__dict__:
            had_profile = True
            old_profile = builtins.__dict__['profile']
        else:
            had_profile = False
            old_profile = None
        builtins.__dict__['profile'] = profile

        try:
            profile.runctx(arg_str, global_ns, local_ns)
            message = ''
        except SystemExit:
            message = "*** SystemExit exception caught in code being profiled."
        except KeyboardInterrupt:
            message = ("*** KeyboardInterrupt exception caught in code being "
                       "profiled.")
        finally:
            if had_profile:
                builtins.__dict__['profile'] = old_profile

        # Trap text output.
        stdout_trap = StringIO()
        show_results(profile, stdout_trap)
        output = stdout_trap.getvalue()
        output = output.rstrip()

        if ipython_version < '0.11':
            page(output, screen_lines=self.shell.rc.screen_length)
        else:
            page(output)
        print(message, )

        text_file = opts.T[0]
        if text_file:
            with open(text_file, 'w') as pfile:
                pfile.write(output)
            print('\n*** Profile printout saved to text file %s. %s' % (
                text_file,
                message))

        return_value = None
        if 'r' in opts:
            return_value = profile

        return return_value

    # a timeit-style %memit magic for IPython
    @line_cell_magic
    def memit(self, line='', cell=None):
        """Measure memory usage of a Python statement

        Usage, in line mode:
          %memit [-r<R>t<T>i<I>] statement

        Usage, in cell mode:
          %%memit [-r<R>t<T>i<I>] setup_code
          code...
          code...

        This function can be used both as a line and cell magic:

        - In line mode you can measure a single-line statement (though multiple
          ones can be chained with using semicolons).

        - In cell mode, the statement in the first line is used as setup code
          (executed but not measured) and the body of the cell is measured.
          The cell body has access to any variables created in the setup code.

        Options:
        -r<R>: repeat the loop iteration <R> times and take the best result.
        Default: 1

        -t<T>: timeout after <T> seconds. Default: None

        -i<I>: Get time information at an interval of I times per second.
            Defaults to 0.1 so that there is ten measurements per second.

        -c: If present, add the memory usage of any children process to the report.

        -o: If present, return a object containing memit run details

        -q: If present, be quiet and do not output a result.

        Examples
        --------
        ::

          In [1]: %memit range(10000)
          peak memory: 21.42 MiB, increment: 0.41 MiB

          In [2]: %memit range(1000000)
          peak memory: 52.10 MiB, increment: 31.08 MiB

          In [3]: %%memit l=range(1000000)
             ...: len(l)
             ...:
          peak memory: 52.14 MiB, increment: 0.08 MiB

        """
        from memory_profiler import memory_usage, _func_exec
        opts, stmt = self.parse_options(line, 'r:t:i:coq', posix=False,
                                        strict=False)

        if cell is None:
            setup = 'pass'
        else:
            setup = stmt
            stmt = cell

        repeat = int(getattr(opts, 'r', 1))
        if repeat < 1:
            repeat == 1
        timeout = int(getattr(opts, 't', 0))
        if timeout <= 0:
            timeout = None
        interval = float(getattr(opts, 'i', 0.1))
        include_children = 'c' in opts
        return_result = 'o' in opts
        quiet = 'q' in opts

        # I've noticed we get less noisier measurements if we run
        # a garbage collection first
        import gc
        gc.collect()

        _func_exec(setup, self.shell.user_ns)

        mem_usage = []
        counter = 0
        baseline = memory_usage()[0]
        while counter < repeat:
            counter += 1
            tmp = memory_usage((_func_exec, (stmt, self.shell.user_ns)),
                               timeout=timeout, interval=interval,
                               max_usage=True, max_iterations=1,
                               include_children=include_children)
            mem_usage.append(tmp)

        result = MemitResult(mem_usage, baseline, repeat, timeout, interval,
                             include_children)

        if not quiet:
            if mem_usage:
                print(result)
            else:
                print('ERROR: could not read memory usage, try with a '
                      'lower interval or more iterations')

        if return_result:
            return result

    @classmethod
    def register_magics(cls, ip):
        from distutils.version import LooseVersion
        import IPython
        ipython_version = LooseVersion(IPython.__version__)

        if ipython_version < '0.13':
            try:
                _register_magic = ip.define_magic
            except AttributeError:  # ipython 0.10
                _register_magic = ip.expose_magic

            _register_magic('mprun', cls.mprun.__func__)
            _register_magic('memit', cls.memit.__func__)
        else:
            ip.register_magics(cls)


# commenting out due to failures with some versions of IPython
# see https://github.com/fabianp/memory_profiler/issues/106
# # Ensuring old interface of magics expose for IPython 0.10
# magic_mprun = MemoryProfilerMagics().mprun.__func__
# magic_memit = MemoryProfilerMagics().memit.__func__


def load_ipython_extension(ip):
    """This is called to load the module as an IPython extension."""

    MemoryProfilerMagics.register_magics(ip)


def profile(func=None, stream=None, precision=1, backend='psutil'):
    """
    Decorator that will run the function and print a line-by-line profile
    """
    backend = choose_backend(backend)
    if backend == 'tracemalloc' and has_tracemalloc:
        if not tracemalloc.is_tracing():
            tracemalloc.start()
    if func is not None:
        get_prof = partial(LineProfiler, backend=backend)
        show_results_bound = partial(
            show_results, stream=stream, precision=precision
        )
        if iscoroutinefunction(func):
            @wraps(wrapped=func)
            @coroutine
            def wrapper(*args, **kwargs):
                prof = get_prof()
                val = yield from prof(func)(*args, **kwargs)
                show_results_bound(prof)
                return val
        else:
            @wraps(wrapped=func)
            def wrapper(*args, **kwargs):
                prof = get_prof()
                val = prof(func)(*args, **kwargs)
                show_results_bound(prof)
                return val

        return wrapper
    else:
        def inner_wrapper(f):
            return profile(f, stream=stream, precision=precision,
                           backend=backend)

        return inner_wrapper


def choose_backend(new_backend=None):
    """
    Function that tries to setup backend, chosen by user, and if failed,
    setup one of the allowable backends
    """

    _backend = 'no_backend'
    all_backends = [
        ('psutil', True),
        ('psutil_pss', True),
        ('psutil_uss', True),
        ('posix', os.name == 'posix'),
        ('tracemalloc', has_tracemalloc),
    ]
    backends_indices = dict((b[0], i) for i, b in enumerate(all_backends))

    if new_backend is not None:
        all_backends.insert(0, all_backends.pop(backends_indices[new_backend]))

    for n_backend, is_available in all_backends:
        if is_available:
            _backend = n_backend
            break
    if _backend != new_backend and new_backend is not None:
        warnings.warn('{0} can not be used, {1} used instead'.format(
            new_backend, _backend))
    return _backend


# Insert in the built-ins to have profile
# globally defined (global variables is not enough
# for all cases, e.g. a script that imports another
# script where @profile is used)
def exec_with_profiler(filename, profiler, backend, passed_args=[]):
    from runpy import run_module
    builtins.__dict__['profile'] = profiler
    ns = dict(_CLEAN_GLOBALS,
              profile=profiler,
             # Make sure the __file__ variable is usable
             # by the script we're profiling
              __file__=filename)
    # Make sure the script's directory in on sys.path
    # credit to line_profiler
    sys.path.insert(0, os.path.dirname(script_filename))

    _backend = choose_backend(backend)
    sys.argv = [filename] + passed_args
    try:
        if _backend == 'tracemalloc' and has_tracemalloc:
            tracemalloc.start()
        with io.open(filename, encoding='utf-8') as f:
            exec(compile(f.read(), filename, 'exec'), ns, ns)
    finally:
        if has_tracemalloc and tracemalloc.is_tracing():
            tracemalloc.stop()


def run_module_with_profiler(module, profiler, backend, passed_args=[]):
    from runpy import run_module
    builtins.__dict__['profile'] = profiler
    ns = dict(_CLEAN_GLOBALS, profile=profiler)
    _backend = choose_backend(backend)
    sys.argv = [module] + passed_args
    if _backend == 'tracemalloc' and has_tracemalloc:
        tracemalloc.start()
    try:
        run_module(module, run_name="__main__", init_globals=ns)
    finally:
        if has_tracemalloc and tracemalloc.is_tracing():
            tracemalloc.stop()


class LogFile(object):
    """File-like object to log text using the `logging` module and the log
    report can be customised."""

    def __init__(self, name=None, reportIncrementFlag=False):
        """
        :param name: name of the logger module
               reportIncrementFlag: This must be set to True if only the steps
               with memory increments are to be reported

        :type self: object
              name: string
              reportIncrementFlag: bool
        """
        self.logger = logging.getLogger(name)
        self.reportIncrementFlag = reportIncrementFlag

    def write(self, msg, level=logging.INFO):
        if self.reportIncrementFlag:
            if "MiB" in msg and float(msg.split("MiB")[1].strip()) > 0:
                self.logger.log(level, msg)
            elif msg.__contains__("Filename:") or msg.__contains__(
                    "Line Contents"):
                self.logger.log(level, msg)
        else:
            self.logger.log(level, msg)

    def flush(self):
        for handler in self.logger.handlers:
            handler.flush()


if __name__ == '__main__':
    from argparse import ArgumentParser, REMAINDER

    parser = ArgumentParser(usage=_CMD_USAGE)
    parser.add_argument('--version', action='version', version=__version__)
    parser.add_argument(
        '--pdb-mmem', dest='max_mem', metavar='MAXMEM',
        type=float, action='store',
        help='step into the debugger when memory exceeds MAXMEM')
    parser.add_argument(
        '--precision', dest='precision', type=int,
        action='store', default=3,
        help='precision of memory output in number of significant digits')
    parser.add_argument('-o', dest='out_filename', type=str,
        action='store', default=None,
        help='path to a file where results will be written')
    parser.add_argument('--timestamp', dest='timestamp', default=False,
        action='store_true',
        help='''print timestamp instead of memory measurement for
        decorated functions''')
    parser.add_argument('--include-children', dest='include_children',
        default=False, action='store_true',
        help='also include memory used by child processes')
    parser.add_argument('--backend', dest='backend', type=str, action='store',
        choices=['tracemalloc', 'psutil', 'psutil_pss', 'psutil_uss', 'posix'], default='psutil',
        help='backend using for getting memory info '
             '(one of the {tracemalloc, psutil, posix, psutil_pss, psutil_uss, posix})')
    parser.add_argument("program", nargs=REMAINDER,
        help='python script or module followed by command line arguments to run')
    args = parser.parse_args()

    if len(args.program) == 0:
        print("A program to run must be provided. Use -h for help")
        sys.exit(1)

    target = args.program[0]
    script_args = args.program[1:]
    _backend = choose_backend(args.backend)
    if args.timestamp:
        prof = TimeStamper(_backend, include_children=args.include_children)
    else:
        prof = LineProfiler(max_mem=args.max_mem, backend=_backend)

    try:
        if args.program[0].endswith('.py'):
            script_filename = _find_script(args.program[0])
            exec_with_profiler(script_filename, prof, args.backend, script_args)
        else:
            run_module_with_profiler(target, prof, args.backend, script_args)
    finally:
        if args.out_filename is not None:
            out_file = open(args.out_filename, "a")
        else:
            out_file = sys.stdout

        if args.timestamp:
            prof.show_results(stream=out_file)
        else:
            show_results(prof, precision=args.precision, stream=out_file)
