"""Profile the memory usage of a Python program"""

__version__ = '0.25'

_CMD_USAGE = "python -m memory_profiler script_file.py"

import time, sys, os, pdb
import warnings
import linecache
import inspect
import subprocess

# TODO: provide alternative when multprocessing is not available
try:
    from multiprocessing import Process, Pipe
except ImportError:
    from multiprocessing.dummy import Process, Pipe


try:
    import psutil

    def _get_memory(pid):
        process = psutil.Process(pid)
        try:
            mem = float(process.get_memory_info()[0]) / (1024 ** 2)
        except psutil.AccessDenied:
            mem = -1
        return mem


except ImportError:

    warnings.warn("psutil module not found. memory_profiler will be slow")

    if os.name == 'posix':
        def _get_memory(pid):
            # ..
            # .. memory usage in MB ..
            # .. this should work on both Mac and Linux ..
            # .. subprocess.check_output appeared in 2.7, using Popen ..
            # .. for backwards compatibility ..
            out = subprocess.Popen(['ps', 'v', '-p', str(pid)],
                  stdout=subprocess.PIPE).communicate()[0].split(b'\n')
            try:
                vsz_index = out[0].split().index(b'RSS')
                return float(out[1].split()[vsz_index]) / 1024
            except:
                return -1
    else:
        raise NotImplementedError('The psutil module is required for non-unix '
                                  'platforms')


class Timer(Process):
    """
    Fetch memory consumption from over a time interval
    """

    def __init__(self, monitor_pid, interval, pipe, *args, **kw):
        self.monitor_pid = monitor_pid
        self.interval = interval
        self.pipe = pipe
        self.cont = True
        super(Timer, self).__init__(*args, **kw)

    def run(self):
        m = _get_memory(self.monitor_pid)
        timings = [m]
        self.pipe.send(0)  # we're ready
        while not self.pipe.poll(self.interval):
            m = _get_memory(self.monitor_pid)
            timings.append(m)
        self.pipe.send(timings)


def memory_usage(proc=-1, interval=0.0, timeout=None):

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

    Returns
    -------
    mem_usage : list of floating-poing values
        memory usage, in MB. It's length is always < timeout / interval
    """
    ret = []

    if timeout is not None:
        max_iter = int(timeout / interval)
    elif isinstance(proc, int):
        # external process and no timeout
        max_iter = 1
    else:
        # for a Python function wait until it finishes
        max_iter = float('inf')

    if hasattr(proc, '__call__'):
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

        aspec = inspect.getargspec(f)
        n_args = len(aspec.args)
        if aspec.defaults is not None:
            n_args -= len(aspec.defaults)
        if n_args != len(args):
            raise ValueError(
            'Function expects %s value(s) but %s where given'
            % (n_args, len(args)))

        child_conn, parent_conn = Pipe()  # this will store Timer's results
        p = Timer(os.getpid(), interval, child_conn)
        p.start()
        parent_conn.recv()  # wait until we start getting memory
        f(*args, **kw)
        parent_conn.send(0)  # finish timing
        ret = parent_conn.recv()
        p.join(5 * interval)
    elif isinstance(proc, subprocess.Popen):
        # external process, launched from Python
        while True:
            ret.append(_get_memory(proc.pid))
            time.sleep(interval)
            if timeout is not None:
                max_iter -= 1
                if max_iter == 0:
                    break
            if proc.poll() is not None:
                break
    else:
        # external process
        if proc == -1:
            proc = os.getpid()
        if max_iter == -1:
            max_iter = 1
        counter = 0
        while counter < max_iter:
            counter += 1
            ret.append(_get_memory(proc))
            time.sleep(interval)
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
        if folder == '':
            continue
        fn = os.path.join(folder, script_name)
        if os.path.isfile(fn):
            return fn

    sys.stderr.write('Could not find script {0}\n'.format(script_name))
    raise SystemExit(1)


class LineProfiler:
    """ A profiler that records the amount of memory for each line """

    def __init__(self, **kw):
        self.functions = list()
        self.code_map = {}
        self.enable_count = 0
        self.max_mem = kw.get('max_mem', None)

    def __call__(self, func):
        self.add_function(func)
        f = self.wrap_function(func)
        f.__module__ = func.__module__
        f.__name__ = func.__name__
        f.__doc__ = func.__doc__
        f.__dict__.update(getattr(func, '__dict__', {}))
        return f

    def add_function(self, func):
        """ Record line profiling information for the given Python function.
        """
        try:
            # func_code does not exist in Python3
            code = func.__code__
        except AttributeError:
            import warnings
            warnings.warn("Could not extract a code object for the object %r"
                          % (func,))
            return
        if code not in self.code_map:
            self.code_map[code] = {}
            self.functions.append(func)

    def wrap_function(self, func):
        """ Wrap a function to profile it.
        """

        def f(*args, **kwds):
            self.enable_by_count()
            try:
                result = func(*args, **kwds)
            finally:
                self.disable_by_count()
            return result
        return f

    def run(self, cmd):
        """ Profile a single executable statment in the main namespace.
        """
        import __main__
        main_dict = __main__.__dict__
        return self.runctx(cmd, main_dict, main_dict)

    def runctx(self, cmd, globals, locals):
        """ Profile a single executable statement in the given namespaces.
        """
        self.enable_by_count()
        try:
            exec(cmd, globals, locals)
        finally:
            self.disable_by_count()
        return self

    def runcall(self, func, *args, **kw):
        """ Profile a single function call.
        """
        # XXX where is this used ? can be removed ?
        self.enable_by_count()
        try:
            return func(*args, **kw)
        finally:
            self.disable_by_count()

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
        if event in ('line', 'return') and frame.f_code in self.code_map:
            lineno = frame.f_lineno
            if event == 'return':
                lineno += 1
            entry = self.code_map[frame.f_code].setdefault(lineno, [])
            entry.append(_get_memory(os.getpid()))

        return self.trace_memory_usage

    def trace_max_mem(self, frame, event, arg):
        # run into PDB as soon as memory is higher than MAX_MEM
        if event in ('line', 'return') and frame.f_code in self.code_map:
            c = _get_memory(os.getpid())
            if c >= self.max_mem:
                t = 'Current memory {0:.2f} MB exceeded the maximum '.format(c) + \
                    'of {0:.2f} MB\n'.format(self.max_mem)
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

        return self.trace_max_mem

    def __enter__(self):
        self.enable_by_count()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.disable_by_count()

    def enable(self):
        if self.max_mem is not None:
            sys.settrace(self.trace_max_mem)
        else:
            sys.settrace(self.trace_memory_usage)

    def disable(self):
        self.last_time = {}
        sys.settrace(None)


def show_results(prof, stream=None, precision=3):
    if stream is None:
        stream = sys.stdout
    template = '{0:>6} {1:>12} {2:>12}   {3:<}'

    for code in prof.code_map:
        lines = prof.code_map[code]
        if not lines:
            # .. measurements are empty ..
            continue
        filename = code.co_filename
        if filename.endswith((".pyc", ".pyo")):
            filename = filename[:-1]
        stream.write('Filename: ' + filename + '\n\n')
        if not os.path.exists(filename):
            stream.write('ERROR: Could not find file ' + filename + '\n')
            if filename.startswith("ipython-input") or filename.startswith("<ipython-input"):
                print("NOTE: %mprun can only be used on functions defined in "
                      "physical files, and not in the IPython environment.")
            continue
        all_lines = linecache.getlines(filename)
        sub_lines = inspect.getblock(all_lines[code.co_firstlineno - 1:])
        linenos = list(range(code.co_firstlineno, code.co_firstlineno +
                        len(sub_lines)))
        lines_normalized = {}

        header = template.format('Line #', 'Mem usage', 'Increment',
                                 'Line Contents')
        stream.write(header + '\n')
        stream.write('=' * len(header) + '\n')
        # move everything one frame up
        keys = sorted(lines.keys())

        k_old = keys[0] - 1
        lines_normalized[keys[0] - 1] = lines[keys[0]]
        for i in range(1, len(lines_normalized[keys[0] - 1])):
            lines_normalized[keys[0] - 1][i] = -1.
        k = keys.pop(0)
        while keys:
            lines_normalized[k] = lines[keys[0]]
            for i in range(len(lines_normalized[k_old]),
                           len(lines_normalized[k])):
                lines_normalized[k][i] = -1.
            k_old = k
            k = keys.pop(0)

        first_line = sorted(lines_normalized.keys())[0]
        mem_old = max(lines_normalized[first_line])
        precision = int(precision)
        template_mem = '{{0:{0}.{1}'.format(precision + 6, precision) + 'f} MB'
        for i, l in enumerate(linenos):
            mem = ''
            inc = ''
            if l in lines_normalized:
                mem = max(lines_normalized[l])
                inc = mem - mem_old
                mem_old = mem
                mem = template_mem.format(mem)
                inc = template_mem.format(inc)
            stream.write(template.format(l, mem, inc, sub_lines[i]))
        stream.write('\n\n')


# A lprun-style %mprun magic for IPython.
def magic_mprun(self, parameter_s=''):
    """ Execute a statement under the line-by-line memory profiler from the
    memory_profilser module.

    Usage:
      %mprun -f func1 -f func2 <statement>

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
    """
    try:
        from io import StringIO
    except ImportError: # Python 3.x
        from io import StringIO

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
    opts, arg_str = self.parse_options(parameter_s, 'rf:T:', list_all=True)
    opts.merge(opts_def)
    global_ns = self.shell.user_global_ns
    local_ns = self.shell.user_ns

    # Get the requested functions.
    funcs = []
    for name in opts.f:
        try:
            funcs.append(eval(name, global_ns, local_ns))
        except Exception as e:
            raise UsageError('Could not find function %r.\n%s: %s' % (name,
                e.__class__.__name__, e))

    profile = LineProfiler()
    for func in funcs:
        profile(func)

    # Add the profiler to the builtins for @profile.
    try:
        import builtins
    except ImportError:  # Python 3x
        import builtins as builtins

    if 'profile' in builtins.__dict__:
        had_profile = True
        old_profile = builtins.__dict__['profile']
    else:
        had_profile = False
        old_profile = None
    builtins.__dict__['profile'] = profile

    try:
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
    print((message,))

    text_file = opts.T[0]
    if text_file:
        with open(text_file, 'w') as pfile:
            pfile.write(output)
        print(('\n*** Profile printout saved to text file %s. %s' % (text_file,
                                                                    message)))

    return_value = None
    if 'r' in opts:
        return_value = profile

    return return_value


def _func_exec(stmt, ns):
    # helper for magic_memit, just a function proxy for the exec
    # statement
    exec(stmt, ns)

# a timeit-style %memit magic for IPython
def magic_memit(self, line=''):
    """Measure memory usage of a Python statement

    Usage, in line mode:
      %memit [-r<R>t<T>] statement

    Options:
    -r<R>: repeat the loop iteration <R> times and take the best result.
    Default: 1

    -t<T>: timeout after <T> seconds. Unused if `-i` is active. Default: None

    Examples
    --------
    ::

      In [1]: import numpy as np

      In [2]: %memit np.zeros(1e7)
      maximum of 1: 76.402344 MB per loop

      In [3]: %memit np.ones(1e6)
      maximum of 1: 7.820312 MB per loop

      In [4]: %memit -r 10 np.empty(1e8)
      maximum of 10: 0.101562 MB per loop

      In [5]: memit -t 3 while True: pass;
      Subprocess timed out.
      Subprocess timed out.
      Subprocess timed out.
      ERROR: all subprocesses exited unsuccessfully. Try again with the `-i`
      option.
      maximum of 1: -inf MB per loop

    """
    opts, stmt = self.parse_options(line, 'r:t:i', posix=False, strict=False)
    repeat = int(getattr(opts, 'r', 1))
    if repeat < 1:
        repeat == 1
    timeout = int(getattr(opts, 't', 0))
    if timeout <= 0:
        timeout = None

    mem_usage = []
    for _ in range(repeat):
        tmp = memory_usage((_func_exec, (stmt, self.shell.user_ns)), timeout=timeout)
        mem_usage.extend(tmp)

    if mem_usage:
        print(('maximum of %d: %f MB per loop' % (repeat, max(mem_usage))))
    else:
        print('ERROR: could not read memory usage, try with a lower interval or more iterations')


def load_ipython_extension(ip):
    """This is called to load the module as an IPython extension."""
    ip.define_magic('mprun', magic_mprun)
    ip.define_magic('memit', magic_memit)


def profile(func, stream=None):
    """
    Decorator that will run the function and print a line-by-line profile
    """
    def wrapper(*args, **kwargs):
        prof = LineProfiler()
        val = prof(func)(*args, **kwargs)
        show_results(prof, stream=stream)
        return val
    return wrapper


if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser(usage=_CMD_USAGE, version=__version__)
    parser.disable_interspersed_args()
    parser.add_option("--pdb-mmem", dest="max_mem", metavar="MAXMEM",
        type="float", action="store",
        help="step into the debugger when memory exceeds MAXMEM")
    parser.add_option('--precision', dest="precision", type="int",
        action="store", default=3,
        help="precision of memory output in number of significant digits")

    if not sys.argv[1:]:
        parser.print_help()
        sys.exit(2)

    (options, args) = parser.parse_args()

    prof = LineProfiler(max_mem=options.max_mem)
    __file__ = _find_script(args[0])
    try:
        if sys.version_info[0] < 3:
            import builtins
            builtins.__dict__['profile'] = prof
            ns = locals()
            ns['profile'] = prof # shadow the profile decorator defined above
            exec(compile(open(__file__).read(), __file__, 'exec'), ns, ns)
        else:
            import builtins
            builtins.__dict__['profile'] = prof
            ns = locals()
            ns['profile'] = prof # shadow the profile decorator defined above
            exec(compile(open(__file__).read(), __file__, 'exec'), ns,
                                                                   globals())
    finally:
        show_results(prof, precision=options.precision)
