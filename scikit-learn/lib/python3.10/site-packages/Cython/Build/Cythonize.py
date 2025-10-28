import concurrent.futures
import os
import shutil
import sys
import tempfile
from collections import defaultdict
from contextlib import contextmanager

from .Dependencies import cythonize, extended_iglob
from ..Utils import is_package_dir
from ..Compiler import Options

try:
    import multiprocessing
    parallel_compiles = int(multiprocessing.cpu_count() * 1.5)
except ImportError:
    multiprocessing = None
    parallel_compiles = 0


def find_package_base(path):
    base_dir, package_path = os.path.split(path)
    while is_package_dir(base_dir):
        base_dir, parent = os.path.split(base_dir)
        package_path = '%s/%s' % (parent, package_path)
    return base_dir, package_path


def cython_compile(path_pattern, options) -> dict:
    all_paths = map(os.path.abspath, extended_iglob(path_pattern))
    ext_modules_by_basedir = _cython_compile_files(all_paths, options)
    _build(list(ext_modules_by_basedir.items()), options.parallel)


def _cython_compile_files(all_paths, options) -> dict:
    ext_modules_to_build = defaultdict(list)

    for path in all_paths:
        if options.build_inplace:
            base_dir = path
            while not os.path.isdir(base_dir) or is_package_dir(base_dir):
                base_dir = os.path.dirname(base_dir)
        else:
            base_dir = None

        if os.path.isdir(path):
            # recursively compiling a package
            paths = [os.path.join(path, '**', '*.{py,pyx}')]
        else:
            # assume it's a file(-like thing)
            paths = [path]

        ext_modules = cythonize(
            paths,
            nthreads=options.parallel,
            exclude_failures=options.keep_going,
            exclude=options.excludes,
            compiler_directives=options.directives,
            compile_time_env=options.compile_time_env,
            force=options.force,
            quiet=options.quiet,
            depfile=options.depfile,
            language=options.language,
            **options.options)

        if ext_modules and options.build:
            ext_modules_to_build[base_dir].extend(ext_modules)

    return dict(ext_modules_to_build)


@contextmanager
def _interruptible_pool(pool_cm):
    with pool_cm as proc_pool:
        try:
            yield proc_pool
        except KeyboardInterrupt:
            proc_pool.terminate_workers()
            proc_pool.shutdown(cancel_futures=True)
            raise


def _build(ext_modules, parallel):
    modcount = sum(len(modules) for _, modules in ext_modules)
    if not modcount:
        return

    serial_execution_mode = modcount == 1 or parallel < 2

    try:
        pool_cm = (
            None if serial_execution_mode
            else concurrent.futures.ProcessPoolExecutor(max_workers=parallel)
        )
    except (OSError, ImportError):
        # `OSError` is a historic exception in `multiprocessing`
        # `ImportError` happens e.g. under pyodide (`ModuleNotFoundError`)
        serial_execution_mode = True

    if serial_execution_mode:
        for ext in ext_modules:
            run_distutils(ext)
        return

    with _interruptible_pool(pool_cm) as proc_pool:
        compiler_tasks = [
            proc_pool.submit(run_distutils, (base_dir, [ext]))
            for base_dir, modules in ext_modules
            for ext in modules
        ]

        concurrent.futures.wait(compiler_tasks, return_when=concurrent.futures.FIRST_EXCEPTION)

        worker_exceptions = []
        for task in compiler_tasks:  # discover any crashes
            try:
                task.result()
            except BaseException as proc_err:  # could be SystemExit
                worker_exceptions.append(proc_err)

        if worker_exceptions:
            exc_msg = 'Compiling Cython modules failed with these errors:\n\n'
            exc_msg += '\n\t* '.join(('', *map(str, worker_exceptions)))
            exc_msg += '\n\n'

            non_base_exceptions = [
                exc for exc in worker_exceptions
                if isinstance(exc, Exception)
            ]
            if sys.version_info[:2] >= (3, 11) and non_base_exceptions:
                raise ExceptionGroup(exc_msg, non_base_exceptions)
            else:
                raise RuntimeError(exc_msg) from worker_exceptions[0]


def run_distutils(args):
    try:
        from distutils.core import setup
    except ImportError:
        try:
            from setuptools import setup
        except ImportError:
            raise ImportError("'distutils' is not available. Please install 'setuptools' for binary builds.")

    base_dir, ext_modules = args
    script_args = ['build_ext', '-i']
    cwd = os.getcwd()
    temp_dir = None
    try:
        if base_dir:
            os.chdir(base_dir)
            temp_dir = tempfile.mkdtemp(dir=base_dir)
            script_args.extend(['--build-temp', temp_dir])
        setup(
            script_name='setup.py',
            script_args=script_args,
            ext_modules=ext_modules,
        )
    finally:
        if base_dir:
            os.chdir(cwd)
            if temp_dir and os.path.isdir(temp_dir):
                shutil.rmtree(temp_dir)


def benchmark(code, setup_code=None, import_module=None, directives=None):
    from Cython.Build.Inline import cymeit

    timings, number = cymeit(code, setup_code, import_module, directives, repeat=9)

    # Based on 'timeit.main()' in CPython 3.13.
    units = {"nsec": 1e-9, "usec": 1e-6, "msec": 1e-3, "sec": 1.0}
    scales = [(scale, unit) for unit, scale in reversed(units.items())]  # biggest first

    def format_time(t):
        for scale, unit in scales:
            if t >= scale:
                break
        else:
            raise RuntimeError("Timing is below nanoseconds: {t:f}")
        return f"{t / scale :.3f} {unit}"

    timings.sort()
    assert len(timings) & 1 == 1  # odd number of timings, for median position
    fastest, median, slowest = timings[0], timings[len(timings) // 2], timings[-1]

    print(f"{number} loops, best of {len(timings)}: {format_time(fastest)} per loop (median: {format_time(median)})")

    if slowest > fastest * 4:
        print(
            "The timings are likely unreliable. "
            f"The worst time ({format_time(slowest)}) was more than four times "
            f"slower than the best time ({format_time(fastest)}).")


def create_args_parser():
    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    from ..Compiler.CmdLine import ParseDirectivesAction, ParseOptionsAction, ParseCompileTimeEnvAction

    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        epilog="""\
Environment variables:
  CYTHON_FORCE_REGEN: if set to 1, forces cythonize to regenerate the output files regardless
        of modification times and changes.
  CYTHON_CACHE_DIR: the base directory containing Cython's caches.
  Environment variables accepted by setuptools are supported to configure the C compiler and build:
  https://setuptools.pypa.io/en/latest/userguide/ext_modules.html#compiler-and-linker-options"""
    )

    parser.add_argument('-X', '--directive', metavar='NAME=VALUE,...',
                      dest='directives', default={}, type=str,
                      action=ParseDirectivesAction,
                      help='set a compiler directive')
    parser.add_argument('-E', '--compile-time-env', metavar='NAME=VALUE,...',
                      dest='compile_time_env', default={}, type=str,
                      action=ParseCompileTimeEnvAction,
                      help='set a compile time environment variable')
    parser.add_argument('-s', '--option', metavar='NAME=VALUE',
                      dest='options', default={}, type=str,
                      action=ParseOptionsAction,
                      help='set a cythonize option')
    parser.add_argument('-2', dest='language_level', action='store_const', const=2, default=None,
                      help='use Python 2 syntax mode by default')
    parser.add_argument('-3', dest='language_level', action='store_const', const=3,
                      help='use Python 3 syntax mode by default')
    parser.add_argument('--3str', dest='language_level', action='store_const', const=3,
                      help='use Python 3 syntax mode by default (deprecated alias for -3)')
    parser.add_argument('-+', '--cplus', dest='language', action='store_const', const='c++', default=None,
                        help='Compile as C++ rather than C')
    parser.add_argument('-a', '--annotate', action='store_const', const='default', dest='annotate',
                      help='Produce a colorized HTML version of the source.')
    parser.add_argument('--annotate-fullc', action='store_const', const='fullc', dest='annotate',
                      help='Produce a colorized HTML version of the source '
                           'which includes entire generated C/C++-code.')
    parser.add_argument('-x', '--exclude', metavar='PATTERN', dest='excludes',
                      action='append', default=[],
                      help='exclude certain file patterns from the compilation')

    parser.add_argument('-b', '--build', dest='build', action='store_true', default=None,
                      help='build extension modules using distutils/setuptools')
    parser.add_argument('-i', '--inplace', dest='build_inplace', action='store_true', default=None,
                      help='build extension modules in place using distutils/setuptools (implies -b)')

    parser.add_argument('--timeit', dest='benchmark', metavar="CODESTRING", type=str, default=None,
                      help="build in place, then compile+run CODESTRING as benchmark in first module's namespace (implies -i)")
    parser.add_argument('--setup', dest='benchmark_setup', metavar="CODESTRING", type=str, default=None,
                      help="use CODESTRING as pre-benchmark setup code for --bench")

    parser.add_argument('-j', '--parallel', dest='parallel', metavar='N',
                      type=int, default=parallel_compiles,
                      help=f'run builds in N parallel jobs (default: {parallel_compiles or 1})')
    parser.add_argument('-f', '--force', dest='force', action='store_true', default=None,
                      help='force recompilation')
    parser.add_argument('-q', '--quiet', dest='quiet', action='store_true', default=None,
                      help='be less verbose during compilation')

    parser.add_argument('--lenient', dest='lenient', action='store_true', default=None,
                      help='increase Python compatibility by ignoring some compile time errors')
    parser.add_argument('-k', '--keep-going', dest='keep_going', action='store_true', default=None,
                      help='compile as much as possible, ignore compilation failures')
    parser.add_argument('--no-docstrings', dest='no_docstrings', action='store_true', default=None,
                      help='strip docstrings')
    parser.add_argument('-M', '--depfile', action='store_true', help='produce depfiles for the sources')
    parser.add_argument('sources', nargs='*')
    return parser


def parse_args_raw(parser, args):
    options, unknown = parser.parse_known_args(args)
    sources = options.sources
    # if positional arguments were interspersed
    # some of them are in unknown
    for option in unknown:
        if option.startswith('-'):
            parser.error("unknown option "+option)
        else:
            sources.append(option)
    del options.sources
    return (options, sources)


def parse_args(args):
    parser = create_args_parser()
    options, args = parse_args_raw(parser, args)

    if options.benchmark is not None:
        options.build_inplace = True
    elif not args:
        parser.error("no source files provided")

    if options.build_inplace:
        options.build = True
    if multiprocessing is None:
        options.parallel = 0
    if options.language_level:
        assert options.language_level in (2, 3, '3str')
        options.options['language_level'] = options.language_level

    if options.lenient:
        # increase Python compatibility by ignoring compile time errors
        Options.error_on_unknown_names = False
        Options.error_on_uninitialized = False

    if options.annotate:
        Options.annotate = options.annotate

    if options.no_docstrings:
        Options.docstrings = False

    return options, args


def main(args=None):
    options, paths = parse_args(args)

    all_paths = []
    for path in paths:
        expanded_path = [os.path.abspath(p) for p in extended_iglob(path)]
        if not expanded_path:
            print("{}: No such file or directory: '{}'".format(sys.argv[0], path), file=sys.stderr)
            sys.exit(1)
        all_paths.extend(expanded_path)

    ext_modules_by_basedir = _cython_compile_files(all_paths, options)

    if ext_modules_by_basedir and options.build:
        _build(list(ext_modules_by_basedir.items()), options.parallel)

    if options.benchmark is not None:
        base_dir = import_module = None
        if ext_modules_by_basedir:
            base_dir, first_extensions = ext_modules_by_basedir.popitem()
            if first_extensions:
                import_module = first_extensions[0].name

        if base_dir is not None:
            sys.path.insert(0, base_dir)

        benchmark(
            options.benchmark, options.benchmark_setup,
            import_module=import_module,
        )

        if base_dir is not None:
            sys.path.remove(base_dir)


if __name__ == '__main__':
    main()
