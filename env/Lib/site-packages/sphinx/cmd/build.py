"""Build documentation from a provided source."""

from __future__ import annotations

import argparse
import bdb
import contextlib
import locale
import multiprocessing
import os
import pdb  # NoQA: T100
import sys
import traceback
from os import path
from typing import TYPE_CHECKING, Any, TextIO

from docutils.utils import SystemMessage

import sphinx.locale
from sphinx import __display_version__
from sphinx.application import Sphinx
from sphinx.errors import SphinxError, SphinxParallelError
from sphinx.locale import __
from sphinx.util._io import TeeStripANSI
from sphinx.util.console import color_terminal, nocolor, red, terminal_safe
from sphinx.util.docutils import docutils_namespace, patch_docutils
from sphinx.util.exceptions import format_exception_cut_frames, save_traceback
from sphinx.util.osutil import ensuredir

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Protocol

    class SupportsWrite(Protocol):
        def write(self, text: str, /) -> int | None: ...  # NoQA: E704


def handle_exception(
    app: Sphinx | None,
    args: Any,
    exception: BaseException,
    stderr: TextIO = sys.stderr,
) -> None:
    if isinstance(exception, bdb.BdbQuit):
        return

    if args.pdb:
        print(
            red(__('Exception occurred while building, starting debugger:')),
            file=stderr,
        )
        traceback.print_exc()
        pdb.post_mortem(sys.exc_info()[2])
    else:
        print(file=stderr)
        if args.verbosity or args.traceback:
            exc = sys.exc_info()[1]
            if isinstance(exc, SphinxParallelError):
                exc_format = '(Error in parallel process)\n' + exc.traceback
                print(exc_format, file=stderr)
            else:
                traceback.print_exc(None, stderr)
                print(file=stderr)
        if isinstance(exception, KeyboardInterrupt):
            print(__('Interrupted!'), file=stderr)
        elif isinstance(exception, SystemMessage):
            print(red(__('reST markup error:')), file=stderr)
            print(terminal_safe(exception.args[0]), file=stderr)
        elif isinstance(exception, SphinxError):
            print(red('%s:' % exception.category), file=stderr)
            print(str(exception), file=stderr)
        elif isinstance(exception, UnicodeError):
            print(red(__('Encoding error:')), file=stderr)
            print(terminal_safe(str(exception)), file=stderr)
            tbpath = save_traceback(app, exception)
            print(
                red(
                    __(
                        'The full traceback has been saved in %s, if you want '
                        'to report the issue to the developers.'
                    )
                    % tbpath
                ),
                file=stderr,
            )
        elif (
            isinstance(exception, RuntimeError)
            and 'recursion depth' in str(exception)
        ):  # fmt: skip
            print(red(__('Recursion error:')), file=stderr)
            print(terminal_safe(str(exception)), file=stderr)
            print(file=stderr)
            print(
                __(
                    'This can happen with very large or deeply nested source '
                    'files. You can carefully increase the default Python '
                    'recursion limit of 1000 in conf.py with e.g.:'
                ),
                file=stderr,
            )
            print('    import sys; sys.setrecursionlimit(1500)', file=stderr)
        else:
            print(red(__('Exception occurred:')), file=stderr)
            print(format_exception_cut_frames().rstrip(), file=stderr)
            tbpath = save_traceback(app, exception)
            print(
                red(
                    __(
                        'The full traceback has been saved in %s, if you '
                        'want to report the issue to the developers.'
                    )
                    % tbpath
                ),
                file=stderr,
            )
            print(
                __(
                    'Please also report this if it was a user error, so '
                    'that a better error message can be provided next time.'
                ),
                file=stderr,
            )
            print(
                __(
                    'A bug report can be filed in the tracker at '
                    '<https://github.com/sphinx-doc/sphinx/issues>. Thanks!'
                ),
                file=stderr,
            )


def jobs_argument(value: str) -> int:
    """
    Special type to handle 'auto' flags passed to 'sphinx-build' via -j flag. Can
    be expanded to handle other special scaling requests, such as setting job count
    to cpu_count.
    """
    if value == 'auto':
        return multiprocessing.cpu_count()
    else:
        jobs = int(value)
        if jobs <= 0:
            raise argparse.ArgumentTypeError(
                __('job number should be a positive number')
            )
        else:
            return jobs


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        usage='%(prog)s [OPTIONS] SOURCEDIR OUTPUTDIR [FILENAMES...]',
        epilog=__('For more information, visit <https://www.sphinx-doc.org/>.'),
        description=__("""
Generate documentation from source files.

sphinx-build generates documentation from the files in SOURCEDIR and places it
in OUTPUTDIR. It looks for 'conf.py' in SOURCEDIR for the configuration
settings. The 'sphinx-quickstart' tool may be used to generate template files,
including 'conf.py'

sphinx-build can create documentation in different formats. A format is
selected by specifying the builder name on the command line; it defaults to
HTML. Builders can also perform other tasks related to documentation
processing.

By default, everything that is outdated is built. Output only for selected
files can be built by specifying individual filenames.
"""),
    )

    parser.add_argument(
        '--version',
        action='version',
        dest='show_version',
        version=f'%(prog)s {__display_version__}',
    )

    parser.add_argument(
        'sourcedir', metavar='SOURCE_DIR', help=__('path to documentation source files')
    )
    parser.add_argument(
        'outputdir', metavar='OUTPUT_DIR', help=__('path to output directory')
    )
    parser.add_argument(
        'filenames',
        nargs='*',
        help=__(
            '(optional) a list of specific files to rebuild. '
            'Ignored if --write-all is specified'
        ),
    )

    group = parser.add_argument_group(__('general options'))
    group.add_argument(
        '--builder',
        '-b',
        metavar='BUILDER',
        dest='builder',
        default='html',
        help=__("builder to use (default: 'html')"),
    )
    group.add_argument(
        '--jobs',
        '-j',
        metavar='N',
        default=1,
        type=jobs_argument,
        dest='jobs',
        help=__(
            'run in parallel with N processes, when possible. '
            "'auto' uses the number of CPU cores"
        ),
    )
    group.add_argument(
        '--write-all',
        '-a',
        action='store_true',
        dest='force_all',
        help=__('write all files (default: only write new and ' 'changed files)'),
    )
    group.add_argument(
        '--fresh-env',
        '-E',
        action='store_true',
        dest='freshenv',
        help=__("don't use a saved environment, always read " 'all files'),
    )

    group = parser.add_argument_group(__('path options'))
    group.add_argument(
        '--doctree-dir',
        '-d',
        metavar='PATH',
        dest='doctreedir',
        help=__(
            'directory for doctree and environment files '
            '(default: OUTPUT_DIR/.doctrees)'
        ),
    )
    group.add_argument(
        '--conf-dir',
        '-c',
        metavar='PATH',
        dest='confdir',
        help=__(
            'directory for the configuration file (conf.py) ' '(default: SOURCE_DIR)'
        ),
    )

    group = parser.add_argument_group('build configuration options')
    group.add_argument(
        '--isolated',
        '-C',
        action='store_true',
        dest='noconfig',
        help=__('use no configuration file, only use settings from -D options'),
    )
    group.add_argument(
        '--define',
        '-D',
        metavar='setting=value',
        action='append',
        dest='define',
        default=[],
        help=__('override a setting in configuration file'),
    )
    group.add_argument(
        '--html-define',
        '-A',
        metavar='name=value',
        action='append',
        dest='htmldefine',
        default=[],
        help=__('pass a value into HTML templates'),
    )
    group.add_argument(
        '--tag',
        '-t',
        metavar='TAG',
        action='append',
        dest='tags',
        default=[],
        help=__('define tag: include "only" blocks with TAG'),
    )
    group.add_argument(
        '--nitpicky',
        '-n',
        action='store_true',
        dest='nitpicky',
        help=__('nitpicky mode: warn about all missing references'),
    )

    group = parser.add_argument_group(__('console output options'))
    group.add_argument(
        '--verbose',
        '-v',
        action='count',
        dest='verbosity',
        default=0,
        help=__('increase verbosity (can be repeated)'),
    )
    group.add_argument(
        '--quiet',
        '-q',
        action='store_true',
        dest='quiet',
        help=__('no output on stdout, just warnings on stderr'),
    )
    group.add_argument(
        '--silent',
        '-Q',
        action='store_true',
        dest='really_quiet',
        help=__('no output at all, not even warnings'),
    )
    group.add_argument(
        '--color',
        action='store_const',
        dest='color',
        const='yes',
        default='auto',
        help=__('do emit colored output (default: auto-detect)'),
    )
    group.add_argument(
        '--no-color',
        '-N',
        action='store_const',
        dest='color',
        const='no',
        help=__('do not emit colored output (default: auto-detect)'),
    )

    group = parser.add_argument_group(__('warning control options'))
    group.add_argument(
        '--warning-file',
        '-w',
        metavar='FILE',
        dest='warnfile',
        help=__('write warnings (and errors) to given file'),
    )
    group.add_argument(
        '--fail-on-warning',
        '-W',
        action='store_true',
        dest='warningiserror',
        help=__('turn warnings into errors'),
    )
    group.add_argument('--keep-going', action='store_true', help=argparse.SUPPRESS)
    group.add_argument(
        '--show-traceback',
        '-T',
        action='store_true',
        dest='traceback',
        help=__('show full traceback on exception'),
    )
    group.add_argument(
        '--pdb', '-P', action='store_true', dest='pdb', help=__('run Pdb on exception')
    )
    group.add_argument(
        '--exception-on-warning',
        action='store_true',
        dest='exception_on_warning',
        help=__('raise an exception on warnings'),
    )

    if parser.prog == '__main__.py':
        parser.prog = 'sphinx-build'

    return parser


def make_main(argv: Sequence[str]) -> int:
    """Sphinx build "make mode" entry."""
    from sphinx.cmd import make_mode

    return make_mode.run_make_mode(argv[1:])


def _parse_arguments(
    parser: argparse.ArgumentParser, argv: Sequence[str]
) -> argparse.Namespace:
    args = parser.parse_args(argv)
    return args


def _parse_confdir(noconfig: bool, confdir: str, sourcedir: str) -> str | None:
    if noconfig:
        return None
    elif not confdir:
        return sourcedir
    return confdir


def _parse_doctreedir(doctreedir: str, outputdir: str) -> str:
    if doctreedir:
        return doctreedir
    return os.path.join(outputdir, '.doctrees')


def _validate_filenames(
    parser: argparse.ArgumentParser,
    force_all: bool,
    filenames: list[str],
) -> None:
    if force_all and filenames:
        parser.error(__('cannot combine -a option and filenames'))


def _validate_colour_support(colour: str) -> None:
    if colour == 'no' or (colour == 'auto' and not color_terminal()):
        nocolor()


def _parse_logging(
    parser: argparse.ArgumentParser,
    quiet: bool,
    really_quiet: bool,
    warnfile: str | None,
) -> tuple[TextIO | None, TextIO | None, TextIO, TextIO | None]:
    status: TextIO | None = sys.stdout
    warning: TextIO | None = sys.stderr
    error = sys.stderr

    if quiet:
        status = None

    if really_quiet:
        status = warning = None

    warnfp = None
    if warning and warnfile:
        try:
            warnfile = path.abspath(warnfile)
            ensuredir(path.dirname(warnfile))
            # the caller is responsible for closing this file descriptor
            warnfp = open(warnfile, 'w', encoding='utf-8')  # NoQA: SIM115
        except Exception as exc:
            parser.error(__('cannot open warning file %r: %s') % (warnfile, exc))
        warning = TeeStripANSI(warning, warnfp)  # type: ignore[assignment]
        error = warning

    return status, warning, error, warnfp


def _parse_confoverrides(
    parser: argparse.ArgumentParser,
    define: list[str],
    htmldefine: list[str],
    nitpicky: bool,
) -> dict[str, Any]:
    confoverrides: dict[str, Any] = {}
    val: Any
    for val in define:
        try:
            key, val = val.split('=', 1)
        except ValueError:
            parser.error(__('-D option argument must be in the form name=value'))
        confoverrides[key] = val

    for val in htmldefine:
        try:
            key, val = val.split('=')
        except ValueError:
            parser.error(__('-A option argument must be in the form name=value'))
        with contextlib.suppress(ValueError):
            val = int(val)

        confoverrides[f'html_context.{key}'] = val

    if nitpicky:
        confoverrides['nitpicky'] = True

    return confoverrides


def build_main(argv: Sequence[str]) -> int:
    """Sphinx build "main" command-line entry."""
    parser = get_parser()
    args = _parse_arguments(parser, argv)
    args.confdir = _parse_confdir(args.noconfig, args.confdir, args.sourcedir)
    args.doctreedir = _parse_doctreedir(args.doctreedir, args.outputdir)
    _validate_filenames(parser, args.force_all, args.filenames)
    _validate_colour_support(args.color)
    args.status, args.warning, args.error, warnfp = _parse_logging(
        parser, args.quiet, args.really_quiet, args.warnfile
    )
    args.confoverrides = _parse_confoverrides(
        parser, args.define, args.htmldefine, args.nitpicky
    )

    app = None
    try:
        confdir = args.confdir or args.sourcedir
        with patch_docutils(confdir), docutils_namespace():
            app = Sphinx(
                srcdir=args.sourcedir,
                confdir=args.confdir,
                outdir=args.outputdir,
                doctreedir=args.doctreedir,
                buildername=args.builder,
                confoverrides=args.confoverrides,
                status=args.status,
                warning=args.warning,
                freshenv=args.freshenv,
                warningiserror=args.warningiserror,
                tags=args.tags,
                verbosity=args.verbosity,
                parallel=args.jobs,
                keep_going=False,
                pdb=args.pdb,
                exception_on_warning=args.exception_on_warning,
            )
            app.build(args.force_all, args.filenames)
            return app.statuscode
    except (Exception, KeyboardInterrupt) as exc:
        handle_exception(app, args, exc, args.error)
        return 2
    finally:
        if warnfp is not None:
            # close the file descriptor for the warnings file opened by Sphinx
            warnfp.close()


def _bug_report_info() -> int:
    from platform import platform, python_implementation

    import docutils
    import jinja2
    import pygments

    print('Please paste all output below into the bug report template\n\n')
    print('```text')
    print(f'Platform:              {sys.platform}; ({platform()})')
    print(f'Python version:        {sys.version})')
    print(f'Python implementation: {python_implementation()}')
    print(f'Sphinx version:        {sphinx.__display_version__}')
    print(f'Docutils version:      {docutils.__version__}')
    print(f'Jinja2 version:        {jinja2.__version__}')
    print(f'Pygments version:      {pygments.__version__}')
    print('```')
    return 0


def main(argv: Sequence[str] = (), /) -> int:
    locale.setlocale(locale.LC_ALL, '')
    sphinx.locale.init_console()

    if not argv:
        argv = sys.argv[1:]

    # Allow calling as 'python -m sphinx build …'
    if argv[:1] == ['build']:
        argv = argv[1:]

    if argv[:1] == ['--bug-report']:
        return _bug_report_info()
    if argv[:1] == ['-M']:
        from sphinx.cmd import make_mode

        return make_mode.run_make_mode(argv[1:])
    else:
        return build_main(argv)


if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1:]))
