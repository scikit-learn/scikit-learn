# -*- coding: utf-8 -*-
"""
    sphinx.cmdline
    ~~~~~~~~~~~~~~

    sphinx-build command-line handling.

    :copyright: Copyright 2007-2018 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""
from __future__ import print_function

import argparse
import multiprocessing
import sys
import traceback
from os import path

from docutils.utils import SystemMessage
from six import text_type, binary_type

from sphinx import __display_version__
from sphinx.application import Sphinx
from sphinx.errors import SphinxError
from sphinx.util import Tee, format_exception_cut_frames, save_traceback
from sphinx.util.console import red, nocolor, color_terminal  # type: ignore
from sphinx.util.docutils import docutils_namespace, patch_docutils
from sphinx.util.osutil import abspath, fs_encoding
from sphinx.util.pycompat import terminal_safe

if False:
    # For type annotation
    from typing import Any, IO, List, Union  # NOQA


def handle_exception(app, args, exception, stderr=sys.stderr):
    # type: (Sphinx, Any, Union[Exception, KeyboardInterrupt], IO) -> None
    if args.pdb:
        import pdb
        print(red('Exception occurred while building, starting debugger:'),
              file=stderr)
        traceback.print_exc()
        pdb.post_mortem(sys.exc_info()[2])
    else:
        print(file=stderr)
        if args.verbosity or args.traceback:
            traceback.print_exc(None, stderr)
            print(file=stderr)
        if isinstance(exception, KeyboardInterrupt):
            print('interrupted!', file=stderr)
        elif isinstance(exception, SystemMessage):
            print(red('reST markup error:'), file=stderr)
            print(terminal_safe(exception.args[0]), file=stderr)
        elif isinstance(exception, SphinxError):
            print(red('%s:' % exception.category), file=stderr)
            print(terminal_safe(text_type(exception)), file=stderr)
        elif isinstance(exception, UnicodeError):
            print(red('Encoding error:'), file=stderr)
            print(terminal_safe(text_type(exception)), file=stderr)
            tbpath = save_traceback(app)
            print(red('The full traceback has been saved in %s, if you want '
                      'to report the issue to the developers.' % tbpath),
                  file=stderr)
        elif isinstance(exception, RuntimeError) and 'recursion depth' in str(exception):
            print(red('Recursion error:'), file=stderr)
            print(terminal_safe(text_type(exception)), file=stderr)
            print(file=stderr)
            print('This can happen with very large or deeply nested source '
                  'files.  You can carefully increase the default Python '
                  'recursion limit of 1000 in conf.py with e.g.:', file=stderr)
            print('    import sys; sys.setrecursionlimit(1500)', file=stderr)
        else:
            print(red('Exception occurred:'), file=stderr)
            print(format_exception_cut_frames().rstrip(), file=stderr)
            tbpath = save_traceback(app)
            print(red('The full traceback has been saved in %s, if you '
                      'want to report the issue to the developers.' % tbpath),
                  file=stderr)
            print('Please also report this if it was a user error, so '
                  'that a better error message can be provided next time.',
                  file=stderr)
            print('A bug report can be filed in the tracker at '
                  '<https://github.com/sphinx-doc/sphinx/issues>. Thanks!',
                  file=stderr)


def jobs_argument(value):
    # type: (str) -> int
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
            raise argparse.ArgumentTypeError('job number should be a positive number')
        else:
            return jobs


def get_parser():
    # type: () -> argparse.ArgumentParser
    parser = argparse.ArgumentParser(
        usage='%(prog)s [OPTIONS] SOURCEDIR OUTPUTDIR [FILENAMES...]',
        epilog='For more information, visit <http://sphinx-doc.org/>.',
        description="""
Generate documentation from source files.

sphinx-build generates documentation from the files in SOURCEDIR and places it
in OUTPUTDIR. It looks for 'conf.py' in SOURCEDIR for the configuration
settings.  The 'sphinx-quickstart' tool may be used to generate template files,
including 'conf.py'

sphinx-build can create documentation in different formats. A format is
selected by specifying the builder name on the command line; it defaults to
HTML. Builders can also perform other tasks related to documentation
processing.

By default, everything that is outdated is built. Output only for selected
files can be built by specifying individual filenames.
""")

    parser.add_argument('--version', action='version', dest='show_version',
                        version='%%(prog)s %s' % __display_version__)

    parser.add_argument('sourcedir',
                        help='path to documentation source files')
    parser.add_argument('outputdir',
                        help='path to output directory')
    parser.add_argument('filenames', nargs='*',
                        help='a list of specific files to rebuild. Ignored '
                        'if -a is specified')

    group = parser.add_argument_group('general options')
    group.add_argument('-b', metavar='BUILDER', dest='builder',
                       default='html',
                       help='builder to use (default: html)')
    group.add_argument('-a', action='store_true', dest='force_all',
                       help='write all files (default: only write new and '
                       'changed files)')
    group.add_argument('-E', action='store_true', dest='freshenv',
                       help='don\'t use a saved environment, always read '
                       'all files')
    group.add_argument('-d', metavar='PATH', dest='doctreedir',
                       help='path for the cached environment and doctree '
                       'files (default: OUTPUTDIR/.doctrees)')
    group.add_argument('-j', metavar='N', default=1, type=jobs_argument, dest='jobs',
                       help='build in parallel with N processes where '
                       'possible (special value "auto" will set N to cpu-count)')
    group = parser.add_argument_group('build configuration options')
    group.add_argument('-c', metavar='PATH', dest='confdir',
                       help='path where configuration file (conf.py) is '
                       'located (default: same as SOURCEDIR)')
    group.add_argument('-C', action='store_true', dest='noconfig',
                       help='use no config file at all, only -D options')
    group.add_argument('-D', metavar='setting=value', action='append',
                       dest='define', default=[],
                       help='override a setting in configuration file')
    group.add_argument('-A', metavar='name=value', action='append',
                       dest='htmldefine', default=[],
                       help='pass a value into HTML templates')
    group.add_argument('-t', metavar='TAG', action='append',
                       dest='tags', default=[],
                       help='define tag: include "only" blocks with TAG')
    group.add_argument('-n', action='store_true', dest='nitpicky',
                       help='nit-picky mode, warn about all missing '
                       'references')

    group = parser.add_argument_group('console output options')
    group.add_argument('-v', action='count', dest='verbosity', default=0,
                       help='increase verbosity (can be repeated)')
    group.add_argument('-q', action='store_true', dest='quiet',
                       help='no output on stdout, just warnings on stderr')
    group.add_argument('-Q', action='store_true', dest='really_quiet',
                       help='no output at all, not even warnings')
    group.add_argument('--color', action='store_const', const='yes',
                       default='auto',
                       help='do emit colored output (default: auto-detect)')
    group.add_argument('-N', '--no-color', dest='color', action='store_const',
                       const='no',
                       help='do not emit colored output (default: '
                       'auto-detect)')
    group.add_argument('-w', metavar='FILE', dest='warnfile',
                       help='write warnings (and errors) to given file')
    group.add_argument('-W', action='store_true', dest='warningiserror',
                       help='turn warnings into errors')
    group.add_argument('-T', action='store_true', dest='traceback',
                       help='show full traceback on exception')
    group.add_argument('-P', action='store_true', dest='pdb',
                       help='run Pdb on exception')

    return parser


def main(argv=sys.argv[1:]):  # type: ignore
    # type: (List[unicode]) -> int

    parser = get_parser()
    args = parser.parse_args(argv)

    # get paths (first and second positional argument)
    try:
        srcdir = abspath(args.sourcedir)
        confdir = abspath(args.confdir or srcdir)
        if args.noconfig:
            confdir = None

        if not path.isdir(srcdir):
            parser.error('cannot find source directory (%s)' % srcdir)
        if not args.noconfig and not path.isfile(path.join(confdir, 'conf.py')):
            parser.error("config directory doesn't contain a conf.py file "
                         "(%s)" % confdir)

        outdir = abspath(args.outputdir)
        if srcdir == outdir:
            parser.error('source directory and destination directory are same')
    except UnicodeError:
        parser.error('multibyte filename not supported on this filesystem '
                     'encoding (%r)' % fs_encoding)

    # handle remaining filename arguments
    filenames = args.filenames
    missing_files = []
    for filename in filenames:
        if not path.isfile(filename):
            missing_files.append(filename)
    if missing_files:
        parser.error('cannot find files %r' % missing_files)

    # likely encoding used for command-line arguments
    try:
        locale = __import__('locale')  # due to submodule of the same name
        likely_encoding = locale.getpreferredencoding()
    except Exception:
        likely_encoding = None

    if args.force_all and filenames:
        parser.error('cannot combine -a option and filenames')

    if args.color == 'no' or (args.color == 'auto' and not color_terminal()):
        nocolor()

    doctreedir = abspath(args.doctreedir or path.join(outdir, '.doctrees'))

    status = sys.stdout
    warning = sys.stderr
    error = sys.stderr

    if args.quiet:
        status = None

    if args.really_quiet:
        status = warning = None

    if warning and args.warnfile:
        try:
            warnfp = open(args.warnfile, 'w')
        except Exception as exc:
            parser.error('cannot open warning file %r: %s' % (
                args.warnfile, exc))
        warning = Tee(warning, warnfp)  # type: ignore
        error = warning

    confoverrides = {}
    for val in args.define:
        try:
            key, val = val.split('=', 1)
        except ValueError:
            parser.error('-D option argument must be in the form name=value')
        if likely_encoding and isinstance(val, binary_type):
            try:
                val = val.decode(likely_encoding)
            except UnicodeError:
                pass
        confoverrides[key] = val

    for val in args.htmldefine:
        try:
            key, val = val.split('=')
        except ValueError:
            parser.error('-A option argument must be in the form name=value')
        try:
            val = int(val)
        except ValueError:
            if likely_encoding and isinstance(val, binary_type):
                try:
                    val = val.decode(likely_encoding)
                except UnicodeError:
                    pass
        confoverrides['html_context.%s' % key] = val

    if args.nitpicky:
        confoverrides['nitpicky'] = True

    app = None
    try:
        with patch_docutils(), docutils_namespace():
            app = Sphinx(srcdir, confdir, outdir, doctreedir, args.builder,
                         confoverrides, status, warning, args.freshenv,
                         args.warningiserror, args.tags, args.verbosity, args.jobs)
            app.build(args.force_all, filenames)
            return app.statuscode
    except (Exception, KeyboardInterrupt) as exc:
        handle_exception(app, args, exc, error)
        return 2
