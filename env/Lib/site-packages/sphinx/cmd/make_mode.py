"""sphinx-build -M command-line handling.

This replaces the old, platform-dependent and once-generated content
of Makefile / make.bat.

This is in its own module so that importing it is fast.  It should not
import the main Sphinx modules (like sphinx.applications, sphinx.builders).
"""

from __future__ import annotations

import os
import subprocess
import sys
from os import path
from typing import TYPE_CHECKING

import sphinx
from sphinx.cmd.build import build_main
from sphinx.util.console import blue, bold, color_terminal, nocolor
from sphinx.util.osutil import rmtree

if sys.version_info >= (3, 11):
    from contextlib import chdir
else:
    from sphinx.util.osutil import _chdir as chdir

if TYPE_CHECKING:
    from collections.abc import Sequence

BUILDERS = [
    ('', 'html', 'to make standalone HTML files'),
    ('', 'dirhtml', 'to make HTML files named index.html in directories'),
    ('', 'singlehtml', 'to make a single large HTML file'),
    ('', 'pickle', 'to make pickle files'),
    ('', 'json', 'to make JSON files'),
    ('', 'htmlhelp', 'to make HTML files and an HTML help project'),
    ('', 'qthelp', 'to make HTML files and a qthelp project'),
    ('', 'devhelp', 'to make HTML files and a Devhelp project'),
    ('', 'epub', 'to make an epub'),
    ('', 'latex', 'to make LaTeX files, you can set PAPER=a4 or PAPER=letter'),
    ('posix', 'latexpdf', 'to make LaTeX and PDF files (default pdflatex)'),
    ('posix', 'latexpdfja', 'to make LaTeX files and run them through platex/dvipdfmx'),
    ('', 'text', 'to make text files'),
    ('', 'man', 'to make manual pages'),
    ('', 'texinfo', 'to make Texinfo files'),
    ('posix', 'info', 'to make Texinfo files and run them through makeinfo'),
    ('', 'gettext', 'to make PO message catalogs'),
    ('', 'changes', 'to make an overview of all changed/added/deprecated items'),
    ('', 'xml', 'to make Docutils-native XML files'),
    ('', 'pseudoxml', 'to make pseudoxml-XML files for display purposes'),
    ('', 'linkcheck', 'to check all external links for integrity'),
    (
        '',
        'doctest',
        'to run all doctests embedded in the documentation ' '(if enabled)',
    ),
    ('', 'coverage', 'to run coverage check of the documentation (if enabled)'),
    ('', 'clean', 'to remove everything in the build directory'),
]


class Make:
    def __init__(self, *, source_dir: str, build_dir: str, opts: Sequence[str]) -> None:
        self.source_dir = source_dir
        self.build_dir = build_dir
        self.opts = [*opts]

    def build_dir_join(self, *comps: str) -> str:
        return path.join(self.build_dir, *comps)

    def build_clean(self) -> int:
        source_dir = path.abspath(self.source_dir)
        build_dir = path.abspath(self.build_dir)
        if not path.exists(self.build_dir):
            return 0
        elif not path.isdir(self.build_dir):
            print('Error: %r is not a directory!' % self.build_dir)
            return 1
        elif source_dir == build_dir:
            print('Error: %r is same as source directory!' % self.build_dir)
            return 1
        elif path.commonpath([source_dir, build_dir]) == build_dir:
            print('Error: %r directory contains source directory!' % self.build_dir)
            return 1
        print('Removing everything under %r...' % self.build_dir)
        for item in os.listdir(self.build_dir):
            rmtree(self.build_dir_join(item))
        return 0

    def build_help(self) -> None:
        if not color_terminal():
            nocolor()

        print(bold('Sphinx v%s' % sphinx.__display_version__))
        print("Please use `make %s' where %s is one of" % ((blue('target'),) * 2))
        for osname, bname, description in BUILDERS:
            if not osname or os.name == osname:
                print(f'  {blue(bname.ljust(10))}  {description}')

    def build_latexpdf(self) -> int:
        if self.run_generic_build('latex') > 0:
            return 1

        # Use $MAKE to determine the make command
        make_fallback = 'make.bat' if sys.platform == 'win32' else 'make'
        makecmd = os.environ.get('MAKE', make_fallback)
        if not makecmd.lower().startswith('make'):
            raise RuntimeError('Invalid $MAKE command: %r' % makecmd)
        try:
            with chdir(self.build_dir_join('latex')):
                if '-Q' in self.opts:
                    with open('__LATEXSTDOUT__', 'w') as outfile:
                        returncode = subprocess.call(
                            [
                                makecmd,
                                'all-pdf',
                                'LATEXOPTS=-halt-on-error',
                            ],
                            stdout=outfile,
                            stderr=subprocess.STDOUT,
                        )
                    if returncode:
                        print(
                            'Latex error: check %s'
                            % self.build_dir_join('latex', '__LATEXSTDOUT__')
                        )
                elif '-q' in self.opts:
                    returncode = subprocess.call(
                        [
                            makecmd,
                            'all-pdf',
                            'LATEXOPTS=-halt-on-error',
                            'LATEXMKOPTS=-silent',
                        ],
                    )
                    if returncode:
                        print(
                            'Latex error: check .log file in %s'
                            % self.build_dir_join('latex')
                        )
                else:
                    returncode = subprocess.call([makecmd, 'all-pdf'])
                return returncode
        except OSError:
            print('Error: Failed to run: %s' % makecmd)
            return 1

    def build_latexpdfja(self) -> int:
        if self.run_generic_build('latex') > 0:
            return 1

        # Use $MAKE to determine the make command
        make_fallback = 'make.bat' if sys.platform == 'win32' else 'make'
        makecmd = os.environ.get('MAKE', make_fallback)
        if not makecmd.lower().startswith('make'):
            raise RuntimeError('Invalid $MAKE command: %r' % makecmd)
        try:
            with chdir(self.build_dir_join('latex')):
                return subprocess.call([makecmd, 'all-pdf'])
        except OSError:
            print('Error: Failed to run: %s' % makecmd)
            return 1

    def build_info(self) -> int:
        if self.run_generic_build('texinfo') > 0:
            return 1

        # Use $MAKE to determine the make command
        makecmd = os.environ.get('MAKE', 'make')
        if not makecmd.lower().startswith('make'):
            raise RuntimeError('Invalid $MAKE command: %r' % makecmd)
        try:
            with chdir(self.build_dir_join('texinfo')):
                return subprocess.call([makecmd, 'info'])
        except OSError:
            print('Error: Failed to run: %s' % makecmd)
            return 1

    def build_gettext(self) -> int:
        dtdir = self.build_dir_join('gettext', '.doctrees')
        if self.run_generic_build('gettext', doctreedir=dtdir) > 0:
            return 1
        return 0

    def run_generic_build(self, builder: str, doctreedir: str | None = None) -> int:
        # compatibility with old Makefile
        paper_size = os.getenv('PAPER', '')
        if paper_size in {'a4', 'letter'}:
            self.opts.extend(['-D', f'latex_elements.papersize={paper_size}paper'])
        if doctreedir is None:
            doctreedir = self.build_dir_join('doctrees')

        args = [
            '--builder',
            builder,
            '--doctree-dir',
            doctreedir,
            self.source_dir,
            self.build_dir_join(builder),
        ]
        return build_main(args + self.opts)


def run_make_mode(args: Sequence[str]) -> int:
    if len(args) < 3:
        print(
            'Error: at least 3 arguments (builder, source '
            'dir, build dir) are required.',
            file=sys.stderr,
        )
        return 1

    builder_name = args[0]
    make = Make(source_dir=args[1], build_dir=args[2], opts=args[3:])
    run_method = f'build_{builder_name}'
    if hasattr(make, run_method):
        return getattr(make, run_method)()
    return make.run_generic_build(builder_name)
