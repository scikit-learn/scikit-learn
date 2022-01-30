"""
    sphinx.cmd.make_mode
    ~~~~~~~~~~~~~~~~~~~~

    sphinx-build -M command-line handling.

    This replaces the old, platform-dependent and once-generated content
    of Makefile / make.bat.

    This is in its own module so that importing it is fast.  It should not
    import the main Sphinx modules (like sphinx.applications, sphinx.builders).

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import os
import subprocess
import sys
from os import path
from typing import List

import sphinx
from sphinx.cmd.build import build_main
from sphinx.util.console import blue, bold, color_terminal, nocolor  # type: ignore
from sphinx.util.osutil import cd, rmtree

BUILDERS = [
    ("",      "html",        "to make standalone HTML files"),
    ("",      "dirhtml",     "to make HTML files named index.html in directories"),
    ("",      "singlehtml",  "to make a single large HTML file"),
    ("",      "pickle",      "to make pickle files"),
    ("",      "json",        "to make JSON files"),
    ("",      "htmlhelp",    "to make HTML files and an HTML help project"),
    ("",      "qthelp",      "to make HTML files and a qthelp project"),
    ("",      "devhelp",     "to make HTML files and a Devhelp project"),
    ("",      "epub",        "to make an epub"),
    ("",      "latex",       "to make LaTeX files, you can set PAPER=a4 or PAPER=letter"),
    ("posix", "latexpdf",    "to make LaTeX and PDF files (default pdflatex)"),
    ("posix", "latexpdfja",  "to make LaTeX files and run them through platex/dvipdfmx"),
    ("",      "text",        "to make text files"),
    ("",      "man",         "to make manual pages"),
    ("",      "texinfo",     "to make Texinfo files"),
    ("posix", "info",        "to make Texinfo files and run them through makeinfo"),
    ("",      "gettext",     "to make PO message catalogs"),
    ("",      "changes",     "to make an overview of all changed/added/deprecated items"),
    ("",      "xml",         "to make Docutils-native XML files"),
    ("",      "pseudoxml",   "to make pseudoxml-XML files for display purposes"),
    ("",      "linkcheck",   "to check all external links for integrity"),
    ("",      "doctest",     "to run all doctests embedded in the documentation "
                             "(if enabled)"),
    ("",      "coverage",    "to run coverage check of the documentation (if enabled)"),
    ("",      "clean",       "to remove everything in the build directory"),
]


class Make:
    def __init__(self, srcdir: str, builddir: str, opts: List[str]) -> None:
        self.srcdir = srcdir
        self.builddir = builddir
        self.opts = opts
        self.makecmd = os.environ.get('MAKE', 'make')  # refer $MAKE to determine make command

    def builddir_join(self, *comps: str) -> str:
        return path.join(self.builddir, *comps)

    def build_clean(self) -> int:
        srcdir = path.abspath(self.srcdir)
        builddir = path.abspath(self.builddir)
        if not path.exists(self.builddir):
            return 0
        elif not path.isdir(self.builddir):
            print("Error: %r is not a directory!" % self.builddir)
            return 1
        elif srcdir == builddir:
            print("Error: %r is same as source directory!" % self.builddir)
            return 1
        elif path.commonpath([srcdir, builddir]) == builddir:
            print("Error: %r directory contains source directory!" % self.builddir)
            return 1
        print("Removing everything under %r..." % self.builddir)
        for item in os.listdir(self.builddir):
            rmtree(self.builddir_join(item))
        return 0

    def build_help(self) -> None:
        if not color_terminal():
            nocolor()

        print(bold("Sphinx v%s" % sphinx.__display_version__))
        print("Please use `make %s' where %s is one of" % ((blue('target'),) * 2))
        for osname, bname, description in BUILDERS:
            if not osname or os.name == osname:
                print('  %s  %s' % (blue(bname.ljust(10)), description))

    def build_latexpdf(self) -> int:
        if self.run_generic_build('latex') > 0:
            return 1

        if sys.platform == 'win32':
            makecmd = os.environ.get('MAKE', 'make.bat')
        else:
            makecmd = self.makecmd
        try:
            with cd(self.builddir_join('latex')):
                return subprocess.call([makecmd, 'all-pdf'])
        except OSError:
            print('Error: Failed to run: %s' % makecmd)
            return 1

    def build_latexpdfja(self) -> int:
        if self.run_generic_build('latex') > 0:
            return 1

        if sys.platform == 'win32':
            makecmd = os.environ.get('MAKE', 'make.bat')
        else:
            makecmd = self.makecmd
        try:
            with cd(self.builddir_join('latex')):
                return subprocess.call([makecmd, 'all-pdf'])
        except OSError:
            print('Error: Failed to run: %s' % makecmd)
            return 1

    def build_info(self) -> int:
        if self.run_generic_build('texinfo') > 0:
            return 1
        try:
            with cd(self.builddir_join('texinfo')):
                return subprocess.call([self.makecmd, 'info'])
        except OSError:
            print('Error: Failed to run: %s' % self.makecmd)
            return 1

    def build_gettext(self) -> int:
        dtdir = self.builddir_join('gettext', '.doctrees')
        if self.run_generic_build('gettext', doctreedir=dtdir) > 0:
            return 1
        return 0

    def run_generic_build(self, builder: str, doctreedir: str = None) -> int:
        # compatibility with old Makefile
        papersize = os.getenv('PAPER', '')
        opts = self.opts
        if papersize in ('a4', 'letter'):
            opts.extend(['-D', 'latex_elements.papersize=' + papersize + 'paper'])
        if doctreedir is None:
            doctreedir = self.builddir_join('doctrees')

        args = ['-b', builder,
                '-d', doctreedir,
                self.srcdir,
                self.builddir_join(builder)]
        return build_main(args + opts)


def run_make_mode(args: List[str]) -> int:
    if len(args) < 3:
        print('Error: at least 3 arguments (builder, source '
              'dir, build dir) are required.', file=sys.stderr)
        return 1
    make = Make(args[1], args[2], args[3:])
    run_method = 'build_' + args[0]
    if hasattr(make, run_method):
        return getattr(make, run_method)()
    return make.run_generic_build(args[0])
