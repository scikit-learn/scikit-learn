"""
    sphinx.setup_command
    ~~~~~~~~~~~~~~~~~~~~

    Setuptools/distutils commands to assist the building of sphinx
    documentation.

    :author: Sebastian Wiesner
    :contact: basti.wiesner@gmx.net
    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import os
import sys
from distutils.cmd import Command
from distutils.errors import DistutilsExecError
from io import StringIO
from typing import Any, Dict

from sphinx.application import Sphinx
from sphinx.cmd.build import handle_exception
from sphinx.util.console import color_terminal, nocolor
from sphinx.util.docutils import docutils_namespace, patch_docutils
from sphinx.util.osutil import abspath


class BuildDoc(Command):
    """
    Distutils command to build Sphinx documentation.

    The Sphinx build can then be triggered from distutils, and some Sphinx
    options can be set in ``setup.py`` or ``setup.cfg`` instead of Sphinx's
    own configuration file.

    For instance, from `setup.py`::

       # this is only necessary when not using setuptools/distribute
       from sphinx.setup_command import BuildDoc
       cmdclass = {'build_sphinx': BuildDoc}

       name = 'My project'
       version = '1.2'
       release = '1.2.0'
       setup(
           name=name,
           author='Bernard Montgomery',
           version=release,
           cmdclass=cmdclass,
           # these are optional and override conf.py settings
           command_options={
               'build_sphinx': {
                   'project': ('setup.py', name),
                   'version': ('setup.py', version),
                   'release': ('setup.py', release)}},
       )

    Or add this section in ``setup.cfg``::

       [build_sphinx]
       project = 'My project'
       version = 1.2
       release = 1.2.0
    """

    description = 'Build Sphinx documentation'
    user_options = [
        ('fresh-env', 'E', 'discard saved environment'),
        ('all-files', 'a', 'build all files'),
        ('source-dir=', 's', 'Source directory'),
        ('build-dir=', None, 'Build directory'),
        ('config-dir=', 'c', 'Location of the configuration directory'),
        ('builder=', 'b', 'The builder (or builders) to use. Can be a comma- '
         'or space-separated list. Defaults to "html"'),
        ('warning-is-error', 'W', 'Turn warning into errors'),
        ('project=', None, 'The documented project\'s name'),
        ('version=', None, 'The short X.Y version'),
        ('release=', None, 'The full version, including alpha/beta/rc tags'),
        ('today=', None, 'How to format the current date, used as the '
         'replacement for |today|'),
        ('link-index', 'i', 'Link index.html to the master doc'),
        ('copyright', None, 'The copyright string'),
        ('pdb', None, 'Start pdb on exception'),
        ('verbosity', 'v', 'increase verbosity (can be repeated)'),
        ('nitpicky', 'n', 'nit-picky mode, warn about all missing references'),
        ('keep-going', None, 'With -W, keep going when getting warnings'),
    ]
    boolean_options = ['fresh-env', 'all-files', 'warning-is-error',
                       'link-index', 'nitpicky']

    def initialize_options(self) -> None:
        self.fresh_env = self.all_files = False
        self.pdb = False
        self.source_dir: str = None
        self.build_dir: str = None
        self.builder = 'html'
        self.warning_is_error = False
        self.project = ''
        self.version = ''
        self.release = ''
        self.today = ''
        self.config_dir: str = None
        self.link_index = False
        self.copyright = ''
        # Link verbosity to distutils' (which uses 1 by default).
        self.verbosity = self.distribution.verbose - 1  # type: ignore
        self.traceback = False
        self.nitpicky = False
        self.keep_going = False

    def _guess_source_dir(self) -> str:
        for guess in ('doc', 'docs'):
            if not os.path.isdir(guess):
                continue
            for root, _dirnames, filenames in os.walk(guess):
                if 'conf.py' in filenames:
                    return root
        return os.curdir

    def finalize_options(self) -> None:
        self.ensure_string_list('builder')

        if self.source_dir is None:
            self.source_dir = self._guess_source_dir()
            self.announce('Using source directory %s' % self.source_dir)

        self.ensure_dirname('source_dir')

        if self.config_dir is None:
            self.config_dir = self.source_dir

        if self.build_dir is None:
            build = self.get_finalized_command('build')
            self.build_dir = os.path.join(abspath(build.build_base), 'sphinx')  # type: ignore

        self.doctree_dir = os.path.join(self.build_dir, 'doctrees')

        self.builder_target_dirs = [
            (builder, os.path.join(self.build_dir, builder))
            for builder in self.builder]

    def run(self) -> None:
        if not color_terminal():
            nocolor()
        if not self.verbose:  # type: ignore
            status_stream = StringIO()
        else:
            status_stream = sys.stdout  # type: ignore
        confoverrides: Dict[str, Any] = {}
        if self.project:
            confoverrides['project'] = self.project
        if self.version:
            confoverrides['version'] = self.version
        if self.release:
            confoverrides['release'] = self.release
        if self.today:
            confoverrides['today'] = self.today
        if self.copyright:
            confoverrides['copyright'] = self.copyright
        if self.nitpicky:
            confoverrides['nitpicky'] = self.nitpicky

        for builder, builder_target_dir in self.builder_target_dirs:
            app = None

            try:
                confdir = self.config_dir or self.source_dir
                with patch_docutils(confdir), docutils_namespace():
                    app = Sphinx(self.source_dir, self.config_dir,
                                 builder_target_dir, self.doctree_dir,
                                 builder, confoverrides, status_stream,
                                 freshenv=self.fresh_env,
                                 warningiserror=self.warning_is_error,
                                 verbosity=self.verbosity, keep_going=self.keep_going)
                    app.build(force_all=self.all_files)
                    if app.statuscode:
                        raise DistutilsExecError(
                            'caused by %s builder.' % app.builder.name)
            except Exception as exc:
                handle_exception(app, self, exc, sys.stderr)
                if not self.pdb:
                    raise SystemExit(1) from exc

            if not self.link_index:
                continue

            src = app.config.root_doc + app.builder.out_suffix  # type: ignore
            dst = app.builder.get_outfilename('index')  # type: ignore
            os.symlink(src, dst)
