"""Build Environment used for isolation during sdist building
"""

import logging
import os
import sys
from distutils.sysconfig import get_python_lib
from sysconfig import get_paths

from pip._internal.utils.misc import call_subprocess
from pip._internal.utils.temp_dir import TempDirectory
from pip._internal.utils.ui import open_spinner

logger = logging.getLogger(__name__)


class BuildEnvironment(object):
    """Creates and manages an isolated environment to install build deps
    """

    def __init__(self):
        self._temp_dir = TempDirectory(kind="build-env")
        self._temp_dir.create()

    @property
    def path(self):
        return self._temp_dir.path

    def __enter__(self):
        self.save_path = os.environ.get('PATH', None)
        self.save_pythonpath = os.environ.get('PYTHONPATH', None)
        self.save_nousersite = os.environ.get('PYTHONNOUSERSITE', None)

        install_scheme = 'nt' if (os.name == 'nt') else 'posix_prefix'
        install_dirs = get_paths(install_scheme, vars={
            'base': self.path,
            'platbase': self.path,
        })

        scripts = install_dirs['scripts']
        if self.save_path:
            os.environ['PATH'] = scripts + os.pathsep + self.save_path
        else:
            os.environ['PATH'] = scripts + os.pathsep + os.defpath

        # Note: prefer distutils' sysconfig to get the
        # library paths so PyPy is correctly supported.
        purelib = get_python_lib(plat_specific=0, prefix=self.path)
        platlib = get_python_lib(plat_specific=1, prefix=self.path)
        if purelib == platlib:
            lib_dirs = purelib
        else:
            lib_dirs = purelib + os.pathsep + platlib
        if self.save_pythonpath:
            os.environ['PYTHONPATH'] = lib_dirs + os.pathsep + \
                self.save_pythonpath
        else:
            os.environ['PYTHONPATH'] = lib_dirs

        os.environ['PYTHONNOUSERSITE'] = '1'

        return self.path

    def __exit__(self, exc_type, exc_val, exc_tb):
        def restore_var(varname, old_value):
            if old_value is None:
                os.environ.pop(varname, None)
            else:
                os.environ[varname] = old_value

        restore_var('PATH', self.save_path)
        restore_var('PYTHONPATH', self.save_pythonpath)
        restore_var('PYTHONNOUSERSITE', self.save_nousersite)

    def cleanup(self):
        self._temp_dir.cleanup()

    def install_requirements(self, finder, requirements, message):
        args = [
            sys.executable, '-m', 'pip', 'install', '--ignore-installed',
            '--no-user', '--prefix', self.path, '--no-warn-script-location',
        ]
        if logger.getEffectiveLevel() <= logging.DEBUG:
            args.append('-v')
        for format_control in ('no_binary', 'only_binary'):
            formats = getattr(finder.format_control, format_control)
            args.extend(('--' + format_control.replace('_', '-'),
                         ','.join(sorted(formats or {':none:'}))))
        if finder.index_urls:
            args.extend(['-i', finder.index_urls[0]])
            for extra_index in finder.index_urls[1:]:
                args.extend(['--extra-index-url', extra_index])
        else:
            args.append('--no-index')
        for link in finder.find_links:
            args.extend(['--find-links', link])
        for _, host, _ in finder.secure_origins:
            args.extend(['--trusted-host', host])
        if finder.allow_all_prereleases:
            args.append('--pre')
        if finder.process_dependency_links:
            args.append('--process-dependency-links')
        args.append('--')
        args.extend(requirements)
        with open_spinner(message) as spinner:
            call_subprocess(args, show_stdout=False, spinner=spinner)


class NoOpBuildEnvironment(BuildEnvironment):
    """A no-op drop-in replacement for BuildEnvironment
    """

    def __init__(self):
        pass

    def __enter__(self):
        pass

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def cleanup(self):
        pass

    def install_requirements(self, finder, requirements, message):
        raise NotImplementedError()
