""":mod:`sassutils.distutils` --- :mod:`setuptools`/:mod:`distutils` integrat\
ion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
~~~

This module provides extensions (and some magical monkey-patches, sorry)
of the standard :mod:`distutils` and :mod:`setuptools` (now it's named
Distribute) for libsass.

To use this, add ``libsass`` into ``setup_requires`` (not ``install_requires``)
option of the :file:`setup.py` script::

    from setuptools import setup

    setup(
        # ...,
        setup_requires=['libsass >= 0.6.0']
    )

It will adds :class:`build_sass` command to the :file:`setup.py` script:

.. sourcecode:: console

   $ python setup.py build_sass

This commands builds Sass/SCSS files to compiled CSS files of the project
and makes the package archive (made by :class:`~distutils.command.sdist.sdist`,
:class:`~distutils.command.bdist.bdist`, and so on) to include these compiled
CSS files.

To set the directory of Sass/SCSS source files and the directory to
store compiled CSS files, specify ``sass_manifests`` option::

    from setuptools import find_packages, setup

    setup(
        name='YourPackage',
        packages=find_packages(),
        sass_manifests={
            'your.webapp': ('static/sass', 'static/css')
        },
        setup_requires=['libsass >= 0.6.0']
    )

The option should be a mapping of package names to pairs of paths, e.g.::

    {
        'package': ('static/sass', 'static/css'),
        'package.name': ('static/scss', 'static')
    }

The option can also be a mapping of package names to manifest dictionaries::

    {
        'package': {
            'sass_path': 'static/sass',
            'css_path': 'static/css',
            'strip_extension': True,
        },
    }

.. versionadded:: 0.15.0
    Added ``strip_extension`` so ``a.scss`` is compiled to ``a.css`` instead
    of ``a.scss.css``.  This option will default to ``True`` in the future.

.. versionadded:: 0.6.0
   Added ``--output-style``/``-s`` option to :class:`build_sass` command.

"""
import functools
import os.path

import distutils.errors
import distutils.log
import distutils.util
from setuptools import Command
from setuptools.command.sdist import sdist

from .builder import Manifest
from sass import OUTPUT_STYLES

__all__ = 'build_sass', 'validate_manifests'


def validate_manifests(dist, attr, value):
    """Verifies that ``value`` is an expected mapping of package to
    :class:`sassutils.builder.Manifest`.

    """
    try:
        Manifest.normalize_manifests(value)
    except TypeError:
        raise distutils.errors.DistutilsSetupError(
            attr + "must be a mapping object like: {'package.name': "
            "sassutils.distutils.Manifest('sass/path')}, or as shorten form: "
            "{'package.name': ('sass/path', 'css/path'}), not " +
            repr(value),
        )


class build_sass(Command):
    """Builds Sass/SCSS files to CSS files."""

    description = __doc__
    user_options = [
        (
            'output-style=', 's',
            'Coding style of the compiled result.  Choose one of ' +
            ', '.join(OUTPUT_STYLES),
        ),
    ]

    def initialize_options(self):
        self.package_dir = None
        self.output_style = 'nested'

    def finalize_options(self):
        self.package_dir = {}
        if self.distribution.package_dir:
            self.package_dir = {}
            for name, path in self.distribution.package_dir.items():
                self.package_dir[name] = distutils.util.convert_path(path)

    def run(self):
        manifests = self.distribution.sass_manifests
        manifests = Manifest.normalize_manifests(manifests)
        self.distribution.sass_manifests = manifests
        package_data = self.distribution.package_data
        data_files = self.distribution.data_files or []
        for package_name, manifest in manifests.items():
            package_dir = self.get_package_dir(package_name)
            distutils.log.info("building '%s' sass", package_name)
            css_files = manifest.build(
                package_dir,
                output_style=self.output_style,
            )
            map(distutils.log.info, css_files)
            package_data.setdefault(package_name, []).extend(css_files)
            data_files.append(
                (
                    package_dir,
                    [os.path.join(package_dir, f) for f in css_files],
                ),
            )
        self.distribution.package_data = package_data
        self.distribution.data_files = data_files
        self.distribution.has_data_files = lambda: True
        # See the below monkey patch (end of this source code).
        self.distribution.compiled_sass_files = data_files

    def get_package_dir(self, package):
        """Returns the directory, relative to the top of the source
        distribution, where package ``package`` should be found
        (at least according to the :attr:`package_dir` option, if any).

        Copied from :meth:`distutils.command.build_py.get_package_dir()`
        method.

        """
        path = package.split('.')
        if not self.package_dir:
            if path:
                return os.path.join(*path)
            return ''
        tail = []
        while path:
            try:
                pdir = self.package_dir['.'.join(path)]
            except KeyError:
                tail.insert(0, path[-1])
                del path[-1]
            else:
                tail.insert(0, pdir)
                return os.path.join(*tail)
        else:
            pdir = self.package_dir.get('')
            if pdir is not None:
                tail.insert(0, pdir)
            if tail:
                return os.path.join(*tail)
            return ''


# Does monkey-patching the setuptools.command.sdist.sdist.check_readme()
# method to include compiled Sass files as data files.
if not hasattr(sdist, '_wrapped_check_readme'):
    @functools.wraps(sdist.check_readme)
    def check_readme(self):
        try:
            files = self.distribution.compiled_sass_files
        except AttributeError:
            pass
        else:
            for _, css_files in files:
                self.filelist.extend(css_files)
        return self._wrapped_check_readme()
    sdist._wrapped_check_readme = sdist.check_readme
    sdist.check_readme = check_readme
