#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Cournapeau David <cournape@gmail.com>
#               2010 Fabian Pedregosa <fabian.pedregosa@inria.fr>

descr = """A set of python modules for machine learning and data mining"""

import sys
import os
import shutil

if sys.version_info[0] < 3:
    import __builtin__ as builtins
else:
    import builtins

# This is a bit (!) hackish: we are setting a global variable so that the main
# sklearn __init__ can detect if it is being loaded by the setup routine, to
# avoid attempting to load components that aren't built yet.
builtins.__SKLEARN_SETUP__ = True

DISTNAME = 'scikit-learn'
DESCRIPTION = 'A set of python modules for machine learning and data mining'
LONG_DESCRIPTION = open('README.rst').read()
MAINTAINER = 'Andreas Mueller'
MAINTAINER_EMAIL = 'amueller@ais.uni-bonn.de'
URL = 'http://scikit-learn.org'
LICENSE = 'new BSD'
DOWNLOAD_URL = 'http://sourceforge.net/projects/scikit-learn/files/'

# We can actually import a restricted version of sklearn that
# does not need the compiled code
import sklearn
VERSION = sklearn.__version__

from numpy.distutils.core import setup

###############################################################################
# Optional setuptools features

# For some commands, use setuptools
if len(set(('develop', 'release', 'bdist_egg', 'bdist_rpm',
           'bdist_wininst', 'install_egg_info', 'build_sphinx',
           'egg_info', 'easy_install', 'upload',
            )).intersection(sys.argv)) > 0:
    extra_setuptools_args = dict(
            zip_safe=False, # the package can run out of an .egg file
            include_package_data=True,
        )
else:
    extra_setuptools_args = dict()

def configuration(parent_package='', top_path=None):
    if os.path.exists('MANIFEST'):
        os.remove('MANIFEST')

    from numpy.distutils.misc_util import Configuration
    config = Configuration(None, parent_package, top_path)

    # Avoid non-useful msg:
    # "Ignoring attempt to set 'name' (from ... "
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)

    config.add_subpackage('sklearn')

    return config


if __name__ == "__main__":

    old_path = os.getcwd()
    local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    # python 3 compatibility stuff.
    # Simplified version of scipy strategy: copy files into
    # build/py3k, and patch them using lib2to3.
    if sys.version_info[0] == 3:
        try:
            import lib2to3cache
        except ImportError:
            pass
        local_path = os.path.join(local_path, 'build', 'py3k')
        if os.path.exists(local_path):
            shutil.rmtree(local_path)
        print("Copying source tree into build/py3k for 2to3 transformation"
              "...")
        shutil.copytree(os.path.join(old_path, 'sklearn'),
                        os.path.join(local_path, 'sklearn'))
        import lib2to3.main
        from io import StringIO
        print("Converting to Python3 via 2to3...")
        _old_stdout = sys.stdout
        try:
            sys.stdout = StringIO()  # supress noisy output
            res = lib2to3.main.main("lib2to3.fixes", ['-x', 'import', '-w', local_path])
        finally:
            sys.stdout = _old_stdout

        if res != 0:
            raise Exception('2to3 failed, exiting ...')

    os.chdir(local_path)
    sys.path.insert(0, local_path)

    setup(configuration=configuration,
          name=DISTNAME,
          maintainer=MAINTAINER,
          maintainer_email=MAINTAINER_EMAIL,
          description=DESCRIPTION,
          license=LICENSE,
          url=URL,
          version=VERSION,
          download_url=DOWNLOAD_URL,
          long_description=LONG_DESCRIPTION,
          classifiers=[
              'Intended Audience :: Science/Research',
              'Intended Audience :: Developers',
              'License :: OSI Approved',
              'Programming Language :: C',
              'Programming Language :: Python',
              'Topic :: Software Development',
              'Topic :: Scientific/Engineering',
              'Operating System :: Microsoft :: Windows',
              'Operating System :: POSIX',
              'Operating System :: Unix',
              'Operating System :: MacOS'
             ],
      **extra_setuptools_args)
