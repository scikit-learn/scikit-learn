#! /usr/bin/env python
# Last Change: Sat Jul 21 03:00 PM 2007 J

# Copyright (C) 2007 Cournapeau David <cournape@gmail.com>
#

descr   = """A set of python modules for machine learning and data mining"""

from os.path import join
import os
import sys

DISTNAME            = 'scikits.learn' 
DESCRIPTION         = 'A set of python modules for machine learning and data mining'
LONG_DESCRIPTION    = descr
MAINTAINER          = 'David Cournapeau',
MAINTAINER_EMAIL    = 'david@ar.media.kyoto-u.ac.jp',
URL                 = 'http://projects.scipy.org/scipy/scikits'
LICENSE             = 'new BSD'
DOWNLOAD_URL        = URL

# The following is more or less random copy/paste from numpy.distutils ...
import setuptools
#from distutils.errors import DistutilsError
#from numpy.distutils.system_info import system_info, NotFoundError, dict_append, so_ext
from numpy.distutils.core import setup, Extension

def configuration(parent_package='',top_path=None, package_name=DISTNAME):
    if os.path.exists('MANIFEST'): os.remove('MANIFEST')
    
    pkg_prefix_dir = os.path.join('scikits', 'learn')
    ## Get the version
    #from scikits.pyaudiolab.info import __version__ as pyaudiolab_version

    from numpy.distutils.misc_util import Configuration
    config = Configuration(package_name,parent_package,top_path,
        #version     = pyaudiolab_version,
        maintainer  = MAINTAINER,
        maintainer_email = MAINTAINER_EMAIL,
        description = DESCRIPTION,
        license = LICENSE,
        url = URL, 
        download_url = DOWNLOAD_URL,
        long_description = LONG_DESCRIPTION)
    config.add_subpackage('scikits.learn.machine')

    return config

if __name__ == "__main__":
    setup(configuration = configuration,
        install_requires='numpy', # can also add version specifiers      
        namespace_packages=['scikits'],
        packages=setuptools.find_packages(),
        include_package_data = True,
        #package_data = {'scikits.pyaudiolab': data_files}, 
        test_suite="tester", # for python setup.py test
        zip_safe=True, # the package can run out of an .egg file
        #FIXME url, download_url, ext_modules
        classifiers = 
            [ 'Development Status :: 4 - Beta',
              'Environment :: Console',
              'Intended Audience :: Developers',
              'Intended Audience :: Science/Research',
              'License :: OSI Approved :: BSD License',
              'Topic :: Scientific/Engineering']
    )
