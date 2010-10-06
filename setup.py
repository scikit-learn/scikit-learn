#! /usr/bin/env python
# Last Change: Sat Jul 21 09:00 PM 2007 J

# Copyright (C) 2007-2009 Cournapeau David <cournape@gmail.com>
#               2010 Fabian Pedregosa <fabian.pedregosa@inria.fr>

descr   = """A set of python modules for machine learning and data mining"""

import os

DISTNAME            = 'scikits.learn' 
DESCRIPTION         = 'A set of python modules for machine learning and data mining'
LONG_DESCRIPTION    = descr
MAINTAINER          = 'Fabian Pedregosa'
MAINTAINER_EMAIL    = 'fabian.pedregosa@inria.fr'
URL                 = 'http://scikit-learn.sourceforge.net'
LICENSE             = 'new BSD'
DOWNLOAD_URL        = 'http://sourceforge.net/projects/scikit-learn/files/'
VERSION             = '0.5-rc'

import setuptools # we are using a setuptools namespace
from numpy.distutils.core import setup

def configuration(parent_package='',top_path=None):
    if os.path.exists('MANIFEST'): os.remove('MANIFEST')

    from numpy.distutils.misc_util import Configuration
    config = Configuration(None, parent_package,top_path,
        namespace_packages=['scikits'])

    config.add_subpackage('scikits.learn')
    config.add_data_files('scikits/__init__.py')

    return config

if __name__ == "__main__":
    setup(configuration = configuration,
        name = DISTNAME,
        maintainer  = MAINTAINER,
        include_package_data = True,
        maintainer_email = MAINTAINER_EMAIL,
        description = DESCRIPTION,
        license = LICENSE,
        url = URL,
        version = VERSION,
        download_url = DOWNLOAD_URL,
        long_description = LONG_DESCRIPTION,
        zip_safe=False, # the package can run out of an .egg file
        classifiers = 
            ['Intended Audience :: Science/Research',
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
             ]
    )
