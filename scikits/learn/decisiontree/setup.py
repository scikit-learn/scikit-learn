from os.path import join
import sys
import numpy

if sys.version_info[0] < 3:
    from ConfigParser import ConfigParser
else:
    from configparser import ConfigParser

import warnings

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info, get_standard_file, \
         BlasNotFoundError

    config = Configuration('decisiontree', parent_package, top_path)

    config.add_subpackage('tests')

    libdecisiontree_sources = ['libdecisiontree.c']
    libdecisiontree_depends = [join('src', 'Histogram.h'),
                               join('src', 'Node.h'),
                               join('src', 'Node.cpp'),
                               join('src', 'Object.h')]

    config.add_extension('libdecisiontree',
                         sources = libdecisiontree_sources,
                         include_dirs = [numpy.get_include(), 'src'],
                         depends = libdecisiontree_depends
                         )

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())

