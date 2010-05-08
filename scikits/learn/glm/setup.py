# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr> 
# License: BSD Style.

# $Id$

import numpy
from os.path import join

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('glm',parent_package,top_path)
    config.add_extension('cd_fast',
                         sources=['src/cd_fast.c'],
                         # libraries=['m'],
                         include_dirs=[numpy.get_include()])
    config.add_data_dir('test')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
