# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr> 
# License: BSD Style.

# $Id$

import numpy
from os.path import join

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('glm',parent_package,top_path)
    config.add_extension('lasso_cd_fast',
                         sources=['src/lasso_cd_fast.c'],
                         # libraries=['m'],
                         include_dirs=[numpy.get_include()])
    config.add_extension('enet_cd_fast',
                         sources=['src/enet_cd_fast.c'],
                         include_dirs=[numpy.get_include()])
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
