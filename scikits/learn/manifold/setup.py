# Author: Jake Vanderplas <vanderplas@astro.washington.edu>
# License: BSD, (C) 2011

import numpy
from os.path import join

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('manifold', parent_package, top_path)
    config.add_extension('shortest_path',
                         sources=['shortest_path.c'],
                         include_dirs=[numpy.get_include()])    
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
