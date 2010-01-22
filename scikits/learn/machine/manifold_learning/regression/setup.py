from os.path import join

import os.path
import numpy

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('regression',parent_package,top_path)
    config.add_subpackage('cluster')
    config.add_subpackage('neighbors')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
