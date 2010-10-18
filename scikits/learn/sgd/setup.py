from os.path import join
import warnings
import numpy
import sys
from ConfigParser import ConfigParser

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info, get_standard_file, BlasNotFoundError
    config = Configuration('sgd', parent_package, top_path)

    site_cfg  = ConfigParser()
    site_cfg.read(get_standard_file('site.cfg'))
    
    

    # add other directories
    config.add_subpackage('tests')
    config.add_subpackage('sparse')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
