from os.path import join
import numpy
from ConfigParser import ConfigParser

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_standard_file

    config = Configuration('sparse', parent_package, top_path)
    site_cfg  = ConfigParser()
    site_cfg.read(get_standard_file('site.cfg'))

    config.add_extension('sgd_fast_sparse',
                         sources=[join('src', 'sgd_fast_sparse.c')],
                         include_dirs=[numpy.get_include()]
                         )

    # add other directories
    # config.add_subpackage('tests')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())


