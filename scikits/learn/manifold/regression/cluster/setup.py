import numpy
from ConfigParser import ConfigParser

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_standard_file

    config = Configuration('cluster',parent_package,top_path)

    site_cfg = ConfigParser()
    site_cfg.read(get_standard_file('site.cfg'))
    if site_cfg.has_section('scikit-learn') and site_cfg.getboolean('scikit-learn', 'use_boost'):
        # build this extension if enabled in site.cfg
        include_dirs=['../../src', '.', numpy.get_include()]
        config.add_extension('_modified_general_clustering',
                             sources=['ModifiedGeneralClustering.i'],
                             include_dirs=include_dirs
                             )

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
