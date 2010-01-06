from os.path import join

import os.path
import numpy

def configuration(parent_package='', top_path=None, package_name='regression'):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(package_name,parent_package,top_path)
    config.add_subpackage('*')
    config.add_extension('neighbors._neighbors',
                         sources=["neighbors/neighbors.cpp"],
                         include_dirs=[os.path.dirname(__file__) + '/..'])
    config.add_extension('cluster._modified_general_clustering',
                         sources=["cluster/ModifiedGeneralClustering.i"],
                         include_dirs=[os.path.dirname(__file__) + '/..', numpy.get_include()])
    config.add_data_dir('tests')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='',
                          package_name='regression').todict())
    #setup(configuration=configuration)
