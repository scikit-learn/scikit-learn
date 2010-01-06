from os.path import join

import os.path
import numpy

def configuration(parent_package='', top_path=None, package_name='manifold_learning'):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(package_name,parent_package,top_path)
    config.add_subpackage('*')
    config.add_library('compression.cost_function._cost_function',
                         sources=["compression/cost_function/cost_function.cpp"],
                         include_dirs=[os.path.dirname(__file__)],)
    config.add_library('regression.neighboors._neighboors',
                         sources=["regression/neighboors/neighboors.cpp"],
                         include_dirs=[os.path.dirname(__file__)])
    config.add_extension('regression.cluster._modified_general_clustering',
                         sources=["regression/cluster/ModifiedGeneralClustering.i"],
                         include_dirs=[os.path.dirname(__file__), numpy.get_include()])
    config.add_data_dir('matrix')
    config.add_data_dir('tests')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='',
                          package_name='manifold_learning').todict())
    #setup(configuration=configuration)
