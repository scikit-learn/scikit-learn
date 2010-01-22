from os.path import join

import os.path
import numpy

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('manifold_learning',parent_package,top_path)
    config.add_subpackage('compression')
    config.add_data_dir('examples')
    config.add_subpackage('projection')
    config.add_subpackage('regression')
    config.add_subpackage('stats')
    include_dirs=['src', numpy.get_include()]
    config.add_extension('regression.cluster._modified_general_clustering',
                         sources=["regression/cluster/ModifiedGeneralClustering.i"],
                         depends=["src/matrix/*.h"],
                         include_dirs=include_dirs
                         )
    config.add_extension('compression.cost_function._cost_function',
                         sources=["compression/cost_function/cost_function.cpp"],
                         depends=["src/matrix/*.h"],
                         include_dirs=include_dirs
                         )

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
