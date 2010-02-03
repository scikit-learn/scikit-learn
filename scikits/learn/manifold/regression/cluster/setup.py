import numpy

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('cluster',parent_package,top_path)

    include_dirs=['../../src', '.', numpy.get_include()]

    config.add_extension('_modified_general_clustering',
                         sources=["ModifiedGeneralClustering.i"],
                         include_dirs=include_dirs
                         )
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
