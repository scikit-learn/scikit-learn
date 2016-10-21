import numpy

from sklearn._build_utils import add_cython_extension


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('sparsetools', parent_package, top_path)

    add_cython_extension(top_path,
                         config,
                         '_traversal',
                         sources=['_traversal.c'],
                         include_dirs=[numpy.get_include()],
                         # libraries=libraries
                         )
    add_cython_extension(top_path,
                         config,
                         '_graph_tools',
                         sources=['_graph_tools.c'],
                         include_dirs=[numpy.get_include()],
                         # libraries=libraries
                         )

    config.add_subpackage('tests')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
