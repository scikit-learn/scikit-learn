from os.path import join
import numpy


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('sparse', parent_package, top_path)

    config.add_extension('cd_fast_sparse',
                         sources=[join('src', 'cd_fast_sparse.c')],
                         include_dirs=[numpy.get_include()]
                         )

    # add other directories
    config.add_subpackage('tests')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
