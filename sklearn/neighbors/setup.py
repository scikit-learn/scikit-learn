import os

from sklearn._build_utils import add_cython_extension


def configuration(parent_package='', top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration('neighbors', parent_package, top_path)
    libraries = []
    if os.name == 'posix':
        libraries.append('m')

    add_cython_extension(top_path,
                         config,
                         'ball_tree',
                         sources=['ball_tree.c'],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    add_cython_extension(top_path,
                         config,
                         'kd_tree',
                         sources=['kd_tree.c'],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    add_cython_extension(top_path,
                         config,
                         'dist_metrics',
                         sources=['dist_metrics.c'],
                         include_dirs=[numpy.get_include(),
                                       os.path.join(numpy.get_include(),
                                                    'numpy')],
                         libraries=libraries)

    add_cython_extension(top_path,
                         config,
                         'typedefs',
                         sources=['typedefs.c'],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    config.add_subpackage('tests')

    return config
