import os


def configuration(parent_package='', top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration('neighbors', parent_package, top_path)
    libraries = []
    if os.name == 'posix':
        libraries.append('m')

    config.add_extension('ball_tree',
                         sources=['ball_tree.pyx'],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    config.add_extension('kd_tree',
                         sources=['kd_tree.pyx'],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    config.add_extension('dist_metrics',
                         sources=['dist_metrics.pyx'],
                         include_dirs=[numpy.get_include(),
                                       os.path.join(numpy.get_include(),
                                                    'numpy')],
                         libraries=libraries)

    config.add_extension('typedefs',
                         sources=['typedefs.pyx'],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)
    config.add_extension("quad_tree",
                         sources=["quad_tree.pyx"],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    config.add_subpackage('tests')

    return config
