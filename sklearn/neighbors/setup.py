import os


def configuration(parent_package='', top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration('neighbors', parent_package, top_path)
    libraries = []
    if os.name == 'posix':
        libraries.append('m')

    config.add_extension('binary_tree',
                         sources=['binary_tree.c'],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    config.add_extension('dist_metrics',
                         sources=['dist_metrics.c'],
                         include_dirs=[numpy.get_include(),
                                       os.path.join(numpy.get_include(),
                                                    'numpy')],
                         libraries=libraries)

    return config
