import os


def configuration(parent_package='', top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration('neural_networks', parent_package, top_path)
    libraries = []
    if os.name == 'posix':
        libraries.append('m')

    config.add_extension('backprop_sgd',
                         sources=['backprop_sgd.c'],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    return config
