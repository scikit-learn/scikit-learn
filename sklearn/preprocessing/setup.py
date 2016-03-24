import os

def configuration(parent_package='', top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration('preprocessing', parent_package, top_path)
    #config.add_subpackage('discretizer')

    libraries = []
    if os.name == 'posix':
        libraries.append('m')
        #cblas_libs.append('m')

    config.add_extension('_discretization', sources=['_discretization.c'],
                         libraries=libraries)

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())

