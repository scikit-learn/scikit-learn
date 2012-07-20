import os
from os.path import join

import numpy
import imp


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('linear_model', parent_package, top_path)

    # get the blas finding routine from sklearn/setup.py
    parent_path = os.path.join(os.path.dirname(__file__), "../")
    sm = imp.find_module("setup", [parent_path])
    setup_mod = imp.load_module("setup_mod", *sm)
    sm[0].close()

    cblas_libs, blas_info = setup_mod.get_blas_info()

    libraries = []
    if os.name == 'posix':
        cblas_libs.append('m')
        libraries.append('m')

    config.add_extension('cd_fast',
         sources=['cd_fast.c'],
         libraries=cblas_libs,
         include_dirs=[join('..', 'src', 'cblas'),
                       numpy.get_include(),
                       blas_info.pop('include_dirs', [])],
         extra_compile_args=blas_info.pop('extra_compile_args', []),
         **blas_info
         )

    config.add_extension('sgd_fast',
         sources=['sgd_fast.c'],
         include_dirs=[numpy.get_include()],
         libraries=libraries,
         )

    # add other directories
    config.add_subpackage('tests')
    config.add_subpackage('sparse')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
