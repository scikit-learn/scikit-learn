import os
from os.path import join

import numpy

from sklearn._build_utils import get_blas_info

from Cython import Tempita


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('linear_model', parent_package, top_path)

    cblas_libs, blas_info = get_blas_info()

    if os.name == 'posix':
        cblas_libs.append('m')

    config.add_extension('cd_fast', sources=['cd_fast.pyx'],
                         libraries=cblas_libs,
                         include_dirs=[join('..', 'src', 'cblas'),
                                       numpy.get_include(),
                                       blas_info.pop('include_dirs', [])],
                         extra_compile_args=blas_info.pop('extra_compile_args',
                                                          []), **blas_info)

    config.add_extension('sgd_fast',
                         sources=['sgd_fast.pyx'],
                         include_dirs=[join('..', 'src', 'cblas'),
                                       numpy.get_include(),
                                       blas_info.pop('include_dirs', [])],
                         libraries=cblas_libs,
                         extra_compile_args=blas_info.pop('extra_compile_args',
                                                          []),
                         **blas_info)

    # generate sag_fast from template
    sag_cython_file = 'sklearn/linear_model/sag_fast.pyx.tp'
    sag_file = sag_cython_file[:-3]

    if not (os.path.exists(sag_file) and
            os.stat(sag_cython_file).st_mtime < os.stat(sag_file).st_mtime):

        with open(sag_cython_file, "r") as f:
            tmpl = f.read()
        tmpl_ = Tempita.sub(tmpl)

        with open(sag_file, "w") as f:
            f.write(tmpl_)

    config.add_extension('sag_fast',
                         sources=['sag_fast.pyx'],
                         include_dirs=numpy.get_include())

    # add other directories
    config.add_subpackage('tests')

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
