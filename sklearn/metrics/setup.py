import os
import os.path

import numpy
from numpy.distutils.misc_util import Configuration

from sklearn._build_utils import get_blas_info


def configuration(parent_package="", top_path=None):
    config = Configuration("metrics", parent_package, top_path)

    cblas_libs, blas_info = get_blas_info()
    if os.name == 'posix':
        cblas_libs.append('m')

    config.add_extension("pairwise_fast",
                         sources=["pairwise_fast.c"],
                         include_dirs=[os.path.join('..', 'src', 'cblas'),
                                       numpy.get_include(),
                                       blas_info.pop('include_dirs', [])],
                         libraries=cblas_libs,
                         extra_compile_args=blas_info.pop('extra_compile_args',
                                                          []),
                         **blas_info)

    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration().todict())
