from __future__ import division, print_function, absolute_import

import sys

if sys.version_info[0] >= 3:
    DEFINE_MACROS = [("SCIPY_PY3K", None)]
else:
    DEFINE_MACROS = []


def configuration(parent_package='', top_path=None):
    from scipy._build_utils.system_info import get_info
    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs
    config = Configuration('cluster', parent_package, top_path)

    blas_opt = get_info('lapack_opt')

    config.add_data_dir('tests')

    config.add_extension('_vq',
        sources=[('_vq.c')],
        include_dirs=[get_numpy_include_dirs()],
        extra_info=blas_opt)

    config.add_extension('_hierarchy',
        sources=[('_hierarchy.c')],
        include_dirs=[get_numpy_include_dirs()])

    config.add_extension('_optimal_leaf_ordering',
        sources=[('_optimal_leaf_ordering.c')],
        include_dirs=[get_numpy_include_dirs()])

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
