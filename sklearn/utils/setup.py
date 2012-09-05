import os
from os.path import join

from sklearn._build_utils import get_blas_info


def configuration(parent_package='', top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration('utils', parent_package, top_path)
    config.add_subpackage('sparsetools')

    cblas_libs, blas_info = get_blas_info()

    libraries = []
    if os.name == 'posix':
        libraries.append('m')
        cblas_libs.append('m')

    config.add_extension('arraybuilder',
         sources=['arraybuilder.c'])

    config.add_extension('sparsefuncs',
         sources=['sparsefuncs.c'],
         libraries=libraries)

    config.add_extension('arrayfuncs',
         sources=['arrayfuncs.c'],
         depends=[join('src', 'cholesky_delete.c')],
         libraries=cblas_libs,
         include_dirs=[join('..', 'src', 'cblas'),
                       numpy.get_include(),
                       blas_info.pop('include_dirs', [])],
         extra_compile_args=blas_info.pop('extra_compile_args', []),
         **blas_info
         )

    config.add_extension(
        'murmurhash',
        sources=['murmurhash.c', join('src', 'MurmurHash3.cpp')],
        include_dirs=['src'])

    config.add_extension('graph_shortest_path',
         sources=['graph_shortest_path.c'],
         include_dirs=[numpy.get_include()])

    config.add_extension('seq_dataset',
         sources=['seq_dataset.c'],
         include_dirs=[numpy.get_include()])

    config.add_extension('weight_vector',
         sources=['weight_vector.c'],
         include_dirs=[numpy.get_include()],
         libraries=libraries)

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
