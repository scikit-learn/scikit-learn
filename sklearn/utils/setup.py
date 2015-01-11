import os
from os.path import join

from sklearn._build_utils import get_blas_info


def configuration(parent_package='', top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration('utils', parent_package, top_path)
    config.add_subpackage('sparsetools')

    cblas_libs, blas_info = get_blas_info()
    cblas_compile_args = blas_info.pop('extra_compile_args', [])
    cblas_includes = [join('..', 'src', 'cblas'),
                      numpy.get_include(),
                      blas_info.pop('include_dirs', [])]

    libraries = []
    if os.name == 'posix':
        libraries.append('m')
        cblas_libs.append('m')

    config.add_extension('sparsefuncs_fast', sources=['sparsefuncs_fast.c'],
                         libraries=libraries)

    config.add_extension('arrayfuncs',
                         sources=['arrayfuncs.c'],
                         depends=[join('src', 'cholesky_delete.h')],
                         libraries=cblas_libs,
                         include_dirs=cblas_includes,
                         extra_compile_args=cblas_compile_args,
                         **blas_info
                         )

    config.add_extension(
        'murmurhash',
        sources=['murmurhash.c', join('src', 'MurmurHash3.cpp')],
        include_dirs=['src'])

    config.add_extension('lgamma',
                         sources=['lgamma.c', join('src', 'gamma.c')],
                         include_dirs=['src'],
                         libraries=libraries)

    config.add_extension('graph_shortest_path',
                         sources=['graph_shortest_path.c'],
                         include_dirs=[numpy.get_include()])

    config.add_extension('fast_dict',
                         sources=['fast_dict.cpp'],
                         language="c++",
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    config.add_extension('seq_dataset',
                         sources=['seq_dataset.c'],
                         include_dirs=[numpy.get_include()])

    config.add_extension('weight_vector',
                         sources=['weight_vector.c'],
                         include_dirs=cblas_includes,
                         libraries=cblas_libs,
                         **blas_info)

    config.add_extension("_random",
                         sources=["_random.c"],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    config.add_extension("_logistic_sigmoid",
                         sources=["_logistic_sigmoid.c"],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
