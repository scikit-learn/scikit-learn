import os
from os.path import join

from sklearn._build_utils import gen_from_templates


def configuration(parent_package='', top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration('utils', parent_package, top_path)

    libraries = []
    if os.name == 'posix':
        libraries.append('m')

    config.add_extension('sparsefuncs_fast',
                         sources=['sparsefuncs_fast.pyx'],
                         libraries=libraries)

    config.add_extension('_cython_blas',
                         sources=['_cython_blas.pyx'],
                         libraries=libraries)

    config.add_extension('arrayfuncs',
                         sources=['arrayfuncs.pyx'],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    config.add_extension('murmurhash',
                         sources=['murmurhash.pyx', join(
                             'src', 'MurmurHash3.cpp')],
                         include_dirs=['src'])

    config.add_extension('graph_shortest_path',
                         sources=['graph_shortest_path.pyx'],
                         include_dirs=[numpy.get_include()])

    config.add_extension('_fast_dict',
                         sources=['_fast_dict.pyx'],
                         language="c++",
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    config.add_extension('_openmp_helpers',
                         sources=['_openmp_helpers.pyx'],
                         libraries=libraries)

    # generate _seq_dataset from template
    templates = ['sklearn/utils/_seq_dataset.pyx.tp',
                 'sklearn/utils/_seq_dataset.pxd.tp']
    gen_from_templates(templates, top_path)

    config.add_extension('_seq_dataset',
                         sources=['_seq_dataset.pyx'],
                         include_dirs=[numpy.get_include()])

    config.add_extension('_weight_vector',
                         sources=['_weight_vector.pyx'],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    config.add_extension("_random",
                         sources=["_random.pyx"],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    config.add_extension("_logistic_sigmoid",
                         sources=["_logistic_sigmoid.pyx"],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    config.add_subpackage('tests')

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
