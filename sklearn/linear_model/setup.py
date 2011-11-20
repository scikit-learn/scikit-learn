from os.path import join
import numpy


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    config = Configuration('linear_model', parent_package, top_path)

    # cd fast needs CBLAS
    blas_info = get_info('blas_opt', 0)
    if (not blas_info) or (
        ('NO_ATLAS_INFO', 1) in blas_info.get('define_macros', [])):
        cblas_libs = ['cblas']
        blas_info.pop('libraries', None)
    else:
        cblas_libs = blas_info.pop('libraries', [])

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
         include_dirs=[numpy.get_include()]
         )
    config.add_extension('sgd_fast_sparse',
         sources=['sgd_fast_sparse.c'],
         include_dirs=[numpy.get_include()]
         )

    # add other directories
    config.add_subpackage('tests')
    config.add_subpackage('sparse')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
