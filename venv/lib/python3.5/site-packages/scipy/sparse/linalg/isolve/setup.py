from __future__ import division, print_function, absolute_import

from os.path import join


def configuration(parent_package='',top_path=None):
    from numpy.distutils.system_info import get_info, NotFoundError
    from numpy.distutils.misc_util import Configuration
    from scipy._build_utils import get_g77_abi_wrappers

    config = Configuration('isolve',parent_package,top_path)

    lapack_opt = get_info('lapack_opt')

    if not lapack_opt:
        raise NotFoundError('no lapack/blas resources found')

    # iterative methods
    methods = ['BiCGREVCOM.f.src',
               'BiCGSTABREVCOM.f.src',
               'CGREVCOM.f.src',
               'CGSREVCOM.f.src',
#               'ChebyREVCOM.f.src',
               'GMRESREVCOM.f.src',
#               'JacobiREVCOM.f.src',
               'QMRREVCOM.f.src',
#               'SORREVCOM.f.src'
               ]

    Util = ['getbreak.f.src']
    sources = Util + methods + ['_iterative.pyf.src']
    sources = [join('iterative', x) for x in sources]
    sources += get_g77_abi_wrappers(lapack_opt)

    config.add_extension('_iterative',
                         sources=sources,
                         extra_info=lapack_opt)

    config.add_data_dir('tests')

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup

    setup(**configuration(top_path='').todict())
