from __future__ import division, print_function, absolute_import

import os.path
from os.path import join

from scipy._build_utils import numpy_nodepr_api

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from scipy._build_utils.system_info import get_info
    config = Configuration('optimize',parent_package, top_path)

    include_dirs = [join(os.path.dirname(__file__), '..', '_lib', 'src')]

    minpack_src = [join('minpack','*f')]
    config.add_library('minpack',sources=minpack_src)
    config.add_extension('_minpack',
                         sources=['_minpackmodule.c'],
                         libraries=['minpack'],
                         depends=(["minpack.h","__minpack.h"]
                                  + minpack_src),
                         include_dirs=include_dirs,
                         **numpy_nodepr_api)

    config.add_library('rectangular_lsap',
                       sources='rectangular_lsap/rectangular_lsap.cpp',
                       headers='rectangular_lsap/rectangular_lsap.h')
    config.add_extension('_lsap_module',
                         sources=['_lsap_module.c'],
                         libraries=['rectangular_lsap'],
                         depends=(['rectangular_lsap/rectangular_lsap.cpp',
                                   'rectangular_lsap/rectangular_lsap.h']),
                         include_dirs=include_dirs,
                         **numpy_nodepr_api)

    rootfind_src = [join('Zeros','*.c')]
    rootfind_hdr = [join('Zeros','zeros.h')]
    config.add_library('rootfind',
                       sources=rootfind_src,
                       headers=rootfind_hdr,
                         **numpy_nodepr_api)

    config.add_extension('_zeros',
                         sources=['zeros.c'],
                         libraries=['rootfind'],
                         depends=(rootfind_src + rootfind_hdr),
                         **numpy_nodepr_api)

    lapack = get_info('lapack_opt')
    if 'define_macros' in numpy_nodepr_api:
        if ('define_macros' in lapack) and (lapack['define_macros'] is not None):
            lapack['define_macros'] = (lapack['define_macros'] +
                                       numpy_nodepr_api['define_macros'])
        else:
            lapack['define_macros'] = numpy_nodepr_api['define_macros']
    sources = ['lbfgsb.pyf', 'lbfgsb.f', 'linpack.f', 'timer.f']
    config.add_extension('_lbfgsb',
                         sources=[join('lbfgsb_src',x) for x in sources],
                         **lapack)

    sources = ['moduleTNC.c','tnc.c']
    config.add_extension('moduleTNC',
                         sources=[join('tnc',x) for x in sources],
                         depends=[join('tnc','tnc.h')],
                         **numpy_nodepr_api)

    config.add_extension('_cobyla',
                         sources=[join('cobyla',x) for x in ['cobyla.pyf',
                                                             'cobyla2.f',
                                                             'trstlp.f']],
                         **numpy_nodepr_api)

    sources = ['minpack2.pyf', 'dcsrch.f', 'dcstep.f']
    config.add_extension('minpack2',
                         sources=[join('minpack2',x) for x in sources],
                         **numpy_nodepr_api)

    sources = ['slsqp.pyf', 'slsqp_optmz.f']
    config.add_extension('_slsqp', sources=[join('slsqp', x) for x in sources],
                         **numpy_nodepr_api)

    config.add_extension('_nnls', sources=[join('nnls', x)
                                          for x in ["nnls.f","nnls.pyf"]],
                         **numpy_nodepr_api)

    config.add_extension('_group_columns', sources=['_group_columns.c'],)

    config.add_extension('_bglu_dense', sources=['_bglu_dense.c'])

    config.add_subpackage('_lsq')

    config.add_subpackage('_trlib')

    config.add_subpackage('_trustregion_constr')

    # cython optimize API for zeros functions
    config.add_subpackage('cython_optimize')
    config.add_data_files('cython_optimize.pxd')
    config.add_data_files(os.path.join('cython_optimize', '*.pxd'))
    config.add_extension(
        'cython_optimize._zeros',
        sources=[os.path.join('cython_optimize', '_zeros.c')])

    config.add_subpackage('_shgo_lib')
    config.add_data_dir('_shgo_lib')

    config.add_data_dir('tests')

    # Add license files
    config.add_data_files('lbfgsb_src/README')

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
