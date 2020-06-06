from __future__ import division, print_function, absolute_import

import os
from os.path import join

from scipy._build_utils import numpy_nodepr_api


def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from scipy._build_utils.system_info import get_info
    config = Configuration('integrate', parent_package, top_path)

    # Get a local copy of lapack_opt_info
    lapack_opt = dict(get_info('lapack_opt',notfound_action=2))
    # Pop off the libraries list so it can be combined with
    # additional required libraries
    lapack_libs = lapack_opt.pop('libraries', [])

    mach_src = [join('mach','*.f')]
    quadpack_src = [join('quadpack', '*.f')]
    lsoda_src = [join('odepack', fn) for fn in [
        'blkdta000.f', 'bnorm.f', 'cfode.f',
        'ewset.f', 'fnorm.f', 'intdy.f',
        'lsoda.f', 'prja.f', 'solsy.f', 'srcma.f',
        'stoda.f', 'vmnorm.f', 'xerrwv.f', 'xsetf.f',
        'xsetun.f']]
    vode_src = [join('odepack', 'vode.f'), join('odepack', 'zvode.f')]
    dop_src = [join('dop','*.f')]
    quadpack_test_src = [join('tests','_test_multivariate.c')]
    odeint_banded_test_src = [join('tests', 'banded5x5.f')]

    config.add_library('mach', sources=mach_src,
                       config_fc={'noopt':(__file__,1)})
    config.add_library('quadpack', sources=quadpack_src)
    config.add_library('lsoda', sources=lsoda_src)
    config.add_library('vode', sources=vode_src)
    config.add_library('dop', sources=dop_src)

    # Extensions
    # quadpack:
    include_dirs = [join(os.path.dirname(__file__), '..', '_lib', 'src')]
    if 'include_dirs' in lapack_opt:
        lapack_opt = dict(lapack_opt)
        include_dirs.extend(lapack_opt.pop('include_dirs'))

    config.add_extension('_quadpack',
                         sources=['_quadpackmodule.c'],
                         libraries=['quadpack', 'mach'] + lapack_libs,
                         depends=(['__quadpack.h']
                                  + quadpack_src + mach_src),
                         include_dirs=include_dirs,
                         **lapack_opt)

    # odepack/lsoda-odeint
    odepack_opts = lapack_opt.copy()
    odepack_opts.update(numpy_nodepr_api)
    config.add_extension('_odepack',
                         sources=['_odepackmodule.c'],
                         libraries=['lsoda', 'mach'] + lapack_libs,
                         depends=(lsoda_src + mach_src),
                         **odepack_opts)

    # vode
    config.add_extension('vode',
                         sources=['vode.pyf'],
                         libraries=['vode'] + lapack_libs,
                         depends=vode_src,
                         **lapack_opt)

    # lsoda
    config.add_extension('lsoda',
                         sources=['lsoda.pyf'],
                         libraries=['lsoda', 'mach'] + lapack_libs,
                         depends=(lsoda_src + mach_src),
                         **lapack_opt)

    # dop
    config.add_extension('_dop',
                         sources=['dop.pyf'],
                         libraries=['dop'],
                         depends=dop_src)

    config.add_extension('_test_multivariate',
                         sources=quadpack_test_src)

    # Fortran+f2py extension module for testing odeint.
    config.add_extension('_test_odeint_banded',
                         sources=odeint_banded_test_src,
                         libraries=['lsoda', 'mach'] + lapack_libs,
                         depends=(lsoda_src + mach_src),
                         **lapack_opt)

    config.add_subpackage('_ivp')

    config.add_data_dir('tests')
    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
