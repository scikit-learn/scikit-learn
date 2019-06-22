from __future__ import division, print_function, absolute_import

import os
import sys
from os.path import join, dirname
from distutils.sysconfig import get_python_inc
import subprocess
import numpy
from numpy.distutils.misc_util import get_numpy_include_dirs

try:
    from numpy.distutils.misc_util import get_info
except ImportError:
    raise ValueError("numpy >= 1.4 is required (detected %s from %s)" %
                     (numpy.__version__, numpy.__file__))


def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from scipy._build_utils.system_info import get_info as get_system_info

    config = Configuration('special', parent_package, top_path)

    define_macros = []
    if sys.platform == 'win32':
        # define_macros.append(('NOINFINITIES',None))
        # define_macros.append(('NONANS',None))
        define_macros.append(('_USE_MATH_DEFINES',None))

    curdir = os.path.abspath(os.path.dirname(__file__))
    inc_dirs = [get_python_inc(), os.path.join(curdir, "c_misc")]
    if inc_dirs[0] != get_python_inc(plat_specific=1):
        inc_dirs.append(get_python_inc(plat_specific=1))
    inc_dirs.insert(0, get_numpy_include_dirs())
    inc_dirs.append(join(dirname(dirname(__file__)), '_lib'))

    # C libraries
    c_misc_src = [join('c_misc','*.c')]
    c_misc_hdr = [join('c_misc','*.h')]
    cephes_src = [join('cephes','*.c')]
    cephes_hdr = [join('cephes', '*.h')]
    config.add_library('sc_c_misc',sources=c_misc_src,
                       include_dirs=[curdir] + inc_dirs,
                       depends=(cephes_hdr + cephes_src
                                + c_misc_hdr + cephes_hdr
                                + ['*.h']),
                       macros=define_macros)
    config.add_library('sc_cephes',sources=cephes_src,
                       include_dirs=[curdir] + inc_dirs,
                       depends=(cephes_hdr + ['*.h']),
                       macros=define_macros)

    # Fortran/C++ libraries
    mach_src = [join('mach','*.f')]
    amos_src = [join('amos','*.f')]
    cdf_src = [join('cdflib','*.f')]
    specfun_src = [join('specfun','*.f')]
    config.add_library('sc_mach',sources=mach_src,
                       config_fc={'noopt':(__file__,1)})
    config.add_library('sc_amos',sources=amos_src)
    config.add_library('sc_cdf',sources=cdf_src)
    config.add_library('sc_specfun',sources=specfun_src)

    # Extension specfun
    config.add_extension('specfun',
                         sources=['specfun.pyf'],
                         f2py_options=['--no-wrap-functions'],
                         depends=specfun_src,
                         define_macros=[],
                         libraries=['sc_specfun'])

    # Extension _ufuncs
    headers = ['*.h', join('c_misc', '*.h'), join('cephes', '*.h')]
    ufuncs_src = ['_ufuncs.c', 'sf_error.c', '_logit.c.src',
                  "amos_wrappers.c", "cdf_wrappers.c", "specfun_wrappers.c"]
    ufuncs_dep = (headers + ufuncs_src + amos_src + c_misc_src + cephes_src
                  + mach_src + cdf_src + specfun_src)
    cfg = dict(get_system_info('lapack_opt'))
    cfg.setdefault('include_dirs', []).extend([curdir] + inc_dirs + [numpy.get_include()])
    cfg.setdefault('libraries', []).extend(['sc_amos','sc_c_misc','sc_cephes','sc_mach',
                                            'sc_cdf', 'sc_specfun'])
    cfg.setdefault('define_macros', []).extend(define_macros)
    config.add_extension('_ufuncs',
                         depends=ufuncs_dep,
                         sources=ufuncs_src,
                         extra_info=get_info("npymath"),
                         **cfg)

    # Extension _ufuncs_cxx
    ufuncs_cxx_src = ['_ufuncs_cxx.cxx', 'sf_error.c',
                      '_faddeeva.cxx', 'Faddeeva.cc',
                      '_wright.cxx', 'wright.cc']
    ufuncs_cxx_dep = (headers + ufuncs_cxx_src + cephes_src
                      + ['*.hh'])
    config.add_extension('_ufuncs_cxx',
                         sources=ufuncs_cxx_src,
                         depends=ufuncs_cxx_dep,
                         include_dirs=[curdir] + inc_dirs,
                         define_macros=define_macros,
                         extra_info=get_info("npymath"))

    cfg = dict(get_system_info('lapack_opt'))
    config.add_extension('_ellip_harm_2',
                         sources=['_ellip_harm_2.c', 'sf_error.c',],
                         **cfg
                         )

    # Cython API
    config.add_data_files('cython_special.pxd')
    
    cython_special_src = ['cython_special.c', 'sf_error.c', '_logit.c.src',
                          "amos_wrappers.c", "cdf_wrappers.c", "specfun_wrappers.c"]
    cython_special_dep = (headers + ufuncs_src + ufuncs_cxx_src + amos_src
                          + c_misc_src + cephes_src + mach_src + cdf_src
                          + specfun_src)
    cfg = dict(get_system_info('lapack_opt'))
    cfg.setdefault('include_dirs', []).extend([curdir] + inc_dirs + [numpy.get_include()])
    cfg.setdefault('libraries', []).extend(['sc_amos','sc_c_misc','sc_cephes','sc_mach',
                                            'sc_cdf', 'sc_specfun'])
    cfg.setdefault('define_macros', []).extend(define_macros)
    config.add_extension('cython_special',
                         depends=cython_special_dep,
                         sources=cython_special_src,
                         extra_info=get_info("npymath"),
                         **cfg)

    # combinatorics
    config.add_extension('_comb',
                         sources=['_comb.c'])

    # testing for _round.h
    config.add_extension('_test_round',
                         sources=['_test_round.c'],
                         depends=['_round.h', 'cephes/dd_idefs.h'],
                         include_dirs=[numpy.get_include()] + inc_dirs,
                         extra_info=get_info('npymath'))

    config.add_data_files('tests/*.py')
    config.add_data_files('tests/data/README')

    # regenerate npz data files
    makenpz = os.path.join(os.path.dirname(__file__),
                           'utils', 'makenpz.py')
    data_dir = os.path.join(os.path.dirname(__file__),
                            'tests', 'data')
    for name in ['boost', 'gsl', 'local']:
        subprocess.check_call([sys.executable, makenpz,
                               '--use-timestamp',
                               os.path.join(data_dir, name)])

    config.add_data_files('tests/data/*.npz')

    config.add_subpackage('_precompute')

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
