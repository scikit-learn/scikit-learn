from __future__ import division, print_function, absolute_import

from os.path import join


def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info

    lapack_opt = get_info('lapack_opt', notfound_action=2)

    config = Configuration('interpolate', parent_package, top_path)

    fitpack_src = [join('fitpack', '*.f')]
    config.add_library('fitpack', sources=fitpack_src)

    config.add_extension('interpnd',
                         sources=['interpnd.c'])

    config.add_extension('_ppoly',
                         sources=['_ppoly.c'],
                         **lapack_opt)

    config.add_extension('_bspl',
                         sources=['_bspl.c'],
                         libraries=['fitpack'],
                         depends=['src/__fitpack.h'] + fitpack_src)

    config.add_extension('_fitpack',
                         sources=['src/_fitpackmodule.c'],
                         libraries=['fitpack'],
                         depends=(['src/__fitpack.h','src/multipack.h']
                                  + fitpack_src)
                         )

    config.add_extension('dfitpack',
                         sources=['src/fitpack.pyf'],
                         libraries=['fitpack'],
                         depends=fitpack_src,
                         )

    config.add_extension('_interpolate',
                         sources=['src/_interpolate.cpp'],
                         include_dirs=['src'],
                         depends=['src/interpolate.h'])

    config.add_data_dir('tests')

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
