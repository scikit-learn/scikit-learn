# Created by Pearu Peterson, August 2002
from __future__ import division, print_function, absolute_import


from os.path import join


def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('fftpack',parent_package, top_path)

    config.add_data_dir('tests')

    dfftpack_src = [join('src/dfftpack','*.f')]
    config.add_library('dfftpack', sources=dfftpack_src)

    fftpack_src = [join('src/fftpack','*.f')]
    config.add_library('fftpack', sources=fftpack_src)

    sources = ['fftpack.pyf','src/zfft.c','src/drfft.c','src/zrfft.c',
               'src/zfftnd.c', 'src/dct.c.src', 'src/dst.c.src']

    config.add_extension('_fftpack',
        sources=sources,
        libraries=['dfftpack', 'fftpack'],
        include_dirs=['src'],
        depends=(dfftpack_src + fftpack_src))

    config.add_extension('convolve',
        sources=['convolve.pyf','src/convolve.c'],
        libraries=['dfftpack'],
        depends=dfftpack_src,
    )
    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
