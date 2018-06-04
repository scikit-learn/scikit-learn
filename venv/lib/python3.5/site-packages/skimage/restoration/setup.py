#!/usr/bin/env python

import os

from skimage._build import cython

base_path = os.path.abspath(os.path.dirname(__file__))


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs

    config = Configuration('restoration', parent_package, top_path)
    config.add_data_dir('tests')

    cython(['_unwrap_1d.pyx'], working_path=base_path)
    cython(['_unwrap_2d.pyx'], working_path=base_path)
    cython(['_unwrap_3d.pyx'], working_path=base_path)
    cython(['_denoise_cy.pyx'], working_path=base_path)
    cython(['_nl_means_denoising.pyx'], working_path=base_path)

    config.add_extension('_unwrap_1d', sources=['_unwrap_1d.c'],
                         include_dirs=[get_numpy_include_dirs()])
    unwrap_sources_2d = ['_unwrap_2d.c', 'unwrap_2d_ljmu.c']
    config.add_extension('_unwrap_2d', sources=unwrap_sources_2d,
                         include_dirs=[get_numpy_include_dirs()])
    unwrap_sources_3d = ['_unwrap_3d.c', 'unwrap_3d_ljmu.c']
    config.add_extension('_unwrap_3d', sources=unwrap_sources_3d,
                         include_dirs=[get_numpy_include_dirs()])
    config.add_extension('_denoise_cy', sources=['_denoise_cy.c'],
        include_dirs=[get_numpy_include_dirs(), '../_shared'])
    config.add_extension('_nl_means_denoising',
                         sources=['_nl_means_denoising.c'],
                         include_dirs=[get_numpy_include_dirs(),
                                       '../_shared'])

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(maintainer='scikit-image Developers',
          author='scikit-image Developers',
          maintainer_email='scikit-image@python.org',
          description='Restoration',
          url='https://github.com/scikit-image/scikit-image',
          license='SciPy License (BSD Style)',
          **(configuration(top_path='').todict())
          )
