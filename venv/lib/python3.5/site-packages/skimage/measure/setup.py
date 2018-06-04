#!/usr/bin/env python

from skimage._build import cython

import os
base_path = os.path.abspath(os.path.dirname(__file__))


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs

    config = Configuration('measure', parent_package, top_path)
    config.add_data_dir('tests')

    cython(['_ccomp.pyx'], working_path=base_path)
    cython(['_find_contours_cy.pyx'], working_path=base_path)
    cython(['_moments_cy.pyx'], working_path=base_path)
    cython(['_marching_cubes_classic_cy.pyx'], working_path=base_path)
    cython(['_marching_cubes_lewiner_cy.pyx'], working_path=base_path)
    cython(['_pnpoly.pyx'], working_path=base_path)

    config.add_extension('_ccomp', sources=['_ccomp.c'],
                         include_dirs=[get_numpy_include_dirs()])
    config.add_extension('_find_contours_cy', sources=['_find_contours_cy.c'],
                         include_dirs=[get_numpy_include_dirs()])
    config.add_extension('_moments_cy', sources=['_moments_cy.c'],
                         include_dirs=[get_numpy_include_dirs()])
    config.add_extension('_marching_cubes_classic_cy',
                         sources=['_marching_cubes_classic_cy.c'],
                         include_dirs=[get_numpy_include_dirs()])
    config.add_extension('_marching_cubes_lewiner_cy',
                         sources=['_marching_cubes_lewiner_cy.c'],
                         include_dirs=[get_numpy_include_dirs()])
    config.add_extension('_pnpoly', sources=['_pnpoly.c'],
                         include_dirs=[get_numpy_include_dirs(), '../_shared'])

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(maintainer='scikit-image Developers',
          maintainer_email='scikit-image@python.org',
          description='Graph-based Image-processing Algorithms',
          url='https://github.com/scikit-image/scikit-image',
          license='Modified BSD',
          **(configuration(top_path='').todict())
          )
