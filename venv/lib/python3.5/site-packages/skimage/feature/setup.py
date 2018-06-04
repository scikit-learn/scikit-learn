#!/usr/bin/env python

import os
from skimage._build import cython

base_path = os.path.abspath(os.path.dirname(__file__))


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs

    config = Configuration('feature', parent_package, top_path)
    config.add_data_dir('tests')

    cython(['corner_cy.pyx'], working_path=base_path)
    cython(['censure_cy.pyx'], working_path=base_path)
    cython(['orb_cy.pyx'], working_path=base_path)
    cython(['brief_cy.pyx'], working_path=base_path)
    cython(['_texture.pyx'], working_path=base_path)
    cython(['_hessian_det_appx.pyx'], working_path=base_path)
    cython(['_hoghistogram.pyx'], working_path=base_path)
    cython(['_haar.pyx'], working_path=base_path)

    config.add_extension('corner_cy', sources=['corner_cy.c'],
                         include_dirs=[get_numpy_include_dirs()])
    config.add_extension('censure_cy', sources=['censure_cy.c'],
                         include_dirs=[get_numpy_include_dirs()])
    config.add_extension('orb_cy', sources=['orb_cy.c'],
                         include_dirs=[get_numpy_include_dirs()])
    config.add_extension('brief_cy', sources=['brief_cy.c'],
                         include_dirs=[get_numpy_include_dirs()])
    config.add_extension('_texture', sources=['_texture.c'],
                         include_dirs=[get_numpy_include_dirs(), '../_shared'])
    config.add_extension('_hessian_det_appx', sources=['_hessian_det_appx.c'],
                         include_dirs=[get_numpy_include_dirs()])
    config.add_extension('_hoghistogram', sources=['_hoghistogram.c'],
                         include_dirs=[get_numpy_include_dirs(), '../_shared'])
    config.add_extension('_haar', sources=['_haar.cpp'],
                         include_dirs=[get_numpy_include_dirs(), '../_shared'],
                         language="c++")

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(maintainer='scikit-image Developers',
          author='scikit-image Developers',
          maintainer_email='scikit-image@python.org',
          description='Features',
          url='https://github.com/scikit-image/scikit-image',
          license='SciPy License (BSD Style)',
          **(configuration(top_path='').todict())
          )
