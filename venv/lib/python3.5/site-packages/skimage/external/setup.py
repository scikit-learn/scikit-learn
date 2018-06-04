#!/usr/bin/env python

import os.path

base_path = os.path.abspath(os.path.dirname(__file__))


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs

    config = Configuration('external', parent_package, top_path)

    config.add_extension('tifffile._tifffile',
                         sources=['tifffile/tifffile.c'],
                         include_dirs=[get_numpy_include_dirs()])
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(maintainer='scikit-image Developers',
          maintainer_email='scikit-image@python.org',
          description='External Libaries',
          url='https://github.com/scikit-image/scikit-image',
          license='Modified BSD',
          **(configuration(top_path='').todict())
          )
