#!/usr/bin/env python

import numpy
import os

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('datasets', parent_package, top_path)
    config.add_data_dir('data')
    config.add_data_dir('descr')
    config.add_data_dir('images')
    config.add_data_dir(os.path.join('tests', 'data'))

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
