# -*- coding: utf-8 -*-

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('joblib', parent_package, top_path)
    config.add_subpackage('externals')

    return config
