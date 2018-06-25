# -*- coding: utf-8 -*-

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('loky', parent_package, top_path)
    config.add_subpackage('backend')

    return config
