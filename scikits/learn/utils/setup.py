def configuration(parent_package='',top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration('utils', parent_package, top_path)

    config.add_subpackage('sparsetools')

    return config

