def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(',parent_package,top_path)
    config.add_subpackage('em')
    config.add_subpackage('manifold')
    config.add_subpackage('svm')
    config.add_subpackage('utils')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
