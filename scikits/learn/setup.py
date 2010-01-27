def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('learn',parent_package,top_path)

    config.add_subpackage('datasets')
    config.add_subpackage('common')
    config.add_subpackage('machine')
    config.add_subpackage('utils')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
