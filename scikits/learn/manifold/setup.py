def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('manifold',parent_package,top_path)
    config.add_subpackage('compression')
    config.add_data_dir('examples')
    config.add_subpackage('projection')
    config.add_subpackage('regression')
    config.add_subpackage('stats')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
