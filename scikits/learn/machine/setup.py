def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('machine',parent_package,top_path)
    config.add_subpackage('em')
    #config.add_subpackage('em2')
    config.add_subpackage('manifold_learning')
    config.add_subpackage('svm')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
