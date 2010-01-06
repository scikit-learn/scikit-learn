def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    print "============================================="
    print "parent package is %s" % parent_package
    print "top path is %s" % top_path
    print "============================================="
    config = Configuration('machine',parent_package,top_path)
    # pyem should be suppresed as some points
    config.add_subpackage('pyem')
    config.add_subpackage('em')
    config.add_subpackage('manifold_learning')
    config.add_subpackage('svm')
    #config.make_svn_version_py()  # installs __svn_version__.py
    #config.make_config_py()
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
