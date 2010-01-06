def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    print "============================================="
    print "parent package is %s" % parent_package
    print "top path is %s" % top_path
    print "============================================="
    config = Configuration('',parent_package,top_path)
    config.add_subpackage('utils')
    config.add_subpackage('datasets')
    config.add_subpackage('machine')
    config.add_subpackage('common')
    #config.make_svn_version_py()  # installs __svn_version__.py
    #config.make_config_py()
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
