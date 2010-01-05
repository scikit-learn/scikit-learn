#!/usr/bin/env python

def configuration(parent_package='',top_path=None):
    # The following two lines with `return config` constitutes a
    # minimal contents of configuration(..) that is suitable for pure
    # Python packages.
    from numpy.distutils.misc_util import Configuration
    config = Configuration( 'mlp', parent_package, top_path )

    # include test scripts from tests
    #config.add_data_dir('tests')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup( **configuration( top_path = '' ).todict() )
