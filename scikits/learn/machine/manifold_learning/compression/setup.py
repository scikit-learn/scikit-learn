from os.path import join

import os.path
import numpy

def configuration(parent_package='', top_path=None, package_name='compression'):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(package_name,parent_package,top_path)
    config.add_subpackage('*')
    config.add_extension('cost_function._cost_function',
                         sources=["cost_function/cost_function.cpp"],
                         include_dirs=[os.path.dirname(os.path.abspath(__file__)) + '/..'],)
    config.add_data_dir('tests')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='',
                          package_name='compression').todict())
    #setup(configuration=configuration)
