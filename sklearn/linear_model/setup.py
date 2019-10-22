import os
import numpy

from .._build_utils import gen_from_templates


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('linear_model', parent_package, top_path)

    libraries = []
    if os.name == 'posix':
        libraries.append('m')

    config.add_extension('cd_fast',
                         sources=['cd_fast.pyx'],
                         include_dirs=numpy.get_include(),
                         libraries=libraries)

    config.add_extension('sgd_fast',
                         sources=['sgd_fast.pyx'],
                         include_dirs=numpy.get_include(),
                         libraries=libraries)

    # generate sag_fast from template
    gen_from_templates('sklearn/linear_model/sag_fast.pyx.tp')

    config.add_extension('sag_fast',
                         sources=['sag_fast.pyx'],
                         include_dirs=numpy.get_include())

    # add other directories
    config.add_subpackage('tests')

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
