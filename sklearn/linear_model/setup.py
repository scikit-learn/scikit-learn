import os

import numpy

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('linear_model', parent_package, top_path)

    libraries = []
    if os.name == 'posix':
        libraries.append('m')

    config.add_extension('_cd_fast',
                         sources=['_cd_fast.pyx'],
                         include_dirs=numpy.get_include(),
                         libraries=libraries)

    config.add_extension('_sgd_fast',
                         sources=['_sgd_fast.pyx'],
                         include_dirs=numpy.get_include(),
                         libraries=libraries)

    # generate sag_fast from template
    sag_cython_file = 'sklearn/linear_model/_sag_fast.pyx.tp'
    sag_file = sag_cython_file.replace('.tp', '')

    if not (os.path.exists(sag_file) and
            os.stat(sag_cython_file).st_mtime < os.stat(sag_file).st_mtime):

        with open(sag_cython_file, "r") as f:
            tmpl = f.read()
        from Cython import Tempita # noqa
        tmpl_ = Tempita.sub(tmpl)

        with open(sag_file, "w") as f:
            f.write(tmpl_)

    config.add_extension('_sag_fast',
                         sources=['_sag_fast.pyx'],
                         include_dirs=numpy.get_include())

    # add other directories
    config.add_subpackage('tests')

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
