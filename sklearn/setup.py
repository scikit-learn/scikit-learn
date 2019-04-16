import os

from sklearn._build_utils import maybe_cythonize_extensions


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    import numpy

    # needs to be called during build otherwise show_version may fail sometimes
    get_info('blas_opt', 0)

    libraries = []
    if os.name == 'posix':
        libraries.append('m')

    config = Configuration('sklearn', parent_package, top_path)

    # submodules with build utilities
    config.add_subpackage('__check_build')
    config.add_subpackage('_build_utils')

    # submodules which do not have their own setup.py
    # we must manually add sub-submodules & tests
    config.add_subpackage('compose')
    config.add_subpackage('compose/tests')
    config.add_subpackage('covariance')
    config.add_subpackage('covariance/tests')
    config.add_subpackage('cross_decomposition')
    config.add_subpackage('cross_decomposition/tests')
    config.add_subpackage('feature_selection')
    config.add_subpackage('feature_selection/tests')
    config.add_subpackage('gaussian_process')
    config.add_subpackage('gaussian_process/tests')
    config.add_subpackage('mixture')
    config.add_subpackage('mixture/tests')
    config.add_subpackage('model_selection')
    config.add_subpackage('model_selection/tests')
    config.add_subpackage('neural_network')
    config.add_subpackage('neural_network/tests')
    config.add_subpackage('preprocessing')
    config.add_subpackage('preprocessing/tests')
    config.add_subpackage('semi_supervised')
    config.add_subpackage('semi_supervised/tests')

    # submodules which have their own setup.py
    config.add_subpackage('cluster')
    config.add_subpackage('datasets')
    config.add_subpackage('decomposition')
    config.add_subpackage('ensemble')
    config.add_subpackage('externals')
    config.add_subpackage('feature_extraction')
    config.add_subpackage('manifold')
    config.add_subpackage('metrics')
    config.add_subpackage('neighbors')
    config.add_subpackage('tree')
    config.add_subpackage('utils')
    config.add_subpackage('svm')
    config.add_subpackage('linear_model')

    # add cython extension module for isotonic regression
    config.add_extension('_isotonic',
                         sources=['_isotonic.pyx'],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries,
                         )

    # add the test directory
    config.add_subpackage('tests')

    maybe_cythonize_extensions(top_path, config)

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
