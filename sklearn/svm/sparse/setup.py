from os.path import join
import numpy


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('sparse', parent_package, top_path)

    libsvm_sparse_sources = ['libsvm.c']
    config.add_extension('libsvm', libraries=['libsvm-skl'],
            sources=libsvm_sparse_sources,
            include_dirs=[numpy.get_include(), join('..', 'src', 'libsvm')],
            depends=[join('..', 'src', 'libsvm', 'svm.h'),
                join('..', 'src', 'libsvm', 'libsvm_sparse_helper.c')],)

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
