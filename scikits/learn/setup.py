from os.path import join
import numpy

# Workaround to enforce building cython extensions while
# maintaining compatibility with NumPy
# Found at http://old.nabble.com/problem-with-numpy.distutils-and-Cython-td25100957.html
# Introduced here by Yaroslav Halchenko <debian@onerussian.com> 2010-04-06
from numpy.distutils.command import build_src
import Cython
import Cython.Compiler.Main
build_src.Pyrex = Cython
build_src.have_pyrex = True


def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('learn',parent_package,top_path)
    config.add_subpackage('em')
    config.add_subpackage('datasets')
    config.add_subpackage('feature_selection')
    config.add_subpackage('glm')
    config.add_subpackage('manifold')
    config.add_subpackage('utils')
    config.add_extension('libsvm',
                         define_macros=[('LIBSVM_EXPORTS', None),
                                        ('LIBSVM_DLL',     None)],
                         sources=[join('src', 'svm.cpp'), 
                                  join('src', 'libsvm.pyx'),
                                  ],
                         include_dirs=[numpy.get_include(),
                                       join('scikits', 'learn', 'src')],
                         depends=[join('src', 'svm.h'),
                                 join('src', 'libsvm_helper.c'),
                                  ])

    config.add_extension('BallTree',
                         sources=[join('src', 'BallTree.cpp')],
                         include_dirs=[numpy.get_include()]
                         )


    return config

    config.add_subpackage('utils')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
