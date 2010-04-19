from os.path import join
import numpy

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
                         sources=[join('src', 'svm.cpp'), 
                                  join('src', 'libsvm.c'),
                                  ],
                         include_dirs=[numpy.get_include()],
                         depends=[join('src', 'svm.h'),
                                 join('src', 'libsvm_helper.c'),
                                  ])
    config.add_extension('liblinear',
                         sources=[join('src', 'linear.cpp'), 
                                  join('src', 'liblinear.c'),
                                  join('src', 'tron.cpp'),
                                  join('src', 'blas', 'daxpy.c'),
                                  join('src', 'blas', 'ddot.c'),
                                  join('src', 'blas', 'dnrm2.c'),
                                  join('src', 'blas', 'dscal.c'),                                  
                                  ],
                         include_dirs=[numpy.get_include()],
                         depends=[join('src', 'linear.h'),
                                  join('src', 'liblinear_helper.c'),
                                  join('src', 'tron.h'),
                                  join('src', 'blas', 'blas.h'),
                                  join('src', 'blas', 'blasp.h'),
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
