from os.path import join

def configuration(parent_package='', top_path=None, package_name='svm'):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(package_name,parent_package,top_path)
    config.add_subpackage('*')
    config.add_extension('libsvm_',
                         define_macros=[('LIBSVM_EXPORTS', None),
                                        ('LIBSVM_DLL', None)],
                         sources=[join('libsvm-2.82', 'svm.cpp')],
                         depends=[join('libsvm-2.82', 'svm.h')])
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='',
                          package_name='scipy.sandbox.svm').todict())
