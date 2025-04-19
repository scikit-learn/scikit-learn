newrand_include = include_directories('src/newrand')
libsvm_include = include_directories('src/libsvm')

libsvm_skl = static_library(
  'libsvm-skl',
  ['src/libsvm/svm.cpp'],
)

py.extension_module(
  '_libsvm',
  ['_libsvm.pyx', utils_cython_tree],
  include_directories: [newrand_include, libsvm_include],
  link_with: libsvm_skl,
  cython_args: cython_args,
  subdir: 'sklearn/mkl/_svm',
  install: true
)
