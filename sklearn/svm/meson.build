newrand_include = include_directories('src/newrand')
libsvm_include = include_directories('src/libsvm')
liblinear_include = include_directories('src/liblinear')

_newrand = py.extension_module(
  '_newrand',
  '_newrand.pyx',
  override_options: ['cython_language=cpp'],
  include_directories: [newrand_include],
  cython_args: cython_args,
  subdir: 'sklearn/svm',
  install: true
)

libsvm_skl = static_library(
  'libsvm-skl',
  ['src/libsvm/libsvm_template.cpp'],
)

py.extension_module(
  '_libsvm',
  ['_libsvm.pyx', utils_cython_tree],
  include_directories: [newrand_include, libsvm_include],
  link_with: libsvm_skl,
  cython_args: cython_args,
  subdir: 'sklearn/svm',
  install: true
)

py.extension_module(
  '_libsvm_sparse',
  ['_libsvm_sparse.pyx', utils_cython_tree],
  include_directories: [newrand_include, libsvm_include],
  link_with: libsvm_skl,
  cython_args: cython_args,
  subdir: 'sklearn/svm',
  install: true
)

liblinear_skl = static_library(
  'liblinear-skl',
  ['src/liblinear/linear.cpp', 'src/liblinear/tron.cpp'],
)

py.extension_module(
  '_liblinear',
  ['_liblinear.pyx', utils_cython_tree],
  include_directories: [newrand_include, liblinear_include],
  link_with: [liblinear_skl],
  cython_args: cython_args,
  subdir: 'sklearn/svm',
  install: true
)
