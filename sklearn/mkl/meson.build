# .pyx is generated, so this is needed to make Cython compilation work
linear_model_cython_tree = [
  fs.copyfile('__init__.py'),
]

py.extension_module(
  '_utils',
  ['_utils.pyx', utils_cython_tree],
  cython_args: cython_args,
  subdir: 'sklearn/mkl',
  install: true
)

subdir('_algo')
subdir('_svm')
