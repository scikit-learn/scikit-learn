py.extension_module(
  '_gradient_boosting',
  ['_gradient_boosting.pyx'] + utils_cython_tree,
  dependencies: [np_dep],
  cython_args: cython_args,
  subdir: 'sklearn/ensemble',
  install: true
)

subdir('_hist_gradient_boosting')
