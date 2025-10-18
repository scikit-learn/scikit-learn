py.extension_module(
  '_gradient_boosting',
  [cython_gen.process('_gradient_boosting.pyx')] + utils_cython_tree,
  dependencies: [np_dep],
  limited_api: '3.12',
  subdir: 'sklearn/ensemble',
  install: true
)

subdir('_hist_gradient_boosting')
