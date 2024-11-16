py.extension_module(
  '_utils',
  ['_utils.pyx', utils_cython_tree],
  cython_args: cython_args,
  subdir: 'sklearn/manifold',
  install: true
)

py.extension_module(
  '_barnes_hut_tsne',
  '_barnes_hut_tsne.pyx',
  dependencies: [np_dep, openmp_dep],
  cython_args: cython_args,
  subdir: 'sklearn/manifold',
  install: true
)
