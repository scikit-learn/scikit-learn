py.extension_module(
  '_online_lda_fast',
  ['_online_lda_fast.pyx', utils_cython_tree],
  cython_args: cython_args,
  subdir: 'sklearn/decomposition',
  install: true
)

py.extension_module(
  '_cdnmf_fast',
  '_cdnmf_fast.pyx',
  dependencies: [np_dep],
  cython_args: cython_args,
  subdir: 'sklearn/decomposition',
  install: true
)
