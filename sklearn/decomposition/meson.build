py.extension_module(
  '_online_lda_fast',
  [cython_gen.process('_online_lda_fast.pyx'), utils_cython_tree],
  limited_api: limited_api,
  subdir: 'sklearn/decomposition',
  install: true
)

py.extension_module(
  '_cdnmf_fast',
  cython_gen.process('_cdnmf_fast.pyx'),
  dependencies: [np_dep],
  limited_api: limited_api,
  subdir: 'sklearn/decomposition',
  install: true
)
