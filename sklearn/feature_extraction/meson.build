py.extension_module(
  '_hashing_fast',
  ['_hashing_fast.pyx', utils_cython_tree],
  dependencies: [np_dep],
  override_options: ['cython_language=cpp'],
  cython_args: cython_args,
  subdir: 'sklearn/feature_extraction',
  install: true
)
