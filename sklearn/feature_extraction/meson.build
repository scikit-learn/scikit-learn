py.extension_module(
  '_hashing_fast',
  [cython_gen_cpp.process('_hashing_fast.pyx'), utils_cython_tree],
  dependencies: [np_dep],
  limited_api: '3.12',
  subdir: 'sklearn/feature_extraction',
  install: true
)
