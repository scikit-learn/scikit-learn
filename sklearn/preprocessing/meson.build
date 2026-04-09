py.extension_module(
  '_csr_polynomial_expansion',
  [cython_gen.process('_csr_polynomial_expansion.pyx'), utils_cython_tree],
  subdir: 'sklearn/preprocessing',
  install: true
)

py.extension_module(
  '_target_encoder_fast',
  [cython_gen_cpp.process('_target_encoder_fast.pyx'), utils_cython_tree],
  subdir: 'sklearn/preprocessing',
  install: true
)
