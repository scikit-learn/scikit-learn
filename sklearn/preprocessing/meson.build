py.extension_module(
  '_csr_polynomial_expansion',
  ['_csr_polynomial_expansion.pyx', utils_cython_tree],
  cython_args: cython_args,
  subdir: 'sklearn/preprocessing',
  install: true
)

py.extension_module(
  '_target_encoder_fast',
  ['_target_encoder_fast.pyx', utils_cython_tree],
  override_options: ['cython_language=cpp'],
  cython_args: cython_args,
  subdir: 'sklearn/preprocessing',
  install: true
)
