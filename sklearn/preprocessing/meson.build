py.extension_module(
  '_csr_polynomial_expansion',
  '_csr_polynomial_expansion.pyx',
  cython_args: cython_args,
  subdir: 'sklearn/preprocessing',
  install: true
)

py.extension_module(
  '_target_encoder_fast',
  '_target_encoder_fast.pyx',
  dependencies: [np_dep],
  override_options: ['cython_language=cpp'],
  cython_args: cython_args,
  subdir: 'sklearn/preprocessing',
  install: true
)
