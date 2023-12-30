_loss = custom_target('_loss',
  output: '_loss.pyx',
  input: '_loss.pyx.tp',
  command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@']
)

loss = py.extension_module(
  '_loss',
  '_loss.pyx',
  cython_args: cython_args,
  install: true,
  subdir: 'sklearn/_loss',
)
