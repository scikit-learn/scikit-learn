fs.copyfile('_loss.pxd')

_loss_pyx = custom_target(
  '_loss_pyx',
  output: '_loss.pyx',
  input: '_loss.pyx.tp',
  command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@']
)

py.extension_module(
  '_loss',
  _loss_pyx,
  cython_args: cython_args,
  install: true,
  subdir: 'sklearn/_loss',
)
