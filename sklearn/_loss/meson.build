# .pyx is generated, so this is needed to make Cython compilation work
_loss_cython_tree = [
  fs.copyfile('_loss.pxd')
]

_loss_pyx = custom_target(
  '_loss_pyx',
  output: '_loss.pyx',
  input: '_loss.pyx.tp',
  command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
  # TODO in principle this should go in py.exension_module below. This is
  # temporary work-around for dependency issue with .pyx.tp files. For more
  # details, see https://github.com/mesonbuild/meson/issues/13212
  depends: _loss_cython_tree,
)

py.extension_module(
  '_loss',
  _loss_pyx,
  dependencies: [openmp_dep],
  cython_args: cython_args,
  install: true,
  subdir: 'sklearn/_loss',
)
