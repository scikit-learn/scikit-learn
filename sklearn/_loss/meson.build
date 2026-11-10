# .pyx is generated, so this is needed to make Cython compilation work
# We add sklearn_root_cython_tree and __init__.py so Cython can detect the
# package hierarchy and set the correct __module__ on extension types.
_loss_cython_tree = [
  sklearn_root_cython_tree,
  fs.copyfile('__init__.py'),
  fs.copyfile('_loss.pxd'),
]

_loss_pyx = custom_target(
  '_loss_pyx',
  output: '_loss.pyx',
  input: '_loss.pyx.tp',
  command: [tempita, '@INPUT@', '-o', '@OUTDIR@'],
  # TODO in principle this should go in py.extension_module below. This is
  # temporary work-around for dependency issue with .pyx.tp files. For more
  # details, see https://github.com/mesonbuild/meson/issues/13212
  depends: _loss_cython_tree,
)

py.extension_module(
  '_loss',
  cython_gen.process(_loss_pyx),
  dependencies: [openmp_dep],
  install: true,
  subdir: 'sklearn/_loss',
)
