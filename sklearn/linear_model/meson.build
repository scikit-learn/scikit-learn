# .pyx is generated, so this is needed to make Cython compilation work
linear_model_cython_tree = [
  fs.copyfile('__init__.py'),
  fs.copyfile('_sgd_fast.pxd'),
]

py.extension_module(
  '_cd_fast',
  ['_cd_fast.pyx', utils_cython_tree],
  cython_args: cython_args,
  subdir: 'sklearn/linear_model',
  install: true
)

name_list = ['_sgd_fast', '_sag_fast']

foreach name: name_list
  pyx = custom_target(
    name + '_pyx',
    output: name + '.pyx',
    input: name + '.pyx.tp',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
    # TODO in principle this should go in py.exension_module below. This is
    # temporary work-around for dependency issue with .pyx.tp files. For more
    # details, see https://github.com/mesonbuild/meson/issues/13212
    depends: [linear_model_cython_tree, utils_cython_tree],
  )
  py.extension_module(
    name,
    pyx,
    cython_args: cython_args,
    subdir: 'sklearn/linear_model',
    install: true
)
endforeach
