fs.copyfile('__init__.py')

py.extension_module(
    '_cd_fast',
    '_cd_fast.pyx',
    dependencies: [np_dep],
    cython_args: cython_args,
    subdir: 'sklearn/linear_model',
    install: true
)

name_list = ['_sgd_fast', '_sag_fast']

fs.copyfile('_sgd_fast.pxd')

foreach name: name_list
  pyx = custom_target(
    name + '_pyx',
    output: name + '.pyx',
    input: name + '.pyx.tp',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@']
  )
  py.extension_module(
    name,
    pyx,
    dependencies: [np_dep],
    cython_args: cython_args,
    subdir: 'sklearn/linear_model',
    install: true
)
endforeach
