# TODO Needed for relative imports in Cython?
fs.copyfile('__init__.py')
fs.copyfile('_classmode.pxd')

name_list_with_pxd = ['_datasets_pair', '_middle_term_computer', '_base', '_argkmin', '_radius_neighbors']

foreach name: name_list_with_pxd
  pxd = custom_target(
    name + '_pxd',
    output: name + '.pxd',
    input: name + '.pxd.tp',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@']
  )
  pyx = custom_target(
    name + '_pyx',
    output: name + '.pyx',
    input: name + '.pyx.tp',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@']
  )
  py.extension_module(
    name,
    [pxd, pyx],
    dependencies: [np_dep, openmp_dep],
    override_options: ['cython_language=cpp'],
    cython_args: cython_args,
    subdir: 'sklearn/metrics/_pairwise_distances_reduction',
    install: true

)
endforeach

name_list_without_pxd = ['_argkmin_classmode', '_radius_neighbors_classmode']

foreach name: name_list_without_pxd
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
    override_options: ['cython_language=cpp'],
    cython_args: cython_args,
    subdir: 'sklearn/metrics/_pairwise_distances_reduction',
    install: true

)
endforeach
