_pairwise_distances_reduction_cython_tree = [
  fs.copyfile('__init__.py'),
  fs.copyfile('_classmode.pxd')
]

_datasets_pair_pxd = custom_target(
    '_datasets_pair_pxd',
    output: '_datasets_pair.pxd',
    input: '_datasets_pair.pxd.tp',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@']
  )
_datasets_pair_pyx = custom_target(
    '_datasets_pair_pyx',
    output: '_datasets_pair.pyx',
    input: '_datasets_pair.pyx.tp',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
  )
py.extension_module(
    name,
    [pxd, pyx, _pairwise_distances_reduction_cython_tree],
    dependencies: [np_dep, openmp_dep],
    override_options: ['cython_language=cpp'],
    cython_args: cython_args,
    subdir: 'sklearn/metrics/_pairwise_distances_reduction',
    install: true
)

name_list_with_pxd = ['_middle_term_computer', '_base', '_argkmin', '_radius_neighbors']

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
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
  )
  py.extension_module(
    name,
    [pxd, pyx, _datasets_pair_pxd, _pairwise_distances_reduction_cython_tree],
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
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
  )
  py.extension_module(
    name,
    [pyx, _pairwise_distances_reduction_cython_tree, _datasets_pair_pxd],
    dependencies: [np_dep],
    override_options: ['cython_language=cpp'],
    cython_args: cython_args,
    subdir: 'sklearn/metrics/_pairwise_distances_reduction',
    install: true
)
endforeach
