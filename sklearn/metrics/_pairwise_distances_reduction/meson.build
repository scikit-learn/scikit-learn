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
_datasets_pair = py.extension_module(
    '_datasets_pair',
    [pxd, pyx, _pairwise_distances_reduction_cython_tree, utils_cython_tree],
    dependencies: [np_dep, openmp_dep],
    override_options: ['cython_language=cpp'],
    cython_args: cython_args,
    subdir: 'sklearn/metrics/_pairwise_distances_reduction',
    install: true
)

_base_pxd = custom_target(
    '_base_pxd',
    output: '_base.pxd',
    input: '_base.pxd.tp',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@']
  )
_base_pyx = custom_target(
    '_base_pyx',
    output: '_base.pyx',
    input: '_base.pyx.tp',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
  )
_base = py.extension_module(
    '_base',
    [pxd, pyx, _pairwise_distances_reduction_cython_tree, utils_cython_tree],
    dependencies: [np_dep, openmp_dep],
    override_options: ['cython_language=cpp'],
    cython_args: cython_args,
    subdir: 'sklearn/metrics/_pairwise_distances_reduction',
    install: true
)

_middle_term_computer_pxd = custom_target(
    '_middle_term_computer_pxd',
    output: '_middle_term_computer.pxd',
    input: '_middle_term_computer.pxd.tp',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@']
  )
_middle_term_computer_pyx = custom_target(
    '_middle_term_computer_pyx',
    output: '_middle_term_computer.pyx',
    input: '_middle_term_computer.pyx.tp',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
  )
_middle_term_computer = py.extension_module(
    '_middle_term_computer',
    [pxd, pyx, _pairwise_distances_reduction_cython_tree, utils_cython_tree],
    dependencies: [np_dep, openmp_dep],
    override_options: ['cython_language=cpp'],
    cython_args: cython_args,
    subdir: 'sklearn/metrics/_pairwise_distances_reduction',
    install: true
)

_argkmin_pxd = custom_target(
    '_argkmin_pxd',
    output: '_argkmin.pxd',
    input: '_argkmin.pxd.tp',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@']
  )
_argkmin_pyx = custom_target(
    '_argkmin_pyx',
    output: '_argkmin.pyx',
    input: '_argkmin.pyx.tp',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
  )
_argkmin = py.extension_module(
    '_argkmin',
    [pxd, pyx, _pairwise_distances_reduction_cython_tree, utils_cython_tree],
    dependencies: [np_dep, openmp_dep],
    override_options: ['cython_language=cpp'],
    cython_args: cython_args,
    subdir: 'sklearn/metrics/_pairwise_distances_reduction',
    install: true
)

_radius_neighbors_pxd = custom_target(
    '_radius_neighbors_pxd',
    output: '_radius_neighbors.pxd',
    input: '_radius_neighbors.pxd.tp',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@']
  )
_radius_neighbors_pyx = custom_target(
    '_radius_neighbors_pyx',
    output: '_radius_neighbors.pyx',
    input: '_radius_neighbors.pyx.tp',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
  )
_radius_neighbors = py.extension_module(
    '_radius_neighbors',
    [pxd, pyx, _pairwise_distances_reduction_cython_tree, utils_cython_tree],
    dependencies: [np_dep, openmp_dep],
    override_options: ['cython_language=cpp'],
    cython_args: cython_args,
    subdir: 'sklearn/metrics/_pairwise_distances_reduction',
    install: true
)

_argkmin_classmode_pyx = custom_target(
  '_argkmin_classmode_pyx',
  output: '_argkmin_classmode.pyx',
  input: '_argkmin_classmode.pyx.tp',
  command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
)
_argkmin_classmode = py.extension_module(
  '_argkmin_classmode',
  [pyx, _pairwise_distances_reduction_cython_tree, _datasets_pair_pxd, _base_pxd, _middle_term_computer_pxd, utils_cython_tree],
  dependencies: [np_dep],
  override_options: ['cython_language=cpp'],
  cython_args: cython_args,
  subdir: 'sklearn/metrics/_pairwise_distances_reduction',
  install: true
)

_radius_neighbors_classmode_pyx = custom_target(
  '_radius_neighbors_classmode_pyx',
  output: '_radius_neighbors_classmode.pyx',
  input: '_radius_neighbors_classmode.pyx.tp',
  command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
)
_radius_neighbors_classmode = py.extension_module(
  '_radius_neighbors_classmode',
  [pyx, _pairwise_distances_reduction_cython_tree, _datasets_pair_pxd, _base_pxd, utils_cython_tree],
  dependencies: [np_dep],
  override_options: ['cython_language=cpp'],
  cython_args: cython_args,
  subdir: 'sklearn/metrics/_pairwise_distances_reduction',
  install: true
)
