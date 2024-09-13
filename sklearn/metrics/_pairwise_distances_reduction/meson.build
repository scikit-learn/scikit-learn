# Note: the dependencies between different Cython files in
# _pairwise_distances_reduction is probably one of the most involved in
# scikit-learn. If you change this file make sure you build from scratch:
# rm -rf build; make dev-meson
# run a command like this:
# ninja -C build/cp312 -t missingdeps
# and make sure that the output is something like:
# No missing dependencies on generated files found.

# _pairwise_distances_reduction is cimported from other subpackages so this is
# needed for the cimport to work
_pairwise_distances_reduction_cython_tree = [
  fs.copyfile('__init__.py'),
  # We are in a sub-module of metrics, so we always need to have
  # sklearn/metrics/__init__.py copied to the build directory to avoid the
  # error:
  # relative cimport beyond main package is not allowed
  metrics_cython_tree
]

_classmode_pxd = fs.copyfile('_classmode.pxd')

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
  # TODO in principle this should go in py.exension_module below. This is
  # temporary work-around for dependency issue with .pyx.tp files. For more
  # details, see https://github.com/mesonbuild/meson/issues/13212
  depends: [_datasets_pair_pxd, _pairwise_distances_reduction_cython_tree, utils_cython_tree],
)
_datasets_pair = py.extension_module(
  '_datasets_pair',
  _datasets_pair_pyx,
  dependencies: [np_dep],
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
  # TODO in principle this should go in py.exension_module below. This is
  # temporary work-around for dependency issue with .pyx.tp files. For more
  # details, see https://github.com/mesonbuild/meson/issues/13212
  depends: [_base_pxd, _pairwise_distances_reduction_cython_tree,
            _datasets_pair_pxd, utils_cython_tree],
)
_base = py.extension_module(
  '_base',
  _base_pyx,
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
  # TODO in principle this should go in py.exension_module below. This is
  # temporary work-around for dependency issue with .pyx.tp files. For more
  # details, see https://github.com/mesonbuild/meson/issues/13212
  depends: [_middle_term_computer_pxd,
            _pairwise_distances_reduction_cython_tree,
            utils_cython_tree],
)
_middle_term_computer = py.extension_module(
  '_middle_term_computer',
  _middle_term_computer_pyx,
  dependencies: [np_dep],
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
    # TODO in principle this should go in py.exension_module below. This is
    # temporary work-around for dependency issue with .pyx.tp files. For more
    # details, see https://github.com/mesonbuild/meson/issues/13212
    depends: [_argkmin_pxd,
              _pairwise_distances_reduction_cython_tree,
              _datasets_pair_pxd, _base_pxd, _middle_term_computer_pxd],
      )
_argkmin = py.extension_module(
    '_argkmin',
    _argkmin_pyx,
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
    # TODO in principle this should go in py.exension_module below. This is
    # temporary work-around for dependency issue with .pyx.tp files. For more
    # details, see https://github.com/mesonbuild/meson/issues/13212
    depends: [_radius_neighbors_pxd,
              _datasets_pair_pxd, _base_pxd, _middle_term_computer_pxd,
              _pairwise_distances_reduction_cython_tree, utils_cython_tree],
)
_radius_neighbors = py.extension_module(
    '_radius_neighbors',
    _radius_neighbors_pyx,
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
  # TODO in principle this should go in py.exension_module below. This is
  # temporary work-around for dependency issue with .pyx.tp files. For more
  # details, see https://github.com/mesonbuild/meson/issues/13212
  depends: [_classmode_pxd,
            _argkmin_pxd, _pairwise_distances_reduction_cython_tree,
            _datasets_pair_pxd, _base_pxd, _middle_term_computer_pxd, utils_cython_tree],
)
_argkmin_classmode = py.extension_module(
  '_argkmin_classmode',
  _argkmin_classmode_pyx,
  dependencies: [np_dep, openmp_dep],
  override_options: ['cython_language=cpp'],
  cython_args: cython_args,
  # XXX: for some reason -fno-sized-deallocation is needed otherwise there is
  # an error with undefined symbol _ZdlPv at import time in manylinux wheels.
  # See https://github.com/scikit-learn/scikit-learn/issues/28596 for more details.
  cpp_args: ['-fno-sized-deallocation'],
  subdir: 'sklearn/metrics/_pairwise_distances_reduction',
  install: true
)

_radius_neighbors_classmode_pyx = custom_target(
  '_radius_neighbors_classmode_pyx',
  output: '_radius_neighbors_classmode.pyx',
  input: '_radius_neighbors_classmode.pyx.tp',
  command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
  # TODO in principle this should go in py.exension_module below. This is
  # temporary work-around for dependency issue with .pyx.tp files. For more
  # details, see https://github.com/mesonbuild/meson/issues/13212
  depends: [_classmode_pxd,
            _middle_term_computer_pxd, _radius_neighbors_pxd,
            _pairwise_distances_reduction_cython_tree,
            _datasets_pair_pxd, _base_pxd, utils_cython_tree],
)
_radius_neighbors_classmode = py.extension_module(
  '_radius_neighbors_classmode',
  _radius_neighbors_classmode_pyx,
  dependencies: [np_dep, openmp_dep],
  override_options: ['cython_language=cpp'],
  cython_args: cython_args,
  subdir: 'sklearn/metrics/_pairwise_distances_reduction',
  install: true
)
