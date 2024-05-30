# Metrics is cimported from other subpackages so this is needed for the cimport
# to work
metrics_cython_tree = [
  fs.copyfile('__init__.py')
]
# Some metrics code cimports code from utils, we may as well copy all the necessary files
metrics_cython_tree += utils_cython_tree

_dist_metrics_pxd = custom_target(
  '_dist_metrics_pxd',
  output: '_dist_metrics.pxd',
  input: '_dist_metrics.pxd.tp',
  command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
  # Need to install the generated pxd because it is needed in other subpackages
  # Cython code, e.g. sklearn.cluster
  install_dir: sklearn_dir / 'metrics',
  install: true,
)
metrics_cython_tree += [_dist_metrics_pxd]

_dist_metrics_pyx = custom_target(
  '_dist_metrics_pyx',
  output: '_dist_metrics.pyx',
  input: '_dist_metrics.pyx.tp',
  command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
  # TODO in principle this should go in py.exension_module below. This is
  # temporary work-around for dependency issue with .pyx.tp files. For more
  # details, see https://github.com/mesonbuild/meson/issues/13212
  depends: metrics_cython_tree,
)

_dist_metrics = py.extension_module(
  '_dist_metrics',
  _dist_metrics_pyx,
  dependencies: [np_dep],
  cython_args: cython_args,
  subdir: 'sklearn/metrics',
  install: true
)

py.extension_module(
  '_pairwise_fast',
  ['_pairwise_fast.pyx', metrics_cython_tree],
  cython_args: cython_args,
  subdir: 'sklearn/metrics',
  install: true
)

subdir('_pairwise_distances_reduction')
subdir('cluster')
