# TODO Needed for relative imports in Cython?
fs.copyfile('__init__.py')

_dist_metrics_pyx = custom_target(
  '_dist_metrics_pyx',
  output: '_dist_metrics.pyx',
  input: '_dist_metrics.pyx.tp',
  command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@']
)

_dist_metrics_pxd = custom_target(
  '_dist_metrics_pxd',
  output: '_dist_metrics.pxd',
  input: '_dist_metrics.pxd.tp',
  command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
  # TODO Need to install because this is needed in other modules i.e. sklearn/cluster
  install_dir: sklearn_dir / 'metrics',
  install: true,
)

_dist_metrics = py.extension_module(
        '_dist_metrics',
        [_dist_metrics_pyx, _dist_metrics_pxd],
        # cython_args: cython_args,
        dependencies: [np_dep],
        cython_args: cython_args,
        subdir: 'sklearn/metrics',
        install: true
)

py.extension_module(
        '_pairwise_fast',
        '_pairwise_fast.pyx',
        # cython_args: cython_args,
        dependencies: [np_dep],
        cython_args: cython_args,
        subdir: 'sklearn/metrics',
        install: true
)

subdir('cluster')
subdir('_pairwise_distances_reduction')
