libs_sources = {
  '_dbscan_inner':
    {'sources': ['_dbscan_inner.pyx'], 'override_options': ['cython_language=cpp'], 'include_np': true},
  '_hierarchical_fast':
    {'sources': ['_hierarchical_fast.pyx'], 'override_options': ['cython_language=cpp'], 'include_np': true},
  '_k_means_common':
    {'sources': ['_k_means_common.pyx'], 'include_np': true},
  '_k_means_lloyd':
    {'sources': ['_k_means_lloyd.pyx'], 'include_np': true},
  '_k_means_elkan':
    {'sources': ['_k_means_elkan.pyx'], 'include_np': true},
  '_k_means_minibatch':
    {'sources': ['_k_means_minibatch.pyx'], 'include_np': true},
}

# TODO this fails the first ninja build but succeeds the second. I am going to
# guess that there is some dependencies that needs to be defined, i.e. this
# depends on sklearn/metrics. Error is:
# /home/lesteve/dev/scikit-learn/sklearn/cluster/_hierarchical_fast.pyx:6:0: 'sklearn/metrics/_dist_metrics.pxd' not found

cython_args = [
    '--include-dir',
    meson.global_build_root(),
]

foreach ext_name, ext_dict : libs_sources
    py.extension_module(
        ext_name,
        ext_dict.get('sources'),
        cython_args: cython_args,
        dependencies: [np_dep],
        # dependencies: ext_dict.get('deps', ''),
        override_options : ext_dict.get('override_options', []),
        subdir: 'sklearn/cluster',
        install: true
    )
endforeach
