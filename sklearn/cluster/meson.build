cluster_extension_metadata = {
  '_dbscan_inner':
    {'sources': ['_dbscan_inner.pyx'], 'override_options': ['cython_language=cpp']},
  '_hierarchical_fast':
    {'sources': ['_hierarchical_fast.pyx', metrics_cython_tree],
     'override_options': ['cython_language=cpp']},
  '_k_means_common':
    {'sources': ['_k_means_common.pyx']},
  '_k_means_lloyd':
    {'sources': ['_k_means_lloyd.pyx']},
  '_k_means_elkan':
    {'sources': ['_k_means_elkan.pyx']},
  '_k_means_minibatch':
    {'sources': ['_k_means_minibatch.pyx']},
}

foreach ext_name, ext_dict : cluster_extension_metadata
  py.extension_module(
    ext_name,
    [ext_dict.get('sources'), utils_cython_tree],
    dependencies: [np_dep, openmp_dep],
    override_options : ext_dict.get('override_options', []),
    cython_args: cython_args,
    subdir: 'sklearn/cluster',
    install: true
  )
endforeach

subdir('_hdbscan')
