cluster_extension_metadata = {
  '_dbscan_inner':
    {'sources': [cython_gen_cpp.process('_dbscan_inner.pyx')]},
  '_hierarchical_fast':
    {'sources': [cython_gen_cpp.process('_hierarchical_fast.pyx'), metrics_cython_tree]},
  '_k_means_common':
    {'sources': [cython_gen.process('_k_means_common.pyx')], 'dependencies': [openmp_dep]},
  '_k_means_lloyd':
    {'sources': [cython_gen.process('_k_means_lloyd.pyx')], 'dependencies': [openmp_dep]},
  '_k_means_elkan':
    {'sources': [cython_gen.process('_k_means_elkan.pyx')], 'dependencies': [openmp_dep]},
  '_k_means_minibatch':
    {'sources': [cython_gen.process('_k_means_minibatch.pyx')], 'dependencies': [openmp_dep]},
}

foreach ext_name, ext_dict : cluster_extension_metadata
  py.extension_module(
    ext_name,
    [ext_dict.get('sources'), utils_cython_tree],
    dependencies: [np_dep] + ext_dict.get('dependencies', []),
    subdir: 'sklearn/cluster',
    install: true
  )
endforeach

subdir('_hdbscan')
