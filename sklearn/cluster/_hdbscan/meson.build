cluster_hdbscan_extension_metadata = {
  '_linkage': {'sources': [cython_gen.process('_linkage.pyx'), metrics_cython_tree]},
  '_reachability': {'sources': [cython_gen.process('_reachability.pyx')]},
  '_tree': {'sources': [cython_gen.process('_tree.pyx')]}
}

foreach ext_name, ext_dict : cluster_hdbscan_extension_metadata
  py.extension_module(
    ext_name,
    ext_dict.get('sources'),
    dependencies: [np_dep],
    subdir: 'sklearn/cluster/_hdbscan',
    install: true
  )
endforeach
