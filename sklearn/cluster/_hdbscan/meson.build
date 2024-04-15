cluster_hdbscan_extension_metadata = {
  '_linkage': {'sources': ['_linkage.pyx', metrics_cython_tree]},
  '_reachability': {'sources': ['_reachability.pyx']},
  '_tree': {'sources': ['_tree.pyx']}
}

foreach ext_name, ext_dict : cluster_hdbscan_extension_metadata
  py.extension_module(
    ext_name,
    ext_dict.get('sources'),
    dependencies: [np_dep],
    cython_args: cython_args,
    subdir: 'sklearn/cluster/_hdbscan',
    install: true
  )
endforeach
