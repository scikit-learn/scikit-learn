libs_sources = {
  '_linkage': {'sources': ['_linkage.pyx']},
  '_reachability': {'sources': ['_reachability.pyx']},
  '_tree': {'sources': ['_tree.pyx']}
}

foreach ext_name, ext_dict : libs_sources
    py.extension_module(
        ext_name,
        ext_dict.get('sources'),
        dependencies: [np_dep],
        cython_args: cython_args,
        subdir: 'sklearn/cluster/_hdbscan',
        install: true
    )
endforeach
