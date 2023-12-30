libs_sources = {
  'linkage': {'sources': ['_linkage.pyx']},
  '_reachability': {'sources': ['_reachability.pyx']},
  '_tree': {'sources': ['_tree.pyx']}
}

# TODO needed to use relative imports from sklearn/metrics otherwise error is:
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
        subdir: 'sklearn/cluster/_hdbscan',
        install: true
    )
endforeach
