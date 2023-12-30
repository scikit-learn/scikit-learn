libs_sources = {
  '_tree': {'sources': ['_tree.pyx'], 'override_options': ['cython_language=cpp'], 'optimization_level': 'O3'},
  '_splitter': {'sources': ['_splitter.pyx'], 'optimization_level': 'O3'},
  '_criterion': {'sources': ['_criterion.pyx'], 'optimization_level': 'O3'},
  '_utils': {'sources': ['_utils.pyx'], 'optimization_level': 'O3'},
}

# TODO how to pass optimization_level 03?
foreach ext_name, ext_dict : libs_sources
    py.extension_module(
        ext_name,
        ext_dict.get('sources'),
        dependencies: [np_dep],
        override_options : ext_dict.get('override_options', []),
        cython_args: cython_args,
        subdir: 'sklearn/tree',
        install: true
    )
endforeach
