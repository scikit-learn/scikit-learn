libs_sources = {
  '_tree': {'sources': ['_tree.pyx'], 'override_options': ['cython_language=cpp']},
  '_splitter': {'sources': ['_splitter.pyx']},
  '_criterion': {'sources': ['_criterion.pyx']},
  '_utils': {'sources': ['_utils.pyx']},
}

foreach ext_name, ext_dict : libs_sources
    py.extension_module(
        ext_name,
        ext_dict.get('sources'),
        dependencies: [np_dep],
        override_options : ext_dict.get('override_options', []),
        cython_args: cython_args,
        c_args: ['-O3'],
        cpp_args: ['-O3'],
        subdir: 'sklearn/tree',
        install: true
    )
endforeach
