tree_extension_metadata = {
  '_tree':
    {'sources': ['_tree.pyx'],
     'override_options': ['cython_language=cpp', 'optimization=3']},
  '_splitter':
    {'sources': ['_splitter.pyx'],
     'override_options': ['optimization=3']},
  '_criterion':
    {'sources': ['_criterion.pyx'],
     'override_options': ['optimization=3']},
  '_utils':
    {'sources': ['_utils.pyx'],
     'override_options': ['optimization=3']},
}

foreach ext_name, ext_dict : tree_extension_metadata
  py.extension_module(
    ext_name,
    [ext_dict.get('sources'), utils_cython_tree],
    dependencies: [np_dep],
    override_options : ext_dict.get('override_options', []),
    cython_args: cython_args,
    subdir: 'sklearn/tree',
    install: true
  )
endforeach
