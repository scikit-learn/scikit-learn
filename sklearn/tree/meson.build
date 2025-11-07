tree_extension_metadata = {
  '_tree':
    {'sources': [cython_gen_cpp.process('_tree.pyx')],
     'override_options': ['optimization=3']},
  '_splitter':
    {'sources': [cython_gen.process('_splitter.pyx')],
     'override_options': ['optimization=3']},
  '_partitioner':
    {'sources': [cython_gen.process('_partitioner.pyx')],
     'override_options': ['optimization=3']},
  '_criterion':
    {'sources': [cython_gen.process('_criterion.pyx')],
     'override_options': ['optimization=3']},
  '_utils':
    {'sources': [cython_gen.process('_utils.pyx')],
     'override_options': ['optimization=3']},
}

foreach ext_name, ext_dict : tree_extension_metadata
  py.extension_module(
    ext_name,
    [ext_dict.get('sources'), utils_cython_tree],
    dependencies: [np_dep],
    override_options : ext_dict.get('override_options', []),
    subdir: 'sklearn/tree',
    install: true
  )
endforeach
