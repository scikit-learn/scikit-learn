# Criterion computations (impurity, sums of weighted y) don't need to handle
# NaNs/Infs or strict IEEE 754 semantics, so relax float math there for speed.
_criterion_c_args = cc.get_supported_arguments('-ffast-math')

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
     'override_options': ['optimization=3'],
     'c_args': _criterion_c_args},
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
    c_args: ext_dict.get('c_args', []),
    subdir: 'sklearn/tree',
    install: true
  )
endforeach
