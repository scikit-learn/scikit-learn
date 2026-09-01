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
    # `override_options` above is not enough: distributions such as
    # conda-forge inject their own `-O2` through CFLAGS/CXXFLAGS, which is
    # appended after Meson's own optimization flag on the compiler command
    # line, so it silently wins (GCC/Clang use the last -O flag). Per-target
    # c_args/cpp_args are added after environment flags by Meson, so this is
    # what actually guarantees -O3 for these hot loops regardless of the
    # build environment. See https://github.com/mesonbuild/meson/issues/15055.
    c_args: ['-O3'],
    cpp_args: ['-O3'],
    subdir: 'sklearn/tree',
    install: true
  )
endforeach
