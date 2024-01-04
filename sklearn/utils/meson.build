utils_cython_tree = [
  fs.copyfile('__init__.py'),
  fs.copyfile('_cython_blas.pxd'),
  fs.copyfile('_heap.pxd'),
  fs.copyfile('_openmp_helpers.pxd'),
  fs.copyfile('_random.pxd'),
  fs.copyfile('_sorting.pxd'),
  fs.copyfile('_typedefs.pxd'),
  fs.copyfile('_vector_sentinel.pxd'),
]

libs_sources = {
  'sparsefuncs_fast': {'sources': ['sparsefuncs_fast.pyx'], 'dependencies': [np_dep]},
  '_cython_blas': {'sources': ['_cython_blas.pyx']},
  'arrayfuncs': {'sources': ['arrayfuncs.pyx']},
  'murmurhash': {
      'sources': ['murmurhash.pyx', 'src' / 'MurmurHash3.cpp'],
      'dependencies': [np_dep]
  },
  '_fast_dict': {'sources': ['_fast_dict.pyx'], 'override_options': ['cython_language=cpp']},
  '_openmp_helpers': {'sources': ['_openmp_helpers.pyx']},
  '_random': {'sources': ['_random.pyx'], 'dependencies': [np_dep]},
  '_typedefs': {'sources': ['_typedefs.pyx']},
  '_heap': {'sources': ['_heap.pyx']},
  '_sorting': {'sources': ['_sorting.pyx']},
  '_vector_sentinel':
    {'sources': ['_vector_sentinel.pyx'], 'override_options': ['cython_language=cpp'], 'dependencies': [np_dep]},
  '_isfinite': {'sources': ['_isfinite.pyx']},
}

foreach ext_name, ext_dict : libs_sources
  py.extension_module(
    ext_name,
    ext_dict.get('sources') + utils_cython_tree,
    dependencies: [np_dep, openmp_dep],
    override_options : ext_dict.get('override_options', []),
    cython_args: cython_args,
    subdir: 'sklearn/utils',
    install: true
  )
endforeach

names = ['_seq_dataset', '_weight_vector']

foreach name: names
  pxd = custom_target(
    name + '_pxd',
    output: name + '.pxd',
    input: name + '.pxd.tp',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
    install: true,
    install_dir: sklearn_dir / 'utils'
  )
  utils_cython_tree += [pxd]

  pyx = custom_target(
    name + '_pyx',
    output: name + '.pyx',
    input: name + '.pyx.tp',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@']
  )
  py.extension_module(
    name,
    [pxd, pyx] + utils_cython_tree,
    dependencies: [np_dep],
    cython_args: cython_args,
    subdir: 'sklearn/utils',
    install: true
   )
endforeach
