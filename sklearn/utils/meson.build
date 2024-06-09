# utils is cimported from other subpackages so this is needed for the cimport
# to work
utils_cython_tree = [
  # We add sklearn_root_cython_tree to make sure sklearn/__init__.py is copied
  # early in the build
  sklearn_root_cython_tree,
  fs.copyfile('__init__.py'),
  fs.copyfile('_cython_blas.pxd'),
  fs.copyfile('_heap.pxd'),
  fs.copyfile('_openmp_helpers.pxd'),
  fs.copyfile('_random.pxd'),
  fs.copyfile('_sorting.pxd'),
  fs.copyfile('_typedefs.pxd'),
  fs.copyfile('_vector_sentinel.pxd'),
]

utils_extension_metadata = {
  'sparsefuncs_fast':
    {'sources': ['sparsefuncs_fast.pyx']},
  '_cython_blas': {'sources': ['_cython_blas.pyx']},
  'arrayfuncs': {'sources': ['arrayfuncs.pyx']},
  'murmurhash': {
      'sources': ['murmurhash.pyx', 'src' / 'MurmurHash3.cpp'],
  },
  '_fast_dict':
    {'sources': ['_fast_dict.pyx'], 'override_options': ['cython_language=cpp']},
  '_openmp_helpers': {'sources': ['_openmp_helpers.pyx'], 'dependencies': [openmp_dep]},
  '_random': {'sources': ['_random.pyx']},
  '_typedefs': {'sources': ['_typedefs.pyx']},
  '_heap': {'sources': ['_heap.pyx']},
  '_sorting': {'sources': ['_sorting.pyx']},
  '_vector_sentinel':
    {'sources': ['_vector_sentinel.pyx'], 'override_options': ['cython_language=cpp'],
     'dependencies': [np_dep]},
  '_isfinite': {'sources': ['_isfinite.pyx']},
}

foreach ext_name, ext_dict : utils_extension_metadata
  py.extension_module(
    ext_name,
    [ext_dict.get('sources'), utils_cython_tree],
    dependencies: ext_dict.get('dependencies', []),
    override_options : ext_dict.get('override_options', []),
    cython_args: cython_args,
    subdir: 'sklearn/utils',
    install: true
  )
endforeach

util_extension_names = ['_seq_dataset', '_weight_vector']

foreach name: util_extension_names
  pxd = custom_target(
    name + '_pxd',
    output: name + '.pxd',
    input: name + '.pxd.tp',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
  )
  utils_cython_tree += [pxd]

  pyx = custom_target(
    name + '_pyx',
    output: name + '.pyx',
    input: name + '.pyx.tp',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
    # TODO in principle this should go in py.exension_module below. This is
    # temporary work-around for dependency issue with .pyx.tp files. For more
    # details, see https://github.com/mesonbuild/meson/issues/13212
    depends: [pxd, utils_cython_tree],
  )
  py.extension_module(
    name,
    pyx,
    cython_args: cython_args,
    subdir: 'sklearn/utils',
    install: true
   )
endforeach
