# Check if scipy's cython_blas exports blas_int (ILP64 support).
scipy_has_blas_int = cython.compiles(
  'from scipy.linalg.cython_blas cimport blas_int',
  name: 'scipy cython_blas blas_int'
)

_blas_int_conf = configuration_data()
if scipy_has_blas_int
  _blas_int_conf.set('BLAS_INT_DEF', 'from scipy.linalg.cython_blas cimport blas_int')
else
  _blas_int_conf.set('BLAS_INT_DEF', 'ctypedef int blas_int')
endif

_blas_int_pxi = configure_file(
  input: '_blas_int.pxi.in',
  output: '_blas_int.pxi',
  configuration: _blas_int_conf,
)

# utils is cimported from other subpackages so this is needed for the cimport
# to work
utils_cython_tree = [
  # We add sklearn_root_cython_tree to make sure sklearn/__init__.py is copied
  # early in the build
  sklearn_root_cython_tree,
  fs.copyfile('__init__.py'),
  fs.copyfile('_bitset.pxd'),
  fs.copyfile('_cython_blas.pxd'),
  fs.copyfile('_heap.pxd'),
  fs.copyfile('_openmp_helpers.pxd'),
  fs.copyfile('_random.pxd'),
  fs.copyfile('_sorting.pxd'),
  fs.copyfile('_typedefs.pxd'),
  fs.copyfile('_vector_sentinel.pxd'),
]

# EXPERIMENTAL: SIMD-accelerated simultaneous_sort (see _sorting.pyx), vendored
# from https://github.com/numpy/x86-simd-sort as a meson subproject, built as a
# static library and linked directly into the _sorting extension module.
#
# x86-simd-sort requires x86_64 (AVX2/AVX-512, runtime-dispatched by the
# library itself); on any other architecture, or if the subproject can't be
# fetched (e.g. no network access) and -Dsimd_sort wasn't force-enabled, we
# fall back to the plain introsort path. `xss_available` and the generated
# _simd_sort_available.pxi are what _sorting.pyx uses to decide, at
# Cython-compile-time, whether the SIMD dispatch code is even emitted.
simd_sort_opt = get_option('simd_sort')

if simd_sort_opt.enabled() and host_machine.cpu_family() != 'x86_64'
  error(
    'simd_sort was force-enabled (-Dsimd_sort=enabled) but the vendored ' +
    'x86-simd-sort subproject only supports x86_64 (host is: ' +
    host_machine.cpu_family() + ').'
  )
endif

xss_available = simd_sort_opt.allowed() and host_machine.cpu_family() == 'x86_64'
if xss_available
  xss_proj = subproject(
    'x86-simd-sort',
    default_options: [
      'lib_type=static',
      'build_tests=false',
      'build_benchmarks=false',
      'use_openmp=false',
    ],
    # required=true (enabled) turns a fetch/build failure into a hard error
    # instead of a silent fallback -- useful to assert in CI/release builds
    # that the accelerated path is really the one being shipped.
    required: simd_sort_opt.enabled(),
  )
  xss_available = xss_proj.found()
endif

message(
  'SIMD-accelerated tree sort (x86-simd-sort): ' +
  (xss_available ? 'enabled' : 'disabled')
)

if xss_available
  xss_dep = xss_proj.get_variable('x86simdsortcpp_dep')
else
  xss_dep = declare_dependency()
endif

_simd_sort_conf = configuration_data()
_simd_sort_conf.set('HAVE_XSS_SIMD_SORT', xss_available ? 'True' : 'False')
_simd_sort_pxi = configure_file(
  input: '_simd_sort_available.pxi.in',
  output: '_simd_sort_available.pxi',
  configuration: _simd_sort_conf,
)

_sorting_sources = [cython_gen.process('_sorting.pyx'), _simd_sort_pxi]
if xss_available
  _sorting_sources += ['src' / 'xss_shim.cpp']
endif

utils_extension_metadata = {
  'sparsefuncs_fast':
    {'sources': [cython_gen.process('sparsefuncs_fast.pyx')]},
  '_cython_blas': {'sources': [cython_gen.process('_cython_blas.pyx'), _blas_int_pxi]},
  'arrayfuncs': {'sources': [cython_gen.process('arrayfuncs.pyx')]},
  '_bitset': {'sources': [cython_gen.process('_bitset.pyx')]},
  'murmurhash': {
      'sources': [cython_gen.process('murmurhash.pyx'), 'src' / 'MurmurHash3.cpp'],
  },
  '_fast_dict':
    {'sources': [cython_gen_cpp.process('_fast_dict.pyx')]},
  '_openmp_helpers': {'sources': [cython_gen.process('_openmp_helpers.pyx')], 'dependencies': [openmp_dep]},
  '_random': {'sources': [cython_gen.process('_random.pyx')]},
  '_typedefs': {'sources': [cython_gen.process('_typedefs.pyx')]},
  '_sorting': {
    'sources': _sorting_sources,
    'dependencies': [xss_dep],
    # x86simdsort.h uses C++17 syntax (e.g. std::function CTAD); sklearn's
    # own project-wide cpp_std default (C++14) doesn't apply to it.
    'override_options': xss_available ? ['cpp_std=c++17'] : [],
  },
  '_vector_sentinel':
    {'sources': [cython_gen_cpp.process('_vector_sentinel.pyx')],
     'dependencies': [np_dep]},
  '_isfinite': {'sources': [cython_gen.process('_isfinite.pyx')]},
}

foreach ext_name, ext_dict : utils_extension_metadata
  py.extension_module(
    ext_name,
    [ext_dict.get('sources'), utils_cython_tree],
    dependencies: ext_dict.get('dependencies', []),
    override_options: ext_dict.get('override_options', []),
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
    command: [tempita, '@INPUT@', '-o', '@OUTDIR@'],
  )
  utils_cython_tree += [pxd]

  pyx = custom_target(
    name + '_pyx',
    output: name + '.pyx',
    input: name + '.pyx.tp',
    command: [tempita, '@INPUT@', '-o', '@OUTDIR@'],
    # TODO in principle this should go in py.extension_module below. This is
    # temporary work-around for dependency issue with .pyx.tp files. For more
    # details, see https://github.com/mesonbuild/meson/issues/13212
    depends: [pxd, utils_cython_tree],
  )
  py.extension_module(
    name,
    cython_gen.process(pyx),
    subdir: 'sklearn/utils',
    install: true
   )
endforeach
