fs = import('fs')

# Platform detection
is_windows = host_machine.system() == 'windows'
is_mingw = is_windows and cc.get_id() == 'gcc'

cython_c_args = []
if is_windows
  # For mingw-w64, link statically against the UCRT.
  gcc_link_args = ['-lucrt', '-static']
  if is_mingw
    add_project_link_arguments(gcc_link_args, language: ['c', 'cpp'])
    # Force gcc to float64 long doubles for compatibility with MSVC
    # builds, for C only.
    add_project_arguments('-mlong-double-64', language: 'c')
    # Make fprintf("%zd") work (see https://github.com/rgommers/scipy/issues/118)
    add_project_arguments('-D__USE_MINGW_ANSI_STDIO=1', language: ['c', 'cpp'])
    # Manual add of MS_WIN64 macro when not using MSVC.
    # https://bugs.python.org/issue28267
    bitness = run_command(
      '_build_utils/gcc_build_bitness.py',
      check: true
    ).stdout().strip()
    if bitness == '64'
      add_project_arguments('-DMS_WIN64', language: ['c', 'cpp'])
    endif
    # Silence warnings emitted by PyOS_snprintf for (%zd), see
    # https://github.com/rgommers/scipy/issues/118.
    # Use as c_args for extensions containing Cython code
    cython_c_args += ['-Wno-format-extra-args', '-Wno-format']
  endif
endif

# TODO Taken from scipy, seems like each project has its own tweak on this
# NumPy include directory - needed in all submodules
# The chdir is needed because within numpy there's an `import signal`
# statement, and we don't want that to pick up scipy's signal module rather
# than the stdlib module. The try-except is needed because when things are
# split across drives on Windows, there is no relative path and an exception
# gets raised. There may be other such cases, so add a catch-all and switch to
# an absolute path. Relative paths are needed when for example a virtualenv is
# placed inside the source tree; Meson rejects absolute paths to places inside
# the source tree.
# For cross-compilation it is often not possible to run the Python interpreter
# in order to retrieve numpy's include directory. It can be specified in the
# cross file instead:
#   [properties]
#   numpy-include-dir = /abspath/to/host-pythons/site-packages/numpy/core/include
#
# This uses the path as is, and avoids running the interpreter.
incdir_numpy = meson.get_external_property('numpy-include-dir', 'not-given')
if incdir_numpy == 'not-given'
  incdir_numpy = run_command(py,
    [
      '-c',
      'import os; os.chdir(".."); import numpy; print(numpy.get_include())'
    ],
    check: true
  ).stdout().strip()
endif

inc_np = include_directories(incdir_numpy)
np_dep = declare_dependency(include_directories: inc_np)

cc = meson.get_compiler('c')

# Deal with M_PI & friends; add `use_math_defines` to c_args
# Cython doesn't always get this correctly itself
# explicitly add the define as a compiler flag for Cython-generated code.
if is_windows
  use_math_defines = ['-D_USE_MATH_DEFINES']
else
  use_math_defines = []
endif

# Don't use the deprecated NumPy C API. Define this to a fixed version instead of
# NPY_API_VERSION in order not to break compilation for released SciPy versions
# when NumPy introduces a new deprecation. Use in a meson.build file::
#
#   py3.extension_module('_name',
#     'source_fname',
#     numpy_nodepr_api)
#
numpy_nodepr_api = '-DNPY_NO_DEPRECATED_API=NPY_1_9_API_VERSION'

# python_sources = [
#   '__init__.py',
#   '__init__.pyi',
#   'conftest.py',
# ]

# py.install_sources(
#   python_sources,
#   pure: false,
#   subdir: 'skimage'
# )

# TODO what is the sklearn equivalent
# cython_cli = find_program('_build_utils/cythoner.py')

# cython_gen = generator(cython_cli,
#   arguments : ['@INPUT@', '@OUTPUT@'],
#   output : '@BASENAME@.c')

# cython_gen_cpp = generator(cython_cli,
#   arguments : ['@INPUT@', '@OUTPUT@', '--cplus'],
#   output : '@BASENAME@.cpp')

c_undefined_ok = ['-Wno-maybe-uninitialized']

# Suppress warning for deprecated Numpy API.
# (Suppress warning messages emitted by #warning directives).
# Replace with numpy_nodepr_api after Cython 3.0 is out
# TODO What to do about this?
_cpp_Wno_cpp = cpp.get_supported_arguments('-Wno-cpp')
cython_c_args += [_cpp_Wno_cpp, use_math_defines]
cython_cpp_args = cython_c_args

extensions = ['_isotonic']

py.extension_module(
  '_isotonic',
  '_isotonic.pyx',
  # TODO what about additional arguments like below
  # install: false,
  #   c_args: [
  #     '-DNPY_NO_DEPRECATED_API=0',  # Cython still uses old NumPy C API
  #     # Require 1.25+ to test datetime additions
  #     '-DNPY_TARGET_VERSION=NPY_2_0_API_VERSION',
  #   ],
  #   include_directories: [npy_include_path],
  cython_args: cython_args,
  install: true,
  subdir: 'sklearn',
)

# TODO Needed for Cython?
__init__py = fs.copyfile('__init__.py')

sklearn_dir = py.get_install_dir() / 'sklearn'

subdir('__check_build')
subdir('_loss')
# utils needs to be quite early since used by other modules
subdir('utils')
# TODO metrics needs to be to be before cluster since it depends on it. Is there a better
# way to do this in meson?
subdir('metrics')
subdir('cluster')
subdir('datasets')
subdir('decomposition')
subdir('ensemble')
subdir('feature_extraction')
subdir('linear_model')
subdir('manifold')
subdir('preprocessing')
subdir('neighbors')
subdir('svm')
subdir('tree')
