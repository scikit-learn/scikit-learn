fs = import('fs')

cython_args = []

# Platform detection
is_windows = host_machine.system() == 'windows'
is_mingw = is_windows and cc.get_id() == 'gcc'

if is_mingw
  # For mingw-w64, link statically against the UCRT.
  gcc_link_args = ['-lucrt', '-static']
  add_project_link_arguments(gcc_link_args, language: ['c', 'cpp'])
  # Force gcc to float64 long doubles for compatibility with MSVC
  # builds, for C only.
  add_project_arguments('-mlong-double-64', language: 'c')
endif

# Adapted from scipy, each project seems to have its own tweaks for this. One
# day using dependency('numpy') will be a thing, see
# https://github.com/mesonbuild/meson/issues/9598.
# NumPy include directory - needed in all submodules
# Relative paths are needed when for example a virtualenv is
# placed inside the source tree; Meson rejects absolute paths to places inside
# the source tree. The try-except is needed because when things are split
# across drives on Windows, there is no relative path and an exception gets
# raised. There may be other such cases, so add a catch-all and switch to
# an absolute path.
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
      '''
import os
import numpy as np
try:
  incdir = os.path.relpath(np.get_include())
except Exception:
  incdir = np.get_include()
print(incdir)
'''
    ],
    check: true
  ).stdout().strip()
endif

inc_np = include_directories(incdir_numpy)
np_dep = declare_dependency(include_directories: inc_np)

openmp_dep = dependency('OpenMP', language: 'c', required: false)

if not openmp_dep.found()
    warning(
'''
                ***********
                * WARNING *
                ***********

It seems that scikit-learn cannot be built with OpenMP.

- Make sure you have followed the installation instructions:

    https://scikit-learn.org/dev/developers/advanced_installation.html

- If your compiler supports OpenMP but you still see this
  message, please submit a bug report at:

    https://github.com/scikit-learn/scikit-learn/issues

- The build will continue with OpenMP-based parallelism
  disabled. Note however that some estimators will run in
  sequential mode instead of leveraging thread-based
  parallelism.

                    ***
''')
endif

# Allow supporting older numpy versions than the version compiled against.
# Set the define to the minimum numpy supported version, see
# sklearn/_min_dependencies.py.
add_project_arguments(
  '-DNPY_TARGET_VERSION=NPY_1_19_API_VERSION', language : ['c', 'cpp'])

# TODO do the equivalent of sklearn/_build_utils/__init__.py that sets
# boundchecks=True if SKLEARN_ENABLE_DEBUG_CYTHON_DIRECTIVES is set to 1.
# There is no access to environment variables inside meson.build files so this
# needs to be a meson option or meson buildtype
scikit_learn_cython_args = [
  '-X language_level=3', '-X boundscheck=False', '-X wraparound=False',
  '-X initializedcheck=False', '-X nonecheck=False', '-X cdivision=True',
  '-X profile=False',
  # Needed for cython imports across subpackages, e.g. cluster pyx that
  # cimports metrics pxd
  '--include-dir', meson.global_build_root(),
]
cython_args += scikit_learn_cython_args

extensions = ['_isotonic']

py.extension_module(
  '_isotonic',
  '_isotonic.pyx',
  cython_args: cython_args,
  install: true,
  subdir: 'sklearn',
)

# Need for Cython cimports across subpackages to work, i.e. avoid errors like
# relative cimport from non-package directory is not allowed
fs.copyfile('__init__.py')

sklearn_dir = py.get_install_dir() / 'sklearn'

# Subpackages are mostly in alphabetical order except to handle Cython
# dependencies across subpackages
subdir('__check_build')
subdir('_loss')
# utils needs to be early since plenty of other modules cimports utils .pxd
subdir('utils')
# metrics needs to be to be before cluster since cluster cimports metrics .pxd
subdir('metrics')
subdir('cluster')
subdir('datasets')
subdir('decomposition')
subdir('ensemble')
subdir('feature_extraction')
subdir('linear_model')
subdir('manifold')
subdir('neighbors')
subdir('preprocessing')
subdir('svm')
subdir('tree')
