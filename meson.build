project(
  'scikit-learn',
  'c', 'cpp', 'cython',
  # TODO should the version be dynamically generated, scikit-image does it but
  # not scipy for some reason
  version: '1.5.0.dev0',
  license: 'BSD-3',
  meson_version: '>= 1.1.0',
  default_options: [
    'buildtype=debugoptimized',
    'c_std=c99',
    'cpp_std=c++14',
  ],
)

cc = meson.get_compiler('c')
cpp = meson.get_compiler('cpp')

# Check compiler is recent enough (see "Toolchain Roadmap" for details)
if cc.get_id() == 'gcc'
  if not cc.version().version_compare('>=8.0')
    error('scikit-image requires GCC >= 8.0')
  endif
elif cc.get_id() == 'msvc'
  if not cc.version().version_compare('>=19.20')
    error('scikit-image requires at least vc142 (default with Visual Studio 2019) ' + \
          'when building with MSVC')
  endif
endif

_global_c_args = cc.get_supported_arguments(
  '-Wno-unused-function',
)
add_project_arguments(_global_c_args, language: ['c', 'cpp'])

# We need -lm for all C code (assuming it uses math functions, which is safe to
# assume for scikit-image). For C++ it isn't needed, because libstdc++/libc++ is
# guaranteed to depend on it. For Fortran code, Meson already adds `-lm`.
m_dep = cc.find_library('m', required : false)
if m_dep.found()
  add_project_link_arguments('-lm', language : 'c')
endif

# Hide symbols when building on Linux with GCC. For Python extension modules,
# we only need `PyInit_*` to be public, anything else may cause problems. So we
# use a linker script to avoid exporting those symbols (this is in addition to
# Meson using `-fvisibility=hidden` for C and `-fvisibility-inlines-hidden` for
# C++ code. See gh-15996 for details.
# TODO: Should I do the same here (need simple link-version-pyinit.map)
# _linker_script = meson.project_source_root() / 'sklearn/_build_utils/link-version-pyinit.map'
# version_link_args = ['-Wl,--version-script=' + _linker_script]
# Note that FreeBSD only accepts version scripts when -shared is passed,
# hence we need to pass that to `cc.links` explicitly (flag is already
# present for `extension_module` invocations).
# if not cc.links('', name: '-Wl,--version-script', args: ['-shared', version_link_args])
#   version_link_args = []
# endif

cython = find_program('cython')
tempita = files('sklearn/_build_utils/tempita.py')

py = import('python').find_installation(pure: false)
py_dep = py.dependency()

subdir('sklearn')

install_subdir('sklearn', install_dir: py.get_install_dir())
