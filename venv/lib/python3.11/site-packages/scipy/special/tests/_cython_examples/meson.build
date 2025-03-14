project('random-build-examples', 'c', 'cpp', 'cython')

fs = import('fs')

py3 = import('python').find_installation(pure: false)

cy = meson.get_compiler('cython')

if not cy.version().version_compare('>=3.0.8')
  error('tests requires Cython >= 3.0.8')
endif

py3.extension_module(
  'extending',
  'extending.pyx',
  install: false,
)

extending_cpp = fs.copyfile('extending.pyx', 'extending_cpp.pyx')
py3.extension_module(
  'extending_cpp',
  extending_cpp,
  install: false,
  override_options : ['cython_language=cpp']
)
