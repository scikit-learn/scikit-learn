cython_program = find_program(cython.cmd_array()[0])

cython_gen = generator(cython_program,
  arguments : cython_args + ['@INPUT@', '--output-file', '@OUTPUT@'],
  output : '@BASENAME@.c',
)

py.extension_module(
  '_check_build',
  cython_gen.process('_check_build.pyx'),
  # cython_args: cython_args,
  install: true,
  subdir: 'sklearn/__check_build',
)
