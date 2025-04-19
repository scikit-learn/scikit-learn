py.extension_module(
  '_utils',
  ['_utils.pyx', utils_cython_tree],
  cython_args: cython_args,
  subdir: 'sklearn/mkl/_algo',
  install: true
)

name_list = [
  '_average_mkl',
  '_simple_mkl',
  '_sum_mkl',
]

foreach name: name_list
  py.extension_module(
    name,
    [name + '.pyx', utils_cython_tree],
    cython_args: cython_args,
    subdir: 'sklearn/mkl/_algo',
    install: true
  )
endforeach
