# .pyx is generated, so this is needed to make Cython compilation work
linear_model_cython_tree = [
  fs.copyfile('__init__.py'),
]

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
