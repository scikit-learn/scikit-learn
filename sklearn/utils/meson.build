pxd_to_copy = ['_openmp_helpers.pxd', '_heap.pxd', '_sorting.pxd', '_vector_sentinel.pxd']
foreach pxd: pxd_to_copy
  fs.copyfile(pxd)
endforeach

py.extension_module(
  'sparsefuncs_fast',
  'sparsefuncs_fast.pyx',
  dependencies: [np_dep],
  subdir: 'sklearn/utils',
  install: true
)

fs.copyfile('_cython_blas.pxd')
py.extension_module(
  '_cython_blas',
  '_cython_blas.pyx',
  dependencies: [np_dep],
  subdir: 'sklearn/utils',
  install: true
)

fs.copyfile('_typedefs.pxd')
py.extension_module(
  '_typedefs',
  '_typedefs.pyx',
  subdir: 'sklearn/utils',
  install: true
)
