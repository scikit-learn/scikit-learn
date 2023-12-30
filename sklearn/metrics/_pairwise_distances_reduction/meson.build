# libs_sources = {
#   '_dataset_pair': [
#         {
#             'sources': ['_datasets_pair.pyx.tp', '_datasets_pair.pxd.tp'],
#             'language': 'c++',
#             'include_np': True,
#             'extra_compile_args': ['-std=c++11'],
#         },
#         {
#             'sources': ['_middle_term_computer.pyx.tp', '_middle_term_computer.pxd.tp'],
#             'language': 'c++',
#             'extra_compile_args': ['-std=c++11'],
#         },
#         {
#             'sources': ['_base.pyx.tp', '_base.pxd.tp'],
#             'language': 'c++',
#             'include_np': True,
#             'extra_compile_args': ['-std=c++11'],
#         },
#         {
#             'sources': ['_argkmin.pyx.tp', '_argkmin.pxd.tp'],
#             'language': 'c++',
#             'include_np': True,
#             'extra_compile_args': ['-std=c++11'],
#         },
#         {
#             'sources': ['_argkmin_classmode.pyx.tp'],
#             'language': 'c++',
#             'include_np': True,
#             'extra_compile_args': ['-std=c++11'],
#         },
#         {
#             'sources': ['_radius_neighbors.pyx.tp', '_radius_neighbors.pxd.tp'],
#             'language': 'c++',
#             'include_np': True,
#             'extra_compile_args': ['-std=c++11'],
#         },
#         {
#             'sources': ['_radius_neighbors_classmode.pyx.tp'],
#             'language': 'c++',
#             'include_np': True,
#             'extra_compile_args': ['-std=c++11'],
#         }
# }

# TODO Needed for relative imports in Cython?
fs.copyfile('__init__.py')

fs.copyfile('_classmode.pxd')

name_list_with_pxd = ['_datasets_pair', '_middle_term_computer', '_base', '_argkmin', '_radius_neighbors']

foreach name: name_list_with_pxd
  pxd = custom_target(
    name + '_pxd',
    output: name + '.pxd',
    input: name + '.pxd.tp',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@']
  )
  pyx = custom_target(
    name + '_pyx',
    output: name + '.pyx',
    input: name + '.pyx.tp',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@']
  )
  py.extension_module(
    name,
    [pxd, pyx],
    dependencies: [np_dep],
    # TODO add std=c++11 to cython_args I have c++14 by default in root meson.build so maybe not?
    override_options: ['cython_language=cpp'],
    subdir: 'sklearn/metrics/_pairwise_distances_reduction',
    install: true

)
endforeach

name_list_without_pxd = ['_argkmin_classmode', '_radius_neighbors_classmode']

foreach name: name_list_without_pxd
  pyx = custom_target(
    name + '_pyx',
    output: name + '.pyx',
    input: name + '.pyx.tp',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@']
  )
  py.extension_module(
    name,
    pyx,
    dependencies: [np_dep],
    # TODO add std=c++11 to cython_args I have c++14 by default in root meson.build so maybe not?
    override_options: ['cython_language=cpp'],
    subdir: 'sklearn/metrics/_pairwise_distances_reduction',
    install: true

)
endforeach
