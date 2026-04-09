_binary_tree_pxi = custom_target(
  '_binary_tree_pxi',
  output: '_binary_tree.pxi',
  input: '_binary_tree.pxi.tp',
  command: [tempita, '@INPUT@', '-o', '@OUTDIR@'],
)

# .pyx is generated so this is needed to make Cython compilation work. The pxi
# file is included avoid "missing dependency paths" with ninja -t missindeps
neighbors_cython_tree = [
  fs.copyfile('__init__.py'),
  fs.copyfile('_partition_nodes.pxd'),
  _binary_tree_pxi,
]

name_list = ['_ball_tree', '_kd_tree']

foreach name: name_list
  pyx = custom_target(
    name + '_pyx',
    output: name + '.pyx',
    input: name + '.pyx.tp',
    command: [tempita, '@INPUT@', '-o', '@OUTDIR@'],
    # TODO in principle this should go in py.exension_module below. This is
    # temporary work-around for dependency issue with .pyx.tp files. For more
    # details, see https://github.com/mesonbuild/meson/issues/13212
    depends: [neighbors_cython_tree, utils_cython_tree, metrics_cython_tree],
  )
  py.extension_module(
    name,
    cython_gen.process(pyx),
    dependencies: [np_dep],
    subdir: 'sklearn/neighbors',
    install: true
)
endforeach

neighbors_extension_metadata = {
  '_partition_nodes':
      {'sources': [cython_gen_cpp.process('_partition_nodes.pyx')],
       'dependencies': [np_dep]},
  '_quad_tree': {'sources': [cython_gen.process('_quad_tree.pyx')], 'dependencies': [np_dep]},
}

foreach ext_name, ext_dict : neighbors_extension_metadata
  py.extension_module(
    ext_name,
    [ext_dict.get('sources'), utils_cython_tree],
    dependencies: ext_dict.get('dependencies'),
    subdir: 'sklearn/neighbors',
    install: true
  )
endforeach
