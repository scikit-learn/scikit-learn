_binary_tree_pxi = custom_target(
  '_binary_tree_pxi',
  output: '_binary_tree.pxi',
  input: '_binary_tree.pxi.tp',
  command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
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
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@']
  )
  py.extension_module(
    name,
    [pyx, neighbors_cython_tree, utils_cython_tree],
    dependencies: [np_dep],
    cython_args: cython_args,
    subdir: 'sklearn/neighbors',
    install: true
)
endforeach

neighbors_extension_metadata = {
  '_partition_nodes':
      {'sources': ['_partition_nodes.pyx'],
       'override_options': ['cython_language=cpp'], 'dependencies': [np_dep]},
  '_quad_tree': {'sources': ['_quad_tree.pyx'], 'dependencies': [np_dep]},
}

foreach ext_name, ext_dict : neighbors_extension_metadata
  py.extension_module(
    ext_name,
    [ext_dict.get('sources'), utils_cython_tree],
    dependencies: ext_dict.get('dependencies'),
    override_options : ext_dict.get('override_options', []),
    cython_args: cython_args,
    subdir: 'sklearn/neighbors',
    install: true
  )
endforeach
