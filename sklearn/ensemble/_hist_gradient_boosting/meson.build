# We add sklearn_root_cython_tree and __init__.py so Cython can detect the
# package hierarchy and set the correct __module__ on extension types.
hist_gradient_boosting_cython_tree = [
  sklearn_root_cython_tree,
  fs.copyfile('__init__.py'),
]

hist_gradient_boosting_extension_metadata = {
  '_gradient_boosting': {'sources': [cython_gen.process('_gradient_boosting.pyx')],
                         'dependencies': [openmp_dep]},
  'histogram': {'sources': [cython_gen.process('histogram.pyx')], 'dependencies': [openmp_dep]},
  'splitting': {'sources': [cython_gen.process('splitting.pyx')], 'dependencies': [openmp_dep]},
  '_binning': {'sources': [cython_gen.process('_binning.pyx')], 'dependencies': [openmp_dep]},
  '_predictor': {'sources': [cython_gen.process('_predictor.pyx')], 'dependencies': [openmp_dep]},
  '_bitset': {'sources': [cython_gen.process('_bitset.pyx')]},
  'common': {'sources': [cython_gen.process('common.pyx')]},
}

foreach ext_name, ext_dict : hist_gradient_boosting_extension_metadata
  py.extension_module(
    ext_name,
    [ext_dict.get('sources'), hist_gradient_boosting_cython_tree],
    dependencies: ext_dict.get('dependencies', []),
    subdir: 'sklearn/ensemble/_hist_gradient_boosting',
    install: true
  )
endforeach
