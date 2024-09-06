hist_gradient_boosting_extension_metadata = {
  '_gradient_boosting': {'sources': ['_gradient_boosting.pyx'], 'dependencies': [openmp_dep]},
  'histogram': {'sources': ['histogram.pyx'], 'dependencies': [openmp_dep]},
  'splitting': {'sources': ['splitting.pyx'], 'dependencies': [openmp_dep]},
  '_binning': {'sources': ['_binning.pyx'], 'dependencies': [openmp_dep]},
  '_predictor': {'sources': ['_predictor.pyx'], 'dependencies': [openmp_dep]},
  '_bitset': {'sources': ['_bitset.pyx']},
  'common': {'sources': ['common.pyx']},
}

foreach ext_name, ext_dict : hist_gradient_boosting_extension_metadata
  py.extension_module(
    ext_name,
    ext_dict.get('sources'),
    dependencies: ext_dict.get('dependencies', []),
    cython_args: cython_args,
    subdir: 'sklearn/ensemble/_hist_gradient_boosting',
    install: true
  )
endforeach
