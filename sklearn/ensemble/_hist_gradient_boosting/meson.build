hist_gradient_boosting_extension_metadata = {
  '_gradient_boosting': {'sources': ['_gradient_boosting.pyx']},
  'histogram': {'sources': ['histogram.pyx']},
  'splitting': {'sources': ['splitting.pyx']},
  '_binning': {'sources': ['_binning.pyx']},
  '_predictor': {'sources': ['_predictor.pyx']},
  '_bitset': {'sources': ['_bitset.pyx']},
  'common': {'sources': ['common.pyx']},
}

foreach ext_name, ext_dict : hist_gradient_boosting_extension_metadata
  py.extension_module(
    ext_name,
    ext_dict.get('sources'),
    dependencies: [openmp_dep],
    cython_args: cython_args,
    subdir: 'sklearn/ensemble/_hist_gradient_boosting',
    install: true
  )
endforeach
