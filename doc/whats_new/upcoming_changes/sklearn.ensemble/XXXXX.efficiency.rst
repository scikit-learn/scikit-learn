- Fitting :class:`ensemble.HistGradientBoostingRegressor` with
  `loss="quantile"` or `loss="absolute_error"` and no `sample_weight` is faster.
  The leaf values (the quantile or median of the residuals in each leaf) are now
  computed in parallel with OpenMP in Cython, instead of in a Python loop over the
  leaves. The fitted model is unchanged. Fit time drops by about 20-45%, most for
  datasets with many samples and few features.
  By :user:`Saswat <saswattulo>`
