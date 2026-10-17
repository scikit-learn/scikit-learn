- Fixed :class:`ensemble.HistGradientBoostingClassifier` and
  :class:`ensemble.HistGradientBoostingRegressor` when a categorical feature is
  not placed before all numerical features: `interaction_cst` was applied to the
  wrong features, :func:`inspection.partial_dependence` with `method="recursion"`
  was computed for the wrong features, and the error for a categorical feature
  with too many categories reported the wrong feature index.
  By :user:`Arthur Lacote <cakedev0>`.
